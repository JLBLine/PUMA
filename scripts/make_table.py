#!/usr/bin/python
import atpy
import numpy as np
import make_table_lib as mkl
import optparse
import copy

parser = optparse.OptionParser()

parser.add_option('-p', '--pref_cats',
	help='Enter names of matched cataloges, in order of preferable position')
	
parser.add_option('-m', '--matched_cats',
	help='Enter names of matched cataloges, in order of match')
	
parser.add_option('-c', '--cat_freqs', 
	help='Enter number of frequencies in each catalogue')
	
parser.add_option('-i', '--input_bayes', 
	help='Enter name of input matched bayes')

parser.add_option('-o', '--output_name', 
	help='Enter name for output catalogue')
	
parser.add_option('-a', '--prob_thresh',
	help='The lower and upper probability thresholds - separate with a comma')

parser.add_option('-e', '--epsilon_thresh', 
	help='Cut-off threshold for the epsilon residuals')

parser.add_option('-x', '--chi_thresh', 
	help='Cut-off threshold for the chi squared residuals')

parser.add_option('-r', '--resolution', 
	help='Resolution of base catalogue in "deg:arcmin:arcsec" ie "00:03:00" for 3arcmins')

parser.add_option('-v', '--verbose',action='store_true', 
	help='Add to have a verbose output of source statistcs')

parser.add_option('-s', '--split',default=0, 
	help='The resolution ("deg:arcmin:arcsec") over which to split combined sources')

options, args = parser.parse_args()

##Set up a bunch of initial parameters-----------------------------------------
##-----------------------------------------------------------------------------
low_prob,high_prob = map(float,options.prob_thresh.split(','))
jstat_thresh = float(options.epsilon_thresh)
chi_thresh = float(options.chi_thresh)

closeness = mkl.dec_to_deg(options.resolution)/2

cat_fs = options.cat_freqs.split(',')
cat_freqs= []
for fs in cat_fs:
	split = fs.split('~')
	split = map(float,split)
	cat_freqs.append(split)

matched_cats = options.matched_cats.split(',')
pref_cats = ['weighted']
for pref in options.pref_cats.split(','): pref_cats.append(pref)
num_freqs = []
for freq in cat_fs: num_freqs.append(len(freq.split('~')))

split = options.split
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------

##Transfer variables to the mkl module
mkl.closeness = closeness
mkl.high_prob = high_prob
mkl.low_prob = low_prob
mkl.chi_thresh = chi_thresh
mkl.jstat_thresh = jstat_thresh
mkl.num_freqs = num_freqs
mkl.split = split

##Lists for stats gathering
sources = []
sources_stats = []

rejected = []
rejected_stats = []

eyeballed = []
eyeballed_stats = []

not_treated = []
not_treated_stats = []

eyeball_outfile = open("%s-eyeball.txt" %options.output_name,'w+')
to_eyeball_comps = []
to_eyeball_accepts = []

accept_outfile = open("%s-accept.txt" %options.output_name,'w+')
to_accept_comps = []
to_accept_accepts = []

def make_entry(match,SI,intercept,SI_err,intercept_err,g_stats,accept_type,low_resids):
	'''Gather all of the information together in to a format that can easily be read in to create
	the output VOTable'''
	source = mkl.get_srcg(match)
	source.SI = SI
	source.intercept = intercept
	source.SI_err = SI_err
	source.intercept_err = intercept_err
	source.low_resids = low_resids
	sources.append(source)
	g_stats.accept_type = accept_type
	sources_stats.append(g_stats)

def make_accept(comp,g_stats,accept_type,accepted_inds):
	'''Write accepted source information to *output_name*-accept.txt to plot 
	the types of matches accepted'''
	a_ind = ""
	if len(accepted_inds)==0: accepted_inds = [0]
	for i in accepted_inds: a_ind+=str(i)+','
	a_ind = a_ind[:-1]
	if comp[0] == 'S':
		pass
	else:
		comp = comp[1:]
	##Stats used by plot_extended.py for plotting
	stat_string = "STATS %d %d %s accept %s\nEND_GROUP\n" %(g_stats.num_matches,g_stats.retained_matches,a_ind,accept_type)
	to_accept_comps.append(comp+stat_string)
	to_accept_accepts.append(stat_string)
	

def make_rejection(comp,g_stats,accept_type,accepted_inds):
	'''Formats info for the rejection group'''
	a_ind = ""
	for i in accepted_inds: a_ind+=str(i)+','
	a_ind = a_ind[:-1]
	if comp[0] == 'S':
		pass
	else:
		comp = comp[1:]
	stat_string = "STATS %d %d %s reject %s\nEND_GROUP\n" %(g_stats.num_matches,g_stats.retained_matches,a_ind,accept_type)
	to_eyeball_comps.append(comp+stat_string)
	to_eyeball_accepts.append(stat_string)
	rejected.append(1)
	g_stats.accept_type = accept_type
	rejected_stats.append(g_stats)
	
def make_eyeballed(comp,g_stats,accept_type,accepted_inds):
	'''Formats the info for the retained for inspection group'''
	a_ind = ""
	for i in accepted_inds: a_ind+=str(i)+','
	a_ind = a_ind[:-1]
	comp = comp[1:]
	stat_string = "STATS %d %d %s eyeball %s\nEND_GROUP\n" %(g_stats.num_matches,g_stats.retained_matches,a_ind,accept_type)
	to_eyeball_comps.append(comp+stat_string)
	to_eyeball_accepts.append(stat_string)
	eyeballed.append(1)
	g_stats.accept_type = accept_type
	eyeballed_stats.append(g_stats)
	
def single_match_test(src_all,comp,accepted_matches,accepted_inds,g_stats,num_matches,repeated_cats,matches):
	'''Takes a combination of sources, one from each catalogue, with positional probabilities,
	and determines whether they are a match or not - Algorithm 2 in the write up'''
	match = accepted_matches[0]
	prob = float(match[-1])
	
	##calculate_resids needs a list of matches - calculate parameters
	jstat_resids,params,bses,chi_resids = mkl.calculate_resids([match])

	##Gather source information
	src_g = mkl.get_srcg(match)
	
	##Play the prob trick again to work out which match has been accepted
	match_probs = [float(m.split()[-1]) for m in matches]
	dom_num = match_probs.index(prob)+1
	match_crit = "Combination (%d)\npossible\n%s repeated cats" %(dom_num,repeated_cats)
	
	##Check to see if all matched sources are within the closeness test - create an
	##error ellipse by combined closeness with base cat error
	##Need to convert closeness in to an RA offset, due to spherical trigonometry
	dr=np.pi/180.0
	delta_RA = np.arccos((np.cos(closeness*dr)-np.sin(src_all.decs[0]*dr)**2)/np.cos(src_all.decs[0]*dr)**2)/dr
	
	##Make a list of the ras and decs of the sources to distance test
	ras_test = [ra for ra in src_g.ras if ra!=-100000.0]
	dec_test = [dec for dec in src_g.decs if dec!=-100000.0]
	
	small_test = []
	for ra,dec in zip(ras_test,dec_test):
		##Even though at same dec, 3arcmis offset in RA isn't neccessarily 3arcmins arcdistance 
		ra_dist = mkl.arcdist(src_all.ras[0],ra,src_all.decs[0],src_all.decs[0])
		dec_dist = src_all.decs[0] - dec
		ra_axis = src_all.rerrs[0] + abs(delta_RA)
		dec_axis = src_all.derrs[0] + closeness

		##Test to see if the source lies with an error ellipse created using semi-major
		##and minor axes defined by the ra and dec error of the base cat + half the resolution
		##of the base cat (closeness)
		ell_test = (ra_dist/ra_axis)**2 + (dec_dist/dec_axis)**2
		if ell_test <= 1:
			small_test.append('yes')
		else:
			##Otherwise, fails
			small_test.append('no')
			#no_names.append(repeat_name)  #Note the name of the sources that are far away
	
	#Fail the positional test if a source is outside of the resolution plus position error
	close_test = 'passed'
	if 'no' in small_test: close_test = 'failed'
	
	##If prob is higher than threshold, ignore position of sources and accept the match
	if prob>high_prob:
		##Accept the source, put it in a way that can be read when constructing the final table
		if chi_resids[0]<=2:
			make_entry(match,params[0][0],params[0][1],bses[0][0],bses[0][1],g_stats,'position',0)
		else:
			make_entry(match,params[0][0],params[0][1],bses[0][0],bses[0][1],g_stats,'position',1)
			
		make_accept(comp,g_stats,'position',accepted_inds)
	else:
		##look to see if all sources are within the resolution of the
		##base catalogue or above some probability theshold, if so check with a spec test else reject them
		if close_test=='passed' or prob>low_prob:  
			##IF below eith threshold, append with the applicable fit label
			if jstat_resids[0]<=jstat_thresh or chi_resids[0]<=chi_thresh:
				if chi_resids[0]<=2:
					make_entry(match,params[0][0],params[0][1],bses[0][0],bses[0][1],g_stats,'spectral',0)
				else:
					make_entry(match,params[0][0],params[0][1],bses[0][0],bses[0][1],g_stats,'spectral',1)
				make_accept(comp,g_stats,'spectral',accepted_inds)
			else:
				g_stats.retained_matches = 1
				##Put accepted inds as [0] just to have something outputted to the investigate text file - 
				##accepted_inds is empty is rejecting at this stage
				make_rejection(comp,g_stats,'spectral',[0])
		else:
			g_stats.retained_matches = 1
			make_rejection(comp,g_stats,'position',[0])

##Open the input text file (output from calculate_bayes.py)
bayes_comp = open(options.input_bayes).read().split('END_GROUP')
del bayes_comp[-1]

for comp in bayes_comp:
	##Get the information into nice usable forms, and get rid of empty/pointless
	##entries
	chunks = comp.split('START_COMP')
	all_info = chunks[0].split('\n')
	
	##FOR SOME REASON CAN'T DO BOTH OF THESE LINES IN THE SAME FOR LOOP?!?!?!
	for entry in all_info:   
		if entry=='': del all_info[all_info.index(entry)]
	for entry in all_info:
		if 'START' in entry: del all_info[all_info.index(entry)]

	matches = chunks[1].split('\n')
	del matches[0],matches[-2:]
	
	##Put the information for every source in the matched group in one source_group() class
	##(see apply_criteria_lib for source_group())
	src_all = mkl.get_allinfo(all_info)
	
	##This line applies positional criteria, and tells us if a simple one catalogue repeat source, returning
	##how many combinations are possible and statistics on them
	repeated_cats,accepted_matches,accepted_inds,accepted_probs,jstats,chi_reds,g_stats = mkl.matches_retained(src_all,matches)
	match_crit = "%d of %d \ncombinations \npossible \n%s repeated cats" %(len(accepted_matches),len(matches),repeated_cats)
	
	##If no combinations are possible, reject all info (goes into the eyeball document)
	if len(accepted_matches)==0:
		cats = src_all.cats
		repeated_inds = [i for i in xrange(len(cats)) if cats.count(cats[i])>1]
		make_rejection(comp,g_stats,'position',repeated_inds)

	##If just one combo positionally possible, do a single combo check
	elif len(accepted_matches)==1:
		single_match_test(src_all,comp,accepted_matches,accepted_inds,g_stats,len(matches),repeated_cats,matches)
		
	##If more than one combination is positionally possible:
	else:
		##Check for a dominant source. The combination must be the only one with high position prob,
		##all others with low positional proability, and must dominate spectrally
		
		num_cat = len(set([cat for cat in src_all.cats if cat!='-100000.0']))
		
		dom_source = mkl.spec_pos_agree(jstats,chi_reds,accepted_probs,num_cat)
		src_g = mkl.get_srcg(accepted_matches[0])
		##If it finds a dominant source, accept it - counts as a spectral match 
		if dom_source!='none':
			jstat_resids,params,bses,chi_resids = mkl.calculate_resids([accepted_matches[dom_source]])
			if chi_resids[0]<=2:
				make_entry(accepted_matches[dom_source],params[0][0],params[0][1],bses[0][0],bses[0][1],g_stats,'spectral',0)
			else:
				make_entry(accepted_matches[dom_source],params[0][0],params[0][1],bses[0][0],bses[0][1],g_stats,'spectral',1)
			make_accept(comp,g_stats,'spectral',accepted_inds)
		##If nothing dominates, send to check if a combined source works
		else:
			comb_crit, comb_source, comb_jstat, comb_chi_red = mkl.combine_flux(src_all,src_g,accepted_inds,'plot=no',len(matches))
				
			if 'Accepted' in comb_crit:
				##If source was combined, add one new source with combine in the g_stat
				if len(comb_source)==1:
					sources.append(comb_source[0])
					g_stats.accept_type = 'combine'
					##Gather stats
					sources_stats.append(g_stats)
					##make_accept makes infomation for plot_extended to use
					make_accept(comp,g_stats,'combine',accepted_inds)
				##Otherwise, add as many components as were made with split function
				else:
					##Will only get here if splitting was turned on by the user
					letters = ['A','B','C','D','E','F']
					for source in comb_source:
						letter = letters[comb_source.index(source)]
						new_gstat = copy.deepcopy(g_stats)
						##Gather stats
						sources.append(source)
						new_gstat.accept_type = 'split%s' %letter
						sources_stats.append(new_gstat)
						##make_accept makes infomation for plot_extended to use
						make_accept(comp,g_stats,'split%s' %letter,accepted_inds)
			else:
				if 'split' in comb_crit:
					make_eyeballed(comp,g_stats,'split',accepted_inds)
				else:
					make_eyeballed(comp,g_stats,'combine',accepted_inds)

##Print information about matches where there was only one source from each catalogue
def print_singles():
	print '\nSINGLE MATCHES++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
	print "............Overall............"
	num_accept = len([1 for g_stats in sources_stats if g_stats.num_matches==1 and g_stats.accept_type!='splitB' and g_stats.accept_type!='splitC' and g_stats.accept_type!='splitD' and g_stats.accept_type!='splitE'and g_stats.accept_type!='splitF'])
	num_reject = len([1 for g_stats in rejected_stats if g_stats.num_matches==1])
	num_eyeball = len([1 for g_stats in eyeballed_stats if g_stats.num_matches==1])
	print "Singles: ", num_accept+num_reject+num_eyeball, " Accepted:", num_accept,' Rejected entirely:', num_reject
	print "...........Overall Positional............."
	print "Singles accepted purely by position: " , len([1 for g_stats in sources_stats if g_stats.num_matches==1 and g_stats.accept_type=='position'])
	print "Singles rejected purely by position: " , len([1 for g_stats in rejected_stats if g_stats.num_matches==1 and g_stats.accept_type=='position'])
	print "...........Overall Spectral............."
	print "Singles accepted at the spectral stage: " , len([1 for g_stats in sources_stats if g_stats.num_matches==1 and g_stats.accept_type=='spectral'])
	print "Singles rejected at the spectral stage: " , len([1 for g_stats in rejected_stats if g_stats.num_matches==1 and g_stats.accept_type=='spectral'])
	print "\n"

##Function to calculate the break down of statistics for a particular number of matches
def print_out(num_matches,description,BIG):
	print '%s MATCHES++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++' %BIG
	print "............Overall............"
	num_accept = len([1 for g_stats in sources_stats if g_stats.num_matches==num_matches and g_stats.accept_type!='splitB' and g_stats.accept_type!='splitC' and g_stats.accept_type!='splitD' and g_stats.accept_type!='splitE'and g_stats.accept_type!='splitF'])
	num_reject = len([1 for g_stats in rejected_stats if g_stats.num_matches==num_matches])
	num_eyeball = len([1 for g_stats in eyeballed_stats if g_stats.num_matches==num_matches])
	print "%s: " %description ,num_accept+num_reject+num_eyeball , " Accepted:", num_accept,' Rejected entirely:', num_reject
	print "\tRetained to investigate: ", num_eyeball
	print "...........Overall Positional............."
	print "%s accepted purely by position: " %description, len([1 for g_stats in sources_stats if g_stats.num_matches==num_matches and g_stats.accept_type=='position'])
	print "%s rejected purely by position: " %description, len([1 for g_stats in rejected_stats if g_stats.num_matches==num_matches and g_stats.accept_type=='position'])
	print "...........Overall Spectral............."
	print "%s accepted at the spectral stage: " %description, len([1 for g_stats in sources_stats if g_stats.num_matches==num_matches and g_stats.accept_type=='spectral'])
	print "%s rejected at the spectral stage: " %description, len([1 for g_stats in rejected_stats if g_stats.num_matches==num_matches and g_stats.accept_type=='spectral'])
	print "...........Overall Combine............."
	print "%s accepted at the combine stage: " %description, len([1 for g_stats in sources_stats if g_stats.num_matches==num_matches and g_stats.accept_type=='combine'])
	print "%s sent to eyeball at the combine stage: " %description, len([1 for g_stats in eyeballed_stats if g_stats.num_matches==num_matches and g_stats.accept_type=='combine'])
	print "...........Overall Split............."
	print "%s accepted at the splitting stage: " %description, len([1 for g_stats in sources_stats if g_stats.num_matches==num_matches and g_stats.accept_type=='splitA'])
	print "%s sent to eyeball at the splitting stage: " %description, len([1 for g_stats in eyeballed_stats if g_stats.num_matches==num_matches and g_stats.accept_type=='split'])
	print "-----------------------------------------------------------------------"
	for retained in reversed(range(1,num_matches+1)):
		print "...........%d Positionally Retained..............." %retained
		accept_num = len([1 for g_stats in sources_stats if g_stats.num_matches==num_matches and g_stats.retained_matches==retained and g_stats.accept_type!='splitB' and g_stats.accept_type!='splitC' and g_stats.accept_type!='splitD' and g_stats.accept_type!='splitE'and g_stats.accept_type!='splitF'])
		reject_num = len([1 for g_stats in rejected_stats if g_stats.num_matches==num_matches and g_stats.retained_matches==retained])
		combined_num = len([1 for g_stats in eyeballed_stats if g_stats.num_matches==num_matches and g_stats.retained_matches==retained])
		print "\ttotal: ", accept_num + reject_num + combined_num
		print "\taccepted purely by position: ", len([1 for g_stats in sources_stats if g_stats.num_matches==num_matches and g_stats.accept_type=='position' and g_stats.retained_matches==retained])
		print "\trejected purely by position: ", len([1 for g_stats in rejected_stats if g_stats.num_matches==num_matches and g_stats.accept_type=='position' and g_stats.retained_matches==retained])
		print "\taccepted at the spectral stage: ", len([1 for g_stats in sources_stats if g_stats.num_matches==num_matches and g_stats.accept_type=='spectral' and g_stats.retained_matches==retained])
		print "\trejected at the spectral stage: ", len([1 for g_stats in rejected_stats if g_stats.num_matches==num_matches and g_stats.accept_type=='spectral' and g_stats.retained_matches==retained])
		print "\taccepted at the combine stage: ", len([1 for g_stats in sources_stats if g_stats.num_matches==num_matches and g_stats.accept_type=='combine' and g_stats.retained_matches==retained])
		print "\tsent to eyeball at the combine stage: ", len([1 for g_stats in eyeballed_stats if g_stats.num_matches==num_matches and g_stats.accept_type=='combine' and g_stats.retained_matches==retained])
		print "\taccepted at the splitting stage: ", len([1 for g_stats in sources_stats if g_stats.num_matches==num_matches and g_stats.accept_type=='splitA' and g_stats.retained_matches==retained])
		print "\tsent to eyeball at the splitting stage: ", len([1 for g_stats in eyeballed_stats if g_stats.num_matches==num_matches and g_stats.accept_type=='split' and g_stats.retained_matches==retained])
	print '\n'

print '-------------------------------------------------------------------------'
print '+++++++++++TOTAL NUMBER OF SOURCES: %d +++++++++++++++++++++++++++++++' %len(bayes_comp)
print '-------------------------------------------------------------------------'
print "All sources accepted: " , len([1 for g_stats in sources_stats if g_stats.accept_type!='splitB' and g_stats.accept_type!='splitC' and g_stats.accept_type!='splitD' and g_stats.accept_type!='splitE'and g_stats.accept_type!='splitF'])
print "\taccepted by position: ", len([1 for g_stats in sources_stats if g_stats.accept_type=='position'])
print "\taccepted by spectral: ", len([1 for g_stats in sources_stats if g_stats.accept_type=='spectral'])
print "\taccepted by combine: ", len([1 for g_stats in sources_stats if g_stats.accept_type=='combine'])
print "\taccepted by splitting: ", len([1 for g_stats in sources_stats if g_stats.accept_type=='splitA'])
print "All sources rejected: " , len(rejected_stats)
print "\trejected by position: ", len([1 for g_stats in rejected_stats if g_stats.accept_type=='position'])
print "\trejected by spectral: ", len([1 for g_stats in rejected_stats if g_stats.accept_type=='spectral'])
print "All sources retained to eyeball: " , len(eyeballed_stats)
#print "\tretained by position: ", len([1 for g_stats in eyeballed_stats if g_stats.accept_type=='position'])
#print "\tretained by spectral: ", len([1 for g_stats in eyeballed_stats if g_stats.accept_type=='spectral'])
print "\tretained by combine: ", len([1 for g_stats in eyeballed_stats if g_stats.accept_type=='combine'])
print "\tretained by splitting: ", len([1 for g_stats in eyeballed_stats if g_stats.accept_type=='split'])

if options.verbose==True:
	print_singles()
	print_out(2,'Doubles','DOUBLE')
	print_out(3,'Triples','TRIPLE')
	print_out(4,'Quadruples','QUADRUPLES')
	##ADD THESE AND ANY MORE TO YOU HEARTS DESIRE
	#print_out(5,'Quintuples','QUINTUPLES')
	#print_out(6,'Sextuples','SEXTUPLES')

###WRITE THE TABLE+++++++++++++++++++++++++++++++++++++++++++++_____________________________------------------------++++++++++++++++++++++++++++++
###+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

t=atpy.Table(masked=True)

updated_ras = []
updated_decs = []
updated_rerrs = []
updated_derrs = []

for source in sources:
	##this is for choosing the most 'reliable' position as defined by
	##the input preference cats
	for p_cat in pref_cats:
		if p_cat in source.cats:
			pos_ind = source.cats.index(p_cat)
			break
	else:
		print 'no catalogue match in positional preferences - OMG'
	updated_ras.append(source.ras[pos_ind])
	updated_decs.append(source.decs[pos_ind])
	updated_rerrs.append(source.rerrs[pos_ind])
	updated_derrs.append(source.derrs[pos_ind])

##Make an 'MWA' name based on position (limited in number of characters due to the RTS)
##should change this really
def make_name(ra,dec,ext):
	ra_str = mkl.deg_to_hour(ra)[:6]
	dec_str = mkl.deg_to_degmins(dec)[:5]
	if ext=='':
		return 'MWA'+ra_str+dec_str
	else:
		return 'MWA'+ra_str+dec_str+ext


original_ras = [source.ras[0] for source in sources]
original_decs = [source.decs[0] for source in sources]
type_matches = [source.accept_type for source in sources_stats]

##Name the sources based on position. If the sources were the result of a split,
##name by the position of the base source, and add the letter that was defined by the
##split letter in the stat report object
names = []
for i in xrange(len(type_matches)):
	if 'split' in type_matches[i]:
		ext = type_matches[i][-1]
		name = make_name(original_ras[i],original_decs[i],ext)
	else:
		name = make_name(updated_ras[i],updated_decs[i],'')
	names.append(name)

t.add_column('Name',names,description='Name based on position of combined source')

##Add lots of data
prim_names = [source.names[0] for source in sources]
t.add_column('%s_name' %matched_cats[0],prim_names,description='Name of %s component' %matched_cats[0])
t.add_column('updated_RA_J2000',np.array(updated_ras),description='Updated Right ascension of source',unit='deg')
t.add_column('updated_DEC_J2000',np.array(updated_decs),description='Updated Declination of source',unit='deg')
t.add_column('updated_RA_err',np.array(updated_rerrs),description='Error on Updated Right ascension of source',unit='deg')
t.add_column('updated_DEC_err',np.array(updated_derrs),description='Error on Updated Declination of source',unit='deg')
t.add_column('original_RA_J2000',np.array(original_ras),description='Original Right ascension of source',unit='deg')
t.add_column('original_DEC_J2000',np.array(original_decs),description='Original Declination of source',unit='deg')
		
##For every catalogue in the match
for cat in xrange(len(num_freqs)):
	##See how many frequencies that source has
	num_freq = num_freqs[cat]
	##For every frequency, make a column of fluxes and flux errors, masking every value with -100000.0
	for freq in xrange(num_freq):
		fluxs = np.array([src.fluxs[cat][freq] for src in sources])
		ferrs = np.array([src.ferrs[cat][freq] for src in sources])
		t.add_column('S_%.1f' %float(cat_freqs[cat][freq]),fluxs,description='Flux at %.1fMHz' %float(cat_freqs[cat][freq]),mask=fluxs==-100000.0, fill=None,unit='Jy')
		t.add_column('e_S_%.1f' %float(cat_freqs[cat][freq]),ferrs,description='Flux error at %.1fMHz' %float(cat_freqs[cat][freq]),mask=ferrs==-100000.0, fill=None,unit='Jy')

##Extrapolate the base catalogue frequency flux density, using the fitted values for late comparison
extrap_freq = cat_freqs[0][0]
		
def extrap(freq,SI,intercept):
	return np.exp((np.log(freq)*SI)+intercept)
	
def extrap_error(freq,ext_flux,int_err,SI_err):
	denom = ext_flux**2
	numer = int_err**2 + ((np.log(freq)**2)*(SI_err)**2)
	return np.sqrt(numer/denom)
	
SIs = np.array([source.SI for source in sources])
SI_errs = np.array([source.SI_err for source in sources])
intercepts = np.array([source.intercept for source in sources])
intercept_errs = np.array([source.intercept_err for source in sources])

##create a mask for all of the errors reported for matched with only two cats - can't have an error on
##a fit to just two data points, so statsmodel spits out nonsenseCatalogue & Frequency (MHz) & Number of sources & Sky Coverage
SI_err_mask = [err==float('inf') or err==float('nan') for err in SI_errs]

t.add_column('SI',SIs,description='Spectral Index of Fit')
t.add_column('e_SI',SI_errs,description='Std error on Spectral Index of Fit',mask=SI_err_mask,fill=None)
t.add_column('Intercept',intercepts,description='Intercept of Fit')
t.add_column('e_Intercept',intercept_errs,description='Std error on Intercept of Fit',mask=SI_err_mask,fill=None)

extrap_base = np.array([extrap(extrap_freq,source.SI,source.intercept) for source in sources])
extrap_base_err = np.array([extrap_error(extrap_freq,extrap_flux,source.intercept_err,source.SI_err) for extrap_flux, source in zip(extrap_base,sources)])

t.add_column('S_%.1f_ext' %extrap_freq,extrap_base,unit='Jy',description='Flux at extrapolted to base catalogue frequency using fitted values')
t.add_column('e_S_%.1f_ext' %extrap_freq,extrap_base_err,unit='Jy',description='Error on flux extrapolated to base frequency using error on fitted values',mask=SI_err_mask,fill=None)

##See how many unique catalogues there are
num_cats = [len(set([cat for cat in source.cats if cat!='-100000.0'])) for source in sources]
t.add_column('Number_cats',np.array(num_cats),description='Number of matched catalogues')

numb_matches = [source.num_matches for source in sources_stats]
t.add_column('Number_matches',np.array(numb_matches),description='Number of possible combinations in the group match')

retained_matches = [source.retained_matches for source in sources_stats]
t.add_column('Retained_matches',np.array(retained_matches),description='Number of retained (by position) combinations in the group match')

type_matches = [source.accept_type for source in sources_stats]
t.add_column('Match_stage',np.array(type_matches),description='The PUMA stage at which a decision was made')

low_resids = [source.low_resids for source in sources]
t.add_column('Low_resids',np.array(low_resids),description='Whether the fitted data had residuals above the threshold of chi_reduced<=2. 1 is above threshold, 0 otherwise')

t.write("%s.vot" %options.output_name,votype='ascii',overwrite=True)

#WRITE THE TO EYEBALL+++++++++++++++++++++++++++++++++++++++++++++_____________________________------------------------++++++++++++++++++++++++++++++
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##Write out two files, one that contains all the information for all sources accepted,
##and one for the sources rejected/eyeballed - useful for later analysis, and used by
##plot_extended.py

##Find the base cat flux, and sort by descending flux so can investigate brightest sources to eyeball first
eyeball_fluxes = [float(comp.split('\n')[1].split()[7]) for comp in to_eyeball_comps]
to_eyeball_comps = [to_eye for flux,to_eye in sorted(zip(eyeball_fluxes,to_eyeball_comps), key=lambda pair: pair[0], reverse=True) ]
for eyeball in to_eyeball_comps:
	eyeball_outfile.write(eyeball)
eyeball_outfile.close()

##Same for the accepted sources
accept_fluxes = [float(comp.split('\n')[1].split()[7]) for comp in to_accept_comps]
to_accept_comps = [to_eye for flux,to_eye in sorted(zip(accept_fluxes,to_accept_comps), key=lambda pair: pair[0], reverse=True) ]
for accept in to_accept_comps:
	accept_outfile.write(accept)
accept_outfile.close()