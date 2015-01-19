#!/usr/bin/python

import atpy
import numpy as np
#from itertools import combinations
import optparse

##Get all of the input variables
parser = optparse.OptionParser()

parser.add_option('-p', '--primary_cat',
	help='Enter primary catalogue prefix')
	
parser.add_option('-g', '--primary_freq',
	help='Enter frequency of primary catalogue')
	
parser.add_option('-m', '--matched_cats',
	help='Enter names of matched cataloges, separated by commas')
	
parser.add_option('-f', '--matched_freqs',
	help='Enter frequencies of matched catalogues, separated by commas. If more than one frequency, separate by ~')
	
parser.add_option('-o', '--out_name', 
	help='Enter name of output matched group file')
	
options, args = parser.parse_args()

##Make all of the information usable
primary_cat = options.primary_cat
primary_freq = options.primary_freq
matched_cats = options.matched_cats.split(',')
matched_freqs = options.matched_freqs.split(',')

num_freqs = [len(primary_cat.split('~'))]
for freq in matched_freqs: num_freqs.append(len(freq.split('~')))
num_freqs = map(int,num_freqs)

out_name = options.out_name

dr = np.pi/180.0

##Opens each indivudual vot match table
def open_table(name):
	data = atpy.Table(name,verbose=False)
	shape = data.shape
	rows = shape[0]
	return data,rows

##A class to store the information for a group of matched sources
class source_info:
	def __init__(self):
		self.cats = []
		self.names = []
		self.ras = []
		self.rerrs = []
		self.decs = []
		self.derrs = []
		self.freqs = []
		self.fluxs = []
		self.ferrs = []
		self.majors = []
		self.minors = []
		self.PAs = []
		self.flags = []
		self.IDs = []

##A class to store the information for a single source
class source_single:
	def __init__(self):
		self.cat = []
		self.name = None
		self.ras = None
		self.rerr = None
		self.dec = None
		self.derr = None
		self.freqs = []
		self.fluxs = []
		self.ferrs = []
		self.major = None
		self.minor = None
		self.PA = None
		self.flag = None
		self.ID = None
		
##The way STILTS is run in cross_match.py means sources in the matched catalogues
##may be matched to more than one base catalogue source. Read through all the matches
##and check for any repeated matched sources; if they exist, chose the match with the 
##smallest angular separation. Delete the other matched row

##Will contain lists of rows to be skipped that correspond to each catalogue in matched_cats
skip_rows = []

for cat in matched_cats:
	match_names = []
	separations = []
	skip_row = []
	
	##Get data
	match_name = 'matched_%s_%s.vot' %(primary_cat,cat)
	data,rows = open_table(match_name)
	for row in data:
		match_names.append(str(row[12]))
		separations.append(float(row[-1]))
	
	##Find repeated names. Should only have one match within each STITLS cross match
	repeat_names = [name for name in set(match_names) if match_names.count(name)>1]
	
	##Find the list index number for each repeated source
	repeat_lists = []
	for name in repeat_names:
		repeat_list = []
		for i in xrange(len(match_names)):
			if match_names[i] == name: repeat_list.append(i)
		repeat_lists.append(repeat_list)
	
	##For each repeated source, retain the smallest separation, and note the other
	##row numbers to skip when reading in the main data
	for ls in repeat_lists:
		seps = [separations[ind] for ind in ls]
		keep_ind = ls[seps.index(min(seps))]
		for ind in ls:
			if ind != keep_ind: skip_row.append(ind)
	skip_rows.append(skip_row)
	
primary_names = []
source_matches = []
scaled_source_nums = []

##Read in the data. First, check each match for a new primary catalogue entry. If new, create a new
##source_info class and append primary cat information
for cat,skip in zip(matched_cats,skip_rows):
	match_name = 'matched_%s_%s.vot' %(primary_cat,cat)
	data,rows = open_table(match_name)
	
	prim_dens = float(data.keywords["%s_nu" %primary_cat])
	match_dens = float(data.keywords["%s_nu" %cat])
	if len(scaled_source_nums)==0:
		scaled_source_nums.append(prim_dens)
		scaled_source_nums.append(match_dens)
	else:
		scaled_source_nums.append(match_dens)
	
	for row in data:
		primary_name = row[0]
		if primary_name not in primary_names: 
			primary_names.append(primary_name)
			src = source_info()
			src.cats.append(primary_cat)
			src.names.append(str(row[0]))
			src.ras.append(str(row[1]))
			src.rerrs.append(str(row[2]))
			src.decs.append(str(row[3]))
			src.derrs.append(str(row[4]))
			src.fluxs.append(str(row[5]))
			src.ferrs.append(str(row[6]))
			src.majors.append(str(row[7]))
			src.minors.append(str(row[8]))
			src.PAs.append(str(row[9]))
			##Sometimes the flag data is empty - need to change to a null
			##result, otherwise goes as a space and later indexing gets ruined
			if str(row[10])=='':
				src.flags.append('-100000.0')
			else:
				src.flags.append(str(row[10]))
			src.IDs.append(str(row[11]))
			src.freqs.append(primary_freq)
			source_matches.append(src)
	
	##Second, get all the matched source data, and append to the appropriate
	##primary catalogue source
	for s_row in xrange(rows):
		##If in the skip list, it's a repeated source, so don't add it to the matched data
		if s_row in skip:
			pass
		else:
			row = data.row(s_row)
			primary_name = row[0]
			ind = primary_names.index(primary_name)
			src = source_matches[ind]
			src.cats.append(cat)
			src.names.append(str(row[12]))
			src.ras.append(str(row[13]))
			src.rerrs.append(str(row[14]))
			src.decs.append(str(row[15]))
			src.derrs.append(str(row[16]))
			src.majors.append(str(row[19]))
			src.minors.append(str(row[20]))
			src.PAs.append(str(row[21]))
			##Sometimes the flag data is empty - need to change to a null
			##result, otherwise goes as a space and later indexing gets ruined
			if str(row[22])=='':
				src.flags.append('-100000.0')
			else:
				src.flags.append(str(row[22]))
			src.IDs.append(str(row[23]))
			ind_freq = matched_cats.index(cat)
			cat_freq = matched_freqs[ind_freq]
			cat_freqs = cat_freq.split('~')
			
			##If the catalogue has more than one frequency, append all entries to a list
			##If just one, create a one entry list
			if len(cat_freqs)>1:
				freqss = []
				fluxss = []
				ferrss = []
				fluxss.append(str(row[17]))
				ferrss.append(str(row[18]))
				freqss.append(str(cat_freqs[0]))
				for i in xrange(len(cat_freqs)-1):
					fluxss.append(str(row[24+(i*1)]))
					ferrss.append(str(row[25+(i*1)]))
					freqss.append(str(cat_freqs[i+1]))
				src.fluxs.append(fluxss)
				src.freqs.append(freqss)
				src.ferrs.append(ferrss)
			else:
				src.fluxs.append(str(row[17]))
				src.ferrs.append(str(row[18]))
				src.freqs.append(str(cat_freq))
			
##Add errors in quadrature
def error(o1,o2):
	o1 = o1*dr
	o2 = o2*dr
	return o1**2+o2**2

##Create possible premutations of sources in a group
def do_match(matches,comp):
	new_match = []
	for obj in matches:
		for source in comp:
			match = [x for x in obj]
			match.append(source)
			new_match.append(match)
	return new_match

###Calculate area of sky bound by 2 points on celestial sphere
#def get_lune(ra1,ra2,dec1,dec2):
	#return abs((ra2*dr-ra1*dr)*(np.sin(dec2*dr)-np.sin(dec1*dr)))
	
##Calculates angular distance between two celestial points - input degrees, output radians
def calc_dist(RA1,RA2,Dec1,Dec2):   
	in1 = (90.0 - Dec1)*dr
	in2 = (90.0 - Dec2)*dr
	RA_d = (RA1 - RA2)*dr
	cosalpha = np.cos(in1)*np.cos(in2) + np.sin(in1)*np.sin(in2)*np.cos(RA_d)
	alpha = np.arccos(cosalpha)
	return alpha

#@profile
def do_bayesian(sources):
	##Calculate the weight of each source based on astrometric precision
	#for source in sources: print source.rerr,source.derr 
	weights = [1/error(source.rerr,source.derr) for source in sources]
	
	##Calculate the bayes factor
	n = len(weights)
	sum1 = sum(weights)
	prod1 = 1
	for i in weights: prod1*=i
	sum2=0
	for i in range(0,n+1):
		for j in range(i+1,n):
			##Calculate the distance between the two sources
			distance = calc_dist(sources[i].ra,sources[j].ra,sources[i].dec,sources[j].dec)
			##Calcualte the second sum
			sum2+=(weights[i]*weights[j]*distance**2)
			
	bayes_factor = (2**(n-1))*(prod1/sum1)*np.exp(-(sum2/(2*sum1)))
	
	##Work out which catalogues are present in this combination
	comps_names = [src.cat for src in sources]
	inds = [all_cats.index(name) for name in comps_names]
	
	##Find out the scaled number of sources for the catalogues present
	source_nums = [scaled_source_nums[i] for i in inds]
	
	#Calcualte prior
	prod_nums = 1
	for i in source_nums: prod_nums *= i
	prior = (source_nums[0]/prod_nums)
	
	##Calculate posterior
	#posterior = (bayes_factor*prior)/(1+(bayes_factor*prior))   ##The approximation to the posterior in the lims of small priors
	posterior = (1 + ((1 - prior)/(bayes_factor*prior)))**-1
	
	return prior, bayes_factor, posterior
	
out_file = open(out_name,'w+')

##Make a list of all catalogues (we needed the primary catalogue name in a separate list
##when reading in data, now we need all names in one list in the order the appear in the
##matched vot tables
all_cats = matched_cats
all_cats.insert(0,primary_cat)

##MAIN LOOP OF THE COOOOOOOODE!!!
for src in source_matches:
	##Write all information for each source in an indivudual line for each group
	##Account for catalogues with more than one frequency
	out_file.write('START_GROUP\n')
	for i in xrange(len(src.names)):
		if type(src.freqs[i])!=list:
			source_string = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" %(src.cats[i],src.names[i],src.ras[i],src.rerrs[i],
			src.decs[i],src.derrs[i],src.freqs[i],src.fluxs[i],src.ferrs[i],src.majors[i],src.minors[i],src.PAs[i],src.flags[i],src.IDs[i])
			out_file.write(source_string)
		else:
			flux_str=''
			for j in xrange(len(src.freqs[i])): flux_str+= "%s %s %s " %(src.freqs[i][j],src.fluxs[i][j],src.ferrs[i][j])
			source_string = "%s %s %s %s %s %s %s%s %s %s %s %s\n" %(src.cats[i],src.names[i],src.ras[i],src.rerrs[i],
			src.decs[i],src.derrs[i],flux_str,src.majors[i],src.minors[i],src.PAs[i],src.flags[i],src.IDs[i])
			out_file.write(source_string)
	
	##Separate the grouped information in to source_single classes and append to cats
	##in a specific order
	cats = [[] for i in all_cats]
	
	for i in xrange(len(src.names)):
		sing_src = source_single()
		sing_src.cat = src.cats[i]
		sing_src.name = src.names[i]
		sing_src.ra = float(src.ras[i])
		sing_src.rerr = float(src.rerrs[i])
		sing_src.dec = float(src.decs[i])
		sing_src.derr = float(src.derrs[i])
		sing_src.fluxs = src.fluxs[i]
		sing_src.ferrs = src.ferrs[i]
		sing_src.major = src.majors[i]
		sing_src.minor = src.minors[i]
		sing_src.PA = src.PAs[i]
		sing_src.flag = src.flags[i]
		sing_src.ID = src.IDs[i]
		sing_src.freqs = src.freqs[i]
		
		ind_name = all_cats.index(src.cats[i])
		cats[ind_name].append(sing_src)
	
	##In the next 5 lines, all possible combinations of the catalogues are created
	matches = [cats[0]]
	del cats[0]
	comps = [cat for cat in cats if len(cat)>0]
	for i in range(0,len(comps)):
		matches = do_match(matches,comps[i])
		
	out_file.write('START_COMP\n')

	##For each combination of sources, calculate a bayesian factor, a prior and a posterior probability
	##Write all source information for that match along with prior, bayes_factor and posterior to a single line
	##If a catalogue hasn't been matched, create a string of -100000.0. This is neccessary for create_table.py
	##and plot_image.py to be able to automatically pick out the correct information.
	for option in matches:
		catss = [source.cat for source in option]
		
		prior,bayes,posterior =  do_bayesian(option)
		match_str=''
		for cat in all_cats:
			if cat in catss:
				src = option[catss.index(cat)]
				if type(src.freqs)!=list:
					
					match_str += "%s %s %s %s %s %s %s %s %s %s %s %s %s %s " %(src.cat,src.name,str(src.ra),str(src.rerr),
					str(src.dec),str(src.derr),src.freqs,src.fluxs,src.ferrs,src.major,src.minor,src.PA,src.flag,src.ID)
				else:
					flux_str=''
					for j in xrange(len(src.freqs)): flux_str+= "%s %s %s " %(src.freqs[j],src.fluxs[j],src.ferrs[j])
					match_str += "%s %s %s %s %s %s %s%s %s %s %s %s " %(src.cat,src.name,str(src.ra),str(src.rerr),
					str(src.dec),str(src.derr),flux_str,src.major,src.minor,src.PA,src.flag,src.ID)
			else:
				if num_freqs[all_cats.index(cat)]==1:
					match_str += "-100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 "
				else:
					extra_str = ""
					for j in xrange(num_freqs[all_cats.index(cat)]-1): extra_str+= "-100000.0 -100000.0 -100000.0 "
					match_str += "-100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 " \
					"-100000.0 -100000.0 -100000.0 -100000.0 %s" %extra_str
		match_str += '%s %s %s\n' %(str(bayes),str(prior),str(posterior))
		out_file.write(match_str)
	out_file.write('END_COMP\n')
	out_file.write('END_GROUP\n')
	
out_file.close()

#def main():
    #do_this_shit()
 
#if __name__ == "__main__":
    #main()
