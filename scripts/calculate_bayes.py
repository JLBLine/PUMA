#!/usr/bin/env python
import numpy as np
import optparse
import make_table_lib as mkl
try:
	import pyfits as fits
except ImportError:
	from astropy.io import fits

import warnings
from astropy.utils.exceptions import AstropyUserWarning
from astropy.io.votable.exceptions import VOTableSpecWarning
warnings.simplefilter('ignore', category=AstropyUserWarning)
warnings.simplefilter('ignore', category=VOTableSpecWarning)
from time import ctime
#from cross_match import get_lune


##Get all of the input variables
parser = optparse.OptionParser()

parser.add_option('-p', '--primary_cat',
	help='Enter primary catalogue prefix')

parser.add_option('-g', '--primary_freq',
	help='Enter frequency of primary catalogue. If more than one frequency, separate by ~')

parser.add_option('-m', '--matched_cats',
	help='Enter names of matched cataloges, separated by commas')

parser.add_option('-f', '--matched_freqs',
	help='Enter frequencies of matched catalogues, separated by commas. If more than one frequency, separate by ~')

parser.add_option('-o', '--out_name',
	help='Enter name of output matched group file')

parser.add_option('-r', '--resolution',
	help='Resolution of base catalogue in "deg:arcmin:arcsec" ie "00:03:00" for 3arcmins')

parser.add_option('-a', '--prob_thresh', default=0.5,
	help='Probability threshold; if all combination of sources are under this value, check for any source from a single matched catalogue being outside resolution + error of base source. Default of 0.5')

options, args = parser.parse_args()

##Make all of the information usable
primary_cat = options.primary_cat
primary_freqs = options.primary_freq.split('~')
matched_cats = options.matched_cats.split(',')
matched_freqs = options.matched_freqs.split(',')

num_freqs = [len(primary_cat.split('~'))]
for freq in matched_freqs: num_freqs.append(len(freq.split('~')))
num_freqs = list(map(int,num_freqs))

closeness = mkl.dec_to_deg(options.resolution)/2

prob_thresh = float(options.prob_thresh)

out_name = options.out_name

dr = np.pi/180.0

def get_lune(ra1,ra2,dec1,dec2):
	'''Calculates the steradian coverage of a lune defined by two RA,Dec
	coords'''
	return abs((ra2*dr-ra1*dr)*(np.sin(dec2*dr)-np.sin(dec1*dr)))

##Opens each indvidual fits match table
def open_table(name):
	table = fits.open(name)
	data = table[1].data
	shape = data.shape
	rows = shape[0]
	##Find the number of sources scaled for to the whole sky
	nu_1 = table[1].header['nu_1']
	nu_2 = table[1].header['nu_2']
	##Find the boundaries of the catalogue coverage
	bound_1 = table[1].header['bound_lims_1']
	bound_2 = table[1].header['bound_lims_2']
	##Get names of columns - this line works for both pyfits and astropy
	colnames = table[1].columns
	colnames = colnames.info(attrib='name',output=False)['name']
	table.close()

	return data,rows,nu_1,nu_2,colnames,bound_1,bound_2

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
		self.close = None

##The way STILTS is run in cross_match.py means sources in the matched catalogues
##may be matched to more than one base catalogue source. Read through all the matches
##and check for any repeated matched sources; if they exist, choose the match with the
##smallest angular separation. Delete the other matched row

##Will contain lists of rows to be skipped that correspond to each catalogue in matched_cats
skip_rows = []
table_data = []

##TODO: Speed this up (limiting loop at the moment)
for cat in matched_cats:
	match_names = []
	separations = []
	skip_row = []
	##Get data
	match_name = 'matched_%s_%s.fits' %(primary_cat,cat)
	data,rows,nu_1,nu_2,colnames,bound_1,bound_2 = open_table(match_name)

	##Find which row entry is the matched catalogue name
	colstart = colnames.index('%s_name' %cat)

	#table_data.append
	for row in data:
		match_names.append(str(row[colstart]))
		separations.append(float(row[-1]))

	##Find repeated names. Should only have one match within each STITLS cross match
	repeat_names = [name for name in set(match_names) if match_names.count(name)>1]

	##Find the list index number for each repeated source
	repeat_lists = []
	for name in repeat_names:
		repeat_list = []
		for i in np.arange(len(match_names)):
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
	table_data.append([data,rows,nu_1,nu_2,colnames,cat,colstart,bound_1,bound_2])

scaled_source_nums = []
cat_bounds = []
source_matches = {}

##Read in the data. First, check each match for a new primary catalogue entry. If new, create a new
##source_info class and append primary cat information
for cat_data,skip in zip(table_data,skip_rows):

	###Find the source densities from the tables that were calcualted by cross_match.py,
	###and add to scaled_source_nums list
	data,rows,prim_dens,match_dens,colnames,cat,colstart,prim_bound,match_bound = cat_data

	if len(scaled_source_nums)==0:
		scaled_source_nums.append(prim_dens)
		scaled_source_nums.append(match_dens)
		cat_bounds.append(prim_bound)
		cat_bounds.append(match_bound)
	else:
		scaled_source_nums.append(match_dens)
		cat_bounds.append(match_bound)
	##If any of the errors are negative or zero, the whole match stops and fails -
	##take an average of rerr, derr and ferr here, and subsitute that in later if
	##needs be
	base_rerr_avg = np.mean(data['%s_e_RAJ2000' %primary_cat])
	base_derr_avg = np.mean(data['%s_e_DEJ2000' %primary_cat])
	base_ferr_avg = np.mean(data['%s_e_S%d' %(primary_cat,float(primary_freqs[0]))])

	cat_rerr_avg = np.mean(data['%s_e_RAJ2000' %cat])
	cat_derr_avg = np.mean(data['%s_e_DEJ2000' %cat])
	cat_ferr_avg = np.mean(data['%s_e_S%d' %(cat,float(matched_freqs[matched_cats.index(cat)].split('~')[0]))])
	##(For all rows of data)
	for s_row in np.arange(rows):
		##If in the skip list, it's a repeated source, so don't add it to the matched data
		if s_row in skip:
			pass
		else:
			row = data[s_row]
			primary_name = row[0]
			if primary_name not in source_matches:
				src = source_info()
				src.cats.append(primary_cat)
				src.names.append(str(row[0]))
				src.ras.append(str(row[1]))
				src.decs.append(str(row[3]))

				##Test the errors for zero or neg values,
				##and insert a proxy error if so
				if float(row[2])<=0.0:
					src.rerrs.append(str(base_rerr_avg))
				else:
					src.rerrs.append(str(row[2]))

				if float(row[4])<=0.0:
					src.derrs.append(str(base_derr_avg))
				else:
					src.derrs.append(str(row[4]))

				src.majors.append(str(row[7]))
				src.minors.append(str(row[8]))
				src.PAs.append(str(row[9]))
				##Sometimes the flag data is empty - need to change to a null
				##result, otherwise goes as a space and later indexing gets ruined
				if str(row[10])=='':
					src.flags.append('-100000.0')
				elif str(row[10])=='--':
					src.flags.append('-100000.0')
				else:
					src.flags.append(str(row[10]))
				if str(row[11])=='':
					src.IDs.append('-100000.0')
				elif str(row[11])=='--':
					src.IDs.append('-100000.0')
				else:
					src.IDs.append(str(row[11]))

				if len(primary_freqs)>1:
					freqss = []
					fluxss = []
					ferrss = []
					freqss.append(str(primary_freqs[0]))
					fluxss.append(str(row[5]))
					if -100000.0<float(row[6])<=0.0:
						ferrss.append(str(cat_ferr_avg))
					else:
						ferrss.append(str(row[6]))
					for i in np.arange(len(primary_freqs)-1):
						freqss.append(str(primary_freqs[i+1]))
						fluxss.append(str(row[12+(2*i)]))
						if -100000.0<float(row[13+(2*i)])<=0.0:
							ferrss.append(str(base_ferr_avg))
						else:
							ferrss.append(str(row[13+(2*i)]))

					src.fluxs.append(fluxss)
					src.freqs.append(freqss)
					src.ferrs.append(ferrss)
				else:
					src.freqs.append(primary_freqs[0])
					src.fluxs.append(str(row[5]))
					if -100000.0<float(row[6])<=0.0:
						src.ferrs.append(str(base_ferr_avg))
					else:
						src.ferrs.append(str(row[6]))
				#source_matches.append(src)
				source_matches[str(row[0])] = src

	##Second, get all the matched source data, and append to the appropriate
	##primary catalogue source - we worked out where the information would start
	##earlier on, at colstart
	for s_row in np.arange(rows):
		##If in the skip list, it's a repeated source, so don't add it to the matched data
		if s_row in skip:
			pass
		else:
			row = data[s_row]
			src = source_matches[row[0]]

			src.cats.append(cat)
			src.names.append(str(row[colstart]))       ##-12
			src.ras.append(str(row[colstart+1]))
			src.decs.append(str(row[colstart+3]))
			##Test the errors for zero or neg values,
			##and insert a proxy error if so
			if float(row[colstart+2])<=0.0:
					src.rerrs.append(str(cat_rerr_avg))
			else:
				src.rerrs.append(str(row[colstart+2]))

			if float(row[colstart+4])<=0.0:
				src.derrs.append(str(cat_derr_avg))
			else:
				src.derrs.append(str(row[colstart+4]))
			src.majors.append(str(row[colstart+7]))
			src.minors.append(str(row[colstart+8]))
			src.PAs.append(str(row[colstart+9]))
			##Sometimes the flag data is empty - need to change to a null
			##result, otherwise goes as a space and later indexing gets ruined
			if str(row[colstart+10])=='':
				src.flags.append('-100000.0')
			elif str(row[colstart+10])=='--':
				src.flags.append('-100000.0')
			else:
				src.flags.append(str(row[colstart+10]))
			if str(row[colstart+11])=='':
				src.IDs.append('-100000.0')
			elif str(row[colstart+11])=='--':
				src.IDs.append('-100000.0')
			else:
				src.IDs.append(str(row[colstart+11]))
			ind_freq = matched_cats.index(cat)
			cat_freq = matched_freqs[ind_freq]
			cat_freqs = cat_freq.split('~')
			##If the catalogue has more than one frequency, append all entries to a list
			##If just one, create a one entry list
			if len(cat_freqs)>1:
				freqss = []
				fluxss = []
				ferrss = []
				freqss.append(str(cat_freqs[0]))
				fluxss.append(str(row[colstart+5]))
				if -100000.0<float(row[colstart+6])<=0.0:
					ferrss.append(str(cat_ferr_avg))
				else:
					ferrss.append(str(row[colstart+6]))
				for i in np.arange(len(cat_freqs)-1):
					freqss.append(str(cat_freqs[i+1]))
					fluxss.append(str(row[colstart+12+(2*i)]))
					if -100000.0<float(row[colstart+13+(2*i)])<=0.0:
						ferrss.append(str(cat_ferr_avg))
					else:
						ferrss.append(str(row[colstart+13+(2*i)]))

				src.fluxs.append(fluxss)
				src.freqs.append(freqss)
				src.ferrs.append(ferrss)
			else:
				src.fluxs.append(str(row[colstart+5]))
				if -100000.0<float(row[colstart+6])<=0.0:
					src.ferrs.append(str(cat_ferr_avg))
				else:
					src.ferrs.append(str(row[colstart+6]))
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

##Calculates angular distance between two celestial points - input degrees, output radians
def calc_dist(RA1,RA2,Dec1,Dec2):
	in1 = (90.0 - Dec1)*dr
	in2 = (90.0 - Dec2)*dr
	RA_d = (RA1 - RA2)*dr
	cosalpha = np.cos(in1)*np.cos(in2) + np.sin(in1)*np.sin(in2)*np.cos(RA_d)
	##Sometimes get floating point errors if two sources at exactly
	##the same position, so account for this:
	if cosalpha>1.0: cosalpha = 1.0
	elif cosalpha<-1.0: cosalpha = -1.0
	alpha = np.arccos(cosalpha)
	return alpha

def get_sky_overlap(catalogue_bounds):
	'''Works out the area of the sky over which the catalogues overlap -
	able to split into multiple lunes of so required by present catalogue'''

	##The prior requires knowing the sky area over which the catalogues overlap
	##Start off setting the region as the base catalogue limits (first catalogue in list catalogue_bounds)
	##The following code has been set up to think of RA as going from -180 to 180 deg - if the
	##base catalogue crosses 180/-180, split overlap region into two - the code can then treat base the same

	##Intially start with sky region of base catalogue
	##Get the limits from cross_match.py
	overlap_ra_lower,overlap_ra_upper,overlap_dec_lower,overlap_dec_upper = list(map(float,catalogue_bounds[0].split(',')))
	#print catalogue_bounds
	#print overlap_ra_lower,overlap_ra_upper,overlap_dec_lower,overlap_dec_upper
	overlap_ra_lowers = []
	overlap_ra_uppers = []

	##From cross_match.py, if the base catalogue crosses the 0/360 deg RA border, it's lower ra
	##bound will be negative - can leave as one region
	if overlap_ra_lower < 0:
		overlap_ra_lowers.append(overlap_ra_lower)
		overlap_ra_uppers.append(overlap_ra_upper)
	##Or, if all RA sit between 0 and 180, can again leave alone
	elif overlap_ra_upper <= 180.0:
		overlap_ra_lowers.append(overlap_ra_lower)
		overlap_ra_uppers.append(overlap_ra_upper)
	##Or, if all sits btw 180-360, set to being negative:
	elif overlap_ra_lower > 180.0 and 180.0 <= overlap_ra_upper < 360.0:
		overlap_ra_lowers.append(float(overlap_ra_lower) - 360.0)
		overlap_ra_uppers.append(float(overlap_ra_upper) - 360.0)

	##If crosses the -180/180 bound, split into two
	else:
		current_ra_lower = float(overlap_ra_lower)
		current_ra_upper = float(overlap_ra_upper)
		overlap_ra_lower = -180
		overlap_ra_upper = current_ra_upper - 360.0
		overlap_ra_lowers.append(overlap_ra_lower)
		overlap_ra_uppers.append(overlap_ra_upper)

		overlap_ra_lower_2 = current_ra_lower
		overlap_ra_upper_2 = 180.0
		overlap_ra_lowers.append(overlap_ra_lower_2)
		overlap_ra_uppers.append(overlap_ra_upper_2)

	##For each matched catalogue
	for bound in catalogue_bounds[1:]:
		##Get limits as defined in cross_match.py
		ra_lower,ra_upper,dec_lower,dec_upper = list(map(float,bound.split(',')))
		##In this case, the matched catalogue crosses both the 0/360
		##and -180/180 lines - split into two for this case
		if ra_lower > ra_upper:
			for ind in np.arange(len(overlap_ra_lowers)):
				##using float here means we make a copy rather than assigning a name to element of list
				overlap_ra_lower,overlap_ra_upper = float(overlap_ra_lowers[ind]),float(overlap_ra_uppers[ind])

				##This overlap region has been killed, so ignore
				if overlap_ra_lower == None:
					pass
				else:

					##Matched cat doesn't cover bound lune at all
					if overlap_ra_lower > ra_upper and overlap_ra_upper < ra_lower:
						overlap_ra_lowers[ind],overlap_ra_uppers[ind] = None,None

					##In this case, there are two regions within the current lune that
					##the matched catalogue covers, with the central region of that lune
					##not covered by the matched catalogue - must create two new lunes
					elif overlap_ra_lower < ra_upper and overlap_ra_upper > ra_lower:
						##First replace current upper overlap with ra_lower cross to create
						##one of the replace lunes
						overlap_ra_uppers[ind] = ra_upper
						##Then, make a new lune bound by ra_lower, and the current overlap_ra_upper
						overlap_ra_lowers.append(ra_lower)
						overlap_ra_uppers.append(overlap_ra_upper)
					##Otherwise, there is only one part of the current lune covered by the matched catalogue;
					##just adjust the edge of the lunes if applicable
					else:
						if overlap_ra_lower < ra_upper < overlap_ra_upper:
							overlap_ra_uppers[ind] = ra_upper
						if overlap_ra_lower < ra_lower < overlap_ra_upper:
							overlap_ra_lowers[ind] = ra_lower

		else:
			##If the matched catalogue spans 0/360 RA border, things
			##are straight forward
			if ra_lower < 0.0 or ra_upper < 180.0:
				for ind in np.arange(len(overlap_ra_lowers)):
					overlap_ra_lower,overlap_ra_upper = overlap_ra_lowers[ind],overlap_ra_uppers[ind]
					##If region already killed, skip
					if overlap_ra_lower == None:
						pass
					else:
						##for these conditions, there is no overlap at all
						if ra_lower > overlap_ra_upper:
							overlap_ra_lowers[ind],overlap_ra_uppers[ind] = None,None
						elif ra_upper < overlap_ra_lower:
							overlap_ra_lowers[ind],overlap_ra_uppers[ind] = None,None
						else:
							if ra_lower > overlap_ra_lower:
								overlap_ra_lowers[ind] = ra_lower
							if ra_upper < overlap_ra_upper:
								overlap_ra_uppers[ind] = ra_upper
			##If the matched catalogue is fully between 180,360, set RAs to negatives
			elif ra_lower > 180.0:
				ra_lower -= 360.0
				ra_upper -= 360.0
				for ind in np.arange(len(overlap_ra_lowers)):
					overlap_ra_lower,overlap_ra_upper = overlap_ra_lowers[ind],overlap_ra_uppers[ind]
					##If region already killed, skip
					if overlap_ra_lower == None:
						pass
					else:
						##No overlap between region and matched catalogue
						if ra_lower > overlap_ra_upper:
							overlap_ra_lowers[ind],overlap_ra_uppers[ind] = None,None
						##No overlap between region and matched catalogue
						elif ra_upper < overlap_ra_lower:
							overlap_ra_lowers[ind],overlap_ra_uppers[ind] = None,None
						##Otherwise, adjust region boundaries appropriately
						else:
							if ra_lower > overlap_ra_lower:
								overlap_ra_lowers[ind] = ra_lower
							if ra_upper < overlap_ra_upper:
								overlap_ra_uppers[ind] = ra_upper

			##Here we may need to do the splitting of the matched catalogue into two ranges
			##Only here because matched catalogue crosses 180/-180 boundary
			else:
				ra_lower_cross = ra_upper - 360.0
				current_ra_upper = float(overlap_ra_upper)

				for ind in np.arange(len(overlap_ra_lowers)):
					##using float here means we make a copy rather than assigning a name to element of list
					overlap_ra_lower,overlap_ra_upper = float(overlap_ra_lowers[ind]),float(overlap_ra_uppers[ind])
					if overlap_ra_lower == None:
						pass
					else:
						##Matched cat doesn't cover bound lune at all
						if overlap_ra_lower > ra_lower_cross and overlap_ra_upper < ra_lower:
							overlap_ra_lowers[ind],overlap_ra_uppers[ind] = None,None

						##In this case, there are two regions within the current lune that
						##the matched catalogue covers, with the central region of that lune
						##not covered by the matched catalogue - must create two new lunes
						elif overlap_ra_lower < ra_lower_cross and overlap_ra_upper > ra_lower:
							##First replace current upper overlap with ra_lower cross to create
							##one of the replace lunes
							overlap_ra_uppers[ind] = ra_lower_cross
							##Then, make a new lune bound by ra_lower, and the current overlap_ra_upper
							overlap_ra_lowers.append(ra_lower)
							overlap_ra_uppers.append(overlap_ra_upper)
						##Otherwise, there is only one part of the current lune covered by the matched catalogue;
						##just adjust the edge of the lunes if applicable
						else:
							if overlap_ra_lower < ra_lower_cross < overlap_ra_upper:
								overlap_ra_uppers[ind] = ra_lower_cross
							if overlap_ra_lower < ra_lower < overlap_ra_upper:
								overlap_ra_lowers[ind] = ra_lower

		##Declination is way easier - doesn't loop!!
		if dec_lower > overlap_dec_lower: overlap_dec_lower = dec_lower
		if dec_upper < overlap_dec_upper: overlap_dec_upper = dec_upper

	#return overlap_ra_lowers, overlap_ra_uppers, overlap_dec_lower, overlap_dec_upper

	sky_area = 0

	for i in np.arange(len(overlap_ra_lowers)):
		if overlap_ra_lowers[i] == None:
			pass
		else:
			sky_area += get_lune(overlap_ra_lowers[i],overlap_ra_uppers[i],overlap_dec_lower,overlap_dec_upper)

	return sky_area

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
	##Find the sky boundaries of all catalogue involved in cross-match
	catalogue_bounds = [cat_bounds[i] for i in inds]

	#Calcualte prior
	prod_nums = 1
	for i in source_nums: prod_nums *= i
	prior = (source_nums[0]/prod_nums)

	sky_area = get_sky_overlap(catalogue_bounds)

	prior *= ((sky_area / 4*np.pi)**(1 - len(source_nums)))

	##Calculate posterior
	#posterior = (bayes_factor*prior)/(1+(bayes_factor*prior))   ##The approximation to the posterior in the lims of small priors -
	##uncomment if you want to get rid of divide by zero messages (possibly less accurate though)
	posterior = (1 + ((1 - prior)/(bayes_factor*prior)))**-1

	return prior, bayes_factor, posterior

##Make a list of all catalogues (we needed the primary catalogue name in a separate list
##when reading in data, now we need all names in one list in the order the appear in the
##matched vot tables
all_cats = matched_cats
all_cats.insert(0,primary_cat)

out_file = open(out_name,'w+')

##MAIN LOOP OF THE COOOOOOOODE!!!
#for src in source_matches:

loop = 0

# for meh,src in source_matches.iteritems():
match_keys = source_matches.keys()
for key in match_keys:
	src = source_matches[key]
	##Find the primary position and errors
	prim_ra,prim_dec,prim_rerr,prim_derr,prim_name = float(src.ras[0]),float(src.decs[0]),float(src.rerrs[0]),float(src.derrs[0]),src.names[0]
	##Separate the grouped information in to source_single classes and append to cats
	##in a specific order
	cats = [[] for i in all_cats]
	for i in np.arange(len(src.names)):
		sing_src = source_single()
		sing_src.cat = src.cats[i]
		sing_src.name = src.names[i]
		ra = float(src.ras[i])
		rerr = float(src.rerrs[i])
		dec = float(src.decs[i])
		derr = float(src.derrs[i])
		sing_src.ra = ra
		sing_src.rerr = rerr
		sing_src.dec = dec
		sing_src.derr = derr

		##Test here to see if source lies within ellipse of resolution plus error
		##Mark down as a flag or not (True = passed test, False = failed test)
		##First off, find the offsets. We will set up delta_RA to give us the equivalent offset in RA that corresponds to the
		##resolution, so we use the offset in RA, not the arcdistance
		ra_dist = ra - prim_ra
		##Code to cope if one source 359.9, other 0.1 etc.
		if abs(ra_dist) > 180.0:
			if ra > 180.0:
				ra -= 360.0
				ra_dist = prim_ra - ra
			else:
				prim_ra -= 360.0
				ra_dist = ra - prim_ra
		dec_dist = prim_dec - dec

		##This basically calculates the distance covered by the resolution in RA
		delta_RA = np.arccos((np.cos(closeness*dr)-np.sin(prim_dec*dr)**2)/np.cos(prim_dec*dr)**2)/dr
		ra_axis = prim_rerr + abs(delta_RA)
		dec_axis = prim_derr + closeness

		##Test to see if the source lies with an error ellipse created using semi-major
		##and minor axes defined by the ra and dec error of the base cat + half the resolution
		##of the base cat (closeness)
		ell_test = (ra_dist/ra_axis)**2 + (dec_dist/dec_axis)**2

		if ell_test <= 1:
			sing_src.close = True
		else:
			sing_src.close = False

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
	#del cats[0]
	comps = [cat for cat in cats[1:] if len(cat)>0]
	for i in range(0,len(comps)):
		matches = do_match(matches,comps[i])

	bayes_infos = []

	##Calculate the bayesian prob
	for option in matches:
		catss = [source.cat for source in option]
		prior,bayes,posterior =  do_bayesian(option)
		bayes_infos.append([prior,bayes,posterior])

	remove_cats = []

	##If all of the probabilities are lower than some threshold (default 0.5),
	##it could be possible that a single source from a single catalogue is messing
	##up all of the combinations - here, check if any single source from a catalogue
	##failed the error ellipse test - flag it for removal in remove_cats, and
	##remove it from cats
	if max([bayes[2] for bayes in bayes_infos]) < prob_thresh:
		##Check that there isn't only one cat matched, otherwise
		##we remove the only match and PUMA will accept the base cat
		##data alone
		num_matched = len([1 for cat in cats[1:] if len(cat)>1])
		##Don't remove a match if there is less 2 matched cats
		if num_matched < 2:
			pass
		else:
			##Don't try to remove base cat, so only [1:]
			for cat in cats[1:]:
				if len(cat) == 1:
					if cat[0].close == False:
						cat_ind = cats.index(cat)
						cats[cat_ind] = []
						remove_cats.append(cat[0].cat)

					##In the next 5 lines, all possible combinations of the catalogues are created
					##but without any flagged catalogues
					matches = [cats[0]]
					#del cats[0]
					comps = [cat for cat in cats[1:] if len(cat)>0]
					for i in range(0,len(comps)):
						matches = do_match(matches,comps[i])

					bayes_infos = []
					##Calculate the bayesian prob
					for option in matches:
						catss = [source.cat for source in option]
						prior,bayes,posterior =  do_bayesian(option)
						bayes_infos.append([prior,bayes,posterior])

	##Write all information for each source in an indivudual line for each group
	##Account for catalogues with more than one frequency
	out_file.write('START_GROUP\n')
	for i in np.arange(len(src.names)):
		if src.cats[i] in remove_cats:
			pass
		else:
			if type(src.freqs[i])!=list:
				#print src.cats[i],src.PAs[i],src.flags[i],src.IDs[i]
				source_string = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" %(src.cats[i],src.names[i],src.ras[i],src.rerrs[i],
				src.decs[i],src.derrs[i],src.freqs[i],src.fluxs[i],src.ferrs[i],src.majors[i],src.minors[i],src.PAs[i],src.flags[i],src.IDs[i])
				out_file.write(source_string)
			else:
				flux_str=''
				#print src.cats[i],src.PAs[i],src.flags[i],src.IDs[i]
				for j in np.arange(len(src.freqs[i])): flux_str+= "%s %s %s " %(src.freqs[i][j],src.fluxs[i][j],src.ferrs[i][j])
				source_string = "%s %s %s %s %s %s %s%s %s %s %s %s\n" %(src.cats[i],src.names[i],src.ras[i],src.rerrs[i],
				src.decs[i],src.derrs[i],flux_str,src.majors[i],src.minors[i],src.PAs[i],src.flags[i],src.IDs[i])
				out_file.write(source_string)

	out_file.write('START_COMP\n')
	##For each combination of sources, calculate a bayesian factor, a prior and a posterior probability
	##Write all source information for that match along with prior, bayes_factor and posterior to a single line
	##If a catalogue hasn't been matched, create a string of -100000.0. This is neccessary for create_table.py
	##and plot_image.py to be able to automatically pick out the correct information.
	for option in np.arange(len(matches)):
		catss = [source.cat for source in matches[option]]
		prior,bayes,posterior =  bayes_infos[option]
		match_str=''
		for cat in all_cats:
			if cat in catss and cat not in remove_cats:
				src = matches[option][catss.index(cat)]
				if type(src.freqs)!=list:
					match_str += "%s %s %s %s %s %s %s %s %s %s %s %s %s %s " %(src.cat,src.name,str(src.ra),str(src.rerr),
					str(src.dec),str(src.derr),src.freqs,src.fluxs,src.ferrs,src.major,src.minor,src.PA,src.flag,src.ID)
				else:
					flux_str=''
					for j in np.arange(len(src.freqs)): flux_str+= "%s %s %s " %(src.freqs[j],src.fluxs[j],src.ferrs[j])
					match_str += "%s %s %s %s %s %s %s%s %s %s %s %s " %(src.cat,src.name,str(src.ra),str(src.rerr),
					str(src.dec),str(src.derr),flux_str,src.major,src.minor,src.PA,src.flag,src.ID)
			else:
				if num_freqs[all_cats.index(cat)]==1:
					match_str += "-100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 "
				else:
					extra_str = ""
					for j in np.arange(num_freqs[all_cats.index(cat)]-1): extra_str+= "-100000.0 -100000.0 -100000.0 "
					match_str += "-100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 -100000.0 " \
					"-100000.0 -100000.0 -100000.0 -100000.0 %s" %extra_str
		match_str += '%s %s %s\n' %(str(bayes),str(prior),str(posterior))
		out_file.write(match_str)
	out_file.write('END_COMP\n')
	out_file.write('END_GROUP\n')

out_file.close()
