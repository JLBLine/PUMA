import atpy
import numpy as np
#import sys
#sys.path.append("/home/jline/Data/tools")
#import tool_lib as tool
import subprocess
#import matplotlib.pyplot as plt
#from matplotlib.patches import Ellipse
import optparse


parser = optparse.OptionParser()

parser.add_option('-a', '--table1',
	help='Enter name of table1')
	
parser.add_option('-b', '--table2',
	help='Enter name of table2')
	
parser.add_option('-c', '--details1', 
	help='Enter table1 details name,RA,RA_error,Dec,Dec_error,freq,flux,flux_error,PA,major,minor,flags,ID,freq2,flux2,flux2_error,flux3,flux3_error etc.' \
	'Any columns that do not exist, input -. Add as many fluxs as required')
	
parser.add_option('-d', '--details2',
	help='Enter table2 details')
	
parser.add_option('-e', '--units1',
	help='Enter table1 units in order of RA,RA_error,Dec,Dec_error,flux,flux_error,PA,major,minor' \
	     'Enter as| for angles: either deg,arcsec,arcmin,min,sec  | for flux: either Jy,mJy      ')
	
parser.add_option('-f', '--units2',
	help='Enter table2 units')
	
parser.add_option('-g', '--prefix1',
	help='Enter name of table1')
	
parser.add_option('-i', '--prefix2',
	help='Enter name of table2')
	
parser.add_option('-m', '--make_table',action='store_true',default=False,
	help='Add to not create tables')

parser.add_option('-s', '--separation',
	help='Enter matching radius in arcsecs')
	
parser.add_option('-w', '--ra_lims1', 
	help='Enter RA lims to calculate source density of catalogue 1')
	
parser.add_option('-x', '--dec_lims1', 
	help='Enter Dec lims to calculate source density of catalogue 1')
	
parser.add_option('-y', '--ra_lims2', 
	help='Enter RA lims to calculate source density of catalogue 2')
	
parser.add_option('-z', '--dec_lims2', 
	help='Enter Dec lims to calculate source density of catalogue 2')

options, args = parser.parse_args()


dr = np.pi/180.0

def get_units(data,detail,unit,unit_type,entries):
	if detail=='-':
		column = np.empty(entries); column.fill(-100000.0)
	else:
		if unit_type=='angle':
			if unit=='deg':
				column = data[detail]
			elif unit=='arcmin':
				column = data[detail]/60.0
			elif unit=='arcsec':
				column = data[detail]/3600.0
			elif unit=='min':
				column = data[detail]*(15.0/60.0)
			elif unit=='sec':
				column = data[detail]*(15.0/3600.0)
			else:
				print 'angle coversion error'
		if unit_type=='flux':
			if unit=='Jy':
				column = data[detail]
			elif unit=='mJy':
				column = data[detail]/1000.0
			else:
				print 'flux conversion error'
	return column


def get_units_blanks(data,detail,unit,unit_type,entries):
	column = np.empty(entries); column.fill(-100000.0)
	if detail=='-':
		return column
	else:
		for i in xrange(entries):
			try:
				entry = float(data[detail][i])
				if unit_type=='angle':
					if unit=='deg':
						column[i] = entry
					elif unit=='arcmin':
						column[i] = entry/60.0
					elif unit=='arcsec':
						column[i] = entry/3600.0
					elif unit=='min':
						column[i] = entry*(15.0/60.0)
					elif unit=='sec':
						column[i] = entry*(15.0/3600.0)
					else:
						print 'angle coversion error'
				if unit_type=='flux':
					if unit=='Jy':
						column[i] = entry
					elif unit=='mJy':
						column[i] = entry/1000.0
					else:
						print 'flux conversion error'
			except ValueError:
				column[i]=-100000.0
				
		return column
			
def make_table(data,details,units,prefix,ra_lims,dec_lims):
	shape = data.shape
	entries = shape[0]
	names = data[details[0]]
	
	RA = get_units(data,details[1],units[0],'angle',entries)
	RA_error = get_units(data,details[2],units[1],'angle',entries)
	Dec = get_units(data,details[3],units[2],'angle',entries)
	Dec_error = get_units(data,details[4],units[3],'angle',entries)
	freq = details[5]
	flux = get_units_blanks(data,details[6],units[4],'flux',entries)
	flux_error = get_units_blanks(data,details[7],units[5],'flux',entries)
	PA = get_units_blanks(data,details[8],units[6],'angle',entries)
	major = get_units_blanks(data,details[9],units[7],'angle',entries)
	minor = get_units_blanks(data,details[10],units[8],'angle',entries)
	if details[11]=='-':
		flags = np.empty(entries); flags.fill(-100000.0)
	else:
		flags = data[details[11]]
	if details[12]=='-':
		ID = np.empty(entries); ID.fill(-100000.0)
	else:
		ID = data[details[12]]

	freqs = []
	fluxs = []
	ferrs = []
	if len(details)>13:
		length = len(details) - 13
		num_of_freqs = length / 3
		for i in xrange(num_of_freqs):
			freqs.append(details[13+(3*i)])
			fluxs.append(get_units_blanks(data,details[13+1+(3*i)],'Jy','flux',entries))
			ferrs.append(get_units_blanks(data,details[13+2+(3*i)],'Jy','flux',entries))
	t=atpy.Table(masked=True)
	
	t.add_column('%s_name' %prefix,names,description='Name from catalogue')
	t.add_column('%s_RAJ2000' %prefix,RA,description='Right Ascension of source (J2000)',unit='deg')
	t.add_column('%s_e_RAJ2000' %prefix,RA_error,description='Error on Right Ascension',unit='deg')
	t.add_column('%s_DEJ2000' %prefix,Dec,description='Declination of source (J2000)',unit='deg')
	t.add_column('%s_e_DEJ2000' %prefix,Dec_error,description='Error on Declination',unit='deg')
	t.add_column('%s_S%d' %(prefix,float(freq)),flux,description='Source flux at %.1fMHz' %float(freq),unit='Jy')
	t.add_column('%s_e_S%d' %(prefix,float(freq)),flux_error,description='Error on flux at %.1fMHz' %float(freq),unit='Jy')
	t.add_column('%s_MajAxis' %prefix,major,description='Fitted major axis',unit='deg')
	t.add_column('%s_MinAxis' %prefix,minor,description='Fitted minor axis',unit='deg')
	t.add_column('%s_PA' %prefix,PA,description='Fitted Position Angle',unit='deg')
	t.add_column('%s_flag' %prefix,flags,description='Any meta flag for inclusion')
	t.add_column('%s_FieldID' %prefix,ID,description='If avaiable, image field ID')
	if len(freqs)>0:
		for i in xrange(len(freqs)):
			t.add_column('%s_Flux_%.1f' %(prefix,float(freqs[i])),fluxs[i],description='Source flux at %.1fMHz' %float(freqs[i]),unit='Jy')
			t.add_column('%s_Flux_%.1f_err' %(prefix,float(freqs[i])),ferrs[i],description='Error on flux at %.1fMHz' %float(freqs[i]),unit='Jy')
			
	def get_lune(ra1,ra2,dec1,dec2):
		return abs((ra2*dr-ra1*dr)*(np.sin(dec2*dr)-np.sin(dec1*dr)))
			
	##If the ra lims do not cross over the RA = 0 line
	if ra_lims[0]<ra_lims[1]:
		sources_in_bounds =  [i for i in xrange(len(RA)) if (RA[i]>=ra_lims[0] and RA[i]<=ra_lims[1]) and (Dec[i]>=dec_lims[0] and Dec[i]<=dec_lims[1])]
		area = get_lune(ra_lims[0],ra_lims[1],dec_lims[0],dec_lims[1])
	##If they do, searching the coords slightly differently and do some rearranging to get the area
	else:
		sources_in_bounds =  [i for i in xrange(len(RA)) if (RA[i]>=ra_lims[0] or RA[i]<=ra_lims[1]) and (Dec[i]>=dec_lims[0] and Dec[i]<=dec_lims[1])]
		extra = 360.0 - ra_lims[0]
		area = get_lune(0,ra_lims[1]+extra,dec_lims[0],dec_lims[1])
		
	scaled_source_density = (4*np.pi*len(sources_in_bounds))/area
			
	t.add_keyword('%s_nu' %prefix,str(scaled_source_density))
	t.write('simple_%s.vot' %(prefix),votype='ascii',overwrite=True)
	
	return scaled_source_density
	
prefix1 = options.prefix1
prefix2 = options.prefix2

if options.make_table==True:
	details1 = options.details1.split(',')
	details2 = options.details2.split(',')
	units1 = options.units1.split(',')
	units2 = options.units2.split(',')
	ra_lims1 = map(int,options.ra_lims1.split(','))
	ra_lims2 = map(int,options.ra_lims2.split(','))
	dec_lims1 = map(int,options.dec_lims1.split(','))
	dec_lims2 = map(int,options.dec_lims2.split(','))
	table1 = options.table1
	table2 = options.table2
	data1 = atpy.Table(table1,verbose=False)
	data2 = atpy.Table(table2,verbose=False)
	
	ss_dens1 = make_table(data1,details1,units1,prefix1,ra_lims1,dec_lims1)
	ss_dens2 = make_table(data2,details2,units2,prefix2,ra_lims2,dec_lims2)
	
else:
	##Read in the tables just to get one lousy number?!?!?!
	simp_1 = atpy.Table("simple_%s.vot" %prefix1)
	simp_2 = atpy.Table("simple_%s.vot" %prefix2)

	ss_dens1 = simp_1.keywords["%s_nu" %prefix1]
	ss_dens2 = simp_2.keywords["%s_nu" %prefix2]
	
##Need to have the directory your command stilts in is added to your PATH, ie on my machine I've added the line
## export PATH=/usr/local/STILTS:$PATH  to my ~/.bashrc
match_string = 'stilts tmatch2 join=1and2 find=all matcher=sky params=%d ' %float(options.separation)
match_string += "in1=simple_%s.vot values1='%s_RAJ2000 %s_DEJ2000' suffix1='' " %(prefix1,prefix1,prefix1)
match_string += "in2=simple_%s.vot values2='%s_RAJ2000 %s_DEJ2000' suffix2='' " %(prefix2,prefix2,prefix2)
match_string += "fixcols=all out=matched_%s_%s.vot" %(prefix1,prefix2)

subprocess.call(match_string,shell=True)

match_table = atpy.Table("matched_%s_%s.vot" %(prefix1,prefix2),verbose=False)

match_table.add_keyword("%s_nu" %prefix1,str(ss_dens1))
match_table.add_keyword("%s_nu" %prefix2,str(ss_dens2))