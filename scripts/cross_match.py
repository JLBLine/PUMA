#!/usr/bin/env python
from numpy import array,arange,sin,cos,pi,empty,where,ma,where,isnan
import subprocess
import optparse
import os
import sys
from astropy.table import Table, Column, MaskedColumn
from astropy.io.votable import parse as vot_parse
try:
    import pyfits as fits
except ImportError:
    from astropy.io import fits
from scipy.stats import mode
from distutils.spawn import find_executable

if find_executable('stilts') == None:
    sys.exit('Could not find `stilts`. It must be installed and on the PATH to use cross_match.py. Exiting now')


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
         'Enter as | for angles: either deg,arcsec,arcmin,min,sec | for flux: either Jy,mJy      ')

parser.add_option('-f', '--units2',
    help='Enter table2 units')

parser.add_option('-g', '--prefix1',
    help='Enter name of table1')

parser.add_option('-i', '--prefix2',
    help='Enter name of table2')

parser.add_option('-m', '--make_table',action='store_true',default=False,
    help='Add to store the simple tables')

parser.add_option('-s', '--separation',
    help='Enter matching radius in arcsecs')

parser.add_option('-v', '--verbose',action='store_true',default=False,
    help='Enter to turn on all of the astropy warnings')

parser.add_option('-w', '--ra_lims1',
    help='Enter RA lims to calculate source density of catalogue 1')

parser.add_option('-x', '--dec_lims1',
    help='Enter Dec lims to calculate source density of catalogue 1')

parser.add_option('-y', '--ra_lims2',
    help='Enter RA lims to calculate source density of catalogue 2')

parser.add_option('-z', '--dec_lims2',
    help='Enter Dec lims to calculate source density of catalogue 2')

options, args = parser.parse_args()

if options.verbose:
    pass
else:
    import warnings
    from astropy.utils.exceptions import AstropyUserWarning
    from astropy.io.votable.exceptions import VOTableSpecWarning
    warnings.simplefilter('ignore', category=AstropyUserWarning)
    warnings.simplefilter('ignore', category=VOTableSpecWarning)

dr = pi/180.0

def get_units(data,detail,unit,unit_type,entries):
    '''A function for either returning an array of -100000.0s of
    a given length, or taking a given array and performing unit conversion,
    given the specified units'''
    if detail=='-':
        column = empty(entries); column.fill(-100000.0)
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
                print('PUMA angle coversion error: you entered %s \nMust enter either deg, arcmin, arcsec, min or sec as units \nNow exiting' %unit)
        if unit_type=='flux':
            if unit=='Jy':
                column = data[detail]
            elif unit=='mJy':
                column = data[detail]/1000.0
            else:
                print('PUMA flux conversion error: you entered \nMust enter eith Jy or mJy as units \nNow exiting' %unit)

    column[where(isnan(column) == True)] = -100000.0
    return column

def get_units_blanks(data,detail,unit,unit_type,entries):
    '''A function for either returning an array of -100000.0s of
    a given length, or taking a given array and performing unit conversion,
    given the specified units. Can handle the input array having a blank
    input. Slower than get_units'''
    column = empty(entries); column.fill(-100000.0)
    if detail=='-':
        return column
    else:
        for i in arange(entries):
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
                        print('PUMA angle coversion error: you entered %s \nMust enter either deg, arcmin, arcsec, min or sec as units \nNow exiting' %unit)
                        sys.exit()
                if unit_type=='flux':
                    if unit=='Jy':
                        column[i] = entry
                    elif unit=='mJy':
                        column[i] = entry/1000.0
                    else:
                        print('PUMA flux conversion error: you entered \nMust enter eith Jy or mJy as units \nNow exiting' %unit)
                        sys.exit()
            except ValueError:
                column[i]=-100000.0
        column[where(isnan(column) == True)] = -100000.0
        return column

def get_lune(ra1,ra2,dec1,dec2):
        '''Calculates the steradian coverage of a lune defined by two RA,Dec
        coords'''
        return abs((ra2*dr-ra1*dr)*(sin(dec2*dr)-sin(dec1*dr)))

def shift_test(ra_range,dec_range):
    '''Finds the RA and Dec bounds of a catalogue. If the catalogue
    crosses the 0/360 but not the -180/180 bound, report the lower RA bound
    as negative. If it crosses both the 0/360 and -180/180, report the larger
    RA value as the lower bound so that we get ra2 -> ra1, instead of ra1 -> ra2'''

    ra_range = array(ra_range)
    min_ra,max_ra = min(ra_range),max(ra_range)
    min_dec,max_dec = min(dec_range),max(dec_range)

    min_shifts = []
    max_shifts = []

    ##Here we iteratively shift the zero reference point of 0-360 deg
    ##In this way, when we set the shift value to an empty patch of sky,
    ##we get the same minimum and maximum shifted RA values. We can
    ##then look for repeated min, max values via the mode of colleted
    ##min, maxes
    shift_array =  arange(0,359.0,0.1)
    for shift in shift_array:
        ##Find the indexes of all extries above shift value
        ##and then shift them by -360
        shift_inds = where(ra_range > shift)
        shift_values = ra_range[shift_inds] - 360

        ##Find all values not for shifting
        remain_inds = where(ra_range <= shift)
        remain_values = ra_range[remain_inds]

        if len(shift_values) > 0 and len(remain_values) == 0:
            min_shift,max_shift = shift_values.min(),shift_values.max()
        elif len(shift_values) == 0 and len(remain_values) > 0:
            min_shift,max_shift = remain_values.min(),remain_values.max()
        elif len(shift_values) > 0 and len(remain_values) > 0:
            min_shift,max_shift = min([shift_values.min(),remain_values.min()]),max([shift_values.max(),remain_values.max()])
        else:
            print('shoulnay happen')
        ##Collect the minimum and maximum for later
        min_shifts.append(min_shift)
        max_shifts.append(max_shift)


    ##Find the mode of the mins and maxes
    missing_inds_min = where(min_shifts == mode(min_shifts)[0][0])
    missing_inds = where(max_shifts == mode(max_shifts)[0][0])
    if len(missing_inds_min) != len(missing_inds):
        print('PUMA warning: mismatched length modes when deriving the bounds of the catalogues')

    ##Use them to find out the bounds of all missing RAs
    missing = shift_array[missing_inds]

    ##If the catalogue spans from 0 - 360.0, there won't be a mode - check using
    ##the length of the mode
    len_mode_min, len_mode_max = mode(min_shifts)[1][0],mode(max_shifts)[1][0]

    lower_ra = None
    upper_ra = None

    ##In this case the distribution covers all RA:
    if len_mode_min == 1 and len_mode_max == 1:
        lower_ra, upper_ra = -180,180.0
    else:
        ##In this case, the data crosses the 0/360 case
        if missing.min() > min_ra and missing.max() < max_ra:
            ##The data also crosses the -180/180 border - in this case, return
            ##the higher value ra as the lower bound - otherwise the bounds
            ##report the missing RAs rather than the RAs we care about
            if 180 >= missing.min() >= 0.0 and 180.0 >= missing.max() >= 0.0:
                lower_ra, upper_ra = missing.max(), missing.min()
            ##Same if both negative - report the larger as the lower bound
            elif missing.min() < 0.0 and missing.max() < 0.0:
                lower_ra, upper_ra = missing.max(), missing.min()
            ##Otherwise we need to make the upper missing RA as the lower bound,
            ##so we don't just pick out the missing RAs. If greater than 180,
            ##shift to make negative
            else:

                if missing.max() > 180.0:
                    lower_ra =  missing.max() - 360.0
                else:
                    lower_ra =  missing.max()

                if missing.min() > 180.0:
                    upper_ra = missing.min() - 360.0
                else:
                    upper_ra = missing.min()
        else:
            ##Other wise just a simple case where no crossing of either
            ##0/360 or -180/180
            lower_ra, upper_ra = min_ra,max_ra
            #print min_ra,max_ra,missing.min(),missing.max()

    return lower_ra,upper_ra,min_dec,max_dec

def make_table(data,details,units,prefix,ra_lims,dec_lims):
    ##Find all of the data in the columns as specified by the
    ##user inputs
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
        flags = empty(entries); flags.fill(-100000.0)
    else:
        flags = data[details[11]]
    if details[12]=='-':
        ID = empty(entries); ID.fill(-100000.0)
    else:
        ID = data[details[12]]
    freqs = []
    fluxs = []
    ferrs = []

    ##This handles the case of more than one frequency in a single catalogue
    ##Assumes the units of all of the fluxes, flux errors are the same
    if len(details)>13:
        length = len(details) - 13
        num_of_freqs = length // 3
        for i in arange(num_of_freqs):
            freqs.append(details[13+(3*i)])
            fluxs.append(get_units_blanks(data,details[13+1+(3*i)],units[4],'flux',entries))
            ferrs.append(get_units_blanks(data,details[13+2+(3*i)],units[5],'flux',entries))

    lower_ra,upper_ra,min_dec,max_dec = shift_test(RA,Dec)

    ##Count the number of sources within the requested user lune to calculate the scaled
    ##source density of the catalogues
    ##If the ra lims do not cross over the RA = 0 line
    if ra_lims[0]<ra_lims[1]:
        sources_in_bounds =  [i for i in arange(len(RA)) if (RA[i]>=ra_lims[0] and RA[i]<=ra_lims[1]) and (Dec[i]>=dec_lims[0] and Dec[i]<=dec_lims[1])]
        area = get_lune(ra_lims[0],ra_lims[1],dec_lims[0],dec_lims[1])
    ##If they do, searching the coords slightly differently and do some rearranging to get the area
    else:
        sources_in_bounds =  [i for i in arange(len(RA)) if (RA[i]>=ra_lims[0] or RA[i]<=ra_lims[1]) and (Dec[i]>=dec_lims[0] and Dec[i]<=dec_lims[1])]
        extra = 360.0 - ra_lims[0]
        area = get_lune(0,ra_lims[1]+extra,dec_lims[0],dec_lims[1])
    scaled_source_density = (4*pi*len(sources_in_bounds))/area

    bound_lims = '%.5f,%.5f,%.5f,%.5f' %(lower_ra,upper_ra,min_dec,max_dec)

    ##Create a new table, and populate with the data in correct units
    t=Table(masked=True,meta={'src_dens':scaled_source_density,'bound_lims':bound_lims})
    t_names = Column(name='%s_name' %prefix,data=names,description='Name from catalogue',dtype=str)
    t_ras = Column(name='%s_RAJ2000' %prefix,data=RA,description='Right Ascension of source (J2000)',unit='deg',dtype=float)
    t_rerrs = Column(name='%s_e_RAJ2000' %prefix,data=RA_error,description='Error on Right Ascension',unit='deg',dtype=float)
    t_decs = Column(name='%s_DEJ2000' %prefix,data=Dec,description='Declination of source (J2000)',unit='deg',dtype=float)
    t_derrs = Column(name='%s_e_DEJ2000' %prefix,data=Dec_error,description='Error on Declination',unit='deg',dtype=float)
    t_fluxes = Column(name='%s_S%d' %(prefix,float(freq)),data=flux,description='Source flux at %.1fMHz' %float(freq),unit='Jy',dtype=float)
    t_ferrs = Column(name='%s_e_S%d' %(prefix,float(freq)),data=flux_error,description='Error on flux at %.1fMHz' %float(freq),unit='Jy',dtype=float)
    t_majors = Column(name='%s_MajAxis' %prefix,data=major,description='Fitted major axis',unit='deg',dtype=float)
    t_minors = Column(name='%s_MinAxis' %prefix,data=minor,description='Fitted minor axis',unit='deg',dtype=float)
    t_PAs = Column(name='%s_PA' %prefix,data=PA,description='Fitted Position Angle',unit='deg',dtype=float)
    mask = []
    for i in flags:
        if type(i)==ma.core.MaskedConstant: mask.append(True)
        else: mask.append(False)
    t_flags = MaskedColumn(name='%s_flag' %prefix,data=flags,description='Any meta flag for inclusion',mask=mask,fill_value='--',dtype=str)
    mask = []
    for i in ID:
        if type(i)==ma.core.MaskedConstant: mask.append(True)
        else: mask.append(False)
    t_fields = MaskedColumn(name='%s_FieldID' %prefix,data=ID,description='If avaiable, image field ID',mask=mask,fill_value='--',dtype=str)

    t.add_columns([t_names,t_ras,t_rerrs,t_decs,t_derrs,t_fluxes,t_ferrs,t_majors,t_minors,t_PAs,t_flags,t_fields])
    #Again, handles multiple frequencies in one catalogues
    if len(freqs)>0:
        for i in arange(len(freqs)):
            t_fextra = Column(name='%s_Flux_%.1f' %(prefix,float(freqs[i])),data=fluxs[i],description='Source flux at %.1fMHz' %float(freqs[i]),unit='Jy',dtype=float)
            t_ferrextra = Column(name='%s_Flux_%.1f_err' %(prefix,float(freqs[i])),data=ferrs[i],description='Error on flux at %.1fMHz' %float(freqs[i]),unit='Jy',dtype=float)
            t.add_columns([t_fextra,t_ferrextra])

    ##Add the source density to the table
    #t.add_keyword('%s_nu' %prefix,str(scaled_source_density))
    t.write('simple_%s.fits' %prefix,overwrite=True,format='fits')

    return scaled_source_density,bound_lims


##Read in all of the user inputs and use them to create the simple tables
prefix1 = options.prefix1
prefix2 = options.prefix2
details1 = options.details1.split(',')
details2 = options.details2.split(',')
units1 = options.units1.split(',')
units2 = options.units2.split(',')
ra_lims1 = list(map(float,options.ra_lims1.split(',')))
ra_lims2 = list(map(float,options.ra_lims2.split(',')))
dec_lims1 = list(map(float,options.dec_lims1.split(',')))
dec_lims2 = list(map(float,options.dec_lims2.split(',')))

def read_table(name):
    '''Read in a table and work out if its .vot or .fits, and grab the
    data accordingly. If not the right format, exit the whole process'''
    if '.vot' in name:
        table = vot_parse(name,pedantic=False).get_first_table()
        data = table.array
        return data
    elif '.fits' or '.FITS' in name:
        table = fits.open(name,pedantic=False)
        return table[1].data
    else:
        sys.exit('Entered table must either be VOTable or FITS')

##Read in the table data
data1 = read_table(options.table1)
data2 = read_table(options.table2)

ss_dens1,bound_lims_1 = make_table(data1,details1,units1,prefix1,ra_lims1,dec_lims1)
ss_dens2,bound_lims_2 = make_table(data2,details2,units2,prefix2,ra_lims2,dec_lims2)

##Using the appropriate user values, run stilts to produce a cross matched
##table with all possible matches and all possible information included
##First create a string with all the commands in
match_string = 'stilts tmatch2 join=1and2 find=all matcher=sky params=%d ' %float(options.separation)
match_string += "in1=simple_%s.fits values1='%s_RAJ2000 %s_DEJ2000' suffix1='' " %(prefix1,prefix1,prefix1)
match_string += "in2=simple_%s.fits values2='%s_RAJ2000 %s_DEJ2000' suffix2='' " %(prefix2,prefix2,prefix2)
match_string += "fixcols=all out=matched_%s_%s.fits" %(prefix1,prefix2)

##Run the stilts match
subprocess.call(match_string,shell=True)

if options.make_table:
    pass
else:
    cwd = os.getcwd()
    subprocess.call('rm %s/simple_%s.fits %s/simple_%s.fits' %(cwd,prefix1,cwd,prefix2),shell=True)

##Add the calculated scaled source density of each catalogue to
###the stilts matched catalogue for use by calculate_bayes.py
stilts_table = fits.open("matched_%s_%s.fits" %(prefix1,prefix2))
match_data = stilts_table[1].data

match_table = Table(data=match_data,meta={"nu_1":ss_dens1,"nu_2":ss_dens2,'bound_lims_1':bound_lims_1,'bound_lims_2':bound_lims_2})
match_table.write("matched_%s_%s.fits" %(prefix1,prefix2),overwrite=True,format='fits')
