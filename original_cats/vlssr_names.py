import atpy
from numpy import *
import subprocess

def deg_to_degmins(x):	    #converts angle degrees form in to dd:mm:ss.ss
	x=float(x)
	deg=abs(x)
	degr=deg-int(deg)
	mins=(degr)*60.00
	secs=(mins-int(mins))*60.0
	if mins!=0:
		if -1e-5<=(secs-60)<1e-5:
			mins=mins+1
			secs=0.0
	if x>=0:
		return '+%02d%02d%02d' %(int(deg),int(mins),int(secs))
	if x<0:
		return '-%02d%02d%02d' %(int(deg),int(mins),int(secs))

def deg_to_hour(x):    #converts angle in degrees in to hh:mm:ss.ss, must input as a string
	x=float(x)
	deg=abs(x)
	hr=deg/15.0
	mins=(hr-int(hr))*60.0
	secs=(mins-int(mins))*60.0
	if mins!=0:
		if -1e-5<=(secs-60)<1e-5:
			mins=mins+1
			secs=0.0
	if x>=0:
		return '%02d%02d%02d' %(int(hr),int(mins),int(secs))
	if x<0:
		return '-%02d:%02d:%08.5f' %(int(hr),int(mins),secs)
	
def hour_to_deg(time):      #converts hh:mm:ss.ss in to degrees, must input as a string
	negtest=time[0]
	time=time.split(':')
	degr=float(time[0])*15.0
	mins=float(time[1])*(15.0/60.0)
	secs=float(time[2])*(15.0/(3600.0))
	if negtest=='-':
		deg=degr-mins-secs
	if negtest!='-':
		deg=degr+mins+secs
	return deg

def dec_to_deg(time):    #converts dd:mm:ss.ss in to degrees, must input as a string
	negtest=time[0]
	time=time.split(':')
	degr=float(time[0])
	mins=float(time[1])*(1.0/60.0)
	secs=float(time[2])*(1.0/(3600.0))
	if negtest=='-':
		deg=degr-mins-secs
	if negtest!='-':
		deg=degr+mins+secs
	return deg

class source_info:
	def __init__(self):
		self.name = None
		self.ra = None
		self.dec = None
		self.rerr = None
		self.derr = None
		self.SI = None
		self.intercept = None
		self.SI_err = None
		self.intercept_err = None
		self.num_match = None
		self.retained_match = None
		self.type_match = None
		self.low_resid = None
		self.f180 = None
		self.f180_ext = None
		self.num_cats = None
		
ras = []
rerrs = []
decs = []
derrs = []
fluxs = []
ferrs = []
majors = []
minors = []
PAs = []
fields = []
flags = []

vlssr_data = open('FullVLSSCatalog.text').read().split('\n')

def fill_full(hh,mm,ss,dd,dm,ds,ori,flux,major,minor,PA,flag,field,X,Y):
	ra = hour_to_deg('%s:%s:%s' %(hh,mm,ss))
	dec = dec_to_deg('%s:%s:%s' %(dd,dm,ds))
	ras.append(ra)
	decs.append(dec)
	fluxs.append(float(flux))
	
	if '<' in major:
		major = major[1:]
	if '<' in minor:
		minor = minor[1:]
	
	majors.append(float(major))
	minors.append(float(minor))
	PAs.append(float(PA))
	fields.append(field)
	flags.append(flag)

for line in vlssr_data:
	if 'RA' in line:
		pass
	elif '#' in line:
		pass
	elif 'NVSS' in line:
		pass
	elif 'deg' in line:
		pass
	else:
		info = line.split()
		if len(info) == 14:
			hh,mm,ss,dd,dm,ds,ori,flux,major,minor,PA,field,X,Y = info
			fill_full(hh,mm,ss,dd,dm,ds,ori,flux,major,minor,PA,'',field,X,Y)
		elif len(info) == 15:
			hh,mm,ss,dd,dm,ds,ori,flux,major,minor,PA,flag,field,X,Y = info
			fill_full(hh,mm,ss,dd,dm,ds,ori,flux,major,minor,PA,flag,field,X,Y)
		else:
			rerrs.append(float(info[0])*(15.0/3600.0))
			derrs.append(float(info[1])*(1.0/3600.0))
			ferrs.append(float(info[3]))

tdata = atpy.Table(masked=True)

names = ['J%s%s' %(deg_to_hour(ras[i]),deg_to_degmins(decs[i])) for i in xrange(len(ras))]
			
tdata.add_column('Name', array(names),description='Name based on J2000 position',dtype=str)
tdata.add_column('RA_J2000',array(ras),description='J2000 Right ascension of source',unit='deg',dtype=float)
tdata.add_column('RA_err',array(rerrs),description='Error on Right ascension of source',unit='deg',dtype=float)
tdata.add_column('DEC_J2000',array(decs),description='J2000 Declination of source',unit='deg',dtype=float)
tdata.add_column('DEC_err',array(derrs),description='Error on Declination of source',unit='deg',dtype=float)
tdata.add_column('Flux',array(fluxs),description='Flux density',unit='Jy',dtype=float)
tdata.add_column('Flux_err',array(ferrs),description='Error on flux density',unit='Jy',dtype=float)
tdata.add_column('major',array(majors),description='Major axes of  gaussian fit',unit='arcsec',dtype=float)
tdata.add_column('minor',array(minors),description='Minor axes of  gaussian fit',unit='arcsec',dtype=float)
tdata.add_column('PA',array(PAs),description='Position angle of gaussian fit',unit='deg',dtype=float)
tdata.add_column('Flags',array(flags),description='If extended, vlssr includes a flag',unit='deg',dtype=str)
tdata.add_column('Field',array(fields),description='The image field ID',unit='deg',dtype=str)
	
tdata.write("vlssr_names.vot",votype='ascii',overwrite=True)
subprocess.call('stilts tcopy ifmt=votable ofmt=fits in=sumss_names.vot out=sumss_names.fits',shell=True)
