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
	
tdata = atpy.Table('vizier_sumss.vot',verbose=False,tid=0)

ras = tdata['_RAJ2000']
decs = tdata['_DEJ2000']

names = ['J%s%s' %(deg_to_hour(ras[i]),deg_to_degmins(decs[i])) for i in xrange(len(ras))]

tdata.add_column('name',names,description='Name based on J2000 position',dtype=str)

tdata.write("sumss_names.vot",votype='ascii',overwrite=True)
subprocess.call('stilts tcopy ifmt=votable ofmt=fits in=sumss_names.vot out=sumss_names.fits',shell=True)