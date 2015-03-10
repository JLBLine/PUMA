#!/usr/bin/python
import atpy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import optparse
from astroML.plotting import hist

parser = optparse.OptionParser()

parser.add_option('-v', '--input_vot', 
	help='Enter name of input .vot table')

options, args = parser.parse_args()

input_vot = options.input_vot

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

tdata = atpy.Table(input_vot,verbose=False)

#names = tdata['mwacs_name']
#ras = tdata['updated_RA_J2000']
#decs = tdata['updated_DEC_J2000']
#rerrs = tdata['updated_RA_err']
#derrs = tdata['updated_DEC_err']
SIs = tdata['SI']
intercepts = tdata['Intercept']
SI_errs = tdata['e_SI']
intercept_errs = tdata['e_Intercept']
num_cats = tdata['Number_cats']
num_matches = tdata['Number_matches']
retained_matches = tdata['Retained_matches']
type_matches = tdata['Type_matches']
low_resids = tdata['Low_resids']
#f180s = tdata['S_180']
#f180_exts = tdata['S_180.0_ext']

def make_source(ind):
	source=source_info()
	#source.name = names[ind]
	#source.ra = ras[ind]
	#source.dec = decs[ind]
	#source.rerr = rerrs[ind]
	#source.derr = derrs[ind]
	source.SI = SIs[ind]
	source.intercept = intercepts[ind]
	source.SI_err = SI_errs[ind] 
	source.intercept_err = intercept_errs[ind]
	source.num_match = num_matches[ind]
	source.retained_match = retained_matches[ind]
	source.type_match = type_matches[ind]
	source.low_resid = low_resids[ind]
	#source.f180 = f180s[ind]
	#source.f180_ext = f180_exts[ind]
	source.num_cats = num_cats[ind]
	return source

##Populate all the sources
sources = []
for i in xrange(len(SIs)): sources.append(make_source(i))

##Plot the SI dist of good vs bad fits-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

from matplotlib import rc
font = {'size': 14}
rc('font', **font)

fig_hist = plt.figure(figsize=(15,10))

colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']

##Plot all of the SIs together
##-----------------------------------------------------------------------------------------------------------------------
ax1 = fig_hist.add_subplot(221)
##Plot the histograms using bayesian blocks 
hist(SIs,bins='blocks',ax=ax1,histtype='step',normed=True,color='k',linewidth=3.0,label='All fits (%d sources)' %len(SIs))

##Compare the good fits to the bad fits
##-----------------------------------------------------------------------------------------------------------------------
ax2 = fig_hist.add_subplot(222)
good_fit_SIs = [source.SI for source in sources if source.low_resid==0]
bad_fit_SIs = [source.SI for source in sources if source.low_resid==1]
hist(good_fit_SIs,bins='blocks',ax=ax2,histtype='step',normed=True,color=colours[0],linewidth=3.0,label='Good fits \n(%d sources)' %len(good_fit_SIs))
hist(bad_fit_SIs,bins='blocks',ax=ax2,histtype='step',normed=True,color=colours[1],linewidth=3.0,label='Bad fits \n(%d sources)' %len(bad_fit_SIs))

##Compare the matches with just one matched frequency to the base catalogue, to those
#with multiple frequencies matched
##-----------------------------------------------------------------------------------------------------------------------
ax3 = fig_hist.add_subplot(223)
two_cat_SIs = [source.SI for source in sources if source.num_cats==2]
nottwo_cat_SIs = [source.SI for source in sources if source.num_cats!=2]
hist(two_cat_SIs,bins='blocks',ax=ax3,histtype='step',normed=True,color=colours[0],linewidth=3.0,label='Only two \ncatalogues \n(%d sources)' %len(two_cat_SIs))
hist(nottwo_cat_SIs,bins='blocks',ax=ax3,histtype='step',normed=True,color=colours[1],linewidth=3.0,label='More than two \ncatalogues \n(%d sources)' %len(nottwo_cat_SIs))

##Compare the SIs for the different types of matches performed by PUMA
##-----------------------------------------------------------------------------------------------------------------------
ax4 = fig_hist.add_subplot(224)
pos_match = [source.SI for source in sources if source.type_match=='position']
spec_match = [source.SI for source in sources if source.type_match=='spectral']
comb_match = [source.SI for source in sources if source.type_match=='combine']
split_match = [source.SI for source in sources if source.type_match=='splitA']
if len(pos_match)!=0: hist(pos_match,bins='blocks',ax=ax4,histtype='step',normed=True,color=colours[0],linewidth=3.0,label='Positional match \n(%d sources)' %len(pos_match))
if len(spec_match)!=0: hist(spec_match,bins='blocks',ax=ax4,histtype='step',normed=True,color=colours[1],linewidth=3.0,label='Spectral match \n(%d sources)' %len(spec_match))
if len(comb_match)!=0: hist(comb_match,bins='blocks',ax=ax4,histtype='step',normed=True,color=colours[2],linewidth=3.0,label='Combined match \n(%d sources)' %len(comb_match))
if len(split_match)!=0: hist(split_match,bins='blocks',ax=ax4,histtype='step',normed=True,color=colours[3],linewidth=3.0,label='Split match \n(%d sources)' %len(split_match))

##Tidy up the labels and add legends
for ax,label in zip([ax1,ax2,ax3,ax4],[r'$(a)$',r'$(b)$',r'$(c)$',r'$(d)$']):
	ax.set_xlabel('Spectral Index', fontsize=20)
	ax.set_ylabel('Number of sources (normed)', fontsize=20)
	ax.legend(loc='best',prop={'size':18})
	ax.axvline(x=-0.8,color='k',linestyle='--')
	ax.text(0.05, 0.95, label, transform=ax.transAxes, fontsize=28,
        verticalalignment='top')#, bbox=dict(boxstyle='none'))
	ax.tick_params(labelsize=16)

plt.show()
