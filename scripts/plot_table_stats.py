#!/usr/bin/python
import atpy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import optparse
from astroML.plotting import hist
from statsmodels.robust.scale import mad
import statsmodels.api as sm

parser = optparse.OptionParser()

parser.add_option('-v', '--input_vot', 
	help='Enter name of input .vot table')

options, args = parser.parse_args()

input_vot = options.input_vot

class source_info:
	def __init__(self):
		self.SI = None
		self.intercept = None
		self.SI_err = None
		self.intercept_err = None
		self.num_match = None
		self.retained_match = None
		self.type_match = None
		self.low_resid = None
		self.num_cats = None

tdata = atpy.Table(input_vot,verbose=False)

SIs = tdata['SI']
intercepts = tdata['Intercept']
SI_errs = tdata['e_SI']
intercept_errs = tdata['e_Intercept']
num_cats = tdata['Number_cats']
num_matches = tdata['Number_matches']
retained_matches = tdata['Retained_matches']
type_matches = tdata['Match_stage']
low_resids = tdata['Low_resids']

def make_source(ind):
	source=source_info()
	source.SI = SIs[ind]
	source.intercept = intercepts[ind]
	source.SI_err = SI_errs[ind] 
	source.intercept_err = intercept_errs[ind]
	source.num_match = num_matches[ind]
	source.retained_match = retained_matches[ind]
	source.type_match = type_matches[ind]
	source.low_resid = low_resids[ind]
	source.num_cats = num_cats[ind]
	return source

##A method to plot the univariate kernal density estimate of a distribution as a filled lineplot
def plot_by_kde(ax,data,colour,linewidth,label,linestyle):
	kde = sm.nonparametric.KDEUnivariate(data)
	kde.fit()
	ax.fill(kde.support, kde.density,  facecolor=colour, alpha=0.2,edgecolor=colour)
	ax.plot(kde.support, kde.density, color=colour,linewidth=linewidth,label=label,linestyle=linestyle)
	
##Populate all the sources
sources = []
for i in xrange(len(SIs)): sources.append(make_source(i))

##Plot the SI dist of good vs bad fits-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
from matplotlib import rc
font = {'size': 14}
rc('font', **font)

fig_hist = plt.figure(figsize=(18,12))

colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']

##Sometimes if the table contains single sources, there is no SI fit, so column contains NaNs
## atpy reads as '--' so need to avoid these
sources_SIs = [source for source in sources if source.SI!='--']
SIs = [float(source.SI) for source in sources_SIs]

##Plot all of the SIs together
##-----------------------------------------------------------------------------------------------------------------------
ax1 = fig_hist.add_subplot(221)
plot_by_kde(ax1,SIs,'k',3.0,'All fits (%d sources)' %len(SIs),'-')

##Compare the good fits to the bad fits
##-----------------------------------------------------------------------------------------------------------------------
ax2 = fig_hist.add_subplot(222)
good_fit_SIs = [float(source.SI) for source in sources_SIs if float(source.low_resid)==0]
bad_fit_SIs = [float(source.SI) for source in sources_SIs if float(source.low_resid)==1]

plot_by_kde(ax2,good_fit_SIs,colours[0],3.0,'$\chi^2_{red}<=2.0$\n(%d sources)' %len(good_fit_SIs),'-')
plot_by_kde(ax2,bad_fit_SIs,colours[3],3.0,'$\chi^2_{red}>2.0$\n(%d sources)' %len(bad_fit_SIs),'--')

##Compare the matches with just one matched frequency to the base catalogue, to those
#with multiple frequencies matched
##-----------------------------------------------------------------------------------------------------------------------
ax3 = fig_hist.add_subplot(223)
two_cat_SIs = [float(source.SI) for source in sources_SIs if source.num_cats==2]
nottwo_cat_SIs = [float(source.SI) for source in sources_SIs if source.num_cats>2]

plot_by_kde(ax3,two_cat_SIs,colours[1],3.0,'Only two \ncatalogues \n(%d sources)' %len(two_cat_SIs),'-')
plot_by_kde(ax3,nottwo_cat_SIs,colours[2],3.0,'More than two \ncatalogues \n(%d sources)' %len(nottwo_cat_SIs),'--')

##Compare the SIs for the different types of matches performed by PUMA
##-----------------------------------------------------------------------------------------------------------------------
ax4 = fig_hist.add_subplot(224)
pos_match = [float(source.SI) for source in sources_SIs if source.type_match=='position']
spec_match = [float(source.SI) for source in sources_SIs if source.type_match=='spectral']
comb_match = [float(source.SI) for source in sources_SIs if source.type_match=='combine']
split_match = [float(source.SI) for source in sources_SIs if source.type_match=='splitA']

if len(pos_match)!=0: plot_by_kde(ax4,pos_match,colours[0],3.0,'Positional match \n(%d sources)' %len(pos_match),'-')
if len(spec_match)!=0: plot_by_kde(ax4,spec_match,colours[1],3.0,'Spectral match \n(%d sources)' %len(spec_match),'--')
if len(comb_match)!=0: plot_by_kde(ax4,comb_match,colours[2],5.0,'Combined match \n(%d sources)' %len(comb_match),':')
if len(split_match)!=0: plot_by_kde(ax4,split_match,colours[3],5.0,'Split match \n(%d sources)' %len(split_match),'-.')

##Tidy up the labels and add legends
for ax,label in zip([ax1,ax2,ax3,ax4],[r'$(a)$',r'$(b)$',r'$(c)$',r'$(d)$']):
	ax.set_xlabel('Spectral Index', fontsize=20)
	ax.set_ylabel('Estimated density', fontsize=20)
	ax.legend(loc='upper right',prop={'size':18})
	#ax.axvline(x=-0.8,color='k',linestyle='--')
	ax.text(0.05, 0.95, label, transform=ax.transAxes, fontsize=28,
        verticalalignment='top')#, bbox=dict(boxstyle='none'))
	ax.tick_params(labelsize=16)
	
plt.tight_layout()
plt.show()
