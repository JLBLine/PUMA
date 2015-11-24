#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import patches
from matplotlib.patches import Ellipse
from itertools import combinations
from astropy.wcs import WCS
from wcsaxes import WCSAxes
import make_table_lib as mkl

dr = np.pi/180.0

##Variables that get used by most of the functions here
global closeness
global high_prob
global low_prob
global chi_thresh
global jstat_thresh
global num_freqs
global split
global matched_cats
global save_plots

##Set up a load of markers, colours and alpha values
markers = ['o','*','s','^','D','8','H','>','<','8','v','d']
marker_sizes = np.array([8,11,8,8,7,8,11,11,11,11,11,11,11]) + 3
marker_colours = ['m', "b", "#e5e500", "r", "g", 'k', "#c0531f",'#660066','#000099','y','#990000','#003300']
ell_colours1 = ["#C370C8", "#698ACD", "#AA6600", "#D84C77", "#5FC34F", 'k', "#D46733",'m','b','y','r','g']
ell_colours2 = ['m', "b", "#8B5A00", "r", "g", 'k', "#c0531f",'#660066','#000099','y','#990000','#003300']
alphas = [0.4,0.4,0.4,0.4,0.4,0.4,0.3,0.45,0.5,0.4,0.5,0.4]


##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
def plot_errors(style,colour,freq,flux,ferr,name,size,ax):
	'''Plots errorbars and markers with no line'''
	ax.errorbar(freq,flux,ferr,
	marker=style,ms=size,mfc=colour,mec=colour,ecolor=colour,markeredgewidth=1,label=name,linestyle='None')
	
def plot_pos(style,colour,ra,dec,rerr,derr,name,size,ax,proj):
	'''Plots a single point with x and y erros bars'''
	if proj==1.0:
		p = ax.errorbar(ra,dec,derr,rerr,marker=style,ms=size,mfc=colour,mec=colour,ecolor=colour,markeredgewidth=1,label=name,linestyle='None')
	else:
		p = ax.errorbar(ra,dec,derr,rerr,marker=style,ms=size,mfc=colour,mec=colour,ecolor=colour,markeredgewidth=1,label=name,linestyle='None',transform=proj)
	return p

def plot_errors_comb(style,colour,freq,flux,ferr,name,size,ax):
	'''Plots errorbars and markers with no line'''
	ax.errorbar(freq,flux,ferr,
	marker=style,ms=size,mfc=colour,mec='k',ecolor=colour,markeredgewidth=1,label=name,linestyle='None')	
	
def plot_pos_comb(style,colour,ra,dec,rerr,derr,name,size,ax,proj):
	'''Plots a single point with x and y erros bars, with a black border around it'''
	if proj==1.0:
		p = ax.errorbar(ra,dec,derr,rerr,marker=style,ms=size,mfc=colour,mec='k',ecolor=colour,markeredgewidth=1.2,label=name,linestyle='None')
	else:
		p = ax.errorbar(ra,dec,derr,rerr,marker=style,ms=size,mfc=colour,mec='k',ecolor=colour,markeredgewidth=1.2,label=name,linestyle='None',transform=proj)
	return p
	
def plt_ell_empty(ra,dec,height,width,PA,ax,colour,colour2,alpha,proj):
	'''Plots an ellipse - either plots on the ax_main which uses a wcs projection
	or on the smaller subplots which don't need transforming'''
	##Position Angle measures angle from direction to NCP towards increasing RA (east)
	##Matplotlib plots the angle from the increasing y-axis toward DECREASING x-axis
	##so have to put in the PA as negative
	if proj==1.0:
		ell = Ellipse([ra,dec],width=width,height=height,angle=-PA,linewidth=1.5)  ##minus???
		ell.set_facecolor('none')
		ell.set_edgecolor(colour2)
		ax.add_artist(ell)
		ell = Ellipse([ra,dec],width=width,height=height,angle=-PA,linewidth=1.5)
		ell.set_facecolor(colour)
		ell.set_alpha(0.3)
		ax.add_artist(ell)
	else:
		ell = Ellipse([ra,dec],width=width,height=height,angle=-PA,transform=proj,linewidth=1.5)  ##minus???
		ell.set_facecolor('none')
		ell.set_edgecolor(colour2)
		ax.add_artist(ell)
		ell = Ellipse([ra,dec],width=width,height=height,angle=-PA,transform=proj,linewidth=1.5)
		ell.set_facecolor(colour)
		ell.set_alpha(0.3)
		ax.add_artist(ell)
	
##POSSIBLE EXTENSION - MAKE GENERIC SO IT CYCLES THROUGH SOME COLOURS, NOT SPECIFIED COLOURS
##FOR A PARTICULAR CATALOGUE
def plot_all(cat,name,ra,rerr,dec,derr,major,minor,PA,ax,proj):
	##Plot the colour by index of catalogue in matched_cats 
	ind = matched_cats.index(cat)
	plot_pos(markers[ind],marker_colours[ind],ra,dec,rerr,derr,name,marker_sizes[ind],ax,proj)
	if cat=='mrc':
		pass
	else:
		##If no minor/major information, don't plot
		if float(minor)!=-100000.0:
			if float(major)!=-100000.0:
				
				plt_ell_empty(ra,dec,float(major),float(minor),float(PA),ax,ell_colours1[ind],ell_colours2[ind],alphas[ind],proj)

##--------------------------------------------------------------------------------------------------------------------
def plot_ind(match,ax,ind_ax,ax_spectral,ra_bottom,ra_top,dec_bottom,dec_top,dom_crit,comb_crit):
	'''Takes a string of information of a particular combination and uses it
	to create a plot of the single combination, and fit and plot a line to
	the spectral information. Returns the positional probability, fitted paramaters
	and the residuals of the fit'''
	
	##Get the information from the particular match given
	info = match.split()
	indexes = [(14+((i-1)*3)) for i in num_freqs]
	starts = [0]
	for i in xrange(len(indexes)-1): starts.append(sum(indexes[:i+1]))
	fluxs = []
	freqs = []
	ferrs = []
	for j in xrange(len(starts)): 
		ind = starts[j]
		cat = info[ind]
		if cat!='-100000.0':
			#num_freq = num_freqs[j]
			freq = num_freqs[j]
			name = info[ind+1]
			ra = float(info[ind+2])
			rerr = float(info[ind+3])
			dec = float(info[ind+4])
			derr = float(info[ind+5])
			nu = float(info[ind+6])
			#flux = float(info[ind+7])
			#ferr = float(info[ind+8])/flux
			for k in xrange(freq):
				if info[7+ind+(3*k)]!='-100000.0':
					if np.isnan(float(info[7+ind+(3*k)])) == False:
						freqs.append(float(info[6+ind+(3*k)]))
						fluxs.append(float(info[7+ind+(3*k)]))
						ferrs.append(float(info[8+ind+(3*k)])/float(info[7+ind+(3*k)]))
			major = info[ind+9+((freq-1)*3)]
			minor = info[ind+10+((freq-1)*3)]
			PA = info[ind+11+((freq-1)*3)]
			#fluxs.append(flux)
			#freqs.append(nu)
			#ferrs.append(ferr)
			##Plot each source on the individual combo plot
			plot_all(cat,name,ra,rerr,dec,derr,major,minor,PA,ax,1.0)
	##Sort the frequencies, fluxes and log them
	log_fluxs = np.log([flux for (freq,flux) in sorted(zip(freqs,fluxs),key=lambda pair: pair[0])])
	sorted_ferrs = np.array([ferr for (freq,ferr) in sorted(zip(freqs,ferrs),key=lambda pair: pair[0])])
	log_freqs = np.log(sorted(freqs))
	ferrs = np.array(ferrs)
	prob = info[-1]
	
	##Fit a line using weighted least squares and plot it
	lin_fit,jstat,bse,chi_red = mkl.fit_line(log_freqs,log_fluxs,sorted_ferrs)
	
	ax.text(0.5,0.925, 'P$_{%d}$=%.2f' %(ind_ax+1,float(prob)), transform=ax.transAxes,verticalalignment='center',horizontalalignment='center',fontsize=16)
	ax.text(0.5,0.06, '$\epsilon_{%d}$=%.2f $\chi_{%d}$=%.1f' %(ind_ax+1,jstat,ind_ax+1,chi_red),
		transform=ax.transAxes,verticalalignment='center',horizontalalignment='center',fontsize=16)
	ax.set_xlim(ra_bottom,ra_top)
	ax.set_ylim(dec_bottom,dec_top)
	
	##Plot RA backwards
	ax.invert_xaxis()
	
	##Plot the fitted spectral line, and return the plot object so we can create a legend from it
	if dom_crit=='No dom. source':
		if 'split' in comb_crit:
			spec_plot = 'na'
		else:
			spec_plot, = ax_spectral.plot(np.exp(log_freqs),np.exp(lin_fit.fittedvalues),linestyle='-',linewidth=1,alpha=0.3)
	else:
		spec_plot, = ax_spectral.plot(np.exp(log_freqs),np.exp(lin_fit.fittedvalues),linestyle='-',linewidth=1,alpha=1)
	
	return prob,jstat,spec_plot,lin_fit.params

def make_left_plots(fig,main_dims,spec_dims,ra_main,dec_main):
	
	##A fits image header with which to create a wcs with
	header = { 'NAXIS'  : 2,             ##Number of data axis
    'NAXIS1' : 10,                  ##Length of X axis
    'CTYPE1' : 'RA---SIN',           ##Projection type of X axis
	'CRVAL1' : ra_main,        ##Central X world coord value
	'CRPIX1' : 5,                    ##Central X Pixel value
	'CUNIT1' : 'deg',                ##Unit of X axes
	'CDELT1' : -0.001*np.cos(dec_main*(np.pi/180.0)),              ##Size of pixel in world co-ord
	'NAXIS2' : 10,                  ##Length of X axis
	'CTYPE2' : 'DEC--SIN',           ##Projection along Y axis
	'CRVAL2' : dec_main,                   ##Central Y world coord value
	'CRPIX2' : 5,                    ##Central Y Pixel value
	'CUNIT2' : 'deg',                ##Unit of Y world coord
	'CDELT2' : +0.001      		     ##Size of pixel in deg
	} 
	
	##Create the ws, and the main axis based on that. Plot top left
	wcs = WCS(header=header)
	ax_main = WCSAxes(fig, main_dims, wcs=wcs)
	fig.add_axes(ax_main)
	tr_fk5 = ax_main.get_transform("fk5")
	
	#ax_main.set_title("All sources within 3'.0")
	ax_main.text(0.01,0.93,"All sources within search area",verticalalignment='bottom',horizontalalignment='left', transform=ax_main.transAxes,fontsize=16)
	
	##Create bottom left plot with log-log axes - set the error bars to plot
	##even if they go off the edge of the plot
	ax_spectral = fig.add_axes(spec_dims)
	ax_spectral.set_xscale("log",nonposx='clip')
	ax_spectral.set_yscale("log",nonposy='clip')
	
	return ax_main,ax_spectral,tr_fk5,wcs

def fill_left_plots(all_info,ra_main,dec_main,ax_main,ax_spectral,tr_fk5,wcs,all_fluxs,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,delta_RA):
	'''Get the information and plot the positions and fluxs on the left hand side plots'''
	ras = []
	decs = []
	all_freqs = []
	all_ferrs = []
	for i in xrange(len(all_info)):
		info=all_info[i].split()
		cat = info[0]
		name = info[1]
		ra = float(info[2])
		ras.append(ra)
		rerr = float(info[3])
		dec = float(info[4])
		decs.append(dec)
		derr = float(info[5])
		major = info[-5]
		minor = info[-4]
		PA = info[-3]
		ID = info[-1]
		##If base catalogue, plot error and mathcing ellipses
		if i==0:
			##Calculate the change in RA for half the user given resolution (and constant dec) using spherical trigo
			error_RA = np.arccos((np.cos(closeness*dr)-np.sin(dec*dr)**2)/np.cos(dec*dr)**2)/dr
			
			##Plot an error ellipse of the base cat error + resolution
			ell = patches.Ellipse((ra_main,dec_main),2*(rerr+error_RA),2*(derr+closeness),angle=0,
				transform=tr_fk5,linestyle='dashed',fc='none',lw=1.1,color='gray')
			ax_main.add_patch(ell)
			##Plot a circle of the match radius
			ell = patches.Ellipse((ra_main,dec_main),2*delta_RA,4*(closeness),angle=0,
				transform=tr_fk5,linestyle='dashdot',fc='none',lw=1.1,color='k')
			ax_main.add_patch(ell)
		
		##Plot positions and elliptical fits
		plot_all(cat,name,ra,rerr,dec,derr,major,minor,PA,ax_main,tr_fk5)
		
		##See if one or more flux for a source, and plot fluxes with errorbars
		cat_ind = matched_cats.index(cat)
		if len(info)==14:
			freq = float(info[6])
			flux = float(info[7])
			ferr = float(info[8])
			plot_errors(markers[cat_ind],marker_colours[cat_ind],freq,flux,ferr,name,marker_sizes[cat_ind],ax_spectral)
			all_fluxs.append(flux)
			all_freqs.append(freq)
			all_ferrs.append(ferr)
		##If the catalogue has more than one frequency:
		else:
			extra = (len(info)-14) / 3
			freqs = []
			fluxs = []
			ferrs = []
			for i in xrange(extra+1):
				if info[7+(3*i)]!='-100000.0' and info[7+(3*i)]!='nan':
					freqs.append(info[6+(3*i)])
					fluxs.append(info[7+(3*i)])
					ferrs.append(info[8+(3*i)])
			
			for i in xrange(len(fluxs)):
				all_fluxs.append(float(fluxs[i]))
				plot_errors(markers[cat_ind],marker_colours[cat_ind],float(freqs[i]),float(fluxs[i]),float(ferrs[i]),'%s-%.1fMHz' %(name,float(freqs[i])),marker_sizes[cat_ind],ax_spectral)
			for freq in freqs: all_freqs.append(float(freq))
			for ferr in ferrs: all_ferrs.append(float(ferr))
				
				
	##Add some labels and coord formatting to ax_main
	ra_ax = ax_main.coords[0]
	dec_ax = ax_main.coords[1]
	ra_ax.set_axislabel('RAJ2000')
	dec_ax.set_axislabel('DECJ2000')
	ra_ax.set_major_formatter('hh:mm:ss')
	dec_ax.set_major_formatter('dd:mm:ss')
	
	##Convert axes limits to ax_main wcs, and apply
	ra_low = wcs.wcs_world2pix(ra_down_lim,dec_main,0)  ##The zero is for the orgin point of the image
	ra_high = wcs.wcs_world2pix(ra_up_lim,dec_main,0)
	dec_low = wcs.wcs_world2pix(ra_main,dec_down_lim,0)
	dec_high = wcs.wcs_world2pix(ra_main,dec_up_lim,0)
	ax_main.set_ylim(dec_low[1],dec_high[1])
	ax_main.set_xlim(ra_high[0],ra_low[0])
	
	return all_freqs

def scale_spectral(all_fluxs,all_freqs,ax_spectral):
	##Make the labels on ax_spectral print in MHz and Jy
	max_lflux = np.log(max(all_fluxs))
	min_lflux = np.log(min(all_fluxs))
	
	freq_ticks = [freq for freq in sorted(set(all_freqs))]
	flux_ticks = np.exp(np.arange(min_lflux,max_lflux+abs(max_lflux-min_lflux)/5,abs(max_lflux-min_lflux)/5))
	ax_spectral.set_xticks(freq_ticks,minor=False)
	ax_spectral.set_xticklabels(freq_ticks)
	ax_spectral.set_yticks(flux_ticks,minor=False)
	ax_spectral.set_yticklabels(['%.3f' %flux for flux in list(flux_ticks)],fontsize=14.0)

	##Set some limits on the spextral axis - need to do it in logspace,
	##and set the edge gaps based on the data
	freq_delta = abs(np.log10(max(all_freqs))-np.log10(min(all_freqs)))/10.0
	flux_delta = abs(np.log10(max(all_fluxs))-np.log10(min(all_fluxs)))/10.0
	
	freq_max = 10**(np.log10(max(all_freqs))+freq_delta)
	freq_min = 10**(np.log10(min(all_freqs))-freq_delta)
	flux_max = 10**(np.log10(max(all_fluxs))+flux_delta)
	flux_min = 10**(np.log10(min(all_fluxs))-flux_delta)
	
	ax_spectral.set_xlim(freq_min,freq_max)
	ax_spectral.set_ylim(flux_min,flux_max)

	##Stick some grid stuff on the log log plot
	ten_steps = []
	for arr in [np.array([1e-5,2e-5,3e-5,4e-5,5e-5,6e-5,7e-5,8e-5,9e-5])*10**x for x in xrange(10)]:
		for i in list(arr): ten_steps.append(i)

	ax_spectral.set_xticks([step for step in ten_steps if step>freq_min and step<freq_max],minor=True)
	ax_spectral.set_yticks([step for step in ten_steps if step>flux_min and step<flux_max],minor=True)

	ax_spectral.xaxis.grid(True, which='minor',linestyle='dashed',alpha=0.1)
	ax_spectral.yaxis.grid(True, which='minor',linestyle='dashed',alpha=0.1)
	ax_spectral.set_xlabel(r'log$_{10}$(Frequency) (MHz)',fontsize=14)
	ax_spectral.set_ylabel(r'log$_{10}$(Flux) (Jy)',fontsize=14)
	

def create_plot(comp,accepted_inds,match_crit,dom_crit,outcome):
	'''The main plotting function that takes the relevant data and plots the outcome'''
	###Split the information up as needed
	chunks = comp.split('START_COMP')
	all_info = chunks[0].split('\n')
	
	##FOR SOME REASON CAN'T DO BOTH OF THESE LINES IN THE SAME FOR LOOP?!?!?!
	for entry in all_info:   
		if entry=='': del all_info[all_info.index(entry)]
	for entry in all_info:
		if 'START' in entry: del all_info[all_info.index(entry)]

	matches = chunks[1].split('\n')
	del matches[0],matches[-2:]

	##See how many matches there are, and set up the number of plots needed. If there are
	##more than 16 matches, only plot the top 15 and print how many more matches there 
	##were - saves heaps of tiny little plots from appearing
	num_matches = len(matches)
	skip_16 = 'no'
	if num_matches==1:
		width = 1
		height = 2
	elif num_matches>16:
		width = 4
		height = 4
		skip_16 = 'yes'
	else:
		width = int(num_matches**0.5)
		height = num_matches/width
		if num_matches%width==0:
			pass
		else:
			height+=1
	
	##Sets up a grid layout for the whole of the figure. We'll use half later on for
	##the individual plots
	gs = gridspec.GridSpec(height,2*width)
	
	##Need a big plot!
	fig = plt.figure(figsize = (18,11))

	##Find out the ra,dec of the base catalogue source
	info=all_info[0].split()
	ra_main = float(info[2])
	dec_main = float(info[4])
	
	##Set up dedicated left plots
	main_dims = [0.16, 0.5, 0.29, 0.35]
	spec_dims = [0.16, 0.1, 0.29, 0.35]
	ax_main,ax_spectral,tr_fk5,wcs = make_left_plots(fig,main_dims,spec_dims,ra_main,dec_main)
	
	##Find the limits out to search area - have to do each edge individual,
	##because of the arcdistance projection malarky
	##Even at the same dec, 3 arcmins apart in RA doesn't translate to 3arcmin arcdist - projection
	##fun. Can use law of cosines on a sphere to work out appropriate delta RA. Use this to define plot
	##limits for a nice looking plot
	delta_RA = np.arccos((np.cos((2*closeness)*dr)-np.sin(dec_main*dr)**2)/np.cos(dec_main*dr)**2)/dr
	plot_lim = (2*closeness) + (0.1/60.0)
	ra_up_lim = ra_main + delta_RA + (0.1/60.0)
	ra_down_lim = ra_main - delta_RA - (0.1/60.0)
	dec_up_lim = dec_main + plot_lim
	dec_down_lim = dec_main - plot_lim
	
	##Plot the individual combination plots - do this first so the error bars go over
	##the top of the line plots
	spec_labels = []
	SIs = []
	for i in xrange(height):
		for j in range(width,2*width):
			if i*j == 21:
				if skip_16=='yes':
					ax = plt.subplot(gs[i,j])
					ax.set_xticklabels([])
					ax.set_yticklabels([])
					ax.text(0.5,0.5,"And %d\nother\nplots" %(num_matches-15),transform=ax.transAxes,verticalalignment='center',horizontalalignment='center',fontsize=16)
			else:
				try:
					ind = (i*width)+(j-width)
					match = matches[ind]
					ax = plt.subplot(gs[i,j])
					ax.set_xticklabels([])
					ax.set_yticklabels([])
					##TODO - if plot centred on or close to RA,Dec = 0,0 then going to get wrapping problems. Should be able to pull the need
					##for a wrap from ra_down_lim,ra_up_lim - one should be <0.0, or >360.0. Need to happen inside plot_ind
					prob,resid,spec_plot,params = plot_ind(match,ax,ind,ax_spectral,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,dom_crit,outcome)
					if spec_plot=='na':
						pass
					else:
						SIs.append([params[0],str(ind+1)])
						spec_labels.append(spec_plot)
				except IndexError:
					pass
	
	#===========================================================#
	##Plot the matching criteria information
	match1 = matches[0].split()
	src_g = mkl.get_srcg(match1)
	
	text_axes = fig.add_axes([0.45,0.5,0.125,0.35])
	text_axes.axis('off')

	##Plot the matching information
	props = dict(boxstyle='round', facecolor='w',lw='1.5')
	text_axes.text(0.5,0.5,'Match Criteria:\n%s\n\nDominace Test:\n%s\n\nOutcome:\n%s' %(match_crit,dom_crit,outcome),
		bbox=props,transform=text_axes.transAxes,verticalalignment='center',horizontalalignment='center',fontsize=16)
	
	all_fluxs = []
	##Fill the left hand plots with information goodness
	all_freqs = fill_left_plots(all_info,ra_main,dec_main,ax_main,ax_spectral,tr_fk5,wcs,all_fluxs,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,delta_RA)
	
	##If no repeated catalogues to combine, skip
	if num_matches==0 or num_matches==1:
		pass
	##Otherwise, plot the combined fluxes
	else:
		##Calculate and plot the combined fluxes of the two sources, even if source one or two has been accepted
		##just as a guide
		src_all = mkl.get_allinfo(all_info)
		
		##If positionally impossible don't plot combined info
		if accepted_inds=='Nope':
			pass
		##If only one position possible, don't plot combined info
		elif len(accepted_inds) == 1:
			pass
		##Otherwise, see what the combined fluxes look like
		else:
			comb_crit, ra_ws, rerr_ws, dec_ws, derr_ws, temp_freqs, comb_freqs, comb_fluxs, comb_ferrs, comb_fit, comb_jstat, comb_chi_red, combined_names, set_freqs, set_fluxs, set_fits = mkl.combine_flux(src_all,src_g,accepted_inds,'plot=yes',len(matches))
		
		##If the criteria sent the double to be combined, actually plot the fitted line
		if dom_crit == 'No dom. source':
			
			for freq,flux in zip(set_freqs,set_fluxs):
				ax_spectral.plot(freq,flux,linestyle='--',linewidth=1,color='r')
			
			split_colors = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']
			for fit in set_fits:
				ind = set_fits.index(fit)
				ax_spectral.plot(set_freqs[ind],set_fluxs[ind],linestyle='--',linewidth=1.0,color=split_colors[ind],alpha=0.7)
				split_p, = ax_spectral.plot(temp_freqs,np.exp(fit.params[1] + np.log(temp_freqs)*fit.params[0]),linestyle='-',linewidth=1.5,color=split_colors[ind])
				spec_labels.append(split_p)
				SIs.append([fit.params[0],'split %d' %(ind+1)])
			
			bright_colours = ['#FF6600','#33FF33','#FF47A3','#00ebb3']
			
			for freq in xrange(len(comb_freqs)):
				plot_errors_comb('*',bright_colours[freq],comb_freqs[freq],comb_fluxs[freq],comb_ferrs[freq],'combo',16,ax_spectral)
			comb_p, = ax_spectral.plot(temp_freqs,np.exp(comb_fit.fittedvalues),linestyle='--',linewidth=1.5,color='k')
			spec_labels.append(comb_p)
			SIs.append([comb_fit.params[0],'comb'])
			
			##Send the combined fluxes to the all_fluxs so that ax_spectral is scaled appropriately
			for flux in comb_fluxs:
				all_fluxs.append(flux)
			
			for pos in xrange(len(ra_ws)):
				patch = plot_pos_comb('*',bright_colours[pos],ra_ws[pos],dec_ws[pos],rerr_ws[pos],derr_ws[pos],combined_names[pos],16,ax_main,ax_main.get_transform("fk5"))
				
	scale_spectral(all_fluxs,all_freqs,ax_spectral)

	##==============================================================

	fig.tight_layout()
	fig.subplots_adjust(bottom=0.1)
	fig.subplots_adjust(left=0.15)
	
	##Make room at the top of the plot for a legend for ax_main, make the legend
	fig.subplots_adjust(top=0.85)
	
	leg_labels = [r'$\alpha_{%s}$ = %.2f' %(SI[1],SI[0]) for SI in SIs]
	main_handles,main_labels = ax_main.get_legend_handles_labels()
	
	main_leg = fig.add_axes([0.05,0.87,0.9,0.05])
	main_leg.axis('off')
	main_leg.legend(main_handles,main_labels,loc='lower center',prop={'size':12},ncol=8) #,bbox_to_anchor=(0,1.02),
	
	spec_leg = fig.add_axes([0.45,0.1,0.125,0.35])
	spec_leg.axis('off')

	##Stop the legend from having so many entries that it goes off the plot	
	if len(spec_labels)>11:
		trim_labels = spec_labels[:10]
		trim_labels.append(spec_labels[-1])
		trim_legs = leg_labels[:10]
		trim_legs.append(leg_labels[-1])
		spec_leg.legend(trim_labels,trim_legs,loc='center',prop={'size':14},fancybox=True)
	else:
		spec_leg.legend(spec_labels,leg_labels,loc='center',prop={'size':14},fancybox=True)
	
	##Create an axes to contain patches for an ellipse legend
	patch_leg = fig.add_axes([0.015,0.1,0.06,0.75])
	patch_leg.set_xticks([])
	patch_leg.set_yticks([])
	patch_leg.set_xticklabels([])
	patch_leg.set_yticklabels([])
	
	##See what catalogues are present in the match
	present_cats = [cat for cat in set(src_g.cats) if cat!='-100000.0']
	##Scale accordingly
	increment = 1.0/(2+len(present_cats))
	ell_positions = np.arange(increment/2,1,increment)
	
	##Find the axes coord transform
	patch_trans = patch_leg.transAxes
	##Plot and name the resolution ellipse
	ell = patches.Ellipse((0.5,ell_positions[-2]),0.9,increment-0.05,angle=0,
		transform=patch_trans, linestyle='dashed',fc='none',lw=1.1,color='gray')
	patch_leg.add_patch(ell)
	patch_leg.text(0.5,ell_positions[-2],'Resolution\n+ error',
		transform=patch_trans,verticalalignment='center',horizontalalignment='center',fontsize=14)
	
	##Plot and name the search ellipse
	ell = patches.Ellipse((0.5,ell_positions[-1]),0.9,increment-0.05,angle=0,
		transform=patch_trans, linestyle='dashdot',fc='none',lw=1.1,color='k')
	patch_leg.add_patch(ell)
	patch_leg.text(0.5,ell_positions[-1],'Search\nradius',
		transform=patch_trans,verticalalignment='center',horizontalalignment='center',fontsize=14)
	
	##Use the same method as plot_all - for some reason was getting transform errors.
	##so do it separately here (sigh)
	for cat in present_cats:
		col_ind = matched_cats.index(cat)
		position_ind = present_cats.index(cat)
		patch_leg.errorbar(0.5,ell_positions[position_ind],0.01,0.075,marker=markers[col_ind],ms=8,mfc=marker_colours[col_ind],
			mec=marker_colours[col_ind],ecolor=marker_colours[col_ind],markeredgewidth=1,label='meh',linestyle='None',transform=patch_trans)
		
		ell = patches.Ellipse((0.5,ell_positions[position_ind]),0.9,increment-0.05,angle=0, transform=patch_trans,
			fc=ell_colours1[col_ind],color=ell_colours2[col_ind],alpha=alphas[col_ind])
		patch_leg.add_patch(ell)
		
		patch_leg.text(0.5,ell_positions[position_ind]-(increment/2-0.04),cat,
			transform=patch_trans,verticalalignment='center',horizontalalignment='center',fontsize=16)
	
	if save_plots:
		plt.savefig('%s-pumaplot.png' %all_info[0].split()[1],bbox_inches='tight',dpi=100)
		plt.close()
			
	else:
		plt.show()