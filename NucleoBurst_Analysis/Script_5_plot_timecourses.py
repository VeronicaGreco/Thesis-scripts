
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import scipy.optimize as optimize

#############################################################
# SETTINGS/PARAMETERS FOR HOW THE GRAPH LOOKS
#############################################################

# Axes titles
matplotlib.rcParams['axes.labelsize']  = 8
# Numbers on each axis
matplotlib.rcParams['ytick.labelsize']  = 8
matplotlib.rcParams['xtick.labelsize']  = 8
# Space between axis and number
matplotlib.rcParams['ytick.major.pad']  = 0.8
matplotlib.rcParams['ytick.minor.pad']  = 0.8
matplotlib.rcParams['xtick.major.pad']  = 1.5
matplotlib.rcParams['xtick.minor.pad']  = 1.5
matplotlib.rcParams['ytick.direction']  = 'out'
matplotlib.rcParams['xtick.direction']  = 'out'
# Lines around the graph
matplotlib.rcParams['axes.spines.left']   = True
matplotlib.rcParams['axes.spines.bottom'] = True
matplotlib.rcParams['axes.spines.top']    = False
matplotlib.rcParams['axes.spines.right']  = False
# Make text editable in Adobe Illustrator
#matplotlib.rcParams['pdf.fonttype']          = 42

# Colour maps to use for the genetic diagrams
# https://personal.sron.nl/~pault/
cmap = {}
cmap['vl_purple'] = (214/255.0, 193/255.0, 222/255.0)
cmap['l_purple']  = (177/255.0, 120/255.0, 166/255.0)
cmap['purple']    = (136/255.0,  46/255.0, 114/255.0)
cmap['blue']      = ( 25/255.0, 101/255.0, 176/255.0)
cmap['l_blue']    = ( 82/255.0, 137/255.0, 199/255.0)
cmap['vl_blue']   = (123/255.0, 175/255.0, 222/255.0)
cmap['green']     = ( 78/255.0, 178/255.0, 101/255.0)
cmap['l_green']   = (144/255.0, 201/255.0, 135/255.0)
cmap['vl_green']  = (202/255.0, 224/255.0, 171/255.0)
cmap['yellow']    = (247/255.0, 238/255.0,  85/255.0)
cmap['vl_orange'] = (246/255.0, 193/255.0,  65/255.0)
cmap['l_orange']  = (241/255.0, 147/255.0,  45/255.0)
cmap['orange']    = (232/255.0,  96/255.0,  28/255.0)
cmap['red']       = (220/255.0,   5/255.0,  12/255.0)
cmap['grey']      = (119/255.0, 119/255.0, 119/255.0)
cmap['vl_grey']   = (230/255.0, 230/255.0, 230/255.0)
cmap['black']     = (0/255.0,     0/255.0,   0/255.0)

#############################################################
# HELPER FUNCTIONS
#############################################################

def fit_fn (x, y_min, y_max, k, n):
	return y_min + ((y_max - y_min)*(np.power(x, n)/(np.power(k, n) + np.power(x, n))))

def fit (x, y):
	bounds = ([0.0, 0.0, 0.0, 0.0], [0.01, 100.0, 1000.0, 1000.0])
	log_params, log_params_cov = optimize.curve_fit(fit_fn, x, y, p0=[0.0, 100.0, 100.0, 10.0], bounds=bounds, method='trf')
	return log_params

#############################################################
# PLOT THE DATA
#############################################################

# Load the data
DATA_PREFIX = './data/'
df_int2_lt = pd.read_csv(DATA_PREFIX+'Integrase2_LT.csv')
df_int2_st = pd.read_csv(DATA_PREFIX+'Integrase2_ST.csv')
df_int3_lt = pd.read_csv(DATA_PREFIX+'Integrase3_LT.csv')
df_int3_st = pd.read_csv(DATA_PREFIX+'Integrase3_ST.csv')
df_int4_st = pd.read_csv(DATA_PREFIX+'Integrase4_ST.csv')
df_int5_st = pd.read_csv(DATA_PREFIX+'Integrase5_ST.csv')
df_int7_st = pd.read_csv(DATA_PREFIX+'Integrase7_ST.csv')
df_int8_lt = pd.read_csv(DATA_PREFIX+'Integrase8_LT.csv')
df_int8_st = pd.read_csv(DATA_PREFIX+'Integrase8_ST.csv')

df_int2 = pd.concat([df_int2_lt, df_int2_st])
df_int3 = pd.concat([df_int3_lt, df_int3_st])
df_int4 = df_int4_st
df_int5 = df_int5_st
df_int7 = df_int7_st
df_int8 = pd.concat([df_int8_lt, df_int8_st])

# Fit to the data
x_data_fit = np.linspace(0, 370, 100)
params_int2 = fit(df_int2['minutes'], df_int2['flipped_percent'])
params_int3 = fit(df_int3['minutes'], df_int3['flipped_percent'])
params_int4 = fit(df_int4['minutes'], df_int4['flipped_percent'])
params_int5 = fit(df_int5['minutes'], df_int5['flipped_percent'])
params_int7 = fit(df_int7['minutes'], df_int7['flipped_percent'])
params_int8 = fit(df_int8['minutes'], df_int8['flipped_percent'])

def plot_timecourse(df, params, line_color, filename, show_y_labels=True):
	# Create the figure
	x_data = df['minutes']
	y_data = df['flipped_percent']
	markersize = 1.5
	fig = plt.figure(figsize=(2.2, 1.35))
	gs = gridspec.GridSpec(1, 1)
	ax = plt.subplot(gs[0])
	# Plot the individual data points (these should in the end nbe the averages and SDs)
	ax.plot(x_data, y_data, 'o', color=(0,0,0), markersize=markersize, zorder=-5)
	# Plot the fitted curves
	ax.plot(x_data_fit, fit_fn(x_data_fit, *params), linestyle='-', color=line_color, lw=1.6, zorder=-10)
	# Plot 100% line
	ax.plot([0,380], [100,100], color=(0,0,0), linestyle='--', lw=0.5, zorder=-1)
	# Plot 50% flipped time
	ax.plot([params[2],params[2]], [0,50], color=(0,0,0), linestyle='--', lw=0.5, zorder=-100)
	# Set the dimensions of the axes
	ax.set_xlim([0, 370])
	ax.set_ylim([-1, 105])
	ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
	ax.set_xticklabels([0, 1, 2, 3, 4, 5, 6])
	ax.set_yticks([0, 25, 50, 75, 100])
	ax.set_yticklabels([0, 25, 50, 75, 100])
	if show_y_labels == False:
		ax.set_yticklabels([])
	plt.subplots_adjust(hspace=.0 , wspace=.00, left=.15, right=.95, top=.95, bottom=.14)
	fig.savefig(filename, transparent=True)
	plt.close('all')

def plot_timecourse_with_st(df, df_st, params, line_color, filename, show_y_labels=True):
	# Create the figure
	x_data = df['minutes']
	y_data = df['flipped_percent']
	x_data_st = df_st['minutes']
	y_data_st = df_st['flipped_percent']
	markersize = 3.0
	markersize_st = 1.5
	fig = plt.figure(figsize=(2.2, 1.35))
	gs = gridspec.GridSpec(1, 1)
	ax = plt.subplot(gs[0])
	# Plot the individual data points (these should in the end nbe the averages and SDs)
	ax.plot(x_data, y_data, 'x', color=(0,0,0), markersize=markersize, lw=0.8, zorder=-5)
	ax.plot(x_data_st, y_data_st, 'o', color=(0,0,0), markersize=markersize_st, zorder=-5)
	# Plot the fitted curves
	ax.plot(x_data_fit, fit_fn(x_data_fit, *params), linestyle='-', color=line_color, lw=1.6, zorder=-10)
	# Plot 100% line
	ax.plot([0,380], [100,100], color=(0,0,0), linestyle='--', lw=0.5, zorder=-1)
	# Plot 50% flipped time
	ax.plot([params[2],params[2]], [0,50], color=(0,0,0), linestyle='--', lw=0.5, zorder=-100)
	# Set the dimensions of the axes
	ax.set_xlim([0, 370])
	ax.set_ylim([-1, 105])
	ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
	ax.set_xticklabels([0, 1, 2, 3, 4, 5, 6])
	ax.set_yticks([0, 25, 50, 75, 100])
	ax.set_yticklabels([0, 25, 50, 75, 100])
	if show_y_labels == False:
		ax.set_yticklabels([])
	plt.subplots_adjust(hspace=.0 , wspace=.00, left=.15, right=.95, top=.95, bottom=.14)
	fig.savefig(filename, transparent=True)
	plt.close('all')

plot_timecourse_with_st(df_int2_lt, df_int2_st, params_int2, cmap['purple'], './plots/int2_timecourse.pdf')
plot_timecourse_with_st(df_int3_lt, df_int3_st, params_int3, cmap['blue'], './plots/int3_timecourse.pdf', show_y_labels=False)
plot_timecourse(df_int4_st, params_int4, cmap['orange'], './plots/int4_timecourse.pdf', show_y_labels=False)
plot_timecourse(df_int5_st, params_int5, cmap['green'], './plots/int5_timecourse.pdf')
plot_timecourse(df_int7_st, params_int7, cmap['grey'], './plots/int7_timecourse.pdf', show_y_labels=False)
plot_timecourse_with_st(df_int8_lt, df_int8_st, params_int8, cmap['red'], './plots/int8_timecourse.pdf', show_y_labels=False)

#############################################################
# SAVE THE FITTED VALUES TO FILE
#############################################################

col_names = ['sample', 'y_min', 'y_max', 'k', 'n']
f = open('timecourse_fits.txt', 'w')
f.write('\t'.join(col_names) + '\n')
f.write('\t'.join(['int2'] + list([str(x) for x in params_int2])) + '\n')
f.write('\t'.join(['int3'] + list([str(x) for x in params_int3])) + '\n')
f.write('\t'.join(['int4'] + list([str(x) for x in params_int4])) + '\n')
f.write('\t'.join(['int5'] + list([str(x) for x in params_int5])) + '\n')
f.write('\t'.join(['int7'] + list([str(x) for x in params_int7])) + '\n')
f.write('\t'.join(['int8'] + list([str(x) for x in params_int8])) + '\n')
f.close()
