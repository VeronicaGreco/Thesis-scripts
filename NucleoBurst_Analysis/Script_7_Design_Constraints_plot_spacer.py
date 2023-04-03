
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
df_int2_spacer = pd.read_csv(DATA_PREFIX+'Integrase2_Spacer.csv')

df_int2_10bp = df_int2_spacer[df_int2_spacer.spacer_length.eq(10)]
df_int2_500bp = df_int2_spacer[df_int2_spacer.spacer_length.eq(500)]
df_int2_1000bp = df_int2_spacer[df_int2_spacer.spacer_length.eq(1000)]

# Fit to the data
x_data_fit = np.linspace(0, 220, 100)
params_int2_10bp = fit(df_int2_10bp['minutes'], df_int2_10bp['flipped_percent'])
params_int2_500bp = fit(df_int2_500bp['minutes'], df_int2_500bp['flipped_percent'])
params_int2_1000bp = fit(df_int2_1000bp['minutes'], df_int2_1000bp['flipped_percent'])

def plot_timecourse(df_list, params_list, line_color_list, point_color_list, filename):
	# Create the figure
	markersize = 1.5
	fig = plt.figure(figsize=(2.3, 1.8))
	gs = gridspec.GridSpec(1, 1)
	ax = plt.subplot(gs[0])
	for idx in range(len(df_list)):
		x_data = df_list[idx]['minutes']
		y_data = df_list[idx]['flipped_percent']
		# Plot the individual data points (these should in the end nbe the averages and SDs)
		ax.plot(x_data, y_data, 'o', color=point_color_list[idx], markersize=markersize, zorder=-5)
		# Plot the fitted curves
		ax.plot(x_data_fit, fit_fn(x_data_fit, *params_list[idx]), linestyle='-', color=line_color_list[idx], lw=1.6, zorder=-10)
	# Plot 100% line
	ax.plot([0,380], [100,100], color=(0,0,0), linestyle='--', lw=0.5, zorder=-1)
	# Set the dimensions of the axes
	ax.set_xlim([60, 220])
	ax.set_ylim([-1, 105])
	ax.set_xticks([60, 120, 180])
	ax.set_xticklabels([1, 2, 3])
	ax.set_yticks([0, 25, 50, 75, 100])
	ax.set_yticklabels([0, 25, 50, 75, 100])
	plt.subplots_adjust(hspace=.0 , wspace=.00, left=.12, right=.95, top=.95, bottom=.12)
	fig.savefig(filename, transparent=True)
	plt.close('all')

df_list = [df_int2_10bp, df_int2_500bp, df_int2_1000bp]
params_list = [params_int2_10bp, params_int2_500bp, params_int2_1000bp]
line_color_list = [(0,0,0), (0.4,0.4,0.4), (0.7,0.7,0.7)]
point_color_list = [(0,0,0), (0.35,0.35,0.35), (0.6,0.6,0.6)]
plot_timecourse(df_list, params_list, line_color_list, point_color_list, './plots/int2_spacer.pdf')

#############################################################
# SAVE THE FITTED VALUES TO FILE
#############################################################

col_names = ['sample', 'y_min', 'y_max', 'k', 'n']
f = open('timecourse_fits_spacer.txt', 'w')
f.write('\t'.join(col_names) + '\n')
f.write('\t'.join(['int2_10bp'] + list([str(x) for x in params_int2_10bp])) + '\n')
f.write('\t'.join(['int2_500bp'] + list([str(x) for x in params_int2_500bp])) + '\n')
f.write('\t'.join(['int2_1000bp'] + list([str(x) for x in params_int2_1000bp])) + '\n')
f.close()
