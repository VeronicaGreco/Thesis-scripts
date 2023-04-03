import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec

# Create the read length data for this plot
# awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' barcode19.fastq > read_length_19.txt

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
matplotlib.rcParams['pdf.fonttype']          = 42 


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
cmap['grey']      = (130/255.0, 130/255.0, 130/255.0)
cmap['vl_grey']   = (230/255.0, 230/255.0, 230/255.0)
cmap['l_grey']   =  (200/255.0, 200/255.0, 200/255.0)
cmap['d_grey']    = ( 50/255.0,  50/255.0,  50/255.0)

#############################################################
# ANALYSIS
#############################################################

num = ['16', '17', '18', '19']
#num = ['16']
col = []

for n in num:
	#f = open('./fastq/read_length_'+n+'.txt', 'r')
	f = open('/Users/fg17411/Desktop/PhD/2-Memory_ReadOut/Memory_analysis/PCR_Calibration_Cleanups_VG/November_2021_results/fastq/barcode'+n+'_lengths.txt', 'r')

	data = [0]*10000

	for l in f:
		s = l.split(' ')
		if(int(s[0]) <= 9999):
			data[int(s[0])] = int(s[1].strip())

	short_frags = sum(data[0:3000])
	long_frags = sum(data[3000:])
	print(n, short_frags, long_frags, float(long_frags)/(short_frags+long_frags))

	# Create the figure
	fig = plt.figure(figsize=(5.0, 5.0))
	gs = gridspec.GridSpec(1, 1)
	ax = plt.subplot(gs[0])

	ax.plot(data, zorder=-5)

	ax.plot([3800, 3800], [0, max(data)*1.1], linestyle="--", zorder=-10)

	# Set the dimensions of the axes
	ax.set_xlim([0, 10000])
	ax.set_ylim([0, max(data)*1.1])
	plt.yscale("log") 
	
	# Sort out the formatting of the plot (fill entire frame)
	plt.subplots_adjust(hspace=.0 , wspace=.00, left=.15, right=.99, top=.95, bottom=.1)
	fig.savefig('dist_'+n+'.pdf', transparent=True)
	plt.close('all')


