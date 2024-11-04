#----------------------------------------------------------------------------------------
#
# This file is part of EFTCAMB.
#
# Copyright (C) 2013-2019 by the EFTCAMB authors
#
# The EFTCAMB code is free software;
# You can use it, redistribute it, and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation;
# either version 3 of the License, or (at your option) any later version.
# The full text of the license can be found in the file eftcamb/LICENSE at
# the top level of the EFTCAMB distribution.
#
#----------------------------------------------------------------------------------------

"""
This python script takes care of producing the speed history plotting.

Developed by Marco Raveri (mraveri@uchicago.edu) for the EFTCAMB code.
"""

# ***************************************************************************************

__version__ = '0.0' #: version of the application

# ***************************************************************************************

"""
import of modules
"""

from import_benchmark import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines    as mlines

# ***************************************************************************************

"""
Hard coded options
"""

x_size        = 1.6180339*10.0      # width of the figure, for every 10 models, in cm
y_size        = 10.0       # height of the figure in cm
main_fontsize = 10.0      # in points
width         = 0.80      # width of the bars
alpha         = 0.3       # alpha for shaded things

# fonts:
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

# spacing:
left    = 0.1
right   = 0.99
wspace  = 0.03
bottom  = 0.23
top     = 0.95

# colors:
colormap = { 0: (203.0/255.0, 15.0/255.0, 40.0/255.0),
                 1: (255.0/255.0, 165.0/255.0, 0.0),
                 2: (42.0/255.0, 46.0/255.0, 139.0/255.0),
                 3: (0.0/255.0, 153.0/255.0, 204.0/255.0),
                 4: (0.0/255.0, 221.0/255.0, 52.0/255.0),
                 5: (0.0, 0.0, 0.0),
                 6: (0.0, 0.75, 0.75),
               }

cont_colormap = 'GnBu'

# ***************************************************************************************

num   = len(sorted_benchmark_results_keys)
eff_x_size = float(max(num,10))/10.0*x_size
# plot:
fig, ax = plt.subplots()
gs = gridspec.GridSpec(1, 1)
ax = plt.subplot(gs[0, 0])
fig.set_size_inches( eff_x_size/2.54, y_size/2.54 )

# get the x values:
x = range(num)

# plot the GR data:
y     = []
y_err = []
for res in sorted_benchmark_results_keys:
    y.append(benchmark_results[res].models_rel_b['1_EFT_GR'][0])
    y_err.append(benchmark_results[res].models_rel_b['1_EFT_GR'][1])

y = np.array(y)
y_err = np.array(y_err)

ax.plot( x, y, color='k', linewidth=1, label='GR' )
ax.plot( x, y+y_err, color='k', linewidth=0.3 )
ax.plot( x, y-y_err, color='k', linewidth=0.3 )
ax.fill_between( x, y-y_err, y+y_err, alpha=0.5, facecolor='k', linewidth=0.0)

# get the number of data lines:
done_keys    = ['1_EFT_GR','base_params']
for res in sorted_benchmark_results_keys:
    keys      = sorted( benchmark_results[res].models_rel_b.keys() )
    plot_keys = [ k for k in keys if k not in done_keys ]
    for k in plot_keys: done_keys.append(k)
number_lines = len(done_keys)

# plot the lines:
done_keys    = ['1_EFT_GR','base_params']
col_ind      = 0
for res in sorted_benchmark_results_keys:
    keys      = sorted( benchmark_results[res].models_rel_b.keys() )
    plot_keys = [ k for k in keys if k not in done_keys ]

    for key in plot_keys:
        x     = []
        y     = []
        y_err = []
        for ind, res_2 in enumerate(sorted_benchmark_results_keys):
            try:
                y.append(benchmark_results[res_2].models_rel_b[key][0])
                y_err.append(benchmark_results[res_2].models_rel_b[key][1])
                x.append(ind)
            except: pass
        x = np.array(x)
        y = np.array(y)
        y_err = np.array(y_err)

        col = plt.get_cmap(cont_colormap)( float(col_ind)/float(number_lines-1) )
        col_ind=col_ind+1
        ax.plot( x, y, color=col, linewidth=0.1 )
        ax.fill_between( x, y-y_err, y+y_err, alpha=0.1, facecolor=col, linewidth=0.0)

    for k in plot_keys: done_keys.append(k)

# set the x lim:
ax.set_xlim( [ np.amin(range(num)), np.amax(range(num)) ] )
# get human readable names:
labels = [ ]
for res in sorted_benchmark_results_keys:
    labels.append( benchmark_results[res].day+'/'+benchmark_results[res].month+'/'+benchmark_results[res].year )
# set x labels:
ax.set_xticks( range(num) )
ax.set_xticklabels( labels, rotation=90, horizontalalignment='center', fontsize=0.8*main_fontsize )
ax.get_xticklabels()[0].set_horizontalalignment('left')
ax.get_xticklabels()[-1].set_horizontalalignment('right')
# axes label:
ax.set_xlabel('day', fontsize=main_fontsize)
ax.set_ylabel('EFTCAMB time / CAMB time', fontsize=main_fontsize)
# set the y labels:
fig.canvas.draw()
ax.tick_params(axis='y', which='both', labelleft='on', labelright='on',labelsize=0.9*main_fontsize)
ax.set_yticklabels( [item.get_text() for item in ax.get_yticklabels()], fontsize=0.9*main_fontsize )
ax.get_yticklabels()[0].set_verticalalignment('bottom')
ax.get_yticklabels()[ len(ax.get_yticklabels())/2 ].set_verticalalignment('bottom')
ax.get_yticklabels()[ len(ax.get_yticklabels())/2-1 ].set_verticalalignment('top')
ax.get_yticklabels()[-1].set_verticalalignment('top')
# legend:
leg_handlers = [ mlines.Line2D([],[],color='k'), mlines.Line2D([],[],color=plt.get_cmap(cont_colormap)(1.0)) ]
names_legend = ['GR', 'other models']

plt.legend( handles=leg_handlers, labels=names_legend, fontsize=0.8*main_fontsize, ncol=2, loc='upper left',borderaxespad=0.0)

# appearence:
gs.update( bottom= bottom, top=top, left=left, right=right, wspace=wspace)
# plot title:
title = 'EFTCAMB performance history'
ax.set_title( title, fontsize=main_fontsize, loc='center')
# save and close figure:
plt.savefig(out_dir+'/pdf/speed_history.pdf')
plt.savefig(out_dir+'/png/speed_history.png')
plt.clf()

exit(0)
