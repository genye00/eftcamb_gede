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
This python script takes care of producing the benchmark speedup results plotting.

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

# ***************************************************************************************

"""
Hard coded options
"""

x_size        = 10.0      # width of the figure, for every 10 models, in cm
y_size        = 10.0       # height of the figure in cm
main_fontsize = 10.0      # in points
width         = 0.80      # width of the bars
alpha         = 0.3       # alpha for shaded things

# fonts:
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

# spacing:
left    = 0.04
right   = 1-left
wspace  = 0.03
bottom  = 0.25
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

# ***************************************************************************************

# plot the single benchmarks:
for ind, res in enumerate(sorted_benchmark_results_keys):
    # do not acess the first element:
    if ind==0: continue
    # process and get the data:
    res_1    = sorted_benchmark_results_keys[ind-1]
    labels_1 = sorted(benchmark_results[res_1].models_rel_b.keys())
    labels_2 = sorted(benchmark_results[res].models_rel_b.keys())
    labels   = [ label for label in labels_1 if label in labels_2 ]
    x        = np.arange( len(labels) )
    y_1      = np.array([ temp[0] for temp in [ benchmark_results[res_1].models_rel_b[key] for key in labels ] ])
    y_2      = np.array([ temp[0] for temp in [ benchmark_results[res].models_rel_b[key] for key in labels ] ])
    y_err_1  = np.array([ temp[1] for temp in [ benchmark_results[res_1].models_rel_b[key] for key in labels ] ])
    y_err_2  = np.array([ temp[1] for temp in [ benchmark_results[res].models_rel_b[key] for key in labels ] ])
    y     = y_1/y_2
    y_err = np.abs(y)*np.sqrt( (y_err_1/y_1)**2 +(y_err_2/y_2)**2  )

    num   = len(x)
    eff_x_size = float(num)/10.0*x_size
    # plot:
    fig, ax = plt.subplots()
    gs = gridspec.GridSpec(1, 1)
    ax = plt.subplot(gs[0, 0])
    fig.set_size_inches( eff_x_size/2.54, y_size/2.54 )
    # plot the data:
    rect_1 = ax.bar( x+0.1, y-1.0 , width, bottom=1.0, yerr=y_err, ecolor='k', color=colormap[0], linewidth=0.8 )
    # plot the upper shaded error bar:
    rect_2 = ax.bar( x+0.1, +y_err, width, bottom=y  , color='k' , alpha=alpha      , linewidth=0   )
    # plot the lower shaded error bar:
    rect_3 = ax.bar( x+0.1, -y_err, width, bottom=y  , color='k' , alpha=alpha      , linewidth=0   )
    # plot the grid:
    ax.grid(True)
    # set the x lim:
    ax.set_xlim( [ np.amin(x), np.amax(x)+1.0 ] )
    # get human readable names:
    labels = [ label.replace( '_', ' ' ) for label in labels ]
    # set x labels:
    ax.set_xticks( np.array( xrange( len(labels) ) ) +0.5 )
    ax.set_xticklabels( labels, rotation=90, horizontalalignment='center', fontsize=0.8*main_fontsize )
    # y-axis positioning and size of the tick labels:
    fig.canvas.draw()
    ax.tick_params(axis='y', which='both', labelleft='on', labelright='on',labelsize=0.9*main_fontsize)
    ax.set_yticklabels( [item.get_text() for item in ax.get_yticklabels()], fontsize=0.9*main_fontsize )
    ax.get_yticklabels()[0].set_verticalalignment('bottom')
    ax.get_yticklabels()[ len(ax.get_yticklabels())/2 ].set_verticalalignment('bottom')
    ax.get_yticklabels()[ len(ax.get_yticklabels())/2-1 ].set_verticalalignment('top')
    ax.get_yticklabels()[-1].set_verticalalignment('top')
    # y-axis label:
    ax.set_ylabel('time old / time new', fontsize=main_fontsize)
    ax.yaxis.set_label_position("left")
    # get the size of the labels:
    fig.canvas.draw()
    # get the default renderer:
    renderer = matplotlib.backend_bases.RendererBase()
    # default dpi is 72.0:
    default_dpi = 72.0
    # get the size of the figure:
    figure_x_size = fig.get_size_inches()[0] #: in inches
    figure_y_size = fig.get_size_inches()[1] #: in inches
    # get the maximum label size:
    max_x_tick_label = [0.0,0.0]
    for xlabel in ax.get_xticklabels():
        x_dimension = xlabel.get_window_extent(renderer).width/default_dpi   #: in inches
        y_dimension = xlabel.get_window_extent(renderer).height/default_dpi  #: in inches
        if x_dimension>max_x_tick_label[0]: max_x_tick_label[0]=x_dimension
        if y_dimension>max_x_tick_label[1]: max_x_tick_label[1]=y_dimension
    bottom_size = max_x_tick_label[1]+ 2.0*matplotlib.rcParams['axes.labelpad']/default_dpi
    # appearence:
    gs.update( bottom= bottom_size/figure_y_size, top=top, left=left, right=right, wspace=wspace)
    # plot title:
    title = 'EFTCAMB speedup of version of '+benchmark_results[res].hour+':'+benchmark_results[res].minute+' of '+benchmark_results[res].day+'/'+benchmark_results[res].month+'/'+benchmark_results[res].year\
        +' w.r.t version of '+benchmark_results[res_1].hour+':'+benchmark_results[res_1].minute+' of '+benchmark_results[res_1].day+'/'+benchmark_results[res_1].month+'/'+benchmark_results[res_1].year
    ax.set_title( title, fontsize=main_fontsize, loc='center')
    # save and close figure:
    plt.savefig(out_dir+'/pdf/'+str(res)+'_speedup.pdf')
    plt.savefig(out_dir+'/png/'+str(res)+'_speedup.png')
    plt.clf()

# speedup of the last wrt the first:

# process and get the data:
res      = sorted_benchmark_results_keys[-1]
res_1    = sorted_benchmark_results_keys[0]
labels_1 = sorted(benchmark_results[res_1].models_rel_b.keys())
labels_2 = sorted(benchmark_results[res].models_rel_b.keys())
labels   = [ label for label in labels_1 if label in labels_2 ]
x        = np.arange( len(labels) )
y_1      = np.array([ temp[0] for temp in [ benchmark_results[res_1].models_rel_b[key] for key in labels ] ])
y_2      = np.array([ temp[0] for temp in [ benchmark_results[res].models_rel_b[key] for key in labels ] ])
y_err_1  = np.array([ temp[1] for temp in [ benchmark_results[res_1].models_rel_b[key] for key in labels ] ])
y_err_2  = np.array([ temp[1] for temp in [ benchmark_results[res].models_rel_b[key] for key in labels ] ])
y        = (y_1/y_2)
y_err    = np.abs(y)*np.sqrt( (y_err_1/y_1)**2 +(y_err_2/y_2)**2  )

num   = len(x)
eff_x_size = float(num)/10.0*x_size
# plot:
fig, ax = plt.subplots()
gs = gridspec.GridSpec(1, 1)
ax = plt.subplot(gs[0, 0])
fig.set_size_inches( eff_x_size/2.54, y_size/2.54 )
# plot the data:
rect_1 = ax.bar( x+0.1, y-1.0 , width, bottom=1.0, yerr=y_err, ecolor='k', color=colormap[0], linewidth=0.8 )
# plot the upper shaded error bar:
rect_2 = ax.bar( x+0.1, +y_err, width, bottom=y  , color='k' , alpha=alpha      , linewidth=0   )
# plot the lower shaded error bar:
rect_3 = ax.bar( x+0.1, -y_err, width, bottom=y  , color='k' , alpha=alpha      , linewidth=0   )
# plot the grid:
ax.grid(True)
# set the x lim:
ax.set_xlim( [ np.amin(x), np.amax(x)+1.0 ] )
# get human readable names:
labels = [ label.replace( '_', ' ' ) for label in labels ]
# set x labels:
ax.set_xticks( np.array( xrange( len(labels) ) ) +0.5 )
ax.set_xticklabels( labels, rotation=90, horizontalalignment='center', fontsize=0.8*main_fontsize )
# y-axis positioning and size of the tick labels:
fig.canvas.draw()
ax.tick_params(axis='y', which='both', labelleft='on', labelright='on',labelsize=0.9*main_fontsize)
ax.set_yticklabels( [item.get_text() for item in ax.get_yticklabels()], fontsize=0.9*main_fontsize )
ax.get_yticklabels()[0].set_verticalalignment('bottom')
ax.get_yticklabels()[ len(ax.get_yticklabels())/2 ].set_verticalalignment('bottom')
ax.get_yticklabels()[ len(ax.get_yticklabels())/2-1 ].set_verticalalignment('top')
ax.get_yticklabels()[-1].set_verticalalignment('top')
# y-axis label:
ax.set_ylabel('time old / time new', fontsize=main_fontsize)
ax.yaxis.set_label_position("left")
# get the size of the labels:
fig.canvas.draw()
# get the default renderer:
renderer = matplotlib.backend_bases.RendererBase()
# default dpi is 72.0:
default_dpi = 72.0
# get the size of the figure:
figure_x_size = fig.get_size_inches()[0] #: in inches
figure_y_size = fig.get_size_inches()[1] #: in inches
# get the maximum label size:
max_x_tick_label = [0.0,0.0]
for xlabel in ax.get_xticklabels():
    x_dimension = xlabel.get_window_extent(renderer).width/default_dpi   #: in inches
    y_dimension = xlabel.get_window_extent(renderer).height/default_dpi  #: in inches
    if x_dimension>max_x_tick_label[0]: max_x_tick_label[0]=x_dimension
    if y_dimension>max_x_tick_label[1]: max_x_tick_label[1]=y_dimension
bottom_size = max_x_tick_label[1]+ 2.0*matplotlib.rcParams['axes.labelpad']/default_dpi
# appearence:
gs.update( bottom= bottom_size/figure_y_size, top=top, left=left, right=right, wspace=wspace)
# plot title:
title = 'EFTCAMB speedup of version of '+benchmark_results[res].hour+':'+benchmark_results[res].minute+' of '+benchmark_results[res].day+'/'+benchmark_results[res].month+'/'+benchmark_results[res].year\
    +' w.r.t version of '+benchmark_results[res_1].hour+':'+benchmark_results[res_1].minute+' of '+benchmark_results[res_1].day+'/'+benchmark_results[res_1].month+'/'+benchmark_results[res_1].year
ax.set_title( title, fontsize=main_fontsize, loc='center')
# save and close figure:
plt.savefig(out_dir+'/pdf/firstlast_speedup.pdf')
plt.savefig(out_dir+'/png/firstlast_speedup.png')
plt.clf()

exit(0)
