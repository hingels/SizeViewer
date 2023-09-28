#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 11:52:38 2023

@author: henryingels
"""

import os
import pandas as pd
import numpy as np

import matplotlib as mpl
resolution = 200
mpl.rcParams['figure.dpi'] = resolution
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.left'] = False
mpl.rcParams['axes.spines.right'] = False
width, height = mpl.rcParamsDefault["figure.figsize"]
from matplotlib import pyplot as plt, cm

x_lim = 250
datafolder = "/Users/henryingels/Desktop/Complete data copy/"
prefix = 'ConstantBinsTable_'
suffix = '.dat'
names = [
    # "230709 1-100000 standard",
    # "230709 1-100000 standard #2",
    "230729 dilut L1-100 contr",
    "230729 dilut L1-100 contr #2",
    "230728 KW+EW",
    "230729 KW+EW re-measure",
    "230729 KW+EW re-measure 2",
    "230729 KW+EW re-measure 3",
    "230729 KW+EW re-measure 4"
]

all_files = {
    filename.removeprefix(prefix).removesuffix(suffix): os.path.join(path, filename)
    for (path, subdirs, files) in os.walk(datafolder) for filename in files
    if filename.startswith(prefix) and filename.endswith(suffix)}
filepaths = [all_files[name] for name in names]


num_of_plots = len(filepaths)
height *= (num_of_plots/3)
height = min(np.floor(65536/resolution), height)
mpl.rcParams["figure.figsize"] = [width, height]
colors = cm.plasma(np.linspace(0, 1, num_of_plots))


fig, axs = plt.subplots(num_of_plots, 1)
fig.subplots_adjust(hspace=-0.05*height)


final_i = num_of_plots - 1
overall_min, overall_max = 0, 0
previous_sizes = None
for i, ax in enumerate(axs):
    data = pd.read_csv(filepaths[i], sep = '\t ', engine = 'python').iloc[:250, :]
    bins = data['CenterBinDiameter_[nm]']
    sizes = data['PSD_corrected_[counts/mL/nm]']
    width = bins[1] - bins[0]
    
    plt.sca(ax)
    plt.bar(bins, sizes, width = width, color = colors[i], alpha = 0.7, align = 'center')
    
    overall_max = max(sizes.max(), overall_max)
    if previous_sizes is not None:
        size_differences = sizes - previous_sizes
        overall_max = max(size_differences.max(), overall_max)
        overall_min = min(size_differences.min(), overall_min)
        plt.bar(bins, size_differences, width = width, color = 'black', alpha = 0.3, align = 'center')
    
    plt.xlim(0, x_lim)
    ax.patch.set_alpha(0)
    
    name = os.path.basename(filepaths[i])
    name = name.removeprefix(prefix)
    name = name.removesuffix(suffix)
    plt.text(x_lim*1.05, 0, name, fontsize=12, transform = ax.transData)
    
    if i == final_i:
        ax.yaxis.get_offset_text().set_x(-0.1)
        plt.xlabel("Diameter (nm)")
        plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))
        ax.spines.left.set_visible(True)
        break
    plt.yticks([])
    plt.xticks([])
    
    previous_sizes = sizes

for i, ax in enumerate(axs):
    plt.sca(ax)
    plt.ylim(overall_min, overall_max)
    if i == final_i:
        plt.axhline(0, color = 'black')
    else:
        ax.spines['bottom'].set_position(('data', 0))

final_ax = ax   # Using the ax from the final (bottom) plot:
final_ax.xaxis.set_tick_params(width = 2)
final_ax.yaxis.set_tick_params(width = 2)
tick_values, tick_labels = plt.xticks()
for tick_value in tick_values:
    display_coords = final_ax.transData.transform([tick_value, overall_min])
    figure_x, figure_y = fig.transFigure.inverted().transform(display_coords)
    
    line = plt.Line2D([figure_x, figure_x], [figure_y, 0.9], lw = 2, color='black', alpha=0.1, transform = fig.transFigure)
    fig.add_artist(line)
    line.set_clip_on(False)
    

plt.text(0, 0.45, "Particle size distribution (counts/mL/nm)", fontsize=12, transform = fig.transFigure, rotation = 'vertical', verticalalignment = 'center')

plt.text(0, 0.95, "Shadows measure difference between a plot and the one above it.", fontsize=12, transform = fig.transFigure, verticalalignment = 'center')
