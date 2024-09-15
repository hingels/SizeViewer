#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 18:01:38 2023

@author: henryingels
"""

import os
import matplotlib as mpl
resolution = 200
mpl.rcParams['figure.dpi'] = resolution
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.right'] = False
width, height = mpl.rcParamsDefault["figure.figsize"]
from matplotlib import pyplot as plt, cm
import pandas as pd
import numpy as np

datafolder = "/Volumes/Lab Drive/ViewSizer 3000/Complete data"
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

overall_max = 0

for i, ax in enumerate(axs):
    data = pd.read_csv(filepaths[i], sep = '\t ', engine = 'python').iloc[:250, :]
    bins = data['CenterBinDiameter_[nm]']
    sizes = data['PSD_corrected_[counts/mL/nm]']
    width = bins[1] - bins[0]

    overall_max = max(sizes.max(), overall_max)
    
    plt.sca(ax)
    plt.bar(bins, sizes, width = width, color = colors[i], alpha = 0.7, align = 'center')
    plt.xlim(0, 250)
    ax.patch.set_alpha(0)
    
    name = os.path.basename(filepaths[i])
    name = name.removeprefix(prefix)
    name = name.removesuffix(suffix)
    plt.text(1.05, 0, name, fontsize=12, transform = ax.transAxes)#plt.gcf().transSubfigure)#transform=ax.figure.transFigure)#plt.gcf().transFigure)
    
    if i == (num_of_plots - 1):
        ax.yaxis.get_offset_text().set_x(-0.1)
        plt.xlabel("Diameter (nm)")
        break
    plt.yticks([])
    plt.xticks([])

for ax in axs:
    plt.sca(ax)
    plt.ylim(0, overall_max)

plt.text(0, 0.45, "Particle size distribution (counts/mL/nm)", fontsize=12, transform = fig.transFigure, rotation = 'vertical', verticalalignment = 'center')