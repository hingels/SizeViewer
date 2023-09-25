#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 18:01:38 2023

@author: henryingels
"""

import os
import glob
import matplotlib as mpl
# resolution = 400
resolution = 200
mpl.rcParams['figure.dpi'] = resolution
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.right'] = False
# width, height = mpl.rcParams["figure.figsize"]
width, height = mpl.rcParamsDefault["figure.figsize"]
from matplotlib import pyplot as plt, cm
import pandas as pd
import numpy as np

datafolder = "/Users/henryingels/Desktop/Complete data copy/"
prefix = 'ConstantBinsTable_'
suffix = '.dat'
# names = [
#     "230701a, pre+96h",
#     "230701a, pre+96h fluor",
#     "230701a, pre+96h fluor +DiO",
#     "230701a, pre+96h fluor +DiO+Triton"
# ]
names = [
    "230709 1-100000 standard",
    "230709 1-100000 standard #2",
    "230729 dilut L1-100 contr",
    "230729 dilut L1-100 contr #2",
    "230728 KW+EW"#,
    # "230729 KW+EW re-measure",
    # "230729 KW+EW re-measure 2",
    # "230729 KW+EW re-measure 3",
    # "230729 KW+EW re-measure 4"
]

# filepaths = [os.path.join(datafolder, path) for path in [
#     "230701a, pre+96h/ConstantBinsTable_230701a, pre+96h.dat",
#     "230701a, pre+96h fluor/ConstantBinsTable_230701a, pre+96h fluor.dat",
#     "230701a, pre+96h fluor +DiO/ConstantBinsTable_230701a, pre+96h fluor +DiO.dat",
#     "230701a, pre+96h fluor +DiO+Triton/ConstantBinsTable_230701a, pre+96h fluor +DiO+Triton.dat"]]

# filepaths = glob.glob(os.path.join(datafolder, f'**/{prefix}*{suffix}'), recursive = True)
all_files = {
    filename.removeprefix(prefix).removesuffix(suffix): os.path.join(path, filename)
    for (path, subdirs, files) in os.walk(datafolder) for filename in files
    if filename.startswith(prefix) and filename.endswith(suffix)}
filepaths = [all_files[name] for name in names]

# filepaths = [os.path.join(datafolder, name, f"{prefix}{name}{suffix}") for name in names]


num_of_plots = len(filepaths)
height *= (num_of_plots/2)
height = min(np.floor(65536/resolution), height)
mpl.rcParams["figure.figsize"] = [width, height]
colors = cm.plasma(np.linspace(0, 1, num_of_plots))


fig, axs = plt.subplots(num_of_plots, 1)
fig.subplots_adjust(hspace=-0.05*height)

overall_max = 0

# axes[1].bar(bins, sizes, width = width)
for i, ax in enumerate(axs):
    data = pd.read_csv(filepaths[i], sep = '\t ', engine = 'python').iloc[:250, :]
    # bins = [0.0, *data['UpperBinDiameter_[nm]']]
    bins = data['CenterBinDiameter_[nm]']
    sizes = data['PSD_corrected_[counts/mL/nm]']
    width = bins[1] - bins[0]

    overall_max = max(sizes.max(), overall_max)
    
    plt.sca(ax)
    plt.bar(bins, sizes, width = width, color = colors[i], alpha = 0.7)
    plt.xlim(0, 250)
    ax.patch.set_alpha(0)
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    
    name = os.path.basename(filepaths[i])
    name = name.removeprefix(prefix)
    name = name.removesuffix(suffix)
    # plt.xlabel(name)
    plt.text(1.05, 0, name, fontsize=12, transform = ax.transAxes)#plt.gcf().transSubfigure)#transform=ax.figure.transFigure)#plt.gcf().transFigure)
    
    if i == (num_of_plots - 1):
        ax.yaxis.get_offset_text().set_x(-0.1)
        plt.xlabel("Diameter (nm)")
        # plt.ylabel("Particle size distribution (counts/mL/nm)", horizontalalignment = 'left', fontsize = 8)
        break
    plt.yticks([])
    plt.xticks([])

for ax in axs:
    plt.sca(ax)
    plt.ylim(0, overall_max)

plt.text(0, 0.45, "Particle size distribution (counts/mL/nm)", fontsize=12, transform = fig.transFigure, rotation = 'vertical', verticalalignment = 'center')