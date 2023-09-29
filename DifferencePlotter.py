#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 11:52:38 2023

@author: henryingels
"""

import os
import xml.etree.ElementTree as ET
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

all_dat_files = {
    filename.removeprefix(prefix).removesuffix(suffix): os.path.join(path, filename)
    for (path, subdirs, files) in os.walk(datafolder) for filename in files
    if filename.startswith(prefix) and filename.endswith(suffix)}
all_xml_files = {
    '_'.join(split_name[:-2]): os.path.join(path, filename)
    for (path, subdirs, files) in os.walk(datafolder) for filename in files
    if filename.endswith('.xml') and len(split_name := filename.split('_')) > 2}
filepaths = [all_dat_files[name] for name in names]
xml_filepaths = [all_xml_files[name] for name in names]


num_of_plots = len(filepaths)
height *= (num_of_plots/3)
height = min(np.floor(65536/resolution), height)
mpl.rcParams["figure.figsize"] = [width, height]
colors = cm.plasma(np.linspace(0, 1, num_of_plots))

grid_color = '0.8'


fig, axs = plt.subplots(num_of_plots, 1)
fig.subplots_adjust(hspace=-0.05*height)
transFigure = fig.transFigure


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
        
    if i == final_i:
        ax.yaxis.get_offset_text().set_x(-0.1)
        plt.xlabel("Diameter (nm)")
        plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))
        ax.spines.left.set_visible(True)
        break
    plt.yticks([])
    plt.xticks([])
    
    previous_sizes = sizes
    
    
    with open(xml_filepaths[i]) as xml_file:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        for setting in root.find('RecordingSettings'):
            print(setting)

origins = []
for i, ax in enumerate(axs):
    plt.sca(ax)
    plt.ylim(overall_min, overall_max)
    if i == final_i:
        plt.axhline(0, color = 'black')
    else:
        ax.spines['bottom'].set_position(('data', 0))
    origin_transDisplay = ax.transData.transform([0, 0])
    origins.append(transFigure.inverted().transform(origin_transDisplay))

final_ax = ax   # Using the ax from the final (bottom) plot:
final_ax.xaxis.set_tick_params(width = 2)
final_ax.yaxis.set_tick_params(width = 2)
tick_values, tick_labels = plt.xticks()
final_i = len(tick_values) - 1

grid_proportion_of_figure = 0.9

right_edge = None
right_edge_figure = None

for i, tick_value in enumerate(tick_values):
    display_coords = final_ax.transData.transform([tick_value, overall_min])
    figure_x, figure_y = transFigure.inverted().transform(display_coords)
    
    line = plt.Line2D([figure_x, figure_x], [figure_y, grid_proportion_of_figure], lw = 1, color = grid_color, transform = transFigure, zorder = 0)
    fig.add_artist(line)
    line.set_clip_on(False)
    
    if i == final_i:
        right_edge, _ = display_coords
        right_edge_figure = figure_x
    

plt.text(0, 0.45, "Particle size distribution (counts/mL/nm)", fontsize=12, transform = transFigure, rotation = 'vertical', verticalalignment = 'center')

plt.text(0, 0.95, "Shadows measure difference between a plot and the one above it.", fontsize=12, transform = transFigure, verticalalignment = 'center')


axis_positions = [origin[1] for origin in origins]
cell_height = axis_positions[0] - axis_positions[1]
table_top = axis_positions[0] + 0.5*cell_height
table_bottom = axis_positions[-1] - 0.5*cell_height

def generate_rows():
    for i, ax in enumerate(axs):
        row = []
        
        name = os.path.basename(filepaths[i])
        name = name.removeprefix(prefix)
        name = name.removesuffix(suffix)
        row.append(name)
        
        row.append('Hello')
        
        row.append('world')
        
        yield row

table_width = 1
margin = 0
display_coords = final_ax.transData.transform([0, overall_min])
edge = right_edge_figure + margin

column_widths = [0.6, 0.2, 0.2]
assert sum(column_widths) == 1

table = plt.table(
    tuple(generate_rows()),
    bbox = mpl.transforms.Bbox([[edge, table_bottom], [edge + table_width, table_top]]),
    transform = transFigure,
    cellLoc = 'left', colWidths = column_widths)
table.auto_set_font_size(False)
table.set_fontsize(12)
for cell in table.get_celld().values():
    cell.set(edgecolor = grid_color)
