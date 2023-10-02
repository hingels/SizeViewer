#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 11:52:38 2023

@author: henryingels
"""

import os
import pandas as pd
import numpy as np
import scipy
from scipy.signal import argrelextrema
from collections import OrderedDict
import matplotlib as mpl
resolution = 200
mpl.rcParams['figure.dpi'] = resolution
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.left'] = False
mpl.rcParams['axes.spines.right'] = False
width, height = mpl.rcParamsDefault["figure.figsize"]
from matplotlib import pyplot as plt, cm
from sample_class import Sample
from settings_classes import Setting, Settings


cumulative_enabled = False
peaks_enabled = True
difference_enabled = True

x_lim = 250
datafolder = "/Volumes/Lab Drive/ViewSizer 3000/Data"
prefix = 'ConstantBinsTable_'
suffix = '.dat'
filenames = [
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
table_width = 1.2
table_left_margin = 0 if difference_enabled else 0.02
minimum_table_right_margin = 0.03
column_names = ["1st\ntreatment\n(µM)", "1st\n4°C\nwait\n(h)", "2nd\ntreatment\n(µM)", "2nd\n4°C\nwait\n(h)", "Experimental\nunit", "Filter", "Power\n(mW)", "Exposure\n(ms)", "Gain\n(dB)", "Video\nduration (s)\nx quantity"]
column_widths = [0.14, 0.07, 0.14, 0.07, 0.19, 0.08, 0.1, 0.13, 0.08, 0.16]

kernel_size = 30
x = np.linspace(0, kernel_size, kernel_size)
kernel_std = 4     # In units of bins
kernel_center = kernel_size/2
gaussian = np.exp(-np.power((x - kernel_center)/kernel_std, 2)/2)/(kernel_std * np.sqrt(2*np.pi))
lowpass_filter = gaussian / gaussian.sum()
filter_description = f"Black lines indicate Gaussian smoothing (a low-pass filter) with $\sigma = {kernel_std}$ bins and convolution kernel of size {kernel_size} bins."

kernel2_size = 20
second_derivative_threshold = -30
maxima_candidate_description = f": Candidate peaks after smoothing, selected using argrelextrema in SciPy {scipy.__version__}."
maxima_description = f": Peaks with under {second_derivative_threshold} counts/mL/nm$^3$ second derivative, computed after smoothing again with simple moving average of size {kernel2_size} bins."

grid_color = '0.8'
rejected_maxima_marker = {'marker': 'o', 'fillstyle': 'none', 'color': '0.5', 'linestyle': 'none'}
maxima_marker = {'marker': 'o', 'fillstyle': 'none', 'color': 'black', 'linestyle': 'none'}


width_sum = sum(column_widths)
table_right_margin = table_width - width_sum
assert table_right_margin >= minimum_table_right_margin, f"table_right_margin = {table_right_margin} < minimum_table_right_margin = {minimum_table_right_margin}."
column_widths = np.append(column_widths, table_right_margin)
column_names.append("")

                    
def generate_samples():
    for folder in os.listdir(datafolder):
        sample = Sample(os.path.join(datafolder, folder), prefix, suffix)
        if sample.name not in filenames: continue
        yield sample.filename, sample
unordered_samples = dict(generate_samples())
samples = [unordered_samples[name] for name in filenames]


num_of_plots = len(samples)
height *= (num_of_plots/3)
height = min(np.floor(65536/resolution), height)
mpl.rcParams["figure.figsize"] = [width, height]
colors = cm.plasma(np.linspace(0, 1, num_of_plots))


fig, axs = plt.subplots(num_of_plots, 1)
fig.subplots_adjust(hspace=-0.05*height)
transFigure = fig.transFigure


final_i = num_of_plots - 1
overall_min, overall_max = 0, 0
cumulative_sums = []
bins = None
previous_sizes = None
for i, ax in enumerate(axs):
    sample = samples[i]
    
    data = pd.read_csv(sample.dat, sep = '\t ', engine = 'python').iloc[:250, :]
    new_bins = data['CenterBinDiameter_[nm]']
    if bins is not None:
        assert np.all(new_bins == bins) == True, 'Unequal sequence of bins between samples detected!'
    bins = new_bins
    sizes = data['PSD_corrected_[counts/mL/nm]']
    width = bins[1] - bins[0]
    
    plt.sca(ax)
    plt.bar(bins, sizes, width = width, color = colors[i], alpha = 0.7, align = 'center')
    
    if peaks_enabled:
        filtered = np.convolve(sizes, lowpass_filter, mode = 'same')
        plt.plot(bins, filtered, '-', color = 'black', linewidth = 0.5)
        maxima_candidates, = argrelextrema(filtered, np.greater)
        
        twice_filtered = np.convolve(filtered, [1/kernel2_size]*kernel2_size, mode = 'same')
        # plt.plot(bins, twice_filtered, linewidth = 0.5)
        derivative = np.gradient(twice_filtered, bins)
        second_derivative = np.gradient(derivative, bins)
        # plt.plot(bins, second_derivative*100, linewidth = 0.5)
        second_deriv_negative, = np.where(second_derivative < second_derivative_threshold)
    
        maxima = np.array([index for index in maxima_candidates if index in second_deriv_negative])
        assert len(maxima) != 0, 'No peaks found. The second derivative threshold may be too high.'
        
        rejected_candidates = np.array([entry for entry in maxima_candidates if entry not in maxima])
        plt.plot(bins[rejected_candidates], filtered[rejected_candidates], **rejected_maxima_marker)
        plt.plot(bins[maxima], filtered[maxima], **maxima_marker)
    
    
    if cumulative_enabled:
        cumulative_sums.append(np.cumsum(sizes))    
        
    
    overall_max = max(sizes.max(), overall_max)
    if previous_sizes is not None and difference_enabled:
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
    

if cumulative_enabled:
    max_of_cumulative_sums = 0
    for cumulative_sum in cumulative_sums:
        this_max = max(cumulative_sum)
        max_of_cumulative_sums = max(max_of_cumulative_sums, this_max)
    cumulative_sum_scaling = overall_max / max_of_cumulative_sums


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
    
    if cumulative_enabled:
        cumulative_sum = cumulative_sums[i]
        plt.plot(bins, cumulative_sum*cumulative_sum_scaling, color = 'red', linewidth = 0.5)

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


text_shift = 0.05

plt.text(0, 0.45, "Particle size distribution (counts/mL/nm)", fontsize=12, transform = transFigure, rotation = 'vertical', verticalalignment = 'center')


text_y = 0 + text_shift

if difference_enabled:
    plt.text(0, text_y, "Shadows show difference between a plot and the one above it.", fontsize=12, transform = transFigure, verticalalignment = 'center')

if peaks_enabled:
    text_y -= 0.02
    plt.text(0, text_y, filter_description, fontsize=12, transform = transFigure, verticalalignment = 'center')

if cumulative_enabled:
    text_y -= 0.02
    plt.text(0, text_y, f"Red lines are cumulative sums of unsmoothed data, scaled by {cumulative_sum_scaling:.3}.", fontsize=12, transform = transFigure, verticalalignment = 'center')


icon_x = 0.01
text_x = 0.02

if peaks_enabled:
    text_y -= 0.02
    rejected_maxima_icon, = plt.plot([icon_x], [text_y], **rejected_maxima_marker, transform = transFigure)
    rejected_maxima_icon.set_clip_on(False)
    plt.text(text_x, text_y, maxima_candidate_description, fontsize=12, transform = transFigure, verticalalignment = 'center')
    
    text_y -= 0.02
    maxima_icon, = plt.plot([icon_x], [text_y], **maxima_marker, transform = transFigure)
    maxima_icon.set_clip_on(False)
    plt.text(text_x, text_y, maxima_description, fontsize=12, transform = transFigure, verticalalignment = 'center')

text_y -= 0.02
plt.text(0, text_y, "Measured at room temperature.", fontsize=12, transform = transFigure, verticalalignment = 'center')

text_y -= 0.04
plt.text(0, text_y, " ", fontsize=12, transform = transFigure, verticalalignment = 'center')


axis_positions = [origin[1] for origin in origins]
cell_height = axis_positions[0] - axis_positions[1]
table_top = axis_positions[0] + 0.5*cell_height
table_bottom = axis_positions[-1] - 0.5*cell_height


red_enabled = Setting('RedLaserEnabled', name = 'Red enabled', datatype = bool)
green_enabled = Setting('GreenLaserEnabled', name = 'Green enabled', datatype = bool)
blue_enabled = Setting('BlueLaserEnabled', name = 'Blue enabled', datatype = bool)
xml_settings = [
    Setting('RedLaserPower', name = 'R', units = 'mW', column = 2, datatype = int, show_name = True, depends_on = red_enabled),
    red_enabled,
    Setting('GreenLaserPower', name = 'G', units = 'mW', column = 2, datatype = int, show_name = True, depends_on = green_enabled),
    green_enabled,
    Setting('BlueLaserPower', name = 'B', units = 'mW', column = 2, datatype = int, show_name = True, depends_on = blue_enabled),
    blue_enabled,
    Setting('Exposure', units = 'ms', column = 3, datatype = int),
    Setting('Gain', units = 'dB', column = 4, datatype = int),
    Setting('FrameRate', name = 'Framerate', units = 'frames/sec', datatype = int),
    Setting('FramesPerVideo', name = 'Frames per video', units = 'frames', datatype = int),
    Setting('NumOfVideos', name = 'Number of videos', datatype = int),
    Setting('StirrerSpeed', name = 'Stirring speed', datatype = int),
    Setting('StirredTime', name = 'Stirred time', datatype = int) ]
md_settings = [
    Setting('experimental_unit', column = 0),
    Setting('treatment'),
    Setting('wait'),
    Setting('filter', column = 1) ]
settings_list = [*xml_settings, *md_settings]
settings = Settings(OrderedDict({setting.tag: setting for setting in settings_list}))

def generate_rows():
    column_quantities = dict()
    def number_of_values(tag, sample):
        if (setting := settings.by_tag(tag)) is None or (value := setting.get_value(sample)) is None:
            return 0
        if type(value) is list:
            return len(value)
        return 1
    def get_multivalued(tag, sample):
        if (setting := settings.by_tag(tag)) is None or (value := setting.get_value(sample)) is None:
            return [None]
        if type(value) is list:
            value.sort()
            value = [subvalue for index, subvalue in value]
        else:
            value = [value]
        return value
    for i in range(num_of_plots):
        sample = samples[i]
        
        settings.read_files(sample)
        
        for tag in ('treatment', 'wait'):
            quantity = number_of_values(tag, sample)
            if tag not in column_quantities:
                column_quantities[tag] = quantity
                continue
            column_quantities[tag] = max(column_quantities[tag], quantity)
        
    for i, ax in enumerate(axs):
        row = []
        
        sample = samples[i]
                
        treatments = get_multivalued('treatment', sample)
        waits = get_multivalued('wait', sample)
        for j in range( max(column_quantities['treatment'], column_quantities['wait']) ):
            if j < len(treatments): row.append(treatments[j])
            else: row.append(None)
            if j < len(waits): row.append(waits[j])
            else: row.append(None)
                    
        columns = list(settings.columns.items())
        columns.sort()
        for j, column in columns:
            content = '\n'.join(
                setting.show_name*f"{setting.name}: " + f"{setting.get_value(sample)}" + setting.show_unit*f" ({setting.units})"
                for setting in column if setting.get_value(sample) is not None )
            row.append(content)
        
        framerate = settings.by_tag('FrameRate').get_value(sample)
        frames_per_video = settings.by_tag('FramesPerVideo').get_value(sample)
        video_duration = frames_per_video / framerate
        if video_duration.is_integer():
            video_duration = int(video_duration)        
        num_of_videos = settings.by_tag('NumOfVideos').get_value(sample)
        row.append(f"{video_duration}x{num_of_videos}")
        
        row.append("")
        
        yield row


display_coords = final_ax.transData.transform([0, overall_min])
edge = right_edge_figure + table_left_margin

table = plt.table(
    tuple(generate_rows()),
    bbox = mpl.transforms.Bbox([[edge, table_bottom], [edge + table_width, table_top]]),
    transform = transFigure,
    cellLoc = 'left', colWidths = column_widths)
table.auto_set_font_size(False)
table.set_fontsize(12)
fig.add_artist(table)


for i, name in enumerate(column_names):
    new_cell = table.add_cell(-1, i, width = column_widths[i], height = 0.1, text = name, loc = 'left')
    new_cell.set_text_props(fontweight = 'bold')
final_column = len(column_widths) - 1
for (row, column), cell in table.get_celld().items():
    if column == final_column:
        cell.set(edgecolor = None)
        continue
    cell.set(edgecolor = grid_color)
