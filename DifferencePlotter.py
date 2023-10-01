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

grid_color = '0.8'

rejected_maxima_marker = {'marker': 'o', 'fillstyle': 'none', 'color': '0.5', 'linestyle': 'none'}
maxima_marker = {'marker': 'o', 'fillstyle': 'none', 'color': 'black', 'linestyle': 'none'}


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


x_lim = 250
datafolder = "/Volumes/Lab Drive/ViewSizer 3000/Complete data"
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


class Sample():
    def __init__(self, folder):
        self.folder = folder
        info_path = None
        xml_path = None
        dat_path = None
        for (path, subdirs, files) in os.walk(folder):
            for filename in files:
                full_path = os.path.join(path, filename)
                if filename == 'info.md':
                    info_path = full_path
                split_name = filename.split('_')
                if filename.endswith('.xml') and len(split_name) > 2:
                    truncated = '_'.join(split_name[:-2])
                    if truncated.endswith('Process') is False and truncated.endswith('Temperature') is False:
                        xml_path = full_path
                if filename.startswith(prefix) and filename.endswith(suffix):
                    dat_path = full_path
        self.xml = xml_path
        self.dat = dat_path
        
        info = dict()
        if info_path is not None:
            with open(info_path) as info_file:
                for line in info_file.readlines():
                    prop, value = line.split('=')
                    value = value.strip()
                    setattr(self, prop, value)
                    info[prop] = value
        self.info = info
        
        filename = os.path.basename(folder).removeprefix(prefix).removesuffix(suffix)
        self.filename = filename
        if hasattr(self, 'name') is False:
            self.name = filename
        
        # self.siblings = []
    # def find_similar(self, samples):
    #     if hasattr(self, 'measurement') is False: return
    #     for sample in samples:
            

def generate_samples():
    for folder in os.listdir(datafolder):
        sample = Sample(os.path.join(datafolder, folder))
        if sample.name not in filenames: continue
        yield sample.filename, sample
unordered_samples = dict(generate_samples())
samples = [unordered_samples[name] for name in filenames]

# # measurement_numbers = {sample.measurement: sample if hasattr(sample, 'measurement') for sample in samples}
# measurement_numbers = dict()
# for sample in samples:
#     if hasattr(sample, 'measurement') is False: continue
#     measurement = sample.measurement
#     if measurement not in measurement_numbers:
#         measurement_numbers[measurement] = []
#     measurement_numbers[measurement].append(sample)
# for measurement_number, shared_samples in measurement_numbers.items():
#     for sample in shared_samples:
#         for other_sample in shared_samples:
            


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
previous_sizes = None
for i, ax in enumerate(axs):
    sample = samples[i]
    
    data = pd.read_csv(sample.dat, sep = '\t ', engine = 'python').iloc[:250, :]
    bins = data['CenterBinDiameter_[nm]']
    sizes = data['PSD_corrected_[counts/mL/nm]']
    width = bins[1] - bins[0]
    
    plt.sca(ax)
    plt.bar(bins, sizes, width = width, color = colors[i], alpha = 0.7, align = 'center')
    
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
    
    
    cumulative = np.cumsum(sizes)
    plt.plot(bins, cumulative/100)
        
    
    
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
    

text_shift = 0.1

plt.text(0, 0.45 + text_shift, "Particle size distribution (counts/mL/nm)", fontsize=12, transform = transFigure, rotation = 'vertical', verticalalignment = 'center')

plt.text(0, 0.95 + text_shift, "Shadows show difference between a plot and the one above it.", fontsize=12, transform = transFigure, verticalalignment = 'center')
plt.text(0, 0.93 + text_shift, filter_description, fontsize=12, transform = transFigure, verticalalignment = 'center')

icon_x = 0.01
text_x = 0.02
rejected_maxima_icon, = plt.plot([icon_x], [0.91 + text_shift], **rejected_maxima_marker, transform = transFigure)
rejected_maxima_icon.set_clip_on(False)
plt.text(text_x, 0.91 + text_shift, maxima_candidate_description, fontsize=12, transform = transFigure, verticalalignment = 'center')
maxima_icon, = plt.plot([icon_x], [0.89 + text_shift], **maxima_marker, transform = transFigure)
maxima_icon.set_clip_on(False)
plt.text(text_x, 0.89 + text_shift, maxima_description, fontsize=12, transform = transFigure, verticalalignment = 'center')
plt.text(0, 0.87 + text_shift, "Measured at room temperature.", fontsize=12, transform = transFigure, verticalalignment = 'center')


axis_positions = [origin[1] for origin in origins]
cell_height = axis_positions[0] - axis_positions[1]
table_top = axis_positions[0] + 0.5*cell_height
table_bottom = axis_positions[-1] - 0.5*cell_height



class Setting():
    def __init__(self, name, units, column = None, sample_values: dict = None, show_unit = False, show_name = False, datatype = str, depends_on = None):
        self.name = name
        self.units = units
        self.column = column
        self.sample_values = dict()
        if sample_values is not None:
            for sample, value in sample_values.items():
                self.set_value(sample, value)
        
        self.show_unit = show_unit
        self.show_name = show_name
        self.datatype = datatype
        self.depends_on = depends_on
    def set_value(self, sample: Sample, value):
        datatype = self.datatype
        self.sample_values[sample] = datatype(value) if datatype is not bool else value.lower() == "true"
    def get_value(self, sample: Sample):
        return self.sample_values[sample]
class Settings():
    def __init__(self, settings_dict):
        self.tags = settings_dict
        columns = []
        for setting in settings_dict.values():
            column = setting.column
            if column is None: continue
            if column > (len(columns) - 1):
                columns.append([])
            columns[column].append(setting)
        self.columns = columns
    def by_tag(self, tag):
        return self.tags[tag]
    def apply_dependencies(self):
        '''
        For any setting that is dependent on another setting, set it to zero (or equivalent) if the dependency has a value of False.
        '''
        for tag in self.tags:
            setting = self.by_tag(tag)
            depends_on = setting.depends_on
            if depends_on is None: continue
            for sample, value in setting.sample_values.items():
                if depends_on.get_value(sample) is False:
                    setting.set_value(sample, setting.datatype(0))


red_enabled = Setting('Red enabled', '', datatype = bool)
green_enabled = Setting('Green enabled', '', datatype = bool)
blue_enabled = Setting('Blue enabled', '', datatype = bool)
settings = Settings(OrderedDict({
    'RedLaserPower': Setting('R', 'mW', column = 0, datatype = int, show_name = True, depends_on = red_enabled),
    'RedLaserEnabled': red_enabled,
    'GreenLaserPower': Setting('G', 'mW', column = 0, datatype = int, show_name = True, depends_on = green_enabled),
    'GreenLaserEnabled': green_enabled,
    'BlueLaserPower': Setting('B', 'mW', column = 0, datatype = int, show_name = True, depends_on = blue_enabled),
    'BlueLaserEnabled': blue_enabled,
    'Exposure': Setting('Exposure', 'ms', column = 1, datatype = int),
    'Gain': Setting('Gain', 'dB', column = 2, datatype = int)}))

# row_groups = []
def generate_rows():
    for i, ax in enumerate(axs):
        row = []
        
        sample = samples[i]
        
        # name = sample.name
        # row.append(name)
        
        if hasattr(sample, 'treatment'):
            treatment1 = sample.treatment
        elif hasattr(sample, 'treatment1'):
            treatment1 = sample.treatment1
        else:
            treatment1 = None
        row.append(treatment1)
        if hasattr(sample, 'wait'):
            wait1 = sample.wait
        elif hasattr(sample, 'wait1'):
            wait1 = sample.wait1
        else:
            wait1 = None
        row.append(wait1)
        
        if hasattr(sample, 'treatment2'):
            treatment2 = sample.treatment2
        else:
            treatment2 = None
        row.append(treatment2)
        if hasattr(sample, 'wait2'):
            wait2 = sample.wait2
        else:
            wait2 = None
        row.append(wait2)
                
        row.append(sample.experimental_unit)
        
        if hasattr(sample, 'filter'):
            filter_used = sample.filter
        else:
            filter_used = None
        row.append(filter_used)
                
        with open(sample.xml) as xml_file:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            for entry in root.find('RecordingSettings'):
                tag = entry.tag
                if tag in settings.tags:
                    setting = settings.by_tag(tag)
                    setting.set_value(sample, entry.text)
        settings.apply_dependencies()
                    
        for column in settings.columns:
            content = '\n'.join(
                setting.show_name*f"{setting.name}: " + f"{setting.get_value(sample)}" + setting.show_unit*f" ({setting.units})"
                for setting in column )
            row.append(content)

        yield row

table_width = 1.5
margin = 0
display_coords = final_ax.transData.transform([0, overall_min])
edge = right_edge_figure + margin

# column_names = ["Sample", "Power\n(mW)"]
# column_names = ["Sample", "Power\n(mW)", "Filter cutoff (nm)"]
column_names = ["1st\ntreatment\n(µM)", "1st\n4°C\nwait\n(h)", "2nd\ntreatment\n(µM)", "2nd\n4°C\nwait\n(h)", "Experimental\nunit", "Filter", "Power\n(mW)", "Exposure\n(ms)", "Gain\n(dB)"]
column_widths = [0.14, 0.07, 0.14, 0.07, 0.19, 0.08, 0.1, 0.13, 0.08]
assert 1 - (width_sum := sum(column_widths)) < 0.001, f"sum(column_widths) = {width_sum} != 1"

table = plt.table(
    tuple(generate_rows()),
    bbox = mpl.transforms.Bbox([[edge, table_bottom], [edge + table_width, table_top]]),
    transform = transFigure,
    cellLoc = 'left', colWidths = column_widths)
table.auto_set_font_size(False)
table.set_fontsize(12)


for i, name in enumerate(column_names):
    new_cell = table.add_cell(-1, i, width = column_widths[i], height = 0.1, text = name, loc = 'left')
    new_cell.set_text_props(fontweight = 'bold')#, rotation = 45)


for cell in table.get_celld().values():
    cell.set(edgecolor = grid_color)

