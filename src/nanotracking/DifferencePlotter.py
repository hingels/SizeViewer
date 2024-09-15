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
from matplotlib import pyplot as plt, cm
from .sample_class import Sample
from .settings_classes import Setting, Settings
from .InfoComparison import compare_info
from .DrawTable import draw_table

volume = 2.3E-06
x_lim = 400
prefix = 'ConstantBinsTable_'
prefix2 = 'Videos_'
suffix = '.dat'
                    
def run_program(datafolder, output_folder, filenames, cumulative_enabled = True, peaks_enabled = False, difference_enabled = False, table_drawing_enabled = True, use_filenames = True, include_treatments = False, include_experimental_unit = False, treatments_and_waits = None, results_column_names = None, column_names = None, column_widths = None, peaks_settings = None, table_settings = None, grid_color = '0.8'):
    width, height = mpl.rcParamsDefault["figure.figsize"]

    if peaks_enabled:
        maxima_marker = peaks_settings['maxima_marker']
        rejected_maxima_marker = peaks_settings['rejected_maxima_marker']
        kernel_size = peaks_settings['kernel_size']
        kernel2_size = peaks_settings['kernel2_size']
        second_derivative_threshold = peaks_settings['second_derivative_threshold']
        x = np.linspace(0, kernel_size, kernel_size)
        kernel_std = 4     # In units of bins
        kernel_center = kernel_size/2
        gaussian = np.exp(-np.power((x - kernel_center)/kernel_std, 2)/2)/(kernel_std * np.sqrt(2*np.pi))
        lowpass_filter = gaussian / gaussian.sum()
        filter_description = f"Black lines indicate Gaussian smoothing (a low-pass filter) with $\sigma = {kernel_std}$ bins and convolution kernel of size {kernel_size} bins."
        maxima_candidate_description = f": Candidate peaks after smoothing, selected using argrelextrema in SciPy {scipy.__version__}."
        maxima_description = f": Peaks with under {second_derivative_threshold} counts/mL/nm$^3$ second derivative, computed after smoothing again with simple moving average of size {kernel2_size} bins."

    def generate_samples():
        for folder in os.listdir(datafolder):
            sample = Sample(os.path.join(datafolder, folder), prefix, suffix, videos_file_prefix = prefix2)
            if sample.filename not in filenames: continue
            yield sample.filename, sample
    if use_filenames:
        unordered_samples = dict(generate_samples())
        samples = [unordered_samples[name] for name in filenames]
    else:
        _, samples = tuple(zip(*generate_samples()))


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
    sums = []
    bins = None
    previous_sizes = None
    for i, ax in enumerate(axs):
        sample = samples[i]
        
        full_data = pd.read_csv(sample.dat, sep = '\t ', engine = 'python')
        data = full_data.iloc[:400, :]
        new_bins = data['/LowerBinDiameter_[nm]']
        
        top_nm = max(data['UpperBinDiameter_[nm]'])
        if top_nm.is_integer():
            top_nm = int(top_nm)
        
        if bins is not None:
            assert np.all(new_bins == bins) == True, 'Unequal sequence of bins between samples detected!'
        bins = new_bins
        sizes = data['PSD_corrected_[counts/mL/nm]']
        width = bins[1] - bins[0]
        bin_centers = bins + width/2
        
        plt.sca(ax)
        plt.bar(bins, sizes, width = width, color = colors[i], alpha = 0.7, align = 'edge')
        
        if peaks_enabled:
            filtered = np.convolve(sizes, lowpass_filter, mode = 'same')
            # plt.plot(bin_centers, filtered, '-', color = 'black', linewidth = 0.5)
            maxima_candidates, = argrelextrema(filtered, np.greater)
            
            twice_filtered = np.convolve(filtered, [1/kernel2_size]*kernel2_size, mode = 'same')
            # plt.plot(bin_centers, twice_filtered, linewidth = 0.5)
            derivative = np.gradient(twice_filtered, bin_centers)
            second_derivative = np.gradient(derivative, bin_centers)
            # plt.plot(bin_centers, second_derivative*100, linewidth = 0.5)
            second_deriv_negative, = np.where(second_derivative < second_derivative_threshold)
        
            maxima = np.array([index for index in maxima_candidates if index in second_deriv_negative])
            assert len(maxima) != 0, 'No peaks found. The second derivative threshold may be too high.'
            
            rejected_candidates = np.array([entry for entry in maxima_candidates if entry not in maxima])
            if len(rejected_candidates) != 0:
                plt.plot(bin_centers[rejected_candidates], filtered[rejected_candidates], **rejected_maxima_marker)
            plt.plot(bin_centers[maxima], filtered[maxima], **maxima_marker)
        
        
        if cumulative_enabled:
            cumulative_sums.append(np.cumsum(sizes)*width)
        
        if table_drawing_enabled:
            sums.append((
                ('All data', np.sum(full_data['PSD_corrected_[counts/mL/nm]'])*width),
                (top_nm, np.sum(sizes*width))
            ))        
        
        overall_max = max(sizes.max(), overall_max)
        if previous_sizes is not None and difference_enabled:
            size_differences = sizes - previous_sizes
            overall_max = max(size_differences.max(), overall_max)
            overall_min = min(size_differences.min(), overall_min)
            # plt.bar(bins, size_differences, width = width, color = 'black', alpha = 0.3, align = 'center')
            plt.bar(bins, size_differences, width = width, color = 'black', alpha = 0.3, align = 'edge')
        
        videos = sample.videos
        all_histograms = np.array([np.histogram(video, bins = bins)[0] for video in videos])
        avg_histogram = np.average(all_histograms, axis = 0)
        total_std = np.std(all_histograms, axis = 0, ddof = 1)
        scale_factor = np.array([sizes[j]/avg if (avg := avg_histogram[j]) != 0 else 0 for j in range(len(sizes)-1)])
        
        error_resizing = 0.1
        errorbars = np.array(list(zip((scale_factor*total_std)*error_resizing, [0]*len(total_std)))).T
        plt.errorbar(bin_centers[:-1], scale_factor*avg_histogram, yerr = errorbars, elinewidth = 1, linestyle = '', marker = '.', ms = 1, alpha = 0.5, color = 'black')            
        
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

    right_edge_figure = None

    for i, tick_value in enumerate(tick_values):
        display_coords = final_ax.transData.transform([tick_value, overall_min])
        figure_x, figure_y = transFigure.inverted().transform(display_coords)
        
        line = plt.Line2D([figure_x, figure_x], [figure_y, grid_proportion_of_figure], lw = 1, color = grid_color, transform = transFigure, zorder = 0)
        fig.add_artist(line)
        line.set_clip_on(False)
        
        if i == final_i:
            right_edge_figure = figure_x



    plt.text(0, 0.45, "Particle size distribution (counts/mL/nm)", fontsize=12, transform = transFigure, rotation = 'vertical', verticalalignment = 'center')

    text_shift = 0.05
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


    previous_setting = Setting('previous', name = 'Previous')
    md_settings = [
        Setting('experimental_unit', name = 'Experimental unit'),#, column = 0),
        Setting('treatment', name = 'Treatment', units = 'µM'),
        Setting('wait', name = 'Wait', units = 'h'),
        Setting('filter', name = 'Filter cut-on', units = 'nm', column = 1),
        previous_setting ]
    red_enabled = Setting('RedLaserEnabled', name = 'Red enabled', datatype = bool)
    green_enabled = Setting('GreenLaserEnabled', name = 'Green enabled', datatype = bool)
    blue_enabled = Setting('BlueLaserEnabled', name = 'Blue enabled', datatype = bool)
    detection_threshold_setting = Setting('DetectionThresholdType', name = 'Detection mode', dependencies_require = 'Manual')
    xml_settings = [
        Setting('RedLaserPower', short_name = '635nm', name = '635nm power', units = 'mW', column = 2, datatype = int, show_name = True, depends_on = red_enabled),
        red_enabled,
        Setting('GreenLaserPower', short_name = '520nm', name = '520nm power', units = 'mW', column = 2, datatype = int, show_name = True, depends_on = green_enabled),
        green_enabled,
        Setting('BlueLaserPower', short_name = '445nm', name = '445nm power', units = 'mW', column = 2, datatype = int, show_name = True, depends_on = blue_enabled),
        blue_enabled,
        Setting('Exposure', units = 'ms', datatype = int),
        Setting('Gain', units = 'dB', datatype = int),
        Setting('MeasurementStartDateTime'),#, column = 5),
        Setting('FrameRate', name = 'Framerate', units = 'fps', datatype = int),
        Setting('FramesPerVideo', name = 'Frames per video', units = 'frames', datatype = int),
        Setting('NumOfVideos', name = 'Number of videos', datatype = int),
        Setting('StirrerSpeed', name = 'Stirring speed', units = 'rpm', datatype = int),
        Setting('StirredTime', name = 'Stirred time', units = 's', datatype = int),
        detection_threshold_setting,
        Setting('DetectionThreshold', name = 'Detection threshold', datatype = float, depends_on = detection_threshold_setting) ]
    settings_list = [*md_settings, *xml_settings]
    settings = Settings(OrderedDict({setting.tag: setting for setting in settings_list}))

    results_for_csv = Setting('_results')
    compare_info(settings, samples, results_for_csv, output_folder)

    if table_drawing_enabled:
        axis_positions = [origin[1] for origin in origins]
        cell_height = axis_positions[0] - axis_positions[1]
        table_top = axis_positions[0] + 0.5*cell_height
        table_bottom = axis_positions[-1] - 0.5*cell_height

        results_for_table = results_for_csv
        
        edges = {'right': right_edge_figure, 'bottom': table_bottom, 'top': table_top}
        draw_table(fig, ax, settings, previous_setting, samples, unordered_samples, sums, treatments_and_waits, results_for_table, results_column_names, column_names, column_widths, edges, table_settings, grid_color, include_treatments, include_experimental_unit)

    fig.savefig(f"{output_folder}/Ridgeline plot.png", dpi = 300, bbox_inches='tight')