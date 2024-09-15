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
num_data_points = 400
prefix = 'ConstantBinsTable_'
prefix2 = 'Videos_'
suffix = '.dat'
grid_proportion_of_figure = 0.9
text_shift = 0.05
                    
class NTA():
    # def __init__(self, datafolder, output_folder, filenames, tmp_filepath = 'tmp.bin'):
    def __init__(self, datafolder, output_folder, filenames):
        self.datafolder, self.output_folder, self.filenames = datafolder, output_folder, filenames
        os.makedirs(output_folder, exist_ok = True)
        self.table_settings, self.peak_settings, self.cumulative_enabled, self.difference_enabled = None, None, False, False
        def generate_samples():
            for folder in os.listdir(datafolder):
                sample = Sample(os.path.join(datafolder, folder), prefix, suffix, videos_file_prefix = prefix2)
                if sample.filename not in filenames: continue
                yield sample.filename, sample
        unordered_samples = dict(generate_samples())
        samples = [unordered_samples[name] for name in filenames]
        num_of_plots = len(samples)
        width, height = mpl.rcParamsDefault["figure.figsize"]
        height *= (num_of_plots/3)
        height = min(np.floor(65536/resolution), height)
        self.figsize = (width, height)
        self.colors = cm.plasma(np.linspace(0, 1, num_of_plots))
        # self.unordered_samples, self.samples, self.num_of_plots, self.tmp_filepath = unordered_samples, samples, num_of_plots, tmp_filepath
        self.unordered_samples, self.samples, self.num_of_plots = unordered_samples, samples, num_of_plots
        self.overall_min, self.overall_max = None, None
        self.maxima, self.rejected_maxima = None, None
        self.tmp_filenames = {
            'bins': os.path.join(output_folder, 'bins'),
            'sizes': os.path.join(output_folder, 'sizes'),
            'filtered_sizes': os.path.join(output_folder, 'filtered_sizes'),
            # 'maxima': os.path.join(output_folder, 'maxima'),
            # 'rejected_maxima': os.path.join(output_folder, 'rejected_maxima'),
            'top_nm': os.path.join(output_folder, 'top_nm'),
            'fulldata_size_sums': os.path.join(output_folder, 'fulldata_size_sums'),
            'cumulative_sums': os.path.join(output_folder, 'cumulative_sums'),
            'cumsum_maxima': os.path.join(output_folder, 'cumsum_maxima')
        }
    def enable_table(self, include_experimental_unit, treatments_and_waits, results_column_names, column_names, column_widths, width, margin_minimum_right, margin_left):
        table_settings = locals(); table_settings.pop('self')
        self.table_settings = table_settings
    def disable_table(self):
        self.table_settings = None
    def enable_peak_detection(self, kernel_size, kernel2_size, kernel_std_in_bins, second_derivative_threshold, maxima_marker, rejected_maxima_marker):
        peak_settings = locals(); peak_settings.pop('self')
        x = np.linspace(0, kernel_size, kernel_size)
        gaussian = np.exp(-np.power((x - kernel_size/2)/kernel_std_in_bins, 2)/2)/(kernel_std_in_bins * np.sqrt(2*np.pi))
        peak_settings['lowpass_filter'] = gaussian / gaussian.sum()
        peak_settings['filter_description'] = f"Black lines indicate Gaussian smoothing (a low-pass filter) with $\sigma = {kernel_std_in_bins}$ bins and convolution kernel of size {kernel_size} bins."
        peak_settings['maxima_candidate_description'] = f": Candidate peaks after smoothing, selected using argrelextrema in SciPy {scipy.__version__}."
        peak_settings['maxima_description'] = f": Peaks with under {second_derivative_threshold} counts/mL/nm$^3$ second derivative, computed after smoothing again with simple moving average of size {kernel2_size} bins."
        self.peak_settings = peak_settings
    def disable_peak_detection(self):
        self.peak_settings = None
    def enable_cumulative(self):
        self.cumulative_enabled = True
    def disable_cumulative(self):
        self.cumulative_enabled = False
    def enable_difference(self):
        self.difference_enabled = True
    def disable_difference(self):
        self.difference_enabled = False
    def compute(self):
        peak_settings = self.peak_settings
        peaks_enabled = (peak_settings is not None)
        if peaks_enabled:
            lowpass_filter, kernel2_size, second_derivative_threshold = peak_settings['lowpass_filter'], peak_settings['kernel2_size'], peak_settings['second_derivative_threshold']
            all_filtered, all_maxima, all_rejected, all_top_nm, fulldata_size_sums = [], [], [], [], []
        cumulative_enabled = self.cumulative_enabled
        if cumulative_enabled:
            cumulative_sums = []
            cumsum_maxima = []
        difference_enabled = self.difference_enabled
        if difference_enabled:
            all_size_differences = []

        overall_min, overall_max = 0, 0
        previous_sizes = None
        bins = None
        all_bins, all_sizes = [], []
        for sample in self.samples:
            full_data = pd.read_csv(sample.dat, sep = '\t ', engine = 'python')
            data = full_data.iloc[:num_data_points, :]
            new_bins = data['/LowerBinDiameter_[nm]']
            
            top_nm = max(data['UpperBinDiameter_[nm]'])
            if top_nm.is_integer():
                top_nm = int(top_nm)
            
            if bins is not None:
                assert np.all(new_bins == bins) == True, 'Unequal sequence of bins between samples detected!'
            bins = new_bins
            sizes = data['PSD_corrected_[counts/mL/nm]']
            width = bins[1] - bins[0]
            
            all_bins.append(bins)
            all_sizes.append(sizes)
            if peaks_enabled:
                bin_centers = bins + width/2
                filtered = np.convolve(sizes, lowpass_filter, mode = 'same')
                maxima_candidates, = argrelextrema(filtered, np.greater)
                twice_filtered = np.convolve(filtered, [1/kernel2_size]*kernel2_size, mode = 'same')
                derivative = np.gradient(twice_filtered, bin_centers)
                second_derivative = np.gradient(derivative, bin_centers)
                second_deriv_negative, = np.where(second_derivative < second_derivative_threshold)
                maxima = np.array([index for index in maxima_candidates if index in second_deriv_negative])
                assert len(maxima) != 0, 'No peaks found. The second derivative threshold may be too high.'
                rejected_candidates = np.array([entry for entry in maxima_candidates if entry not in maxima])
                all_filtered.append(filtered)
                all_maxima.append(maxima)
                all_rejected.append(rejected_candidates)
                all_top_nm.append(top_nm)
                fulldata_size_sums.append(np.sum(full_data['PSD_corrected_[counts/mL/nm]']))
            if cumulative_enabled:
                cumulative_sum = np.cumsum(sizes)*width
                cumulative_sums.append(cumulative_sum)
                cumsum_maxima.append(cumulative_sum.max())
            if difference_enabled and previous_sizes is not None:
                size_differences = sizes - previous_sizes
                all_size_differences.append(size_differences)
                overall_max = max(size_differences.max(), overall_max)
                overall_min = min(size_differences.min(), overall_min)
            overall_max = max(sizes.max(), overall_max)
            overall_min = min(sizes.min(), overall_min)
            previous_sizes = sizes
        # self.bins, self.sizes, self.filtered, self.maxima, self.rejected_candidates = all_bins, all_sizes, all_filtered, all_maxima, all_rejected
        
        # arr = np.vstack([all_bins, all_sizes, all_filtered, all_maxima, all_rejected]).astype('d')
        # arr.tofile(self.tmp_filepath)
        # np.hstack(all_bins).astype('d').tofile('bins')
        # np.hstack(all_sizes).astype('d').tofile('sizes')
        # np.hstack(all_filtered).astype('d').tofile('filtered_sizes')
        # np.hstack(all_maxima).astype('d').tofile('maxima')
        # np.hstack(all_rejected).astype('d').tofile('rejected_maxima')
        tmp_filenames = self.tmp_filenames
        np.save(tmp_filenames['bins'], np.vstack(all_bins))
        np.save(tmp_filenames['sizes'], np.vstack(all_sizes))
        if peaks_enabled:
            np.save(tmp_filenames['filtered_sizes'], np.vstack(all_filtered))
            # np.save(tmp_filenames['maxima'], np.vstack(all_maxima))
            # np.save(tmp_filenames['rejected_maxima'], np.vstack(all_rejected))
            self.maxima = all_maxima
            self.rejected_maxima = all_rejected
            np.save(tmp_filenames['top_nm'], np.vstack(all_top_nm))
            np.save(tmp_filenames['fulldata_size_sums'], fulldata_size_sums)
        if cumulative_enabled:
            np.save(tmp_filenames['cumulative_sums'], np.vstack(cumulative_sums))
            np.save(tmp_filenames['cumsum_maxima'], cumsum_maxima)
        if difference_enabled:
            np.save(tmp_filenames['size_differences'], np.vstack(all_size_differences))
        self.overall_min, self.overall_max = overall_min, overall_max
    # def compare(self, )
    # def plot(self, grid_color = '0.8'):
    #     mpl.rcParams["figure.figsize"] = self.figsize

    def run_program(self, grid_color = '0.8'):
        # num_of_plots, samples, colors, peak_settings, unordered_samples, tmp_filepath = self.num_of_plots, self.samples, self.colors, self.peak_settings, self.unordered_samples, self.tmp_filepath
        num_of_plots, samples, colors, table_settings, peak_settings, unordered_samples, overall_min, overall_max, output_folder = self.num_of_plots, self.samples, self.colors, self.table_settings, self.peak_settings, self.unordered_samples, self.overall_min, self.overall_max, self.output_folder
        # all_bins, all_sizes, all_filtered, all_maxima, all_rejected = self.bins, self.sizes, self.filtered, self.maxima, self.rejected_candidates
        peaks_enabled = (peak_settings is not None)
        table_enabled = (table_settings is not None)
        cumulative_enabled, difference_enabled = self.cumulative_enabled, self.difference_enabled
        (_, height) = self.figsize
        rejected_maxima_marker, maxima_marker, filter_description, maxima_candidate_description, maxima_description = peak_settings['rejected_maxima_marker'], peak_settings['maxima_marker'], peak_settings['filter_description'], peak_settings['maxima_candidate_description'], peak_settings['maxima_description']
        tmp_filenames = self.tmp_filenames
        all_bins, all_sizes = np.load(tmp_filenames['bins']+'.npy'), np.load(tmp_filenames['sizes']+'.npy')
        if peaks_enabled:
            all_filtered = np.load(tmp_filenames['filtered_sizes']+'.npy')
            # all_maxima = np.load(tmp_filenames['maxima']+'.npy')
            # all_rejected = np.load(tmp_filenames['rejected_maxima']+'.npy')
            all_maxima, all_rejected = self.maxima, self.rejected_maxima
            all_top_nm = np.load(tmp_filenames['top_nm']+'.npy')
            fulldata_size_sums = np.load(tmp_filenames['fulldata_size_sums']+'.npy')
        if cumulative_enabled:
            cumulative_sums = np.load(tmp_filenames['cumulative_sums']+'.npy')
            cumsum_maxima = np.load(tmp_filenames['cumsum_maxima']+'.npy')
        if difference_enabled:
            all_size_differences = np.load(tmp_filenames['size_differences']+'.npy')

        fig, axs = plt.subplots(num_of_plots, 1)
        fig.subplots_adjust(hspace=-0.05*height)
        transFigure = fig.transFigure

        final_i = num_of_plots - 1
        sums = []
        bins = None
        for i, ax in enumerate(axs):
            sample = samples[i]
            bins, sizes = all_bins[i], all_sizes[i]
            width = bins[1] - bins[0]
            bin_centers = bins + width/2
            
            plt.sca(ax)
            plt.bar(bins, sizes, width = width, color = colors[i], alpha = 0.7, align = 'edge')
            
            if peaks_enabled:
                filtered, maxima, rejected_candidates = all_filtered[i], all_maxima[i], all_rejected[i]
                if len(rejected_candidates) != 0:
                    plt.plot(bin_centers[rejected_candidates], filtered[rejected_candidates], **rejected_maxima_marker)
                plt.plot(bin_centers[maxima], filtered[maxima], **maxima_marker)
            if table_enabled:
                # full_data = pd.read_csv(sample.dat, sep = '\t ', engine = 'python')
                # sums.append((
                #     ('All data', np.sum(full_data['PSD_corrected_[counts/mL/nm]'])*width),
                #     (top_nm, np.sum(sizes*width))
                # ))
                top_nm = all_top_nm[i]
                fulldata_size_sum = fulldata_size_sums[i]
                sums.append((
                    ('All data', fulldata_size_sum*width),
                    (top_nm, np.sum(sizes*width))
                ))
            
            if difference_enabled and i != 0:
                size_differences = all_size_differences[i]
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
            
        if cumulative_enabled:
            max_of_cumulative_sums = cumsum_maxima.max()
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
        text_y = 0 + text_shift
        if difference_enabled:
            plt.text(0, text_y, "Shadows show difference between a plot and the one above it.", fontsize=12, transform = transFigure, verticalalignment = 'center')
        if peaks_enabled:
            text_y -= 0.02
            plt.text(0, text_y, filter_description, fontsize=12, transform = transFigure, verticalalignment = 'center')
        if cumulative_enabled:
            text_y -= 0.02
            plt.text(0, text_y, f"Red lines are cumulative sums of unsmoothed data, scaled by {cumulative_sum_scaling:.3}.", fontsize=12, transform = transFigure, verticalalignment = 'center')
        if peaks_enabled:
            icon_x = 0.01
            text_x = 0.02

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

        if table_enabled:
            axis_positions = [origin[1] for origin in origins]
            cell_height = axis_positions[0] - axis_positions[1]
            table_top = axis_positions[0] + 0.5*cell_height
            table_bottom = axis_positions[-1] - 0.5*cell_height
            results_for_table = results_for_csv
            edges = {'right': right_edge_figure, 'bottom': table_bottom, 'top': table_top}
            draw_table(fig, ax, settings, previous_setting, samples, unordered_samples, sums, results_for_table, edges, table_settings, grid_color)

        fig.savefig(f"{output_folder}/Ridgeline plot.png", dpi = 300, bbox_inches='tight')