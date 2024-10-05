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
from dataclasses import dataclass
import typing
import matplotlib as mpl
resolution = 200
mpl.rcParams['figure.dpi'] = resolution
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.left'] = False
mpl.rcParams['axes.spines.right'] = False
from matplotlib import pyplot as plt, cm
from .sample_class import Sample
from .settings_classes import Calculation, Setting, Settings
from .InfoComparison import compare_info
from .DrawTable import Table
# from nanotracking import data_handler

volume = 2.3E-06
# x_lim = 400
prefix = 'ConstantBinsTable_'
prefix2 = 'Videos_'
suffix = '.dat'
grid_proportion_of_figure = 0.9
text_shift = 0.05
                    
@dataclass
class NeedRefresh:
    settings: set
    calculations: set
    data: bool
    tabulation: bool
    peaks: bool
    cumulative: bool
    difference: bool

class NTA():
    def __init__(self, datafolder, output_folder, filenames, truncation_size):
        self.datafolder, self.output_folder, self.filenames = datafolder, output_folder, filenames
        os.makedirs(output_folder, exist_ok = True)
        self.table, self.peak_settings = None, None
        self.table_enabled, self.cumulative_enabled, self.difference_enabled = False, False, False
        def generate_samples():
            for folder in os.listdir(datafolder):
                folder_path = os.path.join(datafolder, folder)
                if os.path.isfile(folder_path): continue
                sample = Sample(folder_path, prefix, suffix, videos_file_prefix = prefix2)
                if sample.filename not in filenames: continue
                yield sample.filename, sample
        unordered_samples = dict(generate_samples())
        samples = []
        for i, name in enumerate(filenames):
            sample = unordered_samples[name]
            sample.index = i
            samples.append(sample)
        self.samples = samples
        self.truncation_index = self.find_truncation_index(truncation_size)
        samples_setting = Setting('sample', datatype = Sample)
        for sample in samples:
            samples_setting.set_value(sample, sample)
        self.samples_setting = samples_setting
        num_of_plots = len(samples)
        width, height = mpl.rcParamsDefault["figure.figsize"]
        height *= (num_of_plots/3)
        height = min(np.floor(65536/resolution), height)
        self.figsize = (width, height)
        self.colors = cm.plasma(np.linspace(0, 1, num_of_plots))
        self.unordered_samples, self.num_of_plots = unordered_samples, num_of_plots
        self.overall_min, self.overall_max = None, None
        self.maxima, self.rejected_maxima = None, None
        self.bin_width = None
        self.settings = None
        self.tmp_filenames = {
            'bins': os.path.join(output_folder, 'bins'),
            'sizes': os.path.join(output_folder, 'sizes'),
            'filtered_sizes': os.path.join(output_folder, 'filtered_sizes'),
            'top_nm': os.path.join(output_folder, 'top_nm'),
            'fulldata_size_sums': os.path.join(output_folder, 'fulldata_size_sums'),
            'cumulative_sums': os.path.join(output_folder, 'cumulative_sums'),
            'cumsum_maxima': os.path.join(output_folder, 'cumsum_maxima'),
            'size_differences': os.path.join(output_folder, 'size_differences')
        }
        self.calculations = dict()
        self.need_recompute = True
        self.need_refresh = NeedRefresh(settings = set(), calculations = set(), data = True, tabulation = True, peaks = False, cumulative = False, difference = False)
        self.configure_settings()
    def find_truncation_index(self, truncation_size):
        '''
        Given a maximum particle size (truncation_size), finds the lowest index of the array returned by NTA.data_of_sample(sample, truncated = False) such that all sizes are below truncation_size.
        This index doesn't apply to NTA.bins(), NTA.sizes(), etc.!
        '''
        truncation_index = -1
        for sample in self.samples:
            sample_data = self.data_of_sample(sample, truncated = False)
            _truncation_index = np.argmax(sample_data['UpperBinDiameter_[nm]'] > truncation_size)
            truncation_index = max(truncation_index, _truncation_index)
        return truncation_index
    def add_table(self, width, margin_minimum_right, margin_left):
        assert self.table is None, "Table already exists; must first call NTA.delete_table()."
        table = Table(self, width, margin_minimum_right, margin_left)
        self.table = table
        self.need_recompute = True
        self.need_refresh.tabulation = True
        return table
    def delete_table(self):
        self.table = None
    def enable_table(self):
        self.table_enabled = True
    def disable_table(self):
        self.table_enabled = False
    def get_setting_or_calculation(self, tag):
        output = self.settings.by_tag(tag)
        if output is None: output = self.calculations[tag]
        assert output is not None, f'Could not find tag "{tag}" in settings or calculations.'
        return output

    def new_calculation(self, name, value_callback, *output_names):
        calculation = Calculation(name, value_callback, *output_names, samples = self.samples)
        self.calculations[name] = calculation
        need_refresh = self.need_refresh
        need_refresh.tabulation = True
        # need_refresh.calculations.add(calculation)
        return calculation
            
    def enable_peak_detection(self, kernel_size, kernel2_size, kernel_std_in_bins, second_derivative_threshold, maxima_marker, rejected_maxima_marker):
        peak_settings = locals(); peak_settings.pop('self')
        x = np.linspace(0, kernel_size, kernel_size)
        gaussian = np.exp(-np.power((x - kernel_size/2)/kernel_std_in_bins, 2)/2)/(kernel_std_in_bins * np.sqrt(2*np.pi))
        peak_settings['lowpass_filter'] = gaussian / gaussian.sum()
        peak_settings['filter_description'] = f"Black lines indicate Gaussian smoothing (a low-pass filter) with $\sigma = {kernel_std_in_bins}$ bins and convolution kernel of size {kernel_size} bins."
        peak_settings['maxima_candidate_description'] = f": Candidate peaks after smoothing, selected using argrelextrema in SciPy {scipy.__version__}."
        peak_settings['maxima_description'] = f": Peaks with under {second_derivative_threshold} counts/mL/nm$^3$ second derivative, computed after smoothing again with simple moving average of size {kernel2_size} bins."
        self.peak_settings = peak_settings
        self.need_recompute = True
        self.need_refresh.tabulation = True
        self.need_refresh.peaks = True
    def disable_peak_detection(self):
        self.peak_settings = None
    def enable_cumulative(self):
        self.cumulative_enabled = True
        self.need_recompute = True
        self.need_refresh.tabulation = True
        self.need_refresh.cumulative = True
    def disable_cumulative(self):
        self.cumulative_enabled = False
    def enable_difference(self):
        self.difference_enabled = True
        self.need_recompute = True
        self.need_refresh.tabulation = True
        self.need_refresh.difference = True
    def disable_difference(self):
        self.difference_enabled = False
    def configure_settings(self):
        previous_setting = Setting('previous', name = 'Previous')
        md_settings = [
            Setting('experimental_unit', name = 'Experimental unit'),
            Setting('treatment', name = 'Treatment', units = 'µM'),
            Setting('wait', name = 'Wait', units = 'h'),
            Setting('filter', name = 'Filter cut-on', units = 'nm'),
            previous_setting ]
        red_enabled = Setting('RedLaserEnabled', name = 'Red enabled', datatype = bool)
        green_enabled = Setting('GreenLaserEnabled', name = 'Green enabled', datatype = bool)
        blue_enabled = Setting('BlueLaserEnabled', name = 'Blue enabled', datatype = bool)
        detection_threshold_setting = Setting('DetectionThresholdType', name = 'Detection mode', dependencies_require = 'Manual')
        xml_settings = [
            Setting('RedLaserPower', short_name = '635nm', name = '635nm power', units = 'mW', datatype = int, show_name = True, depends_on = red_enabled),
            red_enabled,
            Setting('GreenLaserPower', short_name = '520nm', name = '520nm power', units = 'mW', datatype = int, show_name = True, depends_on = green_enabled),
            green_enabled,
            Setting('BlueLaserPower', short_name = '445nm', name = '445nm power', units = 'mW', datatype = int, show_name = True, depends_on = blue_enabled),
            blue_enabled,
            Setting('Exposure', units = 'ms', datatype = int),
            Setting('Gain', units = 'dB', datatype = int),
            Setting('MeasurementStartDateTime'),
            Setting('FrameRate', name = 'Framerate', units = 'fps', datatype = int),
            Setting('FramesPerVideo', name = 'Frames per video', units = 'frames', datatype = int),
            Setting('NumOfVideos', name = 'Number of videos', datatype = int),
            Setting('StirrerSpeed', name = 'Stirring speed', units = 'rpm', datatype = int),
            Setting('StirredTime', name = 'Stirred time', units = 's', datatype = int),
            detection_threshold_setting,
            Setting('DetectionThreshold', name = 'Detection threshold', datatype = float, depends_on = detection_threshold_setting) ]
        settings_list = [self.samples_setting, *md_settings, *xml_settings]
        settings = Settings(OrderedDict({setting.tag: setting for setting in settings_list}))
        self.settings = settings
    # def refresh(self):
    #     """
    #     Applies any changes made to this NTA object, such as the addition of a table to the plot.
    #     """
    #     need_refresh = self.need_refresh

    def data_of_sample(self, sample, truncated = True, truncation_size = None):
        """
        Get the size distribution of a sample, using either its name (a string) or the Sample object itself.
        Returns a Pandas DataFrame.
        If truncated is True, only data with sizes below truncation_size will be returned.
        If truncation_size is None, NTA.truncation_size will be used.
        """
        if type(sample) is not Sample:
            assert type(sample) is str, "Argument of data_of_sample must be of type Sample or string."
            sample = self.unordered_samples[sample]
        all_data = pd.read_csv(sample.dat, sep = '\t ', engine = 'python')
        if truncated:
            if truncation_size is None:
                truncation_index = self.truncation_index
            else:
                truncation_index = self.find_truncation_index(truncation_size)
            return all_data.iloc[:truncation_index, :]
        return all_data
    def compute(self, prep_tabulation = True):
        need_refresh = self.need_refresh
        def vstack(arrays):
            if len(arrays) == 0:
                try: return arrays[0]
                except: return arrays # In case "arrays" is an empty list.
            return np.vstack(arrays)
        
        peak_settings = self.peak_settings
        peaks_enabled = (peak_settings is not None)
        refresh_peaks = (peaks_enabled and need_refresh.peaks)
        if refresh_peaks:
            lowpass_filter, kernel2_size, second_derivative_threshold = peak_settings['lowpass_filter'], peak_settings['kernel2_size'], peak_settings['second_derivative_threshold']
            all_filtered, all_maxima, all_rejected  = [], [], []
        refresh_cumulative = (self.cumulative_enabled and need_refresh.cumulative)
        if refresh_cumulative:
            cumulative_sums = []
            cumsum_maxima = []
        refresh_difference = (self.difference_enabled and need_refresh.difference)
        if refresh_difference:
            all_size_differences = []
        refresh_data = need_refresh.data

        overall_min, overall_max = 0, 0
        previous_sizes = None
        last_bins = None
        if refresh_data:
            all_bins, all_sizes = [], []
        else:
            all_bins, all_sizes = self.bins(), self.sizes()
            width = self.bin_width
        data_of_sample = self.data_of_sample
        for i, sample in enumerate(self.samples):
            if not refresh_data:
                bins, sizes = all_bins[i], all_sizes[i]
            else:
                data = data_of_sample(sample, truncated = True)
                bins, sizes = data['/LowerBinDiameter_[nm]'], data['PSD_corrected_[counts/mL/nm]']
                width = bins[1] - bins[0]
                
                if last_bins is not None:
                    assert np.all(bins == last_bins) == True, 'Unequal sequence of bins between samples detected!'
                last_bins = bins
                # data_handler.parse_data(bins.to_numpy(dtype = np.double), sizes.to_numpy(dtype = np.double), sample.filename, self.output_folder, num_data_points)
                # data_handler.parse_data(
                #     bins = bins.to_numpy(dtype = np.double),
                #     sizes = sizes.to_numpy(dtype = np.double),
                #     sample_filename = sample.filename,
                #     outputs_path = self.output_folder,
                #     num_data_points = num_data_points)
                all_bins.append(bins)
                all_sizes.append(sizes)

            # videos = sample.videos
            # all_histograms = np.array([np.histogram(video, bins = bins)[0] for video in videos])
            # avg_histogram = np.average(all_histograms, axis = 0)
            # total_std = np.std(all_histograms, axis = 0, ddof = 1)
            # scale_factor = np.array([sizes[j]/avg if (avg := avg_histogram[j]) != 0 else 0 for j in range(len(sizes)-1)])
            # error_resizing = 0.1
            # total_std *= scale_factor * error_resizing
            # avg_histogram *= scale_factor
            # total_stds.append(total_std)
            # avg_histograms.append(avg_histograms)

            if refresh_peaks:
                assert len(lowpass_filter) <= len(sizes), f"kernel_size={len(lowpass_filter)} is too big, given {len(sizes)=}."
                assert kernel2_size <= len(sizes), f"{kernel2_size=} is too big, given {len(sizes)=}."
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
            if refresh_cumulative:
                cumulative_sum = np.cumsum(sizes)*width
                cumulative_sums.append(cumulative_sum)
                cumsum_maxima.append(cumulative_sum.max())
            if refresh_difference and previous_sizes is not None:
                size_differences = sizes - previous_sizes
                all_size_differences.append(size_differences)
                overall_max = max(size_differences.max(), overall_max)
                overall_min = min(size_differences.min(), overall_min)
            if refresh_data or refresh_difference:
                overall_max = max(sizes.max(), overall_max)
                overall_min = min(sizes.min(), overall_min)
            for setting in tuple(need_refresh.settings):
                value = setting.value_callback(sample)
                setting.set_value(sample, value)
                need_refresh.settings.remove(setting)
            for calculation in tuple(need_refresh.calculations):
                calculation.refresh(sample)
                need_refresh.calculations.remove(calculation)
            previous_sizes = sizes
        
        tmp_filenames = self.tmp_filenames
        if refresh_data:
            np.save(tmp_filenames['bins'], vstack(all_bins))
            np.save(tmp_filenames['sizes'], vstack(all_sizes))
            self.overall_min, self.overall_max = overall_min, overall_max
            self.bin_width = last_bins[1] - last_bins[0]
        # np.save(tmp_filenames['total_stds'], vstack(total_stds))
        # np.save(tmp_filenames['avg_histograms'], vstack(avg_histograms))
        if refresh_peaks:
            np.save(tmp_filenames['filtered_sizes'], vstack(all_filtered))
            self.maxima = all_maxima
            self.rejected_maxima = all_rejected
        if refresh_cumulative:
            np.save(tmp_filenames['cumulative_sums'], vstack(cumulative_sums))
            np.save(tmp_filenames['cumsum_maxima'], cumsum_maxima)
        if refresh_difference:
            np.save(tmp_filenames['size_differences'], vstack(all_size_differences))
        self.need_recompute = False
        if prep_tabulation:
            self.prepare_tabulation()
    def prepare_tabulation(self):
        table, settings, num_of_plots, samples, unordered_samples = self.table, self.settings, self.num_of_plots, self.samples, self.unordered_samples
        table_enabled = self.table_enabled
        if table_enabled:
            self.table.reset_columns() # If prepare_tabulation() has been run before, remove the columns for treatments and waits.
            column_widths = table.column_widths
            column_names = table.column_names
            include_experimental_unit = table.include_experimental_unit
            treatments_and_waits = table.treatments_and_waits
            include_treatments = (treatments_and_waits is not None)
            if include_treatments:
                treatments_waits_columnIndex = treatments_and_waits[0]
            else:
                treatments_waits_columnIndex = -1
        
        def generate_rows():
            column_quantities = dict()
            def number_of_subtags(tag):
                if (setting := settings.by_tag(tag)) is None: return 0
                return max(len(setting.subsettings), 1)
            def get_multivalued(tag, sample):
                if (setting := settings.by_tag(tag)) is None: return []
                if len(setting.subsettings) == 0:
                    value = setting.get_value(sample)
                    if value is None: return []
                    return [value]
                subsettings = list(setting.numbered_subsettings.items())
                subsettings.sort()
                values = [subsetting.get_value(sample) for _, subsetting in subsettings]
                if values[0] is None:
                    values[0] = setting.get_value(sample)
                return values
            for i in range(num_of_plots):
                sample = samples[i]
                settings.read_files(sample)
                settings.parse_time(sample)
                if table_enabled and include_treatments:
                    for tag in ('treatment', 'wait'):
                        quantity = number_of_subtags(tag)
                        if tag not in column_quantities:
                            column_quantities[tag] = quantity
                            continue
                        column_quantities[tag] = max(column_quantities[tag], quantity)
            if table_enabled:
                for i in range(len(column_names) + 1): # +1 accounts for the case where len(column_names) = 1. Still may want to insert treatments_and_waits columns at index 0.
                    if i == treatments_waits_columnIndex:
                        num_of_treatments = column_quantities['treatment']
                        num_of_waits = column_quantities['wait']
                        treatment_column_name, treatment_column_width = treatments_and_waits[1]
                        wait_column_name, wait_column_width = treatments_and_waits[2]
                        index = 0
                        for j in range(max(num_of_treatments, num_of_waits)):
                            if j < num_of_treatments:
                                column_names.insert(i + index, treatment_column_name.format(treatment_number = j + 1))
                                column_widths.insert(i + index, treatment_column_width)
                                index += 1
                            if j < num_of_waits:
                                column_names.insert(i + index, wait_column_name.format(wait_number = j + 1))
                                column_widths.insert(i + index, wait_column_width)
                                index += 1
                            
            for i in range(num_of_plots):
                row = []
                sample = samples[i]
                if table_enabled:
                    if include_treatments:
                        treatments = get_multivalued('treatment', sample)
                        waits = get_multivalued('wait', sample)
                        for j in range( max(column_quantities['treatment'], column_quantities['wait']) ):
                            if j < len(treatments): row.append(treatments[j])
                            elif j < column_quantities['treatment']: row.append(None)
                            if j < len(waits): row.append(waits[j])
                            elif j < column_quantities['wait']: row.append(None)
                    if include_experimental_unit:
                        experimental_unit = settings.by_tag('experimental_unit')
                        text = ''
                        if experimental_unit is not None:
                            value = experimental_unit.get_value(sample)
                            text += value if value is not None else ''
                            if hasattr(experimental_unit, 'age'):
                                age = experimental_unit.age.get_value(sample)
                                text += f"\n{age:.1f} d old" if age is not None else ''
                        row.append(text)
                    columns = list(table.columns_as_Settings_object.column_numbers.items())
                    columns.sort()
                    for j, column in columns:
                        # content = '\n'.join(
                        #     setting.show_name*f"{setting.short_name}: " + f"{setting.get_value(sample)}" + setting.show_unit*f" ({setting.units})"
                        #     for setting in column if setting.get_value(sample) is not None )
                        # row.append(content)
                        assert len(column) == 1, "There can be only one Setting object per column."
                        if column[0].tag.startswith('COLUMN'):
                            group = column[0]
                            grouped_settings = group.subsettings.values()
                            if group.format_callback is not None:
                                row.append(group.format_callback(*(setting.get_value(sample) for setting in grouped_settings)))
                                continue
                            row.append(group.format_string.format(**{setting.tag: setting.get_value(sample) for setting in grouped_settings}))
                        elif column[0].tag.startswith('CALC'):
                            group = column[0]
                            # grouped_settings, values = group.subsettings.values(), group.get_value(sample) # These line up; will use zip below
                            grouped_settings = group.subsettings.values()
                            values = [subsetting.get_value(sample) for subsetting in grouped_settings]
                            if group.format_callback is not None:
                                row.append(group.format_callback(*values))
                                continue
                            row.append(group.format_string.format(**{setting.tag: value for setting, value in zip(grouped_settings, values)}))
                        else:
                            setting = column[0]
                            value = setting.get_value(sample)
                            if setting.format_callback is None:
                                row.append(setting.format_string.format(**{setting.tag: value}))
                            else:
                                row.append(setting.format_callback(value))

                
                # if detection_threshold is None:
                #     row.append(detection_mode)
                # else:
                #     row.append(f"{detection_mode}\n{detection_threshold}")
                # results_for_csv.previous.set_value(sample, previous)
                # results_for_csv.time_since_previous.set_value(sample, time_since_previous)
                # results_for_csv.total_conc.set_value(sample, f"{data_sums[1][1]:.2E}")
                # results_for_csv.total_conc_under_topnm.set_value(sample, f"{data_sums[0][1]:.2E}")
                row.append("")
                yield row
        self.rows = tuple(generate_rows())
        self.need_refresh.tabulation = False
    def compare(self):
        assert self.need_recompute == False, "Must run NTA.compute() first."
        assert self.need_refresh.tabulation == False, "Must run NTA.prepare_tabulation() first."
        compare_info(self.settings, self.samples, self.calculations, self.output_folder)
    def load_tempfile(self, key):
        filename = self.tmp_filenames[key]+'.npy'
        try:
            data = np.load(filename)
        except OSError:
            raise Exception(f"Couldn't find {filename}.")
        return data
    def bins(self, sample = None, truncation_size = None, lower = True, middle = False, upper = False):
        assert (lower + middle + upper) == 1, "Must set either lower, middle, or upper to True, but not more than one."
        
        if truncation_size is not None:
            assert sample is not None, "Must specify a sample when using truncation_size."
            upper_bins = self.bins(sample = sample, truncation_size = None, lower = False, upper = True)
            truncation_index = np.argmax(upper_bins > truncation_size)
        
        all_bins = self.load_tempfile('bins')
        all_bins += (0.5*middle + 1*upper) * self.bin_width
        if sample is not None:
            assert sample.index is not None, "The provided Sample object has no index."
            sample_bins = all_bins[sample.index]
            if truncation_size is not None:
                return sample_bins[:truncation_index]
            return sample_bins
        return all_bins
    def sizes(self, sample = None, truncation_size = None):
        all_sizes = self.load_tempfile('sizes')
        
        if truncation_size is not None:
            assert sample is not None, "Must specify a sample when using truncation_size."
            upper_bins = self.bins(sample = sample, truncation_size = None, lower = False, upper = True)
            truncation_index = np.argmax(upper_bins > truncation_size)
        
        if sample is not None:
            assert sample.index is not None, "The provided Sample object has no index."
            sample_sizes = all_sizes[sample.index]
            if truncation_size is not None:
                return sample_sizes[:truncation_index]
            return sample_sizes
        return all_sizes
    def plot(self, grid_color = '0.8', name = 'Ridgeline plot'):
        assert self.need_recompute == False, "Must run NTA.compute() first."
        num_of_plots, samples, colors, table, peak_settings, overall_min, overall_max, output_folder = self.num_of_plots, self.samples, self.colors, self.table, self.peak_settings, self.overall_min, self.overall_max, self.output_folder
        peaks_enabled = (peak_settings is not None)
        table_enabled = self.table_enabled
        if table_enabled:
            assert self.need_refresh.tabulation == False, "Must run NTA.prepare_tabulation() first."
        cumulative_enabled, difference_enabled = self.cumulative_enabled, self.difference_enabled
        (_, height) = self.figsize
        tmp_filenames = self.tmp_filenames
        all_bins, all_sizes = self.bins(), self.sizes()
        # avg_histograms, total_stds = np.load(tmp_filenames['avg_histograms']+'.npy'), np.load(tmp_filenames['total_stds']+'.npy')
        if peaks_enabled:
            rejected_maxima_marker, maxima_marker, filter_description, maxima_candidate_description, maxima_description = peak_settings['rejected_maxima_marker'], peak_settings['maxima_marker'], peak_settings['filter_description'], peak_settings['maxima_candidate_description'], peak_settings['maxima_description']
            all_filtered = np.load(tmp_filenames['filtered_sizes']+'.npy')
            all_maxima, all_rejected = self.maxima, self.rejected_maxima
        if cumulative_enabled:
            cumulative_sums = np.load(tmp_filenames['cumulative_sums']+'.npy')
            cumsum_maxima = np.load(tmp_filenames['cumsum_maxima']+'.npy')
            max_of_cumulative_sums = cumsum_maxima.max()
            cumulative_sum_scaling = overall_max / max_of_cumulative_sums
        if difference_enabled:
            all_size_differences = np.load(tmp_filenames['size_differences']+'.npy')

        mpl.rcParams["figure.figsize"] = self.figsize
        fig, axs = plt.subplots(num_of_plots, 1, squeeze = False)
        axs = axs[:,0]  # Flatten axs from a 2D array (of size num_of_plots x 1) to a 1D array
        fig.subplots_adjust(hspace=-0.05*height)
        transFigure = fig.transFigure
        transFigure_inverted = transFigure.inverted()

        final_i = num_of_plots - 1
        origins = []
        for i, ax in enumerate(axs):
            sample = samples[i]
            bins, sizes = all_bins[i], all_sizes[i]
            # bins, sizes = data_handler.read_data(sample_filename = sample.filename, outputs_path = output_folder, num_data_points = num_data_points)
            width = bins[1] - bins[0]
            bin_centers = bins + width/2
            # avg_histogram, total_std = avg_histograms[i], total_stds[i]
            
            plt.sca(ax)
            plt.bar(bins, sizes, width = width, color = colors[i], alpha = 0.7, align = 'edge')
            
            if peaks_enabled:
                filtered, maxima, rejected_candidates = all_filtered[i], all_maxima[i], all_rejected[i]
                plt.plot(bins, filtered, linewidth = 0.5, color = 'black')
                if len(rejected_candidates) != 0:
                    plt.plot(bin_centers[rejected_candidates], filtered[rejected_candidates], **rejected_maxima_marker)
                plt.plot(bin_centers[maxima], filtered[maxima], **maxima_marker)
            
            if difference_enabled and i != 0:
                size_differences = all_size_differences[i-1]
                plt.bar(bins, size_differences, width = width, color = 'black', alpha = 0.3, align = 'edge')
            
            videos = sample.videos
            all_histograms = np.array([np.histogram(video, bins = bins)[0] for video in videos])
            avg_histogram = np.average(all_histograms, axis = 0)
            total_std = np.std(all_histograms, axis = 0, ddof = 1)
            scale_factor = np.array([sizes[j]/avg if (avg := avg_histogram[j]) != 0 else 0 for j in range(len(sizes)-1)])
            error_resizing = 0.1
            total_std *= scale_factor * error_resizing
            avg_histogram *= scale_factor
            errorbars = np.array(list(zip(total_std, [0]*len(total_std)))).T
            plt.errorbar(bin_centers[:-1], avg_histogram, yerr = errorbars, elinewidth = 1, linestyle = '', marker = '.', ms = 1, alpha = 0.5, color = 'black')            
            
            # plt.xlim(0, x_lim)
            plt.ylim(overall_min, overall_max)
            ax.patch.set_alpha(0)
                
            if i == final_i:
                ax.yaxis.get_offset_text().set_x(-0.1)
                plt.xlabel("Diameter (nm)")
                plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))
                ax.spines.left.set_visible(True)
            else:
                ax.spines['bottom'].set_position(('data', 0))
                plt.yticks([])
                plt.xticks([])

            origin_transDisplay = ax.transData.transform([0, 0])
            origins.append(transFigure_inverted.transform(origin_transDisplay))
            
            if cumulative_enabled:
                cumulative_sum = cumulative_sums[i]
                plt.plot(bins, cumulative_sum*cumulative_sum_scaling, color = 'red', linewidth = 0.5)

        final_ax = ax   # Using the ax from the final (bottom) plot:
        final_ax.xaxis.set_tick_params(width = 2)
        final_ax.yaxis.set_tick_params(width = 2)
        transData = final_ax.transData
        tick_values, tick_labels = plt.xticks()

        final_i = len(tick_values) - 1
        right_edge_figure = None
        for i, tick_value in enumerate(tick_values):
            display_coords = transData.transform([tick_value, overall_min])
            figure_x, figure_y = transFigure_inverted.transform(display_coords)
            
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

        if table_enabled:
            if len(origins) > 1:
                axis_positions = [origin[1] for origin in origins]
                cell_height = axis_positions[0] - axis_positions[1]
                table_top = axis_positions[0] + 0.5*cell_height
                table_bottom = axis_positions[-1] - 0.5*cell_height
            else:
                cell_height = table_top = 1
                table_bottom = 0
            edges = {'right': right_edge_figure, 'bottom': table_bottom, 'top': table_top}
            table.draw_table(fig, ax, self.rows, edges, grid_color)

        fig.savefig(f"{output_folder}/{name}.png", dpi = 300, bbox_inches='tight')