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
from copy import deepcopy
import typing
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
# from nanotracking import data_handler

volume = 2.3E-06
x_lim = 400
num_data_points = 400
prefix = 'ConstantBinsTable_'
prefix2 = 'Videos_'
suffix = '.dat'
grid_proportion_of_figure = 0.9
text_shift = 0.05
                    
class NTA():
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
        samples = []
        for i, name in enumerate(filenames):
            sample = unordered_samples[name]
            sample.index = i
            samples.append(sample)
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
        self.unordered_samples, self.samples, self.num_of_plots = unordered_samples, samples, num_of_plots
        self.overall_min, self.overall_max = None, None
        self.maxima, self.rejected_maxima = None, None
        self.settings = None
        self.tmp_filenames = {
            'bins': os.path.join(output_folder, 'bins'),
            'sizes': os.path.join(output_folder, 'sizes'),
            'filtered_sizes': os.path.join(output_folder, 'filtered_sizes'),
            'top_nm': os.path.join(output_folder, 'top_nm'),
            'fulldata_size_sums': os.path.join(output_folder, 'fulldata_size_sums'),
            'cumulative_sums': os.path.join(output_folder, 'cumulative_sums'),
            'cumsum_maxima': os.path.join(output_folder, 'cumsum_maxima'),
            'size_differences': os.path.join(output_folder, 'size_differences')#,
            # 'avg_histograms': os.path.join(output_folder, 'avg_histograms'),
            # 'total_stds': os.path.join(output_folder, 'total_stds')
        }
        # self.setting_columnNames = []
        # self.setting_columnWidths = []
        self.results_for_csv = Settings()
        self.result_names = {
            'time': None,
            'concentration': None
        }
        self.sums = []
        self.need_recompute = True
        self.need_reprep_tabulation = True
        self.configure_settings()
    def enable_table(self, width, margin_minimum_right, margin_left):
        table_settings = locals(); table_settings.pop('self')
        table_settings.update({
            'include_experimental_unit': False,
            'treatments_and_waits': None,
            # 'setting_tags': [],
            'columns_as_Settings_object': Settings(),
            'column_names': [],
            'column_widths': [],
            'column_names_without_treatmentsOrWaits': [],
            'column_widths_without_treatmentsOrWaits': []
        })
        self.table_settings = table_settings
        self.need_recompute = True
        self.need_reprep_tabulation = True
    def disable_table(self):
        self.table_settings = None    
    # # def table_add_setting(self, setting_tag, name, width):
    # # def table_add_setting(self, name, width):
    # def table_add_column(self, name, width):
    #     table_settings = self.table_settings
    #     # # table_settings['setting_tags'].append(setting_tag)
    #     # table_settings['column_names'].append(name)
    #     # table_settings['column_widths'].append(width)
    #     # # self.setting_columnNames.append(name)
    #     # # self.setting_columnWidths.append(width)
    #     table_settings['column_names_without_treatmentsOrWaits'].append(name)
    #     table_settings['column_widths_without_treatmentsOrWaits'].append(width)
    def table_add_setting(self, setting: Setting):
        tag = setting.tag
        settings = self.table_settings['columns_as_Settings_object']
        assert tag not in settings.tags, f'Setting with tag "{tag}" already added to table.'
        if setting.column_number is None:
            setting.column_number = len(settings.column_widths)
        settings.add_setting(setting.tag, setting)
    def table_add_new_setting(self, tag, column_name, column_width, column_number = None):
        assert tag not in self.settings.tags, f'Setting with tag "{tag}" already created.'
        if column_number is None:
            column_number = len(self.table_settings['columns_as_Settings_object'].column_widths)
        setting = Setting(tag, column_name = column_name, column_width = column_width, column_number = column_number)
        self.settings.add_setting(tag, setting)
        self.table_add_setting(setting)
    def get_setting_or_result(self, tag):
        output = self.settings.by_tag(tag)
        if output is None: output = self.results_for_csv.by_tag(tag)
        assert output is not None, f'Could not find tag "{tag}" in settings or results_for_csv.'
        return output
    def table_add_settings_by_tag(self, *tags, column_number = None, column_name = None, column_width = None, format_string = None, format_callback = None):
        '''
        Adds multiple Setting objects to the table.
        Example use case: specify column_number to group all specified settings into one column.

        If column_number is not given, the next available column will be used.
        Note: the column_number values of the specified Setting objects will be overwritten!
        
        If neither format_string nor format_callback are given, then for each cell in the column, the settings' individual format_strings will be used on separate lines.
        To use format_string, reference settings' values using their tags in curly braces: for example, format_string = "Red has power {RedLaserPower}."
        To use format_callback, define a function that accepts settings' values as arguments and returns a formatted (value-containing) string.

        If format_callback is given, it will be used instead of format_string.
        '''
        if column_number is None:
            column_number = len(self.table_settings['columns_as_Settings_object'].column_widths)
        settings = [self.get_setting_or_result(tag) for tag in tags]
        if format_string is None:
            format_string = '\n'.join([setting.format_string for setting in settings])
        def prepare_setting(setting):
            setting.column_number = column_number
            # setting.format_string = format_string
            # setting.format_callback = format_callback
            if column_name is not None: setting.column_name = column_name
            if column_width is not None: setting.column_width = column_width
        if len(settings) == 1:
            setting = settings[0]
            prepare_setting(setting)
            setting.set_attributes(format_string = format_string, format_callback = format_callback)
            self.table_add_setting(setting)
            return
        if format_callback is None:
            group_suffix = format_string  # Allows multiple different format_callbacks or format_strings to be used on the same group, without counting as the same group (which would cause an error)
        else:
            group_suffix = format_callback.__name__
        group = Setting('COLUMN_' + '_'.join(tags) + group_suffix, column_number = column_number, column_name = column_name, column_width = column_width, format_string = format_string, format_callback = format_callback)
        for setting in settings:
            prepare_setting(setting)
            group.add_subsetting(setting, setting.tag)
        self.table_add_setting(group)
    def table_add_results(self, results_group, column_number = None, column_name = None, column_width = None, format_string = None, format_callback = None):
        # tags = results_group.subsettings.keys()
        # return self.table_add_settings_by_tag(*tags, column_number = column_number, column_name = column_name, column_width = column_width, format_string = format_string, format_callback = format_callback)
        
        # for setting in results_group.subsettings.values():
        #     self.table_add_setting(setting)

        if format_callback is None:
            group_suffix = format_string  # Allows multiple different format_callbacks or format_strings to be used on the same group, without counting as the same group (which would cause an error)
        else:
            group_suffix = format_callback.__name__
        new_column = deepcopy(results_group)
        new_column.set_attributes(tag = results_group.tag + group_suffix, column_number = column_number, column_name = column_name, column_width = column_width, format_string = format_string, format_callback = format_callback)
        self.table_add_setting(new_column)
    
    def table_add_experimental_unit(self, column_name = "Experimental\nunit", width = 0.3, column_number = None):
        '''
        Adds to the table a column for experimental unit, whose name is given by "experimental_unit=…" in each sample's info.md file.
        '''
        # self.table_add_column(column_name, width)
        experimental_unit = self.settings.by_tag('experimental_unit')
        if column_number is None:
            column_number = len(self.table_settings['columns_as_Settings_object'].column_widths)
        experimental_unit.column_number = column_number
        experimental_unit.column_name = column_name
        experimental_unit.column_width = width
        self.table_add_setting(experimental_unit)
    def table_add_treatments_and_waits(self, treatments_column_name, treatments_width, waits_column_name, waits_width):
        '''
        For each treatment & wait-time listed in samples' info.md files, adds to the table
        (1) a column for the treatment's name, and (2) a column for the time waited after applying the treatment.
        '''
        table_settings = self.table_settings
        start_index = len(table_settings['column_names_without_treatmentsOrWaits'])
        assert table_settings['treatments_and_waits'] is None, "Treatments and waits have already been added to the table."
        table_settings['treatments_and_waits'] = [start_index, (treatments_column_name, treatments_width), (waits_column_name, waits_width)]
        # column_names_without_treatmentsOrWaits = self.column_names_without_treatmentsOrWaits
        # start_index = len(column_names_without_treatmentsOrWaits)
        # assert self.treatments_and_waits is None, "Treatments and waits have already been added to the table."
        # self.treatments_and_waits = [start_index, (treatments_column_name, treatments_width), (waits_column_name, waits_width)]
    def reset_columns(self):
        table_settings = self.table_settings
        table_settings['column_names'] = list(table_settings['columns_as_Settings_object'].column_names.keys())
        table_settings['column_widths'] = table_settings['columns_as_Settings_object'].column_widths.copy()

    # def table_add_time(self, name, width):
    #     # self.table_add_column(name, width)
    #     self.table_add_settings_by_tag('time', column_name = name, column_width = width)
    #     self.result_names['time'] = name
    # def table_add_concentration(self, name, width):
    #     # self.table_add_column(name, width)
    #     self.table_add_new_setting('concentration', name, width)
    #     self.result_names['concentration'] = name
    
    def new_results_group(self, *tags, callback = None):
        assert callback is not None, "Must specify callback."
        group = Setting('RESULTS_' + '_'.join(tags), value_callback = callback)
        for tag in tags:
            group.add_subsetting(Setting(tag), tag)
        self.results_for_csv.add_setting(group.tag, group)
        return group
    # def table_add_results_by_tag(self, *tags, column_number = None, column_name = None, column_width = None, format_string = None, format_callback = None):
    #     '''
    #     Alias for table_add_settings_by_tag with is_result = True.
    #     '''
    #     self.table_add_settings_by_tag(*tags, is_result = True, column_number = column_number, column_name = column_name, column_width = column_width, format_string = format_string, format_callback = format_callback)
            
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
        self.need_reprep_tabulation = True
    def disable_peak_detection(self):
        self.peak_settings = None
    def enable_cumulative(self):
        self.cumulative_enabled = True
        self.need_recompute = True
        self.need_reprep_tabulation = True
    def disable_cumulative(self):
        self.cumulative_enabled = False
    def enable_difference(self):
        self.difference_enabled = True
        self.need_recompute = True
        self.need_reprep_tabulation = True
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
        # self.need_reconfig_settings = False
    def compute(self, prep_tabulation = True):
        def vstack(arrays):
            if len(arrays) == 0:
                try: return arrays[0]
                except: return arrays # In case "arrays" is an empty list.
            return np.vstack(arrays)
        
        peak_settings = self.peak_settings
        peaks_enabled = (peak_settings is not None)
        if peaks_enabled:
            lowpass_filter, kernel2_size, second_derivative_threshold = peak_settings['lowpass_filter'], peak_settings['kernel2_size'], peak_settings['second_derivative_threshold']
            all_filtered, all_maxima, all_rejected, all_top_nm,  = [], [], [], []
        cumulative_enabled = self.cumulative_enabled
        if cumulative_enabled:
            cumulative_sums = []
            cumsum_maxima = []
        difference_enabled = self.difference_enabled
        if difference_enabled:
            all_size_differences = []

        overall_min, overall_max = 0, 0
        previous_sizes = None
        sums = self.sums
        bins = None
        all_bins, all_sizes = [], []
        # total_stds, avg_histograms = [], []
        fulldata_size_sums = []
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
            # data_handler.parse_data(bins.to_numpy(dtype = np.double), sizes.to_numpy(dtype = np.double), sample.filename, self.output_folder, num_data_points)
            # data_handler.parse_data(
            #     bins = bins.to_numpy(dtype = np.double),
            #     sizes = sizes.to_numpy(dtype = np.double),
            #     sample_filename = sample.filename,
            #     outputs_path = self.output_folder,
            #     num_data_points = num_data_points)
            width = bins[1] - bins[0]
            fulldata_size_sum = np.sum(full_data['PSD_corrected_[counts/mL/nm]'])
            
            all_bins.append(bins)
            all_sizes.append(sizes)
            fulldata_size_sums.append(fulldata_size_sum)
            sums.append((
                ('All data', fulldata_size_sum*width),
                (top_nm, np.sum(sizes*width))
            ))

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
        
        tmp_filenames = self.tmp_filenames
        np.save(tmp_filenames['bins'], vstack(all_bins))
        np.save(tmp_filenames['sizes'], vstack(all_sizes))
        np.save(tmp_filenames['fulldata_size_sums'], fulldata_size_sums)
        # np.save(tmp_filenames['total_stds'], vstack(total_stds))
        # np.save(tmp_filenames['avg_histograms'], vstack(avg_histograms))
        if peaks_enabled:
            np.save(tmp_filenames['filtered_sizes'], vstack(all_filtered))
            self.maxima = all_maxima
            self.rejected_maxima = all_rejected
            np.save(tmp_filenames['top_nm'], vstack(all_top_nm))
        if cumulative_enabled:
            np.save(tmp_filenames['cumulative_sums'], vstack(cumulative_sums))
            np.save(tmp_filenames['cumsum_maxima'], cumsum_maxima)
        if difference_enabled:
            np.save(tmp_filenames['size_differences'], vstack(all_size_differences))
        self.overall_min, self.overall_max = overall_min, overall_max
        self.need_recompute = False
        if prep_tabulation:
            self.prepare_tabulation()
    def prepare_tabulation(self):
        # assert self.need_reconfig_settings == False, "Must run NTA.configure_settings() first."
        sums = self.sums
        table_settings, settings, num_of_plots, samples, unordered_samples = self.table_settings, self.settings, self.num_of_plots, self.samples, self.unordered_samples
        table_enabled = (table_settings is not None)
        if table_enabled:
            self.reset_columns() # If prepare_tabulation() has been run before, remove the columns for treatments and waits.
            column_widths = table_settings['column_widths']
            column_names = table_settings['column_names']
            include_experimental_unit = table_settings['include_experimental_unit']
            treatments_and_waits = table_settings['treatments_and_waits']
            include_treatments = (treatments_and_waits is not None)
            if include_treatments:
                treatments_waits_columnIndex = treatments_and_waits[0]
            else:
                treatments_waits_columnIndex = -1
        result_names = self.result_names
        time_enabled, concentration_enabled = (result_names['time'] is not None), (result_names['concentration'] is not None)
        previous_setting = settings.by_tag('previous')
        
        results_for_csv = self.results_for_csv
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
            top_nm = None
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
                data_sums = sums[i]
                assert len(data_sums) == 2
                if top_nm is None:
                    top_nm, _ = data_sums[1]
                assert data_sums[1][0] == top_nm
            if table_enabled:
                # for i, name in enumerate(column_names):
                for i in range(len(column_names) + 1): # +1 accounts for the case where len(column_names) = 1. Still may want to insert treatments_and_waits columns at index 0.
                    if i < len(column_names):
                        name = column_names[i]
                        if '{top_nm}' in name:
                            column_names[i] = name.format(top_nm = top_nm)
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
            
            # results_for_csv.add_subsetting(previous_setting, 'previous')
            # results_for_csv.add_subsetting(Setting("Time since previous (s)"), 'time_since_previous')
            # results_for_csv.add_subsetting(Setting(f"Concentration\n<{top_nm}nm\n(counts/mL)"), 'total_conc_under_topnm')
            # results_for_csv.add_subsetting(Setting("Concentration\n(counts/mL)"), 'total_conc')
                
            time_of_above = None
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
                    # columns = list(settings.column_numbers.items())
                    columns = list(table_settings['columns_as_Settings_object'].column_numbers.items())
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
                            # print(group.subsettings)
                            # print(group.tag)
                            row.append(group.format_string.format(**{setting.tag: setting.get_value(sample) for setting in grouped_settings}))
                        elif column[0].tag.startswith('RESULTS'):
                            group = column[0]
                            grouped_settings = group.subsettings.values()
                            values = group.value_callback(sample)
                            if group.format_callback is not None:
                                row.append(group.format_callback(*values))
                                continue
                            row.append(group.format_string.format(**{setting.tag: value for setting, value in zip(grouped_settings, values)}))
                        else:
                            setting = column[0]
                            if setting.value_callback is None:
                                value = setting.get_value(sample)
                            else:
                                value = setting.value_callback(sample)
                            if setting.format_callback is None:
                                row.append(setting.format_string.format(**{setting.tag: value}))
                            else:
                                row.append(setting.format_callback(value))

                
                # exposure = settings.by_tag('Exposure').get_value(sample)
                # gain = settings.by_tag('Gain').get_value(sample)
                # row.append(f"{exposure} ms,\n{gain} dB")
                # detection_mode = settings.by_tag('DetectionThresholdType').get_value(sample)
                # detection_threshold = settings.by_tag('DetectionThreshold').get_value(sample)
                # if detection_threshold is None:
                #     row.append(detection_mode)
                # else:
                #     row.append(f"{detection_mode}\n{detection_threshold}")
                # framerate = settings.by_tag('FrameRate').get_value(sample)
                # frames_per_video = settings.by_tag('FramesPerVideo').get_value(sample)
                # video_duration = frames_per_video / framerate
                # if video_duration.is_integer():
                #     video_duration = int(video_duration)        
                # num_of_videos = settings.by_tag('NumOfVideos').get_value(sample)
                # row.append(f"{video_duration}x{num_of_videos}")
                # stir_time = settings.by_tag('StirredTime').get_value(sample)
                # stir_rpm = settings.by_tag('StirrerSpeed').get_value(sample)
                # row.append(f"{stir_time}x{stir_rpm}")
                # ID = settings.by_tag('ID').get_value(sample)
                # row.append('\n'.join((ID[0:4], ID[4:8], ID[8:12])))
                
                # previous = settings.by_tag('previous').get_value(sample)
                # results_for_csv.previous.set_value(sample, previous)
                # ID_of_previous = None
                # time = settings.by_tag('time').get_value(sample)
                # time_since_previous = None
                # if previous is not None:
                #     if previous not in unordered_samples:
                #         time_since_previous = '?'
                #     else:
                #         previous_sample = unordered_samples[previous]
                #         ID_of_previous = settings.by_tag('ID').get_value(previous_sample)
                #         time_of_previous = settings.by_tag('time').get_value(previous_sample)
                #         time_since_previous = int((time - time_of_previous).total_seconds())
                # results_for_csv.time_since_previous.set_value(sample, time_since_previous)
                # # # if ID_of_previous is not None:
                # # #     ID_of_previous = '\n'.join((ID_of_previous[0:4], ID_of_previous[4:8], ID_of_previous[8:12]))
                # # # row.append(ID_of_previous)
                # # if time_enabled:
                # #     text = []
                # #     time = settings.by_tag('time').get_value(sample)
                # #     time_since_above = None
                # #     if time_of_above is not None:
                # #         time_since_above = int((time - time_of_above).total_seconds())
                # #         text.append(f"{time_since_above} since above")
                # #     time_of_above = time

                # #     if previous is not None:
                # #         text.append(f"{time_since_previous} since previous")

                # #     # column_names.append(result_names['time'])
                # #     row.append('\n'.join(text))
                # #     text.clear()
                
                # data_sums = sums[i]
                # results_for_csv.total_conc.set_value(sample, f"{data_sums[1][1]:.2E}")
                # results_for_csv.total_conc_under_topnm.set_value(sample, f"{data_sums[0][1]:.2E}")
                # # if concentration_enabled:
                # #     text.append(f"Total: {data_sums[1][1]:.2E}")
                # #     text.append(f"<{top_nm}nm: {data_sums[0][1]:.2E}")
                # #     row.append('\n'.join(text))
                row.append("")
                yield row
        self.rows = tuple(generate_rows())
        self.need_reprep_tabulation = False
    def compare(self):
        assert self.need_recompute == False, "Must run NTA.compute() first."
        assert self.need_reprep_tabulation == False, "Must run NTA.prepare_tabulation() first."
        compare_info(self.settings, self.samples, self.results_for_csv, self.output_folder)
    def plot(self, grid_color = '0.8', name = 'Ridgeline plot'):
        assert self.need_recompute == False, "Must run NTA.compute() first."
        num_of_plots, samples, colors, table_settings, peak_settings, overall_min, overall_max, output_folder = self.num_of_plots, self.samples, self.colors, self.table_settings, self.peak_settings, self.overall_min, self.overall_max, self.output_folder
        peaks_enabled = (peak_settings is not None)
        table_enabled = (table_settings is not None)
        if table_enabled:
            assert self.need_reprep_tabulation == False, "Must run NTA.prepare_tabulation() first."
        cumulative_enabled, difference_enabled = self.cumulative_enabled, self.difference_enabled
        (_, height) = self.figsize
        tmp_filenames = self.tmp_filenames
        all_bins, all_sizes = np.load(tmp_filenames['bins']+'.npy'), np.load(tmp_filenames['sizes']+'.npy')
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
            # print(len(bins))
            # print(bins)
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
            
            plt.xlim(0, x_lim)
            plt.ylim(overall_min, overall_max)
            ax.patch.set_alpha(0)
                
            if i == final_i:
                ax.yaxis.get_offset_text().set_x(-0.1)
                plt.xlabel("Diameter (nm)")
                plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))
                ax.spines.left.set_visible(True)
                # plt.axhline(0, color = 'black')
                # break
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
            draw_table(fig, ax, self.rows, edges, table_settings, grid_color)

        fig.savefig(f"{output_folder}/{name}.png", dpi = 300, bbox_inches='tight')