#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 17:02:47 2023

@author: henryingels
"""

import os
import pandas as pd
import numpy as np

from .settings_classes import Calculation, Setting


def compare_info(settings, samples, calculations, output_folder):
    blank = Setting('')
    
    def generate_entries():
        for tag, setting in settings.tags.items():
            if setting.tag == 'Name': continue
            if setting.hidden is False:
                units = setting.units
                name = setting.name if units == '' else f"{setting.name} ({units})"
                yield name, setting
            for subtag, subsetting in setting.subsettings.items():
                if subsetting.hidden: continue
                units = subsetting.units
                name = subsetting.name if units == '' else f"{subsetting.name} ({units})"
                yield name, subsetting
        yield '', blank
        yield 'CALCULATIONS:', blank
        for name, calculation in calculations.items():
            name = name.replace('\n', ' ')
            units = calculation.units
            if units != '': name += f" ({units})"
            yield name, calculation
    entry_names, entries = zip(*generate_entries())
    def get_values():
        for entry_name, entry in zip(entry_names, entries):
            if type(entry) is Setting:
                yield entry_name, entry, [sample.filename for sample in samples], [entry.get_value(sample) for sample in samples]
                continue
            for format_name in entry.formats.keys():
                yield entry_name, entry, [sample.filename for sample in samples], [entry.apply_format(format_name, sample) for sample in samples]
    entry_names, entries, filenames, values = zip(*get_values())

    same_valued_entries, different_valued_entries = [], []
    for entry_name, entry, sample_filenames, sample_values in zip(entry_names, entries, filenames, values):
        if entry is blank:
            same_valued_entries.append((entry_name, sample_filenames, sample_values))
            different_valued_entries.append((entry_name, sample_filenames, sample_values))
            continue
        sample_values = np.array(sample_values, dtype = object)
        are_same = np.all(sample_values == sample_values[0])
        if are_same:
            same_valued_entries.append((entry_name, sample_filenames, sample_values))
        else:
            different_valued_entries.append((entry_name, sample_filenames, sample_values))
    
    all_csv_dataframe = pd.DataFrame(
        data = (
            pd.Series(sample_values, index = sample_filenames)
            for sample_filenames, sample_values in zip(filenames, values)
        ), index = entry_names
    )
    all_csv_dataframe.to_csv(os.path.join(output_folder, 'all.csv'))
    
    names_of_same, filenames_of_same, values_of_same = zip(*same_valued_entries)
    same_values_csv_dataframe = pd.DataFrame(
        data = (
            pd.Series(sample_values_of_same, index = sample_filenames_of_same)
            for sample_filenames_of_same, sample_values_of_same in zip(filenames_of_same, values_of_same)
        ), index = names_of_same
    )
    same_values_csv_dataframe.to_csv(os.path.join(output_folder, 'same_values.csv'))
    
    names_of_different, filenames_of_different, values_of_different = zip(*different_valued_entries)
    different_values_csv_dataframe = pd.DataFrame(
        data = (
            pd.Series(sample_values_of_different, index = sample_filenames_of_different)
            for sample_filenames_of_different, sample_values_of_different in zip(filenames_of_different, values_of_different)
        ), index = names_of_different
    )
    different_values_csv_dataframe.to_csv(os.path.join(output_folder, 'different_values.csv'))
