#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 17:02:47 2023

@author: henryingels
"""

import os
import pandas as pd
import numpy as np

from .settings_classes import Setting


def compare_info(settings, samples, results_object, output_folder):
    blank = Setting('')
    
    def generate_setting_objects():
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
        yield 'RESULTS:', blank
        for result in results_object.tags.values():
            name = result.name.replace('\n', ' ')
            units = setting.units
            if units != '': name += f" ({units})"
            yield name, result
    all_names, setting_objects = zip(*generate_setting_objects())
    
    same_valued_settings = []
    different_valued_settings = []
    for name, setting in zip(all_names, setting_objects):
        if setting is blank:
            same_valued_settings.append((name, setting))
            different_valued_settings.append((name, setting))
            continue
        sample_values = np.array([setting.get_value(sample) for sample in samples], dtype = object)
        are_same = np.all(sample_values == sample_values[0])
        if are_same:
            same_valued_settings.append((name, setting))
        else:
            different_valued_settings.append((name, setting))    
    
    all_csv_dataframe = pd.DataFrame(
        data = (
            pd.Series((setting.get_value(sample) for sample in samples), index = [sample.filename for sample in samples])
            for setting in setting_objects
        ), index = all_names
    )
    all_csv_dataframe.to_csv(os.path.join(output_folder, 'all.csv'))
    
    names_of_same, settings_of_same = zip(*same_valued_settings)
    same_values_csv_dataframe = pd.DataFrame(
        data = (
            pd.Series((setting.get_value(sample) for sample in samples), index = [sample.filename for sample in samples])
            for setting in settings_of_same
        ), index = names_of_same
    )
    same_values_csv_dataframe.to_csv(os.path.join(output_folder, 'same_values.csv'))
    
    names_of_different, settings_of_different = zip(*different_valued_settings)
    different_values_csv_dataframe = pd.DataFrame(
        data = (
            pd.Series((setting.get_value(sample) for sample in samples), index = [sample.filename for sample in samples])
            for setting in settings_of_different
        ), index = names_of_different
    )
    different_values_csv_dataframe.to_csv(os.path.join(output_folder, 'different_values.csv'))
