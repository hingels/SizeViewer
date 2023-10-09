#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 17:02:47 2023

@author: henryingels
"""

import os
import pandas as pd
import numpy as np

from settings_classes import Setting


def compare_info(settings, samples, results_object, output_folder):
    os.makedirs(output_folder, exist_ok = True)
    
    blank = Setting('')
    
    def generate_setting_objects():
        for tag, setting in settings.tags.items():
            if setting.hidden is False:
                yield tag, setting
            for subtag, subsetting in setting.subsettings.items():
                if subsetting.hidden: continue
                yield f"{tag}.{subtag}", subsetting
        yield '', blank
        yield 'RESULTS:', blank
        for result in results_object.subsettings.values():
            name = result.name.replace('\n', ' ')
            yield name, result
    all_tags, setting_objects = zip(*generate_setting_objects())
    
    same_valued_settings = []
    different_valued_settings = []
    for setting in setting_objects:
        sample_values = np.array([setting.get_value(sample) for sample in samples], dtype = object)
        are_same = np.all(sample_values == sample_values[0])
        if are_same:
            same_valued_settings.append(setting)
        else:
            different_valued_settings.append(setting)
    
    all_csv_dataframe = pd.DataFrame(
        data = (
            pd.Series((setting.get_value(sample) for sample in samples), index = [sample.filename for sample in samples])
            for setting in setting_objects
        ), index = all_tags
    )
    all_csv_dataframe.to_csv(os.path.join(output_folder, 'all.csv'))
    
    same_values_csv_dataframe = pd.DataFrame(
        data = (
            pd.Series((setting.get_value(sample) for sample in samples), index = [sample.filename for sample in samples])
            for setting in same_valued_settings
        ), index = [setting.tag for setting in same_valued_settings]
    )
    same_values_csv_dataframe.to_csv(os.path.join(output_folder, 'same_values.csv'))
    
    different_values_csv_dataframe = pd.DataFrame(
        data = (
            pd.Series((setting.get_value(sample) for sample in samples), index = [sample.filename for sample in samples])
            for setting in different_valued_settings
        ), index = [setting.tag for setting in different_valued_settings]
    )
    different_values_csv_dataframe.to_csv(os.path.join(output_folder, 'different_values.csv'))
