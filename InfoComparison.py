#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 17:02:47 2023

@author: henryingels
"""

import os
import pandas as pd
import numpy as np

from sample_class import Sample
from settings_classes import Setting, Settings


output_folder = "/Users/henryingels/Documents/GitHub/Ridgeline-Plotter/CSV outputs"
datafolder = "/Volumes/Lab Drive/ViewSizer 3000/Data"
prefix = 'ConstantBinsTable_'
suffix = '.dat'
use_filenames = True
filenames = [
    # "230709 1-100000 standard",
    # "230709 1-100000 standard #2",
    "230729 dilut L1-100 contr",
    "230729 dilut L1-100 contr #2",
    "230728 KW+EW",
    "230729 KW+EW re-measure",
    "230729 KW+EW re-measure 2",
    "230729 KW+EW re-measure 3",
    "230729 KW+EW re-measure 4",
    "230701a, pre+96h fluor +DiO+Triton"
]
table_width = 1.2
table_left_margin = 0
minimum_table_right_margin = 0.03
column_names = ["1st\ntreatment\n(µM)", "1st\n4°C\nwait\n(h)", "2nd\ntreatment\n(µM)", "2nd\n4°C\nwait\n(h)", "Experimental\nunit", "Filter", "Power\n(mW)", "Exposure\n(ms)", "Gain\n(dB)", "Video\nduration (s)\nx quantity"]
column_widths = [0.14, 0.07, 0.14, 0.07, 0.19, 0.08, 0.1, 0.13, 0.08, 0.16]


width_sum = sum(column_widths)
table_right_margin = table_width - width_sum
assert table_right_margin >= minimum_table_right_margin, f"table_right_margin = {table_right_margin} < minimum_table_right_margin = {minimum_table_right_margin}."
column_widths = np.append(column_widths, table_right_margin)
column_names.append("")
                    

def generate_samples():
    for folder in os.listdir(datafolder):
        sample = Sample(os.path.join(datafolder, folder), prefix, suffix)
        if use_filenames and sample.filename not in filenames: continue
        yield sample.filename, sample
if use_filenames:
    unordered_samples = dict(generate_samples())
    samples = [unordered_samples[name] for name in filenames]
else:
    _, samples = tuple(zip(*generate_samples()))

settings = Settings()

dependencies = {
    'RedLaserPower': 'RedLaserEnabled',
    'GreenLaserPower': 'GreenLaserEnabled',
    'BlueLaserPower': 'BlueLaserEnabled',
}
for tag, dependency_tag in dependencies.items():
    dependency = Setting(dependency_tag, datatype = bool)
    dependencies[tag] = dependency
    settings.add_setting(dependency_tag, dependency)

for sample in samples:
    settings.read_files(sample, get_all = True, dependencies = dependencies)


os.makedirs(output_folder, exist_ok = True)

def generate_setting_objects():
    for tag, setting in settings.tags.items():
        yield tag, setting
        for subtag, subsetting in setting.subsettings.items():
            yield f"{tag}.{subtag}", subsetting
# all_tags, setting_objects = zip(*settings.tags.items())
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
