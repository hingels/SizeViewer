#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 17:08:22 2023

@author: henryingels
"""

from sample_class import Sample

class Setting():
    def __init__(self, tag, name = None, units = '', column = None, sample_values: dict = None, show_unit = False, show_name = False, datatype = str, depends_on = None):
        self.tag = tag
        if name is None: name = tag
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
    def set_value(self, sample: Sample, value, datatype = None):
        if datatype is None:
            datatype = self.datatype
        self.sample_values[sample] = datatype(value) if datatype is not bool else value.lower() == "true"
    def get_value(self, sample: Sample):
        sample_values = self.sample_values
        if sample not in sample_values: return None
        return sample_values[sample]
class Settings():
    def __init__(self, settings_dict = None):
        if settings_dict is None:
            settings_dict = dict()
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
        tags = self.tags
        if tag not in tags: return None
        return tags[tag]
    def add_setting(self, tag, setting):
        tags = self.tags
        assert tag not in tags
        tags[tag] = setting
        
        columns = self.columns
        column = setting.column
        if column is None: return
        if column > (len(columns) - 1):
            columns.append([])
        columns[column].append(setting)
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