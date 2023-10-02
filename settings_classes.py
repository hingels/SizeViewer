#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 17:08:22 2023

@author: henryingels
"""

from sample_class import Sample
import typing
import xml.etree.ElementTree as ET

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
        
        if value is None:
            converted_value = None
        elif datatype is bool:
            converted_value = (value.lower() == 'true')
        else:
            converted_value = datatype(value)
        
        self.sample_values[sample] = converted_value
    def get_value(self, sample: Sample):
        sample_values = self.sample_values
        if sample not in sample_values: return None
        return sample_values[sample]
    def get_subvalues(self, sample: Sample):
        value = self.get_value(sample)
        if value is None:
            return [None]
        if type(value) is list:
            value.sort()
            value = [subvalue for index, subvalue in value]
        else:
            value = [value]
        print(value)
        return value
class Settings():
    def __init__(self, settings_dict = None):
        if settings_dict is None:
            settings_dict = dict()
        self.tags = dict()
        self.columns = dict()
        for tag, setting in settings_dict.items():
            self.add_setting(tag, setting)
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
        if column not in columns:
            columns[column] = []
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
    def read_files(self, sample: Sample, get_all = False, dependencies: dict = None):
        if dependencies is None: dependencies = dict()
        tags = self.tags
        by_tag = self.by_tag
        add_setting = self.add_setting
        with open(sample.xml) as xml_file:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            for entry in (root.iter() if get_all else root.find('RecordingSettings')):
                tag = entry.tag
                if tag in tags:
                    setting = by_tag(tag)
                    setting.set_value(sample, entry.text)
                    continue
                if tag in dependencies:
                    dependency = dependencies[tag]
                    setting = Setting(tag, depends_on = dependency)
                else:
                    setting = Setting(tag)
                add_setting(tag, setting)
                setting.set_value(sample, entry.text)
        with open(sample.info) as info_file:
            for line in info_file.readlines():
                tag, value = line.split('=')
                value = value.strip()
                tag = tag.split('.')
                is_multivalued = (len(tag) != 1)
                tag_base = tag[0]
                if tag_base in tags:
                    setting = by_tag(tag_base)
                else:
                    if tag_base in dependencies:
                        dependency = dependencies[tag_base]
                        setting = Setting(tag_base, depends_on = dependency)
                    else:
                        setting = Setting(tag_base)
                    add_setting(tag_base, setting)
                if is_multivalued is False:
                    setting.set_value(sample, value)
                    continue
                assert len(tag) == 2
                value_index = tag[1]
                assert value_index.isdigit() and float(value_index).is_integer()
                value_index = int(value_index)
                if sample in setting.sample_values:
                    current_value = setting.get_value(sample)
                    if type(current_value) is not list:
                        current_value = [(1, current_value)]
                    current_value.append((value_index, value))
                else:
                    current_value = [(1, value)]
                setting.set_value(sample, current_value, datatype = list)
        self.apply_dependencies()