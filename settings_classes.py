#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 17:08:22 2023

@author: henryingels
"""

from sample_class import Sample
import typing
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta

class Setting():
    def __init__(self, tag, name = None, units = '', column = None, sample_values: dict = None, show_unit = False, show_name = False, datatype = str, depends_on = None, subsettings = None, hidden = False):
        self.tag = tag
        if name is None: name = tag
        self.name = name
        self.units = units
        self.column = column
        self.sample_values = dict()
        if sample_values is not None:
            for sample, value in sample_values.items():
                self.set_value(sample, value)

        self.subsettings = dict()
        self.numbered_subsettings = dict()
        if subsettings is not None:
            for subtag, subsetting in subsettings.items():
                self.add_subsetting(subtag, subsetting)
        
        self.show_unit = show_unit
        self.show_name = show_name
        self.datatype = datatype
        self.depends_on = depends_on

        self.hidden = hidden
    def add_subsetting(self, subsetting, subtag):
        self.subsettings[subtag] = subsetting
        if type(subtag) is int:
            self.numbered_subsettings[subtag] = subsetting
            return
        assert hasattr(self, subtag) is False
        setattr(self, subtag, subsetting)
    def set_value(self, sample: Sample, value, datatype = None):
        if datatype is None:
            datatype = self.datatype
        
        if value is None:
            converted_value = None
        elif datatype is bool:
            converted_value = (value.lower() == 'true')
        elif datatype is datetime:
            assert type(value) is datetime
            converted_value = value
        elif datatype is timedelta:
            assert type(value) is timedelta
            converted_value = value
        else:
            converted_value = datatype(value)
        
        self.sample_values[sample] = converted_value
    def get_value(self, sample: Sample):
        sample_values = self.sample_values
        if sample not in sample_values: return None
        return sample_values[sample]
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
    def parse_time(self, sample):
        if 'MeasurementStartDateTime' not in self.tags: return
        measurement_time = datetime.strptime(self.by_tag('MeasurementStartDateTime').get_value(sample), '%Y-%m-%d %H:%M:%S')
        # sample.set_time(measurement_time)
        if 'time' not in self.tags:
            time_setting = Setting('time', hidden = True)
            self.add_setting('time', time_setting)
        self.by_tag('time').set_value(sample, measurement_time, datatype = datetime)
        
        if 'experimental_unit' not in self.tags: return
        experimental_unit = self.by_tag('experimental_unit')
        if hasattr(experimental_unit, 'date') is False: return
        experimental_unit_date = datetime.strptime(experimental_unit.date.get_value(sample), '%Y/%m/%d')
        age = (measurement_time - experimental_unit_date).total_seconds() / 86400
        
        if hasattr(experimental_unit, 'age') is False:
            age_subsetting = Setting('age')
            experimental_unit.add_subsetting(age_subsetting, 'age')
        # sample.set_age(f"{age:.1f} days old")
        experimental_unit.age.set_value(sample, f"{age:.1f} d old")
    def read_files(self, sample: Sample, get_all = False, dependencies: dict = None):
        if dependencies is None: dependencies = dict()
        tags = self.tags
        by_tag = self.by_tag
        add_setting = self.add_setting
        with open(sample.xml) as xml_file:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            # for entry in (root.iter() if get_all else root.find('RecordingSettings')):
            for entry in root.iter():
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
                full_tag, value = line.split('=')
                value = value.strip()
                tag_split = full_tag.split('.')
                
                tag_base = tag_split[0]
                if tag_base not in tags:
                    if tag_base in dependencies:                        
                        dependency = dependencies[tag_base]
                        setting = Setting(tag_base, depends_on = dependency)
                    else:
                        setting = Setting(tag_base)
                    add_setting(tag_base, setting)
                else:
                    setting = by_tag(tag_base)
                
                if len(tag_split) == 1:        # If has no subvalues:
                    setting.set_value(sample, value)
                    continue
                
                assert len(tag_split) == 2
                subtag = tag_split[1]
                
                if subtag.isdigit():
                    assert float(subtag).is_integer()
                    subtag = int(subtag)
                
                if subtag not in setting.subsettings:
                    subsetting = Setting(full_tag)
                    setting.add_subsetting(subsetting, subtag)
                else:
                    subsetting = setting.subsettings[subtag]
                
                subsetting.set_value(sample, value)
                
        self.apply_dependencies()
