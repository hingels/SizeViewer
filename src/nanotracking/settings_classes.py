#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 17:08:22 2023

@author: henryingels
"""

from .sample_class import Sample
import typing
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from collections import OrderedDict
import re

def format_string_to_function(format_string, *output_names):
    def format_function(*output_values):
        outputs = dict(zip(output_names, output_values))
        new_format_string = format_string
        for placeholder in re.findall(r'\{(?!\{).*?\}', format_string):
            no_braces = placeholder[1:-1]
            assert no_braces in outputs, f'Placeholder "{no_braces}" was given in format string, but no corresponding output name was given.'
            new_format_string = new_format_string.replace(placeholder, str(outputs[no_braces]))
        return new_format_string
    format_function.__name__ = str(int.from_bytes(format_string.encode(), 'little')) # For compatibility with the hacky line "group_suffix = format.__name__" in DrawTable.py
    return format_function
def format_apply_wordwrap(format_function, letters_per_line, no_hyphens = False):
    if letters_per_line is None: return format_function
    def new_format(*args):
        text = format_function(*args)
        if no_hyphens:
            words = text.split()
            new_text = [[]]
            line_length = 0
            for word in words:
                word_length = len(word)
                if line_length + word_length > letters_per_line:
                    new_text.append([])
                    line_length = 0
                new_text[-1].append(word)
                line_length += word_length
            return '\n'.join((' '.join(line) for line in new_text))
        loops = 0
        for i in range(letters_per_line, len(text), letters_per_line):
            before, after = text[:i+loops], text[i+loops:]
            if before[-1] != ' ' and no_hyphens == False:
                before += '-'
            text = before + '\n' + after
            loops += 1
        return text
    return new_format

class Calculation():
    def __init__(self, name, value_function, *output_names, units = None, samples = None):
        '''
        Defines a calculation with a name, a function determining its value (value_function), and how to name the function's outputs (output_names).
        The function value_function takes a Sample object as its only argument.

        Calculation.refresh() must be called immediately after a new Calculation is initialized.
        To skip this, use the keyword argument "samples" to specify which samples to use when refreshing.
        '''
        if units is None: units = ''
        self.name, self.value_function = name, value_function
        self.output_names, self.sample_output_values = output_names, dict()
        self.units, self.formats = units, dict()
        if samples is not None:
            self.refresh(*samples)
    def add_format(self, name, format):
        if type(format) is str:
            format = format_string_to_function(format, *self.output_names)
        self.formats[name] = format
    def apply_format(self, name, sample):
        output_values = self.sample_output_values[sample]
        format = self.formats[name]
        return format(*output_values)
    def refresh(self, *samples):
        '''
        For each sample specified in samples, recalculates output values using Calculation.value_function.
        '''
        for sample in samples:
            output_values = self.value_function(sample)
            try: iter(output_values)
            except TypeError: output_values = [output_values]
            self.sample_output_values[sample] = output_values
    def representation_as_setting(self, format_name, samples):
        '''
        Returns a Setting object whose subsettings represent the outputs of Calculation.value_function, including their numerical values.
        A new Setting object is created each time this runs!
        '''
        sample_output_values = self.sample_output_values
        settings_representation = Setting(f'CALC_{self.name}_FORMAT_{format_name}')
        for i, output_name in enumerate(self.output_names):
            output = Setting(output_name, datatype = None)
            for sample in samples:
                output.set_value(sample, sample_output_values[sample][i])
            settings_representation.add_subsetting(output_name, output)
        # settings_representation.value_function = self.value_function
        settings_representation.format = self.formats[format_name]
        return settings_representation

class Setting():
    def __init__(self, tag, short_name = None, format = None, value_function = None, name = None, units = '', column_number = None, column_name = None, column_width = None, sample_values: dict = None, show_unit = False, show_name = False, datatype = str, depends_on = None, subsettings = None, hidden = False, dependencies_require = True):
        if name is None: name = tag
        if short_name is None: short_name = name
        if format is None:
            def format_function(value):
                return show_name*f"{short_name}: " + str(value) + show_unit*f" ({units})"
            format = format_function
        if type(format) is str:
            format = format_string_to_function(format, tag)
        self.format = format
        self.tag, self.short_name = tag, short_name
        self.value_function, self.datatype = value_function, datatype
        self.name, self.show_name = name, show_name
        self.units, self.show_unit = units, show_unit
        self.column_number, self.column_name, self.column_width = column_number, column_name, column_width
        self.depends_on, self.dependencies_require = depends_on, dependencies_require
        self.hidden = hidden

        self.sample_values = dict()
        if sample_values is not None:
            for sample, value in sample_values.items():
                self.set_value(sample, value)

        self.subsettings, self.numbered_subsettings = OrderedDict(), OrderedDict()
        if subsettings is not None:
            for subtag, subsetting in subsettings.items():
                self.add_subsetting(subtag, subsetting)
        
    def add_subsetting(self, subtag, subsetting):
        self.subsettings[subtag] = subsetting
        if type(subtag) is int:
            self.numbered_subsettings[subtag] = subsetting
            return
        assert hasattr(self, subtag) is False
        setattr(self, subtag, subsetting)
    def set_value(self, sample: Sample, value, datatype = None):
        if datatype is None:
            datatype = self.datatype
            if datatype is None:
                self.sample_values[sample] = value
                return
        
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
        elif datatype is Sample:
            converted_value = value
        else:
            converted_value = datatype(value)
        
        self.sample_values[sample] = converted_value
    def get_value(self, sample: Sample):
        sample_values = self.sample_values
        if sample not in sample_values: return None
        return sample_values[sample]
    def set_attributes(self, **attrs):
        for name, value in attrs.items():
            self.__setattr__(name, value)
class Settings():
    def __init__(self, settings_dict = None):
        self.tags, self.column_numbers, self.column_names = OrderedDict(), OrderedDict(), OrderedDict()
        self.column_widths = []
        if settings_dict is None:
            settings_dict = OrderedDict()
            return
        for tag, setting in settings_dict.items():
            self.add_setting(tag, setting)
    def by_tag(self, tag):
        tags = self.tags
        if tag not in tags: return None
        return tags[tag]
    def add_setting(self, tag, setting):
        tags, column_numbers, column_names = self.tags, self.column_numbers, self.column_names
        assert tag not in tags
        tags[tag] = setting
        
        self.column_widths.append(setting.column_width)
        column_number, column_name = setting.column_number, setting.column_name
        if column_number is not None:
            if column_number not in column_numbers:
                column_numbers[column_number] = []
            column_numbers[column_number].append(setting)
        if column_name is not None:
            if column_name not in column_names:
                column_names[column_name] = []
            column_names[column_name].append(setting)
    def apply_dependencies(self):
        '''
        For any setting that is dependent on another setting, set it to zero (or equivalent) if the dependency has a value of False.
        '''
        for tag in self.tags:
            setting = self.by_tag(tag)
            depends_on = setting.depends_on
            if depends_on is None: continue
            requirement = depends_on.dependencies_require
            for sample, value in setting.sample_values.items():
                if depends_on.get_value(sample) != requirement:
                    setting.set_value(sample, None)
    def parse_time(self, sample):
        if 'MeasurementStartDateTime' not in self.tags: return
        if 'time' not in self.tags:
            time_setting = Setting('time', hidden = True)
            self.add_setting('time', time_setting)
        measurement_time = datetime.strptime(self.by_tag('MeasurementStartDateTime').get_value(sample), '%Y-%m-%d %H:%M:%S')
        self.by_tag('time').set_value(sample, measurement_time, datatype = datetime)
        
        if 'experimental_unit' not in self.tags: return
        experimental_unit = self.by_tag('experimental_unit')
        if not hasattr(experimental_unit, 'date'): return
        experimental_unit_date = datetime.strptime(experimental_unit.date.get_value(sample), '%Y/%m/%d')
        age = (measurement_time - experimental_unit_date).total_seconds() / 86400
        
        if not hasattr(experimental_unit, 'age'):
            age_subsetting = Setting('age', name = 'Age', units = 'days', datatype = float)
            experimental_unit.add_subsetting('age', age_subsetting)
        experimental_unit.age.set_value(sample, age)
    def read_files(self, sample: Sample, dependencies: dict = None):
        if dependencies is None: dependencies = dict()
        tags = self.tags
        by_tag, add_setting = self.by_tag, self.add_setting
        with open(sample.xml) as xml_file:
            tree = ET.parse(xml_file)
            root = tree.getroot()
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
                    units = setting.units
                else: units = ''
                if subtag not in setting.subsettings:
                    subsetting = Setting(full_tag, name = f"{setting.name}: {subtag}", units = units)
                    setting.add_subsetting(subtag, subsetting)
                else:
                    subsetting = setting.subsettings[subtag]
                
                subsetting.set_value(sample, value)
                
        self.apply_dependencies()
