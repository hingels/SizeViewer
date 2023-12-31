#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 17:05:04 2023

@author: henryingels
"""

import os

class Sample():
    def __init__(self, folder, prefix, suffix):
        self.folder = folder
        info_path = None
        xml_path = None
        dat_path = None
        for (path, subdirs, files) in os.walk(folder):
            for filename in files:
                if filename.startswith('._'): continue
                full_path = os.path.join(path, filename)
                if filename == 'info.md':
                    info_path = full_path
                split_name = filename.split('_')
                if filename.endswith('.xml') and len(split_name) > 2:
                    truncated = '_'.join(split_name[:-2])
                    if truncated.endswith('Process') is False and truncated.endswith('Temperature') is False:
                        xml_path = full_path
                if filename.startswith(prefix) and filename.endswith(suffix):
                    dat_path = full_path
        self.xml = xml_path
        self.dat = dat_path
        self.info = info_path
        
        filename = os.path.basename(folder).removeprefix(prefix).removesuffix(suffix)
        self.filename = filename
        if hasattr(self, 'name') is False:
            self.name = filename