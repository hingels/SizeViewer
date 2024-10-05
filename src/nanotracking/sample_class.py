#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 17:05:04 2023

@author: henryingels
"""

import os
import numpy as np

class Sample():
    def __init__(self, folder, prefix, suffix, videos_file_prefix = None, index = None):
        self.folder, self.index = folder, index
        info_path, xml_path, dat_path = None, None, None
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
                if videos_file_prefix is not None and filename.startswith(videos_file_prefix) and filename.endswith(suffix):
                    rows = [[]]
                    with open(full_path) as datfile:
                        i = 1
                        for line in datfile.readlines():
                            entries = line.split()
                            if len(entries) == 0 or entries[0].isdigit() is False: continue
                            entries = np.array([float(entry) if entry != 'nan' else np.nan for entry in entries], dtype = object)
                            first_entry = entries[0]
                            if first_entry.is_integer():
                                first_entry = int(first_entry)
                                if first_entry == i + 1:
                                    i += 1
                                    rows.append([])
                            rows[-1].extend(entries)
                    particles = []
                    for row in rows:
                        sizes = row[4:]
                        particles.extend(sizes)
                        average = np.average(sizes) if len(sizes) != 0 else 0
                        standard_deviation = np.std(sizes, ddof = 1) if len(sizes)-1 > 0 else 0
                        assert np.isclose(average, row[2]) or np.isnan(row[2])
                        assert np.isclose(standard_deviation, row[3]) or np.isnan(row[3])
                    self.particles = np.array(particles)
                    self.videos = [row[4:] for row in rows]
        
        filename = os.path.basename(folder).removeprefix(prefix).removesuffix(suffix)
        if hasattr(self, 'name') is False:
            self.name = filename
        self.filename, self.xml, self.dat, self.info = filename, xml_path, dat_path, info_path