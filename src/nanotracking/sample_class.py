#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 17:05:04 2023

@author: henryingels
"""

import os
# import pandas as pd
import numpy as np

class Sample():
    def __init__(self, folder, prefix, suffix, videos_file_prefix = None):
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
                if videos_file_prefix is not None and filename.startswith(videos_file_prefix) and filename.endswith(suffix):
                    # df = pd.read_csv(full_path, skiprows = 5, sep = '\t', engine = 'python')
                    # print(df.columns)
                    rows = [[]]
                    with open(full_path) as datfile:
                        i = 1
                        for line in datfile.readlines():
                            entries = line.split()
                            if len(entries) == 0 or entries[0].isdigit() is False: continue
                            entries = np.array([float(entry) if entry != 'nan' else np.nan for entry in entries], dtype = object)
                            # first_entry = float(entries[0])
                            first_entry = entries[0]
                            if first_entry.is_integer():
                                first_entry = int(first_entry)
                                if first_entry == i + 1:
                                    i += 1
                                    rows.append([])
                            rows[-1].extend(entries)
                            # print(line[0], len(line.split()))
                            # if np.almost_equal()    # Last 2 entries are avg and sd?
                    # print(rows)
                    # print(len(rows))
                    particles = []
                    for row in rows:
                        # print(f"Is {np.average([float(thing) for thing in row[4:]])} equal to {row[2]}?")
                        sizes = row[4:]
                        particles.extend(sizes)
                        average = np.average(sizes) if len(sizes) != 0 else 0
                        standard_deviation = np.std(sizes, ddof = 1) if len(sizes)-1 > 0 else 0
                        assert np.isclose(average, row[2]) or np.isnan(row[2])
                        assert np.isclose(standard_deviation, row[3]) or np.isnan(row[3])
                    self.particles = np.array(particles)
                    self.videos = [row[4:] for row in rows]

        self.xml = xml_path
        self.dat = dat_path
        self.info = info_path
        
        filename = os.path.basename(folder).removeprefix(prefix).removesuffix(suffix)
        self.filename = filename
        if hasattr(self, 'name') is False:
            self.name = filename