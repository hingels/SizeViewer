import matplotlib as mpl
from .settings_classes import Setting


def draw_table(fig, ax, settings, previous_setting, samples, unordered_samples, sums, results_for_table, edges, table_settings, grid_color):
    right_edge_figure = edges['right']
    table_bottom = edges['bottom']
    table_top = edges['top']
    
    include_experimental_unit = table_settings['include_experimental_unit']
    treatments_and_waits = table_settings['treatments_and_waits']
    results_column_names = table_settings['results_column_names']
    column_names = table_settings['column_names']
    column_widths = table_settings['column_widths']
    table_width = table_settings['width']
    margin_minimum_right = table_settings['margin_minimum_right']
    margin_left = table_settings['margin_left']

    transFigure = fig.transFigure

    include_treatments = ('_treatments_waits' in column_names)

    def generate_rows():
        num_of_plots = len(samples)
        
        column_quantities = dict()
        def number_of_subtags(tag):
            if (setting := settings.by_tag(tag)) is None:
                return 0
            return max(len(setting.subsettings), 1)
        def get_multivalued(tag, sample):
            if (setting := settings.by_tag(tag)) is None:
                return []
            if len(setting.subsettings) == 0:
                value = setting.get_value(sample)
                if value is None: return []
                return [value]
            
            subsettings = list(setting.numbered_subsettings.items())
            subsettings.sort()
            values = [subsetting.get_value(sample) for _, subsetting in subsettings]
            if values[0] is None:
                values[0] = setting.get_value(sample)
            return values
                
            
        top_nm = None
        for i in range(num_of_plots):
            sample = samples[i]
            
            settings.read_files(sample)
            settings.parse_time(sample)
            
            if include_treatments:
                for tag in ('treatment', 'wait'):
                    quantity = number_of_subtags(tag)
                    if tag not in column_quantities:
                        column_quantities[tag] = quantity
                        continue
                    column_quantities[tag] = max(column_quantities[tag], quantity)
            
            data_sums = sums[i]
            assert len(data_sums) == 2
            if top_nm is None:
                top_nm, _ = data_sums[1]
            assert data_sums[1][0] == top_nm
        
        for i, name in enumerate(column_names):
            if '{top_nm}' in name:
                column_names[i] = name.format(top_nm = top_nm)
            elif name == '_treatments_waits':
                column_names.pop(i)
                num_of_treatments = column_quantities['treatment']
                num_of_waits = column_quantities['wait']
                treatment_column_name, treatment_column_width = treatments_and_waits[0]
                wait_column_name, wait_column_width = treatments_and_waits[1]
                index = 0
                for j in range(max(num_of_treatments, num_of_waits)):
                    if j < num_of_treatments:
                        column_names.insert(i + index, treatment_column_name.format(treatment_number = j + 1))
                        column_widths.insert(i + index, treatment_column_width)
                        index += 1
                    if j < num_of_waits:
                        column_names.insert(i + index, wait_column_name.format(wait_number = j + 1))
                        column_widths.insert(i + index, wait_column_width)
                        index += 1

        for i, name in enumerate(results_column_names):     # This is redundant; should find a better way.
            if '{top_nm}' in name:
                results_column_names[i] = name.format(top_nm = top_nm)
        
        results_for_table.add_subsetting(previous_setting, 'previous')
        results_for_table.add_subsetting(Setting("Time since previous (s)"), 'time_since_previous')
        results_for_table.add_subsetting(Setting(f"Concentration\n<{top_nm}nm\n(counts/mL)"), 'total_conc_under_topnm')
        results_for_table.add_subsetting(Setting("Concentration\n(counts/mL)"), 'total_conc')
            
            
        time_of_above = None
        for i in range(num_of_plots):
            row = []
            
            sample = samples[i]
            
            if include_treatments:
                treatments = get_multivalued('treatment', sample)
                waits = get_multivalued('wait', sample)
                for j in range( max(column_quantities['treatment'], column_quantities['wait']) ):
                    if j < len(treatments): row.append(treatments[j])
                    elif j < column_quantities['treatment']: row.append(None)
                    if j < len(waits): row.append(waits[j])
                    elif j < column_quantities['wait']: row.append(None)
            
            if include_experimental_unit:
                experimental_unit = settings.by_tag('experimental_unit')
                text = ''
                if experimental_unit is not None:
                    value = experimental_unit.get_value(sample)
                    text += value if value is not None else ''
                    if hasattr(experimental_unit, 'age'):
                        age = experimental_unit.age.get_value(sample)
                        text += f"\n{age:.1f} d old" if age is not None else ''
                row.append(text)        
                        
            columns = list(settings.columns.items())
            columns.sort()
            for j, column in columns:
                content = '\n'.join(
                    setting.show_name*f"{setting.short_name}: " + f"{setting.get_value(sample)}" + setting.show_unit*f" ({setting.units})"
                    for setting in column if setting.get_value(sample) is not None )
                row.append(content)
            
            exposure = settings.by_tag('Exposure').get_value(sample)
            gain = settings.by_tag('Gain').get_value(sample)
            row.append(f"{exposure} ms,\n{gain} dB")
            
            detection_mode = settings.by_tag('DetectionThresholdType').get_value(sample)
            detection_threshold = settings.by_tag('DetectionThreshold').get_value(sample)
            if detection_threshold is None:
                row.append(detection_mode)
            else:
                row.append(f"{detection_mode}\n{detection_threshold}")
            
            framerate = settings.by_tag('FrameRate').get_value(sample)
            frames_per_video = settings.by_tag('FramesPerVideo').get_value(sample)
            video_duration = frames_per_video / framerate
            if video_duration.is_integer():
                video_duration = int(video_duration)        
            num_of_videos = settings.by_tag('NumOfVideos').get_value(sample)
            row.append(f"{video_duration}x{num_of_videos}")
            
            stir_time = settings.by_tag('StirredTime').get_value(sample)
            stir_rpm = settings.by_tag('StirrerSpeed').get_value(sample)
            row.append(f"{stir_time}x{stir_rpm}")
            
            ID = settings.by_tag('ID').get_value(sample)
            row.append('\n'.join((ID[0:4], ID[4:8], ID[8:12])))
            
            text = []
            
            time = settings.by_tag('time').get_value(sample)
            time_since_above = None
            if time_of_above is not None:
                time_since_above = int((time - time_of_above).total_seconds())
                text.append(f"{time_since_above} since above")
            time_of_above = time

            previous = settings.by_tag('previous').get_value(sample)
            results_for_table.previous.set_value(sample, previous)
            ID_of_previous = None
            time_since_previous = None
            if previous is not None:
                if previous not in unordered_samples:
                    time_since_previous = '?'
                else:
                    previous_sample = unordered_samples[previous]
                    ID_of_previous = settings.by_tag('ID').get_value(previous_sample)
                    time_of_previous = settings.by_tag('time').get_value(previous_sample)
                    time_since_previous = int((time - time_of_previous).total_seconds())
                text.append(f"{time_since_previous} since previous")
            results_for_table.time_since_previous.set_value(sample, time_since_previous)

            if ID_of_previous is not None:
                ID_of_previous = '\n'.join((ID_of_previous[0:4], ID_of_previous[4:8], ID_of_previous[8:12]))
            row.append(ID_of_previous)
            row.append('\n'.join(text))
            text.clear()
            
            data_sums = sums[i]
            
            text.append(f"Total: {data_sums[1][1]:.2E}")
            results_for_table.total_conc.set_value(sample, f"{data_sums[1][1]:.2E}")
            text.append(f"<{top_nm}nm: {data_sums[0][1]:.2E}")
            results_for_table.total_conc_under_topnm.set_value(sample, f"{data_sums[0][1]:.2E}")
            row.append('\n'.join(text))
            
            row.append("")
            
            yield row
    
    width_sum = sum([col_width for name, col_width in zip(column_names, column_widths) if name != ''])
    margin_right = table_width - width_sum
    assert margin_right >= margin_minimum_right, f"margin_right = {margin_right} < margin_minimum_right = {margin_minimum_right}."
    column_widths.append(margin_right)
    column_names.append("")
    # display_coords = final_ax.transData.transform([0, overall_min])
    edge = right_edge_figure + margin_left
    # table = plt.table(
    table = ax.table(
        tuple(generate_rows()),
        bbox = mpl.transforms.Bbox([[edge, table_bottom], [edge + table_width, table_top]]),
        transform = transFigure,
        cellLoc = 'left', colWidths = column_widths)
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    fig.add_artist(table)
    for i, name in enumerate(column_names):
        new_cell = table.add_cell(-1, i, width = column_widths[i], height = 0.1, text = name, loc = 'left')
        new_cell.set_text_props(fontweight = 'bold')
    final_column = len(column_widths) - 1
    for (row, column), cell in table.get_celld().items():
        if column == final_column:
            cell.set(edgecolor = None)
            continue
        cell.set(edgecolor = grid_color)