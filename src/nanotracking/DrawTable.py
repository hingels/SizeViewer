import matplotlib as mpl
from copy import deepcopy
from .settings_classes import Setting, Settings


class Table():
    def __init__(self,
            nta_obj,
            width, margin_minimum_right, margin_left,
            include_experimental_unit=False,
            treatments_and_waits=None,
            columns_as_Settings_object=None,
            column_names=None,
            column_widths=None,
            column_names_without_treatmentsOrWaits=None,
            column_widths_without_treatmentsOrWaits=None):
        if columns_as_Settings_object is None:
            columns_as_Settings_object = Settings()
        if column_names is None: column_names = []
        if column_widths is None: column_widths = []
        if column_names_without_treatmentsOrWaits is None: column_names_without_treatmentsOrWaits = []
        if column_widths_without_treatmentsOrWaits is None: column_widths_without_treatmentsOrWaits = []
        self.nta_obj = nta_obj
        self.width, self.margin_minimum_right, self.margin_left = width, margin_minimum_right, margin_left
        self.include_experimental_unit, self.treatments_and_waits = include_experimental_unit, treatments_and_waits
        self.columns_as_Settings_object, self.column_names, self.column_widths, self.column_names_without_treatmentsOrWaits, self.column_widths_without_treatmentsOrWaits = columns_as_Settings_object, column_names, column_widths, column_names_without_treatmentsOrWaits, column_widths_without_treatmentsOrWaits
        self.fig, self.ax, self.table_plot = None, None, None
    def table_add_setting(self, setting: Setting):
        tag = setting.tag
        settings = self.columns_as_Settings_object
        assert tag not in settings.tags, f'Setting with tag "{tag}" already added to table.'
        if setting.column_number is None:
            setting.column_number = len(settings.column_widths)
        settings.add_setting(setting.tag, setting)
    def table_add_settings_by_tag(self, *tags, column_number = None, column_name = None, column_width = None, format_string = None, format_callback = None):
        '''
        Adds multiple Setting objects to the table.
        Example use case: specify column_number to group all specified settings into one column.

        If column_number is not given, the next available column will be used.
        Note: the column_number values of the specified Setting objects will be overwritten!
        
        If neither format_string nor format_callback are given, then for each cell in the column, the settings' individual format_strings will be used on separate lines.
        To use format_string, reference settings' values using their tags in curly braces: for example, format_string = "Red has power {RedLaserPower}."
        To use format_callback, define a function that accepts settings' values as arguments and returns a formatted (value-containing) string.

        If format_callback is given, it will be used instead of format_string.
        '''
        if column_number is None:
            column_number = len(self.columns_as_Settings_object.column_widths)
        get_setting_or_result = self.nta_obj.get_setting_or_result
        settings = [get_setting_or_result(tag) for tag in tags]
        if format_string is None:
            format_string = '\n'.join([setting.format_string for setting in settings])
        def prepare_setting(setting):
            setting.column_number = column_number
            if column_name is not None: setting.column_name = column_name
            if column_width is not None: setting.column_width = column_width
        if len(settings) == 1:
            setting = settings[0]
            prepare_setting(setting)
            setting.set_attributes(format_string = format_string, format_callback = format_callback)
            self.table_add_setting(setting)
            return
        if format_callback is None:
            group_suffix = format_string  # Allows multiple different format_callbacks or format_strings to be used on the same group, without counting as the same group (which would cause an error)
        else:
            group_suffix = format_callback.__name__
        group = Setting('COLUMN_' + '_'.join(tags) + group_suffix, column_number = column_number, column_name = column_name, column_width = column_width, format_string = format_string, format_callback = format_callback)
        for setting in settings:
            prepare_setting(setting)
            group.add_subsetting(setting, setting.tag)
        self.table_add_setting(group)
    def table_add_results(self, results_group, column_number = None, column_name = None, column_width = None, format_string = None, format_callback = None):
        if format_callback is None:
            group_suffix = format_string  # Allows multiple different format_callbacks or format_strings to be used on the same group, without counting as the same group (which would cause an error)
        else:
            group_suffix = format_callback.__name__
        new_column = deepcopy(results_group)
        new_column.set_attributes(tag = results_group.tag + group_suffix, column_number = column_number, column_name = column_name, column_width = column_width, format_string = format_string, format_callback = format_callback)
        self.table_add_setting(new_column)
    def table_add_experimental_unit(self, column_name = "Experimental\nunit", width = 0.3, column_number = None):
        '''
        Adds to the table a column for experimental unit, whose name is given by "experimental_unit=…" in each sample's info.md file.
        '''
        experimental_unit = self.nta_obj.settings.by_tag('experimental_unit')
        if column_number is None:
            column_number = len(self.columns_as_Settings_object.column_widths)
        experimental_unit.column_number = column_number
        experimental_unit.column_name = column_name
        experimental_unit.column_width = width
        self.table_add_setting(experimental_unit)
    def table_add_treatments_and_waits(self, treatments_column_name, treatments_width, waits_column_name, waits_width):
        '''
        For each treatment & wait-time listed in samples' info.md files, adds to the table
        (1) a column for the treatment's name, and (2) a column for the time waited after applying the treatment.
        '''
        start_index = len(self.column_names_without_treatmentsOrWaits)
        assert self.treatments_and_waits is None, "Treatments and waits have already been added to the table."
        self.treatments_and_waits = [start_index, (treatments_column_name, treatments_width), (waits_column_name, waits_width)]
    def reset_columns(self):
        self.column_names = list(self.columns_as_Settings_object.column_names.keys())
        self.column_widths = self.columns_as_Settings_object.column_widths.copy()

    def draw_table(self, fig, ax, rows, edges, grid_color):
        right_edge_figure = edges['right']
        table_bottom = edges['bottom']
        table_top = edges['top']
        column_names = self.column_names
        column_widths = self.column_widths
        table_width = self.width
        margin_minimum_right = self.margin_minimum_right
        margin_left = self.margin_left
        transFigure = fig.transFigure
        
        width_sum = sum([col_width for name, col_width in zip(column_names, column_widths) if name != ''])
        margin_right = table_width - width_sum
        assert margin_right >= margin_minimum_right, f"margin_right = {margin_right} < margin_minimum_right = {margin_minimum_right}. Try increasing the table's \"width\" setting."
        column_widths.append(margin_right)
        column_names.append("")
        # display_coords = final_ax.transData.transform([0, overall_min])
        edge = right_edge_figure + margin_left
        table = ax.table(
            rows,
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
        self.fig, self.ax, self.table_plot = fig, ax, table
        return table