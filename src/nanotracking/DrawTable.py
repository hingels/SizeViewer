import matplotlib as mpl
from copy import deepcopy
from .settings_classes import Setting, Settings, format_string_to_function, format_apply_wordwrap


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
    def add_setting(self, setting: Setting):
        tag = setting.tag
        settings = self.columns_as_Settings_object
        assert tag not in settings.tags, f'Setting with tag "{tag}" already added to table.'
        if setting.column_number is None:
            setting.column_number = len(settings.column_widths)
        settings.add_setting(setting.tag, setting)
    def add_settings_by_tag(self, *tags, column_number = None, column_name = None, column_width = None, format = None, letters_per_line = None, no_hyphens = False):
        '''
        Adds multiple Setting objects to the table.
        Example use case: specify column_number to group all specified settings into one column.

        If column_number is not given, the next available column will be used.
        Note: the column_number values of the specified Setting objects will be overwritten!
        
        There are two ways to specify format:
        1. As a format string, which should reference settings' values using their tags in curly braces. For example, format = "Red has power {RedLaserPower}."
        2. As a function, which should accept settings' values as arguments and return a string.
        If format is not given, then for each cell in the column, the settings' individual formats will be used on separate lines.

        If letters_per_line is specified, text will be split into lines with letters_per_line characters each.
        If additionally no_hyphens is True, hyphens will not be used for word wrap; instead, words will be grouped such that each line does not exceed letters_per_line characters.
        '''
        if column_number is None:
            column_number = len(self.columns_as_Settings_object.column_widths)
        def prepare_setting(setting):
            setting.column_number = column_number
            if column_name is not None: setting.column_name = column_name
            if column_width is not None: setting.column_width = column_width
        get_setting_or_calculation = self.nta_obj.get_setting_or_calculation
        settings = [get_setting_or_calculation(tag) for tag in tags]
        if type(format) is str:
            format = format_string_to_function(format, *tags)
        if len(settings) == 1:
            setting = settings[0]
            prepare_setting(setting)
            if format is not None:
                setting.set_attributes(format = format_apply_wordwrap(format, letters_per_line, no_hyphens = no_hyphens))
            self.add_setting(setting)
            return
        if format is None:
            def format_function(*output_values):
                assert len(output_values) == len(settings), f"Number of values given = {len(output_values)}; should be {len(settings)}"
                return '\n'.join([formatted for setting, value in zip(settings, output_values) if (formatted := setting.format(value)) != ''])
            format_function.__name__ = ''.join([setting.tag for setting in settings]) # For compatibility with hacky line "group_suffix = format.__name__" below
            format = format_function
        group_suffix = format.__name__ # TODO: replace with a less hacky solution
        group = Setting('COLUMN_' + '_'.join(tags) + group_suffix, column_number = column_number, column_name = column_name, column_width = column_width, format = format_apply_wordwrap(format, letters_per_line, no_hyphens = no_hyphens))
        for setting in settings:
            prepare_setting(setting)
            group.add_subsetting(setting.tag, setting)
        self.add_setting(group)
    def add_calculation(self, calculation, format_name, column_number = None, column_name = None, column_width = None):
        new_column = calculation.representation_as_setting(format_name, self.nta_obj.samples)
        new_column.set_attributes(column_number = column_number, column_name = column_name, column_width = column_width)
        self.add_setting(new_column)
    def add_experimental_unit(self, column_name = "Experimental\nunit", width = 0.3, column_number = None):
        '''
        Adds to the table a column for experimental unit, whose name is given by "experimental_unit=…" in each sample's info.md file.
        '''
        experimental_unit = self.nta_obj.settings.by_tag('experimental_unit')
        if column_number is None:
            column_number = len(self.columns_as_Settings_object.column_widths)
        experimental_unit.column_number = column_number
        experimental_unit.column_name = column_name
        experimental_unit.column_width = width
        self.add_setting(experimental_unit)
    def add_treatments_and_waits(self, treatments_column_name, treatments_width, waits_column_name, waits_width):
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
        right_edge_figure, table_bottom, table_top = edges['right'], edges['bottom'], edges['top']
        column_names, column_widths = self.column_names, self.column_widths
        table_width, margin_minimum_right, margin_left = self.width, self.margin_minimum_right, self.margin_left
        edge = right_edge_figure + margin_left
        
        width_sum = sum([col_width for name, col_width in zip(column_names, column_widths) if name != ''])
        margin_right = table_width - width_sum
        assert margin_right >= margin_minimum_right, f"margin_right = {margin_right} < margin_minimum_right = {margin_minimum_right}. Try increasing the table's \"width\" setting."
        column_names.append(""); column_widths.append(margin_right)
        table = ax.table(
            rows,
            bbox = mpl.transforms.Bbox([[edge, table_bottom], [edge + table_width, table_top]]),
            transform = fig.transFigure,
            cellLoc = 'left', colWidths = column_widths)
        table.auto_set_font_size(False); table.set_fontsize(12)
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