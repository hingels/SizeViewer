import matplotlib as mpl
from .settings_classes import Setting, Settings


class Table():
    def __init__(self,
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
        self.width, self.margin_minimum_right, self.margin_left = width, margin_minimum_right, margin_left
        self.include_experimental_unit, self.treatments_and_waits = include_experimental_unit, treatments_and_waits
        self.columns_as_Settings_object, self.column_names, self.column_widths, self.column_names_without_treatmentsOrWaits, self.column_widths_without_treatmentsOrWaits = columns_as_Settings_object, column_names, column_widths, column_names_without_treatmentsOrWaits, column_widths_without_treatmentsOrWaits
    def table_add_setting(self, setting: Setting):
        tag = setting.tag
        settings = self.columns_as_Settings_object
        assert tag not in settings.tags, f'Setting with tag "{tag}" already added to table.'
        if setting.column_number is None:
            setting.column_number = len(settings.column_widths)
        settings.add_setting(setting.tag, setting)

def draw_table(fig, ax, rows, edges, table_settings, grid_color):
    right_edge_figure = edges['right']
    table_bottom = edges['bottom']
    table_top = edges['top']
    column_names = table_settings['column_names']
    column_widths = table_settings['column_widths']
    table_width = table_settings['width']
    margin_minimum_right = table_settings['margin_minimum_right']
    margin_left = table_settings['margin_left']
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