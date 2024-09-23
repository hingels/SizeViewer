import matplotlib as mpl
from .settings_classes import Setting


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