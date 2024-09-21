import unittest
from src.nanotracking import DifferencePlotter

class Test_Table(unittest.TestCase):
    def setUp(self):
        nta = DifferencePlotter.NTA(
            datafolder = "tests/Test data",
            output_folder = "tests/Test output",
            filenames = ["1"]
        )
        nta.compute()
        results_column_names = ["Time (s)", "Concentration\n(counts/mL)"]
        self.table_settings = {
            'include_experimental_unit': False,
            'treatments_and_waits': [("Treatment\n{treatment_number}\n(µM)", 0.2), ("4°C\nwait\n{wait_number}\n(h)", 0.07)],
            'results_column_names': results_column_names,
            'column_names': ["_treatments_waits", "Filter\ncut-on\n(nm)", "Power\n(mW)", "Exposure,\ngain", "Detection\nsetting", "Video sec\nx quantity", "Stir sec\nx RPM", "ID", "ID of\nprevious", *results_column_names],
            'column_widths': [0.1, 0.19, 0.14, 0.19, 0.16, 0.12, 0.1, 0.13, 0.33, 0.3],
            'width': 2.2,
            'margin_minimum_right': 0.03,
            'margin_left': 0.2
        }
        self.nta = nta
    def get_num_columns(self):
        table_settings = self.nta.table_settings
        num_column_names = len(table_settings['column_names']) + len(table_settings['results_column_names'])
        assert len(table_settings['column_widths']) == num_column_names, "Unequal numbers of column widths and names."
        return num_column_names
    def test_number_of_columns(self):
        self.nta.enable_table(**self.table_settings)
        self.get_num_columns()

if __name__ == '__main__':
    unittest.main()