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
        self.table_settings = {
            'width': 2.2,
            'margin_minimum_right': 0.03,
            'margin_left': 0.2
        }
        self.nta = nta
    def add_columns(self):
        nta = self.nta
        nta.table_add_treatments_and_waits("Treatment\n{treatment_number}\n(µM)", 0.2, "4°C\nwait\n{wait_number}\n(h)", 0.07)
        names = ["Filter\ncut-on\n(nm)", "Power\n(mW)", "Exposure,\ngain", "Detection\nsetting", "Video sec\nx quantity", "Stir sec\nx RPM", "ID", "ID of\nprevious"]
        widths = [0.1, 0.19, 0.14, 0.19, 0.16, 0.12, 0.1, 0.13]
        assert len(names) == len(widths), f"{len(names)=} does not equal {len(widths)=}"
        for name, width in zip(names, widths):
            nta.table_add_column(name, width)
        nta.results_enable_time("Time (s)")
        nta.table_add_time(width = 0.33)
        nta.results_enable_concentration("Concentration\n(counts/mL)")
        nta.table_add_concentration(width = 0.3)
    def get_num_columns(self):
        table_settings = self.nta.table_settings
        num_column_names = len(table_settings['column_names'])
        assert len(table_settings['column_widths']) == num_column_names, "Unequal numbers of column widths and names."
        return num_column_names
    def test_number_of_columns(self):
        nta = self.nta
        nta.enable_table(**self.table_settings)
        self.add_columns()
        self.get_num_columns()

if __name__ == '__main__':
    unittest.main()