import unittest
from src.nanotracking import DifferencePlotter

class Test_Table(unittest.TestCase):
    def setUp(self):
        nta = DifferencePlotter.NTA(
            datafolder = "tests/Test data",
            output_folder = f"tests/Test output/{self.id()}",
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
    def check_table_persistence(self):
        nta = self.nta
        nta.enable_table(**self.table_settings)
        self.add_columns()
        num_columns = self.get_num_columns()
        nta.compute()
        self.assertEqual(num_columns, self.get_num_columns(), "Column count changed after running NTA.compute().")
        nta.plot(name = "Initial plot")
        self.assertEqual(num_columns, self.get_num_columns(), "Column count changed after running NTA.plot().")
        return num_columns
    def test_persistence_with_peakfinding(self):
        nta = self.nta
        num_columns = self.check_table_persistence()
        nta.enable_peak_detection(
            kernel_size = 30,
            kernel2_size = 20,
            kernel_std_in_bins = 4,
            second_derivative_threshold = -30,
            maxima_marker = {'marker': 'o', 'fillstyle': 'none', 'color': 'black', 'linestyle': 'none'},
            rejected_maxima_marker = {'marker': 'o', 'fillstyle': 'none', 'color': '0.5', 'linestyle': 'none'}
        )
        self.assertEqual(num_columns, self.get_num_columns(), "Column count changed after running NTA.enable_peak_detection().")
        nta.compute()
        self.assertEqual(num_columns, self.get_num_columns(), "Column count changed after running NTA.compute().")
        nta.plot(name = "Final plot")
        self.assertEqual(num_columns, self.get_num_columns(), "Column count changed after running NTA.plot().")

if __name__ == '__main__':
    unittest.main()