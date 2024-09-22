import unittest
from src.nanotracking import DifferencePlotter

class Test_Table(unittest.TestCase):
    def setUp(self):
        filenames = ["1", "1.2", "1.3"]
        nta = DifferencePlotter.NTA(
            datafolder = "tests/Test data",
            output_folder = f"tests/Test output/{self.id()}",
            filenames = filenames
        )
        nta.compute()
        self.assertEqual(nta.num_of_plots, len(filenames), "Number of filenames is not equal to number of plots.")
        self.num_of_plots = nta.num_of_plots
        self.table_options = {
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
        nta.table_add_time("Time (s)", 0.33)
        nta.table_add_concentration("Concentration\n(counts/mL)", 0.3)
    def get_num_columns(self):
        table_settings = self.nta.table_settings
        num_column_names = len(table_settings['column_names_without_treatmentsOrWaits'])
        assert len(table_settings['column_widths_without_treatmentsOrWaits']) == num_column_names, "Unequal numbers of column widths and names."
        return num_column_names

    def test_number_of_columns(self):
        nta = self.nta
        nta.enable_table(**self.table_options)
        self.add_columns()
        self.get_num_columns()
    
    def setup_test_persistence(self):
        nta = self.nta
        nta.enable_table(**self.table_options)
        self.add_columns()
        num_columns = self.get_num_columns()
        nta.compute(prep_tabulation = False)
        self.assertEqual(num_columns, self.get_num_columns(), "Column count changed after running NTA.compute().")
        nta.prepare_tabulation()
        self.assertEqual(num_columns, self.get_num_columns(), "Column count changed after running NTA.prepare_tabulation().")
        self.assertEqual(self.num_of_plots, nta.num_of_plots, "Number of plots has changed.")
        nta.plot(name = "Initial plot")
        self.assertEqual(num_columns, self.get_num_columns(), "Column count changed after running NTA.plot().")
        return num_columns
    def finish_test_persistence(self, num_columns):
        nta = self.nta
        nta.compute(prep_tabulation = False)
        self.assertEqual(num_columns, self.get_num_columns(), "Column count changed after running NTA.compute().")
        nta.prepare_tabulation()
        self.assertEqual(num_columns, self.get_num_columns(), "Column count changed after running NTA.prepare_tabulation().")
        nta.plot(name = "Final plot")
        self.assertEqual(num_columns, self.get_num_columns(), "Column count changed after running NTA.plot().")
    def test_persistence_with_peakfinding(self):
        nta = self.nta
        num_columns = self.setup_test_persistence()
        nta.enable_peak_detection(
            kernel_size = 30,
            kernel2_size = 20,
            kernel_std_in_bins = 4,
            second_derivative_threshold = -30,
            maxima_marker = {'marker': 'o', 'fillstyle': 'none', 'color': 'black', 'linestyle': 'none'},
            rejected_maxima_marker = {'marker': 'o', 'fillstyle': 'none', 'color': '0.5', 'linestyle': 'none'}
        )
        self.assertEqual(num_columns, self.get_num_columns(), "Column count changed after running NTA.enable_peak_detection().")
        self.finish_test_persistence(num_columns)
    def test_persistence_with_cumulative(self):
        nta = self.nta
        num_columns = self.setup_test_persistence()
        nta.enable_cumulative()
        self.assertEqual(num_columns, self.get_num_columns(), "Column count changed after running NTA.enable_cumulative().")
        self.finish_test_persistence(num_columns)
    def test_persistence_with_difference(self):
        nta = self.nta
        num_columns = self.setup_test_persistence()
        nta.enable_difference()
        self.assertEqual(num_columns, self.get_num_columns(), "Column count changed after running NTA.enable_difference().")
        self.finish_test_persistence(num_columns)

if __name__ == '__main__':
    unittest.main()