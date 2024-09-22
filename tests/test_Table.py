import unittest
from src.nanotracking import DifferencePlotter
from src.nanotracking import settings_classes

class Test_Table(unittest.TestCase):
    filenames = ["1", "1.2", "1.3"]
    specifiers = "All measurements"
    def setUp(self):
        filenames = self.filenames
        nta = DifferencePlotter.NTA(
            datafolder = "tests/Test data",
            output_folder = f"tests/Test output/{self.specifiers}/{self.id()}",
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
        # names = ["Filter\ncut-on\n(nm)", "Power\n(mW)", "Exposure,\ngain", "Detection\nsetting", "Video sec\nx quantity", "Stir sec\nx RPM", "ID", "ID of\nprevious"]
        # widths = [0.1, 0.19, 0.14, 0.19, 0.16, 0.12, 0.1, 0.13]
        # assert len(names) == len(widths), f"{len(names)=} does not equal {len(widths)=}"
        nta.table_add_settings_by_tag('filter', column_name = "Filter\ncut-on\n(nm)", column_width = 0.1)
        nta.table_add_settings_by_tag('RedLaserPower', 'GreenLaserPower', 'BlueLaserPower', column_name = "Power\n(mW)", column_width = 0.19)
        nta.table_add_settings_by_tag('Exposure', column_name = "Exposure,\ngain", column_width = 0.14)
        def get_detection_info(threshold_type, threshold):
            if threshold_type == 'Polydisperse': return threshold_type
            return f"{threshold_type}\n{threshold}"
        nta.table_add_settings_by_tag('DetectionThresholdType', 'DetectionThreshold', column_name = "Detection\nsetting", column_width = 0.19, format_callback = get_detection_info)
        
        def get_video_info(framerate, frames_per_video, num_of_videos):
            video_duration = frames_per_video / framerate
            if video_duration.is_integer():
                video_duration = int(video_duration)
            return f"{video_duration}x{num_of_videos}"
        nta.table_add_settings_by_tag('FrameRate', 'FramesPerVideo', 'NumOfVideos', column_name = "Video sec\nx quantity", column_width = 0.16, format_callback = get_video_info)

        def get_stir_info(stir_time, stir_rpm):
            return f"{stir_time}x{stir_rpm}"
        nta.table_add_settings_by_tag('StirredTime', 'StirrerSpeed', column_name = "Stir sec\nx RPM", column_width = 0.12, format_callback = get_stir_info)
        
        def get_ID_info(ID):
            return '\n'.join((ID[0:4], ID[4:8], ID[8:12]))
        nta.table_add_settings_by_tag('ID', column_name = "ID", column_width = 0.1, format_callback = get_ID_info)
        
        def get_previous_ID_info(previous):
            previous_sample = nta.unordered_samples[previous]
            ID_of_previous = nta.settings.by_tag('ID').get_value(previous_sample)
            return '\n'.join((ID_of_previous[0:4], ID_of_previous[4:8], ID_of_previous[8:12]))
        # previous_ID_setting = settings_classes.Setting('previous_ID_setting', column_name = "ID of\nprevious", column_width = 0.13, format_callback = get_previous_ID_info)
        # nta.settings.add_setting(previous_ID_setting.tag, previous_ID_setting)
        nta.table_add_settings_by_tag('previous', column_name = "ID of\nprevious", column_width = 0.13, format_callback = get_previous_ID_info)

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

class Test_Table_OneMeasurement(Test_Table):
    filenames = ["1"]
    specifiers = "One measurement"
class Test_Table_TwoMeasurements(Test_Table):
    filenames = ["1", "1.2"]
    specifiers = "Two measurements"

if __name__ == '__main__':
    unittest.main()