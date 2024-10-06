import unittest
from src.nanotracking import DifferencePlotter
from src.nanotracking import settings_classes
import numpy as np

class Test_Table(unittest.TestCase):
    filenames = ["1-1e5 150nm Nanosphere", "1-1e5 150nm Nanosphere 2", "1-1e5 150nm Nanosphere 32ms", "1-1e5 150nm Nanosphere diff detection setting"]
    specifiers = "All measurements"
    def setUp(self):
        filenames = self.filenames
        nta = DifferencePlotter.NTA(
            datafolder = "tests/Test data",
            output_folder = f"tests/Test output/{self.specifiers}/{self.id()}",
            filenames = filenames,
            truncation_size = 400 # nanometers
        )
        nta.compute()
        self.assertEqual(nta.num_of_plots, len(filenames), "Number of filenames is not equal to number of plots.")
        self.num_of_plots = nta.num_of_plots
        self.table_options = {
            'width': 2.5,
            'margin_minimum_right': 0.03,
            'margin_left': 0.2
        }
        self.nta = nta
        self.table = None
    def add_columns(self):
        nta = self.nta
        table = self.table
        table.add_treatments_and_waits("Treatment\n{treatment_number}\n(µM)", 0.2, "4°C\nwait\n{wait_number}\n(h)", 0.07)
        table.add_settings_by_tag('filter', column_name = "Filter\ncut-on\n(nm)", column_width = 0.1)
        table.add_settings_by_tag('RedLaserPower', 'GreenLaserPower', 'BlueLaserPower', column_name = "Power\n(mW)", column_width = 0.19)
        table.add_settings_by_tag('Exposure', column_name = "Exposure,\ngain", column_width = 0.14)
        def get_detection_info(threshold_type, threshold):
            if threshold is None: return threshold_type
            if threshold_type == 'Polydisperse': return threshold_type
            return f"{threshold_type}\n{threshold}"
        table.add_settings_by_tag('DetectionThresholdType', 'DetectionThreshold', column_name = "Detection\nsetting", column_width = 0.19, format = get_detection_info)
        
        def get_video_info(framerate, frames_per_video, num_of_videos):
            video_duration = frames_per_video / framerate
            if video_duration.is_integer():
                video_duration = int(video_duration)
            return f"{video_duration}x{num_of_videos}"
        table.add_settings_by_tag('FrameRate', 'FramesPerVideo', 'NumOfVideos', column_name = "Video sec\nx quantity", column_width = 0.16, format = get_video_info)

        # def get_stir_info(stir_time, stir_rpm):
        #     return f"{stir_time}x{stir_rpm}"
        # table.add_settings_by_tag('StirredTime', 'StirrerSpeed', column_name = "Stir sec\nx RPM", column_width = 0.12, format = get_stir_info)
        stir_format = '{stir_time}x{stir_rpm}'
        table.add_settings_by_tag('StirredTime', 'StirrerSpeed', column_name = "Stir sec\nx RPM", column_width = 0.12, format = stir_format)
        
        def get_ID_info(ID):
            return '\n'.join((ID[0:4], ID[4:8], ID[8:12]))
        table.add_settings_by_tag('ID', column_name = "ID", column_width = 0.1, format = get_ID_info)
        
        settings = nta.settings
        samples, unordered_samples = nta.samples, nta.unordered_samples
        def get_previous_ID_info(previous):
            if previous is None: return ''
            previous_sample = unordered_samples[previous]
            ID_of_previous = settings.by_tag('ID').get_value(previous_sample)
            return '\n'.join((ID_of_previous[0:4], ID_of_previous[4:8], ID_of_previous[8:12]))
        table.add_settings_by_tag('previous', column_name = "ID of\nprevious", column_width = 0.13, format = get_previous_ID_info)

        times = settings.by_tag('time')
        def value_function(sample):
            truncation_size = 200
            counts = nta.counts(sample = sample)
            truncated_counts = nta.counts(sample = sample, truncation_size = truncation_size)

            size_binwidth = nta.size_binwidth
            truncated_sizes = nta.sizes(sample = sample, truncation_size = truncation_size, lower = False, upper = True)
            top_nm = max(truncated_sizes)
            if top_nm.is_integer():
                top_nm = int(top_nm)
            total_conc = np.sum(counts)*size_binwidth
            total_conc_under_topnm = np.sum(truncated_counts*size_binwidth)
            
            sample_index = sample.index
            time = times.get_value(sample) # times accessed via closure
            time_since_previous = None
            previous = settings.by_tag('previous').get_value(sample)
            if previous is not None:
                if previous not in unordered_samples:
                    time_since_previous = '?'
                else:
                    previous_sample = unordered_samples[previous]
                    time_of_previous = settings.by_tag('time').get_value(previous_sample)
                    time_since_previous = int((time - time_of_previous).total_seconds())
            try: above = samples[sample_index - 1]
            except: above = None
            time_since_above = None
            if above is not None:
                time_of_above = times.get_value(above)
                time_since_above = int((time - time_of_above).total_seconds())
            return previous, time_since_previous, time_since_above, f"{total_conc:.2E}", f"{total_conc_under_topnm:.2E}", top_nm
        calculation = nta.new_calculation(
            'Previous/time/concentrations', value_function,
            'previous', 'time_since_previous', 'time_since_above', 'total_conc', 'total_conc_under_topnm', 'top_nm')
        
        def get_time_info(previous, time_since_previous, time_since_above, total_conc, total_conc_under_topnm, top_nm):
            text = []
            if time_since_above is not None:
                text.append(f"{time_since_above} since above")
            if time_since_previous is not None:
                text.append(f"{time_since_previous} since previous")
            return '\n'.join(text)
        calculation.add_format('Time format', get_time_info)

        # def get_conc_info(previous, time_since_previous, time_since_above, total_conc, total_conc_under_topnm, top_nm):
        #     return f"Total: {total_conc}\n<{top_nm}nm: {total_conc_under_topnm}"
        # calculation.add_format('Concentration format', get_conc_info)
        calculation.add_format('Concentration format', 'Total: {total_conc}\n<{top_nm}nm: {total_conc_under_topnm}')
        
        table.add_calculation(calculation, 'Time format', column_name = "Time (s)", column_width = 0.33)
        table.add_calculation(calculation, 'Concentration format', column_name = "Concentration\n(counts/mL)", column_width = 0.3)

        def get_sample_name(sample):
            return sample.name
        table.add_settings_by_tag('sample', column_name = "Sample name", column_width = 0.25, format = get_sample_name)
            

    def get_num_columns(self):
        table = self.table
        num_column_names = len(table.column_names_without_treatmentsOrWaits)
        assert len(table.column_widths_without_treatmentsOrWaits) == num_column_names, "Unequal numbers of column widths and names."
        return num_column_names

    def test_number_of_columns(self):
        nta = self.nta
        self.table = nta.add_table(**self.table_options)
        nta.enable_table()
        self.add_columns()
        self.get_num_columns()
    
    def setup_test_persistence(self):
        nta = self.nta
        self.table = nta.add_table(**self.table_options)
        nta.enable_table()
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
        nta.compare()
        self.assertEqual(num_columns, self.get_num_columns(), "Column count changed after running NTA.compare().")
    def test_persistence_with_peakfinding(self):
        nta = self.nta
        num_columns = self.setup_test_persistence()
        nta.enable_peak_detection(
            gaussian_width = 30,
            gaussian_std_in_bins = 4,
            moving_avg_width = 20,
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
    filenames = ["1-1e5 150nm Nanosphere"]
    specifiers = "One measurement"
class Test_Table_TwoMeasurements(Test_Table):
    filenames = ["1-1e5 150nm Nanosphere", "1-1e5 150nm Nanosphere 2"]
    specifiers = "Two measurements"

if __name__ == '__main__':
    unittest.main()