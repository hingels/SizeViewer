import unittest
from src.nanotracking import DifferencePlotter

class Test_Basics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.nta = DifferencePlotter.NTA(
            datafolder = "tests/Test data",
            output_folder = "tests/Test output/Basics",
            filenames = ["1-1e5 150nm Nanosphere", "1-1e5 150nm Nanosphere 2", "1-1e5 150nm Nanosphere 32ms", "1-1e5 150nm Nanosphere diff detection setting"],
            truncation_size = 400 # nanometers
        )
    def test_compute(self):
        self.nta.compute()
    def test_compare(self):
        self.nta.compute()
        self.nta.compare()
    def test_plot(self):
        self.nta.plot()

if __name__ == '__main__':
    unittest.main()