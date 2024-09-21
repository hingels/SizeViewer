import unittest
from src.nanotracking import DifferencePlotter

class Test_Basics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.nta = DifferencePlotter.NTA(
            datafolder = "tests/Test data",
            output_folder = "tests/Test output",
            filenames = ["1"]
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