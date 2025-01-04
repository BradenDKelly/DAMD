#import pytest
#from ..src import parameter_reader     #.parameter_reader import read_top_file
import numpy as np
import sys
import os
import unittest

current_dir = os.path.dirname(__file__)
src_dir = os.path.abspath(os.path.join(current_dir, '..', 'src'))
sys.path.append(src_dir)

from parameter_reader import read_top_file

class TestParameterReader(unittest.TestCase):
    def test_read_topology(self):
        print(f"Starting Now")
        # Sample topology input for testing
        top_file = "spce.top"
        expected_masses = {"1": 1.008, "2": 16.00}  # Replace with actual test data
        expected_bonds = [("1", "2")]
          
        parameters = read_top_file(top_file)

        print(parameters)
        assert np.isclose(parameters[2][0].mass, 15.99940)
        assert np.isclose(parameters[2][1].mass, 1.008)
        #assert parameters["masses"] == expected_masses
        #assert parameters["bonds"] == expected_bonds

        

if __name__ == '__main__':
    unittest.main()
