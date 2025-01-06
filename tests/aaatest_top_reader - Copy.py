#import pytest
#from ..src import parameter_reader     #.parameter_reader import read_top_file
import numpy as np
import sys
import os
import unittest

current_dir = os.path.dirname(__file__)
src_dir = os.path.abspath(os.path.join(current_dir, '..', 'src'))
sys.path.append(src_dir)

from topp_reader import read_top_file

class TestParameterReader(unittest.TestCase):
    def test_read_topology(self):
        print(f"Starting Now")
        # Sample topology input for testing
        #top_file = "spce.top"
        #parameters = read_top_file(top_file)
        #atomtypes_list = parameters.atomtypes

        #assert np.isclose(atomtypes_list[0].mass, 15.99940)
        #assert np.isclose(atomtypes_list[1].mass, 1.008)
        
        # Test MEA and tip3p .top file
        top_file = "mea_tip3p.top"
        parameters = read_top_file(top_file)
        atomtypes_list = parameters.atomtypes

        for i, mols in enumerate(parameters.molparams):
            for j, atoms in enumerate(mols.atoms):
               
                print(f"{i}{j}  {atoms.charge}")
        # water
        assert np.isclose(atomtypes_list[0].mass, 15.99940)
        assert np.isclose(atomtypes_list[1].mass, 1.008)
        # mea
        assert np.isclose(atomtypes_list[2].mass, 16.00)    # check first atoms mass
        assert np.isclose(atomtypes_list[12].mass, 1.00)    # check last atoms mass
        assert np.isclose(atomtypes_list[11].charge, 0.345190) # check second last atoms charge
        
        # This tests that #include "mea_added.itp" worked as it should to get parameters of second MEA molecule in .itp file
        assert np.isclose(parameters.molparams[1].atoms[0].mass, 16.00)    # check first atoms mass of second MEA 
        assert np.isclose(parameters.molparams[1].atoms[10].charge, 0.352192) # check last atoms charge in second MEA
        
        assert len(parameters.molparams) == 3
            
        

if __name__ == '__main__':
    unittest.main()
