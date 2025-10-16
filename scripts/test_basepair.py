#!/usr/bin/env python3
"""
Unit tests for the PDB base pair analysis tool.

This module contains comprehensive tests for the PDBBasePairAnalyzer class
and related functionality.
"""

import unittest
import numpy as np
import tempfile
import os
from unittest.mock import Mock, patch, MagicMock
from pdb_basepair import PDBBasePairAnalyzer, analyze_base_pairs


class TestPDBBasePairAnalyzer(unittest.TestCase):
    """Test cases for PDBBasePairAnalyzer class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = PDBBasePairAnalyzer(distance_cutoff=10.0)
        
        # Mock residue data for testing
        self.mock_residues = [
            Mock(resname='G', id=(None, 1)),
            Mock(resname='C', id=(None, 2)),
            Mock(resname='A', id=(None, 3)),
            Mock(resname='U', id=(None, 4)),
        ]
        
        # Mock C4' atoms with coordinates
        self.mock_coords = np.array([
            [0.0, 0.0, 0.0],   # G
            [5.0, 0.0, 0.0],   # C (5 Å from G)
            [0.0, 10.0, 0.0],  # A (10 Å from G)
            [5.0, 10.0, 0.0],  # U (5 Å from A, ~11.2 Å from G)
        ])
        
        # Mock atoms
        for i, residue in enumerate(self.mock_residues):
            atom = Mock()
            atom.coord = self.mock_coords[i]
            residue.__getitem__ = Mock(return_value=atom)
    
    def test_init(self):
        """Test analyzer initialization."""
        analyzer = PDBBasePairAnalyzer(distance_cutoff=12.0, chains=['A', 'B'])
        self.assertEqual(analyzer.distance_cutoff, 12.0)
        self.assertEqual(analyzer.chains, ['A', 'B'])
        self.assertIsNone(analyzer.structure)
        self.assertEqual(len(analyzer.residues), 0)
    
    def test_base_pair_rules(self):
        """Test canonical base pairing rules."""
        # Valid pairs
        self.assertTrue(('G', 'C') in self.analyzer.base_pairs)
        self.assertTrue(('C', 'G') in self.analyzer.base_pairs)
        self.assertTrue(('A', 'U') in self.analyzer.base_pairs)
        self.assertTrue(('U', 'A') in self.analyzer.base_pairs)
        self.assertTrue(('G', 'U') in self.analyzer.base_pairs)
        self.assertTrue(('U', 'G') in self.analyzer.base_pairs)
        self.assertTrue(('A', 'T') in self.analyzer.base_pairs)
        self.assertTrue(('T', 'A') in self.analyzer.base_pairs)
        
        # Invalid pairs
        self.assertFalse(('A', 'C') in self.analyzer.base_pairs)
        self.assertFalse(('G', 'A') in self.analyzer.base_pairs)
        self.assertFalse(('C', 'U') in self.analyzer.base_pairs)
    
    def test_is_nucleic_acid(self):
        """Test nucleic acid residue identification."""
        # Standard bases
        self.assertTrue(self.analyzer._is_nucleic_acid(Mock(resname='A')))
        self.assertTrue(self.analyzer._is_nucleic_acid(Mock(resname='U')))
        self.assertTrue(self.analyzer._is_nucleic_acid(Mock(resname='G')))
        self.assertTrue(self.analyzer._is_nucleic_acid(Mock(resname='C')))
        self.assertTrue(self.analyzer._is_nucleic_acid(Mock(resname='T')))
        
        # DNA bases
        self.assertTrue(self.analyzer._is_nucleic_acid(Mock(resname='DA')))
        self.assertTrue(self.analyzer._is_nucleic_acid(Mock(resname='DT')))
        
        # RNA bases
        self.assertTrue(self.analyzer._is_nucleic_acid(Mock(resname='RA')))
        self.assertTrue(self.analyzer._is_nucleic_acid(Mock(resname='RU')))
        
        # Monophosphates
        self.assertTrue(self.analyzer._is_nucleic_acid(Mock(resname='AMP')))
        self.assertTrue(self.analyzer._is_nucleic_acid(Mock(resname='UMP')))
        
        # Non-nucleic acid
        self.assertFalse(self.analyzer._is_nucleic_acid(Mock(resname='ALA')))
        self.assertFalse(self.analyzer._is_nucleic_acid(Mock(resname='GLY')))
    
    def test_get_c4_atom(self):
        """Test C4' atom extraction."""
        # Mock residue with C4' atom
        residue = Mock()
        c4_atom = Mock()
        residue.__getitem__ = Mock(return_value=c4_atom)
        
        result = self.analyzer._get_c4_atom(residue)
        self.assertEqual(result, c4_atom)
        residue.__getitem__.assert_called_with("C4'")
    
    def test_get_c4_atom_alternative_names(self):
        """Test C4' atom extraction with alternative names."""
        # Test C4* name
        residue = Mock()
        c4_atom = Mock()
        residue.__getitem__ = Mock(side_effect=lambda x: c4_atom if x == "C4*" else None)
        
        result = self.analyzer._get_c4_atom(residue)
        self.assertEqual(result, c4_atom)
    
    def test_get_c4_atom_not_found(self):
        """Test C4' atom extraction when not found."""
        residue = Mock()
        residue.__getitem__ = Mock(return_value=None)
        
        result = self.analyzer._get_c4_atom(residue)
        self.assertIsNone(result)
    
    def test_get_base_type(self):
        """Test base type extraction from residue names."""
        # Standard bases
        self.assertEqual(self.analyzer._get_base_type('A'), 'A')
        self.assertEqual(self.analyzer._get_base_type('U'), 'U')
        self.assertEqual(self.analyzer._get_base_type('G'), 'G')
        self.assertEqual(self.analyzer._get_base_type('C'), 'C')
        self.assertEqual(self.analyzer._get_base_type('T'), 'T')
        
        # DNA bases
        self.assertEqual(self.analyzer._get_base_type('DA'), 'A')
        self.assertEqual(self.analyzer._get_base_type('DT'), 'T')
        
        # RNA bases
        self.assertEqual(self.analyzer._get_base_type('RA'), 'A')
        self.assertEqual(self.analyzer._get_base_type('RU'), 'U')
        
        # Monophosphates
        self.assertEqual(self.analyzer._get_base_type('AMP'), 'A')
        self.assertEqual(self.analyzer._get_base_type('UMP'), 'U')
        
        # Diphosphates
        self.assertEqual(self.analyzer._get_base_type('ADP'), 'A')
        self.assertEqual(self.analyzer._get_base_type('UDP'), 'U')
        
        # Triphosphates
        self.assertEqual(self.analyzer._get_base_type('ATP'), 'A')
        self.assertEqual(self.analyzer._get_base_type('UTP'), 'U')
    
    def test_calculate_distances(self):
        """Test distance calculation between C4' atoms."""
        self.analyzer.coordinates = self.mock_coords
        
        distances = self.analyzer._calculate_distances()
        
        # Check matrix shape
        self.assertEqual(distances.shape, (4, 4))
        
        # Check symmetry
        self.assertTrue(np.allclose(distances, distances.T))
        
        # Check specific distances
        self.assertAlmostEqual(distances[0, 1], 5.0)  # G-C distance
        self.assertAlmostEqual(distances[0, 2], 10.0)  # G-A distance
        self.assertAlmostEqual(distances[1, 3], 10.0)  # C-U distance
    
    def test_get_base_pair_matrix(self):
        """Test base pair matrix generation."""
        # Set up mock data
        self.analyzer.residues = self.mock_residues
        self.analyzer.coordinates = self.mock_coords
        self.analyzer.residue_info = [
            {'residue_name': 'G', 'chain': 'A', 'residue_number': 1},
            {'residue_name': 'C', 'chain': 'A', 'residue_number': 2},
            {'residue_name': 'A', 'chain': 'A', 'residue_number': 3},
            {'residue_name': 'U', 'chain': 'A', 'residue_number': 4},
        ]
        
        matrix = self.analyzer.get_base_pair_matrix()
        
        # Check matrix properties
        self.assertEqual(matrix.shape, (4, 4))
        self.assertTrue(np.allclose(matrix, matrix.T))  # Symmetric
        
        # G-C pair should be found (distance = 5.0 Å)
        self.assertTrue(matrix[0, 1])
        self.assertTrue(matrix[1, 0])
        
        # A-U pair should be found (distance = 5.0 Å)
        self.assertTrue(matrix[2, 3])
        self.assertTrue(matrix[3, 2])
        
        # G-A should not be a pair (distance = 10.0 Å, but A-G is not canonical)
        self.assertFalse(matrix[0, 2])
        self.assertFalse(matrix[2, 0])
    
    def test_get_base_pair_matrix_no_residues(self):
        """Test base pair matrix generation with no residues."""
        with self.assertRaises(ValueError):
            self.analyzer.get_base_pair_matrix()
    
    def test_get_base_pairs(self):
        """Test detailed base pair information extraction."""
        # Set up mock data
        self.analyzer.residues = self.mock_residues
        self.analyzer.coordinates = self.mock_coords
        self.analyzer.residue_info = [
            {'residue_name': 'G', 'chain': 'A', 'residue_number': 1},
            {'residue_name': 'C', 'chain': 'A', 'residue_number': 2},
            {'residue_name': 'A', 'chain': 'A', 'residue_number': 3},
            {'residue_name': 'U', 'chain': 'A', 'residue_number': 4},
        ]
        
        # Generate matrix first
        self.analyzer.get_base_pair_matrix()
        
        pairs = self.analyzer.get_base_pairs()
        
        # Should find G-C and A-U pairs
        self.assertEqual(len(pairs), 2)
        
        # Check pair information
        pair_info = {(p['base1'], p['base2']) for p in pairs}
        self.assertIn(('G', 'C'), pair_info)
        self.assertIn(('A', 'U'), pair_info)
    
    @patch('matplotlib.pyplot.show')
    @patch('matplotlib.pyplot.savefig')
    def test_plot_matrix(self, mock_savefig, mock_show):
        """Test matrix visualization."""
        # Set up mock data
        self.analyzer.residues = self.mock_residues
        self.analyzer.coordinates = self.mock_coords
        self.analyzer.residue_info = [
            {'residue_name': 'G', 'chain': 'A', 'residue_number': 1},
            {'residue_name': 'C', 'chain': 'A', 'residue_number': 2},
        ]
        
        # Generate matrix first
        self.analyzer.get_base_pair_matrix()
        
        # Test plotting
        self.analyzer.plot_matrix(save_path='test_matrix.png')
        
        # Check that savefig was called
        mock_savefig.assert_called_once()
        mock_show.assert_called_once()
    
    def test_print_summary(self):
        """Test summary printing."""
        # Set up mock data
        self.analyzer.residues = self.mock_residues
        self.analyzer.coordinates = self.mock_coords
        self.analyzer.residue_info = [
            {'residue_name': 'G', 'chain': 'A', 'residue_number': 1},
            {'residue_name': 'C', 'chain': 'A', 'residue_number': 2},
        ]
        
        # Generate matrix first
        self.analyzer.get_base_pair_matrix()
        
        # This should not raise an exception
        self.analyzer.print_summary()


class TestConvenienceFunction(unittest.TestCase):
    """Test cases for the convenience function."""
    
    @patch('pdb_basepair.PDBBasePairAnalyzer')
    def test_analyze_base_pairs(self, mock_analyzer_class):
        """Test the convenience function."""
        # Mock the analyzer instance
        mock_analyzer = Mock()
        mock_analyzer_class.return_value = mock_analyzer
        mock_analyzer.residue_info = [{'residue_name': 'G'}]
        
        # Mock the matrix
        mock_matrix = np.array([[False, True], [True, False]])
        mock_analyzer.get_base_pair_matrix.return_value = mock_matrix
        
        # Test the function
        matrix, residue_info = analyze_base_pairs('test.pdb', distance_cutoff=12.0, chains=['A'])
        
        # Check that analyzer was created with correct parameters
        mock_analyzer_class.assert_called_once_with(distance_cutoff=12.0, chains=['A'])
        
        # Check that methods were called
        mock_analyzer.load_pdb.assert_called_once_with('test.pdb')
        mock_analyzer.get_base_pair_matrix.assert_called_once()
        
        # Check return values
        np.testing.assert_array_equal(matrix, mock_matrix)
        self.assertEqual(residue_info, [{'residue_name': 'G'}])


class TestIntegration(unittest.TestCase):
    """Integration tests with mock PDB files."""
    
    def create_mock_pdb_file(self, content):
        """Create a temporary PDB file with given content."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(content)
            return f.name
    
    def tearDown(self):
        """Clean up temporary files."""
        # Clean up any temporary files created during tests
        pass
    
    @patch('pdb_basepair.PDBParser')
    def test_load_pdb_integration(self, mock_parser_class):
        """Test PDB loading with mock parser."""
        # Create mock structure
        mock_structure = Mock()
        mock_model = Mock()
        mock_chain = Mock()
        mock_chain.id = 'A'
        
        # Mock residue with C4' atom
        mock_residue = Mock()
        mock_residue.resname = 'G'
        mock_residue.id = (None, 1)
        mock_c4_atom = Mock()
        mock_c4_atom.coord = np.array([0.0, 0.0, 0.0])
        mock_residue.__getitem__ = Mock(return_value=mock_c4_atom)
        
        mock_chain.__iter__ = Mock(return_value=iter([mock_residue]))
        mock_model.__iter__ = Mock(return_value=iter([mock_chain]))
        mock_structure.__iter__ = Mock(return_value=iter([mock_model]))
        
        mock_parser = Mock()
        mock_parser.get_structure.return_value = mock_structure
        mock_parser_class.return_value = mock_parser
        
        # Test loading
        analyzer = PDBBasePairAnalyzer()
        pdb_file = self.create_mock_pdb_file("ATOM      1  P   G A   1      0.000   0.000   0.000")
        
        try:
            analyzer.load_pdb(pdb_file)
            self.assertEqual(len(analyzer.residues), 1)
            self.assertEqual(len(analyzer.coordinates), 1)
            np.testing.assert_array_equal(analyzer.coordinates[0], [0.0, 0.0, 0.0])
        finally:
            os.unlink(pdb_file)


def run_tests():
    """Run all tests."""
    # Create test suite
    test_suite = unittest.TestSuite()
    
    # Add test cases
    test_suite.addTest(unittest.makeSuite(TestPDBBasePairAnalyzer))
    test_suite.addTest(unittest.makeSuite(TestConvenienceFunction))
    test_suite.addTest(unittest.makeSuite(TestIntegration))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(test_suite)
    
    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_tests()
    exit(0 if success else 1)
