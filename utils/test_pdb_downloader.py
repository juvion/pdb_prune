#!/usr/bin/env python3

import unittest
import os
import shutil
from pdb_downloader import download_rna_pdbs
import logging

class TestPDBDownloader(unittest.TestCase):
    def setUp(self):
        """Set up test environment before each test."""
        self.test_dir = "test_downloads"
        # Disable logging during tests
        logging.disable(logging.CRITICAL)
    
    def tearDown(self):
        """Clean up after each test."""
        # Remove test directory if it exists
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
        # Re-enable logging
        logging.disable(logging.NOTSET)
    
    def test_basic_download(self):
        """Test basic download functionality with default parameters."""
        download_rna_pdbs(
            download_directory=self.test_dir,
            max_entries=5,
            start_index=0
        )
        # Check if directory was created
        self.assertTrue(os.path.exists(self.test_dir))
        # Check if files were downloaded
        pdb_files = [f for f in os.listdir(self.test_dir) if f.endswith('.pdb')]
        self.assertLessEqual(len(pdb_files), 5)
    
    def test_range_download(self):
        """Test downloading from a specific range."""
        download_rna_pdbs(
            download_directory=self.test_dir,
            max_entries=5,
            start_index=100
        )
        pdb_files = [f for f in os.listdir(self.test_dir) if f.endswith('.pdb')]
        self.assertLessEqual(len(pdb_files), 5)
    
    def test_sort_by_resolution(self):
        """Test sorting by resolution."""
        download_rna_pdbs(
            download_directory=self.test_dir,
            max_entries=5,
            sort_by="resolution",
            sort_order="asc"
        )
        pdb_files = [f for f in os.listdir(self.test_dir) if f.endswith('.pdb')]
        self.assertLessEqual(len(pdb_files), 5)
    
    def test_sort_by_sequence_length(self):
        """Test sorting by sequence length."""
        download_rna_pdbs(
            download_directory=self.test_dir,
            max_entries=5,
            sort_by="sequence_length",
            sort_order="desc"
        )
        pdb_files = [f for f in os.listdir(self.test_dir) if f.endswith('.pdb')]
        self.assertLessEqual(len(pdb_files), 5)
    
    def test_invalid_sort_field(self):
        """Test handling of invalid sort field."""
        with self.assertRaises(Exception):
            download_rna_pdbs(
                download_directory=self.test_dir,
                max_entries=5,
                sort_by="invalid_field"
            )
    
    def test_negative_start_index(self):
        """Test handling of negative start index."""
        with self.assertRaises(Exception):
            download_rna_pdbs(
                download_directory=self.test_dir,
                max_entries=5,
                start_index=-1
            )
    
    def test_zero_max_entries(self):
        """Test handling of zero max entries."""
        download_rna_pdbs(
            download_directory=self.test_dir,
            max_entries=0
        )
        pdb_files = [f for f in os.listdir(self.test_dir) if f.endswith('.pdb')]
        self.assertEqual(len(pdb_files), 0)

def run_tests():
    """Run all tests and print results."""
    print("\nRunning PDB Downloader Tests...")
    print("-" * 50)
    
    # Create test suite
    suite = unittest.TestLoader().loadTestsFromTestCase(TestPDBDownloader)
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Print summary
    print("\nTest Summary:")
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print("-" * 50)
    
    return len(result.failures) + len(result.errors) == 0

if __name__ == "__main__":
    success = run_tests()
    exit(0 if success else 1) 