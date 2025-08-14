#!/usr/bin/env python3

from typing import Dict, Any, Optional
from pathlib import Path
import sys
import os

# Add utils to path for importing existing modules
sys.path.append(str(Path(__file__).parent.parent.parent))

from .base_stage import BaseStage, StageResult
from utils.pdb_downloader import download_rna_pdbs

class DownloadStage(BaseStage):
    """Pipeline stage for downloading RNA PDB files."""
    
    def __init__(self, name: str, config_manager, experiment_id: str):
        """Initialize download stage.
        
        Args:
            name: Stage name
            config_manager: Configuration manager instance
            experiment_id: Experiment identifier
        """
        super().__init__(name, config_manager, experiment_id)
        
        # Create download-specific output directory
        self.download_dir = self.output_dir / "raw_pdbs"
        self.download_dir.mkdir(parents=True, exist_ok=True)
    
    def validate_inputs(self, input_data: Any = None) -> bool:
        """Validate download stage inputs.
        
        Args:
            input_data: Input data (not used for download stage)
            
        Returns:
            True if inputs are valid
        """
        # Download stage doesn't require input data
        # Validate configuration parameters
        proc_config = self.config.processing
        
        if proc_config.max_entries <= 0:
            self.logger.error("max_entries must be positive")
            return False
        
        if proc_config.batch_size <= 0:
            self.logger.error("batch_size must be positive")
            return False
        
        return True
    
    def get_expected_outputs(self) -> Dict[str, Any]:
        """Get expected outputs for download stage.
        
        Returns:
            Dictionary describing expected outputs
        """
        return {
            'output_directory': str(self.download_dir),
            'expected_files': {
                'pdb_files': f"Up to {self.config.processing.max_entries} .pdb files",
                'found_pdb_ids.txt': "List of discovered PDB IDs",
                'downloaded_pdb_ids_*.txt': "List of successfully downloaded PDB IDs"
            },
            'file_formats': ['pdb', 'txt']
        }
    
    def execute(self, input_data: Any = None) -> StageResult:
        """Execute the download stage.
        
        Args:
            input_data: Input data (not used for download stage)
            
        Returns:
            StageResult with execution details
        """
        try:
            self.logger.info(f"Starting RNA PDB download to {self.download_dir}")
            
            # Get configuration parameters
            proc_config = self.config.processing
            
            # Check if we should resume from checkpoint
            checkpoint_data = self.load_checkpoint()
            start_from = 0
            search_start_page = 1
            
            if checkpoint_data:
                start_from = checkpoint_data.get('last_downloaded_index', 0)
                search_start_page = checkpoint_data.get('last_search_page', 1)
                self.logger.info(f"Resuming download from index {start_from}, page {search_start_page}")
            
            # Count existing files before download
            existing_files = list(self.download_dir.glob("*.pdb"))
            input_count = len(existing_files)
            
            # Call the existing download function
            download_rna_pdbs(
                download_directory=str(self.download_dir),
                max_entries=proc_config.max_entries,
                start_from=start_from,
                batch_size=proc_config.batch_size,
                max_retries=proc_config.max_retries,
                search_start_page=search_start_page
            )
            
            # Count files after download
            final_files = list(self.download_dir.glob("*.pdb"))
            output_count = len(final_files)
            
            # Create checkpoint with current progress
            self.create_checkpoint({
                'last_downloaded_index': output_count,
                'last_search_page': search_start_page + 1,
                'total_files': output_count,
                'download_directory': str(self.download_dir)
            })
            
            # Validate outputs
            if not self._validate_outputs():
                raise ValueError("Download validation failed")
            
            self.logger.info(f"Download completed. Files: {input_count} -> {output_count}")
            
            return StageResult(
                success=True,
                stage_name=self.name,
                execution_time=0,  # Will be set by base class
                input_count=input_count,
                output_count=output_count,
                metadata={
                    'download_directory': str(self.download_dir),
                    'files_downloaded': output_count - input_count,
                    'total_files': output_count
                }
            )
            
        except Exception as e:
            error_msg = f"Download stage failed: {str(e)}"
            self.logger.error(error_msg)
            
            return StageResult(
                success=False,
                stage_name=self.name,
                execution_time=0,
                input_count=0,
                output_count=0,
                error_message=error_msg
            )
    
    def _validate_outputs(self) -> bool:
        """Validate download outputs.
        
        Returns:
            True if outputs are valid
        """
        # Check if download directory exists
        if not self.download_dir.exists():
            self.logger.error(f"Download directory does not exist: {self.download_dir}")
            return False
        
        # Check for PDB files
        pdb_files = list(self.download_dir.glob("*.pdb"))
        if len(pdb_files) == 0:
            self.logger.warning("No PDB files found after download")
            # This might be valid if resuming and all files already exist
        
        # Check for found_pdb_ids.txt
        found_ids_file = self.download_dir / "found_pdb_ids.txt"
        if not found_ids_file.exists():
            self.logger.error("found_pdb_ids.txt not found")
            return False
        
        # Validate file sizes
        for pdb_file in pdb_files[:10]:  # Check first 10 files
            if pdb_file.stat().st_size == 0:
                self.logger.warning(f"Empty PDB file found: {pdb_file}")
        
        self.logger.info(f"Download validation passed. Found {len(pdb_files)} PDB files")
        return True
    
    def get_output_summary(self) -> Dict[str, Any]:
        """Get summary of download outputs.
        
        Returns:
            Dictionary with output summary
        """
        pdb_files = list(self.download_dir.glob("*.pdb"))
        
        # Read found PDB IDs if available
        found_ids_file = self.download_dir / "found_pdb_ids.txt"
        found_ids_count = 0
        if found_ids_file.exists():
            try:
                with open(found_ids_file, 'r') as f:
                    found_ids_count = len(f.readlines())
            except Exception as e:
                self.logger.warning(f"Could not read found_pdb_ids.txt: {e}")
        
        return {
            'pdb_files_downloaded': len(pdb_files),
            'found_pdb_ids': found_ids_count,
            'download_directory': str(self.download_dir),
            'largest_file_size': max((f.stat().st_size for f in pdb_files), default=0),
            'smallest_file_size': min((f.stat().st_size for f in pdb_files), default=0),
            'total_size_mb': sum(f.stat().st_size for f in pdb_files) / (1024 * 1024)
        }