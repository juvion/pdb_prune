#!/usr/bin/env python3

from typing import Dict, Any, Optional
from pathlib import Path
import sys

# Add utils to path for importing existing modules
sys.path.append(str(Path(__file__).parent.parent.parent))

from .base_stage import BaseStage, StageResult
from utils.extract_loops import LoopExtractor

class LoopExtractionStage(BaseStage):
    """Pipeline stage for extracting RNA loops from NumPy arrays."""
    
    def __init__(self, name: str, config_manager, experiment_id: str):
        """Initialize loop extraction stage.
        
        Args:
            name: Stage name
            config_manager: Configuration manager instance
            experiment_id: Experiment identifier
        """
        super().__init__(name, config_manager, experiment_id)
        
        # Get RNA configuration
        self.rna_config = config_manager.rna
        
        # Create loop extraction-specific output directories
        self.loops_dir = self.output_dir / "extracted_loops"
        self.loops_dir.mkdir(parents=True, exist_ok=True)
    
    def validate_inputs(self, input_data: Any = None) -> bool:
        """Validate loop extraction stage inputs.
        
        Args:
            input_data: Input directory containing NumPy arrays
            
        Returns:
            True if inputs are valid
        """
        if input_data is None:
            self.logger.error("Input data (NumPy arrays directory) is required")
            return False
        
        input_dir = Path(input_data)
        if not input_dir.exists():
            self.logger.error(f"Input directory does not exist: {input_dir}")
            return False
        
        # Check for NumPy files
        npy_files = list(input_dir.glob("*.npy"))
        if len(npy_files) == 0:
            self.logger.error(f"No NumPy files found in input directory: {input_dir}")
            return False
        
        self.logger.info(f"Found {len(npy_files)} NumPy files for loop extraction")
        return True
    
    def get_expected_outputs(self) -> Dict[str, Any]:
        """Get expected outputs for loop extraction stage.
        
        Returns:
            Dictionary describing expected outputs
        """
        return {
            'loops_directory': str(self.loops_dir),
            'expected_files': {
                'loop_arrays': "Extracted RNA loop NumPy arrays"
            },
            'file_formats': ['npy'],
            'extraction_parameters': {
                'min_loop_length': self.rna_config.min_loop_length,
                'max_loop_length': self.rna_config.max_loop_length,
                'distance_cutoff': self.rna_config.distance_cutoff
            }
        }
    
    def execute(self, input_data: Any = None) -> StageResult:
        """Execute the loop extraction stage.
        
        Args:
            input_data: Input directory containing NumPy arrays
            
        Returns:
            StageResult with execution details
        """
        try:
            input_dir = Path(input_data) if input_data else None
            if input_dir is None:
                raise ValueError("Input directory is required")
            
            self.logger.info(f"Starting loop extraction from {input_dir}")
            
            # Count input files
            input_files = list(input_dir.glob("*.npy"))
            input_count = len(input_files)
            
            # Initialize loop extractor with configuration parameters
            extractor = LoopExtractor(
                input_dir=str(input_dir),
                output_dir=str(self.loops_dir),
                min_length=self.rna_config.min_loop_length,
                max_length=self.rna_config.max_loop_length,
                distance_cutoff=self.rna_config.distance_cutoff,
                atom_type=self.rna_config.atom_type
            )
            
            # Process all files
            extractor.process_all_files()
            
            # Count output files
            output_files = list(self.loops_dir.glob("*.npy"))
            output_count = len(output_files)
            
            # Create checkpoint
            self.create_checkpoint({
                'input_directory': str(input_dir),
                'loops_directory': str(self.loops_dir),
                'extracted_loops': output_count,
                'extraction_parameters': {
                    'min_loop_length': self.rna_config.min_loop_length,
                    'max_loop_length': self.rna_config.max_loop_length,
                    'distance_cutoff': self.rna_config.distance_cutoff,
                    'atom_type': self.rna_config.atom_type
                },
                'extraction_stats': extractor.stats if hasattr(extractor, 'stats') else {}
            })
            
            # Validate outputs
            if not self._validate_outputs():
                raise ValueError("Loop extraction validation failed")
            
            self.logger.info(f"Loop extraction completed. Files: {input_count} -> {output_count}")
            
            return StageResult(
                success=True,
                stage_name=self.name,
                execution_time=0,
                input_count=input_count,
                output_count=output_count,
                metadata={
                    'loops_directory': str(self.loops_dir),
                    'extracted_loops': output_count,
                    'extraction_parameters': {
                        'min_loop_length': self.rna_config.min_loop_length,
                        'max_loop_length': self.rna_config.max_loop_length,
                        'distance_cutoff': self.rna_config.distance_cutoff,
                        'atom_type': self.rna_config.atom_type
                    },
                    'extraction_stats': extractor.stats if hasattr(extractor, 'stats') else {}
                }
            )
            
        except Exception as e:
            error_msg = f"Loop extraction stage failed: {str(e)}"
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
        """Validate loop extraction outputs.
        
        Returns:
            True if outputs are valid
        """
        # Check if output directory exists
        if not self.loops_dir.exists():
            self.logger.error(f"Loops directory does not exist: {self.loops_dir}")
            return False
        
        # Check for extracted loop files
        loop_files = list(self.loops_dir.glob("*.npy"))
        
        if len(loop_files) == 0:
            self.logger.warning("No loop files found")
            return True  # This might be valid if no loops were found
        
        # Basic validation of loop files
        try:
            import numpy as np
            sample_file = loop_files[0]
            array = np.load(sample_file)
            self.logger.info(f"Sample loop array shape: {array.shape}")
            
            # Check if loop length is within expected range
            loop_length = array.shape[0] if len(array.shape) > 0 else 0
            if not (self.rna_config.min_loop_length <= loop_length <= self.rna_config.max_loop_length):
                self.logger.warning(f"Sample loop length {loop_length} outside expected range [{self.rna_config.min_loop_length}, {self.rna_config.max_loop_length}]")
            
        except Exception as e:
            self.logger.error(f"Failed to load sample loop file: {e}")
            return False
        
        self.logger.info(f"Loop extraction validation passed. Found {len(loop_files)} loop files")
        return True