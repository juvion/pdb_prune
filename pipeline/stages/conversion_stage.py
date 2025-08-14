#!/usr/bin/env python3

from typing import Dict, Any, Optional
from pathlib import Path
import sys

# Add utils to path for importing existing modules
sys.path.append(str(Path(__file__).parent.parent.parent))

from .base_stage import BaseStage, StageResult
from utils.pdb_to_npy import PDBToNumpyConverter

class ConversionStage(BaseStage):
    """Pipeline stage for converting PDB files to NumPy arrays."""
    
    def __init__(self, name: str, config_manager, experiment_id: str):
        """Initialize conversion stage.
        
        Args:
            name: Stage name
            config_manager: Configuration manager instance
            experiment_id: Experiment identifier
        """
        super().__init__(name, config_manager, experiment_id)
        
        # Create conversion-specific output directories
        self.numpy_arrays_dir = self.output_dir / "numpy_arrays"
        self.numpy_arrays_dir.mkdir(parents=True, exist_ok=True)
    
    def validate_inputs(self, input_data: Any = None) -> bool:
        """Validate conversion stage inputs.
        
        Args:
            input_data: Input directory containing extracted PDB files
            
        Returns:
            True if inputs are valid
        """
        if input_data is None:
            self.logger.error("Input data (extracted PDB directory) is required")
            return False
        
        input_dir = Path(input_data)
        if not input_dir.exists():
            self.logger.error(f"Input directory does not exist: {input_dir}")
            return False
        
        # Check for PDB files
        pdb_files = list(input_dir.glob("*.pdb"))
        if len(pdb_files) == 0:
            self.logger.error(f"No PDB files found in input directory: {input_dir}")
            return False
        
        self.logger.info(f"Found {len(pdb_files)} PDB files for conversion")
        return True
    
    def get_expected_outputs(self) -> Dict[str, Any]:
        """Get expected outputs for conversion stage.
        
        Returns:
            Dictionary describing expected outputs
        """
        return {
            'numpy_arrays_directory': str(self.numpy_arrays_dir),
            'expected_files': {
                'numpy_arrays': "RNA coordinate NumPy arrays"
            },
            'file_formats': ['npy']
        }
    
    def execute(self, input_data: Any = None) -> StageResult:
        """Execute the conversion stage.
        
        Args:
            input_data: Input directory containing extracted PDB files
            
        Returns:
            StageResult with execution details
        """
        try:
            input_dir = Path(input_data) if input_data else None
            if input_dir is None:
                raise ValueError("Input directory is required")
            
            self.logger.info(f"Starting PDB to NumPy conversion from {input_dir}")
            
            # Count input files
            input_files = list(input_dir.glob("*.pdb"))
            input_count = len(input_files)
            
            # Initialize converter
            converter = PDBToNumpyConverter(
                input_dir=str(input_dir),
                output_dir=str(self.numpy_arrays_dir)
            )
            
            # Process all files
            converter.process_all_files()
            
            # Count output files
            output_files = list(self.numpy_arrays_dir.glob("*.npy"))
            output_count = len(output_files)
            
            # Create checkpoint
            self.create_checkpoint({
                'input_directory': str(input_dir),
                'numpy_arrays_directory': str(self.numpy_arrays_dir),
                'converted_files': output_count,
                'conversion_stats': converter.stats if hasattr(converter, 'stats') else {}
            })
            
            # Validate outputs
            if not self._validate_outputs():
                raise ValueError("Conversion validation failed")
            
            self.logger.info(f"Conversion completed. Files: {input_count} -> {output_count}")
            
            return StageResult(
                success=True,
                stage_name=self.name,
                execution_time=0,
                input_count=input_count,
                output_count=output_count,
                metadata={
                    'numpy_arrays_directory': str(self.numpy_arrays_dir),
                    'converted_files': output_count,
                    'conversion_stats': converter.stats if hasattr(converter, 'stats') else {}
                }
            )
            
        except Exception as e:
            error_msg = f"Conversion stage failed: {str(e)}"
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
        """Validate conversion outputs.
        
        Returns:
            True if outputs are valid
        """
        # Check if output directory exists
        if not self.numpy_arrays_dir.exists():
            self.logger.error(f"NumPy arrays directory does not exist: {self.numpy_arrays_dir}")
            return False
        
        # Check for converted files
        numpy_files = list(self.numpy_arrays_dir.glob("*.npy"))
        
        if len(numpy_files) == 0:
            self.logger.warning("No NumPy array files found")
            return False
        
        # Basic validation of NumPy files
        try:
            import numpy as np
            sample_file = numpy_files[0]
            array = np.load(sample_file)
            self.logger.info(f"Sample array shape: {array.shape}")
        except Exception as e:
            self.logger.error(f"Failed to load sample NumPy file: {e}")
            return False
        
        self.logger.info(f"Conversion validation passed. Found {len(numpy_files)} NumPy files")
        return True