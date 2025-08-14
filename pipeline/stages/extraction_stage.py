#!/usr/bin/env python3

from typing import Dict, Any, Optional
from pathlib import Path
import sys

# Add utils to path for importing existing modules
sys.path.append(str(Path(__file__).parent.parent.parent))

from .base_stage import BaseStage, StageResult
from utils.pdb_rna_processor import PDBRNAProcessor

class ExtractionStage(BaseStage):
    """Pipeline stage for extracting RNA chains from PDB files."""
    
    def __init__(self, name: str, config_manager, experiment_id: str):
        """Initialize extraction stage.
        
        Args:
            name: Stage name
            config_manager: Configuration manager instance
            experiment_id: Experiment identifier
        """
        super().__init__(name, config_manager, experiment_id)
        
        # Create extraction-specific output directories
        self.extracted_pdbs_dir = self.output_dir / "extracted_pdbs"
        self.extracted_sequences_dir = self.output_dir / "extracted_sequences"
        self.extracted_pdbs_dir.mkdir(parents=True, exist_ok=True)
        self.extracted_sequences_dir.mkdir(parents=True, exist_ok=True)
    
    def validate_inputs(self, input_data: Any = None) -> bool:
        """Validate extraction stage inputs.
        
        Args:
            input_data: Input directory containing raw PDB files
            
        Returns:
            True if inputs are valid
        """
        if input_data is None:
            self.logger.error("Input data (PDB directory) is required")
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
        
        self.logger.info(f"Found {len(pdb_files)} PDB files for extraction")
        return True
    
    def get_expected_outputs(self) -> Dict[str, Any]:
        """Get expected outputs for extraction stage.
        
        Returns:
            Dictionary describing expected outputs
        """
        return {
            'extracted_pdbs_directory': str(self.extracted_pdbs_dir),
            'extracted_sequences_directory': str(self.extracted_sequences_dir),
            'expected_files': {
                'extracted_pdb_files': "RNA chain PDB files",
                'fasta_files': "RNA sequence FASTA files"
            },
            'file_formats': ['pdb', 'fasta']
        }
    
    def execute(self, input_data: Any = None) -> StageResult:
        """Execute the extraction stage.
        
        Args:
            input_data: Input directory containing raw PDB files
            
        Returns:
            StageResult with execution details
        """
        try:
            input_dir = Path(input_data) if input_data else None
            if input_dir is None:
                raise ValueError("Input directory is required")
            
            self.logger.info(f"Starting RNA extraction from {input_dir}")
            
            # Count input files
            input_files = list(input_dir.glob("*.pdb"))
            input_count = len(input_files)
            
            # Initialize RNA processor
            processor = PDBRNAProcessor(
                original_dir=str(input_dir),
                processed_dir=str(self.extracted_pdbs_dir),
                fasta_dir=str(self.extracted_sequences_dir)
            )
            
            # Process all files
            processor.process_all_files()
            
            # Count output files
            extracted_pdbs = list(self.extracted_pdbs_dir.glob("*.pdb"))
            extracted_fastas = list(self.extracted_sequences_dir.glob("*.fasta"))
            output_count = len(extracted_pdbs)
            
            # Create checkpoint
            self.create_checkpoint({
                'input_directory': str(input_dir),
                'extracted_pdbs': len(extracted_pdbs),
                'extracted_sequences': len(extracted_fastas),
                'processing_stats': processor.stats
            })
            
            # Validate outputs
            if not self._validate_outputs():
                raise ValueError("Extraction validation failed")
            
            self.logger.info(f"Extraction completed. Files: {input_count} -> {output_count}")
            
            return StageResult(
                success=True,
                stage_name=self.name,
                execution_time=0,
                input_count=input_count,
                output_count=output_count,
                metadata={
                    'extracted_pdbs_directory': str(self.extracted_pdbs_dir),
                    'extracted_sequences_directory': str(self.extracted_sequences_dir),
                    'extracted_pdbs': len(extracted_pdbs),
                    'extracted_sequences': len(extracted_fastas),
                    'processing_stats': processor.stats
                }
            )
            
        except Exception as e:
            error_msg = f"Extraction stage failed: {str(e)}"
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
        """Validate extraction outputs.
        
        Returns:
            True if outputs are valid
        """
        # Check if output directories exist
        if not self.extracted_pdbs_dir.exists():
            self.logger.error(f"Extracted PDBs directory does not exist: {self.extracted_pdbs_dir}")
            return False
        
        if not self.extracted_sequences_dir.exists():
            self.logger.error(f"Extracted sequences directory does not exist: {self.extracted_sequences_dir}")
            return False
        
        # Check for extracted files
        extracted_pdbs = list(self.extracted_pdbs_dir.glob("*.pdb"))
        extracted_fastas = list(self.extracted_sequences_dir.glob("*.fasta"))
        
        if len(extracted_pdbs) == 0:
            self.logger.warning("No extracted PDB files found")
        
        if len(extracted_fastas) == 0:
            self.logger.warning("No extracted FASTA files found")
        
        self.logger.info(f"Extraction validation passed. Found {len(extracted_pdbs)} PDBs, {len(extracted_fastas)} FASTAs")
        return True