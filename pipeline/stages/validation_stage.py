#!/usr/bin/env python3

from typing import Dict, Any, Optional
from pathlib import Path
import sys

# Add utils to path for importing existing modules
sys.path.append(str(Path(__file__).parent.parent.parent))

from .base_stage import BaseStage, StageResult
from utils.analyze_pdb_quality import PDBAnalyzer

class ValidationStage(BaseStage):
    """Pipeline stage for validating and analyzing data quality."""
    
    def __init__(self, name: str, config_manager, experiment_id: str):
        """Initialize validation stage.
        
        Args:
            name: Stage name
            config_manager: Configuration manager instance
            experiment_id: Experiment identifier
        """
        super().__init__(name, config_manager, experiment_id)
        
        # Get quality configuration
        self.quality_config = config_manager.quality
        
        # Create validation-specific output directories
        self.validation_reports_dir = self.output_dir / "validation_reports"
        self.quality_plots_dir = self.output_dir / "quality_plots"
        self.validation_reports_dir.mkdir(parents=True, exist_ok=True)
        self.quality_plots_dir.mkdir(parents=True, exist_ok=True)
    
    def validate_inputs(self, input_data: Any = None) -> bool:
        """Validate validation stage inputs.
        
        Args:
            input_data: Dictionary containing directories to validate
            
        Returns:
            True if inputs are valid
        """
        if input_data is None:
            self.logger.error("Input data (directories to validate) is required")
            return False
        
        if not isinstance(input_data, dict):
            self.logger.error("Input data must be a dictionary of directories")
            return False
        
        # Check each directory
        for stage_name, directory in input_data.items():
            dir_path = Path(directory)
            if not dir_path.exists():
                self.logger.error(f"Directory for {stage_name} does not exist: {directory}")
                return False
        
        self.logger.info(f"Found {len(input_data)} directories to validate")
        return True
    
    def get_expected_outputs(self) -> Dict[str, Any]:
        """Get expected outputs for validation stage.
        
        Returns:
            Dictionary describing expected outputs
        """
        return {
            'validation_reports_directory': str(self.validation_reports_dir),
            'quality_plots_directory': str(self.quality_plots_dir),
            'expected_files': {
                'validation_report': "Comprehensive validation report",
                'quality_plots': "Data quality visualization plots",
                'statistics_summary': "Statistical summary of data quality"
            },
            'file_formats': ['json', 'txt', 'png', 'pdf'],
            'quality_thresholds': {
                'max_nan_percentage': self.quality_config.max_nan_percentage,
                'min_residues_per_structure': self.quality_config.min_residues_per_structure,
                'max_residues_per_structure': self.quality_config.max_residues_per_structure
            }
        }
    
    def execute(self, input_data: Any = None) -> StageResult:
        """Execute the validation stage.
        
        Args:
            input_data: Dictionary containing directories to validate
            
        Returns:
            StageResult with validation details
        """
        try:
            if input_data is None:
                raise ValueError("Input directories are required")
            
            self.logger.info("Starting data validation and quality analysis")
            
            validation_results = {}
            total_files_validated = 0
            
            # Validate each stage's outputs
            for stage_name, directory in input_data.items():
                self.logger.info(f"Validating {stage_name} outputs in {directory}")
                
                stage_validation = self._validate_stage_output(stage_name, Path(directory))
                validation_results[stage_name] = stage_validation
                total_files_validated += stage_validation.get('file_count', 0)
            
            # Perform quality analysis on NumPy arrays if available
            quality_analysis = {}
            if 'conversion' in input_data:
                quality_analysis = self._perform_quality_analysis(Path(input_data['conversion']))
            
            # Generate comprehensive validation report
            validation_report = self._generate_validation_report(validation_results, quality_analysis)
            
            # Save validation report
            report_file = self.validation_reports_dir / "validation_report.json"
            import json
            with open(report_file, 'w') as f:
                json.dump(validation_report, f, indent=2)
            
            # Create checkpoint
            self.create_checkpoint({
                'validation_results': validation_results,
                'quality_analysis': quality_analysis,
                'total_files_validated': total_files_validated,
                'validation_report_file': str(report_file)
            })
            
            # Determine overall validation success
            overall_success = self._determine_overall_success(validation_results, quality_analysis)
            
            if not overall_success:
                self.logger.warning("Validation completed with warnings or failures")
            else:
                self.logger.info("Validation completed successfully")
            
            return StageResult(
                success=overall_success,
                stage_name=self.name,
                execution_time=0,
                input_count=len(input_data),
                output_count=total_files_validated,
                metadata={
                    'validation_reports_directory': str(self.validation_reports_dir),
                    'quality_plots_directory': str(self.quality_plots_dir),
                    'validation_results': validation_results,
                    'quality_analysis': quality_analysis,
                    'total_files_validated': total_files_validated,
                    'validation_report_file': str(report_file)
                }
            )
            
        except Exception as e:
            error_msg = f"Validation stage failed: {str(e)}"
            self.logger.error(error_msg)
            
            return StageResult(
                success=False,
                stage_name=self.name,
                execution_time=0,
                input_count=0,
                output_count=0,
                error_message=error_msg
            )
    
    def _validate_stage_output(self, stage_name: str, directory: Path) -> Dict[str, Any]:
        """Validate outputs from a specific stage.
        
        Args:
            stage_name: Name of the stage
            directory: Directory containing stage outputs
            
        Returns:
            Dictionary with validation results
        """
        validation_result = {
            'stage_name': stage_name,
            'directory': str(directory),
            'exists': directory.exists(),
            'file_count': 0,
            'file_types': {},
            'issues': []
        }
        
        if not directory.exists():
            validation_result['issues'].append(f"Directory does not exist: {directory}")
            return validation_result
        
        # Count files by type
        for file_path in directory.rglob('*'):
            if file_path.is_file():
                validation_result['file_count'] += 1
                suffix = file_path.suffix.lower()
                validation_result['file_types'][suffix] = validation_result['file_types'].get(suffix, 0) + 1
        
        # Stage-specific validation
        if stage_name == 'download':
            pdb_files = list(directory.glob('*.pdb'))
            if len(pdb_files) == 0:
                validation_result['issues'].append("No PDB files found")
        
        elif stage_name == 'extraction':
            extracted_pdbs = list(directory.glob('*.pdb'))
            if len(extracted_pdbs) == 0:
                validation_result['issues'].append("No extracted PDB files found")
        
        elif stage_name == 'conversion':
            npy_files = list(directory.glob('*.npy'))
            if len(npy_files) == 0:
                validation_result['issues'].append("No NumPy files found")
        
        elif stage_name == 'loop_extraction':
            loop_files = list(directory.glob('*.npy'))
            # It's okay to have no loops if none were found
            validation_result['loop_count'] = len(loop_files)
        
        return validation_result
    
    def _perform_quality_analysis(self, numpy_dir: Path) -> Dict[str, Any]:
        """Perform quality analysis on NumPy arrays.
        
        Args:
            numpy_dir: Directory containing NumPy arrays
            
        Returns:
            Dictionary with quality analysis results
        """
        try:
            # Initialize analyzer
            analyzer = PDBAnalyzer(str(numpy_dir))
            
            # Analyze all files
            analyzer.analyze_all_files()
            
            # Get statistics
            stats = analyzer.get_statistics()
            
            # Check against quality thresholds
            quality_issues = []
            
            # Check NaN percentage
            if 'nan_percentage' in stats:
                if stats['nan_percentage'] > self.quality_config.max_nan_percentage:
                    quality_issues.append(f"High NaN percentage: {stats['nan_percentage']:.2f}% > {self.quality_config.max_nan_percentage}%")
            
            # Check residue counts
            if 'residue_stats' in stats:
                residue_stats = stats['residue_stats']
                if 'min_residues' in residue_stats:
                    if residue_stats['min_residues'] < self.quality_config.min_residues_per_structure:
                        quality_issues.append(f"Structure with too few residues: {residue_stats['min_residues']} < {self.quality_config.min_residues_per_structure}")
                
                if 'max_residues' in residue_stats:
                    if residue_stats['max_residues'] > self.quality_config.max_residues_per_structure:
                        quality_issues.append(f"Structure with too many residues: {residue_stats['max_residues']} > {self.quality_config.max_residues_per_structure}")
            
            return {
                'statistics': stats,
                'quality_issues': quality_issues,
                'quality_passed': len(quality_issues) == 0
            }
            
        except Exception as e:
            self.logger.error(f"Quality analysis failed: {e}")
            return {
                'statistics': {},
                'quality_issues': [f"Quality analysis failed: {str(e)}"],
                'quality_passed': False
            }
    
    def _generate_validation_report(self, validation_results: Dict, quality_analysis: Dict) -> Dict[str, Any]:
        """Generate comprehensive validation report.
        
        Args:
            validation_results: Results from stage validations
            quality_analysis: Results from quality analysis
            
        Returns:
            Comprehensive validation report
        """
        import datetime
        
        report = {
            'timestamp': datetime.datetime.now().isoformat(),
            'experiment_id': self.experiment_id,
            'validation_summary': {
                'total_stages_validated': len(validation_results),
                'stages_with_issues': sum(1 for result in validation_results.values() if result.get('issues')),
                'total_files_validated': sum(result.get('file_count', 0) for result in validation_results.values())
            },
            'stage_validations': validation_results,
            'quality_analysis': quality_analysis,
            'overall_status': 'PASSED' if self._determine_overall_success(validation_results, quality_analysis) else 'FAILED'
        }
        
        return report
    
    def _determine_overall_success(self, validation_results: Dict, quality_analysis: Dict) -> bool:
        """Determine overall validation success.
        
        Args:
            validation_results: Results from stage validations
            quality_analysis: Results from quality analysis
            
        Returns:
            True if validation passed overall
        """
        # Check if any stage has critical issues
        for result in validation_results.values():
            if not result.get('exists', False):
                return False
            if result.get('file_count', 0) == 0 and result['stage_name'] != 'loop_extraction':
                return False
        
        # Check quality analysis
        if quality_analysis and not quality_analysis.get('quality_passed', True):
            # Allow warnings but not critical failures
            critical_issues = [issue for issue in quality_analysis.get('quality_issues', []) 
                             if 'failed' in issue.lower()]
            if critical_issues:
                return False
        
        return True