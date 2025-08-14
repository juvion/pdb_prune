#!/usr/bin/env python3

from typing import List, Dict, Any, Optional, Type
from pathlib import Path
import logging
import time
import json
from dataclasses import dataclass, asdict
from enum import Enum

from config.config_manager import ConfigManager
from .stages.base_stage import BaseStage, StageResult
from .stages.download_stage import DownloadStage
from .stages.extraction_stage import ExtractionStage
from .stages.conversion_stage import ConversionStage
from .stages.loop_extraction_stage import LoopExtractionStage
from .stages.validation_stage import ValidationStage

class PipelineStatus(Enum):
    """Pipeline execution status."""
    NOT_STARTED = "not_started"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    PAUSED = "paused"

@dataclass
class PipelineResult:
    """Result of pipeline execution."""
    experiment_id: str
    status: PipelineStatus
    total_execution_time: float
    stages_completed: int
    total_stages: int
    stage_results: List[StageResult]
    error_message: Optional[str] = None

class RNAPipeline:
    """Main pipeline orchestrator for RNA structure processing."""
    
    # Define the standard pipeline stages in order
    STANDARD_STAGES = [
        ('download', DownloadStage),
        ('extraction', ExtractionStage), 
        ('conversion', ConversionStage),
        ('loop_extraction', LoopExtractionStage),
        ('validation', ValidationStage)
    ]
    
    def __init__(self, config_manager_or_path = None):
        """Initialize the RNA processing pipeline.
        
        Args:
            config_manager_or_path: Either a ConfigManager instance or path to configuration file
        """
        if isinstance(config_manager_or_path, ConfigManager):
            self.config = config_manager_or_path
        else:
            self.config = ConfigManager(config_manager_or_path)
        self.logger = logging.getLogger("pipeline.orchestrator")
        
        # Validate configuration
        self.config.validate_config()
        self.config.create_directories()
        
        # Pipeline state
        self.current_experiment_id: Optional[str] = None
        self.stages: List[BaseStage] = []
        self.stage_results: List[StageResult] = []
        
    @property
    def experiments_dir(self) -> Path:
        """Get the experiments directory path."""
        return self.config.get_absolute_path(self.config.data.experiments_dir)
        
    def create_experiment(self, experiment_name: str, description: str = "") -> str:
        """Create a new experiment.
        
        Args:
            experiment_name: Name of the experiment
            description: Optional description
            
        Returns:
            Experiment ID
        """
        # Generate experiment ID with timestamp
        timestamp = int(time.time())
        experiment_id = f"{experiment_name}_{timestamp}"
        
        # Create experiment directory
        experiment_dir = self.config.get_absolute_path(
            f"{self.config.data.experiments_dir}/{experiment_id}"
        )
        experiment_dir.mkdir(parents=True, exist_ok=True)
        
        # Save experiment metadata
        metadata = {
            'experiment_id': experiment_id,
            'name': experiment_name,
            'description': description,
            'created_at': timestamp,
            'config': asdict(self.config.data)
        }
        
        metadata_file = experiment_dir / 'experiment_metadata.json'
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        self.logger.info(f"Created experiment: {experiment_id}")
        return experiment_id
    
    def setup_stages(self, experiment_id: str, custom_stages: Optional[List[tuple]] = None) -> None:
        """Setup pipeline stages for an experiment.
        
        Args:
            experiment_id: Experiment identifier
            custom_stages: Optional custom stages list, uses standard if None
        """
        self.current_experiment_id = experiment_id
        self.stages = []
        
        stages_to_use = custom_stages or self.STANDARD_STAGES
        
        for stage_name, stage_class in stages_to_use:
            stage = stage_class(stage_name, self.config, experiment_id)
            self.stages.append(stage)
            self.logger.info(f"Setup stage: {stage_name}")
    
    def run_full_pipeline(self, experiment_name: str, description: str = "") -> PipelineResult:
        """Run the complete pipeline.
        
        Args:
            experiment_name: Name of the experiment
            description: Optional description
            
        Returns:
            PipelineResult with execution details
        """
        start_time = time.time()
        
        try:
            # Create experiment
            experiment_id = self.create_experiment(experiment_name, description)
            
            # Setup stages
            self.setup_stages(experiment_id)
            
            # Run all stages
            return self._execute_stages(start_time)
            
        except Exception as e:
            error_msg = f"Pipeline failed: {str(e)}"
            self.logger.error(error_msg)
            
            return PipelineResult(
                experiment_id=self.current_experiment_id or "unknown",
                status=PipelineStatus.FAILED,
                total_execution_time=time.time() - start_time,
                stages_completed=len(self.stage_results),
                total_stages=len(self.stages),
                stage_results=self.stage_results,
                error_message=error_msg
            )
    
    def run_stage(self, stage_name: str, experiment_id: str, input_data: Any = None) -> StageResult:
        """Run a specific pipeline stage.
        
        Args:
            stage_name: Name of the stage to run
            experiment_id: Experiment identifier
            input_data: Input data for the stage
            
        Returns:
            StageResult with execution details
        """
        # Find stage class
        stage_class = None
        for name, cls in self.STANDARD_STAGES:
            if name == stage_name:
                stage_class = cls
                break
        
        if stage_class is None:
            raise ValueError(f"Unknown stage: {stage_name}")
        
        # Create and run stage
        stage = stage_class(stage_name, self.config, experiment_id)
        return stage.run(input_data)
    
    def resume_pipeline(self, experiment_id: str, from_stage: str) -> PipelineResult:
        """Resume pipeline execution from a specific stage.
        
        Args:
            experiment_id: Experiment identifier
            from_stage: Stage name to resume from
            
        Returns:
            PipelineResult with execution details
        """
        start_time = time.time()
        
        try:
            # Setup stages
            self.setup_stages(experiment_id)
            
            # Find starting stage index
            start_index = 0
            for i, stage in enumerate(self.stages):
                if stage.name == from_stage:
                    start_index = i
                    break
            else:
                raise ValueError(f"Stage '{from_stage}' not found")
            
            # Load previous results
            self._load_previous_results(experiment_id, start_index)
            
            # Execute remaining stages
            return self._execute_stages(start_time, start_index)
            
        except Exception as e:
            error_msg = f"Pipeline resume failed: {str(e)}"
            self.logger.error(error_msg)
            
            return PipelineResult(
                experiment_id=experiment_id,
                status=PipelineStatus.FAILED,
                total_execution_time=time.time() - start_time,
                stages_completed=len(self.stage_results),
                total_stages=len(self.stages),
                stage_results=self.stage_results,
                error_message=error_msg
            )
    
    def _execute_stages(self, start_time: float, start_index: int = 0) -> PipelineResult:
        """Execute pipeline stages.
        
        Args:
            start_time: Pipeline start time
            start_index: Index of first stage to execute
            
        Returns:
            PipelineResult with execution details
        """
        input_data = None
        
        for i in range(start_index, len(self.stages)):
            stage = self.stages[i]
            
            # Execute stage
            result = stage.run(input_data)
            self.stage_results.append(result)
            
            # Check if stage failed
            if not result.success:
                return PipelineResult(
                    experiment_id=self.current_experiment_id,
                    status=PipelineStatus.FAILED,
                    total_execution_time=time.time() - start_time,
                    stages_completed=i,
                    total_stages=len(self.stages),
                    stage_results=self.stage_results,
                    error_message=result.error_message
                )
            
            # Save checkpoint if configured
            if self.config.experiments.save_checkpoints:
                self._save_pipeline_checkpoint()
            
            # Prepare input for next stage
            input_data = self._prepare_next_stage_input(stage, result)
        
        # All stages completed successfully
        return PipelineResult(
            experiment_id=self.current_experiment_id,
            status=PipelineStatus.COMPLETED,
            total_execution_time=time.time() - start_time,
            stages_completed=len(self.stages),
            total_stages=len(self.stages),
            stage_results=self.stage_results
        )
    
    def _prepare_next_stage_input(self, completed_stage: BaseStage, result: StageResult) -> Any:
        """Prepare input data for the next stage.
        
        Args:
            completed_stage: The stage that just completed
            result: Result from the completed stage
            
        Returns:
            Input data for next stage
        """
        # Default implementation - can be customized based on stage types
        return completed_stage.output_dir
    
    def _save_pipeline_checkpoint(self) -> None:
        """Save pipeline checkpoint."""
        if not self.current_experiment_id:
            return
        
        experiment_dir = self.config.get_absolute_path(
            f"{self.config.data.experiments_dir}/{self.current_experiment_id}"
        )
        
        checkpoint_data = {
            'experiment_id': self.current_experiment_id,
            'timestamp': time.time(),
            'stages_completed': len(self.stage_results),
            'total_stages': len(self.stages),
            'stage_results': [asdict(result) for result in self.stage_results]
        }
        
        checkpoint_file = experiment_dir / 'pipeline_checkpoint.json'
        with open(checkpoint_file, 'w') as f:
            json.dump(checkpoint_data, f, indent=2)
    
    def _load_previous_results(self, experiment_id: str, start_index: int) -> None:
        """Load previous stage results.
        
        Args:
            experiment_id: Experiment identifier
            start_index: Index to start loading from
        """
        experiment_dir = self.config.get_absolute_path(
            f"{self.config.data.experiments_dir}/{experiment_id}"
        )
        
        checkpoint_file = experiment_dir / 'pipeline_checkpoint.json'
        if checkpoint_file.exists():
            try:
                with open(checkpoint_file, 'r') as f:
                    checkpoint_data = json.load(f)
                
                # Load stage results up to start_index
                stage_results_data = checkpoint_data.get('stage_results', [])
                for i in range(min(start_index, len(stage_results_data))):
                    result_data = stage_results_data[i]
                    # Convert back to StageResult object
                    result = StageResult(**result_data)
                    self.stage_results.append(result)
                
                self.logger.info(f"Loaded {len(self.stage_results)} previous stage results")
                
            except Exception as e:
                self.logger.warning(f"Failed to load previous results: {e}")
    
    def get_experiment_status(self, experiment_id: str) -> Dict[str, Any]:
        """Get status of an experiment.
        
        Args:
            experiment_id: Experiment identifier
            
        Returns:
            Dictionary with experiment status
        """
        experiment_dir = self.config.get_absolute_path(
            f"{self.config.data.experiments_dir}/{experiment_id}"
        )
        
        if not experiment_dir.exists():
            return {'error': f'Experiment {experiment_id} not found'}
        
        # Load metadata
        metadata_file = experiment_dir / 'experiment_metadata.json'
        metadata = {}
        if metadata_file.exists():
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
        
        # Load checkpoint if exists
        checkpoint_file = experiment_dir / 'pipeline_checkpoint.json'
        checkpoint = {}
        if checkpoint_file.exists():
            with open(checkpoint_file, 'r') as f:
                checkpoint = json.load(f)
        
        return {
            'experiment_id': experiment_id,
            'metadata': metadata,
            'checkpoint': checkpoint,
            'directory': str(experiment_dir)
        }
    
    def list_experiments(self) -> List[Dict[str, Any]]:
        """List all experiments.
        
        Returns:
            List of experiment information
        """
        experiments_dir = self.config.get_absolute_path(self.config.data.experiments_dir)
        experiments = []
        
        if experiments_dir.exists():
            for exp_dir in experiments_dir.iterdir():
                if exp_dir.is_dir():
                    status = self.get_experiment_status(exp_dir.name)
                    if 'error' not in status:
                        experiments.append(status)
        
        return experiments