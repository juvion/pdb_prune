#!/usr/bin/env python3

from abc import ABC, abstractmethod
from typing import Dict, Any, Optional
from pathlib import Path
import logging
import time
from dataclasses import dataclass

@dataclass
class StageResult:
    """Result of a pipeline stage execution."""
    success: bool
    stage_name: str
    execution_time: float
    input_count: int
    output_count: int
    error_message: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None

class BaseStage(ABC):
    """Base class for all pipeline stages."""
    
    def __init__(self, name: str, config_manager, experiment_id: str):
        """Initialize base stage.
        
        Args:
            name: Stage name
            config_manager: Configuration manager instance
            experiment_id: Experiment identifier
        """
        self.name = name
        self.config = config_manager
        self.experiment_id = experiment_id
        self.logger = logging.getLogger(f"pipeline.{name}")
        
        # Create experiment-specific directories
        self.experiment_dir = self.config.get_absolute_path(
            f"{self.config.data.experiments_dir}/{experiment_id}"
        )
        self.experiment_dir.mkdir(parents=True, exist_ok=True)
        
        # Stage-specific output directory
        self.output_dir = self.experiment_dir / name
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    @abstractmethod
    def execute(self, input_data: Any = None) -> StageResult:
        """Execute the pipeline stage.
        
        Args:
            input_data: Input data from previous stage
            
        Returns:
            StageResult with execution details
        """
        pass
    
    @abstractmethod
    def validate_inputs(self, input_data: Any = None) -> bool:
        """Validate stage inputs.
        
        Args:
            input_data: Input data to validate
            
        Returns:
            True if inputs are valid
        """
        pass
    
    @abstractmethod
    def get_expected_outputs(self) -> Dict[str, Any]:
        """Get expected outputs for this stage.
        
        Returns:
            Dictionary describing expected outputs
        """
        pass
    
    def run(self, input_data: Any = None) -> StageResult:
        """Run the stage with error handling and timing.
        
        Args:
            input_data: Input data from previous stage
            
        Returns:
            StageResult with execution details
        """
        start_time = time.time()
        
        try:
            self.logger.info(f"Starting stage: {self.name}")
            
            # Validate inputs
            if not self.validate_inputs(input_data):
                raise ValueError(f"Invalid inputs for stage {self.name}")
            
            # Execute stage
            result = self.execute(input_data)
            
            # Update timing
            result.execution_time = time.time() - start_time
            
            if result.success:
                self.logger.info(
                    f"Stage {self.name} completed successfully in {result.execution_time:.2f}s. "
                    f"Input: {result.input_count}, Output: {result.output_count}"
                )
            else:
                self.logger.error(f"Stage {self.name} failed: {result.error_message}")
            
            return result
            
        except Exception as e:
            execution_time = time.time() - start_time
            error_msg = f"Stage {self.name} failed with exception: {str(e)}"
            self.logger.error(error_msg)
            
            return StageResult(
                success=False,
                stage_name=self.name,
                execution_time=execution_time,
                input_count=0,
                output_count=0,
                error_message=error_msg
            )
    
    def create_checkpoint(self, data: Dict[str, Any]) -> Path:
        """Create a checkpoint file for this stage.
        
        Args:
            data: Data to save in checkpoint
            
        Returns:
            Path to checkpoint file
        """
        import json
        
        checkpoint_file = self.output_dir / f"{self.name}_checkpoint.json"
        
        with open(checkpoint_file, 'w') as f:
            json.dump({
                'stage_name': self.name,
                'experiment_id': self.experiment_id,
                'timestamp': time.time(),
                'data': data
            }, f, indent=2)
        
        self.logger.info(f"Checkpoint created: {checkpoint_file}")
        return checkpoint_file
    
    def load_checkpoint(self) -> Optional[Dict[str, Any]]:
        """Load checkpoint data for this stage.
        
        Returns:
            Checkpoint data if exists, None otherwise
        """
        import json
        
        checkpoint_file = self.output_dir / f"{self.name}_checkpoint.json"
        
        if checkpoint_file.exists():
            try:
                with open(checkpoint_file, 'r') as f:
                    checkpoint_data = json.load(f)
                self.logger.info(f"Checkpoint loaded: {checkpoint_file}")
                return checkpoint_data.get('data')
            except Exception as e:
                self.logger.warning(f"Failed to load checkpoint: {e}")
        
        return None
    
    def cleanup(self) -> None:
        """Clean up temporary files and resources."""
        # Default implementation - can be overridden by subclasses
        temp_files = list(self.output_dir.glob("*.tmp"))
        for temp_file in temp_files:
            try:
                temp_file.unlink()
                self.logger.debug(f"Cleaned up temporary file: {temp_file}")
            except Exception as e:
                self.logger.warning(f"Failed to clean up {temp_file}: {e}")