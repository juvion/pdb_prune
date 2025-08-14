#!/usr/bin/env python3

import yaml
import os
from pathlib import Path
from typing import Dict, Any, Optional
import logging
from dataclasses import dataclass

@dataclass
class DataConfig:
    """Data directory configuration."""
    raw_pdbs_dir: str
    processed_dir: str
    experiments_dir: str
    cache_dir: str
    temp_dir: str

@dataclass
class ProcessingConfig:
    """Processing parameters configuration."""
    max_entries: int
    batch_size: int
    max_retries: int
    timeout_seconds: int
    delay_between_requests: float

@dataclass
class RNAConfig:
    """RNA structure parameters configuration."""
    distance_cutoff: float
    min_loop_length: int
    max_loop_length: int
    atom_types: list
    base_atoms: dict

@dataclass
class QualityConfig:
    """Quality control configuration."""
    max_nan_percentage: float
    min_sequence_length: int
    max_sequence_length: int
    min_residues_per_structure: int
    max_residues_per_structure: int
    resolution_cutoff: float

@dataclass
class LoggingConfig:
    """Logging configuration."""
    level: str
    format: str
    file_rotation: bool
    max_file_size: str
    backup_count: int

@dataclass
class PerformanceConfig:
    """Performance settings configuration."""
    chunk_size: int
    max_memory_usage: str
    parallel_workers: int
    enable_caching: bool

@dataclass
class ExperimentsConfig:
    """Experiment tracking configuration."""
    auto_versioning: bool
    track_metrics: bool
    save_checkpoints: bool
    checkpoint_frequency: int

class ConfigManager:
    """Manages pipeline configuration from YAML files."""
    
    def __init__(self, config_path: Optional[str] = None):
        """Initialize configuration manager.
        
        Args:
            config_path: Path to configuration file. If None, uses default.
        """
        if config_path is None:
            # Default to config file in same directory as this module
            config_path = Path(__file__).parent / "pipeline_config.yaml"
        
        self.config_path = Path(config_path)
        self._config_data = None
        self._load_config()
        self._setup_logging()
    
    def _load_config(self) -> None:
        """Load configuration from YAML file."""
        try:
            with open(self.config_path, 'r') as f:
                self._config_data = yaml.safe_load(f)
        except FileNotFoundError:
            raise FileNotFoundError(f"Configuration file not found: {self.config_path}")
        except yaml.YAMLError as e:
            raise ValueError(f"Invalid YAML configuration: {e}")
    
    def _setup_logging(self) -> None:
        """Setup logging based on configuration."""
        log_config = self.logging
        logging.basicConfig(
            level=getattr(logging, log_config.level),
            format=log_config.format
        )
    
    @property
    def data(self) -> DataConfig:
        """Get data configuration."""
        return DataConfig(**self._config_data['data'])
    
    @property
    def processing(self) -> ProcessingConfig:
        """Get processing configuration."""
        return ProcessingConfig(**self._config_data['processing'])
    
    @property
    def rna(self) -> RNAConfig:
        """Get RNA configuration."""
        return RNAConfig(**self._config_data['rna'])
    
    @property
    def quality(self) -> QualityConfig:
        """Get quality configuration."""
        return QualityConfig(**self._config_data['quality'])
    
    @property
    def logging(self) -> LoggingConfig:
        """Get logging configuration."""
        return LoggingConfig(**self._config_data['logging'])
    
    @property
    def performance(self) -> PerformanceConfig:
        """Get performance configuration."""
        return PerformanceConfig(**self._config_data['performance'])
    
    @property
    def experiments(self) -> ExperimentsConfig:
        """Get experiments configuration."""
        return ExperimentsConfig(**self._config_data['experiments'])
    
    def get_absolute_path(self, relative_path: str) -> Path:
        """Convert relative path to absolute path based on project root.
        
        Args:
            relative_path: Relative path from project root
            
        Returns:
            Absolute path
        """
        # Assume project root is parent of config directory
        project_root = self.config_path.parent.parent
        return project_root / relative_path
    
    def create_directories(self) -> None:
        """Create all configured directories if they don't exist."""
        data_config = self.data
        directories = [
            data_config.raw_pdbs_dir,
            data_config.processed_dir,
            data_config.experiments_dir,
            data_config.cache_dir,
            data_config.temp_dir
        ]
        
        for directory in directories:
            abs_path = self.get_absolute_path(directory)
            abs_path.mkdir(parents=True, exist_ok=True)
            logging.info(f"Created directory: {abs_path}")
    
    def validate_config(self) -> bool:
        """Validate configuration parameters.
        
        Returns:
            True if configuration is valid
            
        Raises:
            ValueError: If configuration is invalid
        """
        # Validate processing parameters
        proc_config = self.processing
        if proc_config.max_entries <= 0:
            raise ValueError("max_entries must be positive")
        if proc_config.batch_size <= 0:
            raise ValueError("batch_size must be positive")
        
        # Validate RNA parameters
        rna_config = self.rna
        if rna_config.distance_cutoff <= 0:
            raise ValueError("distance_cutoff must be positive")
        if rna_config.min_loop_length >= rna_config.max_loop_length:
            raise ValueError("min_loop_length must be less than max_loop_length")
        
        # Validate quality parameters
        quality_config = self.quality
        if not 0 <= quality_config.max_nan_percentage <= 100:
            raise ValueError("max_nan_percentage must be between 0 and 100")
        
        return True
    
    def update_config(self, section: str, key: str, value: Any) -> None:
        """Update configuration value.
        
        Args:
            section: Configuration section name
            key: Configuration key
            value: New value
        """
        if section not in self._config_data:
            raise KeyError(f"Configuration section '{section}' not found")
        
        self._config_data[section][key] = value
        logging.info(f"Updated config: {section}.{key} = {value}")
    
    def save_config(self, output_path: Optional[str] = None) -> None:
        """Save current configuration to file.
        
        Args:
            output_path: Path to save configuration. If None, overwrites current file.
        """
        if output_path is None:
            output_path = self.config_path
        
        with open(output_path, 'w') as f:
            yaml.dump(self._config_data, f, default_flow_style=False, indent=2)
        
        logging.info(f"Configuration saved to: {output_path}")