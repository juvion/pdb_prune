#!/usr/bin/env python3
"""
RNA Pipeline Demonstration Script

This script demonstrates how to use the RNA processing pipeline system
for downloading, extracting, converting, and analyzing RNA structures.
"""

import sys
from pathlib import Path
import logging

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from config.config_manager import ConfigManager
from pipeline.orchestrator import RNAPipeline

def demo_pipeline_usage():
    """Demonstrate basic pipeline usage."""
    
    print("🧬 RNA Pipeline Demonstration")
    print("=" * 50)
    
    # 1. Initialize the pipeline
    print("\n1. Initializing pipeline...")
    config = ConfigManager('config/pipeline_config.yaml')
    pipeline = RNAPipeline(config)
    print(f"✓ Pipeline initialized with {len(pipeline.stages)} stages")
    
    # 2. Create a new experiment
    print("\n2. Creating experiment...")
    experiment_id = pipeline.create_experiment("demo_run")
    print(f"✓ Experiment created: {experiment_id}")
    print(f"  - Directory: {pipeline.experiments_dir}/{experiment_id}")
    
    # 3. Show available stages
    print("\n3. Available pipeline stages:")
    for stage_name, stage_class in pipeline.STANDARD_STAGES:
        print(f"  - {stage_name}: {stage_class.__name__}")
    
    # 4. Setup stages for the experiment
    print("\n4. Setting up stages...")
    pipeline.setup_stages(experiment_id)
    for stage in pipeline.stages:
        print(f"✓ {stage.name}: {stage.__class__.__name__}")
        print(f"  - Output dir: {stage.output_dir}")
    
    # 5. Show configuration
    print("\n5. Pipeline configuration:")
    print(f"  - Max entries: {config.processing.max_entries}")
    print(f"  - Batch size: {config.processing.batch_size}")
    print(f"  - Min loop length: {config.rna.min_loop_length}")
    print(f"  - Max loop length: {config.rna.max_loop_length}")
    print(f"  - Quality threshold: {config.quality.max_nan_percentage}%")
    
    # 6. Show experiment status
    print("\n6. Experiment status:")
    status = pipeline.get_experiment_status(experiment_id)
    print(f"  - Experiment ID: {status['experiment_id']}")
    print(f"  - Directory exists: {Path(status['directory']).exists()}")
    
    # 7. List all experiments
    print("\n7. All experiments:")
    experiments = pipeline.list_experiments()
    print(f"  - Total experiments: {len(experiments)}")
    for exp in experiments[-3:]:  # Show last 3
        print(f"    • {exp['experiment_id']}")
    
    print("\n" + "=" * 50)
    print("🎉 Pipeline demonstration completed!")
    print("\nTo run the actual pipeline:")
    print("  python -m pipeline.orchestrator --experiment demo_run --max-entries 10")
    
    return True

def show_pipeline_help():
    """Show help information for using the pipeline."""
    
    help_text = """
🧬 RNA Pipeline System - Usage Guide
=" * 50

The RNA pipeline processes PDB structures through these stages:

1. DOWNLOAD: Downloads RNA PDB files from RCSB
2. EXTRACTION: Extracts RNA chains and sequences
3. CONVERSION: Converts PDB to NumPy coordinate arrays
4. LOOP_EXTRACTION: Identifies and extracts RNA loops
5. VALIDATION: Performs quality control and generates reports

Basic Usage:
-----------

# Initialize pipeline
from config.config_manager import ConfigManager
from pipeline.orchestrator import RNAPipeline

config = ConfigManager('config/pipeline_config.yaml')
pipeline = RNAPipeline(config)

# Create experiment
experiment_id = pipeline.create_experiment("my_experiment")

# Initialize stages
stages = pipeline.initialize_stages(experiment_id)

# Run individual stages
download_stage = stages['download']
result = download_stage.execute()

Configuration:
-------------
Edit config/pipeline_config.yaml to customize:
- max_entries: Number of PDB files to process
- batch_size: Processing batch size
- min_loop_length: Minimum RNA loop length
- quality thresholds: Data quality criteria

Output Structure:
----------------
data/experiments_data/{experiment_id}/
├── download/           # Downloaded PDB files
├── extraction/         # Extracted RNA data
├── conversion/         # NumPy coordinate arrays
├── loop_extraction/    # Extracted RNA loops
└── validation/         # Quality reports and plots

For more information, see the documentation in docs/
"""
    
    print(help_text)

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--help":
        show_pipeline_help()
    else:
        try:
            success = demo_pipeline_usage()
            sys.exit(0 if success else 1)
        except Exception as e:
            print(f"❌ Demo failed: {e}")
            sys.exit(1)