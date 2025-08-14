# RNA Processing Pipeline System

🧬 **A comprehensive, modular pipeline for processing RNA PDB structures with advanced quality control and loop extraction capabilities.**

## Overview

This pipeline system provides an end-to-end solution for processing RNA structures from the Protein Data Bank (PDB). It's designed for computational biologists and structural bioinformatics researchers who need to:

- Download and process large-scale RNA datasets
- Extract RNA chains and sequences from PDB files
- Convert structural data to NumPy arrays for analysis
- Identify and extract RNA loops and structural motifs
- Perform comprehensive quality control and validation

## Pipeline Architecture

### Core Components

```
pipeline/
├── orchestrator.py          # Main pipeline orchestration
├── stages/                  # Individual processing stages
│   ├── base_stage.py       # Abstract base class
│   ├── download_stage.py   # PDB file downloading
│   ├── extraction_stage.py # RNA chain extraction
│   ├── conversion_stage.py # PDB to NumPy conversion
│   ├── loop_extraction_stage.py # RNA loop identification
│   └── validation_stage.py # Quality control & validation
config/
├── config_manager.py       # Configuration management
└── pipeline_config.yaml    # Pipeline settings
```

### Processing Stages

1. **Download Stage** (`DownloadStage`)
   - Downloads RNA PDB files from RCSB PDB
   - Supports batch processing and resume functionality
   - Handles API rate limiting and error recovery

2. **Extraction Stage** (`ExtractionStage`)
   - Extracts RNA chains from PDB files
   - Generates FASTA sequences
   - Filters by chain type and quality

3. **Conversion Stage** (`ConversionStage`)
   - Converts PDB coordinates to NumPy arrays
   - Standardizes coordinate systems
   - Handles missing atoms and residues

4. **Loop Extraction Stage** (`LoopExtractionStage`)
   - Identifies RNA loops and structural motifs
   - Extracts loop coordinates and sequences
   - Classifies loop types and geometries

5. **Validation Stage** (`ValidationStage`)
   - Performs comprehensive quality control
   - Generates validation reports and plots
   - Checks data integrity and completeness

## Quick Start

### 1. Installation

```bash
# Clone the repository
git clone <repository-url>
cd pdb_prune

# Install dependencies
pip install -r requirements.txt
```

### 2. Configuration

Edit `config/pipeline_config.yaml` to customize processing parameters:

```yaml
processing:
  max_entries: 1000        # Number of PDB files to process
  batch_size: 100          # Batch size for processing
  
rna:
  min_loop_length: 3       # Minimum loop length
  max_loop_length: 20      # Maximum loop length
  distance_cutoff: 8.0     # Distance cutoff for contacts
  
quality:
  max_nan_percentage: 20.0 # Maximum allowed NaN percentage
  min_sequence_length: 10  # Minimum sequence length
  resolution_cutoff: 3.5   # Maximum resolution (Å)
```

### 3. Basic Usage

```python
from config.config_manager import ConfigManager
from pipeline.orchestrator import RNAPipeline

# Initialize pipeline
config = ConfigManager('config/pipeline_config.yaml')
pipeline = RNAPipeline(config)

# Create and run experiment
experiment_id = pipeline.create_experiment("my_rna_analysis")
result = pipeline.run_full_pipeline("my_rna_analysis", "RNA structure analysis")

print(f"Pipeline completed: {result.status}")
print(f"Stages completed: {result.stages_completed}/{result.total_stages}")
```

### 4. Demo Script

Run the demonstration script to explore pipeline features:

```bash
python demo_pipeline.py
```

## Advanced Usage

### Custom Stage Configuration

```python
# Setup pipeline with custom stages
pipeline.setup_stages(experiment_id)

# Run individual stages
download_stage = pipeline.stages[0]
result = download_stage.execute()

if result.success:
    print(f"Downloaded {result.output_count} files")
else:
    print(f"Error: {result.error_message}")
```

### Experiment Management

```python
# List all experiments
experiments = pipeline.list_experiments()
for exp in experiments:
    print(f"Experiment: {exp['experiment_id']}")

# Get experiment status
status = pipeline.get_experiment_status(experiment_id)
print(f"Directory: {status['directory']}")
print(f"Metadata: {status['metadata']}")
```

### Checkpoint and Resume

The pipeline supports automatic checkpointing:

```python
# Stages automatically create checkpoints
stage.create_checkpoint({
    'processed_files': 150,
    'current_batch': 2,
    'timestamp': time.time()
})

# Load checkpoint to resume processing
checkpoint = stage.load_checkpoint()
if checkpoint:
    start_from = checkpoint['processed_files']
```

## Output Structure

Each experiment creates a structured output directory:

```
data/experiments_data/{experiment_id}/
├── experiment_metadata.json    # Experiment information
├── pipeline_checkpoint.json    # Pipeline state
├── download/                   # Downloaded PDB files
│   ├── pdb_files/
│   ├── download_log.json
│   └── download_checkpoint.json
├── extraction/                 # Extracted RNA data
│   ├── rna_chains/
│   ├── sequences/
│   ├── extraction_log.json
│   └── extraction_checkpoint.json
├── conversion/                 # NumPy coordinate arrays
│   ├── numpy_arrays/
│   ├── conversion_log.json
│   └── conversion_checkpoint.json
├── loop_extraction/           # RNA loops and motifs
│   ├── loops/
│   ├── loop_sequences/
│   ├── loop_metadata/
│   └── loop_extraction_checkpoint.json
└── validation/                # Quality reports
    ├── reports/
    ├── plots/
    ├── validation_summary.json
    └── validation_checkpoint.json
```

## Quality Control

The validation stage performs comprehensive quality checks:

- **Coordinate Validation**: Checks for NaN values, outliers, and geometric consistency
- **Sequence Analysis**: Validates RNA sequences and chain completeness
- **Resolution Filtering**: Filters structures by experimental resolution
- **Statistical Analysis**: Generates distribution plots and summary statistics

### Quality Metrics

- NaN percentage in coordinates
- Sequence length distribution
- Resolution distribution
- Loop length statistics
- Structural completeness

## Configuration Reference

### Data Directories
```yaml
data:
  raw_pdbs_dir: "data/raw_data/pdbs"
  processed_dir: "data/processed"
  experiments_dir: "data/experiments_data"
  cache_dir: "data/cache"
  temp_dir: "data/temp"
```

### Processing Parameters
```yaml
processing:
  max_entries: 8904           # Total entries to process
  batch_size: 100             # Batch processing size
  resume_downloads: true      # Resume interrupted downloads
  parallel_processing: false  # Enable parallel processing
```

### RNA Structure Parameters
```yaml
rna:
  distance_cutoff: 8.0        # Contact distance cutoff (Å)
  min_loop_length: 3          # Minimum loop length
  max_loop_length: 20         # Maximum loop length
  include_modified: true      # Include modified nucleotides
```

### Quality Control
```yaml
quality:
  max_nan_percentage: 20.0    # Maximum NaN percentage
  min_sequence_length: 10     # Minimum sequence length
  max_sequence_length: 1000   # Maximum sequence length
  resolution_cutoff: 3.5      # Resolution cutoff (Å)
  min_residues_per_structure: 3
  max_residues_per_structure: 2000
```

## Testing

Run the comprehensive test suite:

```bash
# Test pipeline orchestration
python test_pipeline_orchestration.py

# Test individual components
python -m pytest tests/
```

## Performance Optimization

### Memory Management
- Configurable chunk sizes for large datasets
- Memory usage monitoring and limits
- Efficient NumPy array handling

### Parallel Processing
```yaml
performance:
  chunk_size: 1000
  max_memory_usage: "8GB"
  parallel_workers: 4
  enable_caching: true
```

### Caching
- Automatic caching of intermediate results
- Configurable cache size and retention
- Smart cache invalidation

## Troubleshooting

### Common Issues

1. **Memory Errors**: Reduce `batch_size` and `chunk_size`
2. **Download Failures**: Check internet connection and RCSB PDB availability
3. **Missing Dependencies**: Ensure all required packages are installed
4. **Permission Errors**: Check write permissions for output directories

### Logging

The pipeline provides detailed logging:

```python
# Configure logging level
logging.basicConfig(level=logging.INFO)

# View stage-specific logs
logger = logging.getLogger('pipeline.download')
logger.info("Custom log message")
```

### Debug Mode

```yaml
logging:
  level: "DEBUG"
  file: "pipeline.log"
  format: "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
```

## Contributing

### Adding Custom Stages

1. Inherit from `BaseStage`
2. Implement required abstract methods
3. Add to pipeline configuration

```python
from pipeline.stages.base_stage import BaseStage, StageResult

class CustomStage(BaseStage):
    def execute(self, input_data=None) -> StageResult:
        # Your custom processing logic
        return StageResult(
            success=True,
            stage_name=self.name,
            execution_time=elapsed_time,
            input_count=input_count,
            output_count=output_count
        )
    
    def validate_inputs(self, input_data=None) -> bool:
        # Input validation logic
        return True
    
    def get_expected_outputs(self) -> Dict[str, Any]:
        # Define expected outputs
        return {"custom_output": "description"}
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{rna_pipeline,
  title={RNA Processing Pipeline System},
  author={Your Name},
  year={2024},
  url={https://github.com/your-repo/rna-pipeline}
}
```

## Support

For questions, issues, or contributions:
- Open an issue on GitHub
- Contact the development team
- Check the documentation in `docs/`

---

**Built for computational structural biology research** 🧬