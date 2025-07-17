Here's the revised comprehensive EDA prompt with the nested folder requirement:

---

## **Comprehensive EDA Prompt for RNA Training Data Analysis**

**Objective**: Conduct a thorough Exploratory Data Analysis (EDA) on the RNA training dataset to understand sequence characteristics, coordinate quality, and structural features.

**Code Structure & Modularity**
- Never create a file longer than 500 lines of code. If a file approaches this limit, refactor by splitting it into modules or helper files.
- Organize code into clearly separated modules, grouped by feature or responsibility.

**Task Completion**
- Mark completed tasks in TASK.md immediately after finishing them.
- Add new sub-tasks or TODOs discovered during development to TASK.md under a “Discovered During Work” section.

**Style & Conventions**
- Use Python as the primary language.
- Follow PEP8, use type hints, and format with black.
- Use pydantic for data validation.
- Use FastAPI for APIs and SQLAlchemy or SQLModel for ORM if applicable.
- Write docstrings for every function using the Google style:

```
def example():
    """
    Brief summary.

    Args:
        param1 (type): Description.

    Returns:
        type: Description.
    """
```

**Documentation & Explainability**
- Update README.md when new features are added, dependencies change, or setup steps are modified.
- Comment non-obvious code and ensure everything is understandable to a mid-level developer.
- When writing complex logic, add an inline # Reason: comment explaining the why, not just the what.

**AI Behavior Rules**
- Never assume missing context. Ask questions if uncertain.
- Never hallucinate libraries or functions – only use known, verified Python packages.
- Always confirm file paths and module names exist before referencing them in code or tests.
- Never delete or overwrite existing code unless explicitly instructed to or if part of a task from TASK.md.


### **Step-by-Step Instructions:**

#### **Phase 1: Sequence Analysis (FASTA files in `competition/train/seqs/`)**

1. **Sequence Length Distribution Analysis**
   - Create a script that reads all `.fasta` files from `competition/train/seqs/`
   - Extract sequence lengths and generate comprehensive statistics (min, max, mean, median, std)
   - Create visualizations: histogram with appropriate binning, box plot, and cumulative distribution
   - Save results to `pdb_prune/EDA/sequence_length_analysis/`

2. **Base Composition Analysis**
   - Analyze the frequency of each RNA base (A, U, G, C) across all sequences
   - Calculate base composition percentages per sequence and overall dataset
   - Create visualizations: bar charts for overall base frequencies, heatmap for per-sequence composition
   - Identify any unusual base patterns or biases
   - Save results to `pdb_prune/EDA/base_composition_analysis/`

#### **Phase 2: Coordinate Analysis (NPY files in `competition/train/coords/`)**

3. **NaN Analysis in Coordinate Files**
   - Read all `.npy` files from `competition/train/coords/`
   - Analyze NaN patterns: count NaN values per file, identify positions with NaN values
   - **Important**: Exclude the first phosphorus atom (P) from NaN analysis as per your requirement
   - Create visualizations: histogram of NaN counts per file, heatmap showing NaN positions
   - Calculate statistics: percentage of files with NaN values, average NaN count per file
   - Save results to `pdb_prune/EDA/coordinate_nan_analysis/`

4. **Coordinate Structure Analysis**
   - Analyze the shape and structure of coordinate arrays
   - Calculate statistics on array dimensions (number of residues, atoms per residue)
   - Identify any structural patterns or anomalies
   - Save results to `pdb_prune/EDA/coordinate_structure_analysis/`

#### **Phase 3: PDB Structure Analysis (`competition/official_training_pdbs/`)**

5. **Side Chain Analysis**
   - Parse all `.ent` PDB files in `competition/official_training_pdbs/`
   - Count side chains per PDB structure
   - Analyze side chain distribution and identify crystal structure artifacts
   - Create visualizations: histogram of side chain counts, correlation with sequence length
   - Save results to `pdb_prune/EDA/pdb_sidechain_analysis/`

6. **Cross-Dataset Correlation Analysis**
   - Correlate sequence lengths with coordinate array dimensions
   - Analyze relationship between base composition and coordinate quality (NaN patterns)
   - Identify any systematic patterns across the three data types
   - Save results to `pdb_prune/EDA/cross_correlation_analysis/`

#### **Phase 4: Summary and Reporting**

7. **Generate Comprehensive EDA Report**
   - Create a Jupyter notebook that combines all analyses
   - Include interactive visualizations where appropriate
   - Generate summary statistics and key findings
   - Create a markdown report with insights and recommendations
   - Save all outputs to `pdb_prune/EDA/`

### **Technical Requirements:**

- Use existing utility functions from `utils/` where applicable
- Create new modular functions for each analysis component
- Implement proper error handling for missing or corrupted files
- Use pandas for data manipulation and matplotlib/seaborn for visualizations
- Ensure all outputs are saved with clear naming conventions
- Include progress bars for large file processing operations

### **Expected Outputs:**

1. **CSV files**: Detailed statistics and data tables
2. **PNG/PDF plots**: High-quality visualizations with publication level art standard
3. **JSON files**: Configuration and metadata
4. **Jupyter notebook**: Interactive analysis with explanations
5. **Markdown report**: Executive summary with key findings

### **File Structure for Results:**
```
pdb_prune/EDA/
├── sequence_length_analysis/
├── base_composition_analysis/
├── coordinate_nan_analysis/
├── coordinate_structure_analysis/
├── pdb_sidechain_analysis/
├── cross_correlation_analysis/
└── summary_report/
```

### **Quality Assurance:**
- Validate data integrity at each step
- Include data quality metrics and outlier detection
- Provide clear documentation for each analysis component
- Ensure reproducibility with seed setting for random operations

This prompt provides a comprehensive framework for conducting systematic EDA on your RNA training data, covering all the aspects you mentioned while maintaining clear structure and deliverables, with all results nested within the `pdb_prune/EDA/` folder structure.