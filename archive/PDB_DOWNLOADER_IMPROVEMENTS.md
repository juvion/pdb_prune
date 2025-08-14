# PDB Downloader Improvements

This document describes the improvements made to the PDB downloader to fix connection issues and add resume functionality.

## Issues Fixed

1. **Connection Breaking**: The original downloader would stop when encountering network errors or timeouts
2. **No Resume Capability**: If the download was interrupted, you had to start from the beginning
3. **Limited Error Handling**: Poor handling of API rate limits and network issues
4. **No Progress Tracking**: No way to see what had been found vs. what had been downloaded

## New Features

### 1. Enhanced Error Handling
- **Increased Timeout**: API requests now have a 45-second timeout (up from 30)
- **Exponential Backoff**: Retry delays are capped at 30 seconds to prevent excessive waiting
- **Better Logging**: More detailed logging shows page numbers and progress
- **Graceful Interruption**: Ctrl+C saves progress before exiting

### 2. Progress Saving
- **found_pdb_ids.txt**: Automatically saves all discovered PDB IDs
- **Periodic Saves**: Progress is saved every 5 pages during search
- **Final Save**: Complete list is saved at the end of search phase

### 3. Resume Functionality
- **Search Resume**: New `--search_start_page` parameter to resume API search from specific page
- **Download Resume**: Existing `--start_from` parameter to resume downloads from specific index
- **Auto-detection**: Helper script can automatically detect where to resume

### 4. Safety Limits
- **Max Search Attempts**: Limited to 50 search attempts to prevent infinite loops
- **Increased Delays**: 2-second delays between API calls to be more respectful

## Usage

### Basic Usage (Same as Before)
```bash
python utils/pdb_downloader.py --max_entries 1000 --output_dir my_pdbs
```

### Resume from Specific Search Page
```bash
python utils/pdb_downloader.py --max_entries 7000 --output_dir downloaded_rna_pdbs2 --search_start_page 20
```

### Resume Downloads from Specific Index
```bash
python utils/pdb_downloader.py --max_entries 7000 --output_dir downloaded_rna_pdbs2 --start_from 500
```

### Resume Both Search and Downloads
```bash
python utils/pdb_downloader.py --max_entries 7000 --output_dir downloaded_rna_pdbs2 --search_start_page 20 --start_from 500
```

## Helper Script: resume_pdb_search.py

A new helper script makes resuming easier by automatically analyzing progress.

### Analyze Current Progress
```bash
python resume_pdb_search.py --output_dir downloaded_rna_pdbs2 --analyze_only
```

Output:
```
=== Progress Analysis ===
Found PDB IDs: 1900
Downloaded files: 150
Failed downloads: 5
Suggested search start page: 19
Suggested download start index: 150
========================
```

### Auto-Resume
```bash
python resume_pdb_search.py --output_dir downloaded_rna_pdbs2 --max_entries 7000
```

This automatically detects where to resume and continues the download.

## File Structure

The improved downloader creates these files in the output directory:

```
downloaded_rna_pdbs2/
├── found_pdb_ids.txt          # List of all discovered PDB IDs
├── failed_downloads.txt       # List of failed download attempts
├── pdb1a34.ent               # Downloaded PDB files
├── pdb1a3m.ent
└── ...
```

## Recovery Scenarios

### Scenario 1: Search Phase Interrupted
If the search phase is interrupted (e.g., at page 20 with 1900 IDs found):

```bash
# Resume search from page 20
python utils/pdb_downloader.py --output_dir downloaded_rna_pdbs2 --search_start_page 20 --max_entries 7000
```

### Scenario 2: Download Phase Interrupted
If downloads are interrupted (e.g., 150 files downloaded out of 1900 found):

```bash
# Resume downloads from index 150
python utils/pdb_downloader.py --output_dir downloaded_rna_pdbs2 --start_from 150 --max_entries 7000
```

### Scenario 3: Complete Interruption
If both search and download are interrupted:

```bash
# Let the helper script figure it out
python resume_pdb_search.py --output_dir downloaded_rna_pdbs2 --max_entries 7000
```

## Best Practices

1. **Use the Helper Script**: `resume_pdb_search.py` is the easiest way to resume
2. **Check Progress Regularly**: Use `--analyze_only` to check progress without starting downloads
3. **Backup Progress Files**: Keep copies of `found_pdb_ids.txt` for safety
4. **Monitor Logs**: Watch for repeated failures that might indicate API issues
5. **Be Patient**: The improved delays make the process more reliable but slower

## Troubleshooting

### "Failed to fetch page after X attempts"
- **Cause**: Network issues or API rate limiting
- **Solution**: Wait a few minutes and resume from the failed page

### "Search interrupted by user"
- **Cause**: Ctrl+C pressed during search
- **Solution**: Progress is automatically saved, just resume with suggested parameters

### "Only found X unique PDB entries, less than requested"
- **Cause**: Reached the end of available RNA PDB entries
- **Solution**: This is normal, proceed with available entries

## Technical Details

### Search Algorithm
- Searches RCSB PDB API for RNA-containing entries
- Processes 100 entries per page
- Maintains a set of unique PDB IDs to avoid duplicates
- Saves progress every 5 pages and at completion

### Download Algorithm
- Downloads in configurable batches (default: 50)
- Checks for existing files to avoid re-downloading
- Uses Bio.PDB.PDBList for reliable downloads
- Implements retry logic with exponential backoff

### Resume Logic
- Search resume: Loads existing PDB IDs and continues from specified page
- Download resume: Skips already downloaded files and continues from specified index
- Auto-detection: Analyzes existing files to suggest optimal resume parameters