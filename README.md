# PBLMM

## Peptide-Based Linear Mixed Model for Proteomics Data Analysis

**Version 2.1.1** (2023-10-23)

PBLMM is a Python package that implements a peptide-based linear mixed model approach for analyzing proteomics data. It allows for robust statistical analysis of protein abundance changes across different conditions by properly accounting for peptide-level variation.

## Features

- **Peptide-based statistical modeling** using linear mixed models (LMM)
- **Multiple protein rollup methods** (sum, median, mean)
- **Hypothesis testing** with both LMM and t-test approaches
- **Multiple testing correction** using FDR control
- **Flexible data processing** compatible with various proteomics data formats

## Installation

```bash
pip install git+https://github.com/science64/PBLMM.git
```

## Quick Start

```python
import pandas as pd
from PBLMM import HypothesisTesting, Defaults

# Initialize defaults
defaults = Defaults()

# Load your proteomics data (peptide or PSM level)
data = pd.read_csv("your_proteomics_data.csv")

# Define your experimental conditions
conditions = ["Control", "Treatment", "Control", "Treatment"]

# Define pairs for comparison
pairs = [["Treatment", "Control"]]

# Initialize hypothesis testing
ht = HypothesisTesting(defaults)

# Run peptide-based LMM analysis
results = ht.peptide_based_lmm(
    input_file=data,
    conditions=conditions,
    pairs=pairs
)

# Save results
results.to_csv("pblmm_results.csv")
```

## Understanding the Statistical Output

### P-values

The PBLMM code generates p-values using two different statistical approaches:

#### 1. Linear Mixed Models (LMM)

In the `peptide_based_lmm` method:
- P-values are generated by fitting a linear mixed model to peptide-level data
- The model uses peptide sequences as random effects (`groups='Sequence'`)
- The model formula used is: `"value ~ variable"` where:
  - `value` represents log2-transformed abundance measurements
  - `variable` represents the experimental conditions being compared
- P-values test the significance of the condition effect (null hypothesis: no difference between conditions)

#### 2. Unpaired t-test

In the `ttest` method:
- P-values come from a standard two-sample t-test with equal variance assumption
- For each protein, the test compares abundance values between two specified conditions
- Tests whether the means of the two conditions are significantly different

### Q-values

Q-values are only calculated in the `peptide_based_lmm` method:
- They represent p-values adjusted for multiple hypothesis testing
- Correction is performed using the Benjamini-Hochberg procedure (FDR correction)
- Q-values control the false discovery rate (FDR), which is the expected proportion of false positives among all significant results
- This adjustment is critical when testing many proteins simultaneously to minimize false discoveries

## Advanced Usage

### Custom Column Names

You can customize the column names used for protein accessions, sequences, and abundance values:

```python
defaults = Defaults()
defaults.MasterProteinAccession = "Your Protein ID Column"
defaults.sequence = "Your Sequence Column"
defaults.AbundanceColumn = "Your Abundance Column Prefix"
```

### Technical Replicates and Multiplexing

For experiments with technical replicates or multiple plexes:

```python
results = ht.peptide_based_lmm(
    input_file=data,
    conditions=conditions,
    techreps=["rep1", "rep1", "rep2", "rep2"],
    plexes=["plex1", "plex1", "plex2", "plex2"],
    pairs=pairs
)
```

## Citation

If you use PBLMM in your research, please cite:

```
Original implementation by Kevin Klann
Updated by Süleyman Bozkurt
```

## License

MIT

## Contributors

- Kevin Klann (Original author)
- Süleyman Bozkurt (Maintainer)
