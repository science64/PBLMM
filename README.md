# PBLMM

## Peptide-Based Linear Mixed Model for TMT-based Proteomics Data Analysis

**Version 2.1.3** (2026-02-02)

PBLMM is a Python package that implements a peptide-based linear mixed model approach for analyzing proteomics data. It allows for robust statistical analysis of protein abundance changes across different conditions by properly accounting for peptide-level variation.

## Features

- **Peptide-based statistical modeling** using linear mixed models (LMM)
- **Multiple protein rollup methods** (sum, median, mean)
- **Hypothesis testing** with both LMM and t-test approaches
- **Multiple testing correction** using FDR control
- **Flexible data processing** compatible with various proteomics data formats
- No Support for pandas version 3 (recommended version 2.3.3).

## Compatibility

| Software | Supported versions | Notes |
|---|---|---|
| Proteome Discoverer (TMT multiplex) | up to 3.2 | Validated as of 23.03.2026 |

## Installation

```bash
pip install git+https://github.com/science64/PBLMM.git
pip install pandas==2.3.3
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

## Input Data

Both `peptide_based_lmm` and `ttest` accept the **same input**: a **peptide-level** (or PSM-level) DataFrame where each row is a single peptide/PSM measurement. The file must contain:

| Column | Default name | Description |
|---|---|---|
| Peptide sequence | `Annotated Sequence` | Unique identifier for each peptide |
| Protein accession | `Master Protein Accessions` | Protein the peptide maps to |
| Abundance columns | columns containing `Abundance:` | Raw quantitative values per sample/channel |

> You can override these defaults via the `Defaults` class (see [Custom Column Names](#custom-column-names)).

### How each method uses the peptide input

| | `peptide_based_lmm` | `ttest` |
|---|---|---|
| **Statistics computed at** | Peptide level | Protein level |
| **Peptide → Protein conversion** | Done separately for the abundance output table only; **statistics use each peptide row individually** | Done **before** statistics: peptides are summed per protein (`protein_rollup_sum`), and the t-test is run on the resulting protein-level values |
| **Why peptide level matters for LMM** | Each peptide observation is a data point in the mixed model; peptide sequence is modeled as a random effect, so biological variance at the peptide level directly informs the p-value | Not applicable — peptide identity is discarded after rollup |

In short: **`peptide_based_lmm` keeps every peptide as a separate observation and lets the model account for peptide-to-peptide variability**, while **`ttest` first collapses peptides into a single protein score by summing, then tests those protein scores**.

## Understanding the Statistical Output

### P-values

The PBLMM code generates p-values using two different statistical approaches:

#### 1. Linear Mixed Models (LMM) — `peptide_based_lmm`

- The peptide file is melted into long format; each row is one peptide measurement in one sample
- A linear mixed model is fit per protein: `value ~ variable` with peptide sequence as a random effect (`groups='Sequence'`)
  - `value` — log2-transformed abundance of each individual peptide
  - `variable` — experimental condition label for that sample
- The p-value tests the fixed effect of condition (null hypothesis: no difference between conditions)
- Because the model operates directly on peptide-level data, it properly accounts for the fact that multiple peptides from the same protein are not independent

#### 2. Unpaired t-test — `ttest`

- Peptides are first rolled up to protein level by **summing** all peptide abundances per protein per sample (`protein_rollup_sum`)
- A standard two-sample t-test (equal variance) is then run on these protein-level summed values
- The log2 fold change is calculated from the mean of the summed values across replicates
- Statistics are therefore at the **protein level**, not the peptide level

### Q-values

Q-values are calculated in both methods using the same approach:
- Correction is performed using the Benjamini-Hochberg (BH) procedure (`fdr_bh`) from `statsmodels.stats.multitest.multipletests`
- Q-values control the false discovery rate (FDR): the expected proportion of false positives among all declared significant results
- All per-protein p-values for a given pair are collected first, then FDR correction is applied across all proteins at once
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
