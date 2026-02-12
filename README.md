# MorphGen-Rates Documentation

## Overview

`morphgen-rates` is a Python package that estimates bifurcation and annihilation rates from morphological data obtained from experimentally reconstructed neuron morphologies. These rates are core parameters in mathematical models used to generate realistic neuron morphologies via random-walk–based generation processes.

### What does this package do?

The package analyzes neuronal morphology data to estimate two spatially varying rates that govern dendritic tree generation:

- **Bifurcation rate (β)**: Probability per unit distance that a synthesized dendritic segment splits into two branches
- **Annihilation rate (α)**: Probability per unit distance that a synthesized dendritic segment terminates

The rates vary with distance from the soma and are inferred from:

- **Sholl analysis data**: Mean and standard deviation of dendrite intersection counts across radial distance bins from the soma
- **Bifurcation statistics (Optional)**: Summary statistics describing the number of branching points in the reconstructed dendritic tree

This package implements computational methods from:

- Cavarretta F (2025) *A mathematical model for data-driven synthesis of neuron morphologies based on random walks.* Front. Appl. Math. Stat. 11:1632271. doi: 10.3389/fams.2025.1632271

For questions, issues, or contributions, please contact:

- Author: Francesco Cavarretta
- Email: fcavarretta@ualr.edu

## Installation

### Requirements

- Python >= 3.9
- numpy >= 1.23
- pandas >= 1.5
- pyomo >= 6.6
- A nonlinear optimization solver (typically IPOPT)

### Install the package

```bash
pip install morphgen-rates
```

### Install IPOPT (recommended solver)

The package requires a nonlinear solver. IPOPT is recommended.

On Ubuntu/Debian:

```bash
sudo apt-get install coinor-libipopt-dev
pip install cyipopt
```

On macOS:

```bash
brew install ipopt
pip install cyipopt
```

On Windows:

Follow instructions at:

https://github.com/coin-or/Ipopt

## Data Format Specifications

### Built-in Dataset Format

The package includes a CSV file (`morph_data.csv`) with the following columns:

| Column | Type | Description |
|---|---|---|
| area | str | Brain region (e.g., “aPC”, “Neocortex”, “OB”) |
| neuron_type | str | Cell class (e.g., “PYR”, “SL”, “MITRAL”, “TUFTED”) |
| neuron_name | str | Individual cell identifier |
| section_type | str | Dendrite type ("apical_dendrite" or "basal_dendrite") |
| bifurcation_count | float | Number of bifurcation points |
| total_length | float | Total dendritic length (μm) |
| bin_size | float | Radial bin size for Sholl analysis (μm) |
| Count0 | float | Primary dendrite count (Sholl at r=0) |
| Count1 | float | Sholl intersections at bin 1 |
| Count2 | float | Sholl intersections at bin 2 |
| … | … | … |
| CountN | float | Sholl intersections at bin N |

## Troubleshooting

### Solver Not Found

Error:

```text
RuntimeError: Solver 'ipopt' is not available
```

Solution:

- Install IPOPT solver (see Installation section)

### Data Not Found

Error:

```text
AssertionError: The area ... or neuron class ... are not known
```

Solution:

- Check available combinations using:

```python
import pandas as pd

df = pd.read_csv("morph_data.csv")
print("Available combinations:")
print(df[["area", "neuron_type"]].drop_duplicates())
```

### Numerical Instability

Symptom:

- Very large or very small rate values

Solution:

- Check that max_step_size is appropriate for your data scale
- Ensure Sholl data is smoothly decreasing (no sudden jumps)
- Verify that variance values are non-negative

### Poor Moment Matching in Primary Dendrite Distribution

Symptom:

- Computed distribution doesn't match target mean/std

Solution:

- Increase slack_penalty:

```python
probs = compute_init_number_probs(
    mean_primary_dendrites=3.5,
    sd_primary_dendrites=1.2,
    min_primary_dendrites=1,
    max_primary_dendrites=7,
    slack_penalty=1.0  # Increase from default 0.1
)
```
