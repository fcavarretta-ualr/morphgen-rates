# morphgen-rates

Python package to compute **bifurcation** and **annihilation** rates from morphology-derived data.

## Provenance

This package is based on code originally provided in the repository:

- https://github.com/FrancescoCavarretta/NeuronSynthesis

The original repository was developed for the following publication:

- Cavarretta F (2025) *A mathematical model for data-driven synthesis of neuron morphologies based on random walks.* Front. Appl. Math. Stat. 11:1632271. doi: 10.3389/fams.2025.1632271


## Quick start

```python
from morphgen_rates import compute_rates

# Replace `data` with your input object (events, traces, or morphology structures)
rates = compute_rates(data)

print(rates["bifurcation_rate"])
print(rates["annihilation_rate"])
```


## License

GPL-3.0-or-later

