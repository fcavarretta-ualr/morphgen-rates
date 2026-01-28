````markdown
# morphgen-rates

Python package to compute **bifurcation** and **annihilation** rates from morphology-derived data.

## Install (local development)

From the repository root:

```bash
pip install -e .
````

## Quick start

```python
from morphgen_rates import compute_rates

# Replace `data` with your input object (events, traces, or morphology structures)
rates = compute_rates(data)

print(rates["bifurcation_rate"])
print(rates["annihilation_rate"])
```

## API

### `compute_rates(data, *, dt=None, **kwargs)`

Computes bifurcation and annihilation rates from the provided input.

* `data` is your input dataset (format depends on your project)
* `dt` can be used to convert counts-per-step into a physical/time rate (optional)

Returns a dictionary with (at minimum):

* `bifurcation_rate`
* `annihilation_rate`

## Development

Run tests:

```bash
pytest
```

## License

GPL-3.0-or-later

