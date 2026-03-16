# spatialtuning

Utilities for spatial tuning workflows, including data preparation, analysis, history matching, plotting, and sampling.

## Package Structure

Main package:
- `funs/analysis.py`: analysis helpers and plotting wrappers
- `funs/prep_class.py`: preparation and mask generation workflows
- `funs/hm_class.py`: history matching workflow
- `funs/sampling_functions.py`: hull-based sampling helpers
- `funs/plotting.py`: plotting utilities
- `funs/aux.py`: I/O and metric helpers
- `funs/utils.py`: shared numeric and geometry utilities

Public API is exposed in `funs/__init__.py`.

## Installation

From repository root:

```bash
pip install -e .
```

Install development tools:

```bash
pip install -e .[dev]
```

## Quick Start

```python
from funs import Analysis, HistoryMatching, Prep_Mask_Generation

# Example: initialize analysis
# analysis = Analysis(working_dir, case_name, ppe_para)
```

## Notebooks

Project notebooks are under the repository root and `notebooks/`.
Prefer importing from `funs` rather than using `%run` on module files.

## Development

Formatting and linting configuration is in `pyproject.toml`.

Common checks:

```bash
ruff check .
black --check .
pytest
```
