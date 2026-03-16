# spatialtuning

Utilities for spatial tuning workflows, including data preparation, analysis, history matching, plotting, and sampling.

## Package Structure

Main package:
- `proj2dhullsampler/analysis.py`: analysis helpers and plotting wrappers
- `proj2dhullsampler/prep_class.py`: preparation and mask generation workflows
- `proj2dhullsampler/hm_class.py`: history matching workflow
- `proj2dhullsampler/sampling_functions.py`: hull-based sampling helpers
- `proj2dhullsampler/plotting.py`: plotting utilities
- `proj2dhullsampler/aux.py`: I/O and metric helpers
- `proj2dhullsampler/utils.py`: shared numeric and geometry utilities

Public API is exposed in `proj2dhullsampler/__init__.py`.

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
from proj2dhullsampler import Analysis, HistoryMatching, Prep_Mask_Generation

# Example: initialize analysis
# analysis = Analysis(working_dir, case_name, ppe_para)
```

## Notebooks

Project notebooks are under the repository root and `notebooks/`.
Prefer importing from `proj2dhullsampler` rather than using `%run` on module files.

## Development

Formatting and linting configuration is in `pyproject.toml`.

Common checks:

```bash
ruff check .
black --check .
pytest
```
