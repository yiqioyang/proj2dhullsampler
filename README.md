# proj2dhullsampler

`proj2dhullsampler` is a Python package for parameter-space screening and
history matching of spatial climate diagnostics. It includes utilities for:

- preparing observational and PPE-derived feature tables
- training and applying Gaussian process emulators
- building boolean masks of acceptable simulations
- grouping diagnostics by sensitive parameter pairs
- constructing alpha-shape hulls in normalized parameter space
- drawing new candidate parameter sets from the surviving region

## Repository Layout

```text
proj2dhullsampler/
├── proj2dhullsampler/
│   ├── prep_class.py
│   ├── hm_class.py
│   ├── sampling_functions.py
│   ├── preprocess.py
│   ├── plotting.py
│   ├── aux.py
│   ├── utils.py
│   └── unused_funs.py   # exploratory helpers kept for reference; not imported by the package
├── notebooks/
│   ├── prepare.ipynb
│   └── implementation.ipynb
├── tests/
├── pyproject.toml
└── README.md
```

## Core Workflow

The code is structured around three stages.

### 1. Prepare a case directory

`HistoryMatching` is the package entry point. Its `create_case` method creates
the case folder, writes uniformly sampled normalized parameters, and builds
tabular PPE and observation features. `prepare_case` then trains the emulators
and creates masks for the configured uncertainty thresholds.

Expected inputs are typically:

- `para`: parameter table as a `pandas.DataFrame`
- `tabs`: optional tuple of already-tabulated `(ppe_tab, obs_tab)` data
- `ppe`: PPE outputs as an `xarray.Dataset`
- `obs`: observations as an `xarray.Dataset`
- `obs_dict`: mapping from model variable names to observation variable names
- `lat_bins`: latitude bins for zonal aggregation
- `manul_ppe_info`: table describing manually selected regional averages

Example:

```python
from proj2dhullsampler import HistoryMatching

hm = HistoryMatching("/path/to/work", "case_a")
hm.create_case(
    para=parameter_table,
    tabs=None,
    ppe=ppe_ds,
    obs=obs_ds,
    obs_dict=obs_dict,
    lat_bins=lat_bins,
    manul_ppe_info=manual_regions,
    n_sample=1_000_000,
)

hm.prepare_case(
    {
        "n_cpus": 15,
        "threshold_levels": [2.0, 2.5],
    }
)
```

`prepare_case`'s config only controls `n_cpus` and `threshold_levels`. The
number of sensitive parameters retained per diagnostic (`n_sens_p`, default
`2`) is not exposed through this path; to change it, call
`hm.prep_case.sensitivity_emulation(n_sens_p=..., n_cpus=...)` directly
before calling `hm.load_case()`.

This creates a case directory like:

```text
case_a/
├── sampled_parameters.nc
├── meta.csv
├── tabs/
│   ├── parameters.csv
│   ├── ppe_data.csv
│   └── obs_data.csv
├── y_emu/
├── python_obj/
├── class_obj/
└── output/
```

### 2. Load masks and inspect emulator outputs

For an existing case, `HistoryMatching` loads saved masks, metadata, emulator
inputs, and feature tables.

Example:

```python
from proj2dhullsampler import HistoryMatching

hm = HistoryMatching(
    working_dir="/path/to/work",
    case_name="case_a",
)

hm.load_case()
hm.load_mask(threshold_level=2.0)
hm.visualize_check("PRECT_zonal_-30to30")
```

Typical files expected in the case directory at this stage include:

- `tf_masks_level_<threshold>.csv`
- `meta.csv`
- `sampled_parameters.nc`
- `tabs/parameters.csv`
- `tabs/ppe_data.csv`
- `tabs/obs_data.csv`
- `y_emu/gp_mean_std_<variable>.csv`

### 3. Run history matching and sample new parameters

`HistoryMatching` filters out uninformative diagnostics, groups the remaining
ones by their most sensitive parameter pair, resolves cases where the
surviving regions for co-grouped diagnostics don't overlap, builds alpha-shape
hulls, and draws new samples from the feasible region.

Example (mirrors the order used in `notebooks/apply.ipynb`):

```python
hm.drop_by_name(["local_PRECT_4_7_1_359"])   # drop diagnostics by name prefix
hm.drop_by_n_survive(n_survive=50)           # drop diagnostics that are always/rarely satisfied
hm.remove_var2d_auto(overlapping_threshold=10_000)  # resolve non-overlapping parameter-pair groups
hm.drop_by_nvar_per_pair(n_var_thre=1)       # optional: drop pairs backed by too few diagnostics

hm.prepare_for_sampling(
    shape_alpha=5,
    n_pts=10_000,
    n_threshold=100,
    sample_threshold=100_000,
    max_workers=8,
)
hm.draw(
    n_pts=50_000,
    n_threshold=5_000,
    sample_threshold=100_000_000,
    max_workers=8,
    n_max=1_000,
)
hm.save_samples_specifications(result_name="case_a", top_n=100)
hm.compare_with_original()  # optional: sanity-check sampled vs. original PPE parameter ranges
```

`remove_var2d_auto` calls `group_para_climatology` (grouping diagnostics by
sensitive parameter pair) and `shuffle_vars` (checking overlap within each
group) internally and iterates until no non-overlapping groups remain, or
until `no_iter` iterations are exhausted. Call `group_para_climatology`
directly only if you need the grouping without the overlap-resolution loop.

The final saved outputs are written under `case_a/output/`, including:

- `<result_name>_all_para_realscale.csv`
- `<result_name>_all_para_realscale.nc`
- `<result_name>_topn_para_realscale.csv`
- `<result_name>_topn_para_realscale.nc`
- `<result_name>_specifications.json`
- `<result_name>_dropped_vars.json`

## Public API

The package exposes one public class:

- `HistoryMatching`: create and prepare cases, filter diagnostics, build hulls,
  and generate candidate parameter samples.

The functions in the other package modules support `HistoryMatching` internally
and are not part of the stable top-level API.

## Installation

Python 3.10+ is required.

Install the package in editable mode:

```bash
pip install -e .
```

Install with development dependencies:

```bash
pip install -e .[dev]
```

## Development

Run the basic checks from the repository root:

```bash
ruff check .
black --check .
pytest
```

## Notes

- The documented notebooks are `notebooks/prepare.ipynb` and
  `notebooks/implementation.ipynb`.
- Most geometry operations assume parameters have been normalized to the
  `[0, 1]` range before hull construction and sampling.
