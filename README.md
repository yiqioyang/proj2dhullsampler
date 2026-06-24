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
│   └── utils.py
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

`Prepare_Case` creates a case folder, writes uniformly sampled normalized
parameters, and builds tabular PPE and observation features.

`Prep_Mask_Generation` is kept as a public import alias for older notebooks and
scripts. It points to the same implementation as `Prepare_Case`.

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
from proj2dhullsampler import Prepare_Case

prep = Prepare_Case(
    working_dir="/path/to/work",
    case_name="case_a",
    para=parameter_table,
    tabs=None,
    ppe=ppe_ds,
    obs=obs_ds,
    obs_dict=obs_dict,
    lat_bins=lat_bins,
    manul_ppe_info=manual_regions,
    n_sample=1_000_000,
)

prep.sensitivity_emulation(n_sens_p=2, n_cpus=15)
prep.mask_generation(threshold_level=2.0)
```

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

`HistoryMatching` loads saved masks, metadata, emulator inputs, and feature
tables. `Analysis` is kept as a public import alias to `HistoryMatching` for
older notebooks and scripts.

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

`HistoryMatching` groups diagnostics by parameter pair, builds alpha-shape
hulls, checks overlap, and draws new samples from the feasible region.

Example:

```python
hm.drop_by_n_survive(n_survive=50)
hm.group_para_climatology(overlapping_threshold=10_000)
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
```

The final saved outputs are written under `case_a/output/`, including:

- `<result_name>_all_para_realscale.csv`
- `<result_name>_all_para_realscale.nc`
- `<result_name>_topn_para_realscale.csv`
- `<result_name>_topn_para_realscale.nc`
- `<result_name>_specifications.json`
- `<result_name>_dropped_vars.json`

## Public API

Main classes:

- `Prepare_Case`: prepare case directories, sampled parameters, emulators, and masks
- `HistoryMatching`: filter diagnostics, build hulls, and generate candidate samples
- `Prep_Mask_Generation`: compatibility alias for `Prepare_Case`
- `Analysis`: compatibility alias for `HistoryMatching`

Sampling and plotting helpers:

- `sample_from_hull`
- `sample_from_hulls_n`
- `orchestrate_test`
- `biplot`
- `biplot_original_scale`

Utility helpers:

- `feature_builder`: build PPE and observation feature tables
- `gp_training_application`: fit GP emulators and write mean/std predictions
- `metric_cal_single`: compute weighted average/bias/RMSE style metrics
- `para_csv2nc`: convert sampled parameter CSV files to NetCDF

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
