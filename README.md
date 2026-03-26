# proj2dhullsampler

`proj2dhullsampler` is a Python repository for parameter-space screening and history matching of spatial climate diagnostics. It packages utilities for:

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
│   ├── analysis.py
│   ├── hm_class.py
│   ├── sampling_functions.py
│   ├── plotting.py
│   ├── aux.py
│   └── utils.py
├── notebooks/
├── tests/
├── pyproject.toml
└── README.md
```

## Core Workflow

The code is structured around three stages.

### 1. Prepare a case directory

`Prep_Mask_Generation` creates a case folder, writes uniformly sampled normalized parameters, and builds tabular PPE and observation features.

Expected inputs are typically:

- `ppe`: PPE outputs as an `xarray.Dataset`
- `obs`: observations as an `xarray.Dataset`
- `obs_dict`: mapping from model variable names to observation variable names
- `para`: parameter table as a `pandas.DataFrame`
- `lat_bins`: latitude bins for zonal aggregation
- `manul_ppe_info`: table describing manually selected regional averages
- `added_ppe_obs`: optional additional PPE/obs tabular features

Example:

```python
from proj2dhullsampler import Prep_Mask_Generation

prep = Prep_Mask_Generation(
    working_dir="/path/to/work",
    case_name="case_a",
    ppe=ppe_ds,
    obs=obs_ds,
    obs_dict=obs_dict,
    para=parameter_table,
    lat_bins=lat_bins,
    manul_ppe_info=manual_regions,
    added_ppe_obs=None,
    n_sample=1_000_000,
)
```

This creates a case directory like:

```text
case_a/
├── sampled_parameters.nc
├── tabs/
│   ├── ppe_tab.csv
│   └── obs_tab.csv
├── y_emu/
├── python_obj/
└── class_obj/
```

### 2. Inspect emulator outputs

`Analysis` loads saved masks, metadata, emulator inputs, and feature tables for diagnostics and plotting.

Example:

```python
from proj2dhullsampler import Analysis

analysis = Analysis(
    working_dir="/path/to/work",
    case_name="case_a",
    ppe_para=parameter_table,
    threshold_level=2.0,
)

analysis.plot_onehot()
fig, axes = analysis.plot_by_para("clubb_c1")
analysis.visualize_check("PRECT_zonal_-30to30")
```

Typical files expected in the case directory at this stage include:

- `tf_masks_level_<threshold>.csv`
- `meta.csv`
- `sampled_parameters.nc`
- `tabs/ppe_tab.csv`
- `tabs/obs_tab.csv`
- `y_emu/gp_mean_std_<variable>.csv`

### 3. Run history matching and sample new parameters

`HistoryMatching` groups diagnostics by parameter pair, builds alpha-shape hulls, checks overlap, and draws new samples from the feasible region.

Example:

```python
from proj2dhullsampler import HistoryMatching

hm = HistoryMatching(
    working_dir="/path/to/work",
    case_name="case_a",
    ppe_para=parameter_table,
    threshold_level=2.0,
)

hm.drop_by_n_survive(n_survive=50)
hm.update_meta()
hm.group_para_climatology(overlapping_threshold=10_000)
hm.build_hulls(shape_alpha=5)
hm.orchestrate(n_pts=10_000, n_threshold=100, max_workers=8)
hm.draw(n_pts=50_000, n_threshold=5_000, max_workers=8, n_max=1000)
hm.save_samples(n=100)
hm.write_specifications()
```

The final saved outputs are typically:

- `full_sel_para_realscale.csv`
- `full_sel_para_realscale.nc`
- `sel_para_realscale.csv`
- `sel_para_realscale.nc`
- a JSON specification file written by `write_specifications()`

## Public API

Main classes:

- `Prep_Mask_Generation`: prepare case directories, sampled parameters, and feature tables
- `Analysis`: inspect masks, metadata, and emulator behavior
- `HistoryMatching`: filter diagnostics, build hulls, and generate candidate samples

Sampling and plotting helpers:

- `sample_from_hull`
- `sample_from_hulls_n`
- `orchestrate_test`
- `biplot`
- `biplot_original_scale`
- `plot_histograms_grid_5`

Utility helpers:

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

- The notebooks in [`notebooks/`](/glade/u/home/qingyuany/repos/spatialtuning/notebooks) appear to document the intended workflow and are the best place to see project-specific usage patterns.
- Most geometry operations assume parameters have been normalized to the `[0, 1]` range before hull construction and sampling.
