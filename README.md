# proj2dhullsampler

`proj2dhullsampler` is a Python package for parameter-space screening and
history matching of spatial climate diagnostics. It includes utilities for:

- preparing observational and PPE-derived feature tables
- training and applying Gaussian process emulators based on the most two sensitive parameters
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
│   ├── pipeline.py       # build_case(): shared config-driven setup used by
│   │                     # both run_apply.py and apply.ipynb
│   ├── run_apply.py      # notebook-free, config-driven script
│   ├── sampling_functions.py
│   ├── preprocess.py
│   ├── plotting.py
│   ├── aux.py
│   ├── utils.py
│   └── unused_funs.py   # exploratory helpers kept for reference; not imported by the package
├── application/
│   ├── apply.ipynb        # interactive, notebook-mode walkthrough
│   ├── apply_config.json  # example config for run_apply.py
│   ├── apply_config_annotated.txt # Annotated configuraiton file 
│   └── submit_apply.pbs   # PBS job template for Casper/Derecho
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
        "mode": "notebook",  # or "python"; see "Notebook vs. python mode" below
    }
)
```

`prepare_case`'s config controls `n_cpus`, `threshold_levels`, and (optionally)
`mode`. The number of sensitive parameters retained per diagnostic (`n_sens_p`,
default `2`) is not exposed through this path; to change it, call
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

Example (mirrors the order used in `application/apply.ipynb`):

```python
hm.drop_by_name(["local_PRECT_4_7_1_359"])   # drop diagnostics by name prefix
hm.drop_by_n_survive(n_survive=50)           # drop diagnostics that are always/rarely satisfied
hm.remove_var2d_auto(overlapping_threshold=10_000, added_num=100)  # resolve non-overlapping parameter-pair groups
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
sensitive parameter pair) and `shuffle_vars` (checking pairwise overlap within
each group) internally and iterates until no non-overlapping groups remain.
On each iteration it drops the diagnostic most implicated among variable
pairs whose overlap falls below a pairwise threshold (which starts equal to
`overlapping_threshold`). If a flagged group's *combined* intersection is
below `overlapping_threshold` but every individual pair within it is not
(a 3+-way interaction effect the pairwise check alone can't see), no pair
qualifies that iteration; the pairwise threshold is then relaxed by
`added_num` and the same group is re-evaluated on the next iteration, until
something qualifies. If `no_iter` iterations are exhausted with groups still
unresolved, it raises `ValueError` rather than returning silently. Call
`group_para_climatology` directly only if you need the grouping without the
overlap-resolution loop.

The final saved outputs are written under `case_a/output/`, including:

- `<result_name>_all_para_realscale.csv`
- `<result_name>_all_para_realscale.nc`
- `<result_name>_topn_para_realscale.csv`
- `<result_name>_topn_para_realscale.nc`
- `<result_name>_specifications.json`
- `<result_name>_dropped_vars.json`

## Notebook vs. python mode

`HistoryMatching` accepts a `mode` of `"notebook"` (the default) or `"python"`,
set either at construction (`HistoryMatching(working_dir, case_name, mode=...)`)
or via the `"mode"` key in the dict passed to `prepare_case` (which overrides
whatever was set at construction). It affects the two methods that produce
checkup figures during the workflow, `visualize_check` and
`compare_with_original`:

- `"notebook"`: figures are shown inline with `plt.show()`, as before.
- `"python"`: figures are saved as PNGs to `<case_dir>/diagnostics/` instead of
  being displayed, and the figure `Axes` are closed afterward to avoid leaking
  memory across a long-running script.

`mode` only controls figure display/saving. It doesn't affect any of the
`print()`-based status messages elsewhere in the pipeline (e.g. in
`remove_var2d_auto` or the sensitivity-emulation step); those still go to
stdout regardless of mode.

## Running without a notebook

`proj2dhullsampler/run_apply.py` replays the same workflow as
`application/apply.ipynb` as a plain script, driven entirely by a JSON config
file (see `application/apply_config.json` for a fully worked example,
including data paths, `obs_dict`, `lat_bins`, manually selected regions,
thresholds, and `mode`; see `application/apply_config.annotated.jsonc` for a
field-by-field commented walkthrough of what each key controls — it's a
documentation-only copy, not something `run_apply.py` reads):

```bash
python proj2dhullsampler/run_apply.py --config application/apply_config.json
```

The data-loading and create-or-load-case logic (everything through
`load_mask`) lives in `proj2dhullsampler/pipeline.py`'s `build_case(config,
mode)`, which both `run_apply.py` and `apply.ipynb` call, so the two stay in
sync: a config change only needs to be made once. `run_apply.py` itself
handles only the CLI-specific parts (argument parsing, the batch diagnostics
log, and the drop/sample/save steps that follow `build_case`).

If the case directory doesn't exist yet, `build_case` creates and prepares it
(equivalent to `create_case` + `prepare_case`); if it already exists, it
loads it instead of re-running the (expensive) emulator training. When
`config["mode"]` is `"python"`, `run_apply.py` also:

- switches matplotlib to the non-interactive `Agg` backend, since there's no
  display to show figures on,
- tees all of the pipeline's `print()` output to
  `<case_dir>/diagnostics/run_log.txt` in addition to stdout, so the
  intermediate diagnostic messages survive an unattended/batch run.

### Worker counts

`prepare_case.n_cpus` (GP training parallelism) and `max_workers` (hull
sampling parallelism) in the config are optional. If either is left out,
`run_apply.py` falls back to the number of CPUs actually available to the
process (`len(os.sched_getaffinity(0))`, which reflects a PBS/Slurm job's
cgroup allocation, not just the whole node's core count). An explicit value
in the config always takes precedence. The resolved values are printed (and,
in `"python"` mode, logged to `run_log.txt`) at the start of each run.

### Submitting as a PBS job

`application/submit_apply.pbs` is a template job script for Casper/Derecho:

```bash
qsub application/submit_apply.pbs
```

Fill in `<PROJECT_CODE>` and adjust the `#PBS -q`/`select`/`walltime`
directives for your case size. If you also drop `prepare_case.n_cpus` and
`max_workers` from `apply_config.json`, the script automatically uses
whatever CPU count the job was actually granted (see "Worker counts" above)
instead of requiring the PBS resource request and the JSON config to be kept
in sync by hand.

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

- The application entry points are `application/apply.ipynb` (interactive)
  and `proj2dhullsampler/run_apply.py` (batch/PBS, driven by
  `application/apply_config.json`); both cover the same workflow described
  above via the shared `build_case()` helper.
- Most geometry operations assume parameters have been normalized to the
  `[0, 1]` range before hull construction and sampling.
- `apply_config.json`'s top-level `n_pts`/`n_threshold`/`sample_threshold` are
  passed to both `prepare_for_sampling()` (hull construction/orchestration)
  and `draw()` (final sample generation). `shape_alpha` is not exposed in the
  config and stays at `HistoryMatching.prepare_for_sampling`'s default (`5`)
  for both stages.
- The emulator "validation" error ratio (`validation_error_ratio.csv`, used
  by `drop_by_emulator_performance`/`emultor_error_ratio_threshold`) is
  currently computed by evaluating each GP on the same data it was trained
  on, not a held-out split — treat it as an in-sample fit-quality check
  rather than a generalization estimate.
