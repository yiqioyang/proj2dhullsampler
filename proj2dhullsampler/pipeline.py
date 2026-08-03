"""Shared config-driven setup for the apply pipeline.

build_case() loads the data described by a config dict and returns a ready
HistoryMatching case (created + prepared, or loaded from disk if it already
exists). It is used both by run_apply.py (batch/CLI) and by apply.ipynb
(interactive), so the load/prepare logic only needs to be written once.
"""

import os

import numpy as np
import pandas as pd
import xarray as xr

from proj2dhullsampler import HistoryMatching


def default_worker_count():
    """CPUs actually available to this process.

    Prefers sched_getaffinity, which reflects the cgroup allocation a PBS/Slurm
    job gets (unlike os.cpu_count(), which reports the whole node's core count
    regardless of what the job was actually granted).
    """
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def build_case(config, mode="notebook", on_created=None):
    """Load data per `config` and return a created/prepared or loaded HistoryMatching case.

    `config` has the same shape as apply_config.json (data_paths, lat_bins,
    manual_regions, obs_dict, n_sample, prepare_case, threshold_level, ...).
    `mode` controls HistoryMatching's own notebook/python behavior (e.g.
    plt.show() vs. saving figures to disk) and defaults to "notebook" for
    interactive use; run_apply.py passes config["mode"] through instead.

    `on_created`, if given, is called with the new test_case right after
    create_case() makes the case directory but before prepare_case() runs.
    run_apply.py uses this hook to start teeing stdout to a log file once
    the directory it needs to write into actually exists; interactive
    notebook use has no need for it.
    """
    working_dir = config["working_dir"]
    case_name = config["case_name"]

    paths = config["data_paths"]
    obs_nc = xr.open_dataset(paths["obs_nc"]) if paths["obs_nc"] is not None else None
    obs_tab = (
        pd.read_csv(paths["obs_tab"], index_col=0).squeeze("columns")
        if paths["obs_tab"] is not None
        else None
    )
    para = pd.read_csv(paths["para"], index_col=0)
    ppe_nc = xr.open_dataset(paths["ppe_nc"]) if paths["ppe_nc"] is not None else None
    ppe_tab = (
        pd.read_csv(paths["ppe_tab"], index_col=0) if paths["ppe_tab"] is not None else None
    )

    if config["lat_bins"] is not None:
        lb = config["lat_bins"]
        lat_bins = np.arange(lb["start"], lb["stop"], lb["step"])
    else:
        lat_bins = None

    # Column order matters: local_process() joins row.astype(str) (nm, lat_min, lat_max,
    # lon_min, lon_max) to build diagnostic names like "local_PRECT_4_7_1_359", so the
    # manual_regions entries in the config must list keys in that exact order.
    if config["manual_regions"] is not None:
        manul_ppe_info = pd.DataFrame(
            config["manual_regions"],
            columns=["nm", "lat_min", "lat_max", "lon_min", "lon_max"],
        )
    else:
        manul_ppe_info = None

    obs_dict = config["obs_dict"] if config["obs_dict"] is not None else None

    # n_cpus falls back to the job's actual CPU allocation when left out of the
    # config, so a PBS/Slurm resource request doesn't also have to be
    # duplicated into the JSON by hand. Explicit config values always win.
    default_workers = default_worker_count()
    prepare_case_config = dict(config.get("prepare_case", {}))
    prepare_case_config.setdefault("n_cpus", default_workers)
    print(
        f"Worker count: prepare_case.n_cpus={prepare_case_config['n_cpus']} "
        f"(job allocation reports {default_workers} CPUs available; an explicit "
        f"config value, if set, takes precedence)"
    )

    test_case = HistoryMatching(working_dir, case_name, mode=mode)

    if not test_case.root.exists():
        print("Case directory does not exist yet; creating and preparing it.")
        test_case.create_case(
            para, [ppe_tab, obs_tab], ppe_nc, obs_nc, obs_dict,
            lat_bins, manul_ppe_info, n_sample=config["n_sample"], threshold_scale =config['threshold_level'] 
        )
        if on_created is not None:
            on_created(test_case)

        prepare_case_config["mode"] = mode
        test_case.prepare_case(prepare_case_config)
        test_case.load_mask(threshold_level=config["threshold_level"])

        for yname in config.get("visualize_check_vars", []):
            test_case.visualize_check(yname)
    else:
        print("Case directory already exists; loading it instead of re-preparing.")
        test_case.load_case()
        test_case.load_mask(threshold_level=config["threshold_level"])

    return test_case
