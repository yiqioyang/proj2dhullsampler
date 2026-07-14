"""Config-driven, notebook-free replay of apply.ipynb.

Usage:
    python run_apply.py --config apply_config.json

All paths, thresholds, and the "notebook" vs "python" mode flag live in the
config file instead of being hardcoded, so a run can be reproduced or tuned
without touching this script. See apply_config.json for the expected shape.

In "python" mode, checkup figures that would normally display inline in the
notebook are instead written as PNGs to <working_dir>/<case_name>/diagnostics/,
and all print() output from the pipeline is additionally teed to
<working_dir>/<case_name>/diagnostics/run_log.txt.
"""

import argparse
import json
import os
import sys
from pathlib import Path


def load_config(config_path):
    with open(config_path) as f:
        return json.load(f)


def _default_worker_count():
    """CPUs actually available to this process.

    Prefers sched_getaffinity, which reflects the cgroup allocation a PBS/Slurm
    job gets (unlike os.cpu_count(), which reports the whole node's core count
    regardless of what the job was actually granted).
    """
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


class _Tee:
    """Duplicates writes to multiple streams (used to mirror stdout to a log file)."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)

    def flush(self):
        for s in self.streams:
            s.flush()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Path to a JSON config file")
    args = parser.parse_args()

    config = load_config(args.config)
    mode = config.get("mode", "notebook")
    if mode not in ("notebook", "python"):
        raise ValueError(f"config['mode'] must be 'notebook' or 'python', got {mode!r}")

    working_dir = Path(config["working_dir"])
    case_name = config["case_name"]

    if mode == "python":
        # Non-interactive backend: nothing will be on screen to plt.show() into.
        import matplotlib
        matplotlib.use("Agg")

    _run(config, mode, working_dir, case_name)


def _start_diagnostics_log(case_root, log_file_box):
    """Create <case_root>/diagnostics/run_log.txt and tee stdout into it.

    Must only be called once the case directory actually exists on disk
    (either because it already did, or because create_case() just made it) -
    creating it any earlier would make HistoryMatching.create_case()'s own
    "does this case already exist" check see a directory that we ourselves made.

    The opened handle is appended to log_file_box (rather than returned) so that
    _run()'s finally block can find and close it even if the pipeline raises
    partway through, after the log has already been started.
    """
    diagnostics_dir = case_root / "diagnostics"
    diagnostics_dir.mkdir(parents=True, exist_ok=True)
    log_file = open(diagnostics_dir / "run_log.txt", "w")
    sys.stdout = _Tee(sys.stdout, log_file)
    log_file_box.append(log_file)


def _run(config, mode, working_dir, case_name):
    log_file_box = []
    try:
        _run_pipeline(config, mode, working_dir, case_name, log_file_box)
    finally:
        if log_file_box:
            sys.stdout = sys.stdout.streams[0]
            log_file_box[0].close()


def _run_pipeline(config, mode, working_dir, case_name, log_file_box):
    import numpy as np
    import pandas as pd
    import xarray as xr

    from proj2dhullsampler import HistoryMatching

    case_root = working_dir / case_name
    case_existed = case_root.exists()

    if mode == "python" and case_existed:
        _start_diagnostics_log(case_root, log_file_box)

    paths = config["data_paths"]
    obs_nc = xr.open_dataset(paths["obs_nc"])
    obs_tab = pd.read_csv(paths["obs_tab"], index_col=0).squeeze("columns")
    para = pd.read_csv(paths["para"], index_col=0)
    ppe_nc = xr.open_dataset(paths["ppe_nc"])
    ppe_tab = pd.read_csv(paths["ppe_tab"], index_col=0)

    lb = config["lat_bins"]
    lat_bins = np.arange(lb["start"], lb["stop"], lb["step"])

    # Column order matters: local_process() joins row.astype(str) (nm, lat_min, lat_max,
    # lon_min, lon_max) to build diagnostic names like "local_PRECT_4_7_1_359", so the
    # manual_regions entries in the config must list keys in that exact order.
    manul_ppe_info = pd.DataFrame(
        config["manual_regions"],
        columns=["nm", "lat_min", "lat_max", "lon_min", "lon_max"],
    )

    obs_dict = config["obs_dict"]

    # n_cpus/max_workers fall back to the job's actual CPU allocation when left
    # out of the config, so a PBS/Slurm resource request doesn't also have to be
    # duplicated into the JSON by hand. Explicit config values always win.
    default_workers = _default_worker_count()
    prepare_case_config = dict(config.get("prepare_case", {}))
    prepare_case_config.setdefault("n_cpus", default_workers)
    max_workers = config.get("max_workers", default_workers)
    print(
        f"Worker counts: prepare_case.n_cpus={prepare_case_config['n_cpus']}, "
        f"max_workers={max_workers} (job allocation reports {default_workers} "
        f"CPUs available; explicit config values, if set, take precedence)"
    )

    test_case = HistoryMatching(working_dir, case_name, mode=mode)

    if not case_existed:
        print("Case directory does not exist yet; creating and preparing it.")
        test_case.create_case(
            para, [ppe_tab, obs_tab], ppe_nc, obs_nc, obs_dict,
            lat_bins, manul_ppe_info, n_sample=config["n_sample"],
        )
        # Only safe to start the log now: create_case() just made case_root, and
        # doing this earlier would break its own "does this case exist" check.
        if mode == "python":
            _start_diagnostics_log(case_root, log_file_box)

        prepare_case_config["mode"] = mode
        test_case.prepare_case(prepare_case_config)
        test_case.load_mask(threshold_level=config["threshold_level"])

        for yname in config.get("visualize_check_vars", []):
            test_case.visualize_check(yname)
    else:
        print("Case directory already exists; loading it instead of re-preparing.")
        test_case.load_case()
        test_case.load_mask(threshold_level=config["threshold_level"])

    test_case.drop_by_name(config["vars_to_drop"])
    test_case.drop_by_n_survive(config["n_survive_threshold"])
    test_case.remove_var2d_auto(config["n_survive_threshold_2d"])

    test_case.prepare_for_sampling(max_workers=max_workers)
    test_case.draw(
        n_pts=config["n_pts"],
        n_threshold=config["n_threshold"],
        sample_threshold=config["sample_threshold"],
        max_workers=max_workers,
        n_max=config["n_max"],
    )

    test_case.save_samples_specifications(config["result_name"], top_n=config["top_n"])
    test_case.compare_with_original()

    print("Done.")


if __name__ == "__main__":
    main()
