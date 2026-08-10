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
import sys
from pathlib import Path


def load_config(config_path):
    with open(config_path) as f:
        return json.load(f)


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
    from proj2dhullsampler.pipeline import build_case, default_worker_count

    case_root = working_dir / case_name
    case_existed = case_root.exists()

    if mode == "python" and case_existed:
        _start_diagnostics_log(case_root, log_file_box)

    on_created = None
    if mode == "python":
        # Only safe to start the log once create_case() has made case_root -
        # doing it earlier would make HistoryMatching.create_case()'s own
        # "does this case already exist" check see a directory that we
        # ourselves made.
        on_created = lambda test_case: _start_diagnostics_log(case_root, log_file_box)

    test_case = build_case(config, mode=mode, on_created=on_created)

    # max_workers falls back to the job's actual CPU allocation when left out
    # of the config, so a PBS/Slurm resource request doesn't also have to be
    # duplicated into the JSON by hand. An explicit config value always wins.
    max_workers = config.get("max_workers", default_worker_count())
    print(f"max_workers={max_workers}")

    

    test_case.drop_by_name(config["vars_to_drop"])
    test_case.drop_by_emulator_performance(config["emultor_error_ratio_threshold"])

    test_case.drop_by_n_survive(config["n_survive_threshold"])
    test_case.remove_var2d_auto(config["n_survive_threshold_2d"], added_num = config['added_number_for_pairs'])
    test_case.drop_by_nvar_per_pair(config['n_var_thre'])

    test_case.prepare_for_sampling(
        max_workers=max_workers,
        threshold_ratio_between_para_pairs = config['threshold_ratio_between_para_pairs']
    )
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
