"""Config-driven run of the full pipeline (the entry point used by submit_apply.pbs).

Usage (from application/):
    python ../proj2dhullsampler/run_apply.py --config config_table.json

All paths, thresholds, and the "notebook" vs "python" mode flag live in the
config file instead of being hardcoded, so a run can be reproduced or tuned
without touching this script. See application/config_table.json and
application/config_nc.json for examples, and application/config_annotated.jsonc
for what every key means.

In "python" mode, checkup figures that would normally display inline in the
notebook are instead written as PNGs to <working_dir>/<case_name>/diagnostics/,
and all print() output from the pipeline is additionally teed to
<working_dir>/<case_name>/diagnostics/run_log.txt.
"""

import argparse
import json
import sys
import traceback
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


def _write_constraint_diagnostics(test_case, opts):
    """Per-pair constraint animations, the interlock PDF and the dropped-variable
    pages, into <case>/diagnostics/.

    Called only after the samples have been written - which the interlock figure
    also *requires*, since its right-hand column is the drawn samples themselves
    - and each figure is guarded on its own: these are diagnostics, so a
    plotting failure must never take down a job whose actual results are already
    safely on disk.

    `opts` is the optional "constraint_diagnostics" block of the config; an
    empty dict writes the animations, the interlock PDF and the dropped-variable
    pages with default settings. See application/config_annotated.jsonc for the
    fields.
    """
    from proj2dhullsampler.history_matching_animation import (
        animate_pair_constraints,
        plot_constraint_interlock,
        plot_dropped_constraints,
        save_pair_frames,
    )

    opts = dict(opts or {})

    if not opts.get("enabled", True):
        print("constraint_diagnostics disabled in config; skipping.")
        return

    def _guarded(label, fn):
        print(f"--- constraint diagnostics: {label} ---")
        try:
            fn()
        except Exception:
            # Deliberately swallowed, but never silently: the samples are
            # already saved, so this must not fail the job.
            print(f"WARNING: {label} failed; continuing. Traceback follows.")
            traceback.print_exc(file=sys.stdout)

    def _kw(*names):
        """Config keys -> keyword arguments, skipping anything not set.

        A name may be given as "config_key:argument_name" where the two differ
        (the interlock PDF wants its own point budget, separate from the
        animations').
        """
        out = {}
        for name in names:
            key, _, arg = name.partition(":")
            if key in opts:
                out[arg or key] = opts[key]
        return out

    if opts.get("animations", True):
        _guarded(
            "per-pair animations",
            lambda: animate_pair_constraints(
                test_case,
                mode="python",  # run_apply.py is batch: always write files
                **_kw("fmt", "fps", "dpi", "max_points", "var_order", "show_hull"),
            ),
        )

    if opts.get("frames", False):
        _guarded(
            "per-pair still frames",
            lambda: save_pair_frames(
                test_case, **_kw("dpi", "max_points", "var_order", "show_hull")
            ),
        )

    if opts.get("interlock", True):
        def _interlock():
            import matplotlib.pyplot as plt

            fig = plot_constraint_interlock(
                test_case,
                **_kw(
                    "n_rows",
                    "interlock_max_points:max_points",
                    "panel_size",
                    "show_hull",
                    "sort_by",
                ),
            )
            plt.close(fig)

        _guarded("constraint interlock PDF", _interlock)

    if opts.get("dropped", True):
        _guarded(
            "dropped-variable constraints",
            lambda: plot_dropped_constraints(
                test_case,
                mode="python",  # batch: always write files
                **_kw(
                    "dropped_ncols:ncols",
                    "dropped_rows_per_page:rows_per_page",
                    "dropped_max_vars_per_panel:max_vars_per_panel",
                    "show_hull",
                ),
            ),
        )


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

    # Last: everything above has already been written, so these figures are
    # free to fail without costing the run anything.
    _write_constraint_diagnostics(test_case, config.get("constraint_diagnostics", {}))

    print("Done.")


if __name__ == "__main__":
    main()
