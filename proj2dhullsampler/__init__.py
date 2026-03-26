"""Public API for the proj2dhullsampler package.

Use this module for one-line imports in notebooks/scripts, e.g.:
    from proj2dhullsampler import Analysis, HistoryMatching, Prep_Mask_Generation
"""

from importlib import import_module

_EXPORTS = {
    "Analysis": (".analysis", "Analysis"),
    "HistoryMatching": (".hm_class", "HistoryMatching"),
    "Prep_Mask_Generation": (".prep_class", "Prep_Mask_Generation"),
    "biplot": (".plotting", "biplot"),
    "biplot_original_scale": (".plotting", "biplot_original_scale"),
    "gp_training_application": (".utils", "gp_training_application"),
    "metric_cal_single": (".aux", "metric_cal_single"),
    "orchestrate_test": (".sampling_functions", "orchestrate_test"),
    "para_csv2nc": (".aux", "para_csv2nc"),
    "plot_histograms_grid_5": (".plotting", "plot_histograms_grid_5"),
    "sample_from_hull": (".sampling_functions", "sample_from_hull"),
    "sample_from_hulls_n": (".sampling_functions", "sample_from_hulls_n"),
}

__all__ = list(_EXPORTS)


def __getattr__(name):
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module_name, attr_name = _EXPORTS[name]
    module = import_module(module_name, __name__)
    value = getattr(module, attr_name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(list(globals()) + __all__)
