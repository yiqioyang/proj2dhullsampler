def test_public_api_imports():
    from proj2dhullsampler import HistoryMatching

    assert callable(HistoryMatching)


def test_all_exports_are_importable():
    import proj2dhullsampler

    assert proj2dhullsampler.__all__ == ["HistoryMatching"]
    for name in proj2dhullsampler.__all__:
        assert getattr(proj2dhullsampler, name) is not None


def test_every_submodule_imports_cleanly():
    """Regression test: a bad module-level default (e.g. the removed
    ``np.NAN``, which NumPy >=2.0 no longer defines) raises at import time,
    not when the buggy function is called. Import every submodule directly
    so such errors surface immediately instead of only when a user happens
    to exercise the one broken function."""
    import importlib

    submodules = [
        "aux",
        "hm_class",
        "plotting",
        "prep_class",
        "preprocess",
        "sampling_functions",
        "utils",
    ]
    for name in submodules:
        importlib.import_module(f"proj2dhullsampler.{name}")
