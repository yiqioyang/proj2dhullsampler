"""Public API for the proj2dhullsampler package."""

from importlib import import_module

_EXPORTS = {
    "HistoryMatching": (".hm_class", "HistoryMatching"),
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
