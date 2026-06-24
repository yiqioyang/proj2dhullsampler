def test_public_api_imports():
    from proj2dhullsampler import HistoryMatching

    assert callable(HistoryMatching)


def test_all_exports_are_importable():
    import proj2dhullsampler

    assert proj2dhullsampler.__all__ == ["HistoryMatching"]
    for name in proj2dhullsampler.__all__:
        assert getattr(proj2dhullsampler, name) is not None
