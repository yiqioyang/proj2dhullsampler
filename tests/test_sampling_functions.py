import pandas as pd
from shapely.geometry import Polygon

from proj2dhullsampler.sampling_functions import sample_from_hull


def test_sample_from_hull_keeps_only_points_inside_polygon():
    square = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])

    X = pd.DataFrame(
        {
            "p1": [0.5, 2.0, 0.1],
            "p2": [0.5, 2.0, 0.9],
        }
    )

    result = sample_from_hull(X, ("p1", "p2"), square)

    assert list(result["p1"]) == [0.5, 0.1]
    assert list(result["p2"]) == [0.5, 0.9]


def test_sample_from_hull_returns_empty_when_nothing_matches():
    square = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])

    X = pd.DataFrame({"p1": [5.0], "p2": [5.0]})

    result = sample_from_hull(X, ("p1", "p2"), square)

    assert result.empty
