import pandas as pd

from proj2dhullsampler.hm_class import HistoryMatching, meta_one_hot_shot


def test_meta_one_hot_shot_marks_selected_parameter_indices():
    # columns are diagnostic names; values are the sensitive parameter
    # indices selected for that diagnostic (as produced by
    # Prepare_Case.sensitivity_emulation).
    meta = pd.DataFrame({"y1": [0, 2], "y2": [1, 2]})
    para_nm = ["p0", "p1", "p2", "p3"]

    one_hot = meta_one_hot_shot(meta, para_nm)

    assert list(one_hot.loc["y1"]) == [True, False, True, False]
    assert list(one_hot.loc["y2"]) == [False, True, True, False]


def test_rescale_para_maps_normalized_samples_back_to_original_range():
    hm = HistoryMatching(working_dir="/tmp/unused", case_name="unused")
    hm.ppe_para = pd.DataFrame({"a": [0.0, 10.0], "b": [-5.0, 5.0]})

    sampled = pd.DataFrame({"a": [0.0, 1.0], "b": [0.0, 1.0]})
    rescaled = hm.rescale_para(sampled)

    assert list(rescaled["a"]) == [0.0, 10.0]
    assert list(rescaled["b"]) == [-5.0, 5.0]
