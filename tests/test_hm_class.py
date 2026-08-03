import pandas as pd

from proj2dhullsampler.hm_class import EmulatedDataStorage, HistoryMatching, meta_one_hot_shot


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


def test_remove_var2d_auto_relaxes_threshold_for_3way_interaction(tmp_path, capsys):
    # A, B, C all share the same sensitive parameter pair (px, py). Every
    # *pairwise* overlap among them is healthy (>= 12), but requiring all
    # three to pass together drops to 8 - a 3+-way interaction effect that
    # shuffle_vars's pairwise-only check can't see directly. This used to
    # crash remove_var2d_auto (empty pairwise selection -> IndexError);
    # it should now relax the pairwise threshold via added_num each time
    # nothing qualifies, until a pair does, and resolve the group.
    hm = HistoryMatching(working_dir=tmp_path, case_name="case", mode="notebook")
    (hm.root / "output").mkdir(parents=True)
    hm.dropped_vars = EmulatedDataStorage()
    hm.dropped_vars.nooverlap2d = []
    hm.specifications = EmulatedDataStorage()

    hm.para_nm = ["px", "py"]
    hm.meta = pd.DataFrame({"A": [0, 1], "B": [0, 1], "C": [0, 1]})
    hm.var_nm = ["A", "B", "C"]

    seg1, seg2, seg3, seg4 = 8, 4, 5, 6  # row counts per segment below
    hm.tf_masks = pd.DataFrame({
        "A": [True] * seg1 + [True] * seg2 + [True] * seg3 + [False] * seg4,
        "B": [True] * seg1 + [True] * seg2 + [False] * seg3 + [True] * seg4,
        "C": [True] * seg1 + [False] * seg2 + [True] * seg3 + [True] * seg4,
    })
    # A&B&C = 8 (< threshold); A&B = 12, A&C = 13, B&C = 14 (all >= threshold)

    hm.remove_var2d_auto(overlapping_threshold=10, added_num=3)

    assert "Need to increase threshold" in capsys.readouterr().out
    assert len(hm.dropped_vars.nooverlap2d) == 1
    assert hm.dropped_vars.nooverlap2d[0] in {"A", "B"}
    assert hm.paras_vars_0 == {}
