"""Statistics unit tests on synthetic data (CPU, < 1 min): the Mantel engine, its nulls, the bootstrap and retrieval."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "step2"))
import compare_tiers as ct  # noqa: E402

N, K = 60, 6


def synthetic(shared_unit, shared_class, seed=0):
    """Two 'modalities' for N compounds in K classes; they share per-compound and/or class-level structure."""
    rng = np.random.default_rng(seed)
    labels = np.repeat(np.arange(K), N // K).astype(str)
    cls = rng.normal(size=(K, 8))
    lat = rng.normal(size=(N, 8))
    z1 = (lat if shared_unit else rng.normal(size=(N, 8))) + (cls[labels.astype(int)] if shared_class else 0)
    z2 = (lat if shared_unit else rng.normal(size=(N, 8))) + (cls[labels.astype(int)] if shared_class else 0)
    A = ct.cos(z1 @ rng.normal(size=(8, 30)) + rng.normal(size=(N, 30)))
    B = ct.cos(z2 @ rng.normal(size=(8, 40)) + rng.normal(size=(N, 40)))
    return A, B, labels


def test_rank_fast_path_equals_reranking():
    A, B, _ = synthetic(True, False)
    fast = ct.mantel(A, B, perms=19)["r"]
    slow = ct.mantel(A, B, perms=19, mask=np.ones_like(A, bool))["r"]
    assert abs(fast - slow) < 1e-9


def test_mantel_detects_shared_structure():
    A, B, lab = synthetic(True, False)
    res = ct.mantel(A, B, perms=199)
    assert res["r"] > 0.2 and res["p"] <= 0.01


def test_within_class_null_is_calibrated_when_only_classes_are_shared():
    # plain null is fooled by shared class structure; the within-class null must not be
    A, B, lab = synthetic(False, True, seed=1)
    assert ct.mantel(A, B, perms=199)["p"] <= 0.01
    assert ct.mantel(A, B, perms=199, strata=lab)["p"] > 0.05


def test_within_class_null_still_detects_unit_level_agreement():
    A, B, lab = synthetic(True, True, seed=2)
    assert ct.mantel(A, B, perms=199, strata=lab)["p"] <= 0.01


def test_bootstrap_ci_brackets_estimate():
    A, B, _ = synthetic(True, False, seed=3)
    r = ct.mantel(A, B, perms=9)["r"]
    lo, hi = ct.boot_ci(A, B, n_boot=100)
    assert lo < r < hi


def test_moa_retrieval_beats_label_permutation():
    A, _, lab = synthetic(False, True, seed=4)
    res = ct.nn_retrieval(A, lab, perms=99)
    assert res["mAP"] > res["mAP_null"] and res["p_mAP"] <= 0.05


def test_testability_rule():
    assert ct.testable(np.array(["a"] * 3 + ["b"] * 3 + ["c"] * 3 + ["unclear"] * 5))[0]
    assert not ct.testable(np.array(["a"] * 3 + ["b"] * 3 + ["c"] * 2))[0]
