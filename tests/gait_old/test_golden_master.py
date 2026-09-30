"""
Golden-master regression tests for the deprecated `skdh.gait_old` pipeline.

Every output array of `gait_old.Gait.predict` (all event- and bout-level
endpoints, their asymmetry values, and the bookkeeping columns) is compared against a
committed fixture generated from the unmodified code, at rtol = atol = 1e-10. Key sets,
dtypes, shapes, and NaN placement must also match.

The fixture must only be regenerated with explicit QA approval:

    python tests/gait_old/test_golden_master.py --regenerate
"""

import sys
import warnings
from pathlib import Path

import pytest
from numpy import load, savez_compressed, asarray
from numpy.testing import assert_allclose, assert_array_equal

from skdh.utility.internal import apply_resample
from skdh.gait_old.gait import Gait

RTOL = 1e-10
ATOL = 1e-10
DATA_DIR = Path(__file__).parent / "data"
FIXTURE = DATA_DIR / "golden_master_gait_old.npz"
SEP = "::"

GAIT_KW = dict(
    correct_accel_orient=True,
    use_cwt_scale_relation=True,
    min_bout_time=8.0,
    max_bout_separation_time=0.25,
    max_stride_time=2.25,
    loading_factor=0.2,
    height_factor=0.53,
    prov_leg_length=False,
    filter_order=4,
    filter_cutoff=20.0,
    downsample_aa_filter=True,
    day_window=(0, 24),
)


def _input_50hz():
    data = load(DATA_DIR / "gait_input.npz")
    t50, (acc50,) = apply_resample(
        goal_fs=50.0, time=data["time"], data=(data["accel"],), indices=(), aa_filter=True
    )
    return t50, acc50


def _case_accel_50hz():
    t, acc = _input_50hz()
    return Gait(**GAIT_KW).predict(time=t, accel=acc, fs=50.0, height=1.88)


def _case_accel_50hz_no_height():
    t, acc = _input_50hz()
    return Gait(**GAIT_KW).predict(time=t, accel=acc, fs=50.0, height=None)


def _case_gyro_128hz_turns():
    data = load(DATA_DIR / "gait_input2.npz")
    return Gait(**GAIT_KW).predict(
        time=data["time"],
        accel=data["accel"],
        gyro=data["gyro"],
        fs=128.0,
        height=1.88,
        gait_pred=True,
    )


CASES = {
    "accel_50hz": _case_accel_50hz,
    "accel_50hz_no_height": _case_accel_50hz_no_height,
    "gyro_128hz_turns": _case_gyro_128hz_turns,
}


def _array_outputs(res):
    """Keep every bool/int/float/str/datetime output array (no pickled objects)."""
    out = {}
    for key, value in res.items():
        arr = asarray(value)
        if arr.dtype.kind in "biufUM":
            out[key] = arr
    return out


def run_case(case):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return _array_outputs(CASES[case]())


def regenerate():
    arrays = {}
    for case in CASES:
        for key, arr in run_case(case).items():
            arrays[f"{case}{SEP}{key}"] = arr
    savez_compressed(FIXTURE, **arrays)
    print(f"wrote {len(arrays)} arrays to {FIXTURE}")


@pytest.fixture(scope="module")
def golden():
    with load(FIXTURE, allow_pickle=False) as f:
        return {k: f[k] for k in f.files}


@pytest.mark.parametrize("case", list(CASES))
def test_golden_master(case, golden):
    expected = {
        k.split(SEP, 1)[1]: v for k, v in golden.items() if k.startswith(case + SEP)
    }
    assert expected, f"no golden-master arrays stored for case {case!r}"

    result = run_case(case)

    assert sorted(result) == sorted(expected), "output key set changed"
    for key, exp in expected.items():
        res = result[key]
        assert res.dtype == exp.dtype, f"{key}: dtype {res.dtype} != {exp.dtype}"
        assert res.shape == exp.shape, f"{key}: shape {res.shape} != {exp.shape}"
        if exp.dtype.kind == "f":
            assert_allclose(res, exp, rtol=RTOL, atol=ATOL, equal_nan=True, err_msg=key)
        else:
            assert_array_equal(res, exp, err_msg=key)


if __name__ == "__main__":
    if "--regenerate" in sys.argv:
        regenerate()
    else:
        print(__doc__)
