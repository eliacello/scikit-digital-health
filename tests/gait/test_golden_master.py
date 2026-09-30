"""
Golden-master regression tests for the `skdh.gait` (GaitLumbar) pipeline.

Every output array of `gait.GaitLumbar.predict` (all event- and bout-level endpoints,
including the Model1/Model2 spatial chains and their asymmetry values, plus the
bookkeeping columns) is compared against a committed fixture generated from the
unmodified code, at rtol = atol = 1e-10. Key sets, dtypes, shapes, and NaN placement
must also match.

The fixture must only be regenerated with explicit QA approval:

    python tests/gait/test_golden_master.py --regenerate
"""

import sys
import warnings
from pathlib import Path

import pytest
from numpy import load, savez_compressed, asarray
from numpy.testing import assert_allclose, assert_array_equal

from skdh.utility.internal import apply_resample
from skdh.gait.core import GaitLumbar
from skdh.context import PredictGaitLumbarLgbm

RTOL = 1e-10
ATOL = 1e-10
DATA_DIR = Path(__file__).parent / "data"
FIXTURE = DATA_DIR / "golden_master_gait.npz"
SEP = "::"

APCWT_KW = dict(
    downsample=False,
    height_factor=0.53,
    min_bout_time=8.0,
    max_bout_separation_time=0.5,
    gait_event_method="AP CWT",
    correct_orientation=True,
    filter_cutoff=20.0,
    filter_order=4,
    ic_prom_factor=0.1,
    ic_dist_factor=0.0,
    fc_prom_factor=0.1,
    fc_dist_factor=0.0,
    max_stride_time=lambda x: 2.0 * x + 1.0,
    loading_factor=lambda x: 0.17 * x + 0.05,
)

VCWT_KW = dict(
    downsample=True,
    height_factor=0.53,
    min_bout_time=8.0,
    max_bout_separation_time=0.5,
    gait_event_method="v cwt",
    correct_orientation=True,
    filter_cutoff=20.0,
    filter_order=4,
    use_cwt_scale_relation=True,
    wavelet_scale="default",
    round_scale=True,
    max_stride_time=2.25,
    loading_factor=0.2,
)

VCWT_FIXED_SCALE_KW = dict(VCWT_KW, use_cwt_scale_relation=False, wavelet_scale=8)


def _bouts_50hz():
    data = load(DATA_DIR / "gait_input.npz")
    t50, (acc50,) = apply_resample(
        goal_fs=50.0, time=data["time"], data=(data["accel"],), indices=(), aa_filter=True
    )
    cls = PredictGaitLumbarLgbm()
    cls._in_pipeline = True
    kw, _ = cls.predict(time=t50, accel=acc50, fs=50.0)
    return kw


def _gyro_128hz(gait_kw):
    data = load(DATA_DIR / "gait_input2.npz")
    return GaitLumbar(**gait_kw).predict(
        time=data["time"],
        accel=data["accel"],
        gyro=data["gyro"],
        fs=128.0,
        height=1.88,
        gait_pred=True,
    )


CASES = {
    "apcwt_50hz": lambda: GaitLumbar(**APCWT_KW).predict(height=1.88, **_bouts_50hz()),
    "vcwt_50hz": lambda: GaitLumbar(**VCWT_KW).predict(height=1.88, **_bouts_50hz()),
    "vcwt_50hz_no_height": lambda: GaitLumbar(**VCWT_KW).predict(
        height=None, **_bouts_50hz()
    ),
    "apcwt_128hz_turns": lambda: _gyro_128hz(APCWT_KW),
    "vcwt_128hz_turns": lambda: _gyro_128hz(VCWT_FIXED_SCALE_KW),
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
