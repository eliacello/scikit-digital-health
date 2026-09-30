"""
Cross-module endpoint comparison: `skdh.gait_old.gait_endpoints` vs
`skdh.gait.gait_metrics`.

For every endpoint that exists in both modules, run both implementations on the same
per-step inputs (translated to each module's key convention) and record whether the
outputs are identical. This is a report only; it does not assert anything.

Inputs are captured from each pipeline immediately before endpoint computation, plus
the synthetic unit-test fixture. For the spatial chain (old `StepLength` /
`StrideLength` / `GaitSpeed` vs new `*Model1`), two modes are reported:

- ``shared delta h``: the new module is given the old module's ``delta h`` as
  ``m1 delta h``, so only the endpoint formulas are compared.
- ``native delta h``: each module uses its own delta-h (old: ``get_strides``,
  new: ``DeltaHModel1``). Only possible for the gait_old pipeline input, which is
  the only input carrying a native ``delta h``.

Usage:

    python tests/gait/cross_module_report.py [output.md]
"""

import sys
import json
import warnings
from copy import deepcopy
from pathlib import Path

from numpy import (
    load,
    array,
    arange,
    zeros,
    sin,
    pi,
    asarray,
    array_equal,
    isnan,
    nanmax,
    abs as np_abs,
    allclose,
)

from skdh.utility.internal import apply_resample
from skdh.gait_old import gait_endpoints as old
from skdh.gait_old.gait import Gait
from skdh.gait.gait_metrics import gait_metrics as new
from skdh.gait.core import GaitLumbar
from skdh.context import PredictGaitLumbarLgbm

TESTS = Path(__file__).parents[1]

PAIRS = [
    # (old class, new class)
    ("StrideTime", "StrideTime"),
    ("StanceTime", "StanceTime"),
    ("SwingTime", "SwingTime"),
    ("StepTime", "StepTime"),
    ("InitialDoubleSupport", "InitialDoubleSupport"),
    ("TerminalDoubleSupport", "TerminalDoubleSupport"),
    ("DoubleSupport", "DoubleSupport"),
    ("SingleSupport", "SingleSupport"),
    ("StepLength", "StepLengthModel1"),
    ("StrideLength", "StrideLengthModel1"),
    ("GaitSpeed", "GaitSpeedModel1"),
    ("Cadence", "Cadence"),
    ("IntraStepCovarianceV", "IntraStepCovarianceV"),
    ("IntraStrideCovarianceV", "IntraStrideCovarianceV"),
    ("HarmonicRatioV", "HarmonicRatioV"),
    ("StrideSPARC", "StrideSPARC"),
    ("PhaseCoordinationIndex", "PhaseCoordinationIndex"),
    ("GaitSymmetryIndex", "GaitSymmetryIndex"),
    ("StepRegularityV", "StepRegularityV"),
    ("StrideRegularityV", "StrideRegularityV"),
    ("AutocovarianceSymmetryV", "AutocovarianceSymmetryV"),
    ("RegularityIndexV", "RegularityIndexV"),
]
SPATIAL = {"StepLength", "StrideLength", "GaitSpeed"}
BASE_KEYS = ["IC", "FC", "FC opp foot", "Bout N", "forward cycles"]


class _Capture:
    store = {}

    def predict(self, *args, **kwargs):
        if args:
            fs, leg_length, gait, gait_aux = args
        else:
            fs, leg_length = kwargs["fs"], kwargs["leg_length"]
            gait, gait_aux = kwargs["gait"], kwargs["gait_aux"]
        _Capture.store = dict(
            fs=fs,
            leg_length=leg_length,
            gait=deepcopy(gait),
            gait_aux=deepcopy(gait_aux),
        )


def capture_old_pipeline():
    data = load(TESTS / "gait_old" / "data" / "gait_input.npz")
    t, (a,) = apply_resample(
        goal_fs=50.0, time=data["time"], data=(data["accel"],), indices=(), aa_filter=True
    )
    g = Gait(max_bout_separation_time=0.25)
    g._params = [_Capture]
    g.predict(time=t, accel=a, fs=50.0, height=1.88)
    s = _Capture.store
    gait = {k: s["gait"][k] for k in BASE_KEYS + ["delta h"]}
    aux = dict(
        accel=s["gait_aux"]["accel"],
        axis=s["gait_aux"]["vert axis"],
        idx=s["gait_aux"]["inertial data i"],
    )
    return s["fs"], s["leg_length"], gait, aux


def capture_new_pipeline():
    data = load(TESTS / "gait" / "data" / "gait_input.npz")
    t, (a,) = apply_resample(
        goal_fs=50.0, time=data["time"], data=(data["accel"],), indices=(), aa_filter=True
    )
    cls = PredictGaitLumbarLgbm()
    cls._in_pipeline = True
    kw, _ = cls.predict(time=t, accel=a, fs=50.0)
    g = GaitLumbar()
    g._params = [_Capture]
    g.predict(height=1.88, **kw)
    s = _Capture.store
    gait = {k: s["gait"][k] for k in BASE_KEYS}
    aux = dict(
        accel=s["gait_aux"]["accel"],
        axis=s["gait_aux"]["v axis"],
        idx=s["gait_aux"]["inertial data i"],
    )
    # the new pipeline has no native delta h: use DeltaHModel1 as the shared value
    tmp = to_new(gait, aux, None)
    new.DeltaHModel1().predict(fs=s["fs"], leg_length=s["leg_length"], **tmp)
    gait["delta h"] = tmp["gait"]["m1 delta h"]
    return s["fs"], s["leg_length"], gait, aux


def synthetic():
    t = arange(0, 20, 0.02)
    a = zeros((t.size, 3))
    a[:, 0] = 1 + 0.75 * sin(2 * pi * 1.0 * t)
    a[:, 1] = 0.35 * sin(2 * pi * 1.0 * t)
    a[:, 2] = 0.25 * sin(2 * pi * 2.0 * t + pi / 2)
    gait = {
        "IC": array([50, 100, 150, 200, 250, 300, 350, 400]),
        "FC": array([110, 165, 210, 265, 310, 365, 415, 455]),
        "FC opp foot": array([65, 110, 165, 210, 260, 315, 360, 405]),
        "delta h": array([0.1, 0.2, 0.1, 0.2, 0.2, 0.2, 0.1, 0.1]),
        "Bout N": array([1, 1, 1, 2, 2, 2, 2, 2]),
        "forward cycles": array([2, 1, 0, 2, 2, 2, 1, 0]),
    }
    aux = dict(accel=[a[:199], a], axis=array([0] * 8), idx=array([0, 0, 0, 1, 1, 1, 1, 1]))
    return 50.0, 1.8, gait, aux


def to_old(gait, aux):
    g = {k: asarray(v).copy() for k, v in gait.items()}
    return dict(
        gait=g,
        gait_aux={
            "accel": [x.copy() for x in aux["accel"]],
            "vert axis": asarray(aux["axis"]).copy(),
            "inertial data i": asarray(aux["idx"]).copy(),
        },
    )


def to_new(gait, aux, delta_h_mode):
    g = {k: asarray(v).copy() for k, v in gait.items() if k != "delta h"}
    if delta_h_mode == "shared":
        g["m1 delta h"] = asarray(gait["delta h"]).copy()
    return dict(
        gait=g,
        gait_aux={
            "accel": [x.copy() for x in aux["accel"]],
            "v axis": asarray(aux["axis"]).copy(),
            "inertial data i": asarray(aux["idx"]).copy(),
        },
    )


def compare(a, b):
    a = asarray(a, dtype=float)
    b = asarray(b, dtype=float)
    if a.shape != b.shape:
        return dict(identical=False, close_1e10=False, max_abs_diff=None, n=None,
                    note=f"shape {a.shape} vs {b.shape}")
    nan_match = bool(array_equal(isnan(a), isnan(b)))
    ok = ~(isnan(a) | isnan(b))
    mad = float(nanmax(np_abs(a[ok] - b[ok]))) if ok.any() else 0.0
    return dict(
        identical=bool(array_equal(a, b, equal_nan=True)),
        close_1e10=bool(nan_match and allclose(a, b, rtol=1e-10, atol=1e-10, equal_nan=True)),
        max_abs_diff=mad,
        nan_placement_equal=nan_match,
        n=int(a.size),
    )


def run_pair(old_name, new_name, fs, leg_length, gait, aux, delta_h_mode):
    oc = getattr(old, old_name)()
    nc = getattr(new, new_name)()
    o = to_old(gait, aux)
    n = to_new(gait, aux, delta_h_mode)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        oc.predict(fs, leg_length, o["gait"], o["gait_aux"])
        nc.predict(fs=fs, leg_length=leg_length, gait=n["gait"], gait_aux=n["gait_aux"])
    rows = [(oc.k_, nc.k_, compare(o["gait"][oc.k_], n["gait"][nc.k_]))]
    oa, na = f"{oc.k_} asymmetry", f"{nc.k_} asymmetry"
    if oa in o["gait"] or na in n["gait"]:
        if oa in o["gait"] and na in n["gait"]:
            rows.append((oa, na, compare(o["gait"][oa], n["gait"][na])))
        else:
            rows.append((oa, na, dict(identical=False, close_1e10=False,
                                      note="asymmetry only present in one module")))
    return rows


def main(out=None):
    inputs = {
        "gait_old pipeline (tests/gait_old/data/gait_input.npz)": capture_old_pipeline(),
        "gait pipeline (tests/gait/data/gait_input.npz)": capture_new_pipeline(),
        "synthetic unit-test fixture": synthetic(),
    }
    results = []
    for in_name, (fs, ll, gait, aux) in inputs.items():
        for old_name, new_name in PAIRS:
            modes = ["shared"]
            if old_name in SPATIAL and in_name.startswith("gait_old"):
                modes.append("native")
            for mode in modes:
                for ok, nk, cmp in run_pair(old_name, new_name, fs, ll, gait, aux, mode):
                    results.append(dict(
                        input=in_name,
                        old_class=old_name,
                        new_class=new_name,
                        delta_h=mode if old_name in SPATIAL else "n/a",
                        old_key=ok,
                        new_key=nk,
                        **cmp,
                    ))

    # delta h itself: gait_old get_strides vs DeltaHModel1, on the gait_old pipeline input
    fs, ll, gait, aux = inputs["gait_old pipeline (tests/gait_old/data/gait_input.npz)"]
    n = to_new(gait, aux, None)
    new.DeltaHModel1().predict(fs=fs, leg_length=ll, **n)
    dh = compare(gait["delta h"], n["gait"]["m1 delta h"])
    results.append(dict(
        input="gait_old pipeline (tests/gait_old/data/gait_input.npz)",
        old_class="gait_old.get_strides (internal)", new_class="DeltaHModel1",
        delta_h="native", old_key="delta h", new_key="m1 delta h", **dh,
    ))

    lines = [
        "# Pre-refactor cross-module endpoint comparison",
        "",
        "gait_old.gait_endpoints vs gait.gait_metrics on identical per-step inputs.",
        "`identical` = bit-for-bit equal (NaN-aware); `close` = allclose at rtol=atol=1e-10.",
        "",
        "| input | old endpoint | new endpoint | delta h | old key | new key | identical | close 1e-10 | max abs diff | note |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        mad = r.get("max_abs_diff")
        note = r.get("note", "")
        if r.get("nan_placement_equal") is False:
            note = (note + " NaN placement differs").strip()
        lines.append(
            f"| {r['input']} | {r['old_class']} | {r['new_class']} | {r['delta_h']} | "
            f"`{r['old_key']}` | `{r['new_key']}` | {r['identical']} | {r['close_1e10']} | "
            f"{'' if mad is None else f'{mad:.3e}'} | {note} |"
        )
    text = "\n".join(lines) + "\n"
    if out:
        Path(out).write_text(text)
        Path(out).with_suffix(".json").write_text(json.dumps(results, indent=1))
    print(text)
    return results


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
