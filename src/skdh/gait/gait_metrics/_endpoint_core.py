"""
Shared, key-agnostic gait endpoint implementations.

Single source of truth for the endpoint computations used by both
`skdh.gait.gait_metrics` and the deprecated `skdh.gait_old.gait_endpoints`. Each
function operates only on the arrays passed in and returns the per-step result array;
the endpoint classes in each module are responsible for reading/writing their own
dictionary keys (e.g. ``"stride time"`` vs ``"PARAM:stride time"``, ``"v axis"`` vs
``"vert axis"``).

``mask``/``mask_ofst`` are the boolean index arrays produced by
``GaitEventEndpoint._predict_init``. Per-bout functions return one value per step,
broadcast from the per-bout value via ``inertial_idx``.

Copyright (c) 2023. Pfizer Inc. All rights reserved.
"""

from numpy import (
    zeros,
    full,
    nanmean,
    nanstd,
    nanmedian,
    unique,
    sum,
    sqrt,
    nan,
    nonzero,
    argmin,
    abs,
    round,
    float64,
    fft,
    arange,
    isnan,
    maximum,
    moveaxis,
    ascontiguousarray,
)
from numpy.linalg import norm
from scipy.signal import butter, sosfiltfilt, find_peaks

from skdh.features.lib.extensions.statistics import autocorrelation
from skdh.features.lib.extensions.smoothness import SPARC


def _autocovariancefn(x, max_lag, biased=False, axis=0):
    y = ascontiguousarray(moveaxis(x, axis, -1))

    shape = list(y.shape)
    shape[-1] = max_lag
    ac = zeros(shape, dtype=float64)

    for i in range(min(max_lag, y.shape[-1] - 10)):
        ac[..., i] = autocorrelation(y, i, True)

    # if biased make sure its divided by just y.shape, instead of
    # y.shape - lag
    if biased:
        ac *= (y.shape[-1] - arange(0, max_lag)) / y.shape[-1]

    return moveaxis(ac, axis, -1)


def _nan_like(ic):
    return full(ic.size, nan, dtype=float64)


# ===========================================================
#     GAIT EVENT-LEVEL ENDPOINTS
# ===========================================================
def stride_time(*, fs, ic, mask, mask_ofst):
    """IC(i) to IC(i+2), in seconds. Use with offset 2 masks."""
    out = _nan_like(ic)
    out[mask] = (ic[mask_ofst] - ic[mask]) / fs
    return out


def stance_time(*, fs, ic, fc):
    """IC(i) to FC(i), in seconds."""
    return (fc - ic) / fs


def swing_time(*, fs, ic, fc, mask, mask_ofst):
    """FC(i) to IC(i+2), in seconds. Use with offset 2 masks."""
    out = _nan_like(ic)
    out[mask] = (ic[mask_ofst] - fc[mask]) / fs
    return out


def step_time(*, fs, ic, mask, mask_ofst):
    """IC(i) to IC(i+1), in seconds. Use with offset 1 masks."""
    out = _nan_like(ic)
    out[mask] = (ic[mask_ofst] - ic[mask]) / fs
    return out


def initial_double_support(*, fs, ic, fc_opp_foot):
    """IC(i) to FC of the opposite foot, in seconds."""
    return (fc_opp_foot - ic) / fs


def terminal_double_support(*, fs, ic, fc_opp_foot, mask, mask_ofst):
    """IC(i+1) to FC opp foot(i+1), in seconds. Use with offset 1 masks."""
    out = _nan_like(ic)
    out[mask] = (fc_opp_foot[mask_ofst] - ic[mask_ofst]) / fs
    return out


def double_support(*, initial_double_support, terminal_double_support):
    """Sum of initial and terminal double support."""
    return initial_double_support + terminal_double_support


def single_support(*, fs, ic, fc_opp_foot, mask, mask_ofst):
    """FC opp foot(i) to IC(i+1), in seconds. Use with offset 1 masks."""
    out = _nan_like(ic)
    out[mask] = (ic[mask_ofst] - fc_opp_foot[mask]) / fs
    return out


def inverted_pendulum_step_length(*, leg_length, delta_h):
    """Inverted pendulum model step length: 2 * sqrt(2 * l * h - h**2)."""
    return 2 * sqrt(2 * leg_length * delta_h - delta_h**2)


def stride_length(*, step_length, mask, mask_ofst):
    """Sum of two consecutive step lengths. Use with offset 1 masks."""
    out = _nan_like(step_length)
    out[mask] = step_length[mask_ofst] + step_length[mask]
    return out


def gait_speed(*, stride_length, stride_time):
    """Stride length divided by stride time."""
    return stride_length / stride_time


def cadence(*, step_time):
    """Steps per minute."""
    return 60.0 / step_time


def intra_covariance_v(*, ic, mask, mask_ofst, accel, inertial_idx, v_axis):
    """
    Autocovariance of the vertical acceleration over IC(i) to IC(i) + 2 * lag, at a
    lag of IC(i+k) - IC(i). Offset 1 masks give intra-step covariance, offset 2 masks
    give intra-stride covariance.
    """
    out = _nan_like(ic)

    i1 = ic[mask]
    i2 = ic[mask_ofst]
    i3 = i2 + (i2 - i1)

    for i, idx in enumerate(nonzero(mask)[0]):
        j_ = inertial_idx[idx]
        x = accel[j_][i1[i] : i3[i], v_axis[idx]]

        if (i3[i] - i1[i]) > x.size or x.size == 0:
            out[idx] = nan
        else:
            out[idx] = autocorrelation(x, i2[i] - i1[i], True)

    return out


intra_step_covariance_v = intra_covariance_v
intra_stride_covariance_v = intra_covariance_v


def harmonic_ratio_v(
    *,
    fs,
    ic,
    mask,
    mask_ofst,
    accel,
    inertial_idx,
    v_axis,
    stride_time,
    freq,
    harmonics,
    logger,
):
    """
    Per-stride vertical harmonic ratio. Use with offset 2 masks. `freq` is
    `fft.rfftfreq(1024)` and `harmonics` the harmonic numbers to use.
    """
    out = _nan_like(ic)

    i1 = ic[mask]
    i2 = ic[mask_ofst]

    for i, idx in enumerate(nonzero(mask)[0]):
        va = v_axis[idx]  # shorthand
        F = abs(
            fft.rfft(
                accel[inertial_idx[idx]][i1[i] : i2[i], va],
                n=1024,
            )
        )
        stridef = 1 / stride_time[idx]  # current stride frequency
        # get the indices for the first 20 harmonics
        ix_stridef = argmin(abs(freq * fs - stridef)) * harmonics
        if (ix_stridef < F.size).sum() <= 10:
            logger.warning(
                f"High stride frequency [{stridef:.2f}] results too few harmonics in "
                f"frequency range. Setting to nan"
            )
            out[idx] = nan
            continue
        elif (ix_stridef < F.size).sum() < 20:
            logger.warning(
                f"High stride frequency [{stridef:.2f}] results in use of less than 20 "
                f"harmonics [{(ix_stridef < F.size).sum()}]."
            )
        ix_stridef = ix_stridef[
            ix_stridef < F.size
        ]  # make sure not taking more than possible

        # index 1 is harmonic 2 -> even harmonics / odd harmonics
        out[idx] = sum(F[ix_stridef[1::2]]) / sum(F[ix_stridef[::2]])

    return out


def stride_sparc(*, fs, ic, mask, mask_ofst, accel, inertial_idx):
    """SPARC of the acceleration magnitude less gravity per stride. Offset 2 masks."""
    out = _nan_like(ic)

    i1 = ic[mask]
    i2 = ic[mask_ofst]

    for i, idx in enumerate(nonzero(mask)[0]):
        bout_i = inertial_idx[idx]

        if i2[i] - i1[i] > 0:
            out[idx] = SPARC(
                norm(accel[bout_i][i1[i] : i2[i], :], axis=1) - 1,
                fs,  # fsample
                4,  # padlevel
                10.0,  # fcut
                0.05,  # amplitude threshold
            )
        else:
            out[idx] = nan

    return out


# ===========================================================
#     GAIT BOUT-LEVEL ENDPOINTS
# ===========================================================
def phase_coordination_index(*, step_time, stride_time, n_bouts, inertial_idx):
    """Per-bout phase coordination index, broadcast to each step."""
    pci = zeros(n_bouts, dtype=float64)

    phase = step_time / stride_time  # %, not degrees
    for i in range(n_bouts):
        mask = inertial_idx == i

        psi_abs = nanmean(abs(phase[mask] - 0.5))  # using % not degrees right now
        psi_cv = nanstd(phase[mask], ddof=1) / nanmean(phase[mask])

        pci[i] = 100 * (psi_cv + psi_abs / 0.5)

    return pci[inertial_idx]


def gait_symmetry_index(*, fs, accel, stride_time, inertial_idx):
    """Per-bout gait symmetry index, broadcast to each step."""
    gsi = zeros(len(accel), dtype=float64)

    # setup acceleration filter if its possible to use
    if 0 < (2 * 10 / fs) < 1:
        sos = butter(4, 2 * 10 / fs, btype="low", output="sos")
    else:
        sos = None

    for i, acc in enumerate(accel):
        lag_ = nanmedian(stride_time[inertial_idx == i]) * fs
        if isnan(lag_):  # if only nan values in the bout
            gsi[i] = nan
            continue
        lag = int(round(lag_))
        # GSI uses biased autocovariance
        if sos is not None:
            ac = _autocovariancefn(
                sosfiltfilt(sos, acc, axis=0), int(4.5 * fs), biased=True, axis=0
            )
        else:
            ac = _autocovariancefn(acc, int(4.5 * fs), biased=True, axis=0)

        # C_stride is the sum of 3 axes
        pks, _ = find_peaks(sum(ac, axis=1))
        # find the closest peak to the computed ideal half stride lag
        try:
            t_stride = pks[argmin(abs(pks - lag))]
            idx = int(0.5 * t_stride)

            # maximum ensures no sqrt of negative numbers
            gsi[i] = sqrt(sum(maximum(ac[idx], 0))) / sqrt(3)
        except ValueError:
            gsi[i] = nan

    return gsi[inertial_idx]


def regularity_v(*, fs, accel, cycle_time, inertial_idx, v_axis):
    """
    Per-bout vertical-axis autocovariance at the peak nearest the median cycle lag,
    broadcast to each step. `cycle_time` = step time gives step regularity, stride time
    gives stride regularity.
    """
    reg = full(len(accel), nan, dtype=float64)

    for i in unique(inertial_idx):
        acc = accel[i]
        mask = inertial_idx == i

        va = v_axis[mask][0]
        lag_ = nanmedian(cycle_time[mask]) * fs
        if isnan(lag_):  # if only nan values in the bout
            reg[i] = nan
            continue
        lag = int(round(lag_))
        acf = _autocovariancefn(acc[:, va], int(4.5 * fs), biased=False, axis=0)
        pks, _ = find_peaks(acf)
        try:
            idx = pks[argmin(abs(pks - lag))]
            reg[i] = acf[idx]
        except ValueError:
            reg[i] = nan

    # broadcast regularity into gait for each step
    return reg[inertial_idx]


step_regularity_v = regularity_v
stride_regularity_v = regularity_v


def autocovariance_symmetry(*, step_regularity, stride_regularity):
    """Absolute difference between step and stride regularity."""
    return abs(step_regularity - stride_regularity)


def regularity_index(*, step_regularity, stride_regularity):
    """1 - 2|R_stride - R_step| / (R_step + R_stride)."""
    return 1 - (
        2
        * abs(stride_regularity - step_regularity)
        / (step_regularity + stride_regularity)
    )
