#!/usr/bin/env python3
"""Governed formula layer: bind acoustic validation checks into the Lambda-spine.

Ports the szl.lambda/v1 contract (canonical source:
szl-holdings/szl-lambda-gate reference/szl_lambda_v1.py) as a minimal stdlib
implementation, then aggregates the Aero-FF / Burgers / Duffing check scores
through it.

    Lambda_w(x) = prod_k x_k^w_k = exp(fsum_k(w_k * log x_k)), Lambda = 0 if some x_k == 0

Gate verdicts: GO / NO_GO / ABSTAIN / BLOCK (tau required, no implicit default).
Any failed check is a ZERO_VETO — one failing formula sinks the aggregate, by
design. Lambda is ADVISORY. Lambda uniqueness is Conjecture 1 (machine-checked
false as stated in szl-papers); nothing here depends on it, and this module
never calls Lambda a theorem or a proof.
"""
from __future__ import annotations
import hashlib
import json
import math
import sys
from typing import Any

SCHEMA = "szl.lambda/v1"
WEIGHT_SUM_TOL = 1e-12
TIE_EPS = 1e-9

GO, NO_GO, ABSTAIN, BLOCK = "GO", "NO_GO", "ABSTAIN", "BLOCK"
ZERO_VETO, BELOW_TAU, NUMERIC_TIE = "ZERO_VETO", "BELOW_TAU", "NUMERIC_TIE"


class LambdaV1Error(ValueError):
    def __init__(self, code: str, detail: str = "") -> None:
        self.code = code
        self.detail = detail
        super().__init__(f"{code}: {detail}" if detail else code)


def _is_real(v: Any) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def validate_lambda_inputs(axes: Any, weights: Any) -> None:
    for name, seq in (("axes", axes), ("weights", weights)):
        if not isinstance(seq, (list, tuple)):
            raise LambdaV1Error("LAMBDA_TYPE_INVALID", f"{name} must be a list or tuple")
    if len(axes) == 0 or len(weights) == 0:
        raise LambdaV1Error("LAMBDA_EMPTY", f"len(axes)={len(axes)}")
    if len(axes) != len(weights):
        raise LambdaV1Error("LAMBDA_LENGTH_MISMATCH", f"{len(axes)} != {len(weights)}")
    for name, seq in (("axes", axes), ("weights", weights)):
        for i, v in enumerate(seq):
            if not _is_real(v):
                raise LambdaV1Error("LAMBDA_TYPE_INVALID", f"{name}[{i}] not a real number")
    for i, x in enumerate(axes):
        if isinstance(x, float) and not math.isfinite(x):
            raise LambdaV1Error("LAMBDA_NONFINITE_AXIS", f"axes[{i}]={x!r}")
    for i, x in enumerate(axes):
        if not 0 <= x <= 1:
            raise LambdaV1Error("LAMBDA_AXIS_OUT_OF_RANGE", f"axes[{i}]={x!r} outside [0,1]")
    for i, w in enumerate(weights):
        if isinstance(w, float) and not math.isfinite(w):
            raise LambdaV1Error("LAMBDA_NONFINITE_WEIGHT", f"weights[{i}]={w!r}")
    for i, w in enumerate(weights):
        if w <= 0:
            raise LambdaV1Error("LAMBDA_WEIGHT_NONPOSITIVE", f"weights[{i}]={w!r}")
    # Finite positive weights can still overflow fsum (e.g. [1e308, 1e308]); the
    # reference maps that to LAMBDA_WEIGHT_SUM so the gate BLOCKs, never raises.
    try:
        total = math.fsum(weights)
    except OverflowError:
        raise LambdaV1Error("LAMBDA_WEIGHT_SUM", "sum of weights overflows a float") from None
    if not abs(total - 1.0) <= WEIGHT_SUM_TOL:
        raise LambdaV1Error("LAMBDA_WEIGHT_SUM", f"fsum(weights)={total!r}")


def lambda_v1(axes: list[float], weights: list[float]) -> float:
    """Weighted geometric mean per szl.lambda/v1. Zero axis => 0 (veto)."""
    validate_lambda_inputs(axes, weights)
    if any(x == 0 for x in axes):
        return 0.0
    return math.exp(math.fsum(w * math.log(x) for x, w in zip(axes, weights)))


def gate_v1(axes: list[float], weights: list[float], tau: float) -> tuple[str, str | None]:
    """szl.lambda/v1 gate. Returns (verdict, code). tau is required."""
    if not _is_real(tau) or (isinstance(tau, float) and not math.isfinite(tau)) \
            or not (0 < tau <= 1):
        return BLOCK, "LAMBDA_TAU_INVALID"
    try:
        lam = lambda_v1(axes, weights)
    except LambdaV1Error as e:
        return BLOCK, e.code
    if lam == 0.0:
        return NO_GO, ZERO_VETO
    if abs(math.log(lam) - math.log(tau)) <= TIE_EPS:
        return ABSTAIN, NUMERIC_TIE
    return (GO, None) if math.log(lam) > math.log(tau) else (NO_GO, BELOW_TAU)


# --- Acoustic formula axis scores -------------------------------------------
# Each axis is an honest, documented mapping from a measured check result to
# [0,1]. Binary checks map pass->1.0, fail->0.0 (a fail is a veto). Graded
# checks use an explicit, frozen threshold function. No hidden scalarization.

def acoustic_axes() -> tuple[list[float], list[str]]:
    """Run all acoustic checks and return (scores, names)."""
    sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
    import burgers
    import duffing
    import aero_ff

    scores, names = [], []

    def binary(ok: bool, name: str) -> None:
        scores.append(1.0 if ok else 0.0)
        names.append(name)

    b = burgers.solve_burgers(nx=200, t_final=0.1, cfl=0.5)
    binary(b.analytic_residual < 0.05, "burgers_analytic")
    binary(burgers.refinement_check(), "burgers_refinement")

    traj = duffing.simulate_duffing()
    coeffs = duffing.identify_coefficients(traj)
    binary(abs(coeffs["alpha"] + 1.0) < 0.1 and abs(coeffs["beta"] - 1.0) < 0.1
           and abs(coeffs["delta"] - 0.3) < 0.1, "duffing_recovery")
    binary(duffing.ghost_identifiability(traj), "duffing_ghost_rejection")

    binary(aero_ff.murnaghan_vL(100e9, 80e9, -100e9, -50e9, -30e9, 0.0, 7800) > 0,
           "aero_murnaghan_physical")
    # Three distinct nonzero stress states would expose rank 2 or 3 if A, B, C
    # entered vL separately, so this axis can fail: it checks the documented
    # claim (rank 1, only S = 2A + 10B + 6C identifiable) against the formula.
    ident = aero_ff.murnaghan_identifiability(100e9, 80e9, [0.0, 1e6, 5e6, 2e7], 7800)
    binary(ident["sensitivity_rank"] == 1
           and ident["identifiable_combination"] == "S = 2A + 10B + 6C"
           and ident["individual_ABC_identifiable"] is False, "aero_identifiability_honest")
    ep = aero_ff.exceptional_point_2x2(1.0, 0.5, 0.5)
    binary(ep["is_true_EP"], "aero_ep_defectiveness")
    binary(aero_ff.chern_number_dirac(-1.0, nk=24) == -1, "aero_chern_quantization")

    return scores, names


def governed_validation(tau: float,
                        weights: list[float] | None = None) -> dict:
    """Full governed acoustic validation receipt."""
    scores, names = acoustic_axes()
    if weights is None:
        n = len(scores)
        weights = [1.0 / n] * n
    verdict, code = gate_v1(scores, weights, tau)
    try:
        lam = lambda_v1(scores, weights)
    except LambdaV1Error:
        # gate_v1 already returned BLOCK; preserve the fail-closed receipt
        # instead of raising before the caller sees it. Keep the gate's code:
        # tau is checked first, so with an invalid tau AND invalid weights the
        # contract's code is LAMBDA_TAU_INVALID, not the weights' error.
        if verdict != BLOCK:
            raise
        lam = 0.0
    receipt = {
        "schema": "szl.acoustic.governed/v1",
        "lambda_schema": SCHEMA,
        "axes": dict(zip(names, scores)),
        "weights": weights,
        "lambda": round(lam, 12),
        "tau": tau,
        "verdict": verdict,
        "code": code,
        "lambda_status": "ADVISORY",
        "conjecture_status": "Conjecture 1 (machine-checked false as stated); not a theorem",
        "claims": {
            "device_physical": "HYPOTHESIS_UNVERIFIED",
            "ns_regularity": "NOT_ESTABLISHED",
        },
        "truth_label": "MEASURED_SOFTWARE_SIMULATION",
    }
    return receipt


def selfcheck() -> bool:
    ok = True
    # contract behaviors
    ok &= lambda_v1([1.0, 1.0], [0.5, 0.5]) == 1.0
    ok &= lambda_v1([0.0, 0.5], [0.5, 0.5]) == 0.0
    ok &= abs(lambda_v1([0.25, 1.0], [0.5, 0.5]) - 0.5) < 1e-12
    ok &= gate_v1([0.0, 1.0], [0.5, 0.5], 0.5) == (NO_GO, ZERO_VETO)
    ok &= gate_v1([1.0, 1.0], [0.5, 0.5], 0.5)[0] == GO
    ok &= gate_v1([0.2, 0.2], [0.5, 0.5], 0.5)[0] == NO_GO
    ok &= gate_v1([1.0], [1.0], 1.5) == (BLOCK, "LAMBDA_TAU_INVALID")
    ok &= gate_v1([1.0], [1.0], 0.0) == (BLOCK, "LAMBDA_TAU_INVALID")
    try:
        lambda_v1([1.0], [1.5])  # weights must sum to 1
        ok = False
    except LambdaV1Error as e:
        ok &= e.code == "LAMBDA_WEIGHT_SUM"
    # full pipeline: all checks pass -> GO at a modest tau
    r = governed_validation(tau=0.5)
    ok &= r["verdict"] == GO
    return bool(ok)


def main() -> int:
    ok = selfcheck()
    receipt = governed_validation(tau=0.5)
    receipt["selfcheck_pass"] = ok
    print(json.dumps(receipt, indent=2))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
