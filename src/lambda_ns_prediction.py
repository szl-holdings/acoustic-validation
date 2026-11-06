"""Lambda-scored NS regularization: the Schur-concavity prediction.

Connects the estate's Lean-proved Schur-concavity (the Robin Hood lemma:
equalizing transfers weakly increase Lambda) to the Aero-FF Navier-Stokes
regularization hypothesis, producing a falsifiable prediction:

    Systems with balanced acoustoelastic stress across all Murnaghan
    constants should be more resistant to velocity blowup than systems
    with the same total stress concentrated in one constant.

The Lambda gate's Schur-concavity is PROVED in Lean (Lutar/Lambda/SchurConcave.lean,
0 sorry). The application to NS regularization is a HYPOTHESIS derived from
the proved property — it is not itself proved and does not validate the
Aero-FF physical claims.

SPDX-License-Identifier: Apache-2.0
Source math: szl-holdings/lutar-lean Lutar/Lambda/SchurConcave.lean (PROVEN)
"""
from __future__ import annotations
import json
import math
import sys
from typing import Any

SCHEMA = "szl.acoustic.lambda-scored-ns/v1"


def lambda_aggregate(axes: list[float], weights: list[float] | None = None) -> float:
    """Lambda_w(x) = prod x_k^w_k (szl.lambda/v1 contract, stdlib)."""
    if not axes or any(x < 0 or x > 1 for x in axes):
        raise ValueError("axes must be in [0,1]")
    if any(x == 0 for x in axes):
        return 0.0
    if weights is None:
        weights = [1.0 / len(axes)] * len(axes)
    if abs(math.fsum(weights) - 1.0) > 1e-12:
        raise ValueError("weights must sum to 1")
    return math.exp(math.fsum(w * math.log(x) for x, w in zip(axes, weights)))


def schur_concavity_prediction(stress_distribution: list[float]) -> dict:
    """The Robin Hood lemma applied to stress: balanced > concentrated.

    Given a stress distribution across the Murnaghan-active axes (A, B, C
    contributions), compute Lambda and compare with the equalized distribution.

    Schur-concavity (PROVED in Lean): if x is more balanced than y (x majorized
    by y), then Lambda(x) >= Lambda(y). Equalized stress scores >= concentrated.

    This produces a falsifiable prediction for the NS regularization:
    the regularization mechanism should be MORE effective when the
    acoustoelastic stress is distributed evenly.
    """
    if not stress_distribution or any(s < 0 for s in stress_distribution):
        raise ValueError("stress distribution must be non-negative")
    total = sum(stress_distribution)
    if total == 0:
        return {"prediction": "trivial", "lambda": 1.0,
                "note": "zero stress: no regularization needed"}
    # normalize to [0,1] axes
    max_stress = max(stress_distribution)
    if max_stress == 0:
        axes = [1.0] * len(stress_distribution)
    else:
        axes = [s / max_stress for s in stress_distribution]
    # the equalized distribution (Robin Hood: transfer from rich to poor)
    n = len(axes)
    equalized = [sum(axes) / n] * n
    lam_actual = lambda_aggregate(axes)
    lam_equalized = lambda_aggregate(equalized)
    # Schur-concavity: lam_equalized >= lam_actual (proved)
    schur_gap = lam_equalized - lam_actual
    return {
        "schema": SCHEMA,
        "stress_distribution": stress_distribution,
        "normalized_axes": [round(a, 6) for a in axes],
        "lambda_actual": round(lam_actual, 9),
        "lambda_equalized": round(lam_equalized, 9),
        "schur_gap": round(schur_gap, 9),
        "schur_concavity_holds": lam_equalized >= lam_actual - 1e-12,
        "provenance": "Lean: Lutar/Lambda/SchurConcave.lean (PROVED, 0 sorry)",
        "prediction": (
            f"The NS regularization should be more effective when stress is "
            f"balanced across axes (Lambda_equalized={lam_equalized:.4f}) than "
            f"when concentrated (Lambda_actual={lam_actual:.4f}). "
            f"Gap: {schur_gap:.4f}. This is the Schur-concavity prediction: "
            f"equalizing transfers weakly increase Lambda."
        ),
        "hypothesis_status": "HYPOTHESIS_FROM_PROVED_PROPERTY",
        "note": (
            "The Schur-concavity of Lambda is PROVED in Lean. The application "
            "to NS regularization is a HYPOTHESIS derived from the proved "
            "property. It does not validate the Aero-FF physical claims."
        ),
        "truth_label": "MEASURED_SOFTWARE_SIMULATION",
    }


def murnaghan_balance_index(A: float, B: float, C: float) -> dict:
    """Compute the balance of the third-order Murnaghan constants.

    The Aero-FF stress kick depends on the combination 2A + 10B + 6C.
    The Schur-concavity prediction: a balanced (A, B, C) contribution should
    produce a more stable regularization than an imbalanced one.
    """
    contributions = [abs(2 * A), abs(10 * B), abs(6 * C)]
    total = sum(contributions)
    if total == 0:
        return {"balance_index": 1.0, "note": "zero third-order contributions"}
    # normalized balance: 1.0 = perfectly balanced, 0.0 = maximally concentrated
    n = len(contributions)
    equal = total / n
    variance = sum((c - equal) ** 2 for c in contributions) / n
    max_variance = equal ** 2 * (n - 1)  # when one axis has everything
    balance = 1.0 - (variance / max_variance if max_variance > 0 else 0.0)
    return {
        "contributions": {"2A": contributions[0], "10B": contributions[1], "6C": contributions[2]},
        "balance_index": round(balance, 6),
        "interpretation": (
            f"Balance index {balance:.4f}: "
            f"{'balanced' if balance > 0.7 else 'moderately balanced' if balance > 0.4 else 'concentrated'} "
            f"third-order contributions. Per Schur-concavity, higher balance "
            f"predicts more stable regularization."
        ),
        "schur_prediction": schur_concavity_prediction(contributions),
    }


def selfcheck() -> bool:
    # balanced > concentrated (Schur-concavity)
    balanced = schur_concavity_prediction([1.0, 1.0, 1.0])
    concentrated = schur_concavity_prediction([3.0, 0.0, 0.0])
    ok = balanced["lambda_equalized"] >= concentrated["lambda_actual"]
    # the proved property: equalized >= actual for any distribution
    test_cases = [[1, 2, 3], [5, 1, 1], [0.1, 0.5, 0.9]]
    for tc in test_cases:
        r = schur_concavity_prediction([float(x) for x in tc])
        ok &= r["schur_concavity_holds"]
    # Murnaghan balance
    m = murnaghan_balance_index(-100e9, -50e9, -30e9)
    ok &= 0 <= m["balance_index"] <= 1
    return bool(ok)


if __name__ == "__main__":
    ok = selfcheck()
    result = {
        "schema": SCHEMA,
        "selfcheck_pass": ok,
        "provenance": "Lean: Lutar/Lambda/SchurConcave.lean (PROVED, 0 sorry)",
        "hypothesis": (
            "Schur-concavity of Lambda (PROVED) predicts that balanced "
            "acoustoelastic stress is more resistant to NS velocity blowup "
            "than concentrated stress. This is a falsifiable prediction "
            "derived from the estate's governance math applied to the "
            "Aero-FF hypothesis."
        ),
        "hypothesis_status": "HYPOTHESIS_FROM_PROVED_PROPERTY",
    }
    print(json.dumps(result, indent=1))
    sys.exit(0 if ok else 1)
