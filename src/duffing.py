#!/usr/bin/env python3
"""Duffing oscillator validation with coefficient and ghost-identifiability checks.

Implements the Duffing equation x'' + delta*x' + alpha*x + beta*x^3 = gamma*cos(omega*t)
with parameter identification from synthetic data. Checks:

- Coefficient recovery: identify (alpha, beta, delta) from a trajectory via
  least squares on x'' = -delta*x' - alpha*x - beta*x^3 + gamma*cos(omega*t).
  Signs are handled explicitly: the regression recovers the *negated*
  coefficients, which are then flipped back before reporting.
- Ghost identifiability: a ghost model (sign-flipped cubic stiffness) must be
  *rejected* by the data — its one-step prediction residual must be far worse
  than the recovered true model's. Mirror-initial-condition trajectories are
  not used as evidence.

This is a standalone implementation referencing SZL source equations.
It does not validate physical claims about any device.
"""
from __future__ import annotations
import json
import math
import sys
from typing import Any


def simulate_duffing(
    delta: float = 0.3, alpha: float = -1.0, beta: float = 1.0,
    gamma: float = 0.5, omega: float = 1.2,
    x0: float = 0.1, v0: float = 0.0, dt: float = 0.01, n: int = 1000
) -> list[tuple[float, float]]:
    """Simulate Duffing oscillator with RK4. Returns [(x, v), ...] of length n+1."""
    def deriv(x, v, t):
        return v, -delta * v - alpha * x - beta * x**3 + gamma * math.cos(omega * t)

    x, v = x0, v0
    trajectory = [(x, v)]
    for i in range(n):
        t = i * dt
        k1x, k1v = deriv(x, v, t)
        k2x, k2v = deriv(x + 0.5*dt*k1x, v + 0.5*dt*k1v, t + 0.5*dt)
        k3x, k3v = deriv(x + 0.5*dt*k2x, v + 0.5*dt*k2v, t + 0.5*dt)
        k4x, k4v = deriv(x + dt*k3x, v + dt*k3v, t + dt)
        x += dt/6 * (k1x + 2*k2x + 2*k3x + k4x)
        v += dt/6 * (k1v + 2*k2v + 2*k3v + k4v)
        trajectory.append((x, v))
    return trajectory


def identify_coefficients(trajectory: list[tuple[float, float]],
                          dt: float = 0.01, gamma: float = 0.5,
                          omega: float = 1.2) -> dict:
    """Recover Duffing coefficients from a trajectory using least squares.

    Model: x'' = -alpha*x - beta*x^3 - delta*v + gamma*cos(omega*t)
    Regression on [x, x^3, v] against (accel - forcing) therefore yields
    (-alpha, -beta, -delta); the signs are flipped back here.
    """
    rows, rhs = [], []
    for i in range(1, len(trajectory) - 1):
        x_prev, v_prev = trajectory[i - 1]
        x, v = trajectory[i]
        x_next, v_next = trajectory[i + 1]
        accel = (v_next - v_prev) / (2 * dt)
        t = i * dt
        forcing = gamma * math.cos(omega * t)
        rows.append([x, x**3, v])
        rhs.append(accel - forcing)

    import numpy as np
    A = np.array(rows)
    b = np.array(rhs)
    neg_coeffs, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    return {
        "alpha": round(-float(neg_coeffs[0]), 4),
        "beta": round(-float(neg_coeffs[1]), 4),
        "delta": round(-float(neg_coeffs[2]), 4),
    }


def prediction_residual(trajectory: list[tuple[float, float]],
                        coeffs: dict, dt: float = 0.01,
                        gamma: float = 0.5, omega: float = 1.2) -> float:
    """RMS one-step-acceleration residual of a candidate model against data."""
    sq = 0.0
    count = 0
    for i in range(1, len(trajectory) - 1):
        x_prev, v_prev = trajectory[i - 1]
        x, v = trajectory[i]
        x_next, v_next = trajectory[i + 1]
        accel = (v_next - v_prev) / (2 * dt)
        t = i * dt
        predicted = (-coeffs["alpha"] * x - coeffs["beta"] * x**3
                     - coeffs["delta"] * v + gamma * math.cos(omega * t))
        sq += (accel - predicted) ** 2
        count += 1
    return math.sqrt(sq / max(count, 1))


def ghost_identifiability(trajectory: list[tuple[float, float]] | None = None,
                          ratio_threshold: float = 10.0) -> bool:
    """Reject ghost models: sign-flipped candidates must fit far worse.

    A ghost model (wrong-sign cubic stiffness, or wrong-sign linear
    stiffness) must have a prediction residual at least `ratio_threshold`
    times worse than the recovered model. This tests whether the data can
    distinguish the true dynamics, not merely whether two simulations
    started from different initial conditions differ.
    """
    if trajectory is None:
        trajectory = simulate_duffing()
    recovered = identify_coefficients(trajectory)

    true_residual = prediction_residual(trajectory, recovered)
    ghosts = [
        {**recovered, "beta": -recovered["beta"]},
        {**recovered, "alpha": -recovered["alpha"]},
    ]
    worst_allowed = max(true_residual * ratio_threshold, 1e-9)
    return all(prediction_residual(trajectory, g) > worst_allowed for g in ghosts)


def selfcheck() -> bool:
    """Run all Duffing validation checks."""
    traj = simulate_duffing()
    coeffs = identify_coefficients(traj)
    ghost_ok = ghost_identifiability(traj)
    alpha_ok = abs(coeffs["alpha"] - (-1.0)) < 0.1
    beta_ok = abs(coeffs["beta"] - 1.0) < 0.1
    delta_ok = abs(coeffs["delta"] - 0.3) < 0.1
    return alpha_ok and beta_ok and delta_ok and ghost_ok


def main() -> int:
    ok = selfcheck()
    traj = simulate_duffing()
    coeffs = identify_coefficients(traj)
    result = {
        "schema": "szl.acoustic.duffing/v2",
        "identified_coefficients": coeffs,
        "ghost_distinguishable": ghost_identifiability(traj),
        "selfcheck_pass": ok,
        "truth_label": "MEASURED_SOFTWARE_SIMULATION",
    }
    print(json.dumps(result, indent=2))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
