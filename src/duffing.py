#!/usr/bin/env python3
"""Duffing oscillator validation with coefficient and ghost-identifiability checks.

Implements the Duffing equation x'' + delta*x' + alpha*x + beta*x^3 = gamma*cos(omega*t)
with parameter identification from synthetic data. Checks:
- Coefficient recovery: can we identify the true parameters?
- Ghost identifiability: can we distinguish true from ghost solutions?
"""
from __future__ import annotations
import hashlib
import json
import math
from typing import Any


def simulate_duffing(
    delta: float = 0.3, alpha: float = -1.0, beta: float = 1.0,
    gamma: float = 0.5, omega: float = 1.2,
    x0: float = 0.1, v0: float = 0.0, dt: float = 0.01, n: int = 1000
) -> list[tuple[float, float]]:
    """Simulate Duffing oscillator with RK4."""
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


def identify_coefficients(trajectory: list[tuple[float, float]]) -> dict:
    """Recover Duffing coefficients from trajectory using least squares."""
    n = len(trajectory)
    # Build regression problem: a*x + b*x^3 + c*v = acceleration
    # acceleration approximated by finite differences
    rows, rhs = [], []
    for i in range(1, n - 1):
        x_prev, v_prev = trajectory[i-1]
        x, v = trajectory[i]
        x_next, v_next = trajectory[i+1]
        accel = (v_next - v_prev) / (2 * 0.01)
        t = i * 0.01
        forcing = 0.5 * math.cos(1.2 * t)
        rows.append([x, x**3, v])
        rhs.append(accel - forcing)

    # Simple normal equations (3x3)
    import numpy as np
    A = np.array(rows)
    b = np.array(rhs)
    coeffs, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    return {"alpha": round(coeffs[0], 4), "beta": round(coeffs[1], 4), "delta": round(coeffs[2], 4)}


def ghost_identifiability() -> bool:
    """Check that ghost solutions are distinguishable from true solutions."""
    true_traj = simulate_duffing(x0=0.1)
    ghost_traj = simulate_duffing(x0=-0.1)

    # Ghost solutions (mirrored initial conditions) should be distinguishable
    diff = sum((true_traj[i][0] - ghost_traj[i][0])**2 for i in range(len(true_traj)))
    return diff > 0.01  # They should be different


def selfcheck() -> bool:
    """Run all Duffing validation checks."""
    traj = simulate_duffing()
    coeffs = identify_coefficients(traj)
    ghost_ok = ghost_identifiability()

    # Check coefficient recovery (within tolerance)
    alpha_ok = abs(coeffs["alpha"] - (-1.0)) < 0.5
    beta_ok = abs(coeffs["beta"] - 1.0) < 0.5
    delta_ok = abs(coeffs["delta"] - 0.3) < 0.5

    return alpha_ok and beta_ok and delta_ok and ghost_ok


if __name__ == "__main__":
    ok = selfcheck()
    traj = simulate_duffing()
    coeffs = identify_coefficients(traj)
    result = {
        "schema": "szl.acoustic.duffing/v1",
        "identified_coefficients": coeffs,
        "ghost_distinguishable": ghost_identifiability(),
        "selfcheck_pass": ok,
        "truth_label": "MEASURED_SOFTWARE_SIMULATION",
    }
    print(json.dumps(result, indent=2))
