#!/usr/bin/env python3
"""Burgers equation solver with time-evolved analytic comparison and grid refinement.

Implements the inviscid Burgers equation u_t + u*u_x = 0 with a sin(2*pi*x)
initial condition using the Lax-Friedrichs scheme, with:

- Analytic comparison: for t < t_break = 1/(2*pi) the exact solution satisfies
  the implicit characteristic relation u = sin(2*pi*(x - u*t)), solved per grid
  point by fixed-point iteration. The numerical solution is compared against
  this time-evolved analytic solution, not against the initial condition.
- Grid refinement: at fixed CFL, the L2 error against the analytic solution
  must decrease monotonically as the mesh is refined (first-order method).
- Receipt: results are hashed; sealing into IMMUNE is a separate, explicit step.

This is a standalone implementation referencing SZL source equations.
It does not validate physical claims about any device.
"""
from __future__ import annotations
import hashlib
import json
import math
import sys
from dataclasses import dataclass, field
from typing import Any

T_BREAK = 1.0 / (2.0 * math.pi)


@dataclass
class BurgersResult:
    """Result of a Burgers equation solve with verification metadata."""
    nx: int
    nt: int
    dx: float
    dt: float
    CFL: float
    t_final: float
    solution: list[float] = field(default_factory=list)
    analytic_residual: float = 0.0
    refinement_converged: bool = False
    receipt_hash: str = ""

    def to_receipt(self) -> dict[str, Any]:
        return {
            "schema": "szl.acoustic.burgers/v2",
            "nx": self.nx, "nt": self.nt,
            "dx": self.dx, "dt": self.dt, "CFL": self.CFL, "t_final": self.t_final,
            "analytic_residual": round(self.analytic_residual, 12),
            "refinement_converged": self.refinement_converged,
            "solution_sha256": hashlib.sha256(
                json.dumps(self.solution).encode()).hexdigest(),
            "truth_label": "MEASURED_SOFTWARE_SIMULATION",
        }


def analytic_solution(x: float, t: float, iterations: int = 50) -> float:
    """Exact inviscid Burgers solution for u0 = sin(2*pi*x), valid t < t_break.

    Each characteristic gives the implicit relation u = sin(2*pi*(x - u*t)),
    solved by fixed-point iteration, which converges for t < t_break.
    """
    if t >= T_BREAK:
        raise ValueError(
            f"analytic_solution valid only for t < t_break = {T_BREAK:.6f}; got t={t}")
    u = math.sin(2.0 * math.pi * x)
    for _ in range(iterations):
        u = math.sin(2.0 * math.pi * (x - u * t))
    return u


def solve_burgers(nx: int = 100, t_final: float = 0.1, cfl: float = 0.5) -> BurgersResult:
    """Solve inviscid Burgers with Lax-Friedrichs at fixed CFL.

    Requires t_final < t_break so the analytic comparison is well-defined.
    """
    if t_final >= T_BREAK:
        raise ValueError(f"t_final must be < t_break = {T_BREAK:.6f}")
    if not (0 < cfl <= 1.0):
        raise ValueError("cfl must be in (0, 1]")

    dx = 1.0 / nx
    dt = cfl * dx
    nt = max(1, int(math.ceil(t_final / dt)))
    dt_actual = t_final / nt
    cfl_actual = dt_actual / dx

    u = [math.sin(2.0 * math.pi * i * dx) for i in range(nx)]

    for _ in range(nt):
        u_new = [0.0] * nx
        for i in range(nx):
            left = u[(i - 1) % nx]
            right = u[(i + 1) % nx]
            u_new[i] = 0.5 * (left + right) - 0.5 * cfl_actual * (
                0.5 * right**2 - 0.5 * left**2)
        u = u_new

    # Compare against the time-evolved analytic solution at t_final.
    residual_sq = 0.0
    for i in range(nx):
        exact = analytic_solution(i * dx, t_final)
        residual_sq += (u[i] - exact) ** 2
    residual = math.sqrt(residual_sq / nx)

    return BurgersResult(
        nx=nx, nt=nt, dx=dx, dt=dt_actual, CFL=cfl_actual, t_final=t_final,
        solution=[round(v, 8) for v in u],
        analytic_residual=residual,
        refinement_converged=False,
    )


def refinement_errors(levels: list[int] = (50, 100, 200, 400),
                      t_final: float = 0.1, cfl: float = 0.5) -> list[float]:
    """L2 errors against the analytic solution at fixed t and CFL."""
    return [solve_burgers(nx=nx, t_final=t_final, cfl=cfl).analytic_residual
            for nx in levels]


def refinement_check(levels: list[int] = (50, 100, 200, 400),
                     t_final: float = 0.1, cfl: float = 0.5) -> bool:
    """Verify strict error reduction under mesh refinement at fixed CFL.

    Lax-Friedrichs is first-order, so each doubling of the mesh must reduce
    the analytic error by a clear factor (>= 1.2x here, conservative bound).
    """
    errors = refinement_errors(list(levels), t_final, cfl)
    return all(errors[i + 1] < errors[i] for i in range(len(errors) - 1)) and \
        all(errors[i] / errors[i + 1] >= 1.2 for i in range(len(errors) - 1))


def selfcheck() -> bool:
    """Run all Burgers validation checks."""
    result = solve_burgers()
    ref = refinement_check()
    result.refinement_converged = ref
    # A 100-cell first-order solution at t=0.1 should be well within this bound.
    return result.analytic_residual < 0.05 and ref


def main() -> int:
    ok = selfcheck()
    result = solve_burgers()
    result.refinement_converged = refinement_check()
    receipt = result.to_receipt()
    receipt["selfcheck_pass"] = ok
    print(json.dumps(receipt, indent=2))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
