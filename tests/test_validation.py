import math
import subprocess
import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from burgers import T_BREAK, analytic_solution, refinement_check, refinement_errors, solve_burgers  # noqa: E402
from duffing import ghost_identifiability, identify_coefficients, prediction_residual, simulate_duffing  # noqa: E402


class TestBurgers:
    def test_analytic_solution_matches_initial_condition_at_t0(self):
        for x in [0.0, 0.1, 0.25, 0.5, 0.75]:
            assert abs(analytic_solution(x, 0.0) - math.sin(2 * math.pi * x)) < 1e-9

    def test_analytic_solution_is_time_evolved(self):
        # At t>0 the analytic solution must differ from the initial condition.
        t = 0.1
        diffs = [abs(analytic_solution(x, t) - math.sin(2 * math.pi * x))
                 for x in [0.05, 0.1, 0.15, 0.2]]
        assert max(diffs) > 1e-3

    def test_analytic_rejects_t_past_break(self):
        with pytest.raises(ValueError):
            analytic_solution(0.5, T_BREAK + 0.01)

    def test_numerical_matches_time_evolved_analytic(self):
        result = solve_burgers(nx=400, t_final=0.1, cfl=0.5)
        assert result.analytic_residual < 0.02

    def test_residual_is_not_against_initial_condition(self):
        # A converged solution must differ from the IC: proves the comparison
        # target actually moved in time.
        result = solve_burgers(nx=400, t_final=0.1, cfl=0.5)
        ic_residual = math.sqrt(sum(
            (u - math.sin(2 * math.pi * i * result.dx)) ** 2 for i, u in enumerate(result.solution)
        ) / result.nx)
        assert ic_residual > 0.1
        assert result.analytic_residual < ic_residual

    def test_refinement_error_decreases(self):
        errors = refinement_errors([50, 100, 200, 400], t_final=0.1, cfl=0.5)
        assert all(errors[i + 1] < errors[i] for i in range(len(errors) - 1))
        assert refinement_check()

    def test_invalid_inputs_rejected(self):
        with pytest.raises(ValueError):
            solve_burgers(nx=100, t_final=1.0)
        with pytest.raises(ValueError):
            solve_burgers(nx=100, t_final=0.1, cfl=1.5)

    def test_selfcheck_exit_code_zero(self):
        r = subprocess.run([sys.executable, str(SRC / "burgers.py")], capture_output=True)
        assert r.returncode == 0, r.stdout + r.stderr


class TestDuffing:
    def test_coefficient_recovery_signs(self):
        traj = simulate_duffing()
        coeffs = identify_coefficients(traj)
        assert abs(coeffs["alpha"] - (-1.0)) < 0.1
        assert abs(coeffs["beta"] - 1.0) < 0.1
        assert abs(coeffs["delta"] - 0.3) < 0.1

    def test_recovery_detects_changed_parameters(self):
        # Recovery must track the true parameters, not just return defaults.
        traj = simulate_duffing(delta=0.6, alpha=-2.0, beta=1.5)
        coeffs = identify_coefficients(traj)
        assert abs(coeffs["alpha"] - (-2.0)) < 0.1
        assert abs(coeffs["beta"] - 1.5) < 0.1
        assert abs(coeffs["delta"] - 0.6) < 0.1

    def test_ghost_models_rejected(self):
        traj = simulate_duffing()
        assert ghost_identifiability(traj)

    def test_ghost_rejection_is_not_about_initial_conditions(self):
        # Even with identical initial conditions, a ghost model must fit worse.
        traj = simulate_duffing(x0=0.1)
        recovered = identify_coefficients(traj)
        true_r = prediction_residual(traj, recovered)
        ghost = {**recovered, "beta": -recovered["beta"]}
        assert prediction_residual(traj, ghost) > 10 * true_r

    def test_selfcheck_exit_code_zero(self):
        r = subprocess.run([sys.executable, str(SRC / "duffing.py")], capture_output=True)
        assert r.returncode == 0, r.stdout + r.stderr


class TestFailureInjection:
    def test_burgers_detects_broken_scheme(self):
        # A scheme with the wrong flux sign must fail the analytic check.
        import burgers as B
        orig = B.solve_burgers

        def broken(nx=100, t_final=0.1, cfl=0.5):
            dx = 1.0 / nx
            dt = cfl * dx
            nt = max(1, int(math.ceil(t_final / dt)))
            u = [math.sin(2 * math.pi * i * dx) for i in range(nx)]
            for _ in range(nt):
                u_new = [0.0] * nx
                for i in range(nx):
                    left, right = u[(i - 1) % nx], u[(i + 1) % nx]
                    u_new[i] = 0.5 * (left + right) + 0.5 * cfl * (
                        0.5 * right**2 - 0.5 * left**2)
                u = u_new
            r = orig(nx=nx, t_final=t_final, cfl=cfl)
            r.solution = [round(v, 8) for v in u]
            r.analytic_residual = math.sqrt(sum(
                (u[i] - B.analytic_solution(i * dx, t_final)) ** 2 for i in range(nx)) / nx)
            return r

        B.solve_burgers = broken
        try:
            assert broken().analytic_residual > 0.05
        finally:
            B.solve_burgers = orig

    def test_duffing_detects_wrong_sign_recovery(self):
        # Simulating the pre-fix behavior: reporting negated coefficients
        # must fail the recovery check.
        traj = simulate_duffing()
        coeffs = identify_coefficients(traj)
        wrong = {k: -v for k, v in coeffs.items()}
        assert abs(wrong["alpha"] - (-1.0)) > 0.1
