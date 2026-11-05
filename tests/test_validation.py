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


class TestAeroFF:
    def test_murnaghan_velocity_responds_to_stress(self):
        import sys; sys.path.insert(0, str(SRC))
        from aero_ff import murnaghan_vL
        v0 = murnaghan_vL(100e9, 80e9, -100e9, -50e9, -30e9, 0.0, 7800)
        v1 = murnaghan_vL(100e9, 80e9, -100e9, -50e9, -30e9, 50e6, 7800)
        assert v0 > 0 and v1 > 0 and v1 != v0

    def test_murnaghan_ABC_not_individually_identifiable(self):
        import sys; sys.path.insert(0, str(SRC))
        from aero_ff import murnaghan_identifiability
        r = murnaghan_identifiability(100e9, 80e9, [0.0, 1e6, 5e6], 7800)
        assert r["sensitivity_rank"] == 1
        assert r["individual_ABC_identifiable"] is False

    def test_gain_db_ambiguity_surfaced(self):
        import sys, math; sys.path.insert(0, str(SRC))
        from aero_ff import ep_gain_db
        # G = 2/(1-0.25) = 8/3: amplitude dB = 20log10(G), power dB = 10log10(G)
        assert abs(ep_gain_db(0.5) - 20 * math.log10(8/3)) < 1e-9
        assert abs(ep_gain_db(0.5, power_ratio=True) - 10 * math.log10(8/3)) < 1e-9

    def test_true_EP_requires_defective_eigenvectors(self):
        import sys; sys.path.insert(0, str(SRC))
        from aero_ff import exceptional_point_2x2
        ep = exceptional_point_2x2(1.0, 0.5, 0.5)
        assert ep["is_true_EP"] and ep["defective"] and ep["eigenvectors_parallel"]
        off = exceptional_point_2x2(1.0, 0.5, 0.25)
        assert not off["is_true_EP"] and not off["coalesced"]

    def test_rankine_hugoniot_scaling(self):
        import sys; sys.path.insert(0, str(SRC))
        from aero_ff import rankine_hugoniot
        assert rankine_hugoniot(340.0, 1.4, 50.0) == 400.0

    def test_chern_number_quantization(self):
        import sys; sys.path.insert(0, str(SRC))
        from aero_ff import chern_number_dirac
        assert chern_number_dirac(-1.0, nk=24) == -1
        assert chern_number_dirac(1.0, nk=24) == 1
        assert chern_number_dirac(3.0, nk=24) == 0

    def test_selfcheck_exit_code_zero(self):
        r = subprocess.run([sys.executable, str(SRC / "aero_ff.py")], capture_output=True)
        assert r.returncode == 0, r.stdout + r.stderr

    def test_claims_status_labels_honest(self):
        import sys, json; sys.path.insert(0, str(SRC))
        r = subprocess.run([sys.executable, str(SRC / "aero_ff.py")], capture_output=True)
        d = json.loads(r.stdout)
        assert d["claims_status"]["navier_stokes_regularity_implication"] == "NOT_ESTABLISHED"
        assert d["claims_status"]["host_ai_accuracy_percentages"] == "UNVERIFIED_NO_PROTOCOL"
