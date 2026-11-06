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


class TestFormulaGovernance:
    def test_lambda_v1_contract(self):
        import sys; sys.path.insert(0, str(SRC))
        from formula_governance import lambda_v1, gate_v1, LambdaV1Error, NO_GO, ZERO_VETO, GO, BLOCK
        assert lambda_v1([1.0, 1.0], [0.5, 0.5]) == 1.0
        assert lambda_v1([0.0, 0.5], [0.5, 0.5]) == 0.0
        assert gate_v1([0.0, 1.0], [0.5, 0.5], 0.5) == (NO_GO, ZERO_VETO)
        assert gate_v1([1.0, 1.0], [0.5, 0.5], 0.5)[0] == GO
        assert gate_v1([1.0], [1.0], 1.5) == (BLOCK, "LAMBDA_TAU_INVALID")
        with pytest.raises(LambdaV1Error) as e:
            lambda_v1([1.0], [1.5])
        assert e.value.code == "LAMBDA_WEIGHT_SUM"

    def test_failed_check_is_zero_veto(self):
        import sys; sys.path.insert(0, str(SRC))
        import formula_governance as FG
        original = FG.acoustic_axes
        FG.acoustic_axes = lambda: ([1.0, 0.0, 1.0], ["a", "b", "c"])
        try:
            r = FG.governed_validation(tau=0.5)
            assert r["verdict"] == "NO_GO" and r["code"] == "ZERO_VETO"
        finally:
            FG.acoustic_axes = original

    def test_governed_validation_go_with_honest_labels(self):
        import sys, json, subprocess
        r = subprocess.run([sys.executable, str(SRC / "formula_governance.py")],
                           capture_output=True)
        assert r.returncode == 0
        d = json.loads(r.stdout)
        assert d["verdict"] == "GO"
        assert d["lambda_status"] == "ADVISORY"
        assert d["claims"]["ns_regularity"] == "NOT_ESTABLISHED"
        assert "not a theorem" in d["conjecture_status"]


def test_governed_validation_invalid_weights_fail_closed_receipt():
    """Invalid weights must return a BLOCK receipt, not raise before the
    caller sees the decision (Phase A A1)."""
    import sys; sys.path.insert(0, str(SRC))
    from formula_governance import governed_validation
    r = governed_validation(0.5, weights=[1.5, -0.5])
    assert r["verdict"] == "BLOCK"
    assert r["code"] == "LAMBDA_LENGTH_MISMATCH"
    assert r["lambda"] == 0.0


class TestFailClosedReceiptContract:
    """szl.lambda/v1 (szl-lambda-gate reference/szl_lambda_v1.py): the gate never
    raises on bad input, it BLOCKs; tau is checked first; ports must agree with
    the reference on every error code and verdict."""

    def test_weight_sum_overflow_blocks_instead_of_raising(self):
        import formula_governance as FG
        # Finite, positive weights whose fsum overflows a float: the reference
        # reports LAMBDA_WEIGHT_SUM. An OverflowError must not escape the gate.
        with pytest.raises(FG.LambdaV1Error) as e:
            FG.lambda_v1([0.5, 0.5], [1e308, 1e308])
        assert e.value.code == "LAMBDA_WEIGHT_SUM"
        assert FG.gate_v1([0.5, 0.5], [1e308, 1e308], 0.5) == (FG.BLOCK, "LAMBDA_WEIGHT_SUM")

    def test_governed_validation_weight_overflow_returns_block_receipt(self):
        import formula_governance as FG
        r = FG.governed_validation(0.5, weights=[1e308] * 8)
        assert (r["verdict"], r["code"], r["lambda"]) == ("BLOCK", "LAMBDA_WEIGHT_SUM", 0.0)

    @pytest.mark.parametrize("tau,weights", [
        (1.5, [1.5, -0.5]),           # tau invalid + length mismatch
        (0.0, [0.2] * 8),             # tau invalid + weights sum to 1.6
        (float("nan"), [1e308] * 8),  # tau invalid + weight sum overflows
    ])
    def test_receipt_code_is_the_gate_code_when_tau_is_also_invalid(self, tau, weights):
        import formula_governance as FG
        scores, _ = FG.acoustic_axes()
        assert FG.gate_v1(scores, weights, tau) == (FG.BLOCK, "LAMBDA_TAU_INVALID")
        r = FG.governed_validation(tau, weights=weights)
        assert (r["verdict"], r["code"], r["lambda"]) == ("BLOCK", "LAMBDA_TAU_INVALID", 0.0)

    @pytest.mark.parametrize("axes,weights,tau,expected", [
        # Exact delta = log(Lambda) - log(tau) is +1.29e-4 and -1.36e-4 here,
        # over 1e5 times TIE_EPS, so neither case is a tie. tau is the
        # subnormal float that exp(log Lambda) rounds to: log() of that
        # rounded value would see delta == 0 and wrongly ABSTAIN.
        ([1e-320, 2e-320], [0.5, 0.5], 1.414e-320, ("GO", None)),
        ([3e-321, 7e-322, 1e-319], [0.25, 0.25, 0.5], 1.2045e-320, ("NO_GO", "BELOW_TAU")),
    ])
    def test_gate_compares_the_unrounded_log_lambda(self, axes, weights, tau, expected):
        import formula_governance as FG
        # spec gate/compare: "log space on the unrounded value".
        assert FG.lambda_v1(axes, weights) == tau
        assert FG.gate_v1(axes, weights, tau) == expected
        assert math.exp(FG.log_lambda_v1(axes, weights)) == FG.lambda_v1(axes, weights)
        assert FG.log_lambda_v1([0.0, 1.0], [0.5, 0.5]) == -math.inf


def _separable_rho_vL_sq(lam, mu, A, B, C, sigma):
    """A counterfactual stress law in which A, B and C enter with different
    powers of sigma, so velocity data COULD tell them apart."""
    denom = 3 * lam + 2 * mu
    return (lam + 2 * mu - (7 * lam + 10 * mu) / denom * sigma
            - (2 * A * sigma + 10 * B * sigma ** 2 + 6 * C * sigma ** 3) / denom)


class TestMurnaghanIdentifiabilityIsComputed:
    """The (A, B, C) identifiability result must be computed from the Murnaghan
    formula for the stress states supplied, not asserted as a constant."""

    def test_zero_stress_data_identify_nothing(self):
        from aero_ff import murnaghan_identifiability
        for stresses in ([0.0], [], [0.0, 0.0]):
            r = murnaghan_identifiability(100e9, 80e9, stresses, 7800)
            assert r["sensitivity_rank"] == 0, stresses
            assert r["identifiable_combination"] is None, stresses
            assert r["individual_ABC_identifiable"] is False, stresses

    def test_stressed_data_identify_exactly_S(self):
        from aero_ff import murnaghan_identifiability
        r = murnaghan_identifiability(100e9, 80e9, [0.0, 1e6, 5e6, 2e7], 7800)
        assert r["sensitivity_rank"] == 1
        assert r["sensitivity_direction"] == [1, 5, 3]  # parallel to (2, 10, 6)
        assert r["identifiable_combination"] == "S = 2A + 10B + 6C"
        assert r["individual_ABC_identifiable"] is False

    def test_rank_follows_the_formula(self, monkeypatch):
        # Falsifiability: under a stress law where A, B, C enter separately the
        # same analysis must report it (three distinct nonzero stresses -> rank 3).
        from aero_ff import murnaghan_identifiability
        monkeypatch.setattr("aero_ff.murnaghan_rho_vL_sq", _separable_rho_vL_sq)
        r = murnaghan_identifiability(100e9, 80e9, [1e6, 5e6, 2e7], 7800)
        assert r["sensitivity_rank"] == 3
        assert r["individual_ABC_identifiable"] is True
        assert murnaghan_identifiability(100e9, 80e9, [1e6, 5e6], 7800)["sensitivity_rank"] == 2

    def test_identifiability_axis_can_veto(self, monkeypatch):
        import formula_governance as FG
        monkeypatch.setattr("aero_ff.murnaghan_rho_vL_sq", _separable_rho_vL_sq)
        scores, names = FG.acoustic_axes()
        assert dict(zip(names, scores))["aero_identifiability_honest"] == 0.0
        r = FG.governed_validation(tau=0.5)
        assert (r["verdict"], r["code"]) == ("NO_GO", "ZERO_VETO")
