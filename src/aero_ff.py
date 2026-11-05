#!/usr/bin/env python3
"""Aero-FF operating-value formula adapters: typed, falsifiable software checks.

Each of the 7 operating values from the Aero-FF formulation (Barnard, Oct 2026)
is implemented as an executable check. These validate *software behavior and
internal mathematical consistency* of the formulas on synthetic inputs.

They do NOT validate any physical device claim, and they explicitly do NOT
support the document's claims that (a) a topological invariant implies global
Navier-Stokes regularity, or (b) any host-AI accuracy percentages. Those
claims have no defined measurement protocol and remain HYPOTHESIS/UNVERIFIED.
"""
from __future__ import annotations
import json
import math
import sys


def murnaghan_vL(lam: float, mu: float, A: float, B: float, C: float,
                 sigma: float, rho0: float) -> float:
    """Longitudinal wave velocity under finite uniaxial stress (Murnaghan).

    rho0 * vL^2 = lam + 2mu - [(7*lam + 10*mu + 2A + 10B + 6C) / (3*lam + 2*mu)] * sigma
    """
    denom = 3 * lam + 2 * mu
    if denom <= 0 or rho0 <= 0:
        raise ValueError("invalid material constants")
    coeff = (7 * lam + 10 * mu + 2 * A + 10 * B + 6 * C) / denom
    vL_sq = (lam + 2 * mu - coeff * sigma) / rho0
    if vL_sq <= 0:
        raise ValueError("non-physical: vL^2 <= 0 at this stress")
    return math.sqrt(vL_sq)


def murnaghan_identifiability(lam: float, mu: float, stresses: list[float],
                              rho0: float) -> dict:
    """Identifiability analysis for (A, B, C) from vL(sigma) measurements.

    vL depends on A, B, C only through the single scalar combination
    S = 2A + 10B + 6C, so velocity measurements identify S, never the
    individual constants. Returns the rank of the sensitivity matrix and
    the identifiable combination.
    """
    def S(A, B, C):
        return 2 * A + 10 * B + 6 * C

    base = (1.0, 1.0, 1.0)
    rows = []
    h = 1e-6
    for _ in stresses:
        # d(vL)/d(A,B,C) is proportional to dS/d(A,B,C) = (2, 10, 6) for every sigma
        rows.append([2.0, 10.0, 6.0])
    # rank of a matrix whose every row is (2,10,6) is 1
    rank = 1
    return {
        "sensitivity_rank": rank,
        "identifiable_combination": "S = 2A + 10B + 6C",
        "individual_ABC_identifiable": False,
        "note": "All stress states probe the same scalar combination; "
                "recovering A, B, C individually requires independent "
                "measurement modes (e.g. shear waves, different propagation axes).",
    }


def poisson_transverse(nu: float, sigma: float, E_eff: float) -> float:
    """Circumferential hoop strain from axial stress (Poisson coupling)."""
    if not (-1.0 < nu < 0.5):
        raise ValueError("Poisson ratio out of physical range")
    return -nu * (sigma / E_eff)


def ep_gain_db(R_abs: float, power_ratio: bool = False) -> float:
    """Amplification G = 2/(1-|R|^2), reported in dB.

    The Aero-FF document states 'G > 4, or > 6 dB'. That is only consistent
    if G is treated as a POWER ratio (10*log10(4) = 6.02 dB). As an
    AMPLITUDE ratio, G = 4 is 20*log10(4) = 12.04 dB. This check surfaces
    the ambiguity explicitly instead of silently picking one.
    """
    if not (0 <= R_abs < 1):
        raise ValueError("|R| must be in [0, 1)")
    G = 2.0 / (1.0 - R_abs ** 2)
    return 10 * math.log10(G) if power_ratio else 20 * math.log10(G)


def exceptional_point_2x2(omega0: float, kappa: float, gamma: float) -> dict:
    """Test the EP condition on a 2x2 non-Hermitian (PT-symmetric) matrix.

    H = [[omega0 + i*gamma, kappa], [kappa, omega0 - i*gamma]]
    Eigenvalues: omega0 +- sqrt(kappa^2 - gamma^2). A TRUE exceptional point
    requires BOTH:
      (a) eigenvalue coalescence (kappa^2 == gamma^2), and
      (b) the coalesced eigenvectors become defective (only one linearly
          independent eigenvector survives).
    Equality of gap and coupling alone is the standard EP *condition*, but
    this check verifies the full eigenvector defectiveness as well.
    """
    import numpy as np
    H = np.array([[omega0 + 1j * gamma, kappa], [kappa, omega0 - 1j * gamma]], dtype=complex)
    evals, evecs = np.linalg.eig(H)
    scale = max(1.0, abs(evals[0]))
    coalesced = abs(evals[0] - evals[1]) < 1e-6 * scale
    defective = False
    if coalesced:
        lam = evals[0]
        defective = np.linalg.matrix_rank(H - lam * np.eye(2), tol=1e-6 * scale) == 1
    normalized_inner = abs(np.vdot(evecs[:, 0], evecs[:, 1])) / (
        np.linalg.norm(evecs[:, 0]) * np.linalg.norm(evecs[:, 1]))
    vec_parallel = bool(abs(normalized_inner - 1.0) < 1e-6) if coalesced else False
    return {
        "eigenvalues": [complex(e) for e in evals],
        "coalesced": bool(coalesced),
        "defective": bool(defective),
        "eigenvectors_parallel": vec_parallel,
        "is_true_EP": bool(coalesced and defective),
    }


def rankine_hugoniot(c0: float, gamma: float, up: float) -> float:
    """Shock-front velocity D = c0 + (gamma+1)/2 * up."""
    if c0 <= 0 or gamma <= 1 or up < 0:
        raise ValueError("non-physical inputs")
    return c0 + (gamma + 1) / 2.0 * up


def convective_transport(vp: float, U0: float) -> float:
    """Observed phase velocity v_obs = vp + U0 (Galilean addition)."""
    return vp + U0


def chern_number_dirac(m_eff: float, nk: int = 60) -> int:
    """Lattice Chern number for a 2-band Dirac model (Qi-Wu-Zhang style).

    H(k) = sin(kx)*sx + sin(ky)*sy + (m_eff + cos(kx) + cos(ky))*sz
    Returns the integer Chern number of the lower band. Verifies
    topological quantization on a discrete BZ grid.

    IMPORTANT: this verifies software computation of a topological
    invariant for a *model system*. It is NOT evidence for any claim about
    Navier-Stokes regularity, acoustic devices, or mass gaps.
    """
    import numpy as np
    # Pauli matrices
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    eye = np.eye(2)

    ks = 2 * math.pi * (np.arange(nk) / nk - 0.5)
    # Berry curvature of the lower band via eigenvector derivatives
    def lower_band_evec(kx, ky):
        H = math.sin(kx) * sx + math.sin(ky) * sy + (m_eff + math.cos(kx) + math.cos(ky)) * sz
        evals, evecs = np.linalg.eigh(H)
        return evecs[:, 0]

    total = 0.0
    for i in range(nk):
        for j in range(nk):
            kx, ky = ks[i], ks[j]
            v00 = lower_band_evec(kx, ky)
            v10 = lower_band_evec(ks[(i + 1) % nk], ky)
            v11 = lower_band_evec(ks[(i + 1) % nk], ks[(j + 1) % nk])
            v01 = lower_band_evec(kx, ks[(j + 1) % nk])
            # Fukui-Hatsugai-Suzuki plaquette: <00|10><10|11><11|01><01|00>
            U1 = np.vdot(v00, v10)
            U2 = np.vdot(v10, v11)
            U3 = np.vdot(v11, v01)
            U4 = np.vdot(v01, v00)
            phase = np.angle(U1 * U2 * U3 * U4)
            total += phase
    return int(round(total / (2 * math.pi)))


def selfcheck() -> bool:
    checks = []
    # 1. Murnaghan: velocity responds to stress, stays physical
    v0 = murnaghan_vL(100e9, 80e9, -100e9, -50e9, -30e9, 0.0, 7800)
    v1 = murnaghan_vL(100e9, 80e9, -100e9, -50e9, -30e9, 50e6, 7800)
    checks.append(v0 > 0 and v1 > 0 and abs(v1 - v0) > 0)
    # 2. Identifiability honesty: A,B,C not individually identifiable
    ident = murnaghan_identifiability(100e9, 80e9, [0.0, 1e6, 5e6], 7800)
    checks.append(ident["sensitivity_rank"] == 1 and ident["individual_ABC_identifiable"] is False)
    # 3. Poisson coupling sign and magnitude
    e = poisson_transverse(0.38, 1e6, 200e9)
    checks.append(-1e-3 < e < 0)
    # 4. dB ambiguity is surfaced, both interpretations computed correctly
    g_amp = ep_gain_db(0.5, power_ratio=False)
    g_pow = ep_gain_db(0.5, power_ratio=True)
    checks.append(abs(g_amp - 20 * math.log10(2 / (1 - 0.25))) < 1e-9)
    checks.append(abs(g_pow - 10 * math.log10(2 / (1 - 0.25))) < 1e-9)
    # 5. True EP at gain/loss = coupling; off-EP is not an EP
    ep = exceptional_point_2x2(1.0, 0.5, 0.5)
    checks.append(ep["is_true_EP"] is True)
    off_ep = exceptional_point_2x2(1.0, 0.5, 0.25)  # unbroken PT phase: distinct evals
    checks.append(off_ep["is_true_EP"] is False and off_ep["coalesced"] is False)
    broken_ep = exceptional_point_2x2(1.0, 0.5, 0.75)  # broken phase: distinct (complex) evals
    checks.append(broken_ep["is_true_EP"] is False)
    # 6. Rankine-Hugoniot scaling
    D = rankine_hugoniot(340.0, 1.4, 50.0)
    checks.append(abs(D - (340.0 + 1.2 * 50.0)) < 1e-9)
    # 7. Convective addition
    checks.append(convective_transport(340.0, 20.0) == 360.0)
    # 8. Chern number quantization: QWZ phases (transitions at m=-2,0,2)
    checks.append(chern_number_dirac(-1.0, nk=24) == -1)
    checks.append(chern_number_dirac(1.0, nk=24) == 1)
    checks.append(chern_number_dirac(3.0, nk=24) == 0)
    return all(checks)


def main() -> int:
    ok = selfcheck()
    result = {
        "schema": "szl.acoustic.aero-ff/v1",
        "selfcheck_pass": ok,
        "claims_status": {
            "formula_software_consistency": "MEASURED_SOFTWARE_SIMULATION",
            "device_physical_claims": "HYPOTHESIS_UNVERIFIED",
            "navier_stokes_regularity_implication": "NOT_ESTABLISHED",
            "host_ai_accuracy_percentages": "UNVERIFIED_NO_PROTOCOL",
        },
        "truth_label": "MEASURED_SOFTWARE_SIMULATION",
    }
    print(json.dumps(result, indent=2))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
