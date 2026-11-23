# Leader-Study Upgrade Plan - W27 loop-4 - 2026-11-23

Status: PROPOSE-ONLY DRAFT - DO NOT MERGE. Owner-gated: merging requires governing approval and the estate signature-wall owner unlock.

Repo: acoustic-validation - Acoustic validation framework (software only): Burgers equation, Duffing oscillator and Aero-FF falsifiable tests. No physical claim is validated; this repository has no receipts, runtime or paid API.

Frame: OUROBOROS loop-4 leader study. Governing seat: Kimi kimi-k3.1-pro (json_ok, think off). Advisory seats: Astra/Sol/Claude community rails (provenance unverified, advisory only). Leader evidence web-verified 2026-11-23.

## Leader moats to absorb
- ASME V&V 10/20/40: Verification and validation standards for computational solid mechanics, fluid dynamics, and heat transfer, including code verification, solution verification, and validation against accuracy requirements. Moat: Standardized credibility framework separating numerical correctness from physical validation, with explicit uncertainty and acceptance criteria.
- SciPy and FEniCS: Established scientific-computing baselines for numerical methods, regression testing, and reproducible solver comparisons. Moat: Mature implementations, broad peer scrutiny, and reusable reference behavior for independent cross-checks.
- NVIDIA Modulus: Physics-ML framework for developing and evaluating physics-informed models. Moat: Integrated tooling for physics-informed learning and benchmark-oriented experimentation; relevant only as a possible comparison target, not as evidence of physical validation.
- OpenFOAM V&V practice: Established CFD verification and validation workflows and best practices. Moat: Operational discipline around mesh convergence, numerical uncertainty, benchmark cases, and documented validation limits.
- Burgers and Duffing benchmark suites: Reference problems and expected behaviors for nonlinear wave propagation and nonlinear oscillation. Moat: Falsifiable, widely understood test targets that make solver regressions and misleading comparisons easier to detect.

## Honest gaps (grounded in repo facts)
- The repository currently contains only two self-contained solver areas, Burgers and Duffing; the facts do not establish implemented Aero-FF tests beyond the repository title and stated prototype purpose.
- All verification is on synthetic inputs and software behavior; no physical claim, device claim, exceptional point, quantum-Hall analogue, or AI performance claim is validated.
- No ingestion of Chad's recordings has occurred, and rights and calibration remain unresolved, so there is no real acoustic measurement evidence.
- The controlled measurement campaign and private analysis pipeline are roadmap items, not implemented capabilities; originals-plus-hashes preservation, units/setup capture, spectra, transfer functions, room-effect comparison, and k-Wave propagation are not established as present.
- The README explicitly states that there are no IMMUNE receipt sealing, NEXUS runtime integration, API endpoints, or x402 monetization in this repository despite related topics.
- The existing Burgers analytic comparison is valid only before wave-breaking at t < 1/(2π); the facts do not support a claim of correctness across shocked or post-breaking regimes.
- Grid-refinement convergence at fixed CFL is described, but grid-convergence-index calculation, formal observed-order reporting, and numerical uncertainty budgets are not stated as implemented.
- Duffing coefficient recovery and ghost rejection are verified on synthetic trajectories; robustness to measured noise, missing data, sampling variation, model mismatch, and experimental calibration is not established.
- The facts state that CI and pinned dependencies are roadmap items marked complete, but provide no independent run receipts or test results; this audit therefore cannot claim a fresh execution or passing CI.
- No external benchmark cross-check against SciPy, FEniCS, OpenFOAM, Modulus, or published reference suites is stated as implemented.

## Dosed upgrade plan
- Dose 1 — evidence hygiene: keep the present scope unchanged and add a machine-readable verification report for each selfcheck containing Git commit, dependency versions, solver settings, CFL, mesh sizes, tolerances, errors, exit status, and explicit synthetic-data labeling; do not call it a sealed receipt.
- Dose 2 — deepen current verification: extend Burgers refinement reporting with observed convergence order and GCI where applicable, retain the pre-breaking analytic limit, and add clearly separated post-breaking tests only against an appropriate independent reference rather than the invalid smooth analytic solution.
- Dose 3 — strengthen Duffing falsification: add synthetic noise, sampling-rate, trajectory-length, and parameter-regime sweeps; report recovery bias and uncertainty; retain wrong-sign ghost rejection and failure injection as mandatory CI gates.
- Dose 4 — independent benchmarks: compare both solvers with published reference cases and, where useful, SciPy or another independent implementation; record expected tolerances and discrepancies before considering FEniCS, OpenFOAM, k-Wave, or Modulus.
- Dose 5 — clarify Aero-FF scope: either implement a small, explicitly defined falsifiable Aero-FF software test with inputs, equations, expected outputs, and rejection criteria, or state precisely that Aero-FF remains a target rather than an implemented validated module.
- Dose 6 — prepare, but do not execute, the measurement pipeline: specify consent/rights gates, calibration requirements, uncompressed PCM, gain and units metadata, immutable originals, hashes, and analysis provenance; keep all physical-validation claims blocked until those gates and controlled comparisons are satisfied.
- Dose 7 — defer integrations: propose NEXUS, IMMUNE, API, or x402 work only as separate explicitly scoped projects with their own tests and receipts; do not represent them as capabilities of this repository.
- Risk note: The main risk is credibility leakage: repository topics and roadmap language could be mistaken for implemented Aero-FF validation, acoustic evidence, receipt sealing, runtime integration, or monetization. Numerical overclaiming is also possible if pre-breaking Burgers agreement or synthetic Duffing recovery is generalized beyond its tested regime. Additional dependencies and benchmark frameworks can increase complexity without improving evidence unless introduced incrementally with independent references, explicit tolerances, and preserved negative controls.
- Propose-only note: This is a proposal only. No repository files, tests, CI jobs, external benchmarks, recordings, integrations, receipts, or paid endpoints were executed or changed for this audit. Each dose should land behind review with its own acceptance criteria and a rollback path; physical validation remains out of scope until rights, calibration, controlled measurements, and validation requirements are actually satisfied.

## Advisory notes (unverified rails, non-binding)
- astra: acoustic-validation: Proprietary high-fidelity benchmark datasets for PDE solver verification. ayllu: Formal verification proofs for multi-agent consensus and evidence-binding. aifor-science-protocol: Institutional adoption as the standard for AI-driven research pre-registration.
- sol: acoustic-validation: Proprietary high-fidelity synthetic datasets for benchmark validation. ayllu: Formal verification proofs for multi-agent consensus and evidence-binding. aifor-science-protocol: Standardized, machine-readable registration schema adopted by major journals.
- claude: acoustic-validation: Proprietary benchmark datasets of verified physical ground-truth measurements. ayllu: A standardized, cryptographically verifiable audit trail for multi-agent consensus. aifor-science-protocol: Institutional endorsement and integration into peer-review submission pipelines.

## DO NOT MERGE
This plan is a proposal artifact of the 10-pass OUROBOROS loop (loop-4 of 10). No code, no CI, no behavior changes in this file. A reviewer with authority may convert to a work order later.