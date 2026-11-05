# Acoustic Validation

**SZL × Chad Barnard: falsifiable acoustic framework with governed verification.**

> Turn acoustic and fluid hypotheses into traceable, reproducible tests. Software results only — no physical claim is validated here.

## What this is

A standalone numerical validation prototype for the Aero-FF framework and related acoustic/fluid dynamics claims. It currently contains two self-contained solvers with honest, independent checks:

- **Burgers equation (inviscid, sin initial data)** — Lax-Friedrichs solver compared against the *time-evolved* analytic characteristic solution (valid before wave-breaking at `t < 1/(2π)`), plus grid-refinement convergence at fixed CFL.
- **Duffing oscillator** — RK4 simulation with least-squares coefficient recovery (signs handled explicitly) and a ghost-identifiability check that rejects sign-flipped candidate models from data.

Both selfchecks return a failing exit code when checks fail, so CI can gate on them.

## What is verified vs not

**Verified (software behavior on synthetic inputs):**
- Numerical solutions match the time-evolved analytic solution within expected discretization error.
- Errors decrease monotonically under mesh refinement at fixed CFL.
- Duffing coefficients are recovered from synthetic trajectories within tolerance, including when true parameters are changed.
- Ghost (wrong-sign) models are rejected by prediction residuals, not by initial-condition differences.
- Broken solvers (wrong flux sign, negated coefficients) fail the checks (failure-injection tests).

**NOT claimed:**
- No Navier–Stokes or Yang–Mills proof.
- No validation of any device physical claim, exceptional point, or quantum-Hall analogue.
- No AI performance validation.
- No ingestion of Chad's recordings — rights and calibration are unresolved; see "Evidence requested" below.
- No IMMUNE receipt sealing, no NEXUS runtime integration, no API endpoint, and no x402 monetization are present in this repository. Those are separate systems; any claim that this repo seals receipts or serves a paid API would be false until they are wired up and verified.

## Run it

```
pip install -r requirements.txt pytest
pytest tests/ -v
python src/burgers.py   # exits 0 on pass, 1 on fail
python src/duffing.py    # exits 0 on pass, 1 on fail
```

## Roadmap (in dependency order)

1. ✅ Numerical repairs: time-evolved analytic comparison, sign-correct recovery, meaningful ghost rejection, failing exit codes, pinned deps, CI.
2. Private analysis pipeline: preserve originals + hashes, record units/setup, spectra and transfer functions (SciPy/pyfar), room-effect comparison (pyroomacoustics), k-Wave for defined nonlinear propagation.
3. Controlled measurement campaign spec: uncompressed PCM, documented gain, calibration, known stimulus, measured geometry, repeat trials and controls.
4. SZL evidence layer: bind inputs, parameters, revision, uncertainty and permissions to each result; typed formula adapters with explicit abstention.
5. Only then: a narrow surrogate model (e.g., frequency response) against physics baselines, with session-level holdout and OOD tests.

## Evidence requested from Chad

1. Raw synchronized recordings (500+ trials, including failures and exclusion rules)
2. Sensor types, sample rates, clock sync, calibration records, uncertainty estimates
3. Device geometry, materials, energy-budget measurements
4. Source code, exact reproduction steps, any complete Lean project
5. AI results: dataset provenance, train/test separation, baselines, metrics, scripts

**Rights note:** Chad's manuscripts grant conditional verification permission with ownership and advocacy conditions — not a blanket grant to train on or redistribute his material. No recordings may be ingested into any model until permitted uses are settled in writing.

## Collaborators

- **Stephen P. Lutar Jr.** (SZL Holdings) — validation framework, governed infrastructure
- **Chad Barnard** (Molecular Audio Canada) — acoustic research, device observations
- **Paul Guarino** — Aero-FF framework author

## License

Apache-2.0
