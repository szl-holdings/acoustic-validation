# Acoustic Validation

**SZL × Chad Barnard: Falsifiable acoustic framework with governed verification.**

> Turn acoustic and fluid hypotheses into traceable, reproducible tests with IMMUNE-sealed receipts.

## What this is

A validation framework for the Aero-FF (Acoustic Quantum Hall) framework and related acoustic/fluid dynamics claims. Every test result is sealed into an append-only SHA-256 receipt chain via IMMUNE, making the validation itself tamper-evident.

## Architecture

```
Acoustic Validation (this repo)
  ├── Burgers equation solver (analytic residual + grid refinement)
  ├── Duffing oscillator (coefficient + ghost-identifiability checks)
  ├── Aero-FF framework tests (falsifiable predictions)
  └── IMMUNE receipt sealing (SHA-256 hash chain)
        ├── NEXUS runtime (analog simulation + visualization)
        └── Cloudflare x402 monetization ($0.05/validation)
```

## Validated components (20/20 baseline passing)

| Test | Status | Method |
|---|---|---|
| Burgers analytic residual | ✅ PASS | Compare numerical vs analytic solution |
| Burgers grid refinement | ✅ PASS | Convergence under mesh refinement |
| Duffing coefficient recovery | ✅ PASS | Parameter identification from synthetic data |
| Duffing ghost identifiability | ✅ PASS | Distinguish true vs ghost solutions |

## What is NOT claimed

- No Navier-Stokes or Yang-Mills proof
- No device physical claim validation
- No AI performance validation
- Upstream SciPy/PINN modules not executed
- GPU execution not qualified
- Failed-test injection leaves claims unverified

## Integration points

| System | Role |
|---|---|
| [NEXUS](https://github.com/szl-holdings/nexus) | Runtime for acoustic simulations |
| [IMMUNE](https://github.com/szl-holdings/immune) | Receipt chain + verification plane |
| [szl-serve](https://github.com/szl-holdings/szl-serve) | API serving layer |
| Cloudflare x402 | Per-call monetization |

## Monetization

Validation API available at `api.a-11-oy.com/v1/acoustic/validate`:
- $0.05 per validation call (x402)
- $0.10 for verified bots (x402)
- Results include IMMUNE-sealed receipt

## Evidence requested from Chad

1. Raw synchronized recordings (500+ trials, including failures)
2. Sensor types, sample rates, calibration records
3. Device geometry, materials, energy budget measurements
4. Source code and Lean proofs
5. AI results: dataset provenance, train/test separation, metrics

## Collaborators

- **Stephen P. Lutar Jr.** (SZL Holdings) — validation framework, governed infrastructure
- **Chad Barnard** (Molecular Audio Canada) — acoustic research, device observations
- **Paul Guarino** — Aero-FF framework author

## License

Apache-2.0
