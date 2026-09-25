# Finite local-force scalar preparation

This package certifies a finite rectangular force on the sixteen preparation
ports of the original 64-site scalar model. The bounded software observer
controls local force registers, reads mass/velocity/duration coefficients,
switches each force off, and retains a public certificate. It supplies a
conditional preparation channel for the existing field/record instrument.
The free field Hamiltonian remains active throughout the pulse.

```sh
python3 code/source_scalar_preparation/preparation.py --write
python3 code/source_scalar_preparation/verify_preparation.py
python3 -m pytest -q code/source_scalar_preparation/test_preparation.py
```

The producer writes `preparation_receipt.json`; `--check` checks deterministic
parity. The verifier uses independent exact arithmetic in the `1,sqrt(5)` basis,
checks all 64 action rows, and reconstructs the error budget. It pins the parent
receipts; their own independent verifiers validate the inherited instrument.
The mathematical derivation is in `SOURCE_SCALAR_LOCAL_PREPARATION.tex`.

The pulse is centered on virtual reference time zero. The prepared state is
compared with the free coherent reference only after the pulse ends; it is
not that reference state at the actual pulse midpoint. All later bounded
instrument probabilities inherit the trace-distance bound, including every
correlated readout in the retained 21-slot instrument.

The initial vacuum, quantization, admissible force coupling, relative action
scale, operational time and implementation-error limits are inputs. No Born
bits or hardware execution are synthesized. The ideal added-energy bound
does not bound noisy energy from trace distance alone. The finite control
schedule is a specified operation contract, not an authenticated quantum run.

