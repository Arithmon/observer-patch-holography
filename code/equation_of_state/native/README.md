# Native scalar repair and work-law identification

This package independently verifies a bounded twelve-port software patch with
local scalar state, seam reads and mean-repair writes. Every read/write is retained
in an external diagnostic record. It is an exact primitive-level evidence bundle,
not a thermal bath or a full self-reading physical federation.

The package separates three things:

- Native execution conserves total scalar load and reduces
  `Q(x) = sum_i x_i² / 2` by exactly `(x_i-x_j)²/4` on a repaired seam.
- Recording that decrease in an added accumulator gives exact bookkeeping
  `Q + B = Q_initial`. Neither `Q` nor `B` is identified with physical energy
  or heat.
- Supplied positive work-energy extensions agree with `Q` on every state at
  the reference volume, but have different volume derivatives and candidate
  pressure/density ratios. Their ratios are illustrative assumed work laws,
  not native or measured thermodynamic equations of state. Even restricting
  to strictly positive candidate pressure leaves `1/3` and `1/2` distinct.

The work family, volume, energy identification and derivative convention are
explicit inputs. The package establishes neither thermal equilibrium nor an
extensive material law. It excludes no enriched repair, physical action or
alternative OPH construction. All native thermodynamic fields remain JSON `null`.

Run from the research repository root:

```bash
python3 code/equation_of_state/native/native_eos.py
python3 code/equation_of_state/native/verify_native_eos.py --sim-root ../oph-physics-sim \
  --out code/equation_of_state/native/native_eos_verification.json
python3 -m unittest discover -s code/equation_of_state/native -p 'test_native_eos.py' -v
```

The producer imports the existing simulator's canonical matrix primitive.
Its environment therefore needs the simulator dependencies. The verifier and
tests need only Python's standard library; they import neither producer nor
simulator. Verification reconstructs the thirty seams from exact icosahedral
coordinates in `Q(phi)`, with `phi²=phi+1`, and directly replays every rational
mean. `--sim-root` additionally checks the recorded simulator source bytes.
Source hashes identify actual imported files; the recorded checkout commit alone
does not assert a clean working tree.

`native_eos_receipt.json` retains all 240 attempted updates: two complete sweeps
for each of ascending pulse, descending pulse, constant and disabled-repair
controls. It also includes forty initial/final pressure jets, an energy-offset
control and a fixed direction-weighted scalar moment control. The pressure jets
evaluate the analytic derivative of the declared work family; they do not use
numerical differentiation or represent expansion experiments. The verifier
checks complete scope fields and the tests include coherently rehashed false
pressure, hidden heat and equilibrium promotions, omitted controls, altered
read/write histories and source drift.

The analytic work-extension and quadratic-error derivation is retained in the
producer module docstring; the executable witness remains separate from physical
energy identification.
To inspect the data independently, begin with the JSON `contract`, `claim_boundary`
and `work_extension` objects before reading any numerical ratios.
