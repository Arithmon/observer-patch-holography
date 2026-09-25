# Observer dynamics: executable evidence interface

These tools verify bounded observer-like self-reading patches: local port
loads, repair attempts, boundary readback, accumulated records and public
receipt bundles. The immutable producers, specifications and compact numerical
inputs live in `../../evidence/observer_dynamics_20260925/`.

From the RER root:

```sh
python3 code/observer_dynamics/verify.py
python3 code/observer_dynamics/verify.py --replay-refinement --write-receipt
python3 evidence/exact_federation_tower_20260924/verify_tower.py
python3 evidence/source_net_causal_poset/verify_sampled_manifold.py
python3 evidence/observer_dynamics_20260925/lean/verify.py
```

The first command checks immutable file hashes, recomputes finite-chain and
finite-window receipts, checks the 48 terminal-state digest crosspins, verifies
scale-estimator identities and runs independent and mutation tests. The second
also reconstructs the fresh-noise streams for all 61 retained refinement chains
and writes a separate canonical verification receipt. Neither command modifies
historical evidence. The tower verifier streams load reconstruction in 8 MiB
chunks. The manifold command uses the original independently written verifier
and tests, recovered with hashes matching the sampled-run receipt.

`archive_adapter.py` answers only the historical producers' explicit `git show`
requests for the frozen simulator revision from the archived source blobs. No
sibling simulator checkout, network service or simulator Git history is needed.
It also serializes relative provenance paths with `/` on every operating system;
the frozen producer bytes and numerical comparisons remain unchanged. Install
the repository-root `requirements.txt` for the complete default controls,
including the pinned CAMB interpolation checks. Skipped controls fail verification.
The causal controls use a separate process because they vendor another source
revision. The tests can also be run through this entrypoint without introducing
those modules into an unrelated Python process.

The large terminal arrays and geometry caches remain in their original archive.
Digest crosspins and retained statistical identities do not constitute a fresh
microscopic replay of those arrays. Large transport spectra and boosted-clock
measurements are preserved with their original producers and checks; the
bounded verification reruns their small controls. The sampled manifold record
lacks the original large event samples and joint estimator covariance, so its
large-run uncertainty cannot be independently rebuilt here.

The conditional CMB comparison uses the retained standard-transfer calculation
and supplied source amplitude, tilt, cosmological background and temperature.
It is a reproduction of a conditional observer distribution and measurement
comparison, not a source-dynamics derivation of these inputs. Plot images are
not required by this numerical package. Exact finite-chain identities,
numerical receipts, retained large-run observations and compiled Lean witnesses
remain distinct evidence types.
