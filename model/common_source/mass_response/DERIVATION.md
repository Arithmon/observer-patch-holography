# Structural dispersion and finite interacting charge-transfer gaps

This construction maps native reference log odds to a declared coherent
dimer continuation. The architecture is chosen and its consequences are
computed. The structural gap is not spontaneous mass generation, an observed
Standard Model mass, or a relabeling of native repair relaxation.

The source reference has per-state weights
`w0=2/225`, `w1=1/36`, `w2=2/45`. Define

```text
a = log(w1/w0) = log(25/8),
b = log(w2/w1) = log(8/5),
J = a-b = log(125/64),
Omega = a+b = log(5).
```

`J` is the connected occupation coefficient of that reference. Using the
two log odds as coherent hopping strengths is an additional, explicit
architecture prescription. Native endpoint repair alone does not impose it.
The reference supplies no empirical particle mass. A factorized-reference
control sets `a=b=Omega/2` while preserving the endpoint ratio and clock
frequency; this is a changed reference/control, not a freedom at fixed full
modular action.

## PRIMARY finite interacting ring

The primary four-mode ring has hopping `-a` on `(0,1),(2,3)` and `-b` on
`(1,2),(3,0)`, and interaction

```text
g sum_cycle (n_i-1/2)(n_j-1/2),    a>b>0, g>=0.
```

The analytic formulas in [dispersion.py](dispersion.py) are independent of
the numerical CAR construction. The particle-number sector dimensions are
`1,4,6,4,1`. Empty and full sectors have energy `g`; one- and three-particle
sectors have energies `(-Omega,-J,J,Omega)`. In the half-filled sector the
characteristic polynomial is

```text
lambda^2 (lambda^2+g lambda-4a^2)(lambda^2+g lambda-4b^2).
```

To obtain it directly, use occupation masks `(3,5,6,9,10,12)`. Reflection
combinations split the matrix into two zero lines and the blocks
`[[-g,-2a],[-2a,0]]`, `[[-g,-2b],[-2b,0]]`. The sign of each off-diagonal
entry can also be reversed by a basis phase without changing its spectrum.
The independent test constructs the six-dimensional occupation matrix,
including the boundary-hop CAR signs, before comparing its eigenvalues with
the analytic expressions.

Consequently the full finite-system ground state has two particles and
energy `E2=-(g+sqrt(g^2+16a^2))/2`. The charged addition and removal gaps,
relative to the particle-hole-symmetric energy convention, are

```text
Delta4 = (g+sqrt(g^2+16a^2))/2-Omega.
```

A single-ring odd-charge transition requires a compensating charge/parity
reference. The primary charge-conserving measurement instead compares two
identical copies with occupations `(2,2)` and `(3,1)`. It has frequency
`E3+E1-2E2=2 Delta4`, with the same total charge four on both arms.
Adding a common chemical potential `mu*N` cancels exactly from that
frequency. This global charge cancellation is not a proof of a local
electromagnetic gauge-field realization.
Inferring one charged gap as half the measured pair frequency also uses the
identical-copy and particle-hole-symmetry contracts; it is not a separate
measurement of an infinite-volume charged pole.

The lowest neutral excitation within one half-filled ring is a different
observable:

```text
(sqrt(g^2+16a^2)-sqrt(g^2+16b^2))/2.
```

For `g=0`, `Delta4=J`. For `g=J`, the finite charged gap is approximately
`1.028595953717031`, whereas the single-ring neutral gap is approximately
`1.305496913774672`. The interaction changes the reference and excitation
energies. Neither the free structural gap nor the neutral gap may be silently
substituted for the measured interacting charge-transfer frequency.

This finite ring is a small interacting system. Its exact gap does not
establish the rest mass or dressed pole of an infinite homogeneous theory.
Volume dependence, spatial propagation and the relevant vacuum/clock
identification require their own shared-history argument.

## Separate FREE three-axis extension

Six Hermitian Clifford matrices `Gamma_1,...,Gamma_6` on eight components
satisfy `{Gamma_i,Gamma_j}=2 delta_ij I`. Define the free Bloch operator

```text
H(k) = sum_j [(a-b cos k_j) Gamma_(2j-1)
             + b sin k_j Gamma_(2j)].
```

The minus-intercell convention is explicit. Its minimum is at `k=0`, which
belongs to every periodic finite torus. Anticommutation gives

```text
H(k)^2 = [3J^2 + 4ab sum_j sin(k_j/2)^2] I8.
```

The free positive-band gap is therefore `sqrt(3)|J|`, independent of volume.
All eight components remain present; the construction neither discards
extra species nor identifies physical spin, Weyl chirality or a Standard
Model multiplet. The zero-imbalance control is gapless at the same clock
frequency. The negative bands also require a stated many-body reference;
an empty Fock state is not automatically the ground state of this hopping
Hamiltonian.

For a refining family, set spacing to `1/n`, keep coordinate momentum `p`
fixed and choose

```text
a_n=(n Omega+J)/2, b_n=(n Omega-J)/2, k=p/n.
```

These are explicit architecture and scale prescriptions. The exact square
approaches `3J^2+(Omega^2/4)|p|^2`. With `x_j=p_j/(2n)`, its difference from
that limit is

```text
n^2 Omega^2 sum_j [x_j^2-sin(x_j)^2] + J^2 sum_j sin(x_j)^2.
```

The inequalities `0<=x^2-sin(x)^2<=x^4/3` and `sin(x)^2<=x^2` prove

```text
0 <= E_cont(p)^2-E_n(p)^2
  <= J^2|p|^2/(4n^2) + Omega^2 sum_j p_j^4/(48n^2).
```

For example, the first inequality follows by integrating
`d(x-sin(x)cos(x))/dx=2sin(x)^2<=2x^2` from zero and then integrating once
more. This proves the mathematical all-size bound. Floating numerical
checks evaluate it but are not outward roundoff certificates.

The free extension provides a dispersion and refinement control. It is not
the primary interacting two-copy execution, does not inherit its readout
law, and proves no interacting three-dimensional continuum limit or
physical Lorentzian field identification. In particular, the growing
hopping strengths require accounted native source operations and drive
resources; naming a refinement level does not provide those services.
