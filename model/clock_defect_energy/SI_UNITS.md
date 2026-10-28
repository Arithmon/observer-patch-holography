# Seconds, metres and kilograms for a declared clock

The [SI conversion utility](si_units.py) uses the exact defining constants
`c=299792458 m/s` and `h=6.62607015e-34 J s` from the
[BIPM](https://www.bipm.org/en/measurement-units/si-defining-constants).
Its decimal outputs use floating arithmetic. The input `T` is the declared
seconds per unit of one specified model clock; it is not solver time or host
runtime. This is a conditional unit attachment, not a measured calibration.

If the model's signal speed is `v` model lengths per clock unit, identifying
that signal with light gives `L=c T/v` metres per model length. Under the
additional physical identification `E=hbar omega`, a rest angular frequency
`omega0` per model clock unit gives

```text
E_rest = hbar omega0 / T joules
m_rest = hbar omega0 / (T c^2) kilograms.
```

For the illustrative convention `T=1 second` and `v=1`, a model length is
299792458 metres, and the mass conversion coefficient is approximately
`1.17337e-51 kg` per unit rest angular frequency. This coefficient is not a
particle mass. The clock calibration, coherent rest gap and energy-frequency
identification must come from the same proposed physical process. A positive
repair count or a classical relaxation rate does not supply that gap.

```sh
python -m model.clock_defect_energy.si_units --seconds-per-unit 1 --signal-speed 1
```

Without a supplied rest frequency the command reports no mass. A specified
candidate can additionally pass `--rest-angular-frequency` and an independently
extracted `--kinetic-curvature`. For the convention where angular wavenumber
is measured per model length, a massive relativistic branch must satisfy
`omega0 * omega_second_derivative(0) / v^2 = 1`. Changing the unit attachment
cannot repair a failed dimensionless check. This utility applies no physical
mass label to the screened actively corrected qubit.
