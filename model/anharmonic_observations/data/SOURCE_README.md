# Dataset for Programmable quantum simulation of anharmonic dynamics

The CSV files in this dataset contain the experimental data presented in
"Programmable quantum simulation of anharmonic dynamics" by McGarry et al.
2026.

The file `symmetric_well_characteristic_functions_for_x.csv` contains the
characteristic functions measured to produce the data in Fig. 4 (a).

The file `symmetric_well_characteristic_functions_for_P.csv` contains the data
used to construct the probability distributions shown in Fig. 4 (b-d). This is
the same data that is plotted in Fig. 7.

The file `asymmetric_wells_expectation_x.csv` contains the data plotted in Fig.
5.

The remaining files contain the characteristic functions that are plotted in
Fig. 8.

## File contents

### Symmetric well time series data

File: `symmetric_well_characteristic_functions_for_x.csv`

Figure: 4(a)

This file contains 40 data series, with a header row. Each series 
represents a scan of the imaginary part of a characteristic function, the slope
of which is used to determine the expectation value of the position operator,
as described in the text.

The header row gives the value at which the characteristic function is
measured. The first column gives the time through the evolution for that
series. The remaining columns each hold one data series.

### Symmetric well probability data

File: `symmetric_well_characteristic_functions_for_P.csv`

Figure: 7

This file contains six data series with two header rows. Each series
corresponds to the real or imaginary part of a characteristic function measured
at a given time.

The first header row gives the evolution time at which the characteristic
function is measured.
The second header row denotes if the data is for the real or imaginary part of
the characteristic function.
The first column gives the value at which the characteristic function is
measured.
The remaining columns hold the data series.

### Asymmetric well expectation value data

File: `asymmetric_wells_expectation_x.csv`

Figure 5

This file contains five data series with one header row. Each series
corresponds to the expectation value of position as measured throughout a
simulation of a different anharmonic potential.
The position is given as dimensionless value, as described in the text.

The header row labels the case for the potential a-e, which corresponds to the
potentials described in the paper:
- (a) symmetric well
- (b) asymmetric well (phi=pi/20)
- (c) asymmetric well (phi=pi/20, start on opposite side)
- (d) asymmetric well (phi=pi/10)
- (e) asymmetric well (phi=pi/20, start in centre)


### Symmetric well characteristic function data

Files: all remaining files, labelled `cf_Xms_PPPP.csv`

`X` denotes the time of the evolution that is characterised in milliseconds.
`PPPP` denotes the real (`real`) or imaginary (`imag`) part of the
characteristic function.

Figure 8

These files contain the data for the characteristic function as measured in the
range -5 to +5 with 21 points in both the real and imaginary directions.
