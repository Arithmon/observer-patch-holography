# 22. Why Do All the Clocks Agree?

Put four clocks on a table and look at what each of them is.

The first is the one chapter twenty-one built: a region holds a state, the state determines a flow, and the flow's parameter is the region's own time. It has no dial. It is a family of transformations of the region's questions, labeled by a bare number. The second is a tally. Something happens at the region's edge and a counter goes up by one. The third is a ripple: a wave moving through whatever fills the region, whose phase you can read off with an interference probe. The fourth is light, or whatever plays its part, arriving from somewhere else after a delay that a ruler and a stopwatch would call a distance divided by a speed.

Physics writes one letter for all four. The *t* in the Schrödinger equation, the *t* in the wave equation, the *t* in the light cone and the *t* that thermal physics turns into an inverse temperature are typeset identically and nobody stops to ask why. In a theory that starts with a spacetime, the question never comes up, because the coordinate is handed out before the first equation and every mechanism is required to use it.

This world hands out no coordinate. Each of the four clocks has to be built from the one thing there is, a continuing history of repairs, and each could in principle have come out running at its own rate, in its own units, with no conversion between them. They do not. This chapter builds them one at a time out of the same history and shows that they read one parameter, with fixed conversion factors that the history itself supplies. Along the way it picks up a field, a wave, a distance, a minus sign and a quantum field theory, because those turn out to be what agreeing clocks are made of.

## Counting every touch

Start with the tally, since it is the easiest to build and the easiest to get wrong.

An observer region sits in the repair traffic. Seams around it are visited one at a time, each visit compares the two sides and commits the repair law's answer, and many visits change nothing at all, because the two sides agree. The obvious clock counts changes: tick whenever a unit actually moves across one of the region's seams. It is a bad clock.

Chapter twenty-one's timekeeping observer did something else. It stepped its flow at every repair event that met its region, including the ones that moved nothing across the seam. That looked like bookkeeping. It is the whole difference between a clock and a noise source. Write down how the region's state changes from one tick of the change-only clock to the next and the result has a mode that flips sign at every tick. A transformation with a negative mode has no logarithm a clock can use, which means no continuous family of positive steps runs through it: there is no way to say what half a tick of that clock would be. Count every touch instead, idle visits included, and the tick-to-tick transformation is positive throughout, its logarithm exists, and the clock can be read at any fraction of a tick.

Call the result an **incident clock**: a counter that advances once whenever a repair meets the region, whatever the repair does.

Different regions are touched at different rates. A single port sits on five or six seams. A whole twelve-port carrier is met by thirty-three. A patch of the screen is met by a number that grows with its area. So the raw counts disagree from the first visit. Divide each count by the number of touches its region expects per unit of the history's progress, and the disagreement goes away. The normalized count of every region converges to one parameter, uniformly over any stretch you care to fix, with the probability of missing by more than any given margin falling off as the history gets finer. A port and a carrier read the same time, because the visit schedule is the same schedule for both and each count is a fair sample of it.

The finite version of this was in chapter twenty-one: two clocks with equal average rates whose readings between ticks spread with a variance of sixteen. That spread is real. What the normalization does is make it small against the count as the count grows. A clock is a tally large enough that its fluctuations stop mattering, and the incident clock gets there without anybody synchronizing anything.

## A field out of the tally

Now look at what the tally is counting through.

Every seam visit leaves a trace in the occupations of the two ports it touched. Average those occupations over long stretches of the history and the averages fluctuate, slightly, around the settled value. The fluctuations are not independent from port to port. A unit that left one port has gone somewhere nearby, so neighboring averages are correlated, and the pattern of correlations is a field: a number at every place, with a definite rule for how numbers at different places hang together.

The rule is one operator. Call it A. It is built from the seams, it is positive, and it is local in the strict sense, each port talking only to the ports it shares a seam with. On the screen it has the heat-trace dimension of a surface, two, as it should for a surface. Its inverse is the field's correlation between two places, and as the screen is refined, the correlations on the finer screens converge onto a single continuum operator with a definite domain. Nothing was put in by hand. The occupations were always there; the operator is what their long-run statistics turn out to be.

Then look at the field at two different times. The correlation between the field read at one clock reading and the field read at another decays exponentially in the difference between the readings, at rates set by A, mode by mode. The slow modes forget slowly and the fast ones fast. The clock in that statement is the normalized incident count from the last section. The field's time is the tally's time: the same parameter that counts touches is the parameter in which the field relaxes, with no conversion factor between them other than the normalization already fixed.

That is two clocks agreeing. It is also the first sign of something larger, because A has now appeared twice, once as the field's shape in space and once as the field's rate of forgetting in time.

## A wave on the same count

A field that only relaxes is a heat-like thing. It smooths and forgets. Waves keep their phase, and keeping phase needs something the passive repair traffic does not do.

So give the observer a small quantum memory of its own, kept out of the repair traffic, and couple it to the incident clock in the simplest possible way. Each tick, idle touches included, triggers one fixed step of a walk on that memory. Nothing else touches it between ticks. After N ticks the memory has taken exactly N steps.

The walk is built from the same seams as the field, and its amplitude follows a cosine. In the limit of a fine screen the amplitude is the cosine of the clock reading times the square root of A, with one fixed conversion factor set by how many seams the clock's port sits on. The operator that set the field's shape and its forgetting sets the wave's frequencies. The clock that counted touches drives the oscillation. There is no second clock inside the wave; the wave's time is the tally, read through a cosine.

Because everything here is a finite construction with exact rational odds, it makes predictions down to the last digit. Prepare the memory, run the clock, and read it with an interference probe. The probability of the plus outcome is 1 before any tick, 21/22 after one and 91/110 after two. Over a complete two-visit history, every branch included, it comes to 1297/1320. Change the rule so that only visits that actually move something trigger a step, and the same history gives 7901/7920 instead. Reprepare the memory before every step and it gives 71/72. Those are different numbers from the same repair history, and the difference is how an experimenter tells which clock a region is actually keeping.

## Distance is the first answer

Two regions, some way apart on the screen. Disturb the first and watch the second through the wave.

For a while nothing arrives. The wave's response at the second region is exactly zero, with no small tail leaking ahead, until the clock reaches a definite value. After that the response switches on. The value at which it switches on is the distance between the two regions, measured with the field's own operator as a ruler: the largest amount a function can change from one region to the other while its slope, as A measures slope, nowhere exceeds one.

The relaxing field gives the same number by a different route. Its leakage from one region to the other over a short clock interval falls off like a Gaussian in that same distance, and the exact rate of the fall-off reads the distance back out. The bound is sharp. Nothing arrives before the distance, and something does arrive as soon as the clock passes it.

So distance and duration come out of one apparatus in one unit. A distance is how long the clock runs before a disturbance can be answered, and that makes the speed of the answer one, by construction, for every pair of regions and every observer keeping this clock. Chapter seventeen found one invariant speed from the angle-preserving maps of the sky and noted that the number it takes is a choice of units. This is the same speed from the other end. The cone is the first response of the wave, and the speed of light is the exchange rate between the count the tally keeps and the ruler the field supplies. The meter and the second were tied together in 1983 for exactly this reason: two units for one thing.

## Three directions and one minus sign

The screen supplies two directions. The third comes from the records themselves.

Compare the field with itself some number of repairs later. Those comparisons are lag records, native readouts like any other, and a weighted mixture of them, the longer lags weighted more heavily, gives the field a new coordinate running away from the screen: the deeper the coordinate, the longer the lags it is built from. Independent long stretches of one history then fill that coordinate in, so the field along it fluctuates on its own instead of copying the screen. Glue the radial coordinate onto the screen's field and the joint operator is the screen's A plus the second derivative along the new direction. Its heat-trace dimension is three. The memory of a region sits on its surface, as chapter twenty-four will show, and the depth of that memory in the history is the third direction of the space it describes.

Run the wave in all three directions. Each tick now triggers two fixed pulses, one carrying the screen's frequency and one carrying the radial frequency, and the two pulses act along axes that do not commute. That one fact fixes how the frequencies combine. The joint frequency is the square root of the sum of their squares, the way the hypotenuse of a right triangle combines its sides. Use axes that commute and the frequencies simply add, which is the wrong answer and fails the direct check. The square-root rule is the dispersion law of a wave moving at one speed in three directions, and it has fallen out of the order in which two pulses act on a memory.

Write the wave's action, the quantity whose stationary paths are the wave's motions. The three space directions enter with the coefficients of A and the radial piece, all positive. The clock direction enters with the opposite sign. That sign is no convention. The wave is second order in the clock and its spatial operator is positive, and an action with that structure assigns the clock direction a minus sign against the three others. The coefficients of that action are a metric, and its signature is minus, plus, plus, plus. That is the signature of the spacetime of special relativity, recovered from a tally, a positive operator and a cosine.

The native metric is measurable rather than smooth: it is defined almost everywhere, with coefficients pinned between two positive bounds, and it carries a crease where the original faces of the screen meet. Laboratory physics asks for smooth spacetimes and it gets one. Pick any finite list of measurements, wave readings, clock comparisons, interventions, detector outputs, and any tolerance however small. One smooth static spacetime with one time and three space directions reproduces every joint probability on that list within that tolerance. No finite set of observations tells the crease apart from a smooth spacetime.

## The region's own clock, on the wall

Now go back to the first clock on the table, the flow from chapter twenty-one.

Take one mode of the wave and give it the energy the clock assigns: two levels, one frequency apart. Put the mode into a state and ask the chapter twenty-one question, which flow that state determines. Then ask when that flow is the same motion as the wave's own evolution, run in clock time.

The answer is a theorem with no slack in it. The state's own flow equals the clock evolution, rescaled by a constant, exactly when the state is thermal for the clock's energy: odds dropping off as the exponential of the energy times that constant. The constant is the inverse temperature. The match has to hold axis by axis, so a state that leans the right amount in the wrong direction fails, and two modes leaning by the same amount at different frequencies would need two different constants and so fail to share one. A region in equilibrium at one temperature has one conversion factor for every mode, and that factor is its inverse temperature, which is what chapter twenty-one said a thermometer reads.

That puts four readings on one parameter: a tally that counts touches, a field relaxing in its count, a wave oscillating in its count, and the flow of the region's own state, with conversion factors fixed by the seams and the temperature. The region's own time is the time on the wall, multiplied by the temperature.

## Quanta

One mode in one memory is a single quantum system. A field has many quanta. Take many independent copies of the thermal mode and add up their fluctuations, scaled the way the central limit theorem scales them.

Two of the summed quantities come out special. Their commutator has expectation value exactly i, the defining relation of a position and momentum pair, and the products of their exponentials converge to the expectation values of a thermal free quantum field. The two-point function at two clock readings has a real part that oscillates as the cosine of the frequency times the difference in readings, and an imaginary part that oscillates as the sine. The imaginary part is the commutator. It is the piece that has no classical counterpart, and it runs in the same clock as everything else on the table.

Put the screen modes and the radial modes in together and the joint thermal state of the three-dimensional wave has a modular flow, in chapter twenty-one's sense, that is time translation of that same wave, rescaled by its temperature. Bisognano and Wichmann took a quantum field and its vacuum as given and proved that a region's modular flow is a motion. Chapter twenty-one used their result about an accelerating observer to show that the region's flow is a physical motion. That result needed a quantum field theory to stand on. The field theory is here, built from the repair history and running on the observer's own count.

It does not stop at free fields. Add an interaction, a quartic term coupling the radial modes and the screen's constant mode, and the interacting thermal state can be prepared from native operations, its flow fixes an ideal clock, and a native incident clock driving the same interacting system runs in step with that flow. The particles and forces of part six are interacting fields of this kind.

## Seconds

None of the four clocks has seconds in it. The tally counts touches, the wave oscillates in that count, the region's flow runs in that count multiplied by a temperature, and light covers one unit of the field's distance per unit of count. Every conversion between the four is fixed inside the construction, which is the point of the chapter. The one conversion left is to the clock in a metrology laboratory, and it is a single number: how many touches of some region make one period of the cesium transition chapter twenty-one opened with. Fix that once, and every other duration in this world is quoted in seconds by arithmetic.

Physics wrote one letter for all four clocks and never had to justify it, because its spacetime came with a coordinate and its equations were required to use it. Here the coordinate is an output. It is the one parameter that a tally, a field, a wave, a light cone and a region's own state all turn out to share, and the reason they share it is that all five are readings of the same continuing repair history, taken by observers who count everything that touches them. Appendix D shows how to run the construction and check the numbers in this chapter on a laptop.
