# 21. Why Does Time Feel Like It's Moving?

Augustine was in his forties and a year or two into the bishopric of Hippo Regius, a port on the North African coast that is Annaba in Algeria today, when he wrote the eleventh book of the *Confessions*. The thirteen books took him from about 397 to about 400. Twelve of them are about his life and his God. The eleventh is about time, and in its fourteenth chapter it comes apart in his hands.

"What then is time? If no one asks me, I know; if I want to explain it to a questioner, I do not know."

The line gets quoted as a piece of wisdom, which flatters it. Augustine was after an instrument. He works through the possibilities in the chapters on either side of that sentence and discards each one: the past has gone and cannot be laid alongside a ruler, the future has not arrived and cannot either, and the present has no extent, because anything with extent divides into a part that has gone and a part that has not come. He ends the chapter holding a ruler with nothing in front of it.

What he settles on, having run out of things in the world to measure, is that the measuring happens in the mind: the soul stretched across a memory, an attention and an expectation, so that the quantity being compared against a ruler is an impression rather than an interval of anything. He gives it a name, *distentio animi*. The name has kept philosophers busy for sixteen hundred years. It moves the ruler indoors and leaves the question of what is being laid against it exactly where it was.

Sixteen centuries of instrument-building have not improved his position. In 1967 the thirteenth General Conference on Weights and Measures redefined the second as the duration of 9,192,631,770 periods of the radiation from a particular transition in a cesium-133 atom, which is the most carefully specified duration in the history of the species and is a duration. A cesium clock counts periods and reports a count. Subtract two counts and you have an interval. Every timepiece anybody has built, from a notched candle to the fountain clocks that keep international atomic time to sixteen digits, reports a difference between two events, and none of them has an output for the present moment, because there is no part of the mechanism that could produce one.

Chapter sixteen separated order, direction and duration. Informational order records which commit certified reading a version written by another. The wiring constrains possible conflicts; authenticated provenance fixes the actual edges. Averaging seven inputs into one illustrates a loss from a reduced seam description, while retaining an input record can preserve the distinction elsewhere. Neither this arithmetic nor the dependency graph gives an elapsed time. A clock requires an additional dynamical rate and a calibration.

## An algebra and a state

An observer holds an algebra of questions, which is everything it can ask about its own patch, carrying the brackets that record which pairs can be asked together and which cannot. And it holds a state, which is the odds it puts on every one of them.

That pair fixed the floor under position and momentum, and put the ceiling on separated agreement at 2√2 and not a thousandth higher. It also decided which part of the world is classical, by collecting the questions that commute with everything into a center and letting the records live there.

Take the pair apart looking for a clock.

The algebra lists questions and how they compose. The state assigns their probabilities. Neither object comes with seconds written on it. Together they can define a modular flow, a mathematical family of reversible changes. Interpreting its parameter as a clock reading requires a physical dynamics and an instrument.

Chapter ten’s run makes the distinction concrete. Eighty-two thousand patches worked all 122,880 seams into agreement, but the reported cycles came from the executor. Such bookkeeping does not identify a physical tick. A clock needs a changing degree of freedom, a readable record and a calibration.

## A room left alone

Leave a room alone for long enough and everything in it arrives at one temperature. The coffee, the mug, the table, the air above the table. The description of that room has collapsed to a single number. The collapse discards almost everything that was ever true about it.

A temperature helps determine the direction of heat exchange when two systems meet. Its rate also depends on their materials and contact. The equilibrium description and the dynamics leading to it are different pieces of information.

Run that argument backwards.

Chapter sixteen supplied a reference law preserved by a specified resampling channel. Invariance names an equilibrium only relative to that dynamics. For the modular construction, the question is which reversible flow has the chosen state as its normalized thermal reference.

So work out what it would have to be.

A change that qualifies is a reshuffling of the observer's questions among themselves that respects the algebra: questions go to questions, sums to sums, products to products, and the whole thing is reversible. Line those shuffles up, one for every real number, and the family is what chapter eighteen called a flow, driven by the move it repeats, its generator.

Under the modular theorem’s algebra and faithfulness hypotheses, two demands characterize the normalized flow.

The first is that the state's odds are unchanged by it. Every question is worth after the shuffle exactly what it was worth before. That is what being settled means, written out.

The second is where the non-commuting part earns its keep. In a world where every pair of questions could be answered together, the state's weight for asking A and then B would equal its weight for asking B and then A, because there would be nothing to tell the two orders apart. Chapter nineteen's entire subject is that they differ. So the state carries a gap between the two orders, one gap for every pair of questions. The gap is a fixed feature of the state rather than a nuisance to be argued away. The second demand is that swapping the order costs exactly one unit of the shuffle, the same one unit for every pair of questions, with the amount of shuffling measured along a direction at right angles to the flow's own parameter. That last clause is complex analysis, the one step in the whole construction that cannot be done on a table.

For a faithful state in the modular theorem’s setting, those demands select a unique flow. In the finite matrix example, faithfulness means that every eigenvalue of the density matrix is positive. The resulting parameter is modular time. A nontrivial flow, a suitable prepared readout and a physical calibration are further ingredients of an operational clock.

Minoru Tomita, who lost his hearing at the age of two and published very little in a long career, wrote the construction down in 1967 in a manuscript that the people who received it found close to unreadable. Masamichi Takesaki worked out what it said, passed the results to Jacques Dixmier in the summer of that year, and published a usable version in 1970 as a thin volume of Springer lecture notes. The thing is called **modular flow**, after the algebraic machinery it drops out of rather than after anything it does, and its generator, in chapter eighteen's sense of the word, is called the modular Hamiltonian.

In the finite matrix case the recipe is short: take the matrix logarithm of the strictly positive density matrix and change its sign. This dimensionless operator generates modular flow in a specified sign convention.

$$K = -\log \rho$$

Here rho is the density matrix and K the dimensionless modular Hamiltonian. The logarithm is natural: a doubling contributes 0.693. To apply it, first diagonalize rho, take the logarithm of each positive eigenvalue, and transform back. The formula is a finite matrix realization; general operator algebras use the modular operator.

## One system, two outcomes

That recipe can be run by hand on the smallest system there is.

Take a two-level system with supplied energies zero and E. Put it in a Gibbs state at a supplied temperature T. The upper-to-lower probability ratio is the exponential of minus E divided by T when temperature is written in energy units. These energy and temperature inputs give the state its thermal interpretation.

Write the state as its array. Take the logarithm, which acts on the two diagonal entries separately. Put a minus sign in front. What comes out has E divided by T in the upper slot, zero in the lower, and a constant added to both, contributed by the requirement that the odds add up to one.

The constant does nothing whatsoever. A generator shifted by a constant produces the same flow, since the constant piece commutes with the whole algebra and shuffles nothing. Delete it.

What survives is the energy operator divided by the temperature in energy units. Its modular flow follows the same operator orbits as the supplied Hamiltonian evolution, with the parameter scale and orientation fixed by the temperature, time unit and sign convention.

Put numbers on this diagonal thermal state. Let the lower level have probability 20 in 21 and the upper level 1 in 21. The ratio determines the two diagonal probabilities. Off-diagonal coherences are absent on this specified Gibbs branch.

The lower slot of the generator holds zero, so the whole generator is whatever stands in the upper slot. That is the logarithm of the ratio between the two odds.

$$K = \log \frac{\text{odds of the lower answer}}{\text{odds of the upper answer}} = \log 20 = 2.996$$

Add it up yourself. Twenty is two times ten, so its logarithm is the sum of theirs: write 0.693 for the two, write 2.303 for the ten, and the twenty comes to 2.996, four thousandths short of 3. The check runs from the other end, since the exponent 3 gives 20.0855, a shade over the twenty on the slip.

The dimensionless generator gap is 2.996. On the supplied thermal branch, that equals the energy gap divided by temperature in energy units. It does not determine either quantity separately.

The slip of probabilities therefore determines a modular Hamiltonian. To identify a physical Hamiltonian H, supply a positive inverse temperature beta and the Gibbs relation: K equals beta times H plus log Z times the identity. Scaling H and the temperature by the same positive factor preserves the slip.

## Whose clock

The flow belongs to the algebra–state pair.

Different faithful states can define different modular flows, although some distinct states give the same flow. A tracial state, which weights a finite matrix block uniformly, gives a trivial flow even on a noncommutative algebra. Comparing two physical clocks additionally requires their dynamics, readouts and a relation between their time units.

The two-level example shows the parameter dependence. Keep H fixed and halve the temperature. The probability ratio squares from twenty to one to four hundred to one, and the modular generator gap doubles from 2.996 to 5.99. The two modular parameters then use different scales for the same Hamiltonian evolution. A prediction about two clock instruments requires their physical readout models.

Chapter eighteen approached energy through generators. Here the state supplies a dimensionless modular generator. Joules enter through the separate energy identification and calibration, rather than through a count of unsettled records.

The state dependence here is precise: a normalized modular flow is determined by the specified algebra and faithful state. That mathematical fact does not by itself identify a subjective experience, a laboratory duration or a rate in seconds.

## The center does not move

Records live in the center of the algebra: the questions that commute with everything, which is what makes a record readable without disturbance, copyable to a neighbor, and the same on the second reading as on the first. Run the flow on the center and it does nothing at all. Every element of the center is left exactly where it is by every transformation in the family, at every value of the parameter.

The flow was built to close the gap between asking A and then B and asking B and then A. An element of the center has no gap with anything, because it commutes with everything in the algebra. So there is nothing for the flow to fix about it, and it sits untouched. The price of commuting with everything is having no clock of your own.

Chapter nineteen proved the harsher version, which reaches past this flow to every continuous motion a record layer could have had. The only structure-preserving moves a record algebra has are permutations of its finitely many labels, and a continuous family of permutations would have to jump to get anywhere, so it sits where it started. The classical layer changes one discrete relabeling at a time, or it holds.

These facts locate possible modular motion outside the center. They do not prove a mechanism for experienced time. Even there, a tracial state has trivial flow. To build a clock, one must specify a nontrivial dynamics, a preparation and readout that reveal change, and a physical time scale.

## The flow has no arrow

The flow is reversible. Every transformation in the family has its inverse in the family, which is the transformation at minus the same parameter. That was settled before either demand was written down: what qualified as a change at all was a reshuffling of the questions that loses nothing, and a shuffle that loses nothing can be undone.

This reversible flow does not itself supply a dissipative arrow. Chapter sixteen’s averaging example loses distinctions from the seam description, while a retained record can preserve them elsewhere. The conditional thermal resampling construction gives a separate entropy-production inequality under its stated reference and energy assumptions.

The constructions answer different questions. An authenticated history gives event dependencies. A specified dissipative channel gives an entropy-production direction, while reversible dynamics gives a rate in its own parameter. A physical clock combines a repeatable dynamical process, readable records and a calibration. Recording its readings does not make every commit a logically irreversible erasure.

## The region's own clock

Chapter eighteen defined energy and left an unpaid bill inside the definition.

A region’s dynamics can have an energy generator once its time parameter and energy convention are fixed. A repair-event count alone supplies neither that generator nor an elapsed duration.

For a finite region with a faithful density matrix, the modular construction supplies K = −log rho and its dimensionless parameter. A physical clock adds a changing preparation and readable records. Comparing rates between regions requires a common physical interpretation of those readouts.

On a supplied Gibbs branch, K = beta H + log Z times the identity, with beta = 1/(k_B T) for T in kelvin. The constant disappears from the flow. Identifying H in joules and relating its evolution to seconds require the energy, temperature and time dictionary. Modular reconstruction establishes the relation once those inputs are supplied.

## The ship that light never catches

Now supply a Minkowski quantum field theory in its vacuum state, a wedge algebra satisfying the Bisognano–Wichmann hypotheses, and an ideal observer with constant proper acceleration. This is a physical continuum example of the modular construction.

Chapter seventeen constructed a Lorentz kinematic carrier, with three rotation and three boost directions. In the present supplied spacetime, boosts act geometrically on a wedge, and a uniformly accelerating trajectory follows one boost orbit.

In Minkowski spacetime, an eternally uniformly accelerating observer has a causal horizon. Signals from beyond the relevant horizon never reach that trajectory. This uses the spacetime’s physical causal structure; the finite kinematic carrier by itself does not identify an observer’s physical light cone.

Use the quantum field theory’s wedge algebra and the vacuum restricted to it. Locality, the vacuum and the hypotheses linking the algebra to spacetime are supplied here. The geometric boundary alone does not construct that field theory or state.

For this vacuum–wedge pair, the Bisognano–Wichmann theorem identifies modular flow with Lorentz boosts at the fixed normalization. Along a specified uniformly accelerating orbit, boost parameter is proportional to proper time, with the conversion set by the acceleration and physical units.

Joseph Bisognano and Eyvind Wichmann proved it for quantum field theory in two papers in the *Journal of Mathematical Physics*, "On the duality condition for a Hermitian scalar field" in volume 16 in 1975 and "On the duality condition for quantum fields" in volume 17 in 1976. They were not working on a theory of time. The question in front of them was technical, whether the operators available inside a region exhaust the operators that commute with everything outside it. The flow of the region turned out to be a Lorentz transformation.

Here the algebraic and geometric flows can be compared because the field theory, vacuum, wedge and proper-time convention have all been specified. Their agreement gives a concrete physical realization of modular flow under those hypotheses.

## Warm

The physical identification also gives a thermal detector relation.

The vacuum restricted to the wedge satisfies the thermal KMS condition for boosts. Converting boost parameter to the uniformly accelerated observer’s proper time gives the Unruh temperature. A detector reads that relation through its coupling to the field; it is not a statement that every thermometer or finite observation behaves as an equilibrium bath.

Bill Unruh calculated an ideal detector response. In the stationary, long-duration limit, a ground-state inertial detector in the Minkowski vacuum has zero excitation rate. A uniformly accelerated detector has a thermal excitation-to-de-excitation ratio at the Unruh temperature. Finite switching and the chosen coupling require their own response calculation.

He did that calculation in 1976, in a paper called "Notes on black-hole evaporation" in *Physical Review D*. The temperature is proportional to the acceleration. The constant of proportionality is about four parts in ten to the twenty-first of a degree for every meter per second squared.

Which is why nobody has run into this by accident. Warming yourself by a single degree calls for an acceleration around twenty-five billion billion times the strength of gravity at the Earth's surface, sustained, with you in the vehicle. Every acceleration a human body has ever survived leaves the reading far below what any thermometer built resolves, and the passengers have complaints about the other effects.

The algebra and faithful state define a modular parameter and a normalized KMS relation. The supplied Minkowski field theory and vacuum identify that flow with boosts; acceleration and the physical units convert it to proper time and temperature. The laboratory readings depend on these additional identifications.

In this physical example, duration and temperature are related through one specified dynamics. The general modular theorem supplies the mathematical relation, while a clock and a thermal instrument supply its calibrated readings. What does the temperature reading measure?
