# 25. Why Does Gravity Look Like Geometry?

The crew of the International Space Station float about four hundred kilometers above the ground, where the Earth's gravity is close to nine tenths as strong as it is at the surface. Gravity has kept hold of them. The station and everybody aboard are falling around the planet together, and that is the whole of the trick.

Albert Einstein worked this out in the autumn of 1907. He was examining patents in Bern and owed Johannes Stark a review article on relativity for a yearbook of radioactivity and electronics, and he was assembling it at his desk in the office when the thought arrived that he described fifteen years afterward, in a lecture in Kyoto on 14 December 1922, as the happiest of his life. If a person falls freely, he will not feel his own weight.

Every other influence in physics discriminates. A magnet lifts a steel paperclip and does nothing at all to an aluminum one of the same mass, because the response depends on what the object is made of. Gravity pulls harder on a heavier object by exactly the factor that makes a heavier object harder to shift, the two cancel, and everything falls at one rate, which is why David Scott could stand on the Moon in August 1971 with a hammer in one hand and a falcon feather in the other and drop them together for the camera. Write the description in a frame that falls along with everything else and gravity drops out of it.

The sharpest test of that cancellation flew between 2016 and 2018 in a French satellite called MICROSCOPE, which carried two hollow concentric cylinders, one of a titanium alloy and one of a platinum alloy, around the Earth in an orbit continuously corrected so that nothing but gravity acted on the pair. Electrostatic forces held both cylinders in place. The size of those forces reports any difference in how the two of them fall. Five months of usable free-fall data came out of two and a half years in orbit. The difference in acceleration between titanium and platinum came out at one and a half parts in ten to the fifteenth, with a statistical uncertainty of two and three tenths in the same units: a null result at the level of one part in a thousand million million.

So gravity belongs to the place rather than to the thing. Which leaves two questions. Why does an influence that treats everything identically turn out to be a shape, and what decides how much shape a given lump of matter makes?

## Two walkers on the equator

Two people stand on the equator, a kilometer apart. Both face due north. Both are told to walk straight ahead and never turn, with no course corrections and no navigating.

They collide at the North Pole, about ten thousand kilometers later.

Both of them walked straight the whole way. The gap between them closed from a kilometer to nothing, at an accelerating rate over the last stretch, and no force acted on either walker at any point in the trip. A cord tied between them would have gone slack rather than taut. The convergence is a fact about the surface they crossed. They can establish it without ever seeing the Earth from outside, with nothing more than a tape measure between them. Send a third walker off from a point between them and the same thing happens to each of the three gaps at once. That is what makes the effect a property of the ground rather than an agreement between two people.

The path you get by going straight ahead is called a **geodesic**. On a flat floor it is a straight line, on a sphere it is a great circle, which is why the flight from London to Tokyo goes over the Arctic and why that route looks bent on a wall map. Two geodesics that set off parallel and end up crossing is what curvature means. The walkers measure it from the inside.

Free fall is walking straight ahead. Two ball bearings released side by side above the Earth converge, and the amount they converge by is a fact about where they are rather than about what they are made of, which is what MICROSCOPE measured to one part in a thousand million million. Nothing in the account calls for a sheet with a weight on it, and nothing bends into a further dimension for the bending to happen in.

Release a small ball of particles, all at rest with respect to each other, and let them fall. The ball's volume changes. John Baez and Emory Bunn, in the *American Journal of Physics* in 2005, wrote Einstein's equation as a sentence about that ball: the rate at which it begins to shrink is proportional to its volume times the energy density at its center, plus the pressure in each of the three directions at that point.

## The equation, read rather than solved

The whole of that, in the notation the textbooks use, takes one line.

$$G_{ab} + \Lambda\, g_{ab} = 8\pi G\, T_{ab}$$

That says the shape of the world at a point is fixed by what is sitting at that point.

The letters $a$ and $b$ each stand for one of the four directions, one of time and three of space. Four directions taken two at a time gives sixteen pairs, and swapping the two directions of a pair changes nothing, leaving ten of the sixteen independent, so that single line is ten equations. The object $g_{ab}$ is the **metric**: ten numbers at every point of the world that turn a pair of directions into their contribution to chapter seventeen's interval, which makes the metric the ruler, the thing that says how much time and how much distance separate two nearby events. The object $G_{ab}$ is a particular combination of the rates at which those ten numbers change from point to point, assembled so that its own bookkeeping closes automatically, with nothing leaking in or out. On the right, $T_{ab}$ is the stress: energy density, momentum flow and pressure, everything that is there. And $8\pi G$ carries Newton's gravitational constant, which fixes how much bending a given amount of energy buys.

The one term that spoils the arrangement is $\Lambda$, the cosmological constant, which multiplies the metric itself and has nothing on the other side of the equals sign to account for it.

Einstein reached the equation without that term in November 1915, by requiring that it reduce to the old law of gravity for weak fields, that it not depend on the coordinates anybody chose, and that both sides conserve what has to be conserved. Those requirements narrow the possibilities to one. They do not say why geometry should respond to energy in the first place.

## Four pages in 1995

Ted Jacobson's paper occupies pages 1260 to 1263 of volume 75 of *Physical Review Letters*, titled "Thermodynamics of Spacetime: The Einstein Equation of State".

He needed three things. An observer who accelerates steadily has a horizon behind it, a boundary past which signals never arrive. A horizon carries an entropy equal to a quarter of its area. And an accelerating observer finds empty space warm, at a temperature set by its acceleration. Then Clausius's relation, which is the first thing anybody learns about heat: the heat flowing into a system, divided by the temperature it arrives at, equals the entropy that heat buys.

Jacobson's move was to demand that this hold on every one of those horizons. A black hole has one. So does the air above your kitchen table, and so does every other point of the world, one horizon for every direction an observer passing through that point might accelerate in.

Matter crossing such a horizon carries energy, so it delivers heat. The horizon's entropy is a quarter of its area, so the entropy change is an area change. The temperature is the one the accelerating observer reads. Insist that heat over temperature equals the area change at every point and in every direction at once, and the geometry has no freedom left: it has to distort in exactly the way that keeps the accounts straight everywhere, and that distortion is Einstein's equation, coupling constant included.

Jacobson had to bring the area law and the temperature into the argument. OPH offers candidates for both. Chapter twenty-four counts the capacity of patches on a screen. Identifying that capacity with physical boundary entropy, with the area normalization gravity requires, is a further step. Chapter twenty-one relates modular flow to an accelerating observer under its geometric-flow and temperature-normalization hypotheses. A finite state has a modular flow; establishing that this flow is the physical boost is part of those hypotheses.

The finite repair constructions establish a useful part of the argument. Under the declared schedules and conserved constraints, they converge to states of maximum entropy, whose first-order entropy change vanishes along variations preserving those constraints. Gravity needs a different, physical stationarity statement: the entropy of a small region's contents plus its boundary area divided by four times Newton's constant must be stationary under the admitted state and shape variations at fixed physical volume. A fixed occupation total or a fixed list of patches does not establish that condition. With the area, modular-energy and continuum identifications supplied on the same construction, the following calculation shows what this additional condition requires.

## Where the eight and the pi come from

Take a small ball around a point in the assumed smooth geometry, in the rest frame of one observer, and apply that fixed-volume stationarity condition in the shrinking-ball limit.

The contents side first. The required ball-modular identification weights the physical energy density by a factor that falls from a maximum at the center to zero at the edge: the difference between the squared radius of the ball and the squared distance from the center, divided by twice the radius. This geometric weight is an input to the small-ball theorem, rather than a consequence of having a finite Gibbs state. Integrating it over the ball gives four pi times the fourth power of the radius, divided by fifteen. The integral takes one line, and the fifteen comes out of it.

One factor remains: the geometric modular normalization of two pi. Under the boost and Unruh-temperature identification used in chapter twenty-one, an accelerated observer reads its acceleration divided by two pi, in units where the speed of light, Planck's constant and Boltzmann's constant are all one. The corresponding entropy-energy relation supplies the extra two pi, taking the matter side to eight pi squared times the fourth power of the radius, divided by fifteen. This normalization must belong to the same physical clock and stress as the ball calculation.

The geometry side. Alfred Gray and Lieven Vanhecke, in *Acta Mathematica* in 1979, worked out how the size of a small geodesic sphere depends on the curvature around it. The piece of that needed here says that if you hold the volume of a small ball fixed and raise the curvature, the area of its boundary falls by four pi times the fourth power of the radius, divided by fifteen, times the change in exactly the curvature combination that stands on the left of Einstein's equation. Curving the surroundings costs boundary area.

Set the sum of the two to zero, remembering that the area enters the total entropy divided by four times Newton's constant. At leading order every appearance of the radius cancels. Both fifteens cancel. What is left is eight pi squared on one side against pi over Newton's constant on the other, and their ratio is eight pi times Newton's constant, the coupling in the field equation. Under these entropy and geometric normalizations, the factor eight pi follows without a fit. The calculation does not itself select Newton's constant or establish the physical inputs that gave the two coefficients.

That gives a relation in one observer's rest frame. Suppose the same geometry and stress tensors obey it for observers passing through the same event at every velocity in an open range, with the transformations and remainder bounds required by the theorem. Written out in terms of the velocity, each relation is a polynomial of degree two. A polynomial that vanishes throughout an open range has every coefficient equal to zero. The tensor upgrade recovers all ten components once those common physical inputs are established.

## A geometry that can respond

A wave description on a curved space does not by itself make gravity dynamical. The conditional field construction supplies a useful test. Its static product geometry has no radial null curvature response. The scalar field it describes admits a radial wave carrying positive null stress. Keep the reference and other sources fixed, and the same geometry cannot obey Einstein's balance for both the reference and that excitation.

The test identifies a required response. If the area of a spherical slice is four pi times a radius squared, energy passing along an affinely measured light ray must change the second derivative of that radius by minus four pi times Newton's constant, the radius and the null stress. This is the focusing part of Einstein's equation. Deriving that response from the same repair history is a physical obligation. The entropy and tensor theorems state what follows when their premises are met; the finite field description does not meet them simply by possessing a Lorentzian metric.

## A gas law without molecules

Jacobson ended his abstract with a sentence that has been argued over for thirty years: it may be no more appropriate to canonically quantize the Einstein equation than it would be to quantize the wave equation for sound in air.

Pressure times volume equals a constant times temperature. Émile Clapeyron put it in that form in 1834. It held up in every laboratory in Europe through the rest of the century without containing a single molecule. In 1908 Jean Perrin put grains of gamboge resin in water under a microscope, counted how the population thinned out with height, and got Avogadro's number out of the count, which finished an argument about whether matter came in pieces at all. Wilhelm Ostwald, who had spent a career arguing that atoms were a bookkeeping convenience, gave that position up in the preface to the fourth edition of his own textbook in 1909, citing Perrin by name. The gas law did not change a symbol, because it was never a claim about molecules. It is a claim about a gas that has settled down, and stays true whatever the gas turns out to be made of.

The conditional entropy argument puts Einstein's equation on that shelf. A network satisfying its physical entropy, stress and continuum premises must report that equation. Settling alone does not supply the premises. In this interpretation, quantizing the metric in the hope of finding the atom of space resembles quantizing the pressure of a gas in the hope of finding the molecule. The pressure is real enough to lift a piston, and it belongs to a settled state rather than to a constituent.

Once the smooth metric and physical frame identification are in place, the 1907 thought has its usual precise meaning. A freely falling frame removes the connection at a point, while tidal curvature remains. The finite frame constructions supply algebraic ingredients for this identification; they do not establish physical free fall by themselves. For pressureless matter with nonzero density, covariant stress conservation then makes the flow geodesic. The walkers return as a consequence within the physical geometric description.

## Nine light rays

Everything in the argument crossed a horizon, and horizons are made of light rays.

The horizon probes used in this argument evaluate stress along a light-ray direction. The stress takes two directions and returns a number. Feed it the same null direction twice and you get the number used in the null-balance theorem.

This restriction specifies the data used by the theorem. Light cones define the causal boundary of the physical description, but not every local measurement is a null contraction of stress. A quantity invisible to these horizon probes is undetermined by these probes; a different measurement or an additional physical condition can carry more information.

The conditional scalar-field construction supplies a compatible cone test. Disturb its field in one slab and read a detector in a separate slab. Under the construction's refinement and detector assumptions, the response before that field's light cone arrives shrinks to zero.

Ten numbers describe the stress at a point. Light rays report one number each. How much of the stress do the rays determine?

A symmetric object that reads zero along every light ray has to be a multiple of the metric, and nothing else. So the rays pin down nine of the ten numbers. What they leave undetermined is a multiple of the metric.

For readings from one symmetric tensor, nine directions do it. Take the six rays that travel along one axis at a time, forward and backward along each of the three. Then take three more that share their motion equally among all three axes, each space component equal in size to one over the square root of three, with the signs plus plus plus, plus plus minus, and plus minus plus. The six axis rays never pick up the cross terms, the entries pairing one space direction with a different one, because a ray running along a single axis has no second direction to pair it with. The three tilted rays are there for those, and their three sign patterns read the three cross terms off one another. Nine light rays, nine numbers, one system of linear equations to invert. Its determinant is exactly eight thousand one hundred and ninety-two over twenty-seven, or two to the thirteenth over three cubed. A determinant that misses zero is a system that inverts, and this particular inversion is stable: small errors in the nine readings stay small in the nine numbers recovered from them.

Any nine numbers fit that system, so the fit alone does not show that the readings come from a tensor: a tenth ray, or a hundredth, has to land where the recovered tensor says it will. Given the geometric and stress tensors, balance along the nine rays fixes the field equation up to a multiple of the metric. Conservation makes that multiple constant on each connected region. That term is the tenth number.

## The vacuum and the tenth number

A Lorentz-invariant vacuum contribution to stress is a number times the metric. This is the vacuum-reference condition used here; a general quantum state on a curved background need not have stress of that form.

Contract it with a light ray. The metric evaluated on a light-ray direction twice is the interval of a light ray with itself. That is zero: chapter seventeen defined a light ray as a separation whose interval is exactly zero. So the vacuum's contribution to what a light ray reports is a number times zero.

Exactly zero, for any vacuum energy density whatever, however enormous, and zero by the definition of the probe rather than by a cancellation between large terms.

The single direction the light rays cannot see is exactly the direction the vacuum's energy points in. That is why the cosmological constant survives in the field equation as a term nobody derived: the small-ball calculation that produced everything else on the left-hand side is built out of horizons, horizons are made of rays, and rays are blind to it.

A naive Planck-scale cutoff on quantum zero-point energies gives a vacuum energy density around a hundred and twenty orders of magnitude larger than the one inferred from the sky, whose cosmological constant is about 1.1 times ten to the minus fifty-second per square meter. The null-balance theorem explains a precise limitation: its data cannot fix the metric-proportional term. That limitation does not by itself solve the vacuum-energy problem or prohibit other local measurements and calculations.

A specified horizon-capacity condition offers a way to fix it. On that branch, the smallness of the constant is the largeness of the horizon's record count. Chapter twenty-nine computes candidate counts and compares them with the sky. Selecting the physical capacity condition from the same repair history remains necessary.

## Newton's constant is an area

On the physical area-law branch, boundary entropy is its area divided by four times Newton's constant, in the same natural units. Read that from right to left: Newton's constant is the area divided by four times the entropy. It supplies an area scale for counting distinguishable states. With the measured gravitational conversion, one square meter of boundary has entropy of order ten to the sixty-ninth.

Chapter thirteen counted readings without assigning a width in meters. A specified screen-area dictionary can join the two. Give a cell an area equal to a whole number of one small area unit, and entropy equal to a quarter of that same number. Divide its area by four times its entropy and the number cancels, leaving the area unit. The cancellation is exact. Identifying this dictionary with the physical entropy and area is a separate premise.

Under that identification, Newton's constant is the unit of area. Its square root supplies the corresponding gravitational length.

The measured Planck area is 2.61228 times ten to the minus seventy square meters, with square root 1.6163 times ten to the minus thirty-five meters. Multiply that area by the cube of the speed of light and divide by Planck's constant over two pi, and the measured gravitational constant comes back: 6.67430 times ten to the minus eleven cubic meters per kilogram per second squared. These are conversions of the laboratory value. They do not independently measure a grain in the proposed screen.

The figure carries meters, kilograms and seconds, which people chose. The proposed screen dictionary concerns relations that survive a change of units; its physical attachment still needs an independent area and entropy calibration. The laboratory figure itself is known poorly. Newton's constant is the least precisely known of the fundamental constants by a wide margin, measured to five figures where the electron's magnetic moment is known to twelve, because gravity is feeble and everything in the laboratory has mass. MICROSCOPE measured gravity's indifference to what falls at one part in a thousand million million, and gravity's strength is known to twenty-two parts per million, ten orders of magnitude worse, because a difference between two cylinders nulls out everything the two have in common and a magnitude has to be assembled out of a kilogram, a meter and a second.

## Order, number and curvature

Chapter seventeen constructed a flat spacetime limit under a declared population, read law and counting measure: an order from who read whom, a volume from how many. Extending that route to the curved geometry needed here requires the corresponding order and volume identifications on one common refinement.

Start with the order. Within the distinguishing Lorentzian manifolds covered by the theorems of Hawking, King and McCarthy and of Malament, causal order determines the conformal geometry. It fixes the light cones but leaves a local scale: multiplying the metric by a positive function leaves those cones unchanged. A calibrated volume measure fixes that scale. If record order and counts converge to the required causal order and physical volume on the same smooth limit, they therefore determine its metric. An arbitrary finite log does not meet those conditions merely by having an order and a count; reading curvature also needs the appropriate smooth control.

Spatially flat homogeneous expansion illustrates the division. In conformal coordinates its metric is the flat metric multiplied by the square of a growing scale factor. The cones remain the same, while the four-volume density grows as the fourth power of that factor. If calibrated record counts sample equal small conformal-coordinate four-volumes, their density ratio recovers the fourth power of the scale-factor ratio. For comoving emission and reception, the fourth root is then one plus the redshift. The statement uses that cosmological geometry and count-volume identification; it does not follow from an arbitrary pair of record totals.

Curvature in general is more than one growing number. With count-volume calibration and the small-diamond expansion in place, counts can measure its effect on volume. Boundary area has a separate finite construction: in the declared three-dimensional population and all-neighbor read law, crossing reads approach a surface count with density equal to a quarter of pi times the square of the population density times the fourth power of the reach. This supplies a geometric area proxy for that source. Joining it to physical generalized entropy and to the same curved field history remains a further requirement.

The constructions therefore give concrete candidates for parts of the derivation: order, volume and boundary counts, finite entropy identities, and modular clocks. They do not yet establish that one repair history supplies all the physical identifications together, or that its generalized entropy is stationary at fixed physical volume. If those inputs and the required continuum bounds hold, the small-ball argument yields Einstein's equation. Establishing the responsive geometry is part of completing that physical construction.

Under those premises, the geometric equations close. Ten numbers at every point, a rule fixing how they answer to stress, the factor eight pi from the ball calculation, and a gravitational area scale requiring physical identification. In that description free fall is going straight ahead, and the two walkers on the equator were doing physics the whole way to the pole.

The geometric description is still an effective one. Its field equation relates curvature to energy and momentum; it does not select an internal gauge algebra or matter sector. The port-current constraints and the physical matter lift pose those separate questions in the next chapters.
