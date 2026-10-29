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

Jacobson brought the area law and the temperature in from outside, and the third ingredient, the demand that the books balance on every horizon at once, he postulated. This world supplies all three. Chapter twenty-four counted a region's memory on the screen around it and found a quarter of the area in Planck units. Chapter twenty-one showed that the flow a region's own state determines is the boost of an accelerating observer, running at the temperature that observer reads. That leaves the demand, and it is the interesting part, because here it is no demand at all. It is something repair has finished doing.

## Why every small ball has settled

Chapter fifteen struck the carrier like a bell and listened. The fast band halves in 5.4 passes and the middle band in 6.6, and after a hundred passes the fast band holds under three parts in a million of what it was handed. What outlives them is the slow band, and the slow band is space. So on any scale large compared with a pass of repair, every small region has finished settling: whatever could relax has relaxed, and what is left is the part that carries positions. Gravity lives on that scale. It is a statement about settled regions, which is why it appears when the description is coarse-grained and has nothing to say at the grain itself, where the fast content is in the middle of dying away.

What does a settled region hold? Chapter twenty-three answered that for heat. Among all the descriptions that agree with what has actually been settled, a repair step keeps the one that adds nothing else, the one with the largest entropy. For a small ball the settled facts are the records on its boundary. Everything the rest of the world can learn about the inside arrives through that surface, so the boundary records are the whole of the constraint, and the settled ball holds the largest entropy those records allow.

Count everything a ball holds. Its entropy has two parts: what the contents hold, and what the boundary holds. The second part is chapter twenty-four's count, a quarter of the boundary area in Planck units, because the screen around the ball is made of patches that each hold the same amount of memory. The sum of the two is the **generalized entropy**, and it is the number repair has pushed as high as the boundary records permit.

One more thing has to be held fixed, and chapter seventeen says what. A volume in this world is a number of record events: the order says which events precede which, and the count says how much there is. Comparing the settled ball with a slightly different ball of the same volume compares two ways of arranging the same number of records. That is the only comparison that means anything here, since a ball with more records in it is simply a bigger ball.

At a maximum, a small step in any allowed direction changes the height by nothing at first order. The top of a hill is flat. A ball whose generalized entropy had a slope in some allowed direction would have a direction repair could move it in, and a region with a repair left to make has not settled. So in every settled ball the first-order change of the generalized entropy vanishes, in every direction that keeps the volume fixed. That is the whole of the stationarity gravity needs, and it holds because repair has run.

Jacobson came back to the problem in 2016, in volume 116 of *Physical Review Letters*, with a sharper version of his argument. He asked that the entanglement entropy of every small ball of empty space be at its maximum at fixed volume, and showed that Einstein's equation follows. He had to postulate the maximum and called it a hypothesis about the vacuum. Here it is what the word settled means.

## Heat that buys entropy

Settled is different from empty. Put a little energy into the ball and the entropy of its contents moves. By how much?

Chapter sixteen priced a description against a reference by its relative entropy, the extra questions per state the description costs. It is never negative, and it is zero exactly when the description is the reference. Picture it as a bowl with the settled state at the bottom. Near the bottom of any smooth bowl the height grows as the square of the distance: go a tenth of the way out and you have climbed a hundredth. So a small disturbance of the settled state changes the relative entropy at second order and leaves its first-order change at zero.

Written out, the relative entropy of the disturbed state against the settled one is the change in the region's energy, read on the region's own clock, minus the change in its entropy. Chapter eighteen made that energy the generator of the region's own flow. A difference whose first-order change is zero is a difference of two first-order changes that are equal.

$$\delta S = \delta \langle K \rangle$$

Here S is the entropy of the ball's contents, K is the generator of the ball's own clock, the angle brackets take its average over the ball's state, and the delta marks a first-order change away from the settled state. The entropy the contents gain equals the energy that came in, counted on the clock the ball keeps for itself. That is Clausius's relation, heat over temperature equals entropy, with the temperature folded into the clock. It follows from the bowl having a bottom.

## Where the eight and the pi come from

Take a small ball around a point, in the rest frame of one observer, and do the bookkeeping.

The contents side first. The ball's own clock is a flow that holds the ball's edge fixed and carries its inside toward the future, the counterpart for a ball of chapter twenty-one's boost. It works hardest at the center and not at all at the edge, so it weights the energy density at each point by the difference between the squared radius of the ball and the squared distance from the center, divided by twice the radius. Integrating that weight over the ball gives four pi times the fourth power of the radius, divided by fifteen. The integral takes one line, and the fifteen comes out of it.

One factor remains. The ball's flow closes up on itself after a turn of two pi in its own parameter, which is the same fact that makes an accelerated observer read a temperature equal to its acceleration divided by two pi, in units where the speed of light, Planck's constant and Boltzmann's constant are all one. Reading the flow's parameter as energy brings that two pi in, taking the matter side to eight pi squared times the fourth power of the radius, divided by fifteen.

The geometry side. Alfred Gray and Lieven Vanhecke, in *Acta Mathematica* in 1979, worked out how the size of a small geodesic sphere depends on the curvature around it. The piece of that needed here says that if you hold the volume of a small ball fixed and raise the curvature, the area of its boundary falls by four pi times the fourth power of the radius, divided by fifteen, times the change in exactly the curvature combination that stands on the left of Einstein's equation. Curving the surroundings costs boundary area.

Add the two changes and set the sum to zero, as a settled ball requires, remembering that the area enters the generalized entropy divided by four times Newton's constant. At leading order every appearance of the radius cancels. Both fifteens cancel. What is left is eight pi squared on one side against pi over Newton's constant on the other, and their ratio is eight pi times Newton's constant, the coupling in the field equation. Nothing was fitted to get the eight pi.

That gives the balance in one observer's rest frame. Every observer passing through the same event, at any velocity, has a small ball of its own, and every one of those balls is settled. Written out in terms of the velocity, each balance is a polynomial of degree two. A polynomial that vanishes throughout an open range has every coefficient equal to zero, and the coefficients are the components of the field equation. Every observer's balance holds at once, and the field equation holds as a tensor, up to the one term the light rays below are blind to.

## A geometry that can respond

The balance cannot hold in a geometry that never moves. Send energy across a ball along a light ray and the entropy of the contents moves at first order, so the boundary area has to move by the opposite amount, or the ball would have a slope repair could climb. Something has to give, and the only thing left to give is the shape.

The response is specific. If the area of a spherical slice is four pi times a radius squared, energy passing along an affinely measured light ray changes the second derivative of that radius by minus four pi times Newton's constant, times the radius, times the null stress. Light rays crossing a region that holds energy converge, the way the two walkers did. That is the focusing part of Einstein's equation, and it is what the settled count requires of the geometry.

The balance counts energy and has no other column. The ball's own clock weighs whatever is inside by the energy that clock assigns. Heat counts, light counts, the pressure of a gas counts, and so does the upkeep of a particle. Chapter eighteen read mass as the standing cost of a pattern repair cannot flatten, and chapter twenty-eight shows those patterns are knots in the transport between patches, whose rest energy is the repair their cores cost per tick on the same clock. The area answers to that cost exactly as it answers to light. The geometry has no way to ask whether the energy is titanium or platinum, which is the equivalence principle, and why MICROSCOPE found nothing at one part in a thousand million million.

## A gas law without molecules

Jacobson ended his abstract with a sentence that has been argued over for thirty years: it may be no more appropriate to canonically quantize the Einstein equation than it would be to quantize the wave equation for sound in air.

Pressure times volume equals a constant times temperature. Émile Clapeyron put it in that form in 1834. It held up in every laboratory in Europe through the rest of the century without containing a single molecule. In 1908 Jean Perrin put grains of gamboge resin in water under a microscope, counted how the population thinned out with height, and got Avogadro's number out of the count, which finished an argument about whether matter came in pieces at all. Wilhelm Ostwald, who had spent a career arguing that atoms were a bookkeeping convenience, gave that position up in the preface to the fourth edition of his own textbook in 1909, citing Perrin by name. The gas law did not change a symbol, because it was never a claim about molecules. It is a claim about a gas that has settled down, and stays true whatever the gas turns out to be made of.

Einstein's equation belongs on that shelf. It is the equation of state of settled records, true of any region that has finished settling, whatever its records are written in. Quantizing the metric in the hope of finding the atom of space resembles quantizing the pressure of a gas in the hope of finding the molecule. The pressure is real enough to lift a piston, and it belongs to a settled state rather than to a constituent. The constituents of space are the records, and they were never in the metric.

With the smooth metric in place, the 1907 thought has its usual precise meaning. A freely falling frame removes the connection at a point, while tidal curvature remains. For pressureless matter with nonzero density, conservation of the stress makes the flow geodesic. The walkers return as a consequence.

## Nine light rays

Everything in the argument crossed a horizon, and horizons are made of light rays.

The horizon probes in this argument evaluate stress along a light-ray direction. The stress takes two directions and returns a number. Feed it the same null direction twice and you get the number the balance uses.

The cone those rays travel on is chapter twenty-two's first response. Disturb the field in one region and read a detector in another, and the detector's response is exactly zero until the clock reaches the distance between the two regions, and switches on after that.

Ten numbers describe the stress at a point. Light rays report one number each. How much of the stress do the rays determine?

A symmetric object that reads zero along every light ray has to be a multiple of the metric, and nothing else. So the rays pin down nine of the ten numbers. What they leave undetermined is a multiple of the metric.

For readings from one symmetric tensor, nine directions do it. Take the six rays that travel along one axis at a time, forward and backward along each of the three. Then take three more that share their motion equally among all three axes, each space component equal in size to one over the square root of three, with the signs plus plus plus, plus plus minus, and plus minus plus. The six axis rays never pick up the cross terms, the entries pairing one space direction with a different one, because a ray running along a single axis has no second direction to pair it with. The three tilted rays are there for those, and their three sign patterns read the three cross terms off one another. Nine light rays, nine numbers, one system of linear equations to invert. Its determinant is exactly eight thousand one hundred and ninety-two over twenty-seven, or two to the thirteenth over three cubed. A determinant that misses zero is a system that inverts, and this particular inversion is stable: small errors in the nine readings stay small in the nine numbers recovered from them.

Any nine numbers fit that system, so the fit alone does not show that the readings come from a tensor: a tenth ray, or a hundredth, has to land where the recovered tensor says it will. Given the geometric and stress tensors, balance along the nine rays fixes the field equation up to a multiple of the metric. Conservation makes that multiple constant on each connected region. That term is the tenth number.

## The vacuum and the tenth number

The vacuum is the settled reference, and it looks the same to every observer however fast that observer moves, so its contribution to stress is a number times the metric.

Contract it with a light ray. The metric evaluated on a light-ray direction twice is the interval of a light ray with itself. That is zero: chapter seventeen defined a light ray as a separation whose interval is exactly zero. So the vacuum's contribution to what a light ray reports is a number times zero.

Exactly zero, for any vacuum energy density whatever, however enormous, and zero by the definition of the probe rather than by a cancellation between large terms.

The single direction the light rays cannot see is exactly the direction the vacuum's energy points in. That is why the cosmological constant survives in the field equation as a term nobody derived: the small-ball calculation that produced everything else on the left-hand side is built out of horizons, horizons are made of rays, and rays are blind to it.

A naive Planck-scale cutoff on quantum zero-point energies gives a vacuum energy density around a hundred and twenty orders of magnitude larger than the one inferred from the sky, whose cosmological constant is about 1.1 times ten to the minus fifty-second per square meter. The local balance never sees that number, so nothing ever had to cancel it.

The horizon fixes the tenth number. The largest settled region an observer has is everything inside its cosmological horizon, and that horizon can hold only so many records. The cosmological constant is set by the count: the smallness of the constant is the largeness of the horizon's memory. Chapter twenty-nine counts the records and sets the count beside the sky.

## Newton's constant is an area

Boundary entropy is the area divided by four times Newton's constant, in natural units. Read that from right to left: Newton's constant is the area divided by four times the entropy. It is the area scale for counting distinguishable states. With the measured gravitational conversion, one square meter of boundary holds an entropy of order ten to the sixty-ninth.

Chapter thirteen counted readings without assigning a width in meters, and the screen joins the two. Give a cell an area equal to a whole number of one small area unit, and entropy equal to a quarter of that same number. Divide its area by four times its entropy and the number cancels, leaving the area unit. The cancellation is exact.

So Newton's constant is the unit of area: the piece of screen that holds a quarter of one natural unit of memory. Its square root is the gravitational length.

The measured Planck area is 2.61228 times ten to the minus seventy square meters, with square root 1.6163 times ten to the minus thirty-five meters. Multiply that area by the cube of the speed of light and divide by Planck's constant over two pi, and the measured gravitational constant comes back: 6.67430 times ten to the minus eleven cubic meters per kilogram per second squared. Those are the laboratory's numbers, and in this world they describe one object, the size of the piece of screen that holds a quarter of a natural unit of memory.

The figure carries meters, kilograms and seconds, which people chose. The relation it encodes, area against memory, survives every change of units. The laboratory figure itself is known poorly. Newton's constant is the least precisely known of the fundamental constants by a wide margin, measured to five figures where the electron's magnetic moment is known to twelve, because gravity is feeble and everything in the laboratory has mass. MICROSCOPE measured gravity's indifference to what falls at one part in a thousand million million, and gravity's strength is known to twenty-two parts per million, ten orders of magnitude worse, because a difference between two cylinders nulls out everything the two have in common and a magnitude has to be assembled out of a kilogram, a meter and a second.

## Order, number and curvature

Chapter seventeen got flat spacetime out of an order and a number: who read whom, and how many records there are. Curved spacetime comes out of the same two.

Start with the order. Within the distinguishing Lorentzian manifolds covered by the theorems of Hawking, King and McCarthy and of Malament, causal order determines the conformal geometry. It fixes the light cones but leaves a local scale, because multiplying the metric by a positive function leaves those cones unchanged. A volume measure fixes that scale. Record order converges to the causal order and record counts converge to the volume, on one common refinement, so together they determine the metric, curvature included.

Spatially flat homogeneous expansion shows the division of labor. In conformal coordinates its metric is the flat metric multiplied by the square of a growing scale factor. The cones are the same at every epoch, while the four-volume density grows as the fourth power of that factor. Record counts sampling equal small conformal-coordinate four-volumes have a density ratio equal to the fourth power of the scale-factor ratio. For comoving emission and reception, the fourth root of that ratio is one plus the redshift.

Curvature in general is more than one growing number, and counts measure its effect on volume through the small-diamond expansion. They measure areas as well. In a three-dimensional population with an all-neighbor read law, crossing reads approach a surface count with density equal to a quarter of pi times the square of the population density times the fourth power of the reach. That surface count is the memory chapter twenty-four put on the screen.

Order, volume, boundary area, the first law and the region's own clock all come off one repair history, and the settled balls tie them together. Ten numbers at every point, a rule fixing how they answer to stress, the factor eight pi from the ball calculation, and an area scale that is Newton's constant. In that description free fall is going straight ahead, and the two walkers on the equator were doing physics the whole way to the pole.

The geometric description is an effective one, a summary observers fit to settled records. Its field equation relates curvature to energy and momentum and selects no internal gauge algebra or matter sector. The port-current constraints of the next chapters, and the knots of chapter twenty-eight, answer those questions.
