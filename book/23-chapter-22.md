# 22. Why Does Writing Something Down Cost Heat?

In 1824 a French military engineer of twenty-eight paid a Paris publisher named Bachelier to print six hundred copies of a book about steam engines, and then waited for a scientific community that did not read it. Sadi Carnot's *Reflections on the Motive Power of Fire* runs to a hundred and eighteen pages and asks the question every engine builder in Europe wanted answered, which is how much work you can get out of a fire.

His answer has no pistons in it. The most work any engine whatever can deliver, per unit of heat it swallows, is fixed by two numbers: the temperature of the thing it takes the heat from and the temperature of the thing it dumps the remainder into.

$$\eta = 1 - \frac{T_{\text{cold}}}{T_{\text{hot}}}$$

The Greek letter names the efficiency, the fraction of the heat that comes out as work. The two temperatures are the hot source and the cold sink, both counted from absolute zero. Brass or iron, steam or air, a better valve, a cleverer linkage, a bigger boiler: the ceiling sits where those two temperatures put it and no engineering moves it.

Carnot got that while believing heat was a weightless fluid called caloric, which fell from hot to cold the way water falls down a millrace and did as much work on the way. Caloric does not exist. Every heat engine built since has come in under the number he got out of it.

Cholera reached Paris in the spring of 1832 and killed something near twenty thousand people in the city that year. Carnot caught it on 24 August and was dead within the day, at thirty-six. His belongings were burned as a precaution against contagion, and nearly all of his papers went with them. What survived was a set of working notes his brother Hippolyte published in 1878, in which the fluid had been quietly dropped.

Look at what the surviving formula says. An engine's ceiling depends on the temperatures and on nothing else. That is a peculiar thing for a fact about brass and steam to be. Chapter twenty-one reached the same quantity from the other end, out of an algebra and a state, with no furnace anywhere in the argument. So what is a temperature a count of?

The everyday version sits on your desk. A drive with no moving parts gets warm while it is being written to. Its heat depends on the device and protocol. The universal erasure bound concerns a more specific operation than writing a record.

## The doorkeeper

On 11 December 1867 James Clerk Maxwell wrote to Peter Guthrie Tait describing a way to break the second law of thermodynamics using a very small employee.

Take a box of gas at one uniform temperature, divide it with a wall, and put a door in the wall. Temperature is an average: some molecules in there are moving fast and some slow. Station a being at the door with good eyes and a light touch. When a fast molecule approaches from the left, it opens the door and lets it through to the right. When a slow one approaches from the right, it opens the door and lets it through to the left. The door is weightless and frictionless, so opening and closing it costs nothing.

Wait. The right side heats up and the left side cools down, out of a box that started uniform. Then run Carnot's engine off the difference, and run it again, forever, in a box that nobody has done any work on.

Maxwell called his employee a finite being and later grumbled that it was really more of a valve. William Thomson named it a demon in *Nature* in 1874, meaning the word in its Greek sense, a spirit working quietly in the background.

The demon survived every attempt to kill it for sixty-two years, because every attempt went looking in the wrong place. People examined the door, the hinges, the light the demon needed to see by. Leo Szilard moved the search in 1929, in a paper in *Zeitschrift für Physik* whose title says where he was looking: on the decrease of entropy in a thermodynamic system by the intervention of intelligent beings. He stripped the gas down to a single molecule in a box, with a partition dropped in the middle. The demon looks, records which side the molecule is on, and uses that one recorded fact to let the molecule push the partition outward and lift a weight. The work per cycle comes out as the temperature times the logarithm of two, out of a box at one temperature, from one bit of recorded information.

Rolf Landauer, at IBM, identified the relevant operation in 1961: erasure. Resetting an unknown memory state to a fixed value reduces its entropy. In the standard thermal setup, the memory initially has no correlations with a heat reservoir at temperature T, and the combined evolution preserves information. Entropy lost by the memory must then be accounted for in the reservoir and their correlations. For one equally probable bit, the mean heat delivered to that reservoir obeys

$$E \geq k_B T \ln 2$$

E is the mean heat delivered to the initially thermal reservoir, and T is its temperature. The constant converts kelvin to joules. The logarithm of two is the entropy decrease of an equiprobable bit reset to a definite value, with no usable side record reducing the erasure task.

Put a room temperature into it, three hundred kelvin, and the bill for erasing one bit comes to 2.87 divided by ten to the twenty-first, in joules. A single photon of green light at a wavelength of 550 nanometers carries 3.61 divided by ten to the nineteenth, in joules, which is enough to pay off about a hundred and twenty-six bits.

Charles Bennett finished the demon in a review of the thermodynamics of computation in
1982. The demon's memory is finite. It can run a few cycles free, filling that memory with
records of which side each molecule was on, and then it is full, and to take another
measurement it must reuse or extend that memory. In the ideal thermal cycle, resetting each independent fair bit costs at least the work that bit enabled it to extract. The demon was never
an engine. It was a machine borrowing against its own memory. The loan came due at the
same rate every time.

Antoine Bérut and colleagues put the number on a bench in 2012 and published it in *Nature*. A silica bead about two micrometers across, held in an optical trap shaped into two wells, one well for zero and one for one. Push the bead into a known well and you have erased a bit. Measure the heat that leaves. As the erasure cycles were slowed down, the mean heat per erasure came down onto that value and did not go below it.

Chapter three asked what a physical record costs. Landauer answers the erasure part of that question, once the memory, reservoir and entropy change are specified.

## Four laws, numbered out of order

Thermodynamics has four laws. They were not discovered in the order they are taught.

Start with the one everybody uses without noticing. Nobody has ever established that a bath and a bowl of soup are at the same temperature by putting the soup in the bath. You put a thermometer in the soup, wait, read it, then put the same thermometer in the bath, wait, read it, and compare two numbers taken minutes apart with an instrument that has been in contact with each and never with both at once.

That procedure is a bet on transitivity. If the thermometer settles with the soup, and the same thermometer settles with the bath, then the soup and the bath will settle with each other. The bet is a law: two systems each in equilibrium with a third are in equilibrium with each other. Ralph Fowler and Edward Guggenheim wrote in 1939 that it "could with advantage be known as the zeroth law of thermodynamics", the other three having been numbered decades before anybody noticed this one was doing work underneath them.

Transitivity is what makes temperature a number. A relation that is transitive, and symmetric, and holds between a thing and itself, sorts everything in the world into classes with nothing left over, and classes in a row can be labeled by numbers on a scale. Without transitivity you could say that this is hotter than that, pair by pair, and there would be no scale to write either of them on, no degrees, no thermometer, and no sentence of the form "the water is at 47 degrees".

The first law separates changes of energy into heat and work. Work changes the available energies at fixed probabilities; heat changes probabilities at fixed energies. For a finite update that changes both, comparing only the two endpoints adds the product of those changes. It vanishes for an ordered stroke holding one quantity fixed. In the differential law the product is second order, so the usual two-term formula is exact even when both vary continuously.

The second law arrived in chapter sixteen, in a laboratory in Canberra with a latex bead in water. In the finite model, choose a reference distribution and a stochastic channel that preserves it. Relative entropy to that reference cannot increase under the channel.

The third law is the one popular accounts skip. Walther Nernst put it to the Göttingen Academy in December 1905: as a system is cooled toward absolute zero, its entropy approaches a fixed floor. The floor is set by how many arrangements the system has at its lowest energy. One such arrangement and the floor is zero. Many, and it is the logarithm of how many, and sits there down to the last fraction of a degree.

Ordinary ice does this. Each oxygen atom in ice has four hydrogens near it, two close and two further off, and there are many ways to satisfy that rule across a whole crystal without ever violating it locally, which is chapter six's triangle of frustrated neighbors wearing a different mineral. Linus Pauling counted the arrangements in 1935 and got a floor of 0.805 calories per degree per mole. William Giauque and J. W. Stout measured the heat capacity of ice down to fifteen degrees above absolute zero and reported 0.82, give or take 0.05, in a paper of 1936. That is 3.4 joules per kelvin per mole of entropy that stays in a block of ice at the bottom of the temperature scale.

The other half of the third law is a prohibition. You cannot reach absolute zero in a finite number of steps, however good your refrigerator, because each step takes a fraction of what is left rather than a fixed amount. A sodium gas at the Massachusetts Institute of Technology was cooled in 2003 to 450 picokelvin, give or take 80, which is four hundred and fifty trillionths of a degree above a floor nobody reaches.

## One rule, read twice

Four laws, arrived at over eighty years by people working on engines, gases, chemical affinities and cold. A finite observer model recovers their conditional counterparts from a supplied reference, a compatible resampling rule, and an identification of energy and temperature.

Read it first as a rule about descriptions. Jaynes's instruction from chapter nineteen applies unchanged: among all the descriptions consistent with what has actually been measured, take the one with the largest entropy, the one that adds nothing else. Apply that with one quantity held fixed, the average energy, and the answer is forced. The odds fall off exponentially with energy, cheap arrangements common and expensive ones rare, with a single multiplier in the exponent setting how fast the fall is. Chapter twenty-one ran into that exponential from the far side, in a two-outcome system whose flow came out as the energy divided by the temperature.

The positive multiplier becomes inverse temperature after the energy scale and thermal reference have been identified: beta equals one over k_B T. A probability distribution alone supplies no temperature in kelvin. On that thermal branch, erasing an equiprobable bit has a minimum mean heat cost of k_B T times the logarithm of two, or 2.87 divided by ten to the twenty-first joules at three hundred kelvin.

On a specified spectrum with two distinct energies, equal Gibbs distributions have equal inverse-temperature multipliers. Equality is transitive, which gives the finite zeroth-law identification. A spectrum with only one energy cannot distinguish temperatures from its probabilities. Interpreting equal multipliers as equilibrium between a thermometer, soup and bath additionally requires a physical contact model.

Read the rule the second way, as a rule about transitions. A repair step is handed some things that neighbors have settled between them and some things nobody has settled. What should it do? The same instruction applies: change nothing you do not have to, and assume nothing you have not been given. Leave every settled quantity exactly as it stands, and redraw everything else from the reference, inside the set of arrangements that agree with what is settled.

That instruction picks out one map. Applying it twice does nothing that applying it once did not do, because the second application finds the same settled facts and the same reference. It leaves the reference where it is, since redrawing from the reference cannot move it. And it preserves the average of every quantity that can be read off the settled part, which is the first law in the form chapter nine derived it: a repair that changed a total would have to know the total. The heat and the work come off the same map. Shifting the energies of the arrangements without touching the settled odds is the work channel, redrawing the odds at fixed energies is the heat channel, and the cross term is what a step that does both at once picks up on the way through.

For the declared resampling map, distance to its preserved reference cannot increase. That is the second-law inequality of chapter sixteen. It compares two probability descriptions under one channel; it does not identify the reference with physical truth or identify every source repair with that channel.

Full resampling settles the unresolved probabilities in one step. Local repair can take longer while preserving the same public record. An exact eight-state example keeps one recorded bit fixed and updates two other bits from their conditional probabilities. Their correlations decay toward equilibrium; the recorded bit remains unchanged. This supplies a finite model of memory and relaxation under a specified transition law. Identifying that law and its clock with a physical system requires further evidence.

Both constructions preserve every initially positive state weight. A zero-temperature Gibbs state has zero excited-state weights when excited states exist, so finitely many such steps cannot reach it from full support. Its physical temperature reading uses the specified energy and thermal reference.

Landauer’s bound follows on the identified thermal branch. If a reference-preserving heat stroke lowers entropy by c nats at positive inverse temperature beta, it expels at least c divided by beta of mean energy. Converting that energy decrease into reservoir heat uses the physical heat-stroke model. Setting c to the logarithm of two gives the one-bit bound.

## A one-way loop that obeys the law

Chapter sixteen’s pointwise fluctuation relation uses an extra condition: measured against the reference, the weight a step carries from one state to a second matches the weight it carries back. This is **detailed balance**. The integral identity requires stationarity and its positivity assumptions, not detailed balance.

The second law does not need it. The cheapest way to see that is a machine with three states and a permanent circulation in it.

Label three states one, two and three. From each state, the rule is: stay where you are with probability one half, or step clockwise to the next state with probability one half. There is no counterclockwise move at all. Nothing in this machine can go from state two to state one, ever, except the long way through state three, which takes four steps on average.

Check what it does to the flat description that puts a third of the weight on each state. State one receives half of its own third, from staying, and half of state three's third, from the clockwise step. That is a third. So is every other state, by the same count. The flat description does not move, which makes it the reference.

Detailed balance fails outright. The traffic from state one to state two is a third of a half, one sixth per step. The traffic from state two to state one is zero, because that move does not exist. There is a current going around the loop that never dies down.

Watch the second law hold anyway. Start with all the weight on state one, a description sitting the logarithm of three away from the flat one, 1.585 bits. One step spreads it into halves on states one and two, and the distance falls to 0.585 bits. A second step gives a quarter, a half and a quarter, at 0.085 bits. A third step brings it to 0.024. The descent is exact, it never reverses, and it happens inside a machine with a one-way street in it.

Stationarity is what the second law asks for: that the reference be left alone. Detailed balance is a stronger and separate demand. The extra content it carries buys the detailed fluctuation relations of chapter sixteen and Lars Onsager's paired transport coefficients. A world could have the second law with none of those and this loop is what it would look like.

## The price of settling

Take a memory with two equally likely inputs and reset both to the same output. Its entropy falls from one bit to zero. The same count applies to a reduced seam description only when those are its input probabilities and the operation actually identifies the two inputs.

A comparison or a retained readback need not perform that erasure. Keeping the input in another register preserves the distinction in the combined system. Counting commits alone therefore does not count erased bits.

Bare Shannon entropy and relative entropy are different quantities. The stationary-reference inequality controls the latter. Deterministic settling can reduce the former, but its information loss alone supplies neither a thermal reference nor a heat measurement.

For a concrete thermal example, take the seven equiprobable inputs of chapter sixteen and erase their distinction, with no usable side record. If a physical memory performs that reset with an initially uncorrelated thermal reservoir at three hundred kelvin under the standard Landauer conditions, the mean heat delivered to the reservoir is at least k_B T ln 7: about 8.06 divided by ten to the twenty-first joules. This is a bound per such erasure. A semantic commit supplies neither that entropy decrease nor that reservoir temperature.

## Five entropy readings

The same entropy mathematics appears in five kinds of reading, though the physical dictionaries connecting them have to be supplied and tested.

An engine reads entropy through heat and temperature. A gas reads it through Boltzmann's count of arrangements. A memory reads a Shannon count of questions, with Landauer giving a lower heat cost for physical erasure. A black-hole horizon carries the Bekenstein-Hawking entropy fixed by its area. An erasure loses alternatives from an informational record. These quantities share formal relations, but the finite record count becomes horizon entropy only after a physical carrier, temperature, energy and entropy map identifies the two readings.

A logically irreversible physical write has a Landauer lower cost when its device and thermal environment satisfy the law's premises. The finite observer model counts discarded alternatives and boundary channels without supplying that laboratory realization. A record count comes back in square meters only on a separate horizon-record branch that calibrates the count against physical entropy and area.
