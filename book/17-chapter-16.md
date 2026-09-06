# 16. Why Does Time Run, and Only One Way?

A latex bead a little over six thousandths of a millimeter across sat in a dish of water in a laboratory in Canberra, held in place by a focused laser. A trap of that kind pulls the bead toward the brightest point with a force that grows the further out it strays, which is how a spring behaves. The dish was then dragged sideways past the trap at 1.25 micrometers per second while a quadrant photodiode read the bead's position a thousand times a second, for as long as ten seconds at a stretch.

Dragging a bead through water takes work. The work finishes up in the water as heat. Scored across a few seconds and longer, every one of the five hundred and forty recorded runs did what the second law of thermodynamics requires: the bead lagged behind the trap, the water took up the difference, and the whole business ran the way anybody would expect it to. Score the same recordings over a hundredth of a second instead, and the stretches that run backwards come out about as numerous as the stretches that run forwards. The water hands energy to the bead rather than taking it. The entropy of the arrangement goes down. Widen the window to two seconds and the backwards stretches thin out; widen it past a few seconds and they stop occurring.

G. M. Wang, E. M. Sevick, E. Mittag, D. J. Searles and D. J. Evans published that in *Physical Review Letters* on 29 July 2002, under the title "Experimental Demonstration of Violations of the Second Law of Thermodynamics for Small Systems and Short Time Scales", seventeen words that leave the abstract with nothing to announce.

A rule that fails on short recordings and holds on long ones is a statement about an average. So some quantity is being averaged, the average comes out with a sign, and something has to say where the sign comes from. That is one of the three questions hiding inside the word time, and the other two are older.

## Three jobs, one word

Ask anybody what time is and you get a description with four separate claims folded into it: events fall in a single line, the line runs forward and cannot be run backward, the line is the same for everybody, and it goes at one rate. Chapter seven took away the third one. Two clocks in different places have no shared present in which to disagree.

The other three claims are three different objects with one word stretched over them.

The first is **order**: which of two things happened before which. Order is a relation between events. It has no units, no rate and no duration in it. Chapter seven's recipe carried the whole of it, four steps wired in a line and a fifth that floats.

The second is **direction**: why the relation runs one way round and never the other. That is a property of an operation rather than of a relation, and the one the bead in Canberra was testing.

The third is **flow**: why anything seems to be moving, and why there is a rate at all. Order tells you which of two events came first, where there is a fact about it at all. Direction says why the pair never runs the other way round. Between them they say nothing whatever about passing. You could be handed the complete order of every event in the universe and the direction of every strand in it and have, out of that, no motion, no duration, and nothing that goes by. Chapter seven's recipe card carries an order, and printing it down a card draws an arrow along that order, and the card lies on the counter going at no rate whatever. The answer to flow is a separate answer, arrived at by a different route. It is not the arrow of time experienced from the inside. The state an observer holds determines its own clock, the way a piece of music determines a tempo without a metronome anywhere in the room. Five chapters from here, that stops being an image.

## What gets ordered

Chapter seven taught partial orders on a kitchen counter, before this world contained anything to order, and chapter nine supplied the missing ingredient. Events are authenticated semantic commits. An accepted repair commit is one possible event when its reads and writes are recorded; record, readback and feedback commits can also be events. The archived source-history receipt uses the latter instrumentation carrier and excludes recurrent seam-repair transactions from its physical-event interpretation.

Take the twelve observers of chapter six, wired as an icosahedron, thirty seams between them, every observer with five neighbors. A repair reads the two readings at the ends of one seam, writes the two readings at the ends of that same seam, and touches nothing else in the world.

Sharing an observer makes two repairs possible competitors for the same registers. It does not by itself order them. An informational edge appears when the later accepted commit certifies that it read a register version written by the earlier one. If it merely follows the other repair in an executor's list, there is no edge. Repairs whose certified resources are disjoint have no such fact, and no clock or queue position can add one afterwards.

Two pairs, concretely. A repair on the seam between the first observer and the second may write a version that a later repair between the second and third certifies reading. If it does, the writer mark and read certificate establish the dependency. Against that, a repair between the first and second observers and one between the fourth and fifth touch disjoint registers. Ask which came first and the record system has nothing to interrogate.

Count the possible conflicts. Thirty seams make four hundred and thirty-five pairs of seams. A given seam has two ends, and at each end four other seams arrive, so eight seams share an observer with it and twenty-one do not. That gives one hundred and twenty pairs that may compete for a register and three hundred and fifteen that are structurally disjoint. These are topology counts, not counts of informational precedence pairs. The actual order comes from the versions read in one execution.

That is a great deal less order than a line would carry. A line puts everything in single file and answers all four hundred and thirty-five pairwise questions before a single repair has run. The observers instead produce a web after the fact. Every strand is a certified read-after-write dependency, and unrelated commits remain incomparable. The wiring says where dependencies are possible; the provenance says which ones occurred.

The everyday picture flattens that web into a timeline. The difference shows up in what you are allowed to ask. "Which of these two happened first" has an answer for some pairs and none for others. The wiring constrains which dependencies are possible, while the recorded provenance fixes which comparabilities actually occur.

## Seven pairs go in, one comes out

Chapter nine derived the repair move from three requirements and got the midpoint. Two readings at a seam. Both of them go to the average.

Watch what that does to counts. A seam whose two ends hold 6 and 0 finishes with both ends holding 3. So does a seam holding 5 and 1. So does 4 and 2, and 2 and 4, and 1 and 5, and 0 and 6, and 3 and 3, which was sitting on the answer to begin with. Seven arrangements of two whole numbers go into that move and one arrangement comes out, and every one of the seven is consistent with what the seam says afterwards.

Run it backwards. You have two readings of 3. The question is what they were. Those readings do not say. Without another record, all seven inputs remain compatible with the output. Averaging has removed the distinction from this description of the seam. Whether another part of the system retains it is a separate question.

This repair write leaves the total alone and satisfies the seam while discarding which input it received from the two-reading description. A semantic commit can record such a write, but commits can also retain readback and history. The acceptance rule restricts which repairs are allowed; it does not decide what information survives elsewhere.

Those two restrictions differ. The acceptance test is a rule about available moves. Seven into one is a fact about the averaging map. Relaxing the rule cannot recover an input from the output alone.

A write can be made undoable: retain the input in an initially blank register. The combined output and record then distinguish all seven cases. No erasure is forced by keeping that copy. Reusing finite memory cyclically may require resetting it later, and the cost belongs to that reset under its physical operating conditions. Logical erasure and reversible record retention are different operations.

Order is a relation between two records, checkable by whoever holds both and carrying no units. Direction is a property of a single move, showing up as an arithmetic loss inside that move. Drawing one axis and putting an arrow along it welds the two into a single object. The weld is why people go looking for the arrow in the first instant of the universe.

Boltzmann’s account concerns a different arrow. High-entropy arrangements vastly outnumber low-entropy ones, and an explanation of their prevalence must also account for the initial conditions. David Albert called that initial-condition assumption the past hypothesis in *Time and Chance* in 2000. The seven-into-one map establishes an information loss in a reduced description. Connecting it to sustained physical irreversibility requires a state ensemble, dynamics and an environment; the arithmetic of averaging does not settle the cosmological initial-condition question.

## Counting what got lost

The seven arrangements that collapse into one are a quantity of missing information. Ask how many yes-or-no questions it takes to pin down one item out of a list. Two items, one question. Four items, two questions, since the first splits the four into halves and the second splits a half. Eight items, three. Sixteen, four. Each doubling of the list costs one more question, which is the defining habit of the **logarithm**: it converts multiplying into adding. Counting questions this way makes the logarithm base two. The unit of the answer is the **bit**.

The list is rarely flat. If one item out of the four is likely and the other three are long shots, a well-chosen first question separates the likely one from the rest, and most of the time you are finished after one question rather than two. So the count wanted is an average over the odds. Claude Shannon wrote it down at Bell Labs and published it in the *Bell System Technical Journal* in July and October of 1948, in two halves of a paper called "A Mathematical Theory of Communication".

$$H = -\sum_i p_i \log_2 p_i$$

Here H is the average information in one outcome, measured in bits; the index i runs over the possible outcomes, and p sub i is the probability of outcome i. The minus sign compensates for the negative logarithms of probabilities below one.

Three cases, by hand. A fair coin gives 1 bit. A coin biased nine to one gives 0.469 bits: heads contributes 0.137 and tails 0.332 to the average. A fair eight-sided die gives 3 bits. Assign equal probabilities to the seven inputs at the seam and their entropy is log base two of seven, or 2.807 bits. Mapping that ensemble to one output removes that much entropy from the seam description. Unequal probabilities change the amount, and a retained input record preserves the distinction in the larger system.

Shannon told Myron Tribus that he had gone to John von Neumann for a name, and that von Neumann told him to call it entropy, on the grounds that the same expression appeared in statistical mechanics and that nobody knows what entropy is, so in an argument he would have the advantage. Tribus and Edward McIrvine printed the story in *Scientific American* in September 1971.

The entropy belongs to the description rather than to the thing described. The seam holds 3 and 3 whatever anybody knows. The 2.807 bits are a fact about somebody's odds over the seven arrangements. A different set of odds gives a different number for the same seam. Chapter ten met this quantity under another word: the run that finished with 81,920 records had a spread of 11.3134, which is the same count in a different unit, the way a length is the same length in inches or in centimeters. In bits it is 16.32.

Which is where the schoolroom gloss goes wrong. Repair makes the world tidier and more determinate, and a seam that goes from 6 and 0 to 3 and 3 has had its uncertainty cut, so anybody reading entropy as a measure of how messy the seam is will get the sign backwards. The seam is as messy after the repair as before it: two whole numbers, sitting at two ends of one overlap. What changed is how many questions it takes to say which two.

## Two messages down one channel

The quantity that runs one way is a comparison rather than a count. Take a settled arrangement of the sort the repair law leaves alone, and call the odds it assigns to each state the reference. Take your own odds over states, which need not be the reference. Ask how many extra questions per state your odds cost you compared with the reference, and the answer is a second quantity, the **relative entropy** of your odds against the reference. It is zero exactly when the two agree, positive otherwise, and measures how far your description sits from the settled one.

Pass both descriptions through the same process, whatever the process is, and they get harder to tell apart. Send two distinguishable messages through the same noisy channel and the outputs are closer together than the inputs were; nothing you do to both of them afterwards can make them further apart again. Relative entropy is the number that quantifies that, and can only fall under further processing.

For a stochastic repair channel that preserves a supplied reference, relative entropy to that reference cannot increase. This is the finite second-law inequality. The deterministic averaging rule and the equilibrium resampling channel are separate constructions; applying the inequality requires the stated reference-preservation condition.

## The identity underneath the inequality

The falling is an average of something that does not always fall. Take one repair step. Before it, the world is in some state and your odds give that state a probability, and the reference gives it a probability too. The logarithm of the ratio between them is your surprise at finding that state. After the step, the world is in another state, and there is a second surprise, computed the same way against the same reference. The **entropy produced** by that transition is the first surprise minus the second.

Some transitions have positive entropy production and others negative. For strictly positive normalized initial and reference laws, a reference-preserving stochastic kernel, and a positive output law, the average exponential of minus that quantity is exactly one. Detailed balance is not needed for this integral identity. The positivity assumptions keep the logarithmic ratios defined on every required state.

One further condition is checkable at the seam: measured against the reference, the weight the step carries from one state to a second matches the weight it carries back. A repair in that balance satisfies the following at every pair of states, exactly.

$$p(x)\,K(x,y) = e^{\sigma(x,y)}\,q(y)\,K(y,x)$$

Here p of x is the probability your description gives the state before, K of x y is the chance the repair takes that state to the state after, q of y is the probability of the state after, K of y x is the chance the same repair runs the other way, and sigma is the entropy produced. The equation says the forward weight of any transition exceeds the reverse weight of the same transition by exactly the exponential of the entropy produced, at every pair of states, with nothing approximated. Gavin Crooks published the general relation of this form in *Physical Review E* in 1999, the relation the Canberra bead was measured against.

Two consequences follow. Collect all the transitions that produce the same amount of entropy and add up their weights, and the forward total is that same exponential times the reverse total. Count the entropy produced in bits and the exponential turns into a doubling per bit, which is the same statement in the unit the coins were counted in. A transition producing ten bits runs forwards about a thousand times as often as backwards. One producing forty bits, about a million million times as often. The arrow is a ratio rather than a prohibition, which is why a small enough system watched for a short enough time can be caught running the wrong way, and why nobody has ever caught a cup of coffee doing it.

The mean entropy production equals the fall in relative entropy exactly. Convexity of the exponential gives Jensen’s inequality: its average is at least the exponential of the average exponent. Since the first average is one, mean entropy production is nonnegative. This is an exact inequality between real-valued quantities, with no approximation.

Detailed balance also makes the equilibrium two-time correlation symmetric: measure one quantity before a step and another after it, then swap their positions. The averages agree. This is the finite analogue of the reciprocity associated with Lars Onsager, whose work connected paired transport processes such as heat and electric currents. Reading these finite correlations as physical transport coefficients additionally requires identified currents and a calibrated time scale.

## Faking a past

If the seam no longer distinguishes its inputs, a retained record can preserve access to them. The question is whether a claimed record can be forged.

Set the machine one task. A chain of records claims a history worth twenty bits, meaning twenty yes-or-no questions' worth of detail about what happened. A competing chain claims the same history and has no such history behind it, so the twenty bits it would need to produce are twenty bits it can only guess. A blind guess at twenty bits succeeds one time in two to the twentieth, which is one time in 1,048,576, or 9.5 times in ten million.

The machine ran a hundred thousand attempts. The expected number of successes under that blind-guess model is about a tenth, and the number obtained was zero. The per-attempt probability assumes independent uniform hidden bits and no usable side information. More attempts increase the chance of a success; a retained copy or predictable history would change the guessing problem.

The recorded past is the retained history of commits and their certified read dependencies. Its authentication rests on those records and checks, not on an assumption that every write destroyed its input. Unresolved repairs have no committed outcome in that history yet.

Chapter six left one unit of disagreement in the declared twelve-observer model. One flipped record, thirty seams, nineteen loops, and a leftover of exactly one that no allowed repair removes. Repair order changes where it sits. A retained history can record the writes that moved it and their certified dependencies, whether or not an input remains recoverable from the local readings.

Every strand in that web was laid by one commit reading what another commit wrote. Chapter nine's local rule bounds how many seams that dependency can cross in one model event. This is a graph-local propagation bound. Turning seams and events into meters and seconds requires a physical ruler, a clock calibration and proof that the informational order agrees with physical signal causality. Without those bridges it is not a derivation of the speed of light.
