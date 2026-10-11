# 7. Why Is There No Such Thing As Now?

On 4 October 1971 two men boarded an eastbound commercial flight out of Washington with tickets made out for a party of four. Two of the four were listed as Mr. Clock.

Mr. Clock was four cesium-beam atomic clocks, Hewlett-Packard model 5061A, each about the size of a large suitcase and too delicate for the hold. Joseph Hafele was a physicist at Washington University in St. Louis who had worked out, while preparing a lecture in 1969, that an ordinary airliner and an ordinary atomic clock were enough between them to measure something people had been arguing about since 1905. He spent the following year failing to find anybody who would pay for it. After a talk he gave in 1970, Richard Keating, an astronomer at the United States Naval Observatory, came up and mentioned that he worked with atomic clocks for a living and that his employer owned several. The whole thing cost eight thousand dollars, of which seven thousand six hundred went on economy fares.

They flew east around the world, came back, and flew west from 13 to 17 October. Then they carried the clocks back into the Naval Observatory and compared them against the clocks that had never left the bench.

Every clock in the building had been running the entire time. None of them agreed.

Going east, the flying clocks had lost 59 nanoseconds against the ones on the ground, give or take 10, where the prediction had been a loss of 40 give or take 23. Going west they had gained 273 nanoseconds, give or take 7, against a prediction of 275 give or take 21. Both papers went into *Science* in July 1972, one titled "Around-the-World Atomic Clocks: Predicted Relativistic Time Gains" and the other, published immediately after it by the same two authors, "Around-the-World Atomic Clocks: Observed Relativistic Time Gains".

Look at what those numbers are attached to. Each came out of subtracting one reading from another in one room, at the end, with the clocks side by side.

That is not the only comparison available, and why the others do not help is the whole of the point. Nothing stops you comparing clocks at a distance. Radio does it. A navigation satellite does nothing else all day: the plane could have sent its reading down and the observatory could have written down what arrived.

The trouble is that the message takes time. The reading arrives stale, and correcting for that means knowing how long the trip took. You can time a round trip without difficulty: send a pulse up, wait for it to come back, read your own clock twice. Splitting that total into the trip out and the trip back is where it stops. Calling the two halves equal is the usual move and a perfectly good working rule, but nothing measures it, and checking it would mean comparing two clocks that are far apart, which is the thing you were trying to do.

What did the flying clock read at the instant the bench clock read noon on the second day out? Every route to an answer goes through an assumption somebody made rather than a measurement somebody took.

Everyone carries a picture around without inspecting it. There is a single present moment, everything in the universe is doing something in it, and clocks are devices for finding out what. On that picture the question about noon on the second day has an answer that Hafele and Keating merely failed to measure. It has none. The single present moment is a part of the description rather than a part of the world. The evidence is that no instrument reads it and no procedure produces one.

Agreement between observers has to be built out of messages: one observer sends, another receives, both write. Do messages need a shared present to happen in?

## The messenger

Take away the clock and look at what is left, which is two observers and a channel.

Two armies sit on hills either side of a valley. The enemy is in the valley. Either army alone loses; both together win, provided they attack at the same time. The only way to communicate is to send a messenger down through the valley, where the messenger may be captured.

The first general sends a rider: attack at dawn. The rider may not arrive, so the first general cannot attack on the strength of having sent him. Suppose the rider does arrive. The second general reads the message and cannot attack either, because he knows that the first general does not know the message got through. A general who suspects he may be attacking alone stays put. So the second general sends the rider back with an acknowledgment. That rider may be captured, so the second general cannot rely on the acknowledgment having landed. The first general receives it, and cannot attack, because he knows the second general does not know the acknowledgment arrived.

Acknowledging the acknowledgment moves the problem one rider further down the road and changes nothing else.

E. A. Akkoyunlu, K. Ekanadham and R. V. Huber published the impossibility proof in 1975, in a paper called "Some Constraints and Tradeoffs in the Design of Network Communications". The argument is three lines long. Any protocol that solves the problem delivers some finite number of messages. Consider the last one delivered. Its sender cannot know it arrived, so the sender's decision cannot have depended on it, and the receiver's decision cannot depend on a message that a working protocol has to be prepared to lose. Delete it. Both generals behave exactly as before. Delete the last message of what is left, and keep going, and you reduce any protocol whatever to one that sends no riders at all. Jim Gray named it the Two Generals Paradox in 1978, in a set of lecture notes on database operating systems. The name stuck to a result that was three years old.

The proof kills certainty, for any number of exchanges, on any channel that can lose anything. Anyone who has operated a network knows what that leaves. A participant that has gone quiet and a participant that has stopped existing send you exactly the same thing, which is nothing, and no amount of listening tells the two cases apart. You can wait longer, as every system built on this planet does, which is what a timeout is: a decision to call a slow neighbor dead at a moment somebody chose in advance and wrote into a configuration file.

The channel does something else to messages besides losing them. Send two messages one after the other across a network with more than one route through it, and the second can land first, because the routes are of different lengths and one of them had a queue on it. The receiver writes them down as they arrive. That is a real order, it belongs to the receiver, and the sending order it fails to match is held by somebody else on another hill. There is a fix, and everyone who builds these systems uses it: put a number in each message before sending it. The number works because both parties agreed on the numbering in advance, which makes it a piece of freight the message is carrying rather than a fact the receiver discovered by looking at the sky.

Every requirement chapter three put on an observer survives that. A bounded view, somewhere to write, and a reading that makes a difference, with the record lasting from the writing of it to the reading of it. Read the fifth one again. It asks that the thing outlive its own operations, which is a statement about order. It never mentions a duration, a rate or a second.

## Water, pasta, cheese

A recipe carries the whole distinction. Most people have run one a thousand times without needing a name for it.

Five steps. Fill the pot with water. Bring it to the boil. Put the pasta in. Drain the pasta. Grate the cheese.

Four of those steps are wired together and one is loose. Boiling has to come after filling, since there is nothing to boil otherwise. The pasta goes in after the water boils, and putting it into cold water is a different act with a different result. Draining comes after the pasta goes in. The cheese can be grated at any point, before the pot comes out or while the pasta cooks or the previous evening. No arrangement of the kitchen makes that choice matter.

In how many different orders can this be cooked? The four wired steps have to happen in their one sequence. Grating can be dropped into any of five positions: before filling, between filling and boiling, between boiling and the pasta, between the pasta and draining, or after draining. Five orders, one dinner, and nothing distinguishes them afterward.

Draw it. Four dots in a row with arrows running along them, filling to boiling to pasta to draining, and a fifth dot sitting off to one side with no arrow touching it at all. There is no left-to-right axis in that picture and no place to put one, because putting one in means deciding where the fifth dot goes. Nothing in the kitchen decides that.

The structure underneath is a relation that says, of some pairs of steps, that one must come after the other. The relation has to behave: no step comes after itself, and if grating comes after boiling and boiling comes after filling then grating comes after filling. Beyond that it is permitted to be silent, and about grating and boiling it says nothing at all. A relation like that is a **partial order**. When it is silent about no pair, when every two items are related one way or the other, it is a **total order**.

A printed recipe is a total order. Step one, step two, step three, down the card, and somewhere between the cooking and the cardboard the one choice that never mattered got made and printed as though it had. A card has to be read in some order, being a card, so an order was invented for it.

Leslie Lamport wrote the distributed-systems version in July 1978, in the *Communications of the ACM*, in a paper called "Time, Clocks, and the Ordering of Events in a Distributed System" that runs from page 558 to page 565. He defined one event as happening before another when the first could have influenced the second, showed that this relation is a partial order and nothing stronger, and then gave an algorithm for extending it to a total order so that a computer system could number things. The algorithm needs an arbitrary rule for breaking ties between events that the partial order leaves unrelated. Any rule will serve. Choose a different one and you get a different total order, equally correct, over the same events.

Causality is the partial order. A timeline is a total order. You can have the first without the second. The first is the one the records hand out, built from which read depended on which write. That this order of reads and writes is the causal order of the physical world is a separate identification, which the geometry chapters state where they use it. The second is available to anybody willing to break ties, which can be done more than one way.

There is a physics program built around the same shape. In 1987 Luca Bombelli, Joohan Lee, David Meyer and Rafael Sorkin proposed that spacetime, looked at finely enough, is a partial order in which any two related elements have only finitely many others between them. They called it a **causal set**. Sorkin later put what the program asks of one into five words: order plus number equals geometry.

The order half rests on two theorems from the 1970s, one by Stephen Hawking, Alan King and Patrick McCarthy in 1976 and one by David Malament in 1977. Take a spacetime in which no chain of causes loops back on itself and no two distinct events share both their past and their future, and hand somebody nothing but the list of which events can influence which. From that list alone they can rebuild the light cones at every point, and with the light cones every angle and every path a ray of light can take. What the list cannot tell them is how big anything is. A region and a copy of it blown up to twice the size carry the same list. The number half supplies the size: count the elements in a region and the count is its volume.

The hard part runs the other way. Almost no partial order looks like a spacetime. Daniel Kleitman and Bruce Rothschild proved in 1975 that as the number of elements grows, almost every partial order has three layers, a bottom, a middle and a top, with no chain longer than three. A world like that would be enormously wide and last three ticks. Causal-set theory gets orders that look like spacetime by scattering points at random through a spacetime given in advance, which tests the idea and takes the spacetime as given. A theory that produces its own partial order has to say why the order it produces is one of the rare ones with a continuum in it.

The observers build an order of this kind, and the next section spells out how. Each event gets a **source height**, the number of steps in the longest chain of dependencies ending at it. Height goes up along every dependency and comes from who read what, never from anybody's position in an execution list. It counts steps and carries no unit, so it orders events without measuring any duration between them.

The record grammar is generous. Every finite partial order can be written as a pattern of authenticated reads, so records alone force neither single-file events nor a spacetime. Chapters fifteen and seventeen investigate supplied geometric constructions and continuum limits. Their order/count and calibration premises must hold before read dependencies can be identified with physical light cones, volume or proper time.

## Caused, without earlier

Putting the pasta into unboiled water is a different event from putting it into boiling water. The arrow from boiling to pasta records that dependency and nothing else. It says the second act consumes what the first produced. Nobody consulted a clock to establish it. The arrow would be there in a kitchen with no clock in it.

Carry that across to observers. One record depends on another when its accepted commit certifies a read of a version the other commit actually wrote. The writer mark must name that parent, and the value read must equal the value the parent left. Assemble the authenticated dependencies and their chains and you have the partial order, which is the entire causal structure available. Two records with no chain between them have no fact about which came first. An executor that happened to process one of them first adds nothing, and no clock, retry counter or queue position can supply the missing edge afterwards.

## Eight things called time

Values in a program have types, meaning a declaration of what kind of thing each value is, and compilers convert quietly between types they consider close enough. Write a whole number where a decimal is expected and most languages hand you the decimal without comment. The convenience is real, and so is the class of defect where two quantities that were never the same kind of thing get added together because the machine was being helpful.

Eight typed objects enter this discussion, held apart by exactly that mechanism, run in reverse. They are not eight physical clocks: closure of a complete history and a map from a sequence into spacetime are not clocks at all. The shared six-role distinction groups causal order with its logical labels, and operational clocks with proper duration, while also separating host wall time, solver counters, model coordinates and modular parameters. The more detailed list below keeps the different mathematical objects visible within that distinction.

There is closure of the complete history under its self-reference map. There is the repair order, a list of solver moves and the positions they occupy. There is an observer's record order, which says which records depend on which. There is a modular parameter attached to an algebra-state pair, carrying no unit or duration by itself. There is a worldline, a map from records to events preserving the relevant order. There is a clock readout, a real number attached to each record that increases along that order. There is proper time, the calibrated duration along a physical worldline. And there is a global time function increasing along the model's causal order. Host wall time measures a computation's resource use and belongs to a further level of description.

The eighth needs a distinction. Every finite dependency order admits a numbering that increases along its arrows, and Lamport's construction supplies useful logical clocks. Suitable spacetime models also admit global time coordinates. Neither construction selects a preferred physical clock, a universal simultaneity, or one rate all observers must share. Local records do not supply those stronger claims, and the absence of a preferred clock does not forbid a global ordering used for calculation.

"That repair took about a nanosecond" reads a position in a list as a duration. Accept it once and the positions can be added up, so a stretch of the world's history acquires an elapsed time in units, with no clock anywhere in the reasoning and no calibration performed by anybody. The sentence is humdrum, which is exactly why it gets through. Nothing about it sounds like a mistake.

Crossing between two of the eight costs an explicit map with a name on it, doing the work in the open. Chapter seventeen studies a counting construction: under the stated geometric and sampling assumptions, interval counts estimate volume, and an additional dimension-dependent map relates that volume to proper duration. A count ratio is not automatically a duration ratio. Going from a clock readout to a proper time also takes a calibration.

Take any clock readout and add seventeen to every reading. Every record reads later than the records it depends on, exactly as before, so the shifted readout is a perfectly legal clock for the same history. It is a different clock. Multiply every reading by any positive number instead and the same thing happens. The record order alone fixes neither the origin nor the rate.

## Two meetings

Einstein got to the same place from the opposite direction in 1905. He started from the speed of light being the same for every observer, which is a strong thing to be handed, and derived that two events simultaneous for one observer are not simultaneous for another. The demand was made by a light signal, and simultaneity came apart under it.

Coming the other way, an observer begins with dependencies among records it can check. Those dependencies do not select a shared physical present. A logical numbering or an admissible global coordinate can be constructed with additional choices, but its existence does not make that choice a preferred simultaneity. The missing quantity is the preferred physical present, not every possible common ordering.

Einstein's route needs the speed of light to come out the same for inertial observers, a measured physical input to that argument. A bare grammar of records does not yet derive that fact or relativistic geometry. It supplies dependency orders; a particular physical model must explain signal propagation, operational clocks and their relations. Independent records need no canonical interleaving, even when an executor or coordinate convention chooses one.

Which is why the only numbers Hafele and Keating could bring home were differences taken at a meeting. Two clocks compare where they meet. Each trip gave the flying clocks two meetings with the bench clocks, one on the way out of the Naval Observatory and one on the way back in. Everything between those two events belongs to a route rather than to a moment. The two routes, east and west, gave answers about a third of a microsecond apart.

## Sixteen schedules

Go back to the twelve observers of chapter six, wired as an icosahedron, thirty overlaps between them, and every one of the two thousand and forty-eight states run under sixteen different repair orders.

Calling them schedules was generous. The sixteen are sixteen lists, each naming which repair move occupies which position, with no duration attached to any position and no rate at which positions go by. Two of those lists that differ only by swapping a pair of moves touching no overlap in common are not two histories of the world. They are one partial order, written down twice, and the ordering between those two moves is the numbering on the recipe card.

That is why all sixteen had to be run. Taking turns properly has nothing behind it here, so the agreement the twelve reach cannot be an artifact of a well-behaved order. The same answer comes out of all sixteen lists.

Which leaves the twelve of them in the position the generals were in, with better equipment and no better prospects. Each sees five neighbors and nothing else. Each holds an order over its own records and over nobody else's. A message from a neighbor lands after an unbounded number of intervening moves, out of the sequence it was sent in, and sometimes it does not land at all. A leader would have to be recognized as one by observers with no way of comparing what they hold, and a vote would have to be taken in a present that none of them share.

They converge anyway. Every state, every order, one answer, arrived at under exactly the conditions that a branch of engineering has been building expensive machinery to survive since 1975.
