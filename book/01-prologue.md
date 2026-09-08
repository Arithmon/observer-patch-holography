# No Source Code, No Documentation

In 1632 Galileo published a book in the form of a conversation between three men and gave the best passage in it to the one who agreed with him. Shut yourself below decks on a large ship, Salviati says. Bring flies and butterflies. Bring a bowl with fish in it. Hang a bottle so that it drips into a wide vessel underneath. Watch how everything behaves while the ship is tied up at the dock, and take your time about it.

Then have the ship put to sea at whatever speed you like, provided the motion is smooth and does not rock you about.

The flies do not pile up against the stern wall. The fish swim as easily toward the front of the bowl as toward the back, without having to work harder in either direction. The drops fall into the vessel below and not one of them lands behind it. You can jump the same distance forward as back. Throwing something to a friend takes the same effort whichever way you face. And from none of these effects, Galileo writes, can you determine whether the ship is moving or standing at rest.

That experiment is harder than it looks, and different in kind from what everybody means by testing something. To test a thing from outside you at least get to stand next to it and poke it. Galileo's man is inside the thing he wants to measure. He is a part of it. What he is asking is whether a container gives any sign of itself to what it contains. The answer he gets back is that this one does not.

Negative results of that kind are the most informative things in physics. He has not found that the measurement is delicate. He has found that there is nothing to measure: the information does not exist inside the cabin, no quantity of butterflies will bring it into existence, and the four centuries since have consisted largely of people discovering that the same is true of the cabin they are actually in.

The book was on the Index of forbidden works within a year and its author spent the rest of his life under house arrest, which is a fair measure of how seriously the argument was taken.

Everything that follows is that experiment, run for longer with better equipment.

## Six symptoms

What follows is six behaviors of the machine, for now without explanations which take the rest of the book.

**Galileo's cabin has never been beaten.** Four hundred years of trying, with instruments he could not have imagined. No experiment performed inside a smoothly moving laboratory has ever revealed that it was moving. The nineteenth century made a particular effort, looking for the laboratory's motion through the medium that light was supposed to be waving in, and came up empty to a precision that embarrassed everyone involved. The information is missing with a thoroughness that looks less like a limitation of instruments and more like a decision somebody made.

**Clocks disagree. The disagreement is engineered around.** The satellites your phone uses to find itself carry atomic clocks, wrong in two directions at once. They move quickly, which makes them tick slow by about seven microseconds a day. They sit high up where gravity is weaker, which makes them tick fast by about forty-five. The two effects do not cancel. Without a correction the position your phone reports would drift off by something like ten kilometers a day.

When the first atomic clock of this kind went into orbit in 1977, the correction was built with a switch on it, because not everyone involved was convinced the physicists were right. It launched switched off. They watched the clock run fast by very close to the predicted amount for about three weeks, and then turned it on.

**Asking a question changes the answer.** In 1909 a twenty-three-year-old at Cambridge named Geoffrey Taylor set a needle in front of a gas flame and photographed the pattern of light and dark bands around its shadow. Then he began dimming the flame with sheets of smoked glass. He dimmed it until the light crossing his apparatus was so feeble that particles arrived essentially one at a time, with nothing else present for them to interfere with. To collect enough of them he had to leave the plate in place for three months, which is a long time to wait to find out that nothing has changed. The bands came out looking exactly as they had at full brightness. Whatever each particle was interfering with, it was not the others, because there were no others.

The second half of that result is even stranger. Put a detector in the apparatus that can tell you which side of the needle each particle went, and the bands disappear. Delicacy does not rescue them. Any detector that succeeds in telling you the route destroys the pattern. Any detector that fails to tell you leaves the pattern intact, and how gently it fails makes no difference at all.

**Two things far apart agree too often.** In 2015 a team at Delft trapped a single electron in each of two laboratories 1.3 kilometers apart and measured both. The answers matched more often than they could have if each electron had been carrying its answer around with it. The laboratories were separated widely enough that no signal travelling at the speed of light could have crossed between them while the measurements were happening. That experiment closed the last of the escape routes physicists had been finding in this result for fifty years. The correlation is real, it cannot be used to send a message, and no account in which the particles had properties before anyone looked survives it.

**Falling into a black hole appears to destroy the record of what fell.** On 6 February 1997 Stephen Hawking and Kip Thorne signed a wager with John Preskill. Hawking and Thorne held that information swallowed by a black hole is erased from the world. Preskill held that it survives. The stake was an encyclopedia of the winner's choice, from which, as the bet specifies, information can be recovered at will. In July 2004, at a conference in Dublin, Hawking conceded and handed Preskill a baseball encyclopedia. Thorne did not concede. Hawking remarked afterwards that burning the encyclopedia and handing over the ashes would have made the losing case rather better than the book did.

The bet was worth making because both positions are unacceptable. If the record is erased, quantum mechanics breaks. If it survives, it has to escape from the one place in the universe that nothing escapes from. Three of the most capable physicists alive spent seven years on which of two impossibilities to prefer. At the end of it one of them changed his mind and one of them did not.

**Gravity puts an area ceiling on what a region can hold.** This last one is arithmetic rather than anecdote. A warehouse twice as wide in every direction holds eight times as much. A black-hole horizon behaves differently: double its radius and its entropy grows by four. More general gravitational bounds extend that lesson under their own causal conditions, though an ordinary empty region need not saturate the bound. The observer construction reaches a separate finite boundary count. Calling that count physical entropy or area requires a map to a realized screen and a calibrated exchange rate.

## A universe without a manual

People who take apart hardware for a living meet this problem constantly: a sealed device, no documentation, no source, nobody to ask, and no way to open the case without destroying what is inside. There is a discipline for it. It comes down to a few habits.

Feed it inputs and record what comes out. Change one thing at a time. Watch the timing, because timing leaks structure that the outputs conceal: a system that answers one question quickly and a similar question slowly is telling you that those questions take different routes. Above all, watch what makes it fail. A device gives away more about its architecture in the ten seconds around a crash than in a week of ordinary operation, which is why the six items above are worth more than a century of things that worked.

And a third habit. Never assume a component you have not found on the board. If you have not located the clock, you do not get to reason about what the clock does. If you have not located the memory, you do not get to assume there is a place where the state is being kept.

This book runs that discipline under one aggravating condition, which Galileo's cabin introduced and which never goes away: the engineer is a part. There is no bench to put the machine on.

A reconstruction asks how much physics follows from a smaller starting point. OPH starts with bounded systems that read their boundaries, retain records and repair disagreements with neighbors. Agreement becomes a mathematical constraint. The test is how much familiar physics follows from those constraints, with every extra assumption visible.

The everyday picture of the world is that there is a three-dimensional container, time runs inside it at one rate for everybody, objects occupy positions in it, and physics describes how the objects move.

In 2002 the programmer Joel Spolsky wrote down a rule about descriptions of that shape. All non-trivial abstractions, to some degree, are leaky. An abstraction is a simplified account that lets you work without knowing what is underneath it. It holds until the thing underneath does something the simplification has no words for, at which point it leaks. The container picture of space and time is a **leaky abstraction** and it leaks in at least six places, which you have just finished reading. What comes back through each leak is not noise. It is a wrong answer, delivered confidently, and reproducible to as many decimal places as you care to measure.

## One equation at a time

The framework has a name: Observer Patch Holography. Its case is the reach of one observer architecture. Under explicit probability assumptions, the algebra of readable questions supports quantum probabilities. Complete reversible boundary response and internal transport constrain the local symmetries of forces. A supplied matter law supports controlled motion of a real field in a continuum and, at fixed mesh, an interacting quantum state space. Scalar-density and magnetic readouts use the same field reconstruction in both descriptions; quantum probabilities come from the declared Hilbert representation. Objects that ordinary calculations introduce separately become parts of a connected argument.

The papers state the assumptions behind each result and supply proofs or calculations that can be repeated independently. Physical claims face a further demand: identify what an instrument would measure and what outcome would refute the claim. Appendix C points to that evidence. This book develops the argument in words, explaining what each step rests on.

There are equations in this book, not many, and every one of them is doing work. Each one is introduced by a sentence saying what it is about to say and followed by a sentence saying what it said. Every symbol in it is explained in words on the same page.

An equation is compressed prose. It decompresses under four questions.

What does it count? Every quantity in physics is a count of something: how many, how much, how often, how far. Find out what is being counted and you have read half the sentence.

What is being held fixed? An equation states that certain things change together. Nothing changes together unless something else is pinned in place. The pinned quantity is usually the interesting one, and usually not written down.

What is allowed to vary, and what happens at the edges? Set a term to zero and see what the equation collapses into. The collapsed version is generally a law you learned at school. Watching your old law fall out of the new one is how you check that you have understood either of them.

What do the units say? If one side of the equals sign is a length and the other is a time, then a speed is hiding inside the equals sign, or somebody has made a mistake. In 1999 NASA lost the Mars Climate Orbiter because one team's software reported thrust in pounds and the other team's software read it as newtons. The spacecraft went into the Martian atmosphere instead of around it. Units are the cheapest error detection ever invented.

## One line

There is a single line that working physicists write down when somebody asks what the world is made of. It is short enough to fit on a coffee mug, which is where most of the people who have seen it encountered it. Four or five terms, depending on how you group them. Each term is a compression: unfold one and out comes a force, or a family of particles, or the reason anything has mass.

The Standard Model packages an extraordinary range of measurements in that line. Its symmetry and field choices constrain the terms, while experiments supply parameters. A deeper reconstruction must account for those choices as well as reproduce their successful predictions.

By the last third of this book, you will be able to read that line and ask which parts follow from observer consistency, which require further assumptions, and how the resulting claims can be tested.

That is the promise. This prologue makes no other.
