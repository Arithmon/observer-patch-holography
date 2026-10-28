# Appendix D: Run It Yourself

The prologue complained that the universe ships with no source code and no documentation. The construction this book describes ships with both. The observers, their seams and their repair law, the instruments they ask questions with, the records they keep, their clocks, the fields their clocks read and the waves those fields carry are all a program, written down in full, with every result it reports stored next to the code that produced it and a fingerprint of each result kept so that a replay can be compared against it.

The numbers in chapters nine, fifteen, twenty, twenty-one and twenty-two come out of that program: the coherence a repair leaves behind, the table of angles an observer reads off its own counts, the 245 Delft rounds, the logarithm of five fourths, the variance of sixteen between two event clocks, and the interference odds of 21/22, 91/110 and 1297/1320 for a wave driven by an incident clock. None of them needs to be taken on trust. A laptop checks them.

Get the code:

```sh
git clone https://github.com/FloatingPragma/observer-patch-holography
cd observer-patch-holography
```

The smallest check needs nothing beyond a standard Python of version 3.11 or later:

```sh
python3 -m model.observer native \
  --check model/observer/results/native.json
```

It replays 264 complete histories of one native observer, every branch of every repair included, and checks eighteen facts about them against the stored result. The last line it prints is a short report that ends with the fingerprint of the result, beginning `"all_pass": true, "checks": 18, "histories": 264`. A run of ten thousand sample histories from a fixed starting seed is one more command:

```sh
python3 -m model.observer run \
  --trials 10000 --seed 42
```

The quantum constructions, the fields, the clocks and the waves need a few numerical libraries. Install them into a private environment and run the whole verification:

```sh
python3 -m venv .venv-model
. .venv-model/bin/activate
python3 -m pip install -r model/requirements.txt
python3 -m model.verify
```

That replays every stored result and compares it with the fingerprint it was saved under. A result that comes out different by a single digit fails.

Each construction has its own folder and its own README explaining what it computes and which command checks it. The map of all of them is `model/README.md`. Chapter twenty-two's tally, field, wave and region's clock live in the folders on observer time fields, relative observer geometry and the transfer modular probe, and the formal proofs behind the theorems are the Lean files appendix C describes.

## The continuing C/E campaign

The [C/E campaign](https://github.com/FloatingPragma/observer-patch-holography/tree/main/model/completion) has its own focused replay. It checks the finite-memory observer with counted record exports, the exact joint stationary reference, and the native charge and modular-access certificates:

```sh
python3 -m pytest -q -W error model/completion
python3 tools/build_completion_program.py --check
```

The [completion plan](https://github.com/FloatingPragma/observer-patch-holography/blob/main/docs/COMPLETION_PROGRAM.md) keeps the current work in one place: the extension contract (#1135), observer services (#1136), clock operations (#1138) and the localized defect investigation (#1139). The campaign's README gives each replay and its construction contract.
