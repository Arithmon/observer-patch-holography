# Reverse Engineering Reality

*An Observer-First Theory of Everything.*

Forty-one chapters, a prologue, an epilogue and four appendices, in 48 files.

## Reading order

The two-digit filename prefix is the reading order and the only ordering rule. Sorting
`*.md` by name gives front matter, prologue, chapters 1 to 41, epilogue, appendices A to
D. `tools/build_book_pdf.py` does exactly that and rejects a file with no prefix or a
duplicated one.

The ten parts are a property of the chapter ranges rather than of the directory:

| Part | Chapters | Title |
|---|---|---|
| One | 1 to 6 | The Only Place to Stand |
| Two | 7 to 12 | The Protocol |
| Three | 13 to 16 | The Machine |
| Four | 17 to 22 | What an Observer Finds |
| Five | 23 to 25 | Why Gravity Is Different |
| Six | 26 to 28 | The Parts List |
| Seven | 29 to 31 | The Numbers and the One Line |
| Eight | 32 to 35 | Coming Back Down |
| Nine | 36 to 37 | The Loop Closes |
| Ten | 38 to 41 | Being One of These |

## Building

    python3 tools/build_book_pdf.py

Produces `book/reverse-engineering-reality-book.pdf` on a 6x9 trade trim through pandoc
and tectonic.

## Conventions

Paragraphs are single lines in the source, separated by blank lines, and wrapping is left
to the renderer. Straight quotes, American spelling. Bold marks a key term at first use
and is not used for emphasis. There are no em-dashes anywhere, by rule.


## Interpretation and source scope

The front matter and glossary distinguish complete-history closure, intrinsic local computation, model-specific repair, observer interfaces and calibrated clock readings. Read the maintained [OPH interpretation and working terms](../docs/INTERPRETATION_AND_TERMS.md) alongside the numerical examples. The book writes from inside the theory; the papers and the model documents state which graphs, equations and calibrations are supplied, and the book does not establish a common physically correct microscopic law. Statements about experience and ethics rest on further hypotheses that the papers state.

[Appendix D](47-appendix-d.md#the-continuing-ce-campaign) includes the focused C/E campaign replay and navigation to the current [completion plan](../docs/COMPLETION_PROGRAM.md).
