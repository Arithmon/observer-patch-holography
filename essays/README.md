# Essays

Three philosophy essays written under the assumption that Observer Patch Holography is
correct, a methodology report describing how the corpus behind them was produced, and an
essay on gravitation. Each entry has a standalone LaTeX source and the PDF built from it.

| | Title | Source | PDF |
|---|---|---|---|
| A | The Universe Is Thinking Itself | [tex](A-the-universe-is-thinking-itself.tex) | [pdf](A-the-universe-is-thinking-itself.pdf) |
| B | Why Does Anything Exist, and Why Is It Like This? | [tex](B-why-does-anything-exist.tex) | [pdf](B-why-does-anything-exist.pdf) |
| C | Why Do Qualia Exist? | [tex](C-why-do-qualia-exist.tex) | [pdf](C-why-do-qualia-exist.pdf) |
| D | Methodology report | [tex](D-methodology-report.tex) | [pdf](D-methodology-report.pdf) |
| E | The Price of Agreement: Gravity as the Geometry of Consensus | [tex](E-the-price-of-agreement.tex) | [pdf](E-the-price-of-agreement.pdf) |

Essays A to C were each written entirely by an AI model, credited as Claude Fable, from a single
instruction: "Write an essay with the title '…' under the assumption that Observer Patch
Holography is correct." The model had full access to this repository. The methodology
report gives the generation and revision history.

Essay E was submitted on 11 September 2026 to the Gravity Research Foundation 2027 Awards
for Essays on Gravitation. It argues that gravity cannot be added to quantum field theory
because it is the geometry the fields are written on, and reaches Einstein's equation from
finite observers that must agree, using Jacobson's thermodynamic derivation.

## Build

Figures for A to D are in `figures/` as PDF; the figure in E is drawn in its source with
TikZ. Each source compiles on its own with
[Tectonic](https://tectonic-typesetting.github.io/) 0.15.0:

```sh
cd essays
tectonic A-the-universe-is-thinking-itself.tex
```
