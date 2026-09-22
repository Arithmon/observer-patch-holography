#!/usr/bin/env python3
"""Check Lean comment style and the in-file kernel-axiom audit blocks.

Two independent sections run over the Lean sources under ``Lean/``.

**Section 1, comment prose.**  ``docs/STYLE_GUIDE.md`` binds code comments to
the same two constraints it binds prose to: state the exact research status,
and avoid machine-generated prose mannerisms.  ``tools/check_reader_style.py``
enforces that vocabulary on Markdown and TeX surfaces.  A Lean module docstring
is read by the same audience, so the same vocabulary applies, with one
mechanical difference: a Lean file is mostly code, and code must never trip a
prose rule.  A declaration called ``nowhere_dense``, a hypothesis called
``hstill``, a field called ``already_reduced``: all of those are identifiers,
and none of them is prose.  The scanner therefore masks every non-comment
character before matching, so only ``--``, ``/- -/``, ``/-- -/`` and ``/-! -/``
text can produce a violation.

**Section 2, kernel-axiom audit coverage.**  Most Lean modules end with a block
of ``#print axioms`` commands, so the build log carries the kernel dependency of
every theorem the module states.  A block that names only some of a module's
theorems is worse than no block: it reads as a complete audit.  Section 2 keeps
a registry of modules whose block covers every public theorem and lemma they
state, and fails when one of those blocks loses coverage.  Modules outside the
registry are reported by ``--axiom-coverage-report`` with their unaudited count
and are not failed, so the registry grows by deliberate registration.

Usage::

  python3 tools/check_lean_docstring_style.py
  python3 tools/check_lean_docstring_style.py --axiom-coverage-report
"""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LEAN_ROOT = ROOT / "Lean"


# ---------------------------------------------------------------------------
# Byte-pinned modules
# ---------------------------------------------------------------------------

PINNED_BYTES_REASON = (
    "A committed artifact pins this module's bytes: a replay receipt, an "
    "independent certificate, a frozen-prediction custody manifest, or a "
    "source projection records the file digest, so rewording a comment would "
    "move a digest that a release has published. The prose rules of section 1 "
    "are therefore skipped for these modules, and the wording is corrected "
    "when the owning artifact is regenerated for an independent reason."
)

# Every module whose HEAD bytes a committed digest covers and whose comments
# carry a section-1 hit. Checked against the digest set rather than guessed:
# each entry has its sha256 quoted in at least one tracked receipt, manifest,
# certificate, custody record, or source projection outside ``Lean/``.
PINNED_BYTES_MODULES = [
    "Lean/EventAlgebra/SourceBoundInstrumentInterface.lean",
    "Lean/ObserverPatchHolography/EinsteinBranch/FixedCapacityWLaw.lean",
    "Lean/ObserverPatchHolography/ScalarSeamRepair.lean",
    "Lean/Screen/A2HolonomyBridge.lean",
    "Lean/Screen/A5Commutant.lean",
    "Lean/Screen/A5PortModule.lean",
    "Lean/Screen/ConeCochainBridge.lean",
    "Lean/Screen/DiscreteCoulombGreen.lean",
    "Lean/Screen/KineticFamilyCancellation.lean",
    "Lean/Screen/LocalFaceMaxwellAction.lean",
    "Lean/Screen/OrientedFaceBracketSelector.lean",
    "Lean/Screen/SeamCurrentCarrierQuotient.lean",
    "Lean/Screen/SeamCurrentPhotonLeptonThreshold.lean",
    "Lean/Screen/SeamMaxwellContinuum.lean",
    "Lean/Screen/WhitneyMaxwellDynamics.lean",
    "Lean/Thermodynamics/FiniteConditionalRepair.lean",
    "Lean/Thermodynamics/GreenKubo.lean",
]

# Line-level allowances inside a module that section 1 does check. The intended
# use is a frozen quotation: a sentence whose exact words a committed artifact
# reproduces, sitting in a module that is otherwise free to be reworded. Each
# entry is ``(module path, 1-based line, label)`` and suppresses exactly that
# one match on that one line. The list is expected to stay empty: a module that
# needs more than a line or two belongs in PINNED_BYTES_MODULES instead.
PROSE_EXEMPTIONS: list[tuple[str, int, str]] = []


# ---------------------------------------------------------------------------
# Section 1 vocabulary
# ---------------------------------------------------------------------------

# Progress narration. Comments state what a module proves, not what the project
# did to it. The list mirrors PROGRESS_PATTERNS in tools/check_reader_style.py.
#
# "used to" is deliberately absent, although docs/STYLE_GUIDE.md bans it in
# prose. In Lean comments the phrase is almost always the instrumental sense:
# "the lemma used to close the goal", "the tactic used to discharge the side
# condition". Matching it would flag correct mathematical prose on nearly every
# module, and the historical sense ("this used to be an axiom") is already
# caught by "previously", "formerly", "an earlier version" and "no longer".
PROGRESS_PATTERNS: list[tuple[re.Pattern[str], str]] = [
    (re.compile(r"\bnow\b", re.IGNORECASE), "progress word: now"),
    (re.compile(r"\balready\b", re.IGNORECASE), "progress word: already"),
    (re.compile(r"\bcurrently\b", re.IGNORECASE), "progress word: currently"),
    (re.compile(r"\bpresently\b", re.IGNORECASE), "progress word: presently"),
    (re.compile(r"\blatest\b", re.IGNORECASE), "progress word: latest"),
    (re.compile(r"\brecent(?:ly)?\b", re.IGNORECASE), "progress word: recent/recently"),
    (re.compile(r"\bnewly\b", re.IGNORECASE), "progress word: newly"),
    (re.compile(r"\bpreviously\b", re.IGNORECASE), "progress word: previously"),
    (re.compile(r"\bformerly\b", re.IGNORECASE), "progress word: formerly"),
    (re.compile(r"\bno longer\b", re.IGNORECASE), "progress phrase: no longer"),
    (re.compile(r"\bstill\b", re.IGNORECASE), "progress word: still"),
    (
        re.compile(r"\b(?:has|have|had)?\s*not yet\b", re.IGNORECASE),
        "progress phrase: not yet",
    ),
    (
        re.compile(
            r"\b(?:attachment|bridge|construction|derivation|identification|"
            r"map|obligation|premise|question|receipt|route|selection|theorem|"
            r"verdict|work)s?\s+remains?\s+open\b",
            re.IGNORECASE,
        ),
        "progress phrase: remains open",
    ),
    (re.compile(r"\bso far\b", re.IGNORECASE), "progress phrase: so far"),
    (re.compile(r"\bfuture work\b", re.IGNORECASE), "progress phrase: future work"),
    (re.compile(r"\bnext step\b", re.IGNORECASE), "progress phrase: next step"),
    (re.compile(r"\bcurrent state\b", re.IGNORECASE), "progress phrase: current state"),
    (re.compile(r"\bwhat changed\b", re.IGNORECASE), "progress phrase: what changed"),
    (
        re.compile(r"\bwhat did not change\b", re.IGNORECASE),
        "progress phrase: what did not change",
    ),
    (re.compile(r"\bgoing forward\b", re.IGNORECASE), "progress phrase: going forward"),
    (re.compile(r"\bin the future\b", re.IGNORECASE), "progress phrase: in the future"),
    (re.compile(r"\bwill be added\b", re.IGNORECASE), "progress phrase: will be added"),
    (
        re.compile(r"\ban earlier version\b", re.IGNORECASE),
        "progress phrase: an earlier version",
    ),
    (
        re.compile(r"\bhas been updated\b", re.IGNORECASE),
        "progress phrase: has been updated",
    ),
    (
        re.compile(r"\bsupersed(?:e|es|ed|ing)\b", re.IGNORECASE),
        "progress word: supersede",
    ),
]

# Machine-prose tells. The em dash leads the list: it is the single most common
# tell in generated mathematical prose, and Lean comments collect it faster than
# any other surface because a dash reads as a natural aside between clauses.
AI_TELL_PATTERNS: list[tuple[re.Pattern[str], str]] = [
    (re.compile(r"—"), "em dash"),
    (re.compile(r"\bnot\s+only\b", re.IGNORECASE), "not-only contrast"),
    (re.compile(r"\bnot\s+just\b", re.IGNORECASE), "not-just contrast"),
    (re.compile(r"\bnot\s+merely\b", re.IGNORECASE), "not-merely contrast"),
    (re.compile(r"\bnot\s+simply\b", re.IGNORECASE), "not-simply contrast"),
    # The window stops at sentence punctuation and at a line break, so a
    # negation and a "but" in different wrapped sentences do not pair up.
    (re.compile(r"\bnot\b[^.!?\n]{0,120}\bbut\b", re.IGNORECASE), "not-X-but-Y contrast"),
    (re.compile(r"\bdelve\b", re.IGNORECASE), "AI-tell word: delve"),
    (re.compile(r"\btapestry\b", re.IGNORECASE), "AI-tell word: tapestry"),
    (re.compile(r"\bseamless(?:ly)?\b", re.IGNORECASE), "AI-tell word: seamless"),
    (re.compile(r"\bunpack\b", re.IGNORECASE), "AI-tell word: unpack"),
    (
        re.compile(r"\bIt is important to note\b", re.IGNORECASE),
        "boilerplate intro",
    ),
    (re.compile(r"\bAt its core\b", re.IGNORECASE), "boilerplate intro"),
    (
        re.compile(r"\b(?:In essence|Simply put|Put simply)\b", re.IGNORECASE),
        "boilerplate intro",
    ),
    (
        re.compile(r"\bIt(?: is|'s) worth noting\b", re.IGNORECASE),
        "boilerplate intro",
    ),
    (
        re.compile(
            r"(?:^|[.!?]\s+)(?:Moreover|Furthermore|Crucially|Importantly|Notably|"
            r"In conclusion|Ultimately|That said|With this in mind|"
            r"Building on this idea),",
            re.IGNORECASE | re.MULTILINE,
        ),
        "formulaic sentence opener",
    ),
    (
        re.compile(r"\bthis (?:module|section|file) (?:aims|seeks) to\b", re.IGNORECASE),
        "anthropomorphized document",
    ),
]

# The banned h-word, in any case, assembled so the word itself never appears in
# this file. docs/STYLE_GUIDE.md: the word never appears in any surface.
H_WORD_PATTERN: list[tuple[re.Pattern[str], str]] = [
    (re.compile(r"\b" + "hon" + "est" + r"\w*\b", re.IGNORECASE), "banned h-word"),
]

PROSE_PATTERNS = PROGRESS_PATTERNS + AI_TELL_PATTERNS + H_WORD_PATTERN


# ---------------------------------------------------------------------------
# Comment masking
# ---------------------------------------------------------------------------


def comment_view(text: str) -> str:
    """Return ``text`` with every non-comment character replaced by a space.

    Line breaks survive everywhere, so an offset in the result has the same
    line and column as the same offset in the source. Lean comment forms:

    * ``--`` to end of line;
    * ``/- ... -/`` block comments, which nest;
    * ``/-- ... -/`` declaration docstrings and ``/-! ... -/`` module
      docstrings, which are block comments with a marked opener and nest the
      same way.

    String literals are skipped before comment openers are considered, so a
    ``--`` or a ``/-`` inside a literal never starts a comment. Character
    literals are not treated specially: ``'`` is a legal identifier character
    in Lean, so a general rule would misread ``h'`` and ``σ'``. The narrow
    three-character forms that could matter are handled explicitly.
    """

    length = len(text)
    out = ["\n" if char == "\n" else " " for char in text]
    index = 0

    def keep(start: int, stop: int) -> None:
        for cursor in range(start, min(stop, length)):
            out[cursor] = text[cursor]

    while index < length:
        char = text[index]

        # Comment openers first: a block comment swallows string syntax.
        if char == "/" and text.startswith("/-", index):
            start = index
            depth = 0
            while index < length:
                if text.startswith("/-", index):
                    depth += 1
                    index += 2
                    continue
                if text.startswith("-/", index):
                    depth -= 1
                    index += 2
                    if depth == 0:
                        break
                    continue
                index += 1
            keep(start, index)
            continue

        if char == "-" and text.startswith("--", index):
            start = index
            while index < length and text[index] != "\n":
                index += 1
            keep(start, index)
            continue

        if char == '"':
            index += 1
            while index < length:
                if text[index] == "\\":
                    index += 2
                    continue
                if text[index] == '"':
                    index += 1
                    break
                index += 1
            continue

        # ``'-'`` and ``'"'`` are the only character literals whose payload
        # could otherwise open a comment or a string.
        if char == "'" and text[index : index + 3] in ("'-'", "'\"'"):
            index += 3
            continue

        index += 1

    return "".join(out)


def code_view(text: str) -> str:
    """Return ``text`` with every comment character replaced by a space.

    The exact complement of :func:`comment_view`, so section 2 reads
    declarations from code alone. Without it a docstring sentence such as "the
    general matrix theorem in this module" parses as a declaration named
    ``in``, and the coverage count is wrong in both directions.
    """

    comments = comment_view(text)
    out = []
    for source_char, comment_char in zip(text, comments):
        if source_char == "\n":
            out.append("\n")
        elif comment_char == " ":
            out.append(source_char)
        else:
            out.append(" ")
    return "".join(out)


# ---------------------------------------------------------------------------
# Section 2: kernel-axiom audit coverage
# ---------------------------------------------------------------------------

# Public declarations only, matching tools/check_lean_theorem_count.py: a
# leading ``private`` keeps a declaration out, same-line attributes and the
# ``protected``/``noncomputable`` modifiers do not.
DECLARATION_RE = re.compile(
    r"^[ \t]*(?:@\[[^\]\n]*\][ \t]*)*(?:(?:protected|noncomputable)[ \t]+)*"
    r"(?:theorem|lemma)[ \t]+([^\s({\[:]+)",
    re.MULTILINE,
)

PRINT_AXIOMS_RE = re.compile(r"^[ \t]*#print[ \t]+axioms[ \t]+(\S+)", re.MULTILINE)

# Modules whose ``#print axioms`` block names every public theorem and lemma in
# the file. Registration is deliberate: a module joins the list when its block
# is completed, and section 2 then keeps it complete. Two module families are
# deliberately absent although their blocks are large: Thermodynamics/GreenKubo
# and Thermodynamics/FiniteConditionalRepair appear in PINNED_BYTES_MODULES, so
# their blocks cannot grow without moving a published digest.
AXIOM_BLOCK_COMPLETE = [
    "Lean/Dynamics/ProtectedCharge.lean",
    "Lean/Dynamics/StoneConverse.lean",
    "Lean/Dynamics/WardLimitManifest.lean",
    "Lean/EventAlgebra/Basic.lean",
    "Lean/EventAlgebra/ExpectationBound.lean",
    "Lean/EventAlgebra/FiniteBornFrame.lean",
    "Lean/EventAlgebra/Lueders.lean",
    "Lean/EventAlgebra/NoBroadcastingAdapter.lean",
    "Lean/EventAlgebra/OperationalAdditivityBoundary.lean",
    "Lean/EventAlgebra/OperationalPhaseAttainment.lean",
    "Lean/EventAlgebra/PartitionAverage.lean",
    "Lean/EventAlgebra/PartitionPinching.lean",
    "Lean/EventAlgebra/PartitionPinchingCP.lean",
    "Lean/EventAlgebra/PhaseInstrumentDetermination.lean",
    "Lean/EventAlgebra/QuantumAdequacySurface.lean",
    "Lean/EventAlgebra/SchroedingerFrameFlow.lean",
    "Lean/EventAlgebra/SlotLocalitySurface.lean",
    "Lean/EventAlgebra/SourceBoundInstrumentInterface.lean",
    "Lean/EventAlgebra/StateExpectation.lean",
    "Lean/EventAlgebra/Superselection.lean",
    "Lean/EventAlgebra/Tsirelson.lean",
    "Lean/EventAlgebra/TsirelsonSaturation.lean",
    "Lean/Geometry/CommonWorldJointAction.lean",
    "Lean/Geometry/CommonWorldKinematicsWitness.lean",
    "Lean/Geometry/IntegerKCombInvariance.lean",
    "Lean/Geometry/InverseSquareShellLaw.lean",
    "Lean/Geometry/ScreenCarrierMapCandidate.lean",
    "Lean/Geometry/SourceOrderEinsteinComposition.lean",
    "Lean/Geometry/SourceOrderFrameCompatibilityPacket.lean",
    "Lean/InformationProjection/PathGibbs.lean",
    "Lean/ObservableNormalForms/ObservableNormalForms/AxiomAudit.lean",
    "Lean/ObserverPatchHolography/BridgeEquivalence.lean",
    "Lean/ObserverPatchHolography/CapacityClosurePrinciple.lean",
    "Lean/ObserverPatchHolography/CapacityFixedPoint.lean",
    "Lean/ObserverPatchHolography/ClebschRatio.lean",
    "Lean/ObserverPatchHolography/ClosurePreflight.lean",
    "Lean/ObserverPatchHolography/CoreAxioms.lean",
    "Lean/ObserverPatchHolography/EinsteinBranch/AxiomAudit.lean",
    "Lean/ObserverPatchHolography/EinsteinBranch/CollarPremiseDerivation.lean",
    "Lean/ObserverPatchHolography/EinsteinBranch/CollarScaleSaturation.lean",
    "Lean/ObserverPatchHolography/EinsteinBranch/Composition.lean",
    "Lean/ObserverPatchHolography/EinsteinBranch/Consensus.lean",
    "Lean/ObserverPatchHolography/EinsteinBranch/DeepProfileClosure.lean",
    "Lean/ObserverPatchHolography/EinsteinBranch/Entropy.lean",
    "Lean/ObserverPatchHolography/EinsteinBranch/FiniteCapGeneratorSplit.lean",
    "Lean/ObserverPatchHolography/EinsteinBranch/FixedCapacityWLaw.lean",
    "Lean/ObserverPatchHolography/EinsteinBranch/Manifest.lean",
    "Lean/ObserverPatchHolography/EinsteinBranch/PerCutCollarComposition.lean",
    "Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean",
    "Lean/ObserverPatchHolography/EinsteinBranch/SmallBall.lean",
    "Lean/ObserverPatchHolography/EinsteinBranch/SparseRecordDefectWitness.lean",
    "Lean/ObserverPatchHolography/EinsteinBranch/Tensor.lean",
    "Lean/ObserverPatchHolography/IcosahedralAntibridge.lean",
    "Lean/ObserverPatchHolography/IcosahedralAxisNoGo.lean",
    "Lean/ObserverPatchHolography/IcosahedralOrbitStabilizer.lean",
    "Lean/ObserverPatchHolography/Provenance/FiniteCausetCompiler.lean",
    "Lean/ObserverPatchHolography/RepairGapChain.lean",
    "Lean/ObserverPatchHolography/SeedPi.lean",
    "Lean/ObserverPatchHolography/YangMillsGap.lean",
    "Lean/ObserverPatchHolography/YangMillsProp81.lean",
    "Lean/QFT/CPRestrictionNet.lean",
    "Lean/QFT/EffectiveQuantumComparison.lean",
    "Lean/QFT/FiniteUnitaryScatteringNoGo.lean",
    "Lean/QFT/InhabitedTowerRow.lean",
    "Lean/QFT/InheritanceMatrix.lean",
    "Lean/QFT/InheritanceMatrixGuards.lean",
    "Lean/QFT/JointInstance.lean",
    "Lean/QFT/PathTimeSliceInterface.lean",
    "Lean/QFT/ScalarRegionalTimeSlice.lean",
    "Lean/QFT/SectorInheritance.lean",
    "Lean/QFT/SourceLinkSquare.lean",
    "Lean/QFT/StructuralInheritance.lean",
    "Lean/QFT/StructuralNetAdequacySurface.lean",
    "Lean/QFT/TowerAnchoredCorrelation.lean",
    "Lean/Screen/A2HolonomyBridge.lean",
    "Lean/Screen/A5AngularBands.lean",
    "Lean/Screen/A5AngularMultiplets.lean",
    "Lean/Screen/A5IncidenceResponse.lean",
    "Lean/Screen/A5OrbitRaySeparation.lean",
    "Lean/Screen/A5PrimitivePortPrediction.lean",
    "Lean/Screen/ActionTimeGram.lean",
    "Lean/Screen/BipoSHFrameInvariant.lean",
    "Lean/Screen/BipoSHTransferInvariant.lean",
    "Lean/Screen/CartanScalarReduction.lean",
    "Lean/Screen/CertifiedScaledStepInstrument.lean",
    "Lean/Screen/CommonEWOrderUnit.lean",
    "Lean/Screen/DiscreteGauss.lean",
    "Lean/Screen/DispersionArmingInterface.lean",
    "Lean/Screen/ElectroweakPoleScaleQuotient.lean",
    "Lean/Screen/ExteriorComponentBridge.lean",
    "Lean/Screen/GlobalFormCharacterDescent.lean",
    "Lean/Screen/KineticFamilyCancellation.lean",
    "Lean/Screen/KineticFormDichotomy.lean",
    "Lean/Screen/LayeredDiscreteGauss.lean",
    "Lean/Screen/LightSignalAdequacySurface.lean",
    "Lean/Screen/LightSignalMaxwellComposition.lean",
    "Lean/Screen/MatterGrammarIndexBridge.lean",
    "Lean/Screen/PhysicalA5ForcingNoGo.lean",
    "Lean/Screen/PortGramRepairBand.lean",
    "Lean/Screen/PrimitiveHopSelection.lean",
    "Lean/Screen/PrimitivePortDualMeasure.lean",
    "Lean/Screen/PrimitivePortMetricAttachment.lean",
    "Lean/Screen/PrimitivePortOperatorSelectionBoundary.lean",
    "Lean/Screen/PrimitivePortScaleBoundary.lean",
    "Lean/Screen/RGRepresentationFrontier.lean",
    "Lean/Screen/RegionalContinuity.lean",
    "Lean/Screen/S2DesignSignature.lean",
    "Lean/Screen/SMStructureAdequacySurface.lean",
    "Lean/Screen/SMStructureComposition.lean",
    "Lean/Screen/SeamCurrentAuxiliaryOscillatorLift.lean",
    "Lean/Screen/SeamCurrentEdge30Remainder.lean",
    "Lean/Screen/SeamCurrentPhysicalMetricAttachment.lean",
    "Lean/Screen/SeamU1HolonomyClassification.lean",
    "Lean/Screen/SourceActionTime.lean",
    "Lean/Screen/TopThreeKernelFix.lean",
    "Lean/Screen/UnitSplit12.lean",
    "Lean/Screen/VolumeReadoutBridge.lean",
    "Lean/Thermodynamics/CoherentRefinementFamily.lean",
    "Lean/Thermodynamics/FirstLawIdentity.lean",
    "Lean/Thermodynamics/FluctuationTheorems.lean",
    "Lean/Thermodynamics/FourLawAdequacySurface.lean",
    "Lean/Thermodynamics/HorizonThermalitySurface.lean",
    "Lean/Thermodynamics/PhysicalCalibrationImport.lean",
    "Lean/Time/ProperTimeCalibration.lean",
    "Lean/Time/SourceCountClockEnclosure.lean",
    "Lean/Tower/ConsensusTower.lean",
    "Lean/Tower/CumulativeCapacityEndpoint.lean",
    "Lean/Tower/EventGeometryReadout.lean",
    "Lean/Tower/FixedFederationEndpoint.lean",
    "Lean/Tower/FixedFederationExecutionEndpoint.lean",
    "Lean/Tower/FixedFederationFanoutEndpoint.lean",
    "Lean/Variational/EnrichmentCharacterization.lean",
    "Lean/Variational/FiniteHistoryBridge.lean",
    "Lean/Variational/MechanicsAdequacySurface.lean",
    "Lean/Variational/ModeExtremalEnrichment.lean",
    "Lean/Variational/SourceToHamiltonianComposed.lean",
    "Lean/Variational/StationarySaddleCoverage.lean",
    "Lean/Variational/TranslationInvariantComposedInstance.lean",
]


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------


def lean_sources() -> list[Path]:
    """Lean sources under ``Lean/`` that git tracks.

    Counting tracked files keeps the local gate equal to the CI gate, the same
    reason tools/check_lean_theorem_count.py gives. Outside a git checkout the
    walk falls back to the filesystem.
    """

    try:
        listing = subprocess.run(
            ["git", "-C", str(ROOT), "ls-files", "-z", "--", "Lean"],
            check=True,
            capture_output=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        listing = b""
    if listing:
        paths = [
            ROOT / entry.decode("utf-8")
            for entry in listing.split(b"\0")
            if entry.endswith(b".lean")
        ]
    else:
        paths = list(LEAN_ROOT.rglob("*.lean"))
    return sorted(
        path
        for path in paths
        if path.is_file()
        and not any(
            part.startswith(".") or part == "lake-packages" for part in path.parts
        )
    )


def line_col(text: str, index: int) -> tuple[int, int]:
    line = text.count("\n", 0, index) + 1
    line_start = text.rfind("\n", 0, index) + 1
    return line, index - line_start + 1


# ---------------------------------------------------------------------------
# Section 1
# ---------------------------------------------------------------------------


def prose_issues(rel: str, text: str) -> list[str]:
    """Section-1 violations in one module's comments."""

    view = comment_view(text)
    exempt = {(path, line, label) for path, line, label in PROSE_EXEMPTIONS}
    issues: list[tuple[int, int, str]] = []
    for pattern, label in PROSE_PATTERNS:
        for match in pattern.finditer(view):
            line, col = line_col(view, match.start())
            if (rel, line, label) in exempt:
                continue
            issues.append((line, col, label))
    return [f"{rel}:{line}:{col}: {label}" for line, col, label in sorted(issues)]


# ---------------------------------------------------------------------------
# Section 2
# ---------------------------------------------------------------------------


def audit_block(text: str) -> tuple[list[str], set[str]]:
    """Public declaration names and audited names, read from code alone.

    A ``#print axioms`` argument may be fully qualified or bare, so names are
    compared on their final component. Comparing full names would demand a
    namespace resolver, and the final component is unique inside one module in
    every case the corpus contains.
    """

    code = code_view(text)
    declared = [
        match.group(1).split(".")[-1] for match in DECLARATION_RE.finditer(code)
    ]
    audited = {name.split(".")[-1] for name in PRINT_AXIOMS_RE.findall(code)}
    return declared, audited


def unaudited_names(text: str) -> list[str]:
    declared, audited = audit_block(text)
    seen: set[str] = set()
    out: list[str] = []
    for name in declared:
        if name in audited or name in seen:
            continue
        seen.add(name)
        out.append(name)
    return out


def coverage_issues(rel: str, text: str) -> list[str]:
    _, audited = audit_block(text)
    if not audited:
        return [
            f"{rel}:1:1: registered for full kernel-axiom coverage but carries "
            "no #print axioms block"
        ]
    missing = unaudited_names(text)
    if not missing:
        return []
    return [
        f"{rel}:1:1: kernel-axiom block omits {len(missing)} public "
        f"declaration(s): {', '.join(missing)}"
    ]


def coverage_report(paths: list[Path]) -> int:
    rows: list[tuple[str, int, int, bool]] = []
    for path in paths:
        text = path.read_text(encoding="utf-8", errors="ignore")
        declared, audited = audit_block(text)
        if not audited:
            continue
        missing = unaudited_names(text)
        rel = path.relative_to(ROOT).as_posix()
        rows.append((rel, len(declared), len(missing), rel in set(AXIOM_BLOCK_COMPLETE)))
    registered = sum(1 for row in rows if row[3])
    complete = sum(1 for row in rows if row[2] == 0)
    print(
        f"modules carrying a kernel-axiom block: {len(rows)}; "
        f"complete: {complete}; registered: {registered}"
    )
    for rel, declared, missing, is_registered in rows:
        mark = "registered" if is_registered else "unregistered"
        print(f"- {rel}: {declared} public declaration(s), {missing} unaudited [{mark}]")
    return 0


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--axiom-coverage-report",
        action="store_true",
        help="list every module carrying a kernel-axiom block with its unaudited count",
    )
    parser.add_argument(
        "--pinned-residue",
        action="store_true",
        help=(
            "list the section-1 hits deferred inside the byte-pinned modules; "
            "they are corrected when the owning artifact is regenerated"
        ),
    )
    args = parser.parse_args(argv)

    paths = lean_sources()
    if args.axiom_coverage_report:
        return coverage_report(paths)

    pinned = set(PINNED_BYTES_MODULES)
    registry = set(AXIOM_BLOCK_COMPLETE)
    issues: list[str] = []
    residue: list[str] = []
    seen_registry: set[str] = set()

    for path in paths:
        rel = path.relative_to(ROOT).as_posix()
        text = path.read_text(encoding="utf-8", errors="ignore")
        if rel in pinned:
            residue.extend(prose_issues(rel, text))
        else:
            issues.extend(prose_issues(rel, text))
        if rel in registry:
            seen_registry.add(rel)
            issues.extend(coverage_issues(rel, text))

    for rel in sorted(registry - seen_registry):
        issues.append(f"{rel}:1:1: registered for full kernel-axiom coverage but absent")
    for rel in sorted(pinned):
        if not (ROOT / rel).is_file():
            issues.append(
                f"{rel}:1:1: listed as byte-pinned but absent. "
                f"{PINNED_BYTES_REASON}"
            )

    if args.pinned_residue:
        print(
            f"section-1 hits deferred inside {len(pinned)} byte-pinned modules: "
            f"{len(residue)}"
        )
        for hit in residue:
            print(f"- {hit}")
        return 0

    if issues:
        print("lean docstring style check failed:")
        for issue in issues:
            print(f"- {issue}")
        return 1

    print(
        f"lean docstring style check OK "
        f"({len(paths)} modules, {len(pinned)} byte-pinned carrying "
        f"{len(residue)} deferred section-1 hit(s), "
        f"{len(registry)} with a complete kernel-axiom block)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
