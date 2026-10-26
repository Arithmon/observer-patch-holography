"""Historical lane issues, their superseding lanes, and the open issues that own open work.

Every register gate, consuming lane, owning issue and blocking issue names a
GitHub issue. The lane issues of the retired V3 and V4 programmes are all
closed; a register keeps such a number as labelled history and names the open
issue that owns the remaining work through ``live_owner``. Generated pages
render a historical number as "historical lane #n" beside its live owner.

``SUPERSEDED_LANE_SUCCESSORS`` records lanes that a later lane replaced while
the programme still ran; a claim gate names the successor lane, never the
superseded one. ``LIVE_OWNERS`` maps every historical lane that a register
still cites to the open issue that owns its remaining work. Add a lane here when
it closes, and update ``OPEN_ISSUES`` when the open physics issues change.
"""

from __future__ import annotations

# The open physics issues, by number, with their titles.
OPEN_ISSUES: dict[int, str] = {
    1025: "Define and test a discriminating consequence of the repair model",
    1026: "Construct one common microscopic repair model",
    1033: "Ongoing audit and repair of existing evidence",
    1076: "Specify the minimal repair law and test dynamical alternatives",
    1077: "Complete one source-coded preparation and retained-record history",
    1078: "Search for programme-wide OPH obstructions and stopping criteria",
}

# The document that states the open construction work in prose.
LIVE_OWNER_DOCUMENT = "model/CONSTRUCTION.md#remaining-construction-work"

# Lanes replaced by later lanes while the programme ran (historical record).
SUPERSEDED_LANE_SUCCESSORS: dict[int, tuple[int, ...]] = {
    # Spacetime adequacy, superseded in the 2026-09-10 V4 review.
    728: (740, 777),
    # Mechanics and thermodynamics, closed as bounded conditional milestones.
    731: (739,),
    732: (739,),
    # Electromagnetism, continued by the same-history Maxwell episode.
    733: (754,),
    # Standard Model structure, continued by common-world integration.
    734: (740,),
    # Simulation instruments, continued by the support-wiring readout.
    737: (776,),
    # Conserved current on a declared flat family.
    780: (735, 740, 745),
}

# Backward-compatible name used by the validators.
CLOSED_LANE_SUCCESSORS = SUPERSEDED_LANE_SUCCESSORS

# Every closed lane or issue that a register cites, mapped to the open issue
# that owns its remaining work. Physical attachments, sector completions and
# the constants belong to the common model (#1026); discriminating tests and
# freezes to #1025; repair-law, scheduler and provenance derivations to #1076;
# source-coded histories and their readouts to #1077.
LIVE_OWNERS: dict[int, int] = {
    314: 1026,  # physical matter and Spin action proof packet
    425: 1026,  # hadron backend
    506: 1025,  # independent same-scheme alpha/HVP falsification
    545: 1026,  # fine-structure endpoint reduction
    546: 1026,  # charged-lepton Yukawa landing
    554: 1026,  # discrete structural selectors
    565: 1026,  # twelve-port unit splitting and pairing
    566: 1026,  # twelve-port current from source tomography
    607: 1025,  # frozen-prediction register establishment
    618: 1026,  # gravity premise-elimination ladder
    625: 1026,  # integer load fibre and quadratic readback
    628: 1026,  # operational integer-counting readback
    643: 1025,  # A5 angular test
    649: 1026,  # N reserve attachment
    655: 1025,  # FZ-11 physical bridge
    664: 1025,  # carrier scale and positive lower bound
    666: 1025,  # seam-current discriminator
    667: 1025,  # FZ-12 comparison
    688: 1026,  # thermodynamic physical receipts
    692: 1026,  # finite causal observer net
    708: 1025,  # H0 preflight
    712: 1026,  # source-linked observer quantum-event tower
    726: 1026,  # observation ledger establishment
    727: 1026,  # premise register establishment
    728: 1026,  # spacetime adequacy
    729: 1026,  # gravitation
    730: 1026,  # quantum law and physical readout
    731: 1026,  # mechanics
    732: 1026,  # thermodynamics
    733: 1026,  # electromagnetism
    734: 1026,  # Standard Model structure
    735: 1026,  # matter-action coefficients
    736: 1026,  # masses and constants: P, N and the particle spectrum
    737: 1077,  # simulation instruments
    739: 1026,  # premise discharge and axiomatization queue
    740: 1026,  # one fixed model and cross-observable comparisons
    742: 1026,  # cosmology and astrophysics adequacy
    743: 1026,  # interactions, running and scattering
    744: 1026,  # QCD, hadrons, nuclear and atomic physics
    745: 1026,  # EWSB, weak interactions, flavour and neutrinos
    750: 1076,  # scheduler and source hypotheses from the repair process
    751: 1026,  # dark sector
    754: 1077,  # finite observer Maxwell episode on one history
    763: 1076,  # event precedence from source provenance
    776: 1077,  # support-wiring provenance readout
    777: 1077,  # the read law as a production
    778: 1077,  # population selection of committed record states
    779: 1076,  # one source law and its geometry response
    780: 1026,  # current readout on the sourced causal poset
    782: 1026,  # count-volume and ordering-fraction limits
    1087: 1078,  # audit sub-issues of the obstruction search
    1088: 1078,
    1089: 1078,
    1090: 1078,
}


def successors(lane: int) -> tuple[int, ...]:
    """Return the successor lanes of a superseded lane, or the lane itself."""
    return SUPERSEDED_LANE_SUCCESSORS.get(lane, (lane,))


def is_open(issue: int) -> bool:
    return issue in OPEN_ISSUES


def live_owner(issue: int) -> int:
    """Return the open issue that owns the work a cited issue scoped."""
    if issue in OPEN_ISSUES:
        return issue
    if issue in LIVE_OWNERS:
        return LIVE_OWNERS[issue]
    raise KeyError(f"issue #{issue} is neither open nor a registered historical lane")


def live_owners(issues) -> tuple[int, ...]:
    """Return the distinct live owners of several cited issues, in ascending order."""
    return tuple(sorted({live_owner(issue) for issue in issues}))


def is_known(issue: int) -> bool:
    return issue in OPEN_ISSUES or issue in LIVE_OWNERS


ISSUE_URL = "https://github.com/FloatingPragma/observer-patch-holography/issues"


def issue_link(issue: int) -> str:
    return f"[#{issue}]({ISSUE_URL}/{issue})"


def owner_text(issues) -> str:
    """Render live owners as links, with cited historical lanes as labelled history.

    A historical number is printed as plain text, never as a link, so a reader
    never mistakes a closed lane for a live owner.
    """
    issues = list(issues)
    if not issues:
        return ""
    rendered = ", ".join(issue_link(owner) for owner in live_owners(issues))
    historical = [issue for issue in issues if not is_open(issue)]
    if historical:
        noun = "lanes" if len(historical) > 1 else "lane"
        listed = ", ".join(f"#{issue}" for issue in historical)
        rendered += f" (historical {noun} {listed}, closed)"
    return rendered
