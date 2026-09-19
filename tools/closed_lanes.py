"""Closed GitHub lane issues and the live issues that own their open work.

Each successor list follows the closing banner of the lane issue. Claim
registry gates and gravity-ladder owners name the issues that own open work,
so the validators reject these numbers there. The observation ledger and the
instrument register keep a closed lane as the composition lane of a row and
render its successors beside it.

Add a lane here when it closes, with the successors its closing banner names.
"""

from __future__ import annotations

CLOSED_LANE_SUCCESSORS: dict[int, tuple[int, ...]] = {
    # Spacetime adequacy, superseded in the 2026-09-10 V4 review: common-world
    # integration and the V4 read-law production own its open work.
    728: (740, 777),
    # Mechanics and thermodynamics, closed as bounded conditional milestones:
    # source-side premise work sits in the premise discharge queue.
    731: (739,),
    732: (739,),
    # Electromagnetism, closed as a bounded finite milestone: the
    # same-history Maxwell episode continues it.
    733: (754,),
    # Standard Model structure, closed as a bounded milestone: physical
    # global form and attachment sit in common-world integration.
    734: (740,),
    # Simulation instruments, superseded: the remaining instrument question
    # is the support-wiring readout.
    737: (776,),
}


def successors(lane: int) -> tuple[int, ...]:
    """Return the live owners of a closed lane, or the lane itself if open."""
    return CLOSED_LANE_SUCCESSORS.get(lane, (lane,))
