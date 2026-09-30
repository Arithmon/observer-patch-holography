"""Finite spherical support and explicit fibre refinement inventories."""

from source_selection_model.geometry import tower


def candidate():
    rows = []
    previous = None
    for stage in tower(4):
        stage['faces'] = [list(face) for face in stage['faces']]
        level, count = stage['level'], stage['vertices']
        fibres = [2+level]+[1]*(count-1)
        carriers = [[v, slot] for v in range(count) for slot in range(fibres[v])]
        if previous is None:
            mapping = None
        else:
            old = previous['carriers']
            mapping = [old.index([stage['coarsen'][v], slot] if [v, slot] in old
                                 else [stage['coarsen'][v], 0]) for v, slot in carriers]
        rows.append(dict(support=stage, fibres=fibres, carriers=carriers, coarsen=mapping,
                         section=[carriers.index([v, 0]) for v in range(count)],
                         process_links=[[i, j] for i, (v, _) in enumerate(carriers)
                                       for j, (w, _) in enumerate(carriers) if i < j and v == w]))
        previous = rows[-1]
    return rows
