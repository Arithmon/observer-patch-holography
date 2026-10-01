"""Channel dilations with public outcomes distinct from private Kraus labels."""

import numpy as np
from m1_source_realization.model import code_transfer_kraus


def dilate(groups):
    """Return V indexed as (private branch, public outcome, output; input).

    This is a channel representation, not a new reversible full-M6 service.
    Only dilations inside accessible proper-code registers are synthesized.
    """
    if type(groups) is not list or not 0 < len(groups) <= 16:
        raise ValueError('bounded nonempty public outcome list required')
    matrices = []
    shape = None
    for group in groups:
        if type(group) is not list or not group:
            raise ValueError('nonempty Kraus group required')
        converted = []
        for raw in group:
            matrix = np.asarray(raw)
            if (matrix.ndim != 2 or matrix.dtype.kind not in 'fciu'
                    or not all(0 < d <= 16 for d in matrix.shape)
                    or not np.isfinite(matrix).all()):
                raise ValueError('finite numeric Kraus matrix required')
            if shape is None:
                shape = matrix.shape
            if matrix.shape != shape:
                raise ValueError('common input/output dimensions required')
            converted.append(matrix.astype(complex))
        matrices.append(converted)
    count = sum(map(len, matrices))
    if count > 32:
        raise ValueError('bounded environment required')
    complete = sum(k.conj().T @ k for group in matrices for k in group)
    if np.max(np.abs(complete-np.eye(shape[1]))) > 1e-12:
        raise ValueError('complete TP instrument required')
    v = np.zeros((count, len(groups), shape[0], shape[1]), complex)
    branch = 0
    for outcome, group in enumerate(matrices):
        for k in group:
            v[branch, outcome] = k
            branch += 1
    return v.reshape(-1, shape[1])


def candidate():
    rows = []
    for d in (2, 4):
        ks = code_transfer_kraus(d)
        v = dilate([[ks[0]], ks[1:]])
        rows.append(dict(code_dimension=d, environment=7-d,
                         real=v.real.tolist(), imag=v.imag.tolist()))
    return rows
