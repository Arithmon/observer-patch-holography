"""Outcome-blind nomination; authenticity of metadata needs external evidence."""
from datetime import datetime
import re


def timestamp(value):
    if type(value) is not str or re.fullmatch(r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z', value) is None:
        raise ValueError('exact UTC timestamp required')
    return datetime.strptime(value, '%Y-%m-%dT%H:%M:%SZ')


def nominate(records, registration_utc):
    start = timestamp(registration_utc)
    if type(records) is not list or len(records) > 100:
        raise ValueError('bounded metadata inventory required')
    seen, candidates, excluded = set(), [], []
    for row in records:
        if type(row) is not dict or row.keys() != {'id', 'first_public_utc', 'kind', 'observable', 'exposed_by'}:
            raise ValueError('metadata-only schema; outcomes and precision must not enter selection')
        name = row['id']
        if type(name) is not str or re.fullmatch(r'[a-z0-9][a-z0-9./_-]{0,119}', name) is None or name in seen:
            raise ValueError('unique canonical release identifier')
        seen.add(name)
        time = timestamp(row['first_public_utc'])
        if type(row['exposed_by']) is not list or not all(type(s) is str and 0 < len(s) < 120 for s in row['exposed_by']):
            raise ValueError('explicit exposure inventory required')
        if row['kind'] not in ('dedicated_measurement', 'world_average', 'theory', 'revision'):
            raise ValueError('declared document kind')
        if row['observable'] not in ('charged_tau_mass', 'other'):
            raise ValueError('declared observable')
        if time <= start or row['kind'] != 'dedicated_measurement' or row['observable'] != 'charged_tau_mass':
            excluded.append(name)
        else:
            candidates.append((time, name, row))
    candidates.sort(key=lambda x: (x[0], x[1]))
    if not candidates:
        return dict(state='NOT_READY', selected=None, reason='no post-registration dedicated tau release', excluded=sorted(excluded))
    winner = candidates[0][2]
    return dict(state='NOT_READY' if winner['exposed_by'] else 'NOMINATED_NOT_ADMITTED',
                selected=winner['id'], reason='exposed first candidate' if winner['exposed_by'] else
                'requires source-byte, dataset-overlap and uncertainty admission before evaluation',
                excluded=sorted(excluded))
