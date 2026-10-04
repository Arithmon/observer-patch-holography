"""Explicit, bounded metadata discovery. Never downloads outcome documents.

This command is optional and is never run by the offline verifier or CI.
It writes a new file, preserving the checked-in discovery snapshot.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen


QUERIES = (
    'date > 2026-07-28 and title tau and title mass',
    'date > 2026-07-28 and (title tau or title \u03c4 or title lepton) and title mass',
    'date > 2023-01-01 and date < 2024-01-01 and title tau and title mass',
    'date > 2023-01-01 and date < 2024-01-01 and title lepton and title mass',
)
FIELDS = ('titles', 'arxiv_eprints', 'dois', 'preprint_date', 'document_type', 'collaborations')


def fetch():
    queries = []
    for query in QUERIES:
        endpoint = 'https://inspirehep.net/api/literature?'+urlencode(
            dict(q=query, fields=','.join(FIELDS), size=25, sort='mostrecent'))
        with urlopen(endpoint, timeout=30) as response:
            raw = response.read(2_000_001)
        if len(raw) > 2_000_000:
            raise ValueError('metadata response exceeds resource contract')
        data = json.loads(raw.decode('utf-8'))
        rows = [dict(id=hit['id'], created=hit['created'],
                     **{k: hit['metadata'].get(k) for k in FIELDS}) for hit in data['hits']['hits']]
        queries.append(dict(query=query, endpoint=endpoint, total=data['hits']['total'], rows=rows))
    return dict(retrieved_utc=datetime.now(timezone.utc).isoformat(), queries=queries,
                limit=25, outcomes_requested=False,
                scope='Metadata index only. ASCII tau misses the historical Belle II title; the broader lepton query is the positive control. Index date may be later than first public release.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path, help='new snapshot path; never overwrite an existing file')
    args = parser.parse_args()
    packet = fetch()
    with args.output.open('x', encoding='utf-8', newline='\n') as output:
        json.dump(packet, output, indent=2, sort_keys=True, ensure_ascii=True)
        output.write('\n')
    print(json.dumps([(q['query'], q['total']) for q in packet['queries']], ensure_ascii=True))
