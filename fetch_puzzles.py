#!/usr/bin/env python3
"""Pull Eterna puzzle definitions into a CSV that check_design.py can score.

  fetch_puzzles.py 3058243 2447780 ... --out puzzles.csv
  fetch_puzzles.py --lab 16022932 --out puzzles.csv

Reads the public Eterna API (no login needed); Python 3 standard library only.
Columns written: puzzle_id, title, folder, structure, constraints, locks, begin_sequence.
"""
import argparse, csv, json, sys, time, urllib.request

API = 'https://eternagame.org/get/?type={}&nid={}'


def get(kind, nid):
    with urllib.request.urlopen(API.format(kind, nid), timeout=60) as r:
        return json.load(r)['data']


def lab_puzzle_ids(lab):
    """Every puzzle in a lab/project, in the order the lab page lists them."""
    rounds = get('project', lab)['lab']['puzzles']
    if isinstance(rounds, str):                           # some labs return this JSON-encoded
        rounds = json.loads(rounds)
    return [str(p['nid']) for rnd in rounds for p in rnd['puzzles']]


def puzzle_row(pid):
    p = get('puzzle', pid)['puzzle']
    struct = p['secstruct']
    return {
        'puzzle_id': str(p['id']),
        'title': p['title'],
        # Vienna/Vienna1 -> ViennaRNA 1.8.5; Vienna2 -> ViennaRNA 2.x.
        # The API leaves this null for puzzles on the default engine, which is Vienna 1.
        'folder': p['folder'] or 'Vienna',
        'structure': struct,
        'constraints': p['constraints'] or '',
        # 'x' marks a base the player may not change; it keeps its begin_sequence letter
        'locks': p['locks'] or 'o' * len(struct),
        'begin_sequence': p['beginseq'] or 'A' * len(struct),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('ids', nargs='*', help='Eterna puzzle ids')
    ap.add_argument('--lab', help='lab/project id: fetch every puzzle in it')
    ap.add_argument('--out', help='output CSV (default: stdout)')
    a = ap.parse_args()

    ids = list(a.ids) + (lab_puzzle_ids(a.lab) if a.lab else [])
    if not ids:
        ap.error('give puzzle ids, --lab, or both')

    rows = []
    for i, pid in enumerate(ids):
        rows.append(puzzle_row(pid))
        print(f'  {i + 1}/{len(ids)} {pid}', file=sys.stderr)
        time.sleep(0.2)                                   # be kind to the server

    out = open(a.out, 'w', newline='') if a.out else sys.stdout
    w = csv.DictWriter(out, list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
    if a.out:
        out.close()
        print(f'wrote {a.out} ({len(rows)} puzzles)', file=sys.stderr)


if __name__ == '__main__':
    main()
