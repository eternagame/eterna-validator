#!/usr/bin/env python3
"""Single source of truth for whether a design solves an Eterna puzzle.

A design passes when, as in the Eterna game (EternaJS src/eterna/constraints):
  SHAPE          stock RNAfold MFE == target   (ViennaRNA 1.8.5 for folder Vienna/Vienna1,
                                                ViennaRNA 2.x for Vienna2; `RNAfold -T 37`)
  locks          locked positions ('x') equal begin_sequence
  GC / AUMAX / GUMAX       pairs of that type in the MFE  <= value
  GCMIN / AU / GU / PAIRS  pairs of that type in the MFE  >= value
  CONSECUTIVE_X  longest run of base X  < value
Any other constraint name is an error, never silently ignored.

Binaries: $RNAFOLD_V1 (default ~/Dropbox/src/ViennaRNA-1.8.5/Progs/RNAfold) and $RNAFOLD_V2
(default `RNAfold` on PATH).

Usage:
  check_design.py PUZZLES.csv DESIGNS.tsv [--target-only] [--unique-mfe] [--out RESULTS.tsv]
PUZZLES.csv needs puzzle_id, folder, structure, constraints, locks, begin_sequence
(data/eterna_hard_all.csv or data/eterna100_player_constraints.csv).
DESIGNS.tsv has columns puzzle_id, sequence (plus anything else, passed through).
--target-only checks SHAPE only (the original Eterna100 benchmark condition).
--unique-mfe additionally demands the target be the ONLY structure at the minimum free energy.
"""
import argparse, csv, os, re, subprocess, sys, tempfile
from pathlib import Path

RNAFOLD = {
    1: os.environ.get('RNAFOLD_V1', str(Path.home() / 'Dropbox/src/ViennaRNA-1.8.5/Progs/RNAfold')),
    2: os.environ.get('RNAFOLD_V2', 'RNAfold'),
}
# --unique-mfe only. Both must be ViennaRNA 2.x, even for a Vienna 1 puzzle: 1.8.5's subopt.c
# forces dangles=2 while its RNAfold uses dangles=1, so 1.8.5's own RNAsubopt enumerates the
# wrong window and silently misses structures that tie under the Vienna 1 model.
RNASUBOPT = os.environ.get('RNASUBOPT', 'RNAsubopt')
RNAEVAL = os.environ.get('RNAEVAL', 'RNAeval')
TURNER1999 = os.environ.get('TURNER1999_PAR')          # misc/rna_turner1999.par from a 2.x source tree
MAX_PAIRS = {'GC': 'GC', 'AUMAX': 'AU', 'GUMAX': 'GU'}
MIN_PAIRS = {'GCMIN': 'GC', 'AU': 'AU', 'GU': 'GU', 'PAIRS': 'ANY'}
TIE = 0.005                                            # kcal/mol; energies are reported to the cent
WINDOWS = (0.3, 1.0, 2.0, 4.0)                         # widened until the enumeration is certified
MAX_CANDS = 400_000                                    # give up rather than exhaust memory
ENERGY = re.compile(r'^([.()]+)\s+\(\s*(-?[\d.]+)\)\s*$', re.M)


def engine(folder):
    if folder in ('Vienna', 'Vienna1'):
        return 1
    if folder == 'Vienna2':
        return 2
    raise ValueError(f'unsupported folder {folder!r}')


def fold(seqs, version):
    """MFE structures for seqs with stock RNAfold -T 37, one process per call."""
    if not seqs:
        return []
    nops = '-noPS' if version == 1 else '--noPS'
    with tempfile.TemporaryDirectory() as tmp:  # 1.8.5 can still drop rna.ps in cwd
        out = subprocess.run([RNAFOLD[version], '-T', '37', nops], input='\n'.join(seqs) + '\n',
                             capture_output=True, text=True, cwd=tmp, check=True,
                             timeout=60 + len(seqs)).stdout
    structs = [m.group(1) for m in re.finditer(r'^([.()]+)\s+\(\s*-?[\d.]+\)\s*$', out, re.M)]
    if len(structs) != len(seqs):
        raise RuntimeError(f'RNAfold returned {len(structs)} structures for {len(seqs)} sequences')
    return structs


def model_flags(version, dangles):
    """Flags that put a ViennaRNA 2.x binary on the model this puzzle's engine actually uses."""
    if version == 2:
        return [f'-d{dangles}', '-T', '37']
    if not TURNER1999:
        sys.exit('--unique-mfe on a Vienna 1 puzzle needs TURNER1999_PAR set to '
                 'misc/rna_turner1999.par from a ViennaRNA 2.x source tree')
    return [f'-d{dangles}', '-T', '37', '-P', TURNER1999]


def reference_energy(seq, version):
    """The MFE energy the benchmark's own engine reports, which the enumeration must reach."""
    nops = '-noPS' if version == 1 else '--noPS'
    with tempfile.TemporaryDirectory() as tmp:
        out = subprocess.run([RNAFOLD[version], '-T', '37', nops], input=seq + '\n',
                             capture_output=True, text=True, cwd=tmp, check=True).stdout
    m = list(ENERGY.finditer(out))
    if not m:
        raise RuntimeError(f'could not parse RNAfold output for {seq[:20]}...')
    return float(m[-1].group(2))


def cooptimal_count(seq, version):
    """(n_tied, certified). n_tied is meaningful only when certified is True.

    RNAsubopt enumerates a window around the MFE and RNAeval re-scores every candidate on the
    puzzle's own model, so ties are counted exactly. Generation always uses dangles=2: RNAsubopt
    does not correctly support dangles=1 -- it warns "backtracking failed in repeat" and silently
    omits structures, sometimes the MFE itself. The window is widened until the re-scored minimum
    independently reaches reference_energy; that agreement is the completeness certificate. Never
    inject the known MFE into the candidate list, which would mask an incomplete enumeration.
    """
    e_ref = reference_energy(seq, version)
    score_flags = model_flags(version, 1 if version == 1 else 2)
    for w in WINDOWS:
        with tempfile.TemporaryDirectory() as tmp:
            gen = subprocess.run([RNASUBOPT, '-e', str(w)] + model_flags(version, 2),
                                 input=seq + '\n', capture_output=True, text=True,
                                 cwd=tmp, check=True).stdout
            cands = [l.split()[0] for l in gen.splitlines()[1:] if l.strip()]
            cands = [c for c in cands if len(c) == len(seq) and not set(c) - set('.()')]
            if not cands or len(cands) > MAX_CANDS:
                return len(cands), False
            out = subprocess.run([RNAEVAL] + score_flags,
                                 input=''.join(f'{seq}\n{c}\n' for c in cands),
                                 capture_output=True, text=True, cwd=tmp, check=True).stdout
        es = [(m.group(1), float(m.group(2))) for m in ENERGY.finditer(out)]
        if len(es) != len(cands):
            raise RuntimeError(f'{RNAEVAL} returned {len(es)} energies for {len(cands)} structures')
        emin = min(e for _, e in es)
        if abs(emin - e_ref) < TIE:
            return len({s for s, e in es if abs(e - emin) < TIE}), True
        if emin < e_ref - TIE:
            raise RuntimeError(f'{RNASUBOPT} found {emin} below the reference MFE {e_ref}')
    return -1, False


def pair_counts(seq, struct):
    stack, counts = [], {'GC': 0, 'AU': 0, 'GU': 0, 'ANY': 0}
    for j, c in enumerate(struct):
        if c == '(':
            stack.append(j)
        elif c == ')':
            i = stack.pop()
            key = {'CG': 'GC', 'AU': 'AU', 'GU': 'GU'}.get(''.join(sorted(seq[i] + seq[j])), 'other')
            counts[key] = counts.get(key, 0) + 1
            counts['ANY'] += 1
    return counts


def longest_run(seq, base):
    return max((len(m) for m in re.findall(f'{base}+', seq)), default=0)


def parse_constraints(constraints):
    """Strict parser shared with the solver wrappers: [(name, int value)]; raises on anything odd."""
    c = constraints.split(',') if constraints else []
    if len(c) % 2:
        raise ValueError(f'odd number of tokens in constraint string {constraints!r}')
    out = []
    for name, value in zip(c[0::2], c[1::2]):
        if not (name == 'SHAPE' or name in MAX_PAIRS or name in MIN_PAIRS
                or re.fullmatch(r'CONSECUTIVE_[ACGU]', name)):
            raise ValueError(f'unknown constraint {name!r}')
        if not re.fullmatch(r'\d+', value):
            raise ValueError(f'bad value {value!r} for {name}')
        if name == 'SHAPE' and value != '0':
            raise ValueError(f'multistate SHAPE,{value} not supported')
        out.append((name, int(value)))
    return out


def validate_puzzle(puzzle):
    """Raise ValueError unless the puzzle record is well formed."""
    t, locks, begin = puzzle['structure'], puzzle['locks'], puzzle['begin_sequence']
    pid = puzzle['puzzle_id']
    if not (len(t) == len(locks) == len(begin)):
        raise ValueError(f'{pid}: structure/locks/begin_sequence lengths {len(t)}/{len(locks)}/{len(begin)}')
    if set(t) - set('.()') or set(locks) - set('ox') or set(begin) - set('ACGU'):
        raise ValueError(f'{pid}: bad alphabet in structure, locks or begin_sequence')
    depth = 0
    for c in t:
        depth += {'(': 1, ')': -1}.get(c, 0)
        if depth < 0:
            break
    if depth:
        raise ValueError(f'{pid}: unbalanced target structure')
    engine(puzzle['folder'])
    parse_constraints(puzzle['constraints'])


def failures(puzzle, seq, mfe, target_only=False):
    """List of human-readable reasons the design fails; empty list means solved."""
    validate_puzzle(puzzle)
    target = puzzle['structure']
    if len(seq) != len(target) or set(seq) - set('ACGU'):
        return [f'bad sequence (length {len(seq)} vs {len(target)}, or non-ACGU)']
    bad = [] if mfe == target else ['SHAPE: MFE != target']
    if target_only:
        return bad
    for i, (lock, b0, b) in enumerate(zip(puzzle['locks'], puzzle['begin_sequence'], seq)):
        if lock == 'x' and b != b0:
            bad.append(f'lock: position {i + 1} is {b}, must be {b0}')
    counts = pair_counts(seq, mfe)
    for name, v in parse_constraints(puzzle['constraints']):
        if name == 'SHAPE':
            pass
        elif name in MAX_PAIRS:
            if counts[MAX_PAIRS[name]] > v:
                bad.append(f'{name}: {counts[MAX_PAIRS[name]]} > {v}')
        elif name in MIN_PAIRS:
            if counts[MIN_PAIRS[name]] < v:
                bad.append(f'{name}: {counts[MIN_PAIRS[name]]} < {v}')
        elif name.startswith('CONSECUTIVE_'):
            run = longest_run(seq, name[-1])
            if run >= v:
                bad.append(f'{name}: run of {run} >= {v}')
    return bad


def check(puzzles, designs, target_only=False, unique_mfe=False):
    """puzzles: {puzzle_id: row}; designs: list of dicts with puzzle_id, sequence.
    Returns designs with mfe, solved, reasons added (and n_mfe_tied with unique_mfe)."""
    for v in (1, 2):
        todo = [d for d in designs if engine(puzzles[d['puzzle_id']]['folder']) == v]
        seqs = [d['sequence'].strip() for d in todo]  # no case/T->U normalization: ACGU only
        ok = [i for i, s in enumerate(seqs) if s and not set(s) - set('ACGU')]
        mfes = dict(zip(ok, fold([seqs[i] for i in ok], v)))
        for i, (d, seq) in enumerate(zip(todo, seqs)):
            mfe = mfes.get(i, '')
            reasons = failures(puzzles[d['puzzle_id']], seq, mfe, target_only)
            if unique_mfe:
                # Only worth asking once the design already folds to the target.
                n, cert = cooptimal_count(seq, v) if not reasons else ('', True)
                d['n_mfe_tied'] = n if cert else ''
                if not cert:
                    reasons.append('UNIQUE_MFE: undetermined, enumeration never reached the MFE')
                elif n and n > 1:
                    reasons.append(f'UNIQUE_MFE: {n} structures tie at the MFE')
            d.update(mfe=mfe, solved=int(not reasons), reasons='; '.join(reasons))
    return designs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('puzzles')
    ap.add_argument('designs')
    ap.add_argument('--target-only', action='store_true')
    ap.add_argument('--unique-mfe', action='store_true',
                    help='also require the target to be the only structure at the MFE')
    ap.add_argument('--out')
    a = ap.parse_args()
    puzzles = {r['puzzle_id']: r for r in csv.DictReader(open(a.puzzles))}
    designs = list(csv.DictReader(open(a.designs), delimiter='\t'))
    check(puzzles, designs, a.target_only, a.unique_mfe)
    out = open(a.out, 'w', newline='') if a.out else sys.stdout
    w = csv.DictWriter(out, list(designs[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(designs)
    solved = {d['puzzle_id'] for d in designs if d['solved']}
    print(f'{sum(d["solved"] for d in designs)}/{len(designs)} designs pass; '
          f'{len(solved)}/{len({d["puzzle_id"] for d in designs})} puzzles solved', file=sys.stderr)


if __name__ == '__main__':
    main()
