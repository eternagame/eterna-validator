# Eterna Validator

Check a solution to an RNA sequence against an Eterna puzzle offline.

Used by [Eterna Hard](https://github.com/eternagame/eterna-hard) and
[Eterna100](https://github.com/eternagame/eterna100-benchmarking).

## What counts as solved

A design passes when, as in the Eterna game (EternaJS `src/eterna/constraints`):

- the folding engine's MFE structure equals the target structure;
- every locked base (`x` in `locks`) still matches `begin_sequence`;
- every pair-type limit holds: `GC`/`AUMAX`/`GUMAX` are maxima, `GCMIN`/`AU`/`GU`/`PAIRS` are minima;
- `CONSECUTIVE_X` holds: the longest run of base `X` is shorter than the stated value.

Pairs are counted on the MFE structure. An unrecognized constraint is an error.

## Install ViennaRNA

Eterna puzzles are folded by one of two engines, named in each puzzle's `folder` column:
`Vienna` (or `Vienna1`) means **ViennaRNA 1.8.5**, `Vienna2` means **ViennaRNA 2.x**.
Install whichever your puzzles need.

### ViennaRNA 1.8.5

```sh
curl -fsSLO https://www.tbi.univie.ac.at/RNA/packages/source/ViennaRNA-1.8.5.tar.gz
tar xzf ViennaRNA-1.8.5.tar.gz && cd ViennaRNA-1.8.5
./configure --without-perl --without-kinfold --without-forester \
  CFLAGS="-O2 -std=gnu89 -fcommon -Wno-implicit-int -Wno-implicit-function-declaration"
make
```

The binary lands at `Progs/RNAfold`. The `CFLAGS` are needed because this 2009 release predates
several C standard changes; without them it will not build on a current compiler.

SHA-256 of the tarball:
`f4e2d94beaf77165e8321758e4ab0ad1c5d49879cefa12e48b07d09ed2d0ecf9`

### ViennaRNA 2.x

Package versions starting at ViennaRNA 2.1.9 (matching eternagame.org) up to 2.6.4 have been checked.

```sh
conda install -c bioconda viennarna     # or build from source as above
```

### Point the validator at them

```sh
export RNAFOLD_V1=/path/to/ViennaRNA-1.8.5/Progs/RNAfold
export RNAFOLD_V2=/path/to/RNAfold      # only needed for Vienna2 puzzles
```

`RNAFOLD_V2` defaults to `RNAfold` on your `PATH`.

## Pull puzzle definitions

`fetch_puzzles.py` reads Eterna's public API. No login is needed.

```sh
python3 fetch_puzzles.py 3058243 3308813 --out puzzles.csv   # by puzzle id
python3 fetch_puzzles.py --lab 16022932 --out puzzles.csv    # every puzzle in a lab
```

It writes `puzzle_id, title, folder, structure, constraints, locks, begin_sequence`, which are the columns
the validator needs. The CSVs shipped with [Eterna Hard](https://github.com/eternagame/eterna-hard) already carry these columns, so you can
pass one of those instead.

A puzzle's `folder` (i.e., its folding engine) is null in the API when it uses the default engine; the script records that as `Vienna`.

## Test a design

```sh
python3 check_design.py puzzles.csv designs.tsv
```

`designs.tsv` is tab-separated with columns `puzzle_id` and `sequence`; any other columns are
passed through. Results go to stdout (or `--out results.tsv`) with `mfe`, `solved` and `reasons`
added, and a summary line on stderr:

```
5/5 designs pass; 5/5 puzzles solved
```

`--target-only` checks the structure alone and ignores the player constraints, which is the
original Eterna100 benchmark condition.

## Worked example

Three Eterna100-V1 puzzles, scored against the sample solutions published in
[eterna100-benchmarking](https://github.com/eternagame/eterna100-benchmarking). They cover the
three cases: no player constraints, a pair-type minimum, and locked bases.

The files this example writes are checked in under `example/`, so you can compare your own run
against them.

```sh
export RNAFOLD_V1=/path/to/ViennaRNA-1.8.5/Progs/RNAfold
curl -fsSLO https://raw.githubusercontent.com/eternagame/eterna100-benchmarking/master/data/eterna100_puzzles.tsv
python3 fetch_puzzles.py 20111 215430 446042 --out puzzles.csv
```

```
  1/3 20111
  2/3 215430
  3/3 446042
wrote puzzles.csv (3 puzzles)
```

`puzzles.csv`, with the long sequence fields cut short here:

```
puzzle_id,title,folder,structure,constraints,locks,begin_sequence
20111,Simple Hairpin,Vienna,(((((......))))),"SHAPE,0",oooooooooooooooo,AAAAAAAAAAAAAAAA
215430,Arabidopsis Thaliana 6 RNA...,Vienna,(((((((....(((...,"SHAPE,0,GU,2",oooooooo...,AAAAAAAA...
446042,Six legd turtle 2,Vienna,.....((...((....,"SHAPE,0,CONSECUTIVE_G,4,CONSECUTIVE_C,4",xxoooo...,GGAAAA...
```

So puzzle 215430 needs at least 2 GU pairs, and 446042 allows runs of at most 3 G and 3 C and
locks 22 of its 98 bases. Take those three sample solutions (column 9 holds the V1/Vienna1 ones):

```sh
awk -F'\t' 'BEGIN {OFS="\t"; print "puzzle_id", "sequence"}
            $5 == "20111" || $5 == "215430" || $5 == "446042" {print $5, $9}' \
  eterna100_puzzles.tsv > designs.tsv
python3 check_design.py puzzles.csv designs.tsv --out results.tsv
```

```
3/3 designs pass; 3/3 puzzles solved
```

`results.tsv` adds `mfe`, `solved` and `reasons` to each input row (the `sequence` column is
dropped here for width; `reasons` is empty on a pass):

```
puzzle_id  mfe                                                   solved  reasons
20111      (((((......)))))                                      1
215430     (((((((....(((...........)))((((((((..(((((((((((...  1
446042     .....((...((....))...((....))...((....))...((....))...  1
```

All three pass with the constraints enforced, so `--target-only` gives the same 3/3 here.

### A design that fails

The two conditions do come apart. Eterna100-V2 puzzle 6158733 ("Mat - Martian 2", Vienna2) asks
for at least 12 GU pairs and at most 4 GC pairs. The design below folds to the target but uses
4 GU and 43 GC, so it is not a legal in-game solution. (It is the solution
[eterna100-benchmarking](https://github.com/eternagame/eterna100-benchmarking) published for this
puzzle before its sample solutions were replaced with constraint-satisfying ones; the solution
shipped there now satisfies both conditions.)

```sh
export RNAFOLD_V2=/path/to/RNAfold          # ViennaRNA 2.x, for a Vienna2 puzzle
python3 fetch_puzzles.py 6158733 --out martian.csv
printf 'puzzle_id\tsequence\n6158733\tGUGGGAAACCACAAAGGUAUCGAAGCACGUGGAAAAAACACGUGCCGAUACCAAAGGCAUCCCGAAAAAGAUCUAGAGGGGGAGGAAAAAACUCCUUCCUAGAUCAAAGCUAAGCGUGUAGCGAAAAAAGCUAUACAGCUUAGCAAAAAAGGGAUGUCAAAGUUGUGCCGUGGCGGAAAAAACGCCACGAAGCACAACAAACAACGAAAGUUG\n' > martian.tsv
python3 check_design.py martian.csv martian.tsv --out martian_results.tsv
```

```
0/1 designs pass; 0/1 puzzles solved
```

The `reasons` column says which constraints broke and by how much, `<` and `>` pointing at the
offending side of each limit:

```
GU: 4 < 12; GC: 43 > 4
```

Re-run that with `--target-only` and it passes `1/1`, because the MFE structure does equal the
target. Across the whole benchmark the two conditions now agree: every disclosed sample solution
published in [eterna100-benchmarking](https://github.com/eternagame/eterna100-benchmarking) passes
both ways — 99/99 for Eterna100-V1 and 99/99 for Eterna100-V2 under `--target-only`, and the same
with their player constraints enforced (the 100th puzzle is undisclosed).

## Scope

Single-state puzzles only. 

More advanced folding engines (EternaFold) are scored differently; validator will be made available at a later date.

Switch (multi-state) puzzles are also scored differently. The game applies a ligand binding bonus during folding, which needs a patched engine or a constrained-folding workaround, which will be made available at a later date.

## Requirements

Python 3, standard library only, plus the ViennaRNA binaries above.
