# affine-gap-aligner

[![CI](https://github.com/Xianxin-Long/affine-gap-aligner/actions/workflows/ci.yml/badge.svg)](https://github.com/Xianxin-Long/affine-gap-aligner/actions/workflows/ci.yml)
[![C11](https://img.shields.io/badge/C-11-blue.svg)](https://en.cppreference.com/w/c/11)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An exact pairwise global sequence aligner written in C11. It implements a
Gotoh-style three-state dynamic program for affine gap penalties, reconstructs
the optimal alignment, and reports alignment statistics through human-readable
or JSON output.

This project began as a dynamic-programming course implementation and was
rebuilt as a tested command-line tool and reusable C library.

## Highlights

- Exact global alignment in \(O(mn)\) time
- Affine gap model \(g(k)=\sigma+(k-1)\varepsilon\)
- Packed one-byte-per-cell traceback with \(O(n)\) score storage
- Direct sequence input or two-record FASTA input
- Deterministic traceback without floating-point equality checks
- Text and machine-readable JSON output
- CTest suite, strict compiler warnings, AddressSanitizer/UBSan, and
  Linux/macOS continuous integration

## Build

Requirements: a C11 compiler and CMake 3.16 or newer.

```bash
git clone https://github.com/Xianxin-Long/affine-gap-aligner.git
cd affine-gap-aligner
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

The executable is written to `build/affine-align`.

A Makefile is also provided:

```bash
make
make test
```

## Quick start

Align two sequences supplied on the command line:

```bash
./build/affine-align \
  --seq-a GATTACA \
  --seq-b GCATGCU \
  --match 2 \
  --mismatch 2 \
  --gap-open 2 \
  --gap-extend 1
```

```text
Affine Global Alignment
=======================
Score          : 0
Aligned length : 8
Identity       : 4/8 (50.0%)
Mismatches     : 2
Gap characters : 2
Gap openings   : 2

SeqA      1  G-ATTACA  7
             | | |.|.
SeqB      1  GCA-TGCU  7
```

Or read exactly two records from a FASTA file:

```bash
./build/affine-align --fasta examples/pair.fasta
```

Use `--format json` for structured output:

```bash
./build/affine-align --fasta examples/pair.fasta --format json
```

## Scoring convention

Scores are maximized. Penalties are passed as non-negative magnitudes and are
subtracted internally:

- match: `+match`
- mismatch: `-mismatch`
- gap of length \(k\): `-(gap_open + (k - 1) * gap_extend)`

Thus a one-character gap pays the gap-open penalty exactly once. See
[the algorithm note](docs/algorithm.md) for the recurrences, initialization,
traceback representation, complexity analysis, and correctness sketch.

## Command-line interface

```text
Usage:
  affine-align --seq-a SEQUENCE --seq-b SEQUENCE [options]
  affine-align --fasta FILE [options]

Input:
  -a, --seq-a SEQUENCE       First sequence
  -b, --seq-b SEQUENCE       Second sequence
  -f, --fasta FILE           FASTA file containing exactly two records

Scoring:
  -m, --match VALUE          Match score (default: 2)
  -x, -M, --mismatch VALUE   Mismatch penalty (default: 1)
  -o, --gap-open VALUE       Gap-open penalty (default: 3)
  -e, --gap-extend VALUE     Gap-extension penalty (default: 1)

Output:
      --format text|json     Output format (default: text)
  -w, --width INTEGER        Alignment columns per block (default: 60)
```

Input symbols are normalized to uppercase. The CLI accepts alphabetic IUPAC
symbols, so it can handle nucleotide and protein sequences, including
ambiguous residue codes. The current scoring model uses a uniform match score
and mismatch penalty rather than a substitution matrix.

## Repository layout

```text
include/affine_align.h   Public library API
src/affine_align.c      Three-state DP and packed traceback
src/main.c              CLI, FASTA parser, and output formatting
tests/                  Unit and invariant tests
docs/algorithm.md       Mathematical formulation and correctness sketch
examples/pair.fasta     Reproducible example input
```

The tests include an independent exhaustive-search oracle. For every pair of
binary-alphabet sequences up to length three under two scoring schemes, the
dynamic-programming score must equal the best score over all possible
alignment paths.

## License

[MIT](LICENSE)
