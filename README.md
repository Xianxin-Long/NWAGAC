# Affine Gap Aligner

A small command-line program written in C for global pairwise sequence alignment using the Needleman-Wunsch algorithm with an affine gap penalty.

The program aligns two nucleotide or protein sequences, reports the optimal alignment score, alignment length, number of identities, and number of gap positions, and prints the aligned sequences.

## Algorithm

The implementation uses three dynamic-programming matrices to distinguish matches/mismatches from gaps in either sequence.

For a gap of length $k$, the penalty is

$$
g(k) = o + (k - 1)e.
$$

where:

- $o$ is the gap-opening penalty;
- $e$ is the gap-extension penalty.

The algorithm performs global alignment, so both sequences are aligned from beginning to end.

For input sequences of lengths $m$ and $n$, each dynamic-programming matrix has dimensions $(m+1) \times (n+1)$. The algorithm therefore runs in $O(mn)$ time and uses $O(mn)$ memory.

## Requirements

- A C11-compatible compiler such as GCC or Clang
- A POSIX-like environment providing `getopt`

## Build

Compile the program with GCC:

```bash
gcc -std=c11 -O2 -Wall -Wextra -Wpedantic needle.c -o needle
```

## Usage

```bash
./needle -a <seq1> -b <seq2> -m <match_score> -M <mismatch_penalty> -o <gap_open> -e <gap_extend>
```

All six arguments are required.

| Option | Description |
| --- | --- |
| `-a` | First sequence |
| `-b` | Second sequence |
| `-m` | Match score |
| `-M` | Mismatch penalty |
| `-o` | Gap-opening penalty |
| `-e` | Gap-extension penalty |
| `-h` | Show help |

Sequence input is case-insensitive and is converted to uppercase before alignment.

## Example

```bash
./needle -a ATGC -b ATCG -m 2 -M 1 -o 1 -e 0.5
```

Example output:

```text
########################################
# Program: needle
# ...
########################################
#=======================================
# Length: 4
# Identity: 2
# Gaps: 0
# Score: 2.0
#=======================================
Alignments:
SeqA  1  ATGC  4
         ||  
SeqB  1  ATCG  4
```

The exact date line is generated when the program runs.

## Slides

The course presentation for this project is available here:

[Pairwise Global Alignment Using Affine Gap Penalties](slides.pdf)

## Notes

The current implementation allocates the three dynamic-programming matrices on the stack. Although the theoretical memory complexity is $O(mn)$, very long input sequences can exceed the available stack size, so the program is intended for short to moderate sequence lengths.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
