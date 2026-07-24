#ifndef AFFINE_ALIGN_H
#define AFFINE_ALIGN_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double match_score;
    double mismatch_penalty;
    double gap_open_penalty;
    double gap_extend_penalty;
} aga_scoring_t;

typedef struct {
    char *aligned_a;
    char *aligned_b;
    size_t length;
    double score;
} aga_alignment_t;

typedef struct {
    size_t identities;
    size_t mismatches;
    size_t gap_characters;
    size_t gap_openings;
    double identity_fraction;
} aga_metrics_t;

typedef enum {
    AGA_OK = 0,
    AGA_INVALID_ARGUMENT,
    AGA_SIZE_OVERFLOW,
    AGA_OUT_OF_MEMORY,
    AGA_INTERNAL_ERROR
} aga_status_t;

/*
 * Compute an exact global alignment with affine gap penalties.
 *
 * A gap of length k costs:
 *
 *     gap_open_penalty + (k - 1) * gap_extend_penalty
 *
 * `sequence_a` and `sequence_b` must be NUL-terminated strings. The library
 * treats their bytes as case-sensitive symbols; input normalization belongs to
 * the caller. On success, `out` owns two heap-allocated aligned strings that
 * must be released with aga_alignment_free().
 */
aga_status_t aga_global_align(const char *sequence_a,
                              const char *sequence_b,
                              const aga_scoring_t *scoring,
                              aga_alignment_t *out);

void aga_alignment_free(aga_alignment_t *alignment);

aga_metrics_t aga_alignment_metrics(const aga_alignment_t *alignment);

/*
 * Recalculate the score of an alignment. Returns NaN when the alignment is
 * malformed (different string lengths or a column containing two gaps).
 */
double aga_rescore_alignment(const aga_alignment_t *alignment,
                             const aga_scoring_t *scoring);

const char *aga_status_string(aga_status_t status);

#ifdef __cplusplus
}
#endif

#endif
