#include "affine_align.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TOLERANCE 1e-9

static int failures = 0;

#define CHECK(condition)                                                    \
    do {                                                                    \
        if (!(condition)) {                                                 \
            fprintf(stderr,                                                 \
                    "FAIL %s:%d: %s\n",                                     \
                    __FILE__,                                               \
                    __LINE__,                                               \
                    #condition);                                            \
            ++failures;                                                     \
        }                                                                   \
    } while (0)

static char *remove_gaps(const char *alignment)
{
    const size_t length = strlen(alignment);
    char *sequence = malloc(length + 1);
    size_t output = 0;

    if (sequence == NULL) {
        return NULL;
    }
    for (size_t index = 0; index < length; ++index) {
        if (alignment[index] != '-') {
            sequence[output++] = alignment[index];
        }
    }
    sequence[output] = '\0';
    return sequence;
}

static void enumerate_alignment_scores(const char *sequence_a,
                                       const char *sequence_b,
                                       size_t index_a,
                                       size_t index_b,
                                       char previous_gap,
                                       double score,
                                       const aga_scoring_t *scoring,
                                       double *best_score)
{
    const size_t length_a = strlen(sequence_a);
    const size_t length_b = strlen(sequence_b);

    if (index_a == length_a && index_b == length_b) {
        if (score > *best_score) {
            *best_score = score;
        }
        return;
    }

    if (index_a < length_a && index_b < length_b) {
        const double substitution =
            sequence_a[index_a] == sequence_b[index_b]
                ? scoring->match_score
                : -scoring->mismatch_penalty;
        enumerate_alignment_scores(sequence_a,
                                   sequence_b,
                                   index_a + 1,
                                   index_b + 1,
                                   '\0',
                                   score + substitution,
                                   scoring,
                                   best_score);
    }

    if (index_a < length_a) {
        const double penalty =
            previous_gap == 'B'
                ? scoring->gap_extend_penalty
                : scoring->gap_open_penalty;
        enumerate_alignment_scores(sequence_a,
                                   sequence_b,
                                   index_a + 1,
                                   index_b,
                                   'B',
                                   score - penalty,
                                   scoring,
                                   best_score);
    }

    if (index_b < length_b) {
        const double penalty =
            previous_gap == 'A'
                ? scoring->gap_extend_penalty
                : scoring->gap_open_penalty;
        enumerate_alignment_scores(sequence_a,
                                   sequence_b,
                                   index_a,
                                   index_b + 1,
                                   'A',
                                   score - penalty,
                                   scoring,
                                   best_score);
    }
}

static double brute_force_score(const char *sequence_a,
                                const char *sequence_b,
                                const aga_scoring_t *scoring)
{
    double best_score = -INFINITY;
    enumerate_alignment_scores(sequence_a,
                               sequence_b,
                               0,
                               0,
                               '\0',
                               0.0,
                               scoring,
                               &best_score);
    return best_score;
}

static void binary_sequence(unsigned int bits,
                            size_t length,
                            char *sequence)
{
    for (size_t index = 0; index < length; ++index) {
        sequence[index] =
            (bits & (1U << index)) == 0U ? 'A' : 'B';
    }
    sequence[length] = '\0';
}

static void check_invariants(const char *sequence_a,
                             const char *sequence_b,
                             const aga_scoring_t *scoring,
                             const aga_alignment_t *alignment)
{
    char *recovered_a = remove_gaps(alignment->aligned_a);
    char *recovered_b = remove_gaps(alignment->aligned_b);

    CHECK(recovered_a != NULL);
    CHECK(recovered_b != NULL);
    if (recovered_a != NULL) {
        CHECK(strcmp(recovered_a, sequence_a) == 0);
    }
    if (recovered_b != NULL) {
        CHECK(strcmp(recovered_b, sequence_b) == 0);
    }
    CHECK(strlen(alignment->aligned_a) == alignment->length);
    CHECK(strlen(alignment->aligned_b) == alignment->length);

    const double rescored = aga_rescore_alignment(alignment, scoring);
    CHECK(isfinite(rescored));
    CHECK(fabs(rescored - alignment->score) < TOLERANCE);

    free(recovered_a);
    free(recovered_b);
}

static aga_alignment_t align_checked(const char *sequence_a,
                                     const char *sequence_b,
                                     const aga_scoring_t *scoring)
{
    aga_alignment_t alignment = {0};
    const aga_status_t status =
        aga_global_align(sequence_a, sequence_b, scoring, &alignment);

    CHECK(status == AGA_OK);
    if (status == AGA_OK) {
        check_invariants(sequence_a, sequence_b, scoring, &alignment);
    }
    return alignment;
}

static void test_empty_sequences(void)
{
    const aga_scoring_t scoring = {2.0, 1.0, 3.0, 1.0};
    aga_alignment_t alignment = align_checked("", "", &scoring);

    CHECK(alignment.length == 0);
    CHECK(fabs(alignment.score) < TOLERANCE);
    CHECK(strcmp(alignment.aligned_a, "") == 0);
    CHECK(strcmp(alignment.aligned_b, "") == 0);
    aga_alignment_free(&alignment);
}

static void test_empty_against_nonempty(void)
{
    const aga_scoring_t scoring = {2.0, 1.0, 3.0, 0.5};
    aga_alignment_t alignment = align_checked("", "ACGT", &scoring);

    CHECK(fabs(alignment.score - (-4.5)) < TOLERANCE);
    CHECK(strcmp(alignment.aligned_a, "----") == 0);
    CHECK(strcmp(alignment.aligned_b, "ACGT") == 0);
    aga_alignment_free(&alignment);
}

static void test_identical_sequences(void)
{
    const aga_scoring_t scoring = {2.0, 1.0, 3.0, 1.0};
    aga_alignment_t alignment = align_checked("ACGT", "ACGT", &scoring);
    const aga_metrics_t metrics = aga_alignment_metrics(&alignment);

    CHECK(fabs(alignment.score - 8.0) < TOLERANCE);
    CHECK(strcmp(alignment.aligned_a, "ACGT") == 0);
    CHECK(strcmp(alignment.aligned_b, "ACGT") == 0);
    CHECK(metrics.identities == 4);
    CHECK(metrics.mismatches == 0);
    CHECK(metrics.gap_characters == 0);
    CHECK(metrics.gap_openings == 0);
    aga_alignment_free(&alignment);
}

static void test_internal_gap(void)
{
    const aga_scoring_t scoring = {2.0, 1.0, 3.0, 1.0};
    aga_alignment_t alignment = align_checked("ACGT", "AGT", &scoring);
    const aga_metrics_t metrics = aga_alignment_metrics(&alignment);

    CHECK(fabs(alignment.score - 3.0) < TOLERANCE);
    CHECK(strcmp(alignment.aligned_a, "ACGT") == 0);
    CHECK(strcmp(alignment.aligned_b, "A-GT") == 0);
    CHECK(metrics.identities == 3);
    CHECK(metrics.gap_characters == 1);
    CHECK(metrics.gap_openings == 1);
    aga_alignment_free(&alignment);
}

static void test_gap_switch_when_cheaper_than_mismatch(void)
{
    const aga_scoring_t scoring = {2.0, 10.0, 1.0, 0.5};
    aga_alignment_t alignment = align_checked("A", "B", &scoring);
    const aga_metrics_t metrics = aga_alignment_metrics(&alignment);

    CHECK(fabs(alignment.score - (-2.0)) < TOLERANCE);
    CHECK(alignment.length == 2);
    CHECK(metrics.gap_characters == 2);
    CHECK(metrics.gap_openings == 2);
    aga_alignment_free(&alignment);
}

static void test_fractional_penalties(void)
{
    const aga_scoring_t scoring = {5.0, 4.0, 10.0, 0.5};
    aga_alignment_t alignment = align_checked(
        "ATGAGTCTCTCTGATAAGGACAAGGCTGCTGTGAAAGCCCTATGG",
        "CTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAG",
        &scoring);
    const aga_metrics_t metrics = aga_alignment_metrics(&alignment);

    CHECK(fabs(alignment.score - 41.5) < TOLERANCE);
    CHECK(alignment.length == 57);
    CHECK(metrics.identities == 28);
    CHECK(metrics.mismatches == 5);
    CHECK(metrics.gap_characters == 24);
    CHECK(metrics.gap_openings == 7);
    CHECK(fabs(alignment.score -
               aga_rescore_alignment(&alignment, &scoring)) < TOLERANCE);
    aga_alignment_free(&alignment);
}

static void test_symmetry(void)
{
    const aga_scoring_t scoring = {3.0, 2.0, 4.0, 0.75};
    aga_alignment_t forward =
        align_checked("GATTACA", "GCATGCU", &scoring);
    aga_alignment_t reverse =
        align_checked("GCATGCU", "GATTACA", &scoring);

    CHECK(fabs(forward.score - reverse.score) < TOLERANCE);
    aga_alignment_free(&forward);
    aga_alignment_free(&reverse);
}

static void test_moderately_long_input(void)
{
    const size_t length = 512;
    char *sequence_a = malloc(length + 1);
    char *sequence_b = malloc(length);
    const aga_scoring_t scoring = {1.0, 1.0, 2.0, 0.25};

    CHECK(sequence_a != NULL);
    CHECK(sequence_b != NULL);
    if (sequence_a == NULL || sequence_b == NULL) {
        free(sequence_a);
        free(sequence_b);
        return;
    }

    memset(sequence_a, 'A', length);
    sequence_a[length] = '\0';
    memset(sequence_b, 'A', length - 1);
    sequence_b[length - 1] = '\0';

    aga_alignment_t alignment =
        align_checked(sequence_a, sequence_b, &scoring);
    CHECK(alignment.length == length);
    CHECK(fabs(alignment.score - 509.0) < TOLERANCE);

    aga_alignment_free(&alignment);
    free(sequence_a);
    free(sequence_b);
}

static void test_invalid_scoring(void)
{
    const aga_scoring_t invalid = {2.0, -1.0, 3.0, 1.0};
    aga_alignment_t alignment = {0};

    CHECK(aga_global_align("A", "A", &invalid, &alignment) ==
          AGA_INVALID_ARGUMENT);
}

static void test_against_exhaustive_oracle(void)
{
    const aga_scoring_t scoring_schemes[] = {
        {2.0, 1.0, 3.0, 1.0},
        {1.5, 4.0, 0.75, 0.25}};
    char sequence_a[4];
    char sequence_b[4];

    for (size_t scheme = 0;
         scheme < sizeof(scoring_schemes) / sizeof(scoring_schemes[0]);
         ++scheme) {
        for (size_t length_a = 0; length_a <= 3; ++length_a) {
            for (size_t length_b = 0; length_b <= 3; ++length_b) {
                const unsigned int count_a = 1U << length_a;
                const unsigned int count_b = 1U << length_b;

                for (unsigned int bits_a = 0; bits_a < count_a; ++bits_a) {
                    binary_sequence(bits_a, length_a, sequence_a);
                    for (unsigned int bits_b = 0;
                         bits_b < count_b;
                         ++bits_b) {
                        binary_sequence(bits_b, length_b, sequence_b);
                        const double expected =
                            brute_force_score(sequence_a,
                                              sequence_b,
                                              &scoring_schemes[scheme]);
                        aga_alignment_t alignment =
                            align_checked(sequence_a,
                                          sequence_b,
                                          &scoring_schemes[scheme]);

                        CHECK(fabs(alignment.score - expected) < TOLERANCE);
                        aga_alignment_free(&alignment);
                    }
                }
            }
        }
    }
}

int main(void)
{
    test_empty_sequences();
    test_empty_against_nonempty();
    test_identical_sequences();
    test_internal_gap();
    test_gap_switch_when_cheaper_than_mismatch();
    test_fractional_penalties();
    test_symmetry();
    test_moderately_long_input();
    test_invalid_scoring();
    test_against_exhaustive_oracle();

    if (failures != 0) {
        fprintf(stderr, "%d test assertion(s) failed\n", failures);
        return EXIT_FAILURE;
    }

    printf("All affine alignment tests passed.\n");
    return EXIT_SUCCESS;
}
