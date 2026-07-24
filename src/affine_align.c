#include "affine_align.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef enum {
    STATE_MATCH = 0,
    STATE_GAP_B = 1,
    STATE_GAP_A = 2,
    STATE_NONE = 3
} state_t;

enum {
    TRACE_MATCH_SHIFT = 0,
    TRACE_GAP_B_SHIFT = 2,
    TRACE_GAP_A_SHIFT = 4,
    TRACE_MASK = 0x03
};

typedef struct {
    double score;
    state_t state;
} candidate_t;

static candidate_t best_of_three(candidate_t first,
                                 candidate_t second,
                                 candidate_t third)
{
    candidate_t best = first;

    if (second.score > best.score) {
        best = second;
    }
    if (third.score > best.score) {
        best = third;
    }
    return best;
}

static void set_predecessor(uint8_t *trace,
                            size_t cell,
                            unsigned int shift,
                            state_t predecessor)
{
    const uint8_t mask = (uint8_t)(TRACE_MASK << shift);
    trace[cell] =
        (uint8_t)((trace[cell] & (uint8_t)~mask) |
                  ((uint8_t)predecessor << shift));
}

static state_t get_predecessor(const uint8_t *trace,
                               size_t cell,
                               unsigned int shift)
{
    return (state_t)((trace[cell] >> shift) & TRACE_MASK);
}

static void reverse_string(char *text, size_t length)
{
    for (size_t left = 0, right = length == 0 ? 0 : length - 1;
         left < right;
         ++left, --right) {
        const char temporary = text[left];
        text[left] = text[right];
        text[right] = temporary;
    }
}

static int scoring_is_valid(const aga_scoring_t *scoring)
{
    return scoring != NULL &&
           isfinite(scoring->match_score) &&
           isfinite(scoring->mismatch_penalty) &&
           isfinite(scoring->gap_open_penalty) &&
           isfinite(scoring->gap_extend_penalty) &&
           scoring->mismatch_penalty >= 0.0 &&
           scoring->gap_open_penalty >= 0.0 &&
           scoring->gap_extend_penalty >= 0.0;
}

aga_status_t aga_global_align(const char *sequence_a,
                              const char *sequence_b,
                              const aga_scoring_t *scoring,
                              aga_alignment_t *out)
{
    if (sequence_a == NULL || sequence_b == NULL || out == NULL ||
        !scoring_is_valid(scoring)) {
        return AGA_INVALID_ARGUMENT;
    }

    out->aligned_a = NULL;
    out->aligned_b = NULL;
    out->length = 0;
    out->score = 0.0;

    const size_t length_a = strlen(sequence_a);
    const size_t length_b = strlen(sequence_b);

    if (length_a == SIZE_MAX || length_b == SIZE_MAX) {
        return AGA_SIZE_OVERFLOW;
    }

    const size_t rows = length_a + 1;
    const size_t columns = length_b + 1;

    if (rows > SIZE_MAX / columns) {
        return AGA_SIZE_OVERFLOW;
    }
    const size_t cells = rows * columns;

    if (columns > SIZE_MAX / (6 * sizeof(double)) ||
        length_a > SIZE_MAX - length_b ||
        length_a + length_b == SIZE_MAX) {
        return AGA_SIZE_OVERFLOW;
    }

    uint8_t *trace = malloc(cells * sizeof(*trace));
    double *score_rows = malloc(6 * columns * sizeof(*score_rows));
    const size_t maximum_alignment_length = length_a + length_b;
    char *aligned_a = malloc(maximum_alignment_length + 1);
    char *aligned_b = malloc(maximum_alignment_length + 1);

    if (trace == NULL || score_rows == NULL ||
        aligned_a == NULL || aligned_b == NULL) {
        free(trace);
        free(score_rows);
        free(aligned_a);
        free(aligned_b);
        return AGA_OUT_OF_MEMORY;
    }

    memset(trace, 0xff, cells * sizeof(*trace));

    double *match_previous = score_rows;
    double *gap_b_previous = score_rows + columns;
    double *gap_a_previous = score_rows + 2 * columns;
    double *match_current = score_rows + 3 * columns;
    double *gap_b_current = score_rows + 4 * columns;
    double *gap_a_current = score_rows + 5 * columns;

    match_previous[0] = 0.0;
    gap_b_previous[0] = -INFINITY;
    gap_a_previous[0] = -INFINITY;

    for (size_t column = 1; column <= length_b; ++column) {
        match_previous[column] = -INFINITY;
        gap_b_previous[column] = -INFINITY;
        gap_a_previous[column] =
            column == 1
                ? -scoring->gap_open_penalty
                : gap_a_previous[column - 1] -
                      scoring->gap_extend_penalty;

        set_predecessor(trace,
                        column,
                        TRACE_GAP_A_SHIFT,
                        column == 1 ? STATE_MATCH : STATE_GAP_A);
    }

    for (size_t row = 1; row <= length_a; ++row) {
        const size_t first_cell = row * columns;
        match_current[0] = -INFINITY;
        gap_a_current[0] = -INFINITY;
        gap_b_current[0] =
            row == 1
                ? -scoring->gap_open_penalty
                : gap_b_previous[0] -
                      scoring->gap_extend_penalty;

        set_predecessor(trace,
                        first_cell,
                        TRACE_GAP_B_SHIFT,
                        row == 1 ? STATE_MATCH : STATE_GAP_B);

        for (size_t column = 1; column <= length_b; ++column) {
            const size_t cell = first_cell + column;
            const double substitution =
                sequence_a[row - 1] == sequence_b[column - 1]
                    ? scoring->match_score
                    : -scoring->mismatch_penalty;

            /*
             * Tie-breaking is deterministic:
             *   M: diagonal M, X, Y
             *   X/Y: extend the current gap, then open from M, then switch
             * Final state: M, X, Y
             */
            const candidate_t match_best = best_of_three(
                (candidate_t){match_previous[column - 1], STATE_MATCH},
                (candidate_t){gap_b_previous[column - 1], STATE_GAP_B},
                (candidate_t){gap_a_previous[column - 1], STATE_GAP_A});
            match_current[column] = match_best.score + substitution;
            set_predecessor(trace,
                            cell,
                            TRACE_MATCH_SHIFT,
                            match_best.state);

            const candidate_t gap_b_best = best_of_three(
                (candidate_t){
                    gap_b_previous[column] -
                        scoring->gap_extend_penalty,
                    STATE_GAP_B},
                (candidate_t){
                    match_previous[column] -
                        scoring->gap_open_penalty,
                    STATE_MATCH},
                (candidate_t){
                    gap_a_previous[column] -
                        scoring->gap_open_penalty,
                    STATE_GAP_A});
            gap_b_current[column] = gap_b_best.score;
            set_predecessor(trace,
                            cell,
                            TRACE_GAP_B_SHIFT,
                            gap_b_best.state);

            const candidate_t gap_a_best = best_of_three(
                (candidate_t){
                    gap_a_current[column - 1] -
                        scoring->gap_extend_penalty,
                    STATE_GAP_A},
                (candidate_t){
                    match_current[column - 1] -
                        scoring->gap_open_penalty,
                    STATE_MATCH},
                (candidate_t){
                    gap_b_current[column - 1] -
                        scoring->gap_open_penalty,
                    STATE_GAP_B});
            gap_a_current[column] = gap_a_best.score;
            set_predecessor(trace,
                            cell,
                            TRACE_GAP_A_SHIFT,
                            gap_a_best.state);
        }

        double *temporary = match_previous;
        match_previous = match_current;
        match_current = temporary;

        temporary = gap_b_previous;
        gap_b_previous = gap_b_current;
        gap_b_current = temporary;

        temporary = gap_a_previous;
        gap_a_previous = gap_a_current;
        gap_a_current = temporary;
    }

    const candidate_t final = best_of_three(
        (candidate_t){match_previous[length_b], STATE_MATCH},
        (candidate_t){gap_b_previous[length_b], STATE_GAP_B},
        (candidate_t){gap_a_previous[length_b], STATE_GAP_A});

    size_t row = length_a;
    size_t column = length_b;
    size_t alignment_length = 0;
    state_t state = final.state;
    aga_status_t status = AGA_OK;

    while (row > 0 || column > 0) {
        const size_t cell = row * columns + column;
        state_t predecessor = STATE_NONE;

        if (state == STATE_MATCH && row > 0 && column > 0) {
            predecessor =
                get_predecessor(trace, cell, TRACE_MATCH_SHIFT);
            aligned_a[alignment_length] = sequence_a[row - 1];
            aligned_b[alignment_length] = sequence_b[column - 1];
            --row;
            --column;
        } else if (state == STATE_GAP_B && row > 0) {
            predecessor =
                get_predecessor(trace, cell, TRACE_GAP_B_SHIFT);
            aligned_a[alignment_length] = sequence_a[row - 1];
            aligned_b[alignment_length] = '-';
            --row;
        } else if (state == STATE_GAP_A && column > 0) {
            predecessor =
                get_predecessor(trace, cell, TRACE_GAP_A_SHIFT);
            aligned_a[alignment_length] = '-';
            aligned_b[alignment_length] = sequence_b[column - 1];
            --column;
        } else {
            status = AGA_INTERNAL_ERROR;
            break;
        }

        ++alignment_length;
        state = predecessor;
    }

    if (status == AGA_OK) {
        reverse_string(aligned_a, alignment_length);
        reverse_string(aligned_b, alignment_length);
        aligned_a[alignment_length] = '\0';
        aligned_b[alignment_length] = '\0';

        out->aligned_a = aligned_a;
        out->aligned_b = aligned_b;
        out->length = alignment_length;
        out->score = final.score;
    } else {
        free(aligned_a);
        free(aligned_b);
    }

    free(trace);
    free(score_rows);
    return status;
}

void aga_alignment_free(aga_alignment_t *alignment)
{
    if (alignment == NULL) {
        return;
    }

    free(alignment->aligned_a);
    free(alignment->aligned_b);
    alignment->aligned_a = NULL;
    alignment->aligned_b = NULL;
    alignment->length = 0;
    alignment->score = 0.0;
}

aga_metrics_t aga_alignment_metrics(const aga_alignment_t *alignment)
{
    aga_metrics_t metrics = {0, 0, 0, 0, 0.0};

    if (alignment == NULL || alignment->aligned_a == NULL ||
        alignment->aligned_b == NULL) {
        return metrics;
    }

    char previous_gap = '\0';
    for (size_t index = 0; index < alignment->length; ++index) {
        const char symbol_a = alignment->aligned_a[index];
        const char symbol_b = alignment->aligned_b[index];
        char current_gap = '\0';

        if (symbol_a == '-') {
            current_gap = 'A';
            ++metrics.gap_characters;
        } else if (symbol_b == '-') {
            current_gap = 'B';
            ++metrics.gap_characters;
        } else if (symbol_a == symbol_b) {
            ++metrics.identities;
        } else {
            ++metrics.mismatches;
        }

        if (current_gap != '\0' && current_gap != previous_gap) {
            ++metrics.gap_openings;
        }
        previous_gap = current_gap;
    }

    if (alignment->length > 0) {
        metrics.identity_fraction =
            (double)metrics.identities / (double)alignment->length;
    }
    return metrics;
}

double aga_rescore_alignment(const aga_alignment_t *alignment,
                             const aga_scoring_t *scoring)
{
    if (alignment == NULL || alignment->aligned_a == NULL ||
        alignment->aligned_b == NULL || !scoring_is_valid(scoring) ||
        strlen(alignment->aligned_a) != alignment->length ||
        strlen(alignment->aligned_b) != alignment->length) {
        return NAN;
    }

    double score = 0.0;
    char previous_gap = '\0';

    for (size_t index = 0; index < alignment->length; ++index) {
        const char symbol_a = alignment->aligned_a[index];
        const char symbol_b = alignment->aligned_b[index];
        char current_gap = '\0';

        if (symbol_a == '-' && symbol_b == '-') {
            return NAN;
        }

        if (symbol_a == '-') {
            current_gap = 'A';
        } else if (symbol_b == '-') {
            current_gap = 'B';
        }

        if (current_gap != '\0') {
            score -= current_gap == previous_gap
                         ? scoring->gap_extend_penalty
                         : scoring->gap_open_penalty;
        } else {
            score += symbol_a == symbol_b
                         ? scoring->match_score
                         : -scoring->mismatch_penalty;
        }
        previous_gap = current_gap;
    }

    return score;
}

const char *aga_status_string(aga_status_t status)
{
    switch (status) {
    case AGA_OK:
        return "success";
    case AGA_INVALID_ARGUMENT:
        return "invalid argument";
    case AGA_SIZE_OVERFLOW:
        return "input is too large";
    case AGA_OUT_OF_MEMORY:
        return "out of memory";
    case AGA_INTERNAL_ERROR:
        return "internal traceback error";
    default:
        return "unknown error";
    }
}
