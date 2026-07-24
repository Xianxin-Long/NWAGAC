#include "affine_align.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PROGRAM_VERSION "1.0.0"
#define DEFAULT_LINE_WIDTH 60

typedef enum {
    OUTPUT_TEXT,
    OUTPUT_JSON
} output_format_t;

typedef struct {
    const char *sequence_a_argument;
    const char *sequence_b_argument;
    const char *fasta_path;
    aga_scoring_t scoring;
    size_t line_width;
    output_format_t output_format;
} options_t;

static void print_usage(FILE *stream)
{
    fprintf(stream,
            "Usage:\n"
            "  affine-align --seq-a SEQUENCE --seq-b SEQUENCE [options]\n"
            "  affine-align --fasta FILE [options]\n\n"
            "Input:\n"
            "  -a, --seq-a SEQUENCE       First sequence\n"
            "  -b, --seq-b SEQUENCE       Second sequence\n"
            "  -f, --fasta FILE           FASTA file containing exactly two records\n\n"
            "Scoring (penalties are supplied as non-negative values):\n"
            "  -m, --match VALUE          Match score (default: 2)\n"
            "  -x, -M, --mismatch VALUE   Mismatch penalty (default: 1)\n"
            "  -o, --gap-open VALUE       Gap-open penalty (default: 3)\n"
            "  -e, --gap-extend VALUE     Gap-extension penalty (default: 1)\n\n"
            "Output:\n"
            "      --format text|json     Output format (default: text)\n"
            "  -w, --width INTEGER        Alignment columns per block (default: 60)\n"
            "  -h, --help                 Show this help message\n"
            "  -v, --version              Show program version\n\n"
            "A gap of length k costs gap_open + (k - 1) * gap_extend.\n");
}

static int option_takes_value(const char *argument)
{
    return strcmp(argument, "-a") == 0 ||
           strcmp(argument, "--seq-a") == 0 ||
           strcmp(argument, "-b") == 0 ||
           strcmp(argument, "--seq-b") == 0 ||
           strcmp(argument, "-f") == 0 ||
           strcmp(argument, "--fasta") == 0 ||
           strcmp(argument, "-m") == 0 ||
           strcmp(argument, "--match") == 0 ||
           strcmp(argument, "-x") == 0 ||
           strcmp(argument, "-M") == 0 ||
           strcmp(argument, "--mismatch") == 0 ||
           strcmp(argument, "-o") == 0 ||
           strcmp(argument, "--gap-open") == 0 ||
           strcmp(argument, "-e") == 0 ||
           strcmp(argument, "--gap-extend") == 0 ||
           strcmp(argument, "-w") == 0 ||
           strcmp(argument, "--width") == 0 ||
           strcmp(argument, "--format") == 0;
}

static int parse_double(const char *text, double *value)
{
    char *end = NULL;
    errno = 0;
    const double parsed = strtod(text, &end);

    if (errno == ERANGE || end == text || *end != '\0' ||
        !isfinite(parsed)) {
        return 0;
    }
    *value = parsed;
    return 1;
}

static int parse_size(const char *text, size_t *value)
{
    char *end = NULL;
    errno = 0;
    const unsigned long long parsed = strtoull(text, &end, 10);

    if (errno == ERANGE || end == text || *end != '\0' ||
        parsed == 0 || parsed > (unsigned long long)SIZE_MAX ||
        parsed > (unsigned long long)INT_MAX) {
        return 0;
    }
    *value = (size_t)parsed;
    return 1;
}

static int parse_options(int argc, char **argv, options_t *options)
{
    *options = (options_t){
        .sequence_a_argument = NULL,
        .sequence_b_argument = NULL,
        .fasta_path = NULL,
        .scoring = {2.0, 1.0, 3.0, 1.0},
        .line_width = DEFAULT_LINE_WIDTH,
        .output_format = OUTPUT_TEXT};

    for (int index = 1; index < argc; ++index) {
        const char *argument = argv[index];

        if (strcmp(argument, "-h") == 0 ||
            strcmp(argument, "--help") == 0) {
            print_usage(stdout);
            return 1;
        }
        if (strcmp(argument, "-v") == 0 ||
            strcmp(argument, "--version") == 0) {
            printf("affine-align %s\n", PROGRAM_VERSION);
            return 1;
        }

        if (!option_takes_value(argument)) {
            fprintf(stderr, "error: unknown option '%s'\n\n", argument);
            print_usage(stderr);
            return -1;
        }
        if (index + 1 >= argc) {
            fprintf(stderr, "error: option '%s' requires a value\n",
                    argument);
            return -1;
        }

        const char *value = argv[++index];
        if (strcmp(argument, "-a") == 0 ||
            strcmp(argument, "--seq-a") == 0) {
            options->sequence_a_argument = value;
        } else if (strcmp(argument, "-b") == 0 ||
                   strcmp(argument, "--seq-b") == 0) {
            options->sequence_b_argument = value;
        } else if (strcmp(argument, "-f") == 0 ||
                   strcmp(argument, "--fasta") == 0) {
            options->fasta_path = value;
        } else if (strcmp(argument, "-m") == 0 ||
                   strcmp(argument, "--match") == 0) {
            if (!parse_double(value, &options->scoring.match_score)) {
                fprintf(stderr, "error: invalid match score '%s'\n", value);
                return -1;
            }
        } else if (strcmp(argument, "-x") == 0 ||
                   strcmp(argument, "-M") == 0 ||
                   strcmp(argument, "--mismatch") == 0) {
            if (!parse_double(value,
                              &options->scoring.mismatch_penalty)) {
                fprintf(stderr,
                        "error: invalid mismatch penalty '%s'\n",
                        value);
                return -1;
            }
        } else if (strcmp(argument, "-o") == 0 ||
                   strcmp(argument, "--gap-open") == 0) {
            if (!parse_double(value,
                              &options->scoring.gap_open_penalty)) {
                fprintf(stderr,
                        "error: invalid gap-open penalty '%s'\n",
                        value);
                return -1;
            }
        } else if (strcmp(argument, "-e") == 0 ||
                   strcmp(argument, "--gap-extend") == 0) {
            if (!parse_double(value,
                              &options->scoring.gap_extend_penalty)) {
                fprintf(stderr,
                        "error: invalid gap-extension penalty '%s'\n",
                        value);
                return -1;
            }
        } else if (strcmp(argument, "-w") == 0 ||
                   strcmp(argument, "--width") == 0) {
            if (!parse_size(value, &options->line_width)) {
                fprintf(stderr, "error: invalid line width '%s'\n", value);
                return -1;
            }
        } else if (strcmp(argument, "--format") == 0) {
            if (strcmp(value, "text") == 0) {
                options->output_format = OUTPUT_TEXT;
            } else if (strcmp(value, "json") == 0) {
                options->output_format = OUTPUT_JSON;
            } else {
                fprintf(stderr,
                        "error: output format must be 'text' or 'json'\n");
                return -1;
            }
        }
    }

    if (options->scoring.mismatch_penalty < 0.0 ||
        options->scoring.gap_open_penalty < 0.0 ||
        options->scoring.gap_extend_penalty < 0.0) {
        fprintf(stderr, "error: penalties must be non-negative\n");
        return -1;
    }

    const int has_direct_input =
        options->sequence_a_argument != NULL ||
        options->sequence_b_argument != NULL;

    if (options->fasta_path != NULL && has_direct_input) {
        fprintf(stderr,
                "error: use either --fasta or --seq-a/--seq-b\n");
        return -1;
    }
    if (options->fasta_path == NULL &&
        (options->sequence_a_argument == NULL ||
         options->sequence_b_argument == NULL)) {
        fprintf(stderr,
                "error: provide --fasta or both --seq-a and --seq-b\n\n");
        print_usage(stderr);
        return -1;
    }

    return 0;
}

static int append_symbol(char **sequence,
                         size_t *length,
                         size_t *capacity,
                         char symbol)
{
    if (*length + 1 >= *capacity) {
        const size_t new_capacity =
            *capacity == 0 ? 128 : *capacity * 2;
        if (new_capacity <= *capacity) {
            return 0;
        }

        char *resized = realloc(*sequence, new_capacity);
        if (resized == NULL) {
            return 0;
        }
        *sequence = resized;
        *capacity = new_capacity;
    }

    (*sequence)[(*length)++] = symbol;
    (*sequence)[*length] = '\0';
    return 1;
}

static int read_two_record_fasta(const char *path,
                                 char **sequence_a,
                                 char **sequence_b)
{
    FILE *stream = fopen(path, "r");
    if (stream == NULL) {
        fprintf(stderr, "error: cannot open FASTA file '%s'\n", path);
        return 0;
    }

    char *sequences[2] = {NULL, NULL};
    size_t lengths[2] = {0, 0};
    size_t capacities[2] = {0, 0};
    int record = -1;
    int at_line_start = 1;
    int in_header = 0;
    int valid = 1;
    int character = 0;

    while ((character = fgetc(stream)) != EOF) {
        if (at_line_start && character == '>') {
            ++record;
            if (record >= 2) {
                fprintf(stderr,
                        "error: FASTA file must contain exactly two records\n");
                valid = 0;
                break;
            }
            in_header = 1;
            at_line_start = 0;
            continue;
        }

        if (character == '\n' || character == '\r') {
            at_line_start = 1;
            in_header = 0;
            continue;
        }

        if (in_header) {
            at_line_start = 0;
            continue;
        }

        if (isspace((unsigned char)character)) {
            continue;
        }
        at_line_start = 0;

        if (record < 0) {
            fprintf(stderr,
                    "error: FASTA sequence appears before the first header\n");
            valid = 0;
            break;
        }
        if (!isalpha((unsigned char)character)) {
            fprintf(stderr,
                    "error: invalid FASTA sequence character '%c'\n",
                    character);
            valid = 0;
            break;
        }
        if (!append_symbol(&sequences[record],
                           &lengths[record],
                           &capacities[record],
                           (char)toupper((unsigned char)character))) {
            fprintf(stderr, "error: out of memory while reading FASTA\n");
            valid = 0;
            break;
        }
    }

    if (ferror(stream)) {
        fprintf(stderr, "error: failed while reading FASTA file '%s'\n",
                path);
        valid = 0;
    }
    fclose(stream);

    if (valid && (record != 1 || lengths[0] == 0 || lengths[1] == 0)) {
        fprintf(stderr,
                "error: FASTA file must contain exactly two non-empty records\n");
        valid = 0;
    }

    if (!valid) {
        free(sequences[0]);
        free(sequences[1]);
        return 0;
    }

    *sequence_a = sequences[0];
    *sequence_b = sequences[1];
    return 1;
}

static char *normalize_sequence(const char *input, const char *label)
{
    const size_t length = strlen(input);
    char *normalized = malloc(length + 1);

    if (normalized == NULL) {
        fprintf(stderr, "error: out of memory\n");
        return NULL;
    }

    for (size_t index = 0; index < length; ++index) {
        const unsigned char symbol = (unsigned char)input[index];
        if (!isalpha(symbol)) {
            fprintf(stderr,
                    "error: invalid character '%c' in %s\n",
                    input[index],
                    label);
            free(normalized);
            return NULL;
        }
        normalized[index] = (char)toupper(symbol);
    }
    normalized[length] = '\0';
    return normalized;
}

static size_t count_residues(const char *text,
                             size_t offset,
                             size_t length)
{
    size_t count = 0;
    for (size_t index = offset; index < offset + length; ++index) {
        if (text[index] != '-') {
            ++count;
        }
    }
    return count;
}

static void print_text_alignment(const aga_alignment_t *alignment,
                                 size_t line_width)
{
    const aga_metrics_t metrics = aga_alignment_metrics(alignment);

    printf("Affine Global Alignment\n");
    printf("=======================\n");
    printf("Score          : %.6g\n", alignment->score);
    printf("Aligned length : %zu\n", alignment->length);
    printf("Identity       : %zu/%zu (%.1f%%)\n",
           metrics.identities,
           alignment->length,
           100.0 * metrics.identity_fraction);
    printf("Mismatches     : %zu\n", metrics.mismatches);
    printf("Gap characters : %zu\n", metrics.gap_characters);
    printf("Gap openings   : %zu\n\n", metrics.gap_openings);

    size_t position_a = 0;
    size_t position_b = 0;

    for (size_t offset = 0; offset < alignment->length;
         offset += line_width) {
        const size_t remaining = alignment->length - offset;
        const size_t block_length =
            remaining < line_width ? remaining : line_width;
        const size_t residues_a =
            count_residues(alignment->aligned_a, offset, block_length);
        const size_t residues_b =
            count_residues(alignment->aligned_b, offset, block_length);
        const size_t start_a = residues_a == 0 ? position_a : position_a + 1;
        const size_t start_b = residues_b == 0 ? position_b : position_b + 1;

        printf("SeqA %6zu  %.*s  %zu\n",
               start_a,
               (int)block_length,
               alignment->aligned_a + offset,
               position_a + residues_a);
        printf("             ");
        for (size_t index = 0; index < block_length; ++index) {
            const char symbol_a = alignment->aligned_a[offset + index];
            const char symbol_b = alignment->aligned_b[offset + index];

            if (symbol_a == '-' || symbol_b == '-') {
                putchar(' ');
            } else if (symbol_a == symbol_b) {
                putchar('|');
            } else {
                putchar('.');
            }
        }
        putchar('\n');
        printf("SeqB %6zu  %.*s  %zu\n\n",
               start_b,
               (int)block_length,
               alignment->aligned_b + offset,
               position_b + residues_b);

        position_a += residues_a;
        position_b += residues_b;
    }
}

static void print_json_alignment(const aga_alignment_t *alignment)
{
    const aga_metrics_t metrics = aga_alignment_metrics(alignment);

    printf("{\n");
    printf("  \"score\": %.17g,\n", alignment->score);
    printf("  \"aligned_length\": %zu,\n", alignment->length);
    printf("  \"identities\": %zu,\n", metrics.identities);
    printf("  \"identity_fraction\": %.17g,\n", metrics.identity_fraction);
    printf("  \"mismatches\": %zu,\n", metrics.mismatches);
    printf("  \"gap_characters\": %zu,\n", metrics.gap_characters);
    printf("  \"gap_openings\": %zu,\n", metrics.gap_openings);
    printf("  \"aligned_a\": \"%s\",\n", alignment->aligned_a);
    printf("  \"aligned_b\": \"%s\"\n", alignment->aligned_b);
    printf("}\n");
}

int main(int argc, char **argv)
{
    options_t options;
    const int parse_result = parse_options(argc, argv, &options);
    if (parse_result != 0) {
        return parse_result > 0 ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    char *sequence_a = NULL;
    char *sequence_b = NULL;

    if (options.fasta_path != NULL) {
        if (!read_two_record_fasta(options.fasta_path,
                                   &sequence_a,
                                   &sequence_b)) {
            return EXIT_FAILURE;
        }
    } else {
        sequence_a =
            normalize_sequence(options.sequence_a_argument, "sequence A");
        sequence_b =
            normalize_sequence(options.sequence_b_argument, "sequence B");
        if (sequence_a == NULL || sequence_b == NULL) {
            free(sequence_a);
            free(sequence_b);
            return EXIT_FAILURE;
        }
    }

    aga_alignment_t alignment;
    const aga_status_t status =
        aga_global_align(sequence_a,
                         sequence_b,
                         &options.scoring,
                         &alignment);
    free(sequence_a);
    free(sequence_b);

    if (status != AGA_OK) {
        fprintf(stderr, "error: alignment failed: %s\n",
                aga_status_string(status));
        return EXIT_FAILURE;
    }

    if (options.output_format == OUTPUT_JSON) {
        print_json_alignment(&alignment);
    } else {
        print_text_alignment(&alignment, options.line_width);
    }

    aga_alignment_free(&alignment);
    return EXIT_SUCCESS;
}
