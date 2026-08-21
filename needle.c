/*
* Author: Xianxin Long
* Date: 2025-12-21
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>

#define max(a, b) ((a) > (b) ? (a) : (b))
#define MAX_LEN 10000 // the maximum input sequence length
#define INF 1e10

// Define AlignmentResult so Needleman_Wunsch can return the score and both aligned sequences
typedef struct {
    float score;
    char x[MAX_LEN * 2];
    char y[MAX_LEN * 2];
} AlignmentResult;

/*
* Function: Needleman_Wunsch
* Description: Standard implementation of Needleman Wunsch Algorithm for 2 sequence alignment using affine gap penalty function
* Input: v - the first sequence
*        w - the second sequence
* Output: None
* Returns: AlignmentResult
*/
AlignmentResult Needleman_Wunsch(char v[], char w[], float omega, float mu, float sigma, float epsilon){
    AlignmentResult result;
    int m = strlen(v);
    int n = strlen(w);
    int i, j;

    // Matrix initialization
    float lower[m + 1][n + 1], middle[m + 1][n + 1], upper[m + 1][n + 1];
    lower[0][0] = middle[0][0] = upper[0][0] = 0;
    for (i = 1; i <= m; i = i + 1){
        lower[i][0] = -sigma - (i-1) * epsilon;
        middle[i][0] = -INF;
        upper[i][0] = -INF;
    }
    for (j = 1; j <= n; j = j + 1){
        lower[0][j] = -INF;
        middle[0][j] = -INF;
        upper[0][j] = -sigma - (j-1) * epsilon;
    }

    // Fill the matrix
    for (i = 1; i <= m; i = i + 1){
        for (j = 1; j <= n; j = j + 1){
            lower[i][j] = max(
                lower[i-1][j] - epsilon,
                middle[i-1][j] - sigma 
            );
            upper[i][j] = max(
                upper[i][j-1] - epsilon,
                middle[i][j-1] - sigma
            );
            float score = (v[i-1] == w[j-1]) ? omega : -mu;
            middle[i][j] = max( 
                max(
                lower[i][j],
                middle[i-1][j-1] + score),
                upper[i][j]
            );
        }
    }
    
    // Calculate the NW score
    int ii = m, jj = n, state = 1, index = 0;
    char x[MAX_LEN * 2] = "", y[MAX_LEN * 2] = "";
    float NW_score = max(
        max(
        lower[m][n],
        middle[m][n]),
        upper[m][n]
    );
    result.score = NW_score;

    // Traceback
    if (lower[m][n] == NW_score){
        state = 0;
    }
    else if (middle[m][n] == NW_score){
        state = 1;
    }
    else if (upper[m][n] == NW_score){
        state = 2;
    }
    while (ii > 0 || jj > 0){
        // gap in w
        if (state == 0){
            x[index] = v[ii-1];
            y[index] = '-';
            if (lower[ii][jj] == lower[ii-1][jj] - epsilon){
                state = 0;
            }
            else if (lower[ii][jj] == middle[ii-1][jj] - sigma){
                state = 1;
            }
            ii = ii - 1;
        }
        // match
        else if (state == 1){
            if (middle[ii][jj] == lower[ii][jj]){
                state = 0;
                continue;
            }
            else if (middle[ii][jj] == upper[ii][jj]){
                state = 2;
                continue;
            }
            else{
                ii = ii - 1;
                jj = jj - 1;
                x[index] = v[ii];
                y[index] = w[jj];
                state = 1;
            }
        }
        // gap in v
        else if (state == 2) {
            x[index] = '-';
            y[index] = w[jj-1];
            if (upper[ii][jj] == upper[ii][jj-1] - epsilon){
                state = 2;
            }
            else if (upper[ii][jj] == middle[ii][jj-1] - sigma){
                state = 1;
            }
            jj = jj - 1;
        }
        index = index + 1;
    }

    // Reverse the strings
    int len = strlen(x);
    for (i = 0; i < len; i++) {
        result.x[i] = x[len - 1 - i];
        result.y[i] = y[len - 1 - i];
    }
    result.x[len] = '\0';
    result.y[len] = '\0';

    return result;
}

/*
* Function: alignmentIdentities
* Description: Calculate the identities in the alignment
* Input: result - the AlignmentResult
*        length - the alignment length
* Output: None
* Returns: The identical lengths of the two sequences in the alignment
*/
int alignmentIdentities(AlignmentResult result, int length) {
    int count = 0;
    for (int i = 0; i < length; i = i + 1) {
        if (result.x[i] == result.y[i]) {
            count = count + 1;
        }
    }
    return count;
}

/*
* Function: alignmentGaps
* Description: Calculate the gaps in the alignment
* Input: result - the AlignmentResult
*        length - the alignment length
* Output: None
* Returns: The gaps lengths of the two sequences in the alignment
*/
int alignmentGaps(AlignmentResult result, int length) {
    int count = 0;
    for (int i = 0; i < length; i = i + 1) {
        if (result.x[i] == '-' || result.y[i] == '-') {
            count = count + 1;
        }
    }
    return count;
}

/*
* Function: print_help
* Description: Output help information, including parameters and usage.
* Input: None
* Output: help information
* Returns: None
*/
void print_help(){
    fprintf(stderr, "Usage: needle -a <seq1> -b <seq2> -m <match_score> -M <mismatch_penalty> -o <gap_open> -e <gap_extend>\n");
    fprintf(stderr, "Example: needle -a ATGC -b ATCG -m 2 -M 1 -o 1 -e 0.5\n");
    return;
}

int main(int argc, char* argv[]){
    char seq1[MAX_LEN] = "";
    char seq2[MAX_LEN] = "";
    float match = 1, mismatch = 2, gap_open= 1, gap_extend = 1;

    // Parse the command line parameters
    int opt, a_flag = 0, b_flag = 0, m_flag = 0, M_flag = 0, o_flag = 0, e_flag = 0;
    while ((opt = getopt(argc, argv, "a:b:m:M:o:e:h")) != -1){
        switch (opt) {
            case 'a':
                strncpy(seq1, optarg, MAX_LEN - 1);
                seq1[MAX_LEN - 1] = '\0';
                a_flag = 1;
                break;
            case 'b':
                strncpy(seq2, optarg, MAX_LEN - 1);
                seq2[MAX_LEN - 1] = '\0';
                b_flag = 1;
                break;
            case 'm':
                match = atof(optarg);
                m_flag = 1;
                break;
            case 'M':
                mismatch = atof(optarg);
                M_flag = 1;
                break;
            case 'o':
                gap_open = atof(optarg);
                o_flag = 1;
                break;
            case 'e':
                gap_extend = atof(optarg);
                e_flag = 1;
                break;
            case 'h':
                print_help();
                return 0;
            default:
                fprintf(stderr, "Unknown option: %c\n", optopt);
                print_help();
                return 1;
        }
    }

    // Validate parameters
    int valid = 1;
    if (!a_flag || !b_flag || !m_flag || !M_flag || !o_flag || !e_flag){
        fprintf(stderr, "Error: Missing required parameters\n");
        if (!a_flag){fprintf(stderr, "Missing: -a <seq1>\n");}
        if (!b_flag){fprintf(stderr, "Missing: -b <seq2>\n");}
        if (!m_flag){fprintf(stderr, "Missing: -m <match_score>\n");}
        if (!M_flag){fprintf(stderr, "Missing: -M <mismatch_penalty>\n");}
        if (!o_flag){fprintf(stderr, "Missing: -o <gap_open>\n");}
        if (!e_flag){fprintf(stderr, "Missing: -e <gap_extend>\n");}
        valid = 0;
        print_help();
    }
    for (int i = 0; seq1[i] != '\0'; i = i + 1){
        if (seq1[i] > 'z' || seq1[i] < 'A' || (seq1[i] > 'Z' && seq1[i] < 'a')) {
            fprintf(stderr, "Error: Invalid character '%c' in sequence 1\n", seq1[i]);
            valid = 0;
        }
    }
    for (int i = 0; seq2[i] != '\0'; i = i + 1){
        if (seq2[i] > 'z' || seq2[i] < 'A' || (seq2[i] > 'Z' && seq2[i] < 'a')) {
            fprintf(stderr, "Error: Invalid character '%c' in sequence 2\n", seq2[i]);
            valid = 0;
        }
    }

    if (!valid) {
        return 1;
    }

    // Convert to uppercase letters
    for (int i = 0; seq1[i] != '\0'; i++) {
        if (seq1[i] >= 'a' && seq1[i] <= 'z') {
            seq1[i] = seq1[i] - 'a' + 'A';
        }
    }
    for (int i = 0; seq2[i] != '\0'; i++) {
        if (seq2[i] >= 'a' && seq2[i] <= 'z') {
            seq2[i] = seq2[i] - 'a' + 'A';
        }
    }

    // Align
    AlignmentResult result = Needleman_Wunsch(seq1, seq2, match, mismatch, gap_open, gap_extend);

    // Visualization
    printf("########################################\n");
    printf("# Program: needle\n");
    time_t t = time(NULL);
    printf("# Date: %s", ctime(&t)); 
    printf("# Commandline: needle -a %s -b %s -m %.1f -M %.1f -o %.1f -e %.1f\n", seq1, seq2, match, mismatch, gap_open, gap_extend);
    printf("########################################\n");

    printf("#=======================================\n");
    int length = strlen(result.x);
    printf("# Length: %d\n", length);
    printf("# Identity: %d\n", alignmentIdentities(result, length));
    printf("# Gaps: %d\n", alignmentGaps(result, length));
    printf("# Score: %.1f\n", result.score);
    printf("#=======================================\n");

    printf("Alignments:\n");
    printf("SeqA  1  %s  %zu\n", result.x, strlen(seq1));
    printf("         ");
    for (int i = 0; i < length; i = i + 1) {
        if (result.x[i] == result.y[i]) {
            printf("|");
        }
        else{
            printf(" ");
        }
    }
    printf("\n");
    printf("SeqB  1  %s  %zu\n", result.y, strlen(seq2));
    return 0;
}

/*
needle -a ATGAGTCTCTCTGATAAGGACAAGGCTGCTGTGAAAGCCCTATGG -b CTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAG -m 5 -M 4 -o 10 -e 0.5
*/
