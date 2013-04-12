#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define M 64
#define N 64
#define MAX(x,y) (((x) > (y)) ? (x) : (y))

int* allocate_dp_mat() {
    return (int*) calloc(N*M, sizeof(int));
}

int compute_align_score(int* dp_mat, char* tSeq, char* qSeq) {
    int ipenalty   = -1;
    int dpenalty   = -1;
    int match      =  2;
    int mpenalty   = -2;
    int best_score = 0;
    int iscore     = 0;
    int dscore     = 0;
    int mscore     = 0;
    int i,j;

    memset(dp_mat, 0, M*N*sizeof(int));

    for (i = 1; i < strlen(tSeq) + 1; i++) {
	for (j = 1; j < strlen(qSeq) + 1; j++) {
	    iscore = dp_mat[i*M + j-1] + ipenalty;
	    dscore = dp_mat[(i-1)*M + j] + dpenalty;
	    mscore = dp_mat[(i-1)*M + j-1] + ((tSeq[i-1] == qSeq[j-1]) ? match : mpenalty);
	    dp_mat[i*M + j] = MAX(MAX(0, iscore), MAX(dscore, mscore));
 	    if (dp_mat[i*M + j] >= best_score) 
		best_score = dp_mat[i*M + j];
	}
    }
    return best_score;
}

void compute_align_scores(int* scores, int n, int* dp_mat, char* tSeq, 
                          char** qSeqs) {
    int i = 0;
    for (i; i < n; i++) {
        scores[i] = compute_align_score(dp_mat, tSeq, qSeqs[i]);
    }
}


void print_dp_mat(int* dp_mat, char* tSeq, char* qSeq) {
    int i,j;
    for (j = 0; j < strlen(qSeq) + 1; j++) {
    	for (i = 0; i < strlen(tSeq) + 1; i++) {
    	    printf("%d ", dp_mat[i*M + j]);
    	}
    	printf("\n");
    }
}
