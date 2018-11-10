#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#define MASTER 0
unsigned int matrix_checksum(int N, double *M);

/*
 * This function is to check if has correct number of arguments provided
 */
void checkArgc(int argc) {
    if (argc != 2) {
        fprintf(stderr, "Usage: ./mmm_mpi N\n");
        exit(1);
    }
}

/*
 * This function is to check if the arguments provided are correct, if correct,
 * read the arguments, otherwise report error
 */
void readArgv(char **argv, int *n) {
    *n = atoi(argv[1]);
    if (*n <= 0 ) {
        fprintf(stderr, "Error: wrong matrix order (0 < N)\n");
        exit(1);
    }
}

/*
 * This function is to initialize the matrices
 */
void initMatrix(int n, double **A, double **B) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            (*A)[i * n + j] = i + j;
            (*B)[i * n + j] = i + j * 2;
        }
    }
}


// ikj
void matrixMultiplicationIKJ(int n, int taskid, int numtasks, double *A, double *B, double **C) {
    int i, j, k;
    for (i = taskid * n/numtasks; i < (taskid + 1) * n/numtasks; i++) {
        for (k = 0; k < n; k++) {
            double r = A[i * n + k];
            for (j = 0; j < n; j++) {
                (*C)[i * n + j] += r * B[k * n + j];
                    //(*C)[i*n+j]=(*C)[i*n+j]+r*B[k*n+j];
            }
        }
    }
    
}


// display the time and checksum
void displayResult(double time, int N, double *A, double *B, double *C) {
    printf("Running time: %f secs\n", time);
    printf("A: %u\n", matrix_checksum(N, A));
    printf("B: %u\n", matrix_checksum(N, B));
    printf("C: %u\n", matrix_checksum(N, C));
}

//void srandom (unsigned seed);

int main(int argc, char **argv) {
    int n,numtasks,taskid;
    struct timespec t1, t2;
    double time_pass, time_sec, time_nsec;
    checkArgc(argc);
    readArgv(argv, &n);
    double *A = malloc(n * n * sizeof(double));
    double *B = malloc(n * n * sizeof(double));
    double *C = malloc(n * n * sizeof(double));
    initMatrix(n, &A, &B);
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  
    clock_gettime(CLOCK_MONOTONIC, &t1);
   
    matrixMultiplicationIKJ(n, taskid, numtasks, A, B, &C);
    
    clock_gettime(CLOCK_MONOTONIC, &t2);
    
    time_sec = (double)(t2.tv_sec - t1.tv_sec);
    time_nsec = (double)(t2.tv_nsec - t1.tv_nsec);
    time_pass = time_sec + (time_nsec / 1000000000);
    displayResult(time_pass, n, A, B, C);
    
    free(A);
    free(B);
    free(C);
    MPI_Finalize();
    return 0;
}
