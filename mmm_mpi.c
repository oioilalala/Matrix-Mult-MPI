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
void matrixMultiplicationIKJ(int n, double *A, double *B, double **C, int myrow) {
 //   int i, j, k;
    
    for (int i = 0; i < myrow;
         i++) {
      for (int k = 0; k < n; k++) {
        double r = A[i * n + k];
        for (int j = 0; j < n; j++) {
          (*C)[i * n + j] += r * B[k * n + j];
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


int main(int argc, char **argv) {
    
    int n, numtasks, taskid, rowpertask, leftover, tid, offset, myrow;
  
    struct timespec t1, t2;
    double time_pass, time_sec, time_nsec;
    
    MPI_Init(&argc,&argv);
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
     
    if (taskid == MASTER) { 
        
        checkArgc(argc);
        readArgv(argv, &n);
  //      tasksz = n * n / (numtasks - 1);
  //      allsz = n * n;
     
        double *A = malloc(n * n * sizeof(double));
        double *B = malloc(n * n * sizeof(double));
        double *C = malloc(n * n * sizeof(double));
        
        initMatrix(n, &A, &B);
        
        rowpertask = n / (numtasks - 1);
        leftover = n % (numtasks - 1);
        offset = 0;
    
        clock_gettime(CLOCK_MONOTONIC, &t1);
        
        for (tid = 1; tid < numtasks; tid++) {
            myrow = (tid <= leftover) ? rowpertask + 1 : rowpertask;
            MPI_Send(&n, 1, MPI_INT, tid, 1, MPI_COMM_WORLD);
           // MPI_Send(&tasksz, 1, MPI_INT, tid, 1, MPI_COMM_WORLD);
            MPI_Send(&offset, 1, MPI_INT, tid, 1, MPI_COMM_WORLD);
            MPI_Send(&myrow, 1, MPI_INT, tid, 1, MPI_COMM_WORLD);
            MPI_Send(&B[0], n * n, MPI_DOUBLE, tid, 1, MPI_COMM_WORLD);
            MPI_Send(&A[offset], myrow * n, MPI_DOUBLE, tid, 1, MPI_COMM_WORLD);
            offset += myrow;
        }
     
        for (tid = 1; tid < numtasks; tid++) {
            MPI_Recv(&offset, 1, MPI_INT, tid, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&myrow, 1, MPI_INT, tid, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&C[offset], myrow * n, MPI_DOUBLE, tid, 2, MPI_COMM_WORLD, &status);
            
        }
        
        clock_gettime(CLOCK_MONOTONIC, &t2);
        time_sec = (double)(t2.tv_sec - t1.tv_sec);
        time_nsec = (double)(t2.tv_nsec - t1.tv_nsec);
        time_pass = time_sec + (time_nsec / 1000000000);
        displayResult(time_pass, n, A, B, C);
    
        free(A);
        free(B);
        free(C);
        
       
    }else {
     
        MPI_Recv(&n, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&offset, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&myrow, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
        
        double *A = malloc(n * myrow * sizeof(double));
        double *B = malloc(n * n * sizeof(double));
        double *C = malloc(n * myrow * sizeof(double));
        
        MPI_Recv(&B[0], n * n, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&A[0], n * myrow, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD, &status);
        matrixMultiplicationIKJ(n,A,B,&C,myrow);
        MPI_Send(&offset, 1, MPI_INT, MASTER, 2, MPI_COMM_WORLD);
        MPI_Send(&myrow, 1, MPI_INT, MASTER, 2, MPI_COMM_WORLD);
        MPI_Send(&C[0], n * myrow, MPI_DOUBLE, taskid, 2, MPI_COMM_WORLD);
        
    }
    MPI_Finalize();
     
    return 0;
}
