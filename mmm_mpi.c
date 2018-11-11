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
void matrixMultiplicationIKJ(int n, double *A, double *B, double **C, int taskid, int numtasks) {
 //   int i, j, k;
    
    for (int i = taskid * n / numtasks; i < (taskid + 1) * n / numtasks;
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
    int n, numtasks, taskid, tasksz, tid, allsz;
    struct timespec t1, t2;
    double time_pass, time_sec, time_nsec;
    
      printf("0");
    MPI_Init(&argc,&argv);
        printf("0");
    MPI_Status status;
        printf("0");
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
        printf("0");
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
        printf("0");
   
      printf("1");
     
    clock_gettime(CLOCK_MONOTONIC, &t1);
    
    if (taskid == MASTER) {
        checkArgc(argc);
        readArgv(argv, &n);
        tasksz = n * n / numtasks;
        allsz = n * n;
        
          printf("2");
     
        double *A = malloc(n * n * sizeof(double));
        double *B = malloc(n * n * sizeof(double));
        double *C = malloc(n * n * sizeof(double));
        
        initMatrix(n, &A, &B);
          printf("3");
     
        clock_gettime(CLOCK_MONOTONIC, &t1);
        for (tid = 1; tid < numtasks; tid++) {
           
            MPI_Send(&n, 1, MPI_INT, tid, 1, MPI_COMM_WORLD);
            MPI_Send(&tasksz, 1, MPI_INT, tid, 1, MPI_COMM_WORLD);
            MPI_Send(&B[0], allsz, MPI_DOUBLE, tid, 1, MPI_COMM_WORLD);
            MPI_Send(&A[tid * tasksz], tasksz, MPI_DOUBLE, tid, 1, MPI_COMM_WORLD);
        }
          printf("4");
     
        for (tid = 1; tid < numtasks; tid++) {
            MPI_Recv(&C[tid * tasksz], tasksz, MPI_DOUBLE, tid, 2, MPI_COMM_WORLD, &status);
            int sender = status.MPI_SOURCE;

            printf("%f returned from process %d\n", C[tid * tasksz], sender);
     
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
        MPI_Recv(&tasksz, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
        
        double *A = malloc(n * n * sizeof(double) / tasksz);
        double *B = malloc(n * n * sizeof(double));
        double *C = malloc(n * n * sizeof(double) / tasksz);
        
        MPI_Recv(&B[0], n * n, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&A[taskid * tasksz], tasksz, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD, &status);
        matrixMultiplicationIKJ(n,A,B,&C,taskid,numtasks);
        MPI_Send(&C[0], tasksz, MPI_DOUBLE, taskid, 2, MPI_COMM_WORLD);
        
    }
    MPI_Finalize();
    //signals
    //MPI reduce??
    
    return 0;
}
