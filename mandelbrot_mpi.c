#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <complex.h>

#define AXIS_LEN 1024
#define BILLION 1e9
#define MASTER 0

unsigned int matrix_checksum(int N, void *M, unsigned int size);


int iteratePoint(double x, double y, int cutoff){
    double complex c = x + y * I;
    double complex result = 0;
    int i;

    for(i = 0; (i < cutoff) && (cabs(result) < 2.0); i++){
        result = cpow(result, 2) + c;
    }

    return i;
}

void getInput(int argc, char *argv[], double *x_cen, double *y_cen, int *zoom,
                int *cutoff){
    *x_cen = atof(argv[1]);
    *y_cen = atof(argv[2]);
    *zoom = atoi(argv[3]);
    *cutoff = atoi(argv[4]);

    if(5 != argc){
        fprintf(stderr, "Usage: %s xcenter ycenter zoom cutoff", argv[0]);
        exit(1);
    }
    // Get and validate arguments
    if(10 < *x_cen || -10 > *x_cen){
        fprintf(stderr, "Error: wrong x-center (-10.000000 <= N <= 10.000000)\n");
        exit(1);
    }
    if(10 < *y_cen || -10 > *y_cen) {
        fprintf(stderr, "Error: wrong y-center (-10.000000 <= N <= 10.000000)\n");
        exit(1);
    }
    if(100 < *zoom || 0 > *zoom){
        fprintf(stderr, "Error: wrong zoom (0 <= N <= 100)\n");
        exit(1);
    }
    if(1000 < *cutoff || 50 > *cutoff){
        fprintf(stderr, "Error: wrong cutoff (50 <= N <= 1000)");
        exit(1);
    }
}

int main(int argc, char*argv[]){
    double x_center, y_center, dist;
    int zoom, cutoff, *plane, comm_size, comm_rank, section;
    struct timespec start, finish;

    // Initializing mpi
    MPI_Init(&argc,&argv);
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&comm_rank);
    if (comm_size != 2) {
        if (!comm_rank)
            MPI_Abort(MPI_COMM_WORLD, 1);
    }


    if(MASTER == comm_rank){        // Master code
        getInput(argc, argv, &x_center, &y_center, &zoom, &cutoff);

        dist = pow(2, -1 * zoom);

        section = AXIS_LEN / (comm_size - 1);
        if (0 < (AXIS_LEN % (comm_size - 1))){      // Remainder exists, so round up 1
            section += 1;
        }

        plane = malloc(sizeof(int) * AXIS_LEN * AXIS_LEN);

        clock_gettime(CLOCK_REALTIME, &start);
        for(int tid = 1; tid < comm_size; tid++){
            MPI_Send(&x_center, 1, MPI_DOUBLE, tid, 1, MPI_COMM_WORLD);
            MPI_Send(&y_center, 1, MPI_DOUBLE, tid, 1, MPI_COMM_WORLD);
            MPI_Send(&dist, 1, MPI_DOUBLE, tid, 1, MPI_COMM_WORLD);
            MPI_Send(&cutoff, 1, MPI_DOUBLE, tid, 1, MPI_COMM_WORLD);
            MPI_Send(&section, 1, MPI_INT, tid, 1, MPI_COMM_WORLD);
            MPI_Send(&plane, AXIS_LEN * AXIS_LEN, MPI_INT, tid, 1, MPI_COMM_WORLD);
        }

        for(int tid = 1; tid < comm_size; tid++){
            if(AXIS_LEN < tid * section){           // Only populates within the dimensions of the matrix
                MPI_Recv(&plane[(tid - 1) * section], section - ((tid * section) - AXIS_LEN),
                         MPI_INT, tid, 1, MPI_COMM_WORLD, &status);
            }
            else {
                MPI_Recv(&plane[(tid - 1) * section], section, MPI_INT, tid, 1,
                         MPI_COMM_WORLD, &status);
            }
        }
        clock_gettime(CLOCK_REALTIME, &finish);

        double time = (finish.tv_sec - start.tv_sec) +
                      (finish.tv_nsec - start.tv_nsec) / BILLION;
        printf("Runtime: %f\n", time);
        printf("M: %u\n", matrix_checksum(AXIS_LEN, plane, sizeof(int)));

        free(plane);
    }
    else{       // Slave code
        MPI_Recv(&x_center, 1, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD, &status);         // Receive parameters from master
        MPI_Recv(&y_center, 1, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&dist, 1, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&cutoff, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&section, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&plane, AXIS_LEN * AXIS_LEN, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);

        int start = (comm_rank - 1) * section,
            end = (start + section) - 1;

        for(int i = start; i < end; i++){
            for(int j = 0; j < AXIS_LEN; j++){
                plane[i * AXIS_LEN + j] = iteratePoint(dist * j + x_center, dist * i + y_center, cutoff);
            }
        }
        MPI_Send(&plane[start], AXIS_LEN * section, MPI_INT, MASTER, 2, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}
