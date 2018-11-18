#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <complex.h>

#define RES 1024
#define BILLION 1e9
#define MASTER 0


unsigned int matrix_checksum(int N, void *M, unsigned int size);


int iteratePoint(double x, double y, int cutoff){
    double complex c = x + y * I;
    double complex result = 0;
    int i = 0;

    while((i < cutoff) && (cabs(result) <= 2.0)){
        result = result * result + c;
        i++;
    }

    return i;
}

void writeImage(int *pgmdata, double x, double y,
                int zoom, int cutoff){
    char buf[0x40];
    int temp = 0;
    snprintf(buf, sizeof(buf), "mandel_%f_%f_%d_%d.pgm",x,y,zoom,cutoff);
    FILE *outPGM;
    outPGM = fopen(buf,"w");
    fprintf(outPGM, "P2\n");
    fprintf(outPGM, "1024 1024\n");
    fprintf(outPGM, "%d\n",cutoff);
    for (int i=0; i<1024; i++){
        for (int j=0; j<1024; j++){
            temp = pgmdata[i*1024+j];
            fprintf(outPGM, "%d ", temp);
        }
        fprintf(outPGM, "\n");
    }
    fclose(outPGM);

}


void getInput(int argc, char *argv[], double *x_cen, double *y_cen, int *zoom,
              int *cutoff){
    *x_cen = atof(argv[1]);
    *y_cen = atof(argv[2]);
    *zoom = atoi(argv[3]);
    *cutoff = atoi(argv[4]);

    // Get and validate arguments
    if (argc != 5) {
        fprintf(stderr, "Usage: %s xcenter ycenter zoom cutoff", argv[0]);
        exit(1);
    }
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
    int zoom, cutoff, comm_size, comm_rank, next;
    struct timespec start, finish;

    // Initializing mpi
    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&comm_rank);

    if(MASTER == comm_rank){        // Master code
        int working_tasks,next,source_id;
        int tid;
        working_tasks = 0;
        source_id =0 ;
        next = 0;                   // Row
        getInput(argc, argv, &x_center, &y_center, &zoom, &cutoff);

        int *plane = malloc(sizeof(int) * RES * RES);
        dist = pow(2, -1 * zoom);
        double x = x_center;
        double y = y_center;

        x_center -= 512 * dist;
        y_center -= 512 * dist;

        clock_gettime(CLOCK_REALTIME, &start);

        for (tid=1;tid<comm_size;tid++) {                                       /* initial distribution of work */
            MPI_Send(&x_center, 1, MPI_DOUBLE, tid, 1, MPI_COMM_WORLD);
            MPI_Send(&y_center, 1, MPI_DOUBLE, tid, 1, MPI_COMM_WORLD);
            MPI_Send(&dist, 1, MPI_DOUBLE, tid, 1, MPI_COMM_WORLD);
            MPI_Send(&cutoff, 1, MPI_INT, tid, 1, MPI_COMM_WORLD);
            MPI_Send(&next, 1, MPI_INT, tid, 1, MPI_COMM_WORLD);
            next+=2;
            working_tasks++;
        }


        while(working_tasks>0){                                                 /* dynamic distribution */
            MPI_Recv(&next, 1, MPI_INT, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);
            source_id = status.MPI_SOURCE;
            MPI_Recv(&plane[next * RES], RES, MPI_INT, source_id, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&plane[(next+1) * RES], RES, MPI_INT, source_id, 4, MPI_COMM_WORLD, &status);
            working_tasks--;

            if(next < RES-2) {
                next+=2;
                MPI_Send(&next, 1, MPI_INT, source_id, 1, MPI_COMM_WORLD);
                working_tasks++;
            }  else {
                MPI_Send(&next, 0, MPI_INT, source_id, 3, MPI_COMM_WORLD);
            }
        }

        clock_gettime(CLOCK_REALTIME, &finish);

        double time = (finish.tv_sec - start.tv_sec) +
                      (finish.tv_nsec - start.tv_nsec) / BILLION;
        printf("Runtime: %f\n", time);
        printf("M: %u\n", matrix_checksum(RES, plane, sizeof(int)));
        writeImage(plane, x, y, zoom, cutoff);

        free(plane);

    }
    else{       // Slave code
        int *plane = malloc(sizeof(int) * RES * 2);

        MPI_Recv(&x_center, 1, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&y_center, 1, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&dist, 1, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&cutoff, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);

        while ( (MPI_Recv(&next, 1, MPI_INT, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &status)==MPI_SUCCESS) && (status.MPI_TAG == 1)) {
            for (int i = 0; i < RES; i++){
                plane[i]=iteratePoint(dist * next + x_center, dist * i + y_center, cutoff);
                plane[i+RES]=iteratePoint(dist * (next+1) + x_center, dist * i + y_center, cutoff);
            }

            MPI_Send(&next, 1, MPI_INT, MASTER, 2, MPI_COMM_WORLD);
            MPI_Send(&plane[0], RES, MPI_INT, MASTER, 2, MPI_COMM_WORLD);
            MPI_Send(&plane[RES], RES, MPI_INT, MASTER, 4, MPI_COMM_WORLD);

        }

        free(plane);
    }

    MPI_Finalize();

    return 0;
}
