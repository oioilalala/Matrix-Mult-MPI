#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>

unsigned int matrix_checksum(int N, void *M, unsigned int size);

#define AXIS_LEN 1024
#define BILLION 1e9

int iteratePoint(double x, double y, int cutoff){
    double complex c = x + y * I;
    double complex result = 0;
    int i;
    for(i = 0; (i < cutoff) && (cabs(result) < 2.0); i++){
        result = cpow(result, 2) + c;
    }

    return i;
}

void calcPlane(int *plane, int x_center, int y_center, int zoom, int cutoff){
    struct timespec start, finish;
    double dist = pow(2, -1 * zoom);

    // Iterate through the 1024x1024 points in the map, passing x and y
    // to iteratePoint() to determine Mandelbrot membership for each
    clock_gettime(CLOCK_REALTIME, &start);
    for(int i = 0; i < AXIS_LEN; i++){
        for(int j = 0; j < AXIS_LEN; j++){
            plane[i * AXIS_LEN + j] = iteratePoint(dist * j, dist * i, cutoff);
        }
    }
    clock_gettime(CLOCK_REALTIME, &finish);

    double time = (finish.tv_sec - start.tv_sec) +
                  (finish.tv_nsec - start.tv_nsec) / BILLION;
    printf("Runtime: %f\n", time);
    printf("M: %u\n", matrix_checksum(AXIS_LEN, plane, sizeof(int)));
}

int main(int argc, char*argv[]){
    double x_center, y_center;
    int zoom, cutoff, *plane;

    if(5 != argc){
        fprintf(stderr, "Usage: %s xcenter ycenter zoom cutoff", argv[0]);
        return 1;
    }

    // Get and validate arguments
    x_center = atof(argv[1]);
    y_center = atof(argv[2]);
    zoom = atoi(argv[3]);
    cutoff = atoi(argv[4]);

    if(10 < x_center || -10 > x_center){
        fprintf(stderr, "Error: wrong x-center (-10.000000 <= N <= 10.000000)\n");
        return 1;
    }
    if(10 < y_center || -10 > y_center) {
        fprintf(stderr, "Error: wrong y-center (-10.000000 <= N <= 10.000000)\n");
        return 1;
    }
    if(100 < zoom || 0 > zoom){
        fprintf(stderr, "Error: wrong zoom (0 <= N <= 100)\n");
        return 1;
    }
    if(1000 < cutoff || 50 > cutoff){
        fprintf(stderr, "Error: wrong cutoff (50 <= N <= 1000)");
        return 1;
    }

    plane = malloc(sizeof(int) * AXIS_LEN * AXIS_LEN);

    calcPlane(plane, x_center, y_center, zoom, cutoff);
}