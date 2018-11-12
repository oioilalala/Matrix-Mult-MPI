#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define AXIS_LEN 1024

int main(int argc, char*argv[]){
    double x_center, y_center;
    int zoom, cutoff, *plane;

    // Get and validate arguments
    x_center = atof(argv[1]);
    y_center = atof(argv[2]);
    zoom = atoi(argv[3]);
    cutoff = atoi(argv[4]);

    if(5 != argc){
        fprintf(stderr, "Usage: %s xcenter ycenter zoom cutoff", argv[0]);
        return 1;
    }
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

}