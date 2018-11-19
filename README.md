# Project 3
#### Oi Lam Sou and Keith Luong

### Phases
1. Program #1: Matrix multiplication using MPI
2. Program #2: Mandelbrot set

### Program #1: Matrix multiplication using MPI
The performance of matrix multiplication this time is optimized by MPI. The idea is quite similar to theading, except that it is done so through message passing. This main characteristics of MPI allows the use of remote processors. As a result, the program becomes more scalable as it's no longer bounded by a local host. In `mmm_mpi.c`, we initialized MPI by   `MPI_Init(&argc,&argv);`, then we obtain the size of the communicator created and the rank of a task by `MPI_Comm_size(MPI_COMM_WORLD,&numtasks);` and `MPI_Comm_rank(MPI_COMM_WORLD,&taskid);` respectively. Task with rank `0` will be the master while tasks with ranks larger than `0` will be the workers/slaves. The distribution of workload is static in this program. Each worker will be assigned the same amount of workload `rowpertask = n / (numtasks - 1)` at once. If the matrix order is not divisible by number of tasks, the remainder will be made up by adding one more row. 

#### Master task
`MASTER` is responsible for allocating memory space for all matrices and initializing `A` and `B` so that it can send them to the workers. It starts by sending all parameters and input needed for computation to workers, including the matrix order `n`, `offset` for `A` and `C`, segment of rows `myrow`, the whole matrix `B` and the respective segment of matrix `A`. Messages from `MASTER` to workers are tagged with `1`.
```c
...
MPI_Send(&B[0], n * n, MPI_DOUBLE, tid, 1, MPI_COMM_WORLD);
MPI_Send(&A[offset], myrow * n, MPI_DOUBLE, tid, 1, MPI_COMM_WORLD);
...
```
Then, `MASTER` will receive the results from workers including a segment of `C`, size of the segment (number of rows) `myrow` and `offset` in order to put the segment of `C` in the right place:
```c
...
MPI_Recv(&offset, 1, MPI_INT, tid, 2, MPI_COMM_WORLD, &status);
MPI_Recv(&myrow, 1, MPI_INT, tid, 2, MPI_COMM_WORLD, &status);
MPI_Recv(&C[offset], myrow * n, MPI_DOUBLE, tid, 2, MPI_COMM_WORLD, &status);
...
```
Messages from workers to `MASTER` is tagged as `2`.
#### Worker tasks
A worker is responsible for the actual computation of matrix multiplication. First it receives the parameters needed for computation: `n` (matrix order), `offset` (for `A` and `C`) and `myrow` (number of rows). Using this information, worker is able to allocate appropriate memory space for the matrices. Then, it will receive the whole `B` and the portion of `A` in order to calculate its portion of `C`. Finally, it will send the calculated portion of `C` and its respective parameters to `MASTER` so that the portion of `C` will be put in the correct place. 

MPI is terminated by `MPI_Finalize();`. 

#### Timing measurement
`mpirun -n N ./mandelbrot_mpi -1 -0.2 9 256`

| Number of tasks  | Running time #1 | Running time #2 | Running time #3 | Mean runtime    |
| ---------------- |:---------------:|:---------------:|:---------------:|:---------------:| 
| 2 (serial)       | 0.065661        | 0.066007        | 0.076338        | 0.069335        |
| 4                | 0.037143        | 0.036668        | 0.025774        | 0.033195        |
| 8                | 0.032724        | 0.034662        | 0.032600        | 0.033329        |
| 16               | 0.066835        | 0.050057        | 0.051531        | 0.056141        |
| 32               | 0.117624        | 0.109797        | 0.118193        | 0.115205        |
| 64               | 0.436475        | 0.405924        | 0.382945        | 0.408448        |

### Program #2: Mandelbrot Set



#### Timing Measurement
`mpirun -n N ./mandelbrot_mpi -1 -0.2 9 256`

| Number of tasks  | Running time #1 | Running time #2 | Running time #3 | Mean runtime    |
| ---------------- |:---------------:|:---------------:|:---------------:|:---------------:| 
| 2 (serial)       | 0.847501        | 0.844279        | 0.844985        | 0.845588        |
| 4                | 0.307523        | 0.303908        | 0.304338        | 0.305256        |
| 8                | 0.253007        | 0.251352        | 0.251570        | 0.251976        |
| 16               | 0.260951        | 0.262446        | 0.246086        | 0.256494        |
| 32               | 0.245656        | 0.236120        | 0.240654        | 0.240810        |
| 64               | 0.295983        | 0.250135        | 0.294961        | 0.280360        |


### How we tested the project
We ran the test scripts provided to see if our output matches the expected output.

### Sources used
1. Project 1 and 2
