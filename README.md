# Algorithm Conjugate Gradient

Conjugate Gradient Method is an algorithm for the numerical solution of particular systems of linear equations, namely those whose matrix is positive-definite

## Sequential Implementation

    gradiente.c

### Build

    gcc -o gradiente hb_io.c gradiente.c

### Run

    ./gradiente

## Parallel Implementation - MPI

    gradiente2_mpi.c

### Build

    mpicc hb_io.c gradiente2_mpi.c -o gradiente_mpi

### Run

     mpirun -np 1 ./gradiente_mpi



# Comments

- File name of the test matrix

	"bcsstruc2.data"
