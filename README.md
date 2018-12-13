# Gradiente Conjugado

## Sequencial

    gradiente.c

### Compilando

    gcc -o gradiente hb_io.c gradiente.c

### Executando

    ./gradiente

## Paralelo

    gradiente2_mpi.c

### Compilando

    mpicc hb_io.c gradiente2_mpi.c -o gradiente_mpi

### Executando

     mpirun -np 1 ./gradiente_mpi



# Observações

- Nome do arquivo da matriz de testes

	"bcsstruc2.data"
