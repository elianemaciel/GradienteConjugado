# Gradiente Conjugado

## Sequencial

    gradiente.c

### Compilando

    gcc -o gradiente hb_io.c gradiente.c

### Executando

    ./gradiente

## Sequencial em funções

    gradiente2.c

### Compilando

    gcc -o gradiente hb_io.c gradiente2.c

### Executando

    ./gradiente2

## Paralelo

    gradiente2_mpi.c

### Compilando

    mpicc hb_io.c gradiente2_mpi.c -o gradiente_mpi

### Executando

    ./gradiente2_mpi

## Paralelo OpenMP (sem utilizar)

    gcc -o gradiente hb_io.c gradiente2.c -fopenmp

# Executando

    ./gradiente


# Observações

- Nome do arquivo da matriz de testes

	"bcsstruc2.data"

- Nome do arquivo vetor

    "vetor.txt"

# Perfilação de Código

- Resultados:

	saida.txt
