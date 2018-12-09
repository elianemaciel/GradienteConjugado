# Compilando

gcc -o gradiente hb_io.c gradiente.c

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



# Paralelo (não UTILIZADO)

mpicc hb_io.c gradiente2.c -o gradiente

mpirun -np 2 ./gradiente

# Paralelo OpenMP

gcc -o gradiente hb_io.c gradiente2.c -fopenmp

# Executando

./gradiente
