# Compilando

gcc -o gradiente hb_io.c gradiente1.c

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



# Paralelo

mpicc hb_io.c gradiente2.c -o gradiente

mpirun -np 2 ./gradiente
