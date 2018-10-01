#include <stdio.h>
#include <stdlib.h>
#include "mmio.c"

int main(int argc, char *argv[]){

	int ret_code;
	MM_typecode matcode;
	FILE *f;
	int M, N, nz;
	int i, *I, *J;
	double *val;

	// if (argc < 2){
	// 	fprintf(stderr, "%s < arquivo matriz >\n", argv[0]);
	// 	return 0;
	// }

	f = fopen("boing.rsa", "r");

	if ( f == NULL ){
		printf("Erro ao abrir o arquivo\n");
		return 0;
	}

	// Le os Banner do arquivo

	 mm_read_banner(f, &matcode);

	// Le as dimensoes da matriz

	mm_read_mtx_crd_size(f, &M, &N, &nz);

	// Aloca memoria para armazenar a matriz

	I = (int *) malloc(nz * sizeof(int));
  J = (int *) malloc(nz * sizeof(int));
  val = (double *) malloc(nz * sizeof(double));

	// Le a matriz

	for (i=0; i<nz; i++){
		fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        	I[i]--;  /* Ajusta indices de 1 para 0 */
        	J[i]--;
    	}

	// Escreve a matriz na tela

	mm_write_banner(stdout, matcode);
	mm_write_mtx_crd_size(stdout, M, N, nz);
  for (i=0; i<nz; i++){
        	fprintf(stdout, "%d %d %20.19g\n", I[i]+1, J[i]+1, val[i]);
	}

	return 0;
}
