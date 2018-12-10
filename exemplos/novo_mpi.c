#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>

/******************************************************************************************************************************************/
void abre_arquivo(FILE **arquivo, char *nome_arquivo){
	if ((*arquivo=fopen(nome_arquivo,"r")) == NULL){
	    printf("Erro na abertura do arquivo %s \n",nome_arquivo);
	    exit(0);
	}
	else{
		printf("Arquivo %s aberto com sucesso\n",nome_arquivo);
	}
}

/******************************************************************************************************************************************/
void fecha_arquivo(FILE **aux){
     fclose(*aux);
}

/******************************************************************************************************************************************/
int leitura_vetor_valor(FILE **arq, float *vetor, int posicao){
	char linha[500], aux[40];  
    register int i=0, p=0; 	
	int  encontrouNumero = 0;
    
	fgets(linha, 500, *arq); // leitura da linha
        
    while(linha[i] != '\0'){
    	if(linha[i] != ' ' && linha[i] != '\n'){
	 		encontrouNumero = 1;
			aux[p++] = linha[i];	
		}
		else{
			if(encontrouNumero){
				encontrouNumero = 0;
				aux[p]='\0';
				p=0;
				
				vetor[posicao++] = atof(aux);
			}				
		}
		i++;
	}

	return posicao;
}

/******************************************************************************************************************************************/
int leitura_vetor(FILE **arq, int *vetor, int posicao){
    char linha[500], aux[40];  
    register int i=0, p=0;
    int encontrouNumero = 0;
    
    fgets(linha, 500, *arq); // leitura da linha
    
    while(linha[i] != '\0'){
    	if(linha[i] != ' ' && linha[i] != '\n'){
	 		encontrouNumero = 1;
			aux[p++] = linha[i];	
		}
		else{
			if(encontrouNumero){
				encontrouNumero = 0;
				aux[p]='\0';
				p=0;
				
				vetor[posicao++] = atoi(aux);
			}				
		}
		i++;
	}
    	
	return posicao;
}

/******************************************************************************************************************************************/
void leitura_cabecalho(FILE **arq, int *tam_linha, int *tam_coluna, int *tam_valor, int *N, int *nao_zeros){
    char linha[500], aux[40];  
    register int i=0, p=0;
    
    fgets(linha,500,*arq); //ignora a primeira linha
    fgets(linha,500,*arq); //leitura da segunda linha

    while(linha[i] == ' ')
    	i++;
    
    while(linha[i] != ' ')
    	i++;

    while(linha[i] == ' ')
    	i++;

    while(linha[i] != ' ')
    	aux[p++] = linha[i++];

    aux[p] = '\0';
    p = 0;
    *tam_coluna = atoi(aux);


    while(linha[i] == ' ')
    	i++;

    while(linha[i] != ' ')
    	aux[p++] = linha[i++];

    aux[p] = '\0';
    p=0;
    *tam_linha = atoi(aux);
    

    while(linha[i] == ' ')
    	i++;

    while(linha[i] != ' ')
    	aux[p++] = linha[i++];
    
    aux[p] = '\0';
    p=0;
    *tam_valor = atoi(aux);    

    fgets(linha,500,*arq); //leitura da terceira linha

    i=0;
    while(linha[i] != ' ')
    	i++;

    while(linha[i] == ' ')
    	i++;

    while(linha[i] != ' ')
    	i++;

    while(linha[i] == ' ')
    	i++;

    while(linha[i] != ' ')
    	aux[p++] = linha[i++];

    aux[p] = '\0';
    p = 0;
    *N = atoi(aux);   

    while(linha[i] == ' ')
    	i++;

    while(linha[i] != ' ')
    	aux[p++] = linha[i++];

    aux[p] = '\0';
    p = 0;
    *nao_zeros = atoi(aux);

 	fgets(linha,500,*arq); //ignora a quarta linha
}

/******************************************************************************************************************************************/
void inicializa_vetor(float *vetor, int N, int valor){
	register int i;
	for(i=0; i<N; i++){
 		vetor[i] = valor;
	}
}

/******************************************************************************************************************************************/
void escreve_vetor(int tam_coluna, int tam_linha, int tam_valor, int N, int nao_zeros, int quant_vetor_coluna, int *colptr, int *linhas, float *valores, float *vetX){
	register int i;
	
	printf("\nCOLUNA: %d - LINHA: %d - VALOR: %d\n", tam_coluna, tam_linha, tam_valor);
    printf("DIMENSAO: %d - VALORES: %d\n", N, nao_zeros);
			
	printf("\nVETOR COLUNA:\n");
    for(i=0; i<quant_vetor_coluna; i++)
		printf("%d  ", colptr[i]);
				
	printf("\n\nVETOR LINHA:\n");
	for(i=0; i<nao_zeros; i++)
		printf("%d  ", linhas[i]);
		
	printf("\n\nVETOR VALOR:\n");
	for(i=0; i<nao_zeros; i++)
		printf("%f  ", valores[i]);
		
	printf("\n\nVETOR RESPOSTA:\n");
	for(i=0; i<N; i++){
		printf("%f\n", vetX[i]);
	}	
}

/******************************************************************************************************************************************/
void Mutiplicacao_Matriz_Vetor(int *colptr, int *linhas, float *valores, float *vet, float *vet_res, int N, int nao_zeros, int id, int np){
	register int i=0, j=0;
	
	for(i=id; i<N; i+=np){				
		for(j=(colptr[i]-1); j<(colptr[i+1]-1); j++){
			vet_res[i] += valores[j] * vet[linhas[j]-1];	
				
			if(i != (linhas[j]-1)){
				vet_res[linhas[j]-1] += valores[j] * vet[i];	
			}		
		}		
   	}
}

/******************************************************************************************************************************************/
void Subtracao_Vetor_Vetor(float *vet_1, float *vet_2, float *vet_res, int N){
	register int i;
	
	for(i=0; i<N; i++){
		vet_res[i] = vet_1[i] - vet_2[i];
	}
}

/******************************************************************************************************************************************/
void Atibuicao_Vetor(float *vet_1, float *vet_res, int N){
	register int i;
	
	for(i=0; i<N; i++){
		vet_res[i] = vet_1[i];
	}
}

/******************************************************************************************************************************************/
float Multiplicacao_Vetor_Vetor_Soma(float valor, float *vet_1, float *vet_2, int N){
	register int i;
	
	for(i=0; i<N; i++){
		valor = valor + (vet_1[i] * vet_2[i]);
	}
	
	return valor;
}

/******************************************************************************************************************************************/
void gradiente(int *colptr, int *linhas, float *valores, float *vetB, float *vetX, int N, int nao_zeros, int quant_vetor_coluna, int id, int np){
	register int i;
	
	float sigma_novo = 0;
	float sigma = 0;
	float alpha_linha = 0;
	float alpha = 0;
	float sigma_velho = 0;
	float sigma_linha = 0;	
	float beta = 0;

	float imax = 1000;
	float erro = 0.00001;
	register int iteracao = 1;

	float *vet_Multiplicacao = NULL;
	float *vetR = NULL;
	float *vetD = NULL;
	float *vetQ = NULL;
	float *vet_Aux = NULL;
	float *vet_Soma = NULL;
	
	vet_Multiplicacao = (float*) malloc(N * sizeof(float));
	vetR = (float*) malloc(N * sizeof(float));
	vetD = (float*) malloc(N * sizeof(float));
	vetQ = (float*) malloc(N * sizeof(float));
	vet_Aux = (float*) malloc(N * sizeof(float));
	vet_Soma = (float*) malloc(N * sizeof(float));
	
	inicializa_vetor(vet_Multiplicacao, N, 0);
	inicializa_vetor(vetR, N, 0);
	inicializa_vetor(vetD, N, 0);
	inicializa_vetor(vetQ, N, 0);
	inicializa_vetor(vet_Aux, N, 0);
	inicializa_vetor(vet_Soma, N, 0);
	
	//Matriz Vetor
    Mutiplicacao_Matriz_Vetor(colptr, linhas, valores, vetX, vet_Multiplicacao, N, nao_zeros, id, np);	

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(vet_Multiplicacao, vet_Soma, N, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	Atibuicao_Vetor(vet_Soma, vet_Multiplicacao, N);
	MPI_Bcast(vet_Multiplicacao,N,MPI_FLOAT,0,MPI_COMM_WORLD);

	Subtracao_Vetor_Vetor(vetB, vet_Multiplicacao, vetR, N);
	Atibuicao_Vetor(vetR, vetD, N);

	sigma_novo = Multiplicacao_Vetor_Vetor_Soma(sigma_novo, vetR, vetR, N);
	sigma = sigma_novo;

	while(iteracao < imax && sigma_novo > (erro*erro) * sigma){
		inicializa_vetor(vetQ, N, 0);

		Mutiplicacao_Matriz_Vetor(colptr, linhas, valores, vetD, vetQ, N, nao_zeros, id, np);
		MPI_Barrier(MPI_COMM_WORLD);
		inicializa_vetor(vet_Soma, N, 0);
		MPI_Reduce(vetQ, vet_Soma, N, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		Atibuicao_Vetor(vet_Soma, vetQ, N);
		MPI_Bcast(vetQ,N,MPI_FLOAT,0,MPI_COMM_WORLD);

		alpha_linha = 0;
		alpha_linha = Multiplicacao_Vetor_Vetor_Soma(alpha_linha, vetD, vetQ, N);	
		
		alpha = 0;
		alpha = sigma_novo/alpha_linha;	
		
		for(i=0; i<N; i++){
			vetX[i] = vetX[i] + (alpha * vetD[i]);
		}
		
		if(iteracao % 50 == 0){
			inicializa_vetor(vet_Aux, N, 0);

			Mutiplicacao_Matriz_Vetor(colptr, linhas, valores, vetX, vet_Aux, N, nao_zeros, id, np);	
			MPI_Barrier(MPI_COMM_WORLD);
			inicializa_vetor(vet_Soma, N, 0);
			MPI_Reduce(vet_Aux, vet_Soma, N, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
			Atibuicao_Vetor(vet_Soma, vet_Aux, N);
			MPI_Bcast(vet_Aux,N,MPI_FLOAT,0,MPI_COMM_WORLD);

			Subtracao_Vetor_Vetor(vetB, vet_Aux, vetR, N);
		}
		else{
			for(i=0; i<N; i++){
				vetR[i] = vetR[i] - (alpha * vetQ[i]);
			}
		}

		sigma_velho = sigma_novo;
		sigma_linha = 0;	

		sigma_linha = Multiplicacao_Vetor_Vetor_Soma(sigma_linha, vetR, vetR, N);	
		sigma_novo = sigma_linha;
		
		beta = sigma_novo/sigma_velho;
		
		for(i=0; i<N; i++){
			vetD[i] = vetR[i] + (beta * vetD[i]);
		}
		
		iteracao++;
	}
	
	free(vet_Multiplicacao);
	free(vetR);
	free(vetD);
	free(vetQ);
	free(vet_Aux);
	free(vet_Soma);	
}

/******************************************************************************************************************************************/
int main(int argc, char **argv){
	char nome_arq[20] = "bcsstruc2.data";
    //char nome_arq[20] = "bcsstk01.rsa";
    //char nome_arq[20] = "matrizMenor.rsa";
    FILE *arq;
	
    int pos = 0, i,id,np;	
    int tam_linha;  // quantidade de linhas que tem o vetor linha
    int tam_coluna; // quantidade de linhas que tem o vetor coluna
    int tam_valor;  // quantidade de linhas que tem o vetor valores
    int N;          // tamanho da matriz
    int nao_zeros;  // valores nao zero
	
    int quant_vetor_coluna; // quantidade de itens que possue o vetor coluna
    	
    int *colptr = NULL;    // vetor de colunas 
    int *linhas = NULL;    // vetor de linha 
    float *valores = NULL; // vetor de valores


	MPI_Init(&argc, &argv);//bloqueia e so sai quando todos chegarem nesse ponto
	MPI_Comm_rank(MPI_COMM_WORLD, &id);//MPI_COMM_WORLD eh o grupo principal
	MPI_Comm_size(MPI_COMM_WORLD, &np);
             
	abre_arquivo(&arq,nome_arq);
	leitura_cabecalho(&arq, &tam_linha, &tam_coluna, &tam_valor, &N, &nao_zeros);
   
    colptr = (int*) malloc(tam_coluna * 16 * sizeof(int));
    linhas = (int*) malloc(nao_zeros * sizeof(int));
    valores = (float*) malloc(nao_zeros * sizeof(float));
    
    for(i=0; i<tam_coluna; i++)
  		pos = leitura_vetor(&arq, colptr, pos);
	quant_vetor_coluna = pos;
 	 	
 	pos = 0;
    for(i=0; i<tam_linha; i++)
    	pos = leitura_vetor(&arq, linhas, pos);
	
	pos = 0;	
    for(i=0; i<tam_valor; i++)
    	pos = leitura_vetor_valor(&arq, valores, pos);

	
	float *vetB = NULL;
	float *vetX = NULL;
	
	vetB = (float *)malloc(N * sizeof(float));
	vetX = (float *)malloc(N * sizeof(float));
	
	inicializa_vetor(vetB, N, 1);
	inicializa_vetor(vetX, N, 1);
	
	gradiente(colptr, linhas, valores, vetB, vetX, N, nao_zeros, quant_vetor_coluna, id, np);
	
	if(id == 0){
		escreve_vetor(tam_coluna, tam_linha, tam_valor, N, nao_zeros, quant_vetor_coluna, colptr, linhas, valores, vetX);
	}

    free(colptr);
    free(linhas);
    free(valores);
	
    free(vetB);
    free(vetX);
    MPI_Finalize();
    fecha_arquivo(&arq);
}