float* multiplica_matriz_vetor_esparso(struct pointer pointerParam, float *vetX,int id, int nproc){
	float vlNodo = 0.0, *vetR = NULL, *vetor = NULL;
	float *val = NULL;
	int *col = NULL; int *ptr = NULL;
	int i,j,x,y, resto = 0, restoAux = 0, qtde_elementos, ini_linha,fim_linha, N, auxiliar;

	/*sincroniza*/
	MPI_Barrier(MPI_COMM_WORLD);

	//seta os vetores
	val = (float*) pointerParam.vet_value;
	ptr = (int*)   pointerParam.vet_row;
	col = (int*)   pointerParam.vet_col;

	N = pointerParam.tot_linhas;

	//inicializa vetor r
	if (id == 0)
		vetR = (float *)malloc(N * sizeof(float));
	else
		vetR = (float *)malloc(nproc / N * sizeof(float));

	qtde_elementos = N/nproc;
	ini_linha = qtde_elementos * id;

	if (id==0){
		resto = N % nproc;
		qtde_elementos += resto;
	}
	else
		ini_linha += N % nproc;

	x = 0;
	for(i = ini_linha; i < ini_linha + qtde_elementos; i++){
		vlNodo = 0.0;

		for(j = ptr[i] - 1; j < ptr[i+1] - 1; j++){
			vlNodo += val[j] * vetX[col[j] - 1];
		}

		vetR[x] = vlNodo;
		x++;
	}

	MPI_Gather(&vetR[resto], N/nproc, MPI_FLOAT, &vetR[resto], N/nproc, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(vetR,N,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	//essa logica esta funcionando
	/*if (id == 0 && nproc == 1){
		for(i = ini_linha; i < N - 1; i++){
			vlNodo = vetR[N - i - 1];

			for(j = ptr[i]; j < ptr[i+1] - 1; j++){
				vlNodo += val[j] * vetX[N - col[j]];
			}
			vetR[N - i - 1] = vlNodo;
		}
	}
	else{*/
		qtde_elementos = (N - 1) /nproc;
		restoAux	   = (N - 1) % nproc;
		auxiliar	   = qtde_elementos * id;

		if (id == 0)
			qtde_elementos += restoAux;
		else
			auxiliar += restoAux;

		ini_linha = (N - 1 - (qtde_elementos * id));
		ini_linha--;

		if (id > 0){
			ini_linha -= restoAux;
		}
		fim_linha = ini_linha - (qtde_elementos - 1);
		x = auxiliar + 1;
		y = 0;

		//printf("id - %d ini_linha %d fim_linha %d x - %d\n",id,ini_linha,fim_linha,x);

		for(i = ini_linha; i >= fim_linha;i--){
			vlNodo = vetR[x];

			for(j = ptr[i]; j < ptr[i+1] - 1; j++){
				vlNodo += val[j] * vetX[N - col[j]];
			}

			if (id == 0)
				vetR[x] = vlNodo;
			else{
				vetR[y] = vlNodo;
				y++;
			}

			x++;
		}
	/*}*/

	MPI_Barrier(MPI_COMM_WORLD);

	if (nproc > 1){
		MPI_Gather(&vetR[resto], N/nproc, MPI_FLOAT, &vetR[resto], N/nproc, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(vetR,N,MPI_FLOAT,0,MPI_COMM_WORLD);
	}

	return vetR;
}
