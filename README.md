 matriz A é uma matriz de banda, simétrica, definida positiva, e altamente esparsa que deve ser gerada utilizando a função definida abaixo.

O vetor de termos independentes b é dado por

```bash
$ b[i] = f( i * π / n ) para todo i = 0:n-1,
```
onde

```bash
$ f(x) = 4π² ( sin(2πx) + sin(2π(π-x)) )
```

A norma euclidiana do resíduo (||r||) deve ser calculada a cada iteração;

O critério de parada deve ser o número máximo de iterações OU quando o erro aproximado considerando a norma do resíduo for menor do que uma tolerância especificada (o que vier primeiro);

#Execução do Programa

O pacote de software a ser construído deve gerar um executável chamado cgSolver , que deve ser invocado da seguinte forma:

``` bash
$ cgSolver n k -i <maxIter> -t <tolerancia> -o <arquivo_saida> 
```
Onde:

n: (n>0) parâmetro obrigatório definindo a dimensão do Sistema Linear.
k: parâmetro obrigatório definindo o número de bandas da matriz A. k=0 indica uma matriz diagonal. k=1 é uma matriz tridiagonal. k=2 é uma matriz pentadiagonal, e assim por diante.
-i maxIter: parâmetro opcional definindo o número máximo de iterações a serem executadas. Caso não seja definido, utilizar o valor n.
-t tolerancia: parâmetro opcional definindo o erro aproximado absoluto máximo, considerando a norma Euclidiana do resíduo.
-o arquivo_saida: parâmetro obrigatório no qual arquivo_saida é o caminho completo para o arquivo que vai conter a solução. Esta solução deve estar formatada de acordo com as regras abaixo:
No início do arquivo, deve constar sob a forma de comentários (precedido de #):
O tempo de execução de cada iteração do Método de CG;
O tempo para o cálculo do resíduo;
O valor da norma Euclidiana do resíduo a cada iteração;
O erro aproximado absoluto a cada iteração.
Em seguida o tamanho do vetor solução (n) deve ser escrito em uma linha, e na linha seguinte os valores da solução (x), separados por espaços, com printf("%.14g").
 

```bash
$ ###########
$ # Tempo Método CG: <min> <media> <max>
$ # Tempo Resíduo: <min> <media> <max>
$ #
$ # Norma Euclidiana do Residuo e Erro aproximado
$ # i=1: <norma> <erro>
$ # i=2: <norma> <erro>
$ # i=3: <norma> <erro>
$ # ...
$ ###########
$ N
$ x_1 x_2 ... x_n
```
Os Tempos devem ser calculados em mili-segundos, utilizando-se a função especificada aqui. Onde <min> <media> e <max> são, respectivamente, o menor tempo, a média de todas iterações e o maior tempo.
Tempo Método: tempo calculado a partir do início da iteração do método até a obtenção do vetor solução daquela iteração.
Tempo Resíduo: Tempo para calcular o resíduo.
 

#MENSAGENS DE ERRO:

Em caso de erros, uma mensagem explicando o ocorrido deve ser impressa em stderr e a execução do programa deve ser encerrada com código diferente de 0.

Produto a ser Entregue
O trabalho deve ser desenvolvido por um grupo composto por no máximo DOIS alunos regularmente matriculados na disciplina. O grupo NÃO PODE SER ALTERADO na próxima parte do trabalho.

Cada grupo deve entregar um pacote de software completo contendo os fontes em linguagem C. O pacote deve ser arquivado e compactado com tar(1) e gzip(1) em um arquivo chamado login1.tar.gz (se grupo com 1 membro) ou login1-login2.tar.gz (se grupo com 2 membros), onde login1 e login2 são os logins dos alunos que compõem o grupo.

O pacote deve ter a seguinte estrutura de diretório e arquivos:

./login1-login2/: diretório principal;
./login1-login2/LEIAME;
./login1-login2/Makefile;
Note que a extração dos arquivos de login1-login2.tar.gz deve criar o diretório login1-login2 contendo todos os arquivo acima. Os arquivos fonte também devem estar contidos no diretório, ou em algum sub-diretório, desde que o Makefile funcione.


#LEIAME
O pacote deve conter um arquivo de documentação em texto plano (ASCII). Este arquivo, chamado LEIAME, deve conter as seguintes informações:

autoria do software, isto é, nome e RA dos membros do grupo;
lista dos arquivos e diretórios contidos no pacote e sua descrição (breve);
um capítulo especial descrevendo os algoritmos e as estruturas de dados utilizadas, as alternativas de implementação consideradas e/ou experimentadas e os motivos que o levaram a optar pela versão entregue, as dificuldades encontradas e as maneiras pelas quais foram contornadas.
bugs conhecidos;

#Makefile
O arquivo Makefile deve possuir as regras necessárias para compilar os módulos individualmente e gerar o programa executável. As seguintes regras devem existir OBRIGATORIAMENTE:

all: compila e produz um executável chamado cgSolver no diretório login1-login2/;
clean: remove todos os arquivos temporários e os arquivos gerados pelo Makefile (*.o, executável, etc.).
 

#Entrega
O prazo final para a entrega deste trabalho é dia 25 de setembro de 2016, 23:59:59h, IMPRETERIVELMENTE.

O trabalho deve ser enviado como anexo ao email do professor com o Assunto (Subject): CI164 - Trabalho 1.
No corpo da mensagem DEVE CONSTAR OBRIGATORIAMENTE o Nome e Números de Registro Acadêmico (RA) dos membros do grupo;
O grupo deverá considerar o trabalho como entregue SOMENTE APÓS RECEBER DO PROFESSOR UMA MENSAGEM DE CONFIRMAÇÃO DE RECEBIMENTO dentro de 24 horas desde o envio do trabalho;
 

#Critérios de Avaliação
APENAS OS TRABALHOS QUE FUNCIONAREM SERÃO CORRIGIDOS. Se o trabalho não compilar ou acusar falha de segmentação (Segmentation fault) prematura durante os testes realizados pelo professor (sem que qualquer operação se efetue a contento), trará para o grupo NOTA 0 (ZERO). Também receberão NOTA 0 (ZERO) trabalhos plagiados de qualquer fonte, e/ou com códigos similares. Além disso, apenas trabalhos entregues no prazo marcado receberão nota.

Os itens de avaliação do trabalho e respectivas pontuações são:

Qualidade da documentação: arquivo LEIAME (10 pontos)
Entrada e saída: funcionamento do programa de acordo com a especificação no que tange execução, entrada e saída (10 pontos)
Funcionamento: corretude das respostas nos testes executados (60 pontos)
Eficiência das estruturas de dados utilizadas, desde que devidamente justificadas no arquivo LEIAME (20 pontos)
Defesa: A defesa do trabalho será oral, e definirá a nota individual de cada membro da equipe, de acordo com seu conhecimento a respeito do trabalho.

 

#Função para Medição de Tempo
O tempo de execução deve ser medido em mili-segundos, considerando tempo de relógio, utilizando a função especificada abaixo:

```bash
$ #include <time.h>
$
$ double timestamp(void){
$    struct timeval tp;
$    gettimeofday(&tp, NULL);
$    return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
$ }
```
onde o tempo decorrido é medido pela diferença do "timestamp" medidos antes e depois da região de interesse.

 

#Função para geração de sistemas lineares
A matriz A do SL a ser resolvido deve ser gerada utilizando a função abaixo.

ATENÇÃO: A matriz deve ser simétrica, isto é, os valores nas diagonais superiores e inferiores devem ser idênticos.

``` bash
$ #include <stdlib.h>

$ /***********************
$  * N: tamanho do sistema linear
$ * k: numero da diagonal, 0 = diagonal principal, 1 = acima/abaixo da diagonal, 2 = ...
$ * kMax: numero de bandas do sistema linear
$ * diag: vetor para armazenar os valores da diagonal. Deve ser alocado por quem chama a função.
$ ***********************/
$ int generateRandomDiagonal( unsigned int N, unsigned int k, unsigned int kMax, double *diag )
$ {
$  if ( !diag || N < 3 || kMax > N/2 || k < 0 || k > kMax )
$    return (-1);
$
$  /* garante valor dominante para diagonal principal */
$  double fator = (k == 0) ? ((double)(kMax-1)) : (0.0);
$
$  double invRandMax = 1.0 / (double)RAND_MAX;
$
$  for (int i=0; i < N-k; ++i)
$  {
$    diag[i] = fator + (double)rand() * invRandMax;
$  }
$
$  return (0);
}

```
Antes da primeira chamada da função "generateRandomDiagonal()", e somente uma vez em todo código, você deve inicializar a sequência de números aleatóreos chamando a função:
```bash
$ srand( 20162 );
```
"# GradienteConjugado" 
