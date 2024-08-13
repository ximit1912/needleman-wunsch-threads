#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <pthread.h>

#define A 0 // representa uma base Adenina
#define T 1 // representa uma base Timina
#define G 2 // representa uma base Guanina
#define C 3 // representa uma base Citosina
#define X 4 // representa um gap

#define sair 11

#define maxSeq 10000 // tamanho maximo de bases em uma sequencia genomica
#define maxThreads 1000 // tamanho maximo de threads (k)
/////////////////////////////////////////////////////////////////////////////



/***** VETORES GLOBAIS *****/
/* mapaBases mapeia indices em caracteres que representam as bases, sendo 0='A',
1='T', 2='G', 3='C' e 4='-' representando gap */
char mapaBases[5]={'A','T','G','C','-'};


/* seqMaior e seqMenor representam as duas sequencias de bases de entrada, a
   serem comparadas, inicializadas conforme segue. Elas conterao os indices aos
   inves dos proprios caracteres. seqMenor deve ser menor ou igual a seqMaior. */
int  seqMaior[maxSeq]={A,A,C,T,T,A},
     seqMenor[maxSeq]={A,C,T,T,G,A};


/* alinhaGMaior representava a sequencia maior ja alinhada, assim como alinhaGMenor,
   ambas obtidas no traceback. As duas juntas, pareadas, formam o alinhamento
   global. Tal alinhamento global pode ser obtido de duas formas: a partir do
   primeiro maior escore ou a partir do ultimo maior escore 

    int  alinhaGMaior[maxSeq],
         alinhaGMenor[maxSeq];  
*/
/////////////////////////////////////////////////////////////////////////////




/***** VARIÁVEIS GLOBAIS *****/
int tamSeqMaior=6,  /* tamanho da sequencia maior, inicializado como 6 */
    tamSeqMenor=6,  /* tamanho da sequencia menor, inicializado como 6 */
    penalGap=0,     /* penalidade de gap, a ser descontada no escore acumulado
                       quando um gap eh encontrado */
    grauMuta=0,     /* porcentagem maxima de mutacao na geracao aleatoria da
                       sequencia menor, a qual eh copiada da maior e sofre algumas
                       trocas de bases */
    k, AlinhamentosParciais, /* quantidade de threads, deve ser menor ou igual a tamSeqMenor, 
                                pois k threads irão preencher as linhas da matriz de escore de
                                forma igualmente distribuídas, além de que farão o alinhamento/traceback
                                da seqMenor na maior */
    tamAlinha[maxThreads];      /* tamanho do alinhamento global obtido */
/////////////////////////////////////////////////////////////////////////////




/***** VARIÁVEIS GLOBAIS AUXILIARES *****/
int indRef=-1,  // indice da sequencia maior a partir do qual extrai a sequencia menor, no caso de geracao aleatoria
    nTrocas=-1, // quantidade de trocas na geracao automatica da sequencia menor, a partir de um segmento da sequencia maior
    linPMaior, colPMaior, PMaior, // suporte para deteccao do primeiro maior escore
    linUMaior, colUMaior, UMaior; // suporte para deteccao do ultimo maior escore
/////////////////////////////////////////////////////////////////////////////




/***** STRUCT's PARAMETRO DAS THREADS *****/
typedef struct{
  int id;
  int *aux;
}tPack;

typedef struct{
  int id, tbLin, tbCol, pos, tamParcial;
  int (*auxMatriz)[10000]; /* apontará para a matriz de alinhamentos que será criada na função traceback */
}tPackAlinhamento;
 ////////////////////////////////////////////////////////////////////////////




/***** MATRIZ DE ESCORES *****/
int matrizEscores[maxSeq+1][maxSeq+1];
/////////////////////////////////////////////////////////////////////////////

int alinhaParciais[maxThreads*2][maxSeq]; /* matriz contendo o par do alinhamento globlal
                                         mais ateh k-1 pares de alinhamentos parciais  */


/***** MATRIZ DE PESOS *****/
/*  
    matrizPesos contem os pesos do pareamento de bases. Estruturada e inicializada
    conforme segue, onde cada linha ou coluna se refere a uma das bases A, T, G
    ou C. Considera-se a primeira dimensao da matriz como linhas e a segunda como
    colunas. Na configuracao default, o peso de bases iguais eh 1 e o peso de bases
    diferentes eh 0, mas pode ser alterado. Outra configuracao usual eh 2 para
    bases iguais e -1 para bases diferentes 

                                        0 1 2 3
                                        A T G C
                                    0 A 1 0 0 0
                                    1 T 0 1 0 0
                                    2 G 0 0 1 0
                                    3 C 0 0 0 1 
*/
int matrizPesos[4][4]={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
/////////////////////////////////////////////////////////////////////////////




/******************** INICIO DAS FUNCOES *********************/
/* leitura do tamanho da sequencia maior */
int leTamMaior(void)
{
  bool aux_err=false;

  printf("\nLeitura do Tamanho da Sequencia Maior:");
  do
  { if(aux_err==true)
        printf("\n**ERRO** Voce digitou um valor incorreto!");
    printf("\nDigite 0 < valor < %d = ", maxSeq);
    scanf("%d", &tamSeqMaior);

    aux_err=true;
  } while ((tamSeqMaior<1)||(tamSeqMaior>maxSeq));
}

/* leitura do tamanho da sequencia menor */
int leTamMenor(void)
{
  bool aux_err=false;

  printf("\nLeitura do Tamanho da Sequencia Menor:");
  do
  { if(aux_err==true)
        printf("\n**ERRO** Voce digitou um valor incorreto!");
    printf("\nDigite 0 < valor <= %d = ", tamSeqMaior);
    scanf("%d", &tamSeqMenor);

    aux_err=true;
  } while ((tamSeqMenor<1)||(tamSeqMenor>tamSeqMaior));
}


/* leitura do valor da penalidade de gap */
int lePenalidade(void)
{   int penal;

  printf("\nLeitura da Penalidade de Gap:");
  do
  {
    printf("\nDigite valor >= 0 = ");
    scanf("%d", &penal);
  } while (penal<0);

  return penal;
}


/* leitura da matriz de pesos */
void leMatrizPesos()
{ int i,j;

  printf("\nLeitura da Matriz de Pesos:\n");
  for (i=0; i<4; i++)
  {
    for (j=0; j<4; j++)
    {
      printf("Digite valor %c x %c = ", mapaBases[i], mapaBases[j]);
      scanf("%d",&(matrizPesos[i][j]));
    }
    printf("\n");
  }
}

/* mostra da matriz de pesos */
void mostraMatrizPesos(void)
{ int i,j;

  printf("\nMatriz de Pesos Atual:");
  printf("\n%4c%4c%4c%4c%4c\n",' ','A','T','G','C');
  for (i=0; i<4; i++)
  {
    printf("%4c",mapaBases[i]);
    for (j=0; j<4; j++)
      printf("%4d",matrizPesos[i][j]);
    printf("\n");
  }
}


/* leitura da porcentagem maxima (grau) de mutacao aleatoria */
int leGrauMutacao(void)
{ int prob;

  printf("\nLeitura da Porcentagem Maxima de Mutacao Aleatoria:\n");
  do
  { printf("\nDigite 0 <= valor <= 100 = ");
    scanf("%d", &prob);
  } while ((prob<0)||(prob>100));

  return prob;
}
/////////////////////////////////////////////////////////////////////////////




/***** SEQUENCIAS *****/
/* leitura manual das sequencias de entrada seqMaior e seqMenor */
void leSequencias(void)
{ int i, erro;
  char seqMaiorAux[maxSeq], seqMenorAux[maxSeq];

  indRef=-1;
  nTrocas=-1;
  printf("\nLeitura das Sequencias:\n");

  /* lendo a sequencia maior */
  do
  {
    printf("\nPara a Sequencia Maior,");
    printf("\nDigite apenas caracteres 'A', 'T', 'G' e 'C'");
    do
    { printf("\n> ");
      fgets(seqMaiorAux,maxSeq,stdin);
      tamSeqMaior=strlen(seqMaiorAux)-1; /* remove o enter */
    } while (tamSeqMaior<1);
    printf("\ntamSeqMaior = %d\n",tamSeqMaior);
    i=0;
    erro=0;
    do
    {
      switch (seqMaiorAux[i])
      {
        case 'A': seqMaior[i]=A;
                  break;
        case 'T': seqMaior[i]=T;
                  break;
        case 'G': seqMaior[i]=G;
                  break;
        case 'C': seqMaior[i]=C;
                  break;
        default: erro=1;  /* nao eh permitido qquer outro caractere */
      }
      i++;
    } while ((erro==0)&&(i<tamSeqMaior));
  }while (erro==1);

  /* lendo a sequencia menor */
  do
  {
    printf("\nPara a Sequencia Menor, ");
    printf("\nDigite apenas caracteres 'A', 'T', 'G' e 'C'");
    do
    { printf("\n> ");
      fgets(seqMenorAux,maxSeq,stdin);
      tamSeqMenor=strlen(seqMenorAux)-1; /* remove o enter */
    } while ((tamSeqMenor<1)||(tamSeqMenor>tamSeqMaior));
    printf("\ntamSeqMenor = %d\n",tamSeqMenor);

    i=0;
    erro=0;
    do
    {
      switch (seqMenorAux[i])
      {
        case 'A': seqMenor[i]=A;
                  break;
        case 'T': seqMenor[i]=T;
                  break;
        case 'G': seqMenor[i]=G;
                  break;
        case 'C': seqMenor[i]=C;
                  break;
        default: erro=1;
      }
      i++;
    } while ((erro==0)&&(i<tamSeqMenor));
  }while (erro==1);
}

/* leitura por arquivo das sequencias de entrada seqMaior e seqMenor */
void leSequenciasArquivo()
{  int i, erro1, erro2;
   char seqMaiorAux[maxSeq], seqMenorAux[maxSeq];
   FILE *arq;

   indRef=-1;
   nTrocas=-1;
   printf("\nLeitura das Sequencias no arquivo, DIRETRIZES:\n - o nome do arquivo deve ser sequencias.txt\n - as sequencias deverao conter apenas os caracteres maiusculos A,T,G,C\n - a sequencia maior deve estar na primeira linha e a menor na 2 linha\n");

   /* abre o arquivo para leitura */ 
   arq=fopen("sequencias.txt", "r");
   if (arq==NULL) /* verifica se foi aberto corretamente */
   {
      fprintf(stderr,"\nArquivo sequencias.txt nao existe!");  
      exit(EXIT_FAILURE);      
   }
   else
   {
      do
       { 
         /* verifica se há algum erro no modo como as sequencias devem estar escritas no arquivo  */
         if(erro1==1 || erro2==1)
         { printf("\n*********************************************************************");
           printf("\nVoce inseriu as sequencias de forma errada no arquivo sequencias.txt!");
           printf("\nVerifique se ha alguma incongruencia, como: "); 
           printf("\n - um espaco antes de quebra de linha\n - algum caracter em minusculo");
           printf("\nApos verificar e corrijir, salve o arquivo novamente com o mesmo nome...");
           printf("\n*********************************************************************\n");
      
           fclose(arq);                            /*  fecha o arquivo  */
           system("pause");                        /*  pausa o código   */
           
           arq=fopen("sequencias.txt", "r");       /* abre novamente */
           if (arq==NULL) 
           {
               fprintf(stderr,"\n\nArquivo sequencias.txt nao existe!");  
               exit(EXIT_FAILURE);      
           } 
        }                                         
      

         /* lê as 2 sequencias no arquivo e atribui para as
            variáveis auxiliares seqMaiorAux e seqMenorAux */ 
         fscanf(arq, "%s", seqMaiorAux);
         fscanf(arq, "%s", seqMenorAux);

         /* atribui os seus tamanhos às variáveis referentes */
         tamSeqMaior=strlen(seqMaiorAux);
         tamSeqMenor=strlen(seqMenorAux);

         printf("\nTamanho da sequencia maior = %i", tamSeqMaior);
         printf("\nTamanho da sequencia menor = %i", tamSeqMenor);

         i=0;
         erro1=0;
         do
         {
           switch (seqMaiorAux[i])
           {
             case 'A': seqMaior[i]=A;   
                       break;
             case 'T': seqMaior[i]=T;
                       break;
             case 'G': seqMaior[i]=G;
                       break;
             case 'C': seqMaior[i]=C;
                       break;
             default:  erro1=1;  /* nao eh permitido qualquer outro caractere */
            }
            i++;
         } while(erro1==0 && i<tamSeqMaior);

         i=0; 
         erro2=0; 
         do
         {
           switch (seqMenorAux[i])
           {
             case 'A': seqMenor[i]=A;
                       break;
             case 'T': seqMenor[i]=T;
                       break;
             case 'G': seqMenor[i]=G;
                       break;
             case 'C': seqMenor[i]=C;
                       break;
             default:  erro2=1;   /* nao eh permitido qualquer outro caractere */
            }
            i++;  
         } while(erro2==0 && i<tamSeqMenor);
     } while(erro1==1 || erro2==1);
    
     

     if (fclose(arq))   
        fprintf(stderr,"\n\nArquivo sequencias.txt nao foi fechado corretamente!");
     else 
        printf("\n\nArquivo sequencias.txt foi fechado corretamente!");
   }
}

/* geracao das sequencias aleatorias, conforme tamanho e grau de mutação */
void geraSequencias(void)
{   int i, dif, probAux;
    char base;

    printf("\nGeracao Aleatoria das Sequencias:\n");

    /* gerando a sequencia maior */
    for (i=0; i<tamSeqMaior; i++)
      {
        base=rand()%4; /* produz valores de 0 a 3 */
        seqMaior[i]= base;
      }

    dif=tamSeqMaior-tamSeqMenor; /* diferenca entre os tamanhos das sequencias */

    indRef=0;
    if (dif>0)
      indRef=rand()%dif; /* produz um indice aleatorio para indexar a sequencia maior,
                         para a partir dele extrair a primeira versao da sequencia
                         menor */

    /* gerando a sequencia menor a partir da maior. Copia trecho da sequencia
       maior, a partir de um indice aleatorio que nao ultrapasse os limites do
       vetor maior */
    for (i=0; i<tamSeqMenor; i++)
        seqMenor[i]=seqMaior[indRef+i];

    /* causa mutacoes aleatorias na sequencia menor para gerar "gaps",
       sobre cada base, de acordo com o grau (porcentagem) informado.
       A mutacao causa a troca da base original por outra base aleatoria
       necessariamente diferente. Gera-se uma probabilidade aleatoria
       ateh 100 e se ela estiver dentro do grau informado, a mutacao
       ocorre na base, caso contrario, mantem a base copiada. */

    i=0;
    nTrocas=0;
    while ((i<tamSeqMenor)&&(nTrocas<((grauMuta*tamSeqMenor)/100)))
    {
      probAux=rand()%100+1;

      if (probAux<=grauMuta)
      {
        seqMenor[i]=(seqMenor[i]+(rand()%3)+1)%4;
        nTrocas++;
      }
      i++;
    }

    printf("\nSequencias Geradas: Dif = %d, IndRef = %d, NTrocas = %d\n ",dif, indRef, nTrocas);
}


/* mostra das sequencias seqMaior e seqMenor */
void mostraSequencias(void)
{   int i, nTrocas;

  /* imprime a sequência maior utilizando o vetor de mapeamento */
  printf("\nSequencias Atuais:\n");
  printf("\nSequencia Maior, Tam = %d\n", tamSeqMaior);
  for (i=0; i<tamSeqMaior; i++)
    printf("%c",mapaBases[seqMaior[i]]);
  printf("\n");

  /* trata o índice de referencia (indRef) */
  if(indRef == -1)
  { printf("Nao existe Indice de Referencia ('^')= %d\n", indRef);
    indRef = 0;
  }else
      { /* imprime embaixo aonde está o índice de referência calculado a partir da diferença de tamanho das sequências */
        for (i=0; i<tamSeqMaior; i++)
          if (i!=indRef)
            printf(" ");
          else 
            printf("^");

        printf("\nIndice de Referencia ('^')= %d\n", indRef);
      }
  /* imprime a sequência menor utilizando o vetor de mapeamento */
  printf("\nSequencia Menor, Tam = %d\n", tamSeqMenor);
  for (i=0; i<tamSeqMenor; i++)
    printf("%c",mapaBases[seqMenor[i]]);
  printf("\n");

  

  nTrocas=0;


  /* imprime embaixo de cada base da sequência menor se houve correspondência (+) ou se houve mutação (-) */
  for (i=0; i<tamSeqMenor; i++)
      if (seqMenor[i]!=seqMaior[indRef+i])
      {
           printf("-");
           nTrocas++;
      }
      else 
           printf("+");
  printf("\nQuantidade de trocas ('-') = %d\n", nTrocas);

}
/////////////////////////////////////////////////////////////////////////////

 

/***** MATRIZ DE ESCORES *****/
/* grava a matriz de escores bem formatada em um arquivo.txt */
void gravarMatrizArquivo()
{
  int i, lin, col;
  FILE *arq;

  time_t mytime;
  mytime = time(NULL);
  struct tm tm = *localtime(&mytime);

  arq=fopen("matrizes.txt", "a");
   if (arq==NULL) /* verifica se foi aberto corretamente */
   {
      fprintf(stderr,"\nFalha na gravação da matriz no arquivo!");  
      exit(EXIT_FAILURE);      
   }
   else
   {
    fprintf(arq, "*****Matriz de escores Atual (%d threads)*****\n", k);
    fprintf(arq, "***** [Gerada em %d/%d/%d -> %d:%d:%d] *****\n\n", tm.tm_mday, tm.tm_mon + 1, tm.tm_year + 1900, tm.tm_hour, tm.tm_min, tm.tm_sec);
    fprintf(arq, "%4c%4c",' ',' ');
    for (i=0; i<=tamSeqMaior; i++)
      fprintf(arq, "%4d",i);
    fprintf(arq, "\n");

    fprintf(arq, "%4c%4c%4c",' ',' ','-');
    for (i=0; i<tamSeqMaior; i++)
      fprintf(arq, "%4c",mapaBases[(seqMaior[i])]);
    fprintf(arq, "\n");

    fprintf(arq, "%4c%4c",'0','-');
    for (col=0; col<=tamSeqMaior; col++)
      fprintf(arq, "%4d",matrizEscores[0][col]);
    fprintf(arq, "\n");

    for (lin=1;lin<=tamSeqMenor;lin++)
    {
      fprintf(arq, "%4d%4c",lin,mapaBases[(seqMenor[lin-1])]);
      for (col=0;col<=tamSeqMaior;col++)
      {
        fprintf(arq, "%4d",matrizEscores[lin][col]);
      }
      fprintf(arq, "\n");
    }

    nTrocas=0;
    for (i=0; i<tamSeqMenor; i++)
      if (seqMenor[i]!=seqMaior[indRef+i])
        nTrocas++;

    // grava abaixo a qnt de trocas, o primeiro maior escore e o último maior escore
    fprintf(arq, "\nQuantidade de trocas  = %d", nTrocas);
    fprintf(arq, "\nPrimeiro Maior escore = %d na celula [%d,%d]", PMaior, linPMaior, colPMaior);
    fprintf(arq, "\nUltimo Maior escore   = %d na celula [%d,%d]\n", UMaior, linUMaior, colUMaior);

    fprintf(arq, "__________________________________________________________________\n");
    
    fclose(arq);
   }
}


/* função que cada thread fará */
void *geraComThreads(void *ptr)
{ int lin, col, peso,
      escoreDiagAnterior,      /* escore da diagonal anterior da matriz de escores */
      escoreLinAnterior,       /* escore da linha anterior da matriz de escores */
      escoreColAnterior;       /* escore da coluna anterior da matriz de escores */

  tPack *t = (void *) ptr;


  /* calculando os demais escores, percorrendo todas as posicoes
     da matriz, linha por linha, coluna por coluna, aplicando
     a seguinte formula:
                             /  f(lin-1,col-1)+matrizPesos[lin,col]
     f(lin,col) = m�ximo de {  f(lin,col-1)-penalGap
                             \  f(lin-1,col)-penalGap               
  */
  for (lin=t->id; lin<=tamSeqMenor;)
  // thread preencherá de k em k linhas, começando pela linha referente ao seu id
  {
    // printf("\n***Preenchimento das colunas (c) da linha %d(i) pela thread %d(t) iniciada!***\n", lin, t->id);
    /* 
         REMOVER COMO COMENTÁRIO PARA IMPRIMIR NA TELA QUAL ELEMENTO
         DA MATRIZ DE ESCORES ESTÁ SENDO PREENCHIDO 
      */
    for (col=1; col<=tamSeqMaior; col++)
    {
      peso=matrizPesos[(seqMenor[lin-1])][(seqMaior[col-1])];
      escoreDiagAnterior = matrizEscores[lin-1][col-1]+peso;
      escoreLinAnterior = matrizEscores[lin][col-1]-penalGap;
      while(t->aux[lin-1] < col); 
      /* 
         faz aguardar enquanto a thread anterior
         ainda não preencheu o elemento de cima
         necessário para a função de pontuação ser realizada
      */
      escoreColAnterior = matrizEscores[lin-1][col]-penalGap; 
      
      if ((escoreDiagAnterior>escoreLinAnterior)&&(escoreDiagAnterior>escoreColAnterior)) // trocado >= por >
      {
        matrizEscores[lin][col]=escoreDiagAnterior;
      }
      else if (escoreLinAnterior>escoreColAnterior)
           {
              matrizEscores[lin][col]=escoreLinAnterior;
           }
           else
           { 
              matrizEscores[lin][col]=escoreColAnterior;
           }
    
      // printf("(i,j)=(%d,%d)->t:%d].", col, lin, t->id);   
      /* 
         REMOVER COMO COMENTÁRIO PARA IMPRIMIR NA TELA QUAL ELEMENTO
         DA MATRIZ DE ESCORES ESTÁ SENDO PREENCHIDO 
      */
      t->aux[lin] = col; // incrementa a quantidade de colunas preenchidas pela thread
    }
      // printf("\n");

      lin = lin + k;
  }
}

/* faz o preenchimento da matriz de escores */
void geraMatrizEscores(void)
{ int i, lin, col, 
      iret[k]; // k auxiliares de retorno das threads
  
  tPack t[k];  // k structs que cada thread conterá como parâmetros
  pthread_t thread[k];  // k threads

  int vetorAuxiliar[tamSeqMenor+1];  /* vetor auxiliar com tamSeqMenor+1 posicoes, para cada linha */
  vetorAuxiliar[0] = tamSeqMaior;    /* ele conterá a qnt de colunas já prenchidas servindo para   */
                                     /* sincronizar as threads na leitura do escore da coluna anterior */

  printf("\nGeracao da Matriz de escores:\n");
  /*  A matriz sera gerada considerando que ela representa o cruzamento
      da seqMenor[] associada as linhas e a seqMaior[] associada as
      colunas. */
 
  /* inicializando a linha de penalidades/gaps */
  for (col=0; col<=tamSeqMaior; col++)
    matrizEscores[0][col]=-1*(col*penalGap);

  /* inicializando a coluna de penalidades/gaps */
  for (lin=0; lin<=tamSeqMenor; lin++)
    matrizEscores[lin][0]=-1*(lin*penalGap);


  /* inicialização do vetor auxiliar e das k threads */
  for(i=1; i<=tamSeqMenor; i++)
  {
    vetorAuxiliar[i] = 0;
  } 

  for(i=0; i<k; i++)
  {
    t[i].id = i+1;
    t[i].aux = vetorAuxiliar;
  }

  for(i=0; i<k; i++)
  {
    iret[i] = pthread_create(&thread[i], NULL, geraComThreads, (void*) &t[i]);
  }

  for(i=0; i<k; i++)
  {
    pthread_join(thread[i], NULL);
  }
  
  /* localiza o primeiro e o ultimo maior escores e suas posicoes. */
  linPMaior=1;
  colPMaior=1;
  PMaior=matrizEscores[1][1];

  linUMaior=1;
  colUMaior=1;
  UMaior=matrizEscores[1][1];

  for (lin=1; lin<=tamSeqMenor; lin++)
  {
    for (col=1; col<=tamSeqMaior; col++)
    {
      if (PMaior<matrizEscores[lin][col])
      {
        linPMaior=lin;
        colPMaior=col;
        PMaior=matrizEscores[lin][col];
      }
      if (UMaior<=matrizEscores[lin][col])
      {
        linUMaior=lin;
        colUMaior=col;
        UMaior=matrizEscores[lin][col];
      }
    }
  }
  printf("\nMatriz de escores Gerada.");
  printf("\nPrimeiro Maior escore = %d na celula [%d,%d]", PMaior, linPMaior, colPMaior);
  printf("\nUltimo Maior escore = %d na celula [%d,%d]", UMaior, linUMaior, colUMaior);


  gravarMatrizArquivo();
}


/* imprime a matriz de escores */
void mostraMatrizEscores()
{ int i, lin, col;

  printf("\nMatriz de escores Atual:\n");

  printf("%4c%4c",' ',' ');
  for (i=0; i<=tamSeqMaior; i++)
    printf("%4d",i);
  printf("\n");

  printf("%4c%4c%4c",' ',' ','-');
  for (i=0; i<tamSeqMaior; i++)
    printf("%4c",mapaBases[(seqMaior[i])]);
  printf("\n");

  printf("%4c%4c",'0','-');
  for (col=0; col<=tamSeqMaior; col++)
    printf("%4d",matrizEscores[0][col]);
  printf("\n");

  for (lin=1;lin<=tamSeqMenor;lin++)
  {
    printf("%4d%4c",lin,mapaBases[(seqMenor[lin-1])]);
    for (col=0;col<=tamSeqMaior;col++)
    {
      printf("%4d",matrizEscores[lin][col]);
    }
    printf("\n");
  }
}
/////////////////////////////////////////////////////////////////////////////

void *traceBackComThreads(void *ptr)
{
  tPackAlinhamento *t = (void *) ptr;

  printf("\n*Thread %d criada pra fazer o alinhamento %d*", t->id, AlinhamentosParciais);

  int aux, i, peso,
      escoreDiagAnterior,      /* escore da diagonal anterior da matriz de escores */
      escoreLinAnterior,       /* escore da linha anterior da matriz de escores */
      escoreColAnterior;       /* escore da coluna anterior da matriz de escores */

  do
  {
    // a seguir verifica o escore do elemento [tbLin, tbCol]

    peso=matrizPesos[(seqMenor[t->tbLin-1])][(seqMaior[t->tbCol-1])];
    escoreDiagAnterior = matrizEscores[t->tbLin-1][t->tbCol-1]+peso;
    escoreLinAnterior = matrizEscores[t->tbLin][t->tbCol-1]-penalGap;
    escoreColAnterior = matrizEscores[t->tbLin-1][t->tbCol]-penalGap;

      if ((escoreDiagAnterior>escoreLinAnterior)&&(escoreDiagAnterior>escoreColAnterior)) 
      {
        if (seqMenor[t->tbLin-1]!=seqMaior[t->tbCol-1])
        {   /* O escore da diagonal venceu, mas os elementos correspondentes entre
               as sequencias menor e maior sao diferentes. Nesse caso, surge um
               gap duplo */

            printf("\nALERTA no TraceBack: Pos = %d Lin = %d e Col = %d\n", t->pos, t->tbLin, t->tbCol);

            t->auxMatriz[(t->id-1)*2][t->pos]=seqMaior[t->tbCol-1]; //alinhaGMaior[pos]=seqMaior[tbCol-1];
            t->auxMatriz[(t->id-1)*2 + 1][t->pos]=X;  // alinhaGMenor[pos]=X;
            t->tbCol--;
            t->pos++;

            t->auxMatriz[(t->id-1)*2][t->pos]=X; //alinhaGMaior[pos]=X;
            t->auxMatriz[(t->id-1)*2 + 1][t->pos]=seqMenor[t->tbLin-1];  //alinhaGMenor[pos]=seqMenor[tbLin-1];
            t->tbLin--;
            t->pos++;

        }
        else {
            t->auxMatriz[(t->id-1)*2][t->pos]=seqMaior[t->tbCol-1]; //alinhaGMaior[pos]=seqMaior[tbCol-1];
            t->tbCol--;
            t->auxMatriz[(t->id-1)*2 + 1][t->pos]=seqMenor[t->tbLin-1]; //alinhaGMenor[pos]=seqMenor[tbLin-1];
            t->tbLin--;
            t->pos++;
        }
      }
      else if (escoreLinAnterior>escoreColAnterior)
           {
              t->auxMatriz[(t->id-1)*2][t->pos]=seqMaior[t->tbCol-1]; //alinhaGMaior[pos]=seqMaior[tbCol-1];
              t->auxMatriz[(t->id-1)*2 + 1][t->pos]=X; //alinhaGMenor[pos]=X;
              t->tbCol--; 
              t->pos++; 
           }
           else if (escoreLinAnterior<escoreColAnterior)
                {
                    t->auxMatriz[(t->id-1)*2][t->pos]=X; //alinhaGMaior[pos]=X;
                    t->auxMatriz[(t->id-1)*2 + 1][t->pos]=seqMenor[t->tbLin-1]; //alinhaGMenor[pos]=seqMenor[tbLin-1];
                    t->tbLin--;
                    t->pos++;
                } 
                else if (escoreLinAnterior==escoreColAnterior)
                      /* 
                         se for igual vai continuar fazendo o alinhamento global na horizontal com a 
                         thread 1 mesmo e além disso, criará mais uma thread para continuar mais um 
                         alinhamento parcial utilizando indo na vertical.
                      */
                     {
                        if(AlinhamentosParciais<k)
                        {
                        /* 
                          verificar se  já foram as k threads disparadas, neste caso não serão feitos
                          mais alinhamentos parciais pois só serão mostrado até k alinhamentos
                        */
                          // (escoreLinAnterior > escoreColAnterior) 
                          t->auxMatriz[(t->id-1)*2][t->pos]=seqMaior[t->tbCol-1]; //alinhaGMaior[pos]=seqMaior[tbCol-1];
                          t->auxMatriz[(t->id-1)*2 + 1][t->pos]=X; //alinhaGMenor[pos]=X;
                          t->tbCol--; 
                          t->pos++; 

                           // (escoreLinAnterior < escoreColAnterior)
                          t->auxMatriz[(t->id)*2][t->pos]=X; //alinhaGMaior[pos]=X;
                          t->auxMatriz[(t->id)*2 + 1][t->pos]=seqMenor[t->tbLin-1]; //alinhaGMenor[pos]=seqMenor[tbLin-1];
                          
                          tPackAlinhamento tParcial;
                          
                          tParcial.tbLin=t->tbLin-1; 
                          tParcial.pos=t->pos+1;
                          tParcial.tbCol=t->tbCol; 
                          tParcial.id=t->id+1; // atribui o id da thread nova
                          tParcial.tamParcial=0;
                          tParcial.auxMatriz=alinhaParciais;  ////////////////////

                          AlinhamentosParciais++;

                          pthread_t threadParcial;

                          int iret = pthread_create(&threadParcial, NULL, traceBackComThreads, (void *) &tParcial);
                          pthread_join(threadParcial, NULL);
                        }
                        else
                         /* 
                            se for a thread de numero k, irá fazer o último alinhamento parcial
                            utilizando (escoreLinAnterior > escoreColAnterior) como desempate, 
                            e então, para de criar threads 
                         */
                        { 
                          t->auxMatriz[(t->id-1)*2][t->pos]=seqMaior[t->tbCol-1]; //alinhaGMaior[pos]=seqMaior[tbCol-1];
                          t->auxMatriz[(t->id-1)*2 + 1][t->pos]=X; //alinhaGMenor[pos]=X;
                          t->tbCol--; 
                          t->pos++; 
                        }
                     }
    t->tamParcial++;
  } while ((t->tbLin!=0)&&(t->tbCol!=0));


  /* descarrega o restante de gaps da linha 0, se for o caso */
  while (t->tbLin>0)
  {
    t->auxMatriz[(t->id-1)*2][t->pos]=X; //alinhaGMaior[pos]=X;
    t->auxMatriz[(t->id-1)*2 + 1][t->pos]=seqMenor[t->tbLin-1];  //alinhaGMenor[pos]=seqMenor[tbLin-1];
    t->tbLin--;
    t->pos++;
  }

  /* descarrega o restante de gaps da coluna 0, se for o caso */
  while (t->tbCol>0)
  {
    t->auxMatriz[(t->id-1)*2][t->pos]=seqMaior[t->tbCol-1]; //alinhaGMaior[pos]=seqMaior[tbCol-1];
    t->auxMatriz[(t->id-1)*2 + 1][t->pos]=X;  // alinhaGMenor[pos]=X;
    t->tbCol--;
    t->pos++;
  }

  tamAlinha[t->id-1]=t->tamParcial;

  /* o alinhamento foi feito de tras para frente e deve ser
     invertido, conforme segue */
  for (i=0;i<(tamAlinha[t->id-1]/2);i++)
  {
    aux=t->auxMatriz[(t->id-1)*2 + 1][t->pos]; // aux=alinhaGMenor[i];
    t->auxMatriz[(t->id-1)*2 + 1][t->pos]=*t->auxMatriz[tamAlinha[t->id-1]-i-1]; // alinhaGMenor[i]=alinhaGMenor[tamAlinha-i-1];
    t->auxMatriz[tamAlinha[t->id-1]-i-1][t->pos]=aux; // alinhaGMenor[tamAlinha-i-1]=aux;

    aux=t->auxMatriz[(t->id-1)*2][t->pos]; // aux=alinhaGMaior[i];
    t->auxMatriz[(t->id-1)*2][t->pos]=*t->auxMatriz[tamAlinha[t->id-1]-i-1]; // alinhaGMaior[i]=alinhaGMaior[tamAlinha-i-1];
    t->auxMatriz[tamAlinha[t->id-1]-i-1][t->pos]=aux; // alinhaGMaior[tamAlinha-i-1]=aux;
  }


}


/***** TRACEBACK *****/
/* 
  gera o alinhamento global por meio do percurso de retorno na Matriz de escores,
  em duas formas possiveis, conforme o parametro tipo: 1) a partir da celula do
  primeiro maior escore [linPMaior,colPMaior] ou 2) a partir da celula de ultimo
  maior escore [linUMaior,colUMaior], ambos em direcao a celula inicial [0,0],
  uma celula por vez. O caminho de retorno deve ser feito seguindo o mesmo caminho
  inverso que gerou a celular final a partir da celula inicial. O alinhamento
  global eh composto por duas sequencias alinhaGMenor e alinhaGMaior.

  Note que outros alinhamentos globais podem ser encontrados se o percurso do
  traceback for ramificado em celulas que foram "escoreadas/pontuadas" por meio
  de uma estrategia de desempate, pelo fato de ter havido empate na pontuacao
  dos nohs vizinhos. Alem disso, alinhamentos parciais tambem podem ser obtidos
  com traceback iniciado a partir de qualquer c�lula 
*/
void traceBack(int tipo)
{ 
  AlinhamentosParciais=0;
 
  /* inicialização da primeira thread (contendo o alinhamento global) */
  tPackAlinhamento t;
  pthread_t thread;

  if (tipo==1)
  { printf("\nGeracao do Primeiro Maior Alinhamento Global:\n");
    t.tbLin=linPMaior;
    t.tbCol=colPMaior;
  }
  else {
    printf("\nGeracao do Ultimo Maior Alinhamento Global:\n");
    t.tbLin=linUMaior;
    t.tbCol=colUMaior;
  }
  t.id=1;
  t.tamParcial=0;
  t.auxMatriz=alinhaParciais; //////////////////
  t.pos=0;


  int iret = pthread_create(&thread, NULL, traceBackComThreads, (void *) &t);
  pthread_join(thread, NULL);

  printf("\nAlinhamento Global Gerado.");
}


/* mostra os alinhamentos */
void mostraAlinhamentoGlobal(void)
{  
   int i, j;
   
   printf("\nForam obtidos %d alinhamentos parciais, alem do global logo abaixo:", AlinhamentosParciais);

   for(j=0; j<=AlinhamentosParciais; j++)
   {
      if(j==0)
        printf("\nAlinhamento Global Obtido - Tamanho = %d:\n", tamAlinha[j]);
      else 
        printf("\nAlinhamento Parcial Obtido - Tamanho = %d:\n", tamAlinha[j]);
      
      // Imprimindo a sequência alinhada da seqMaior
      for (i=tamAlinha[0]-1; i>=0; i--)
      {
        if(i < tamAlinha[0] - (tamAlinha[j]))
          printf("=");
        else
          printf("%c",mapaBases[alinhaParciais[j*2][i]]);
      }
      printf("\n");

      
      // Imprimindo a sequência alinhada da seqMenor
      for (i=tamAlinha[0]-1; i>=0; i--)
      {
        if(i < tamAlinha[0] - (tamAlinha[j]))
          printf("=");
        else
          printf("%c",mapaBases[alinhaParciais[j*2 + 1][i]]);
      }
      printf("\n");
   }
}
/////////////////////////////////////////////////////////////////////////////




/***** MENU *****/
/* menu de opcoes fornecido para o usuario */
int menuOpcao(void)
{ int op;
  char enter;

  do
  {
    printf("\nMenu de Opcoes:");
    printf("\n<01> Ler Matriz de Pesos");
    printf("\n<02> Mostrar Matriz de Pesos");
    printf("\n<03> Ler Penalidade de Gap");
    printf("\n<04> Mostrar Penalidade");
    printf("\n<05> Definir Sequencias Genomicas");
    printf("\n<06> Mostrar Sequencias");
    printf("\n<07> Gerar Matriz de Escores");
    printf("\n<08> Mostrar Matriz de Escores");
    printf("\n<09> Gerar Alinhamento Global");
    printf("\n<10> Mostrar Alinhamento Global");
    printf("\n<11> Sair");
    printf("\nDigite a opcao => ");
    scanf("%d",&op);
    scanf("%c",&enter);
  } while ((op<1)||(op>sair));

  return (op);
}

/* trata a opcao fornecida pelo usuario, executando o modulo pertinente */
void trataOpcao(int op)
{ int resp;
  char enter;

  switch (op)
  {
    case 1: leMatrizPesos();
            break;
    case 2: mostraMatrizPesos();
            break;
    case 3: penalGap=lePenalidade();
            break;
    case 4: printf("\nPenalidade = %d",penalGap);
            break;
    case 5: printf("\nDeseja Definicao: <1>MANUAL, <2>ALEATORIA ou <3>POR ARQUIVO?  = ");
            scanf("%d",&resp);
            scanf("%c",&enter); /* remove o enter */
            if (resp==1)
            {
              leSequencias();
            }
            else if(resp==2)
            { leTamMaior();
              leTamMenor();
              grauMuta=leGrauMutacao();
              geraSequencias();
            }
            else
            {
              leSequenciasArquivo();
            }
            break;
    case 6: mostraSequencias();
            break;
    case 7: do
            {
              printf("\nCom k <= tamSeqMenor (MAX 1000). Quantas threads voce deseja utilizar? k = ");
              scanf("%d",&k);
              scanf("%c",&enter); /* remove o enter */
            } while(k > tamSeqMenor);
            geraMatrizEscores();
            break;
    case 8: mostraMatrizEscores();
            break;
    case 9: printf("\nDeseja: <1> Primeiro Maior ou <2> Ultimo Maior? = ");
            scanf("%d",&resp);
            scanf("%c",&enter); /* remove o enter */
            traceBack(resp);
            break;
    case 10: mostraAlinhamentoGlobal();
            break;
  }
}
/////////////////////////////////////////////////////////////////////////////




/***** programa principal *****/
void main(void)
{ int opcao;

  srand(time(NULL));

  do
  { printf("\n\n\n---------------------------------------");
    printf("\n Programa Needleman-Wunsch com threads ");
    printf("\n---------------------------------------\n");
    opcao=menuOpcao();
    trataOpcao(opcao);

  } while (opcao!=sair);

}