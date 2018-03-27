#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "mpi.h"
//#include <omp.h>

double startTime, finalTime;
time_t timer1, timer2;
char buffer1[25], buffer2[25];
struct tm* tm_info;

//-------------------Functions-------------------
double vZ(int t, double X, double gama,  double vez, double H, double I, double gama_wpow);

double rT(double fx, double fy, double fz);
double vT(double dx, double dy, double dz);

double dZ (int t, double X, double gama, double vez, double G, double H, double I, double gama_wpow);

double brute_G (double x0, double yl0, double X, double vex, double vey, double vex3);
double brute_H (double z0, double gama, double vex, double vexgama, double gama_wpow);
double brute_I (double zl0, double gama, double X, double vez, double gama_wpow);

double x=0, y=0, z=0, xl0=0, yl0=0, zl0=0;
int Alt= 220;
double w = 398600.4418/sqrt((6378.0 + 220)*(6378.0 + 220)*(6378.0 + 220));
//otimizacao ---------------------
double ww;
//--------------------------------

double getRealTime(){
    struct timeval tm;
    gettimeofday(&tm, NULL);
    return ((double)tm.tv_sec + (double)tm.tv_usec/1000000.0);
}

static int const N = 20;
double nave = 0;

/*
 * main
 */
void main(int argc, char** argv) {

    if (argc < 2) {
        printf("Usage : Rendezvous <Número de linhas que serão lidas do arquivo>\n");
        exit(EXIT_FAILURE);
    }

    int Tmax = 86400;
    //float z_out[43200][10];

    //otimizacao ----------------------
    ww = w*w;
    //Start time
    time(&timer1);
    tm_info = localtime(&timer1);
    strftime(buffer1, 25, "%d/%m/%Y %H:%M:%S", tm_info);
    startTime = getRealTime(); 

    //---------------------------------

    int NPI = atoi(argv[1]); // numero de posicoes iniciais
    FILE *arq, *out;
    char url[] = "in.dat";
    arq = fopen(url, "r");
    out = fopen("serial-out.txt", "w");
    double var1;
    int time_start, time_end;

    int  meu_id;                                      /*identificador do processo*/
    int  numero_processos;                            /*numero de processos*/
    int  origem;                                      /*rank de quem envia*/
    int  destino;                                     /*rank de quem recebe*/
    int  tag=0;                                       /*identificador de mensagens*/
    double mensagem[100];                             /*mensagem (integer))*/
    char nome_host [MPI_MAX_PROCESSOR_NAME];          /*Comprimento maximo do nome do host retornado por MPI_Get_processor_name*/
    int  tamanho_nome;                                /*Comprimento (em caracteres) do nome*/
    double hora_inicio, hora_fim;                     /*Hora de inicio e fim do processamento*/
    double hora_inicio_total, hora_fim_total;         /*Hora de inicio e fim do processamento total*/
    int contador;                                     /*contador*/
    MPI_Status status;

    /*Inicializa MPI*/
    MPI_Init(&argc, &argv);
    /*Descobre quem sou eu -- rank*/
    MPI_Comm_rank(MPI_COMM_WORLD,&meu_id);
    /*Descobre quanto processos existem*/
    MPI_Comm_size(MPI_COMM_WORLD,&numero_processos);
    /*Descobre meu nome*/
    MPI_Get_processor_name(nome_host,&tamanho_nome);
    
    printf("executando...\n");

    for(int np = 1; np <= NPI; np++) {
        
        time_start = 0;
        time_end = 0;

        if(arq == NULL) {
            printf("Erro, nao foi possivel abrir o arquivo\n");
            exit(EXIT_FAILURE);
        } else {
            fscanf(arq,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &var1, &var1, &var1, &x, &y, &z, &var1, &xl0, &yl0, &zl0, &var1, &var1, &var1, &var1, &var1, &var1, &var1, &var1, &var1);
        }

        if(meu_id==0) { /*se eu sou o root*/
         
            printf("\nMestre:%s  Processo %d\n", nome_host,meu_id);

            for(destino=1;destino<numero_processos;destino++){  /*Inicio do for*/

                if(86400%(numero_processos - 1) == 0) {
                    time_start = time_end;
                    time_end = time_end + 86400/(numero_processos - 1);
                } else {
                    if(destino == (numero_processos - 1)) {
                        time_start = time_end;
                        time_end = 86400;
                    } else {
                        time_start = time_end;
                        time_end = time_end + 86400/(numero_processos - 1);
                    }
                }

                mensagem[1]=t_start;
                mensagem[2]=t_end;

                MPI_Send(mensagem,100,MPI_INT,destino,tag,MPI_COMM_WORLD);
                printf("Enviei trabalho para o processo:%d de %d \n",destino, numero_processos);

            }  /*Final do for*/

            for(destino=1;destino<numero_processos;destino++){  /*Inicio do for*/
                
                MPI_Recv(mensagem,100,MPI_DOUBLE,source,tag,MPI_COMM_WORLD,&status);   /*Recebe mensagem*/
                MPI_Recv(mensagem1,100,MPI_CHAR,destino,tag,MPI_COMM_WORLD,&status); /*Recebe mensagem - tempo de processamento*/

                printf("\n%s \n",mensagem1);                                         /*Mostra mensagem - tempo de processamento*/
                printf("Resultado do processamento: %f \n\n",mensagem[2]);
            }
            time_end = 0;
            time_start = 0;

            /*Descobre hora final- processamento total*/
            hora_fim_total = MPI_Wtime();
            printf("Mestre de novo:(%s)\n", nome_host);
            printf("Tempo de processamento total do cluster (computo+comunicação): %f Milisegundos !\n\n",((hora_fim_total-hora_inicio_total)*1000)); /*Calcula - tempo de processamento total*/

    }    /*Final do if*/

    else {  /*se não sou o root*/
        hora_inicio = MPI_Wtime(); /*MPI_Wtime - retorna um número de ponto flutuante (tempo global em segundos),*/
        origem=0;
        MPI_Recv(mensagem,100,MPI_DOUBLE,origem,tag,MPI_COMM_WORLD,&status);   /*Recebe mensagem*/

           //mensagem[2]=sqrt(mensagem[1]);
           //#pragma omp parallel for
        for(double Ve = 0.5; Ve<=5; Ve+=0.5) {
            //printf("Ve %d\n", Ve);
            double vex, vey, vez;
            vex = vey = vez =Ve*Ve/3;

            //otimizacao -------------------------------
            //brute_A
            double vex2_w = (2*vex)/w;
            //brute_B
            double vey_w = vey/w;
            //brute_G
            double vex3 = vex*3;
            //dx vx
            double vey2_w = vex2_w;
            double vex4 = vex*4;
            //------------------------------------------

            //#pragma omp parallel for
            for(int aux = -14; aux<=2; aux++){
                //printf("Gama %d\n", aux);
                double gama = pow(10, aux);

                //otimizacao -------------------------------
                //brute_A
                double gama_w = gama/w;
                //brute_B
                double gamavex_ww = (gama*vex)/ww;
                //brute_H
                double vexgama = vex*gama;
                //dy vy A
                double gamavey_ww = gamavex_ww;
                //vx vy vz dx dy dz B H I
                double gama_wpow = (gama/w)*(gama/w);
                
                for(int Xaux=1; Xaux<=100; Xaux++) {
                    //printf("X %d\n", Xaux);
                    double X = Xaux;                   

                    double G = brute_G (x, yl0, X, vex, vey, vex3);
                    double H = brute_H (z, gama, vex, vexgama, gama_wpow);
                    double I = brute_I (zl0, gama, X, vez, gama_wpow);
                
                    for(int t = mensagem[1]; t <= mensagem[2]; t++) {
                        //printf("t %d\n", t);
                        double fz = dZ(t, X, gama, vez, G, H, I, gama_wpow);
                
                    }
                }
            }
        }
    }

    hora_fim = MPI_Wtime(); /*MPI_Wtime - retorna um número de ponto flutuante (tempo global em segundos),*/
    MPI_Send(mensagem,100,MPI_INT,origem,tag,MPI_COMM_WORLD);            /*Envia mensagem para origem*/

    /*Inicio: calculo e envio do tempo de processamento*/
    sprintf(mensagem1,"Trabalhador:%s\nTempo de processamento foi de %f Milisegundos !",nome_host,((hora_fim-hora_inicio)*1000)); /*Calcula - tempo de processamento*/
    MPI_Send(mensagem1,strlen(mensagem1)+1,MPI_CHAR,origem,tag,MPI_COMM_WORLD); /*Envia mensagem - tempo de processamento*/
    /*Final: calculo e envio do tempo de processamento*/

    }   /*Final do else*/

        
    time(&timer2);
    tm_info = localtime(&timer2);
    strftime(buffer2, 25, "%d/%m/%Y %H:%M:%S", tm_info);
    
    finalTime = getRealTime();

    fprintf(out,"\nTempo de Execucao: %f segundos\n",finalTime-startTime);
    fclose(out);
    printf("concluido!\n");
}

double brute_G (double x0, double yl0, double X, double vex, double vey, double vex3) {
    double result = 0;
    double sum = 0;
    double aux;

    result= 2*yl0/w + x0 + (2*vey*(log((X+1)/X)))/w;
    ////#pragma omp parallel for reduction(+:sum) private(aux)
    for (int n = 1; n <= N; n++) {
        aux = vex3/(n*n*pow(X,n)*w);

        if (n%2 == 0) {
            aux = -aux;
        }
        sum +=aux;
    }

    result-=sum;

    return result;
}


double brute_H (double z0, double gama, double vex, double vexgama, double gama_wpow) {
    double result = 0;
    double sum = 0;
    double aux;

    result = z0;
    //Calculo do somatorio
    ////#pragma omp parallel for reduction(+:sum) private(aux)
    for (int n = 1; n <= N; n++) {
        //aux = ((vex*gama)/(pow(gama,n)*(ww)))/(1+((n*gama)/w)*((n*gama)/w));
        aux = ((vexgama)/(pow(gama,n)*(ww)))/(1+(n*n*gama_wpow));
        
        if (n%2 == 0) {
            aux = -aux;
        }
        sum += aux;
    }
    result += sum;

    return result;
}


double brute_I (double zl0, double gama, double X, double vez, double gama_wpow) {
    double result = 0;
    double sum = 0;
    double aux;

    result = zl0/w - (vez/w)*(log((X+1)/X));

    //Calculo do somatorio
    ////#pragma omp parallel for reduction(+:sum) private(aux)
    for (int n = 1; n <= N; n++) {
        aux = ((vez)/(n*n*pow(X,n)*w))/(1+(n*n*gama_wpow));
        if (n%2 == 0) {
            aux = -aux;
        }
        sum += aux;
    }

    result += sum;

    return result;
}


/* @author Weverson, Iago
 * vetor Z da distancia
 */
double dZ (int t, double X, double gama, double vez, double G, double H, double I, double gama_wpow) {
    //otimizacao
    double wt = w*t;
    double gamat = gama*t;

    double resultJn = 0;
    double result1 = H * cos(wt) + I * sin(wt);
    double result2 = 0;

    ////#pragma omp parallel for reduction(+:result2)
    for (int n = 1; n <= N; n++) {
        //brute_J
        resultJn = vez/(n*pow(X,n)*w)/(1+(n*n*gama_wpow));
        if (n%2 == 0) {
            resultJn = -resultJn;
        }
        //brute_J

        result2 += resultJn * pow(M_E, -(n * gamat));
    }
    return result1 - result2;
}

/* @author Weverson
 * Modificada por Gledson
 * Velocidade
 */
double rT(double fx, double fy, double fz) {
    return sqrt(fx*fx + fy*fy + fz*fz);
}

/* @author Weverson
 * Modificada por Gledson
 * Posicao
 */
double vT(double dx, double dy, double dz) {
    return sqrt(dx*dx + dy*dy + dz*dz);
}
