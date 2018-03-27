#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

double startTime, finalTime;
time_t timer1, timer2;
char buffer1[25], buffer2[25];
struct tm* tm_info;

//-------------------Functions-------------------
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
void main(int argc, char *argv[]) {

    if (argc < 3) {
        printf("Usage : Rendezvous <Número de linhas que serão lidas do arquivo> <Numero de cores>\n");
        exit(EXIT_FAILURE);
    }
    //otimizacao ----------------------
    ww = w*w;

    //---------------------------------
    //Start time
    time(&timer1);
    tm_info = localtime(&timer1);
    strftime(buffer1, 25, "%d/%m/%Y %H:%M:%S", tm_info);
    startTime = getRealTime();

    int Tmax = 86400;
    int NPI = atoi(argv[1]); // numero de posicoes iniciais
    int num_cores = atoi(argv[2]);

    FILE *arq, *out;
    char url[] = "in.dat";
    arq = fopen(url, "r");
    out = fopen("parallel-out.txt", "w");
    double var1;

    omp_set_num_threads(num_cores);
    
    int tid = omp_get_max_threads();
    printf("Ccongigurações definicas \n");
    printf("Numero de cores: %d\n", num_cores);
    printf("Numero de maximo threads: %d\n", tid);

    printf("executando...\n");

    for(int np = 1; np <= NPI; np++) {
        //printf("Problema %d\n", np);
        if(arq == NULL) {
            printf("Erro, nao foi possivel abrir o arquivo\n");
            exit(EXIT_FAILURE);
        } else {
            fscanf(arq,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &var1, &var1, &var1, &x, &y, &z, &var1, &xl0, &yl0, &zl0, &var1, &var1, &var1, &var1, &var1, &var1, &var1, &var1, &var1);
        }
        
        #pragma omp parallel for
        for(int VeAux = 1; VeAux<=10; VeAux++) {

            double Ve =VeAux;
            Ve = Ve/2;
            double vex, vey, vez;
            vex = vey = vez =Ve*Ve/3;

            //otimizacao -------------------------------
            
            double vex2_w = (2*vex)/w;
            
            double vey_w = vey/w;
            
            double vex3 = vex*3;
            
            double vey2_w = vex2_w;
            double vex4 = vex*4;
            //------------------------------------------

            for(int aux = -14; aux<=2; aux++){
                
                double gama = pow(10, aux);

                //otimizacao -------------------------------
                double gama_w = gama/w;

                double gamavex_ww = (gama*vex)/ww;
 
                double vexgama = vex*gama;
 
                double gamavey_ww = gamavex_ww;

                double gama_wpow = (gama/w)*(gama/w);
                
                for(int Xaux=1; Xaux<=100; Xaux++) {
                    
                    double X = Xaux;                   

                    double G = brute_G (x, yl0, X, vex, vey, vex3);
                    double H = brute_H (z, gama, vex, vexgama, gama_wpow);
                    double I = brute_I (zl0, gama, X, vez, gama_wpow);
                
                    for(int t = 0; t <= Tmax; t++) {
                        
                        double fz = dZ(t, X, gama, vez, G, H, I, gama_wpow);
                
                    }
                }
            }
        }
    }
    time(&timer2);
    tm_info = localtime(&timer2);
    strftime(buffer2, 25, "%d/%m/%Y %H:%M:%S", tm_info);
    
    finalTime = getRealTime();

    fprintf(out,"\nTempo de execucao: %f segundos\n",finalTime-startTime);
    fclose(out);
    printf("concluido!\n");
}

double brute_G (double x0, double yl0, double X, double vex, double vey, double vex3) {
    double result = 0;
    double sum = 0;
    double aux;

    result= 2*yl0/w + x0 + (2*vey*(log((X+1)/X)))/w;
    //#pragma omp parallel for reduction(+:sum) private(aux)
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
    for (int n = 1; n <= N; n++) {
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


double dZ (int t, double X, double gama, double vez, double G, double H, double I, double gama_wpow) {
    double wt = w*t;
    double gamat = gama*t;

    double resultJn = 0;
    double result1 = H * cos(wt) + I * sin(wt);
    double result2 = 0;

    for (int n = 1; n <= N; n++) {
        //brute_J
        resultJn = vez/(n*pow(X,n)*w)/(1+(n*n*gama_wpow));
        if (n%2 == 0) {
            resultJn = -resultJn;
        }

        result2 += resultJn * pow(M_E, -(n * gamat));
    }
    return result1 - result2;
}
