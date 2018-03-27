//gcc /cygdrive/c/Windows/System32/OpenCL.DLL dZ.c device_info.c wtime.c -o dZ

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <sys/types.h>
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#include <unistd.h>
#else
#include <CL/cl.h>
#endif

#include "err_code.h"

//pick up device type from compiler command line or from
//the default type
#ifndef DEVICE
#define DEVICE CL_DEVICE_TYPE_DEFAULT
#endif


extern double wtime();       // returns time since some fixed past point (wtime.c)
extern int output_device_info(cl_device_id );


//------------------------------------------------------------------------------

#define TOL    (0.001)   // tolerance used in floating point comparisons
#define LENGTH (134217728)    // length of return_vectors a, b, and c

//------------------------------------------------------------------------------

//Start time
double startTime, finalTime;
time_t timer1, timer2;
char buffer1[25], buffer2[25];
struct tm* tm_info;

char * getKernelSource(char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        fprintf(stderr, "Error: Could not open kernel source file\n");
        exit(EXIT_FAILURE);
    }

    fseek(file, 0, SEEK_END);
    int len = ftell(file) + 1;
    rewind(file);

    char *source = (char *)calloc(sizeof(char), len);
    if (!source)
    {
        fprintf(stderr, "Error: Could not allocate memory for source string\n");
        exit(EXIT_FAILURE);
    }
    fread(source, sizeof(char), len, file);
    fclose(file);
    return source;
}

double getRealTime(){
    struct timeval tm;
    gettimeofday(&tm, NULL);
    return ((double)tm.tv_sec + (double)tm.tv_usec/1000000.0);
}

int main(int argc, char** argv) {

    if (argc < 2) {
        printf("Usage : Rendezvous <Número de linhas que serão lidas do arquivo>\n");
        exit(EXIT_FAILURE);
    }

    int          err;               
    int count = 86400;
    float w = 398600.4418/sqrt((6378.0 + 220)*(6378.0 + 220)*(6378.0 + 220));
    int i = 0;
    
    float*       return_vector = (float*) calloc(10, sizeof(float));      
    float*       vector = (float*) calloc(10, sizeof(float));       

    unsigned int correct;           
    size_t global = count;
    
    //Ler arquivo de entrada
    int NPI = atoi(argv[1]);
    FILE *arq, *out;
    char url[] = "in.dat";
    arq = fopen(url, "r");
    out = fopen("serial-out.txt", "w");
    double var1;
    
    // Carregando o codigo fonte do kernel
    char *KernelSource = getKernelSource("dZ_cl.cl");                 

    cl_device_id     device_id;     
    cl_context       context;       
    cl_command_queue commands;      
    cl_program       program;       
    cl_kernel        ko_dZ;       
    
    cl_mem d_return_vector;                     
    cl_mem d_vector;                    
    cl_uint numPlatforms;

    // Buscando platadormas
    err = clGetPlatformIDs(0, NULL, &numPlatforms);
    checkError(err, "Finding platforms");
    if (numPlatforms == 0)
    {
        printf("Found 0 platforms!\n");
        return EXIT_FAILURE;
    }

    cl_platform_id Platform[numPlatforms];
    err = clGetPlatformIDs(numPlatforms, Platform, NULL);
    checkError(err, "Getting platforms");

    for (i = numPlatforms - 1; i < numPlatforms; i++)
    {
        err = clGetDeviceIDs(Platform[i], DEVICE, 1, &device_id, NULL);
        if (err == CL_SUCCESS)
        {
            break;
        }
    }

    if (device_id == NULL)
        checkError(err, "Finding a device");

    err = output_device_info(device_id);
    checkError(err, "Printing device output");

    // Criando contexto
    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
    checkError(err, "Creating context");

    // Criando as command queues
    commands = clCreateCommandQueue(context, device_id, 0, &err);
    checkError(err, "Creating command queue");

    program = clCreateProgramWithSource(context, 1, (const char **) & KernelSource, NULL, &err);
    checkError(err, "Creating program");

    // Compilando os programas
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);

    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048];

        printf("Error: Failed to build program executable!\n%s\n", err_code(err));
        clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("Buffer: %s\n", buffer);
        return EXIT_FAILURE;
    }

    // Usa o programa para criar o kernel
    ko_dZ = clCreateKernel(program, "Kernel_dZ", &err);
    checkError(err, "Creating kernel");

    //Iniciar timer
    time(&timer1);
    tm_info = localtime(&timer1);
    strftime(buffer1, 25, "%d/%m/%Y %H:%M:%S", tm_info);
    startTime = getRealTime();

    for(int np = 1; np <= NPI; np++) {
        
        vector[0] = w;
        //Leitura do arquivo de entrada
        if(arq == NULL) {
            printf("Erro, nao foi possivel abrir o arquivo\n");
            exit(EXIT_FAILURE);
        } else {
            fscanf(arq,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &var1, &var1, &var1, &vector[1], &vector[2], &vector[3], &var1, &vector[4], &vector[5], &vector[6], &var1, &var1, &var1, &var1, &var1, &var1, &var1, &var1, &var1);
        }
     
        // Criando os buffers de entrada e saida na memória do device
        d_return_vector  = clCreateBuffer(context,  CL_MEM_READ_ONLY,  sizeof(float) * 10, NULL, &err);
        checkError(err, "Creating buffer d_return_vector");

        d_vector  = clCreateBuffer(context,  CL_MEM_READ_WRITE, sizeof(float) * 10, NULL, &err);
        checkError(err, "Creating buffer d_vector");

        //Copiando valores no host para o buffer na memoria do device
        err = clEnqueueWriteBuffer(commands, d_vector, CL_TRUE, 0, sizeof(float) * 10, vector, 0, NULL, NULL);
        checkError(err, "Coping vector to buffer d_vector");

        // Definindo os argumentos para o kernel
        err |= clSetKernelArg(ko_dZ, 0, sizeof(cl_mem), &d_vector);
        err |= clSetKernelArg(ko_dZ, 1, sizeof(cl_mem), &d_return_vector);
        checkError(err, "Setting kernel arguments"); 

        // Executa o kernel
        err = clEnqueueNDRangeKernel(commands, ko_dZ, 1, NULL, &global, NULL, 0, NULL, NULL);
        checkError(err, "Enqueueing kernel");

        // Espera os comandos concluirem para parar o timer
        err = clFinish(commands);
        checkError(err, "Waiting for kernel to finish");

        // Pegar o resultado da GPU
        err = clEnqueueReadBuffer( commands, d_return_vector, CL_TRUE, 0, sizeof(float) * 10, return_vector, 0, NULL, NULL );  
        if (err != CL_SUCCESS)
        {
            printf("Error: Failed to read output array!\n%s\n", err_code(err));
            exit(1);
        }
    }
    // Espera os comandos concluirem para parar o timer
    time(&timer2);
    tm_info = localtime(&timer2);
    strftime(buffer2, 25, "%d/%m/%Y %H:%M:%S", tm_info);
    
    finalTime = getRealTime();

    printf("\nTempo de execucao: %f segundos\n",finalTime-startTime);

    for(i = 0; i < 10; i++)
    {
        float tmp = vector[i];
        //printf("%.50f\n", tmp);
    }
 
    clReleaseMemObject(d_return_vector);
    clReleaseMemObject(d_vector);
    clReleaseProgram(program);
    clReleaseKernel(ko_dZ);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);

    //free(number_of_lines);
    free(return_vector);
    free(vector);

    return 0;
}