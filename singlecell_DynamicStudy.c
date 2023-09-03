//
//  main.c
//  RK2OPT_multicell
//
//  Created by Pilar Gil
//

/* PROPIEDADES DEL PROGRAMA:
 Evolución temporal de N agentes sin interacción.
 Integración de la ecuación diferencial con RK2, interpolación con punto medio.
 Arrays temporales de longitud reducida para optimizar memoria.
 Implementación de cálculo del parámetro de orden R y la frecuencia de oscilación.
 Cálculo de la solución estacionaria y estudio de estabilidad de la misma.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

// - NUMEROS ALEATORIOS -- START

#define PI 3.14159265
#define NormaRanu (2.3283063671E-10F)
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;

double RapuanoRandom(void)
{
    double x;
    ig1=ind_ran-24;
    ig2=ind_ran-55;
    ig3=ind_ran-61;

    irr[ind_ran]=(irr[ig1]+irr[ig2]);
    ir1=(irr[ind_ran]^irr[ig3]);
    ind_ran++;

    x=(ir1*NormaRanu);
    return x;
}

void IniRapu(void)
{
    int i;
    double maxx;
    srand(time(NULL));
    maxx=4294967296;    //2^32
    for (i=0; i<256; i++)
        irr[i]=maxx*(rand()/((double)RAND_MAX+1));

    ind_ran=ig1=ig2=ig3=0;
    for (i=0; i<1024; i++)
        RapuanoRandom();
}

// - NUMEROS ALEATORIOS -- END

// PARÁMETROS GLOBALES:
int N_AGENTS = 50; // number of cells
double tf = 1000; // tomamos como tiempo inicial t0 = 0
double h = 0.01;

/** FUNCIÓN A RESOLVER:
    f(t, y(t), y(t-tau))
    *X: Array que contiene la evolucion temporal de la funcion incognita, X[time][cell]
    i: iteración de RK2 (sirve para acceder a la solución en el tiempo desfasado tau
**/
void f(double output[N_AGENTS], double x[][N_AGENTS], double y[N_AGENTS], int i, double tau, double nHill, double K, double alpha, double Km, double V, double gam)
{
    int n_tau = tau/h;
    int index_delayed = i-n_tau;
    double aux1, aux2;
    for(int j=0; j<N_AGENTS; j++){
        //output[j] = -(x[index_delayed][j])*0.5; // ECUACIÓN SENCILLA
        output[j] = alpha/(1+pow(x[index_delayed][j]/K,nHill))-V*y[j]/(Km+y[j])-gam*y[j];
    }
}

/** RK
    Algoritmo Runge-Kutta de orden 2 para sistemas generales de N ecuaciones con N incognitas
    *f:     Puntero a la funcion vectorial
    t0:     Tiempo inicial
    tf:      Tiempo final
    x0:    Condición inicial
    N:     Dimension del sistema (numero de ecuaciones y de incognitas); dimension de los arrays x0, t, x, y
    X[time][cell]:     Array que contiene la evolucion temporal de la funcion incognita
**/
void RK2(double *t, double X[][N_AGENTS], double tau, double nHill, double K, double alpha, double Km, double V, double gam)
{
    int n_tau = tau/h; //Numero de puntos almacenados con las condiciones iniciales

    double k1[N_AGENTS];
    double k2[N_AGENTS];
    double kAux[N_AGENTS];

    int i, j;

    for (i=n_tau; i<2*n_tau; i++) // Usar el punto i-esimo para calcular el (i+1)-esimo
    {
        
        f(k1, X, X[i], i, tau, nHill, K, alpha, Km, V, gam); // k1 = f(t_i,y_i)
        
        for (j=0; j<N_AGENTS; j++){
            kAux[j] = X[i][j] + k1[j]*h;
        }
        f(k2, X, kAux, i+1, tau, nHill, K, alpha, Km, V, gam); // k2 = f(t_i+h,y_i+h*k1)
        
        for (j=0; j<N_AGENTS; j++){
            X[i+1][j] = X[i][j] + (k1[j]+k2[j])*h*0.5;
        }
    }
}

double f_stat(double x, double tau, double nHill, double K, double alpha, double Km, double V, double gam){ // FUNCIÓN f_stat DE LA QUE HALLAR RAÍZ (solución estacionaria del sistema)
    return alpha/(1+pow(x/K,nHill))-V*x/(Km+x)-gam*x;
}

double g_stat(double x, double tau, double nHill, double K, double alpha, double Km, double V, double gam){ // DERIVADA DE f_stat
    return -alpha*nHill*pow(x/K,nHill)/(x*pow(1+pow(x/K,nHill),2))-V*Km/pow(Km+x,2)-gam;
}

// FUNCIÓN CARACTERÍSTICA DEL SISTEMA: f =  f[0] (parte real) + i*f[1] (parte imaginaria)
void f_carac(double f[2], double x[2], double a, double b, double tau){ // COMP 2 de FUNCIÓN f
    f[0] = x[0]-b*exp(-x[0]*tau)*cos(x[1]*tau)-a;
    f[1] = x[1]+b*exp(-x[0]*tau)*sin(x[1]*tau);
}

// JACOBIANO DE LA FUNCIÓN f: matriz 2x2
void Jf_carac(double Jf[2][2], double x[2], double a, double b, double tau){ // DERIVADA DE f
    Jf[0][0] = 1+b*tau*exp(-x[0]*tau)*cos(x[1]*tau);
    Jf[0][1] = b*tau*exp(-x[0]*tau)*sin(x[1]*tau);
    Jf[1][0] = -b*tau*exp(-x[0]*tau)*sin(x[1]*tau);
    Jf[1][1] = 1+b*tau*exp(-x[0]*tau)*cos(x[1]*tau);
}

int main()
{
    IniRapu();
    
    // PARÁMETROS DE LA ECUACIÓN
    double tau = 3.0;
    double nHill = 2.0;
    double K = 0.15;
    double alpha = 1.0;
    double Km = 50;
    double V = 1.0;
    double gam = 1.0; //Normalización implícita respecto al eje temporal
    
    // variables para integración
    int n_tau = tau/h; //Numero de puntos almacenados con las condiciones iniciales
    int n_time = tf/h; //Numero de puntos a calcular sin contar con las condiciones iniciales
    int n_total = n_tau + n_time; //Numero total de puntos que almacenamos
    int iter_max = ceil(tf/tau); //Numero de iteraciones de RK2 (por intervalos de tiempo de longitud tau) necesarias para resolver el sistema
    
    double t[2*n_tau]; // almacenamos solo el intervalo actual y el anterior
    double x_cells[2*n_tau+1][N_AGENTS];
    printf("%.3d \n",n_total);
    
    //INICIALIZACIÓN DEL SISTEMA: entre -tau y 0, cada célula toma valor constante x0 asignado aleatoriamente
    double aux=0;
    for (int j=0; j<N_AGENTS; j++){
        aux = 0.5*RapuanoRandom(); //0.2 + 0.02*RapuanoRandom(); // condiciones iniciales próximas a solución estacionaria
        for (int i=0; i<=n_tau; i++){
            x_cells[i][j]=aux;
        }
    }
    
    //FICHERO con EVOLUCIÓN TEMPORAL de CONCENTRACIONES
    FILE *g;
    g=fopen("evol_temp.txt","w");
    fprintf(g,"#VARIABLES DEL SISTEMA:\n");
    fprintf(g,"#tau = %.3f\t nHill = %.3f\t alpha = %.3f\t K = %.3f\t V = %.3f\t Km = %.3f\t gamma = %.3f\t\n", tau, nHill, alpha, K, V, Km, gam);
    fprintf(g,"#t\t");
    for(int i=0; i<N_AGENTS; i++){
        fprintf(g,"x[%d]\t",i);
    }
    fprintf(g,"\n");
    
    // Pasar concentraciones de condiciones iniciales a fichero
    for(int i=0; i<n_tau; i++){
        t[i]=i*h-tau;
        fprintf(g,"%.3f\t",t[i]);
        for(int j=0; j<N_AGENTS; j++){
            fprintf(g,"%.3f\t",x_cells[i][j]);
        }
        fprintf(g,"\n");
    }
    
    //FICHERO con EVOLUCIÓN TEMPORAL de PARÁMETRO DE ORDEN
    FILE *g_order;
    g_order=fopen("order_R.txt","w");
    fprintf(g_order,"#VARIABLES DEL SISTEMA:\n");
    fprintf(g_order,"#tau = %.3f\t nHill = %.3f\t alpha = %.3f\t K = %.3f\t V = %.3f\t Km = %.3f\t gamma = %.3f\t\n", tau, nHill, alpha, K, V, Km, gam);
    fprintf(g_order,"#t_order\tR\n");
    
    int iter_term = ceil(0.1*tf/tau); // TERMALIZACIÓN: iteración del total de simulacion a partir de la cual calculamos R
    
    // VARIABLES PARÁMETRO DE ORDEN
    double numt=0, denomt=0;
    double M=0, M_sumt=0, M2_sumt=0; // En M_sumt almacenamos la suma de M para el promedio temporal. Íden con M2_sumt y la suma de cuadrados.
    double x_sumt[N_AGENTS], x2_sumt[N_AGENTS]; // En x_sumt almacenamos la suma de x para el promedio temporal. Íden con x2_sumt y la suma de cuadrados.
    for(int j=0; j<N_AGENTS; j++){
        x_sumt[j] = 0;
        x2_sumt[j] = 0;
    }
    
    //------------- BÚSQUEDA SOLUCIÓN ESTACIONARIA ------------- START
    
    float x0 = 0.3, x1, f0, f1, g0;
    float tol = 10e-8; // tolerance
    int itmax = 10e2; // maximum number of iterations
    float x_stat;
    
    printf("\nSearch for value of stable concentration.\n");
    printf("it\t x0\t f0\n");
    
    for (int i=0; i<itmax; i++){
        f0 = f_stat(x0, tau, nHill, K, alpha, Km, V, gam);
        g0 = g_stat(x0, tau, nHill, K, alpha, Km, V, gam);
        
        if(g0 == 0.0){
            printf("Mathematical Error.\n");
            return 0;
        }
        
        x1 = x0 - f0/g0; // ITERACIÓN k-ésima
        
        printf("%d\t\t%f\t%f\n",i,x0,f0);
        x0 = x1;
        
        f1 = f_stat(x1, tau, nHill, K, alpha, Km, V, gam);
        if (fabs(f1-f0)<tol){
            printf("Solution reached.");
            break;
        }
    }
    
    x_stat = x1;
    printf("\nStable solution for concentration x=%f\n", x1);
    
    //------------- BÚSQUEDA SOLUCIÓN ESTACIONARIA ------------- END
    
    //------------- CÁLCULO DE LA SOLUCIÓN NUMÉRICA ------------- START
    
    //INTEGRACIÓN DE LA ED POR INTERVALOS DE TIEMPO DE LONGITUD tau
    int steps_t=1; // Número de puntos usados para calcular el promedio temporal (contador de iteraciones global)
    int n_cycles=0;
    double x_max=0, x_min=1, amplitude=0;
    
    for (int cont_iter=0; cont_iter<iter_max; cont_iter++)
    {
        // Actualizar array de tiempo en la iteracion cont_iter
        for (int i=0; i<=2*n_tau; i++){
            t[i]=cont_iter*tau+(i*h-tau);
        }
        
        // Ejecutar RK2 para cada tramo para calcular x_cells en [cont_iter*tau, (cont_iter+1)*tau]
        RK2(t, x_cells, tau, nHill, K, alpha, Km, V, gam);
        
        // Pasar concentraciones a fichero
        for(int i=n_tau; i<=2*n_tau; i++){
            fprintf(g,"%.3f\t",t[i]);
            for(int j=0; j<N_AGENTS; j++){
                fprintf(g,"%.3f\t",x_cells[i][j]);
            }
            fprintf(g,"\n");
        }
        
        if(cont_iter>=iter_term){
            for (int i=n_tau; i<2*n_tau; i++){
                //printf("steps_t: %d\t", steps_t);
                
                //NUMERADOR
                for(int j=0; j<N_AGENTS; j++){
                    M += x_cells[i][j]; // Cálculo de media de concentraciones en tiempo t[i]
                }
                M = M/N_AGENTS;
                M_sumt += M;
                M2_sumt += pow(M,2);
                numt = M2_sumt/steps_t-pow(M_sumt/steps_t,2);
                
                //DENOMINADOR
                for(int j=0; j<N_AGENTS; j++){
                    x_sumt[j] += x_cells[i][j];
                    x2_sumt[j] += pow(x_cells[i][j],2);
                    denomt += x2_sumt[j]/steps_t - pow(x_sumt[j]/steps_t,2);
                }
                denomt = denomt/N_AGENTS;
                
                // Pasar parámetro de orden a fichero
                if (i>n_tau){
                    fprintf(g_order,"%.3f\t%.3f\n", t[i], numt/denomt);
                    //printf("t[i]=%.2f\tR=%.3f\n", t[i], numt/denomt);
                }
                
                M = 0;
                denomt = 0;
                
                //CÁLCULO DE LA FRECUENCIA PROMEDIO DE OSCILACIÓN
                if (t[i]>700 && t[i]<800){
                    steps_t++;
                    for(int j=0; j<N_AGENTS; j++){
                        if (x_cells[i][j]<x_stat && x_cells[i+1][j]>=x_stat){
                            n_cycles++;
                        }
                        if (x_cells[i+1][j]>x_max)
                            x_max = x_cells[i+1][j];
                        if (x_cells[i+1][j]<x_min)
                            x_min = x_cells[i+1][j];
                    }
                }
                amplitude = x_max-x_min;
            }
        }
        
        // Actualizar array de concentraciones
        for(int i=0; i<=n_tau; i++){ //MENOR O IGUAL O MENOR ESTRICTO (?)
            for(int j=0; j<N_AGENTS; j++){
                x_cells[i][j] = x_cells[i+n_tau][j];
            }
        }
    }
    
    double omega = n_cycles/(steps_t*h*N_AGENTS);
    
    printf("\nAmplitud oscilación = %.3f", amplitude);
    //if (amplitude<0.05)
    //    omega=0;
    
    printf("\nFrecuencia promedio = %.3f\n\n", omega);
    fprintf(g,"#Frecuencia promedio = %.3f\n", omega);
    
    //------------- CÁLCULO DE LA SOLUCIÓN NUMÉRICA ------------- END
    
    //------------- ESTUDIO ESTABILIDAD ------------- START
    
    double u0[2]={-0.3, -0.5}, u1[2]={0};
    double f[2] = {0}, f_test[2] = {0};
    double Jf[2][2] = {{0, 0},{0, 0}}, Jf_inv[2][2] = {{0, 0},{0, 0}};
    double detJf=0;
    
    // Cálculo de coeficientes de la ecuación característica
    double a = - gam - V*Km/pow(Km+x_stat,2);
    double b = -alpha*nHill*pow(x_stat/K,nHill)/(x_stat*(1+pow(pow(x_stat/K,nHill),2)));
    
    printf("\nSearch for eigenvalue of characteristic equation.\n");
    printf("\nit\t x0\t f0\n");
    
    for (int i=0; i<itmax; i++){
        
        f_carac(f, u0, a, b, tau);
        Jf_carac(Jf, u0, a, b, tau);
        detJf = Jf[0][0]*Jf[1][1]-Jf[0][1]*Jf[1][0];
        
        if(detJf == 0.0){
            printf("Mathematical Error.\n");
            return 0;
        }
        
        // CALCULO DE LA INVERSA DEL JACOBIANO (componente a componente)
        Jf_inv[0][0] = Jf[1][1]/detJf;
        Jf_inv[0][1] = -Jf[0][1]/detJf;
        Jf_inv[1][0] = -Jf[1][0]/detJf;
        Jf_inv[1][1] = Jf[0][0]/detJf;
        
        // ITERACIÓN k-ésima
        u1[0] = u0[0] -(Jf_inv[0][0]*f[0] + Jf_inv[0][1]*f[1]);
        u1[1] = u0[1] -(Jf_inv[1][0]*f[0] + Jf_inv[1][1]*f[1]);
        
        printf("%d\t%f + i%f\t%f + i%f\n",i,u0[0],u0[1],f[0],f[1]);
        u0[0] = u1[0]; u0[1] = u1[1];
        
        f_carac(f_test, u1, a, b, tau);
        if (sqrt(pow(f[0]-f_test[0],2)+pow(f[1]-f_test[1],2)) < tol){
            printf("Solution reached for a = %.3f, b = %.3f.\n\n", a, b);
            break;
        }
    }
    
    printf("Root from characteristic eq. is: %f + i*%f\n\n", u1[0],u1[1]);
    
    //------------- ESTUDIO ESTABILIDAD ------------- END
    
    //------------- CHECK SOLUTION ------------- START
    
    double Eval[2]={0};
    Eval[0] = u1[0]-b*exp(-tau*u1[0])*cos(tau*u1[1]) - a;
    Eval[1] = u1[1]+b*exp(-tau*u1[0])*sin(tau*u1[1]);
    // f_carac(Eval[2], u1[2], a, b, tau)
    printf("COMPROBATION: %f + i*%f\n\n", Eval[0],Eval[1]);
    
    //------------- CHECK SOLUTION ------------- END
    
    printf("All good.\n");
    
    fclose(g);
    fclose(g_order);

    return 0;
}


