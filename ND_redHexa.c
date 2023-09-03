//
//  main.c
//  NotchDelay_redHexa
//
//  Created by Pilar Gil
//

/* PROPIEDADES DEL PROGRAMA:
 Evolución temporal de una red hexagonal de células con interacción propuesta por Tomka.
 En principio cada célula tiene un delay efectivo diferente.
 Integración de la ecuación diferencial con RK2, interpolación con punto medio.
 Arrays temporales de longitud reducida para optimizar memoria.
 Implementación de cálculo del parámetro de orden R.
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

#define Lx 20 // Longitud x de la red en número de células
#define Ly 20 // Longitud y de la red en número de células
#define N 400 // Número total de células (N=Lx*Ly)

// PARÁMETROS GLOBALES:
int N_AGENTS = N; // number of cells
double tf = 1500; // tomamos como tiempo inicial t0 = 0
double h = 0.01;

// ALGORITMO PARA DETERMINAR LOS CLOSEST NEIGHBORS DE i-ésimo ELTO. EN RED HEXAGONAL
void CondContorno(int*v1, int*v2, int*v3, int*v4, int*v5, int*v6, int i)
{
    int lx,ly;
    ly= (int)(i/Lx);
    lx= i-(ly*Lx);
    
    //CASO 1a: Célula en una columna interior par
    if (lx>0 && lx<(Lx-1) && ly>0 && ly<(Ly-1) && (int)(lx%2)==0)
    {
        *v1=i+1;
        *v2=i-1;
        *v3=i+Lx;
        *v4=i-Lx;
        *v5=i-(Lx+1);
        *v6=i-(Lx-1);
    }
    
    //CASO 1b: Célula en una columna interior impar
    else if (lx>0 && lx<(Lx-1) && ly>0 && ly<(Ly-1)&& (int)(lx%2)!=0)
    {
        *v1=i+1;
        *v2=i-1;
        *v3=i+Lx;
        *v4=i-Lx;
        *v5=i+(Lx+1);
        *v6=i+(Lx-1);
    }
    
    //CASO 2a: Célula en el borde dcho. o izq. de la red (exceptuando esquinas)
    else if ((lx==0 && i!=0 && i!=(N-Lx)) || (lx==(Lx-1) && i!=(Lx-1) && i!=N-1))
    {
        *v1= i+1;
        *v2= i-1;
        *v3= i+Lx;
        *v4= i-Lx;
        *v5=i-(Lx-1);
        *v6=i+(Lx-1);
    }
    
    //CASO 2b: Célula en el borde inferior columna par (exceptuando esquinas)
    else if ((ly==0 && i!=0 && i!=Lx-1 && (int)(lx%2)==0))
    {
        *v1= i+1;
        *v2= i-1;
        *v3= i+Lx;
        *v4= i+(N-Lx);
        *v5=i+(N-Lx+1);
        *v6=i+(N-Lx-1);
    }
    
    //CASO 2c: Célula en el borde inferior columna impar (exceptuando esquinas)
    else if ((ly==0 && i!=0 && i!=Lx-1 && (int)(lx%2)!=0))
    {
        *v1= i+1;
        *v2= i-1;
        *v3= i+Lx;
        *v4= i+(N-Lx);
        *v5=i+(Lx+1);
        *v6=i+(Lx-1);
    }

    //CASO 2d: Célula en el borde superior columna par (exceptuando esquinas)
    else if ((ly==(Ly-1) && i!=(N-Ly) && i!=(N-1) && (int)(lx%2)==0))
    {
        *v1= i+1;
        *v2= i-1;
        *v3= i-Lx;
        *v4= i-(N-Lx);
        *v5=i-(Lx+1);
        *v6=i-(Lx-1);
    }
    
    //CASO 2e: Célula en el borde superior columna impar (exceptuando esquinas)
    else if ((ly==(Ly-1) && i!=(N-Ly) && i!=(N-1) && (int)(lx%2)!=0)){
        *v1= i+1;
        *v2= i-1;
        *v3= i-Lx;
        *v4= i-(N-Lx);
        *v5=i-(N-Lx+1);
        *v6=i-(N-Lx-1);
    }

    //CASO 3: CÉLULAS EN ESQUINAS
    else if (i==0) // Célula en la esquina inferior izq.
    {
        *v1= i+1;
        *v2= i+(N-1);
        *v3= i+Lx;
        *v4= i+(N-Lx);
        *v5=i+(Lx-1);
        *v6=i+(N-Lx+1);
    }
    else if ( i==(Lx-1)) // Célula en la esquina inferior dcha.
    {
        *v1= i+1;
        *v2= i-1;
        *v3= i+Lx;
        *v4= i+(Lx-1);
        *v5=i-(Lx-1);
        *v6=i+(N-Lx);
    }
    else if (i==(N-Lx)) // Célula en la esquina superior izda.
    {
        *v1= i+1;
        *v2= i-1;
        *v3= i-Lx;
        *v4= i-(Lx-1);
        *v5=i+(Lx-1);
        *v6=i-(N-Lx);
    }
    else if (i==(N-1)) // Célula en la esquina superior dcha.
    {
        *v1= i-1;
        *v2= i-(N-Lx);
        *v3= i-Lx;
        *v4= i-(N-1);
        *v5=i-(Lx-1);
        *v6=i-(N-Lx+1);
    }
}

/** FUNCIÓN A RESOLVER:
    f(t, y(t), y(t-tau))
    *X: Array que contiene la evolucion temporal de la funcion incognita, X[time][cell]
    i: iteración de RK2 (sirve para acceder a la solución en el tiempo desfasado tau
**/
void f(double output[N_AGENTS], double x[][N_AGENTS], double *y, int i, double *tau, double tau_c, double nHill, double K, double alpha, double Km, double V, double gam, double eps, double nHill_c, double K_c)
{
    int v[6]={0}; // Array para almacenar posiciones de primeros vecinos
    int index_delayed=0, index_couple=0; // Posiciones en el array temporal correspondientes a los delays
    double x_ext = 0; // Concentración promedio de los primeros vecinos
    
    for(int j=0; j<N_AGENTS; j++){
        CondContorno(&v[1],&v[2],&v[3],&v[4],&v[5],&v[6], j); // Calculamos los vecinos próximos del agente j
        index_delayed = i-round(tau[j]/h);
        index_couple = i-round(tau_c/h);
        for(int k=0; k<6; k++)
            x_ext += x[index_couple][v[k]];
        x_ext = x_ext/6.;
        output[j] = alpha/(1+pow(x[index_delayed][j]/K,nHill))*(1-eps+eps/K_c/(1+pow(x_ext/K_c,nHill_c)))-V*y[j]/(Km+y[j])-gam*y[j];
        x_ext = 0;
    }
}

/** RK
    Algoritmo Runge-Kutta de orden 2 para sistemas generales de N ecuaciones con N incognitas
    X[time][cell]:     Array que contiene la evolucion temporal de la funcion incognita
**/
void RK2(double *t, double X[][N_AGENTS], double *tau, double tau_c, double tau_max, double nHill, double K, double alpha, double Km, double V, double gam, double eps, double nHill_c, double K_c)
{
    int n_tau = tau_max/h; //Numero de puntos almacenados con las condiciones iniciales
   
    double k1[N_AGENTS];
    double k2[N_AGENTS];
    double kAux[N_AGENTS];
    
    int i;
    for (i=n_tau; i<2*n_tau; i++) // Usar el punto i-esimo para calcular el (i+1)-esimo
    {
        f(k1, X, X[i], i, tau, tau_c, nHill, K, alpha, Km, V, gam, eps, nHill_c, K_c); // k1 = f(t_i,y_i)
        
        for (int j=0; j<N_AGENTS; j++){
            kAux[j] = X[i][j] + k1[j]*h;
        }
        f(k2, X, kAux, i+1, tau, tau_c, nHill, K, alpha, Km, V, gam, eps, nHill_c, K_c); // k2 = f(t_i+h,y_i+h*k1)
        
        for (int j=0; j<N_AGENTS; j++){
            X[i+1][j] = X[i][j] + (k1[j]+k2[j])*h*0.5;
        }
    }
}

int main()
{
    IniRapu();
    
    // PARÁMETROS ECUACIÓN CÉLULA INDIVIDUAL
    double nHill = 2.0;
    double K = 0.15;
    double alpha = 1.0;
    double Km = 50;
    double V = 1.0;
    double gam = 1.0; // Normalización implícita respecto al eje temporal
    
    // Asignación de los delays de cada célula
    double* tau = malloc(N_AGENTS * sizeof(double));
    if (tau == NULL) {
        printf("Memory allocation failed. Exiting...\n");
    }
    for (int j=0; j<N_AGENTS; j++){
        tau[j] = round((3.0 + 1.0*RapuanoRandom())*100.0)/100.0; // precisión máxima de cada delay: 2 decimales
    }
    
    // PARÁMETROS INTERACCIÓN
    double eps = 1.0; // Interacción de la interacción (entre 0 y 1)
    double tau_c = 7.00; // Delay intercelular (siempre mayor que cualquier otro tau)
    double nHill_c = 1.0;
    double K_c = 0.15;
    
    //Búsqueda del tau máximo
    double tau_max = tau_c;
    for (int i=0; i<N_AGENTS; i++){
        if (tau[i]>tau_max)
            tau_max = tau[i];
    }
    
    // Variables para integración
    int n_tau = tau_max/h; // Numero de puntos almacenados con las condiciones iniciales
    int iter_max = ceil(tf/tau_max); // Numero de iteraciones de RK2 (por intervalos de tiempo de longitud tau_max) necesarias para resolver el sistema
    printf("\n\niter_max = %d", iter_max);
    printf("\n2*ntau = %d", 2*n_tau);
    
    double t[2*n_tau];
    double x_cells[2*n_tau+1][N_AGENTS];
    
    //INICIALIZACIÓN DEL SISTEMA: entre -tau y 0, cada célula toma valor constante x0 asignado aleatoriamente
    double aux=0;
    for (int j=0; j<N_AGENTS; j++){
        aux = 0.7*RapuanoRandom(); // condiciones iniciales aleatorias
        for (int i=0; i<=n_tau; i++){
            x_cells[i][j]=aux;
        }
    }
    
    //FICHERO con EVOLUCIÓN TEMPORAL de CONCENTRACIONES
    FILE *g;
    g=fopen("evol_temp.txt","w");
    fprintf(g,"#VARIABLES DEL SISTEMA:\n");
    fprintf(g,"#nHill = %.3f\t alpha = %.3f\t K = %.3f\t V = %.3f\t Km = %.3f\t gamma = %.3f\t\n", nHill, alpha, K, V, Km, gam);
    fprintf(g,"#t\t");
    for(int i=0; i<N_AGENTS; i++){
        fprintf(g,"x[%d]\t",i);
    }
    fprintf(g,"\n");
    
    // Pasar concentraciones de condiciones iniciales a fichero
    for(int i=0; i<n_tau; i++){
        t[i]=i*h-tau_c;
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
    fprintf(g_order,"#nHill = %.3f\t alpha = %.3f\t K = %.3f\t V = %.3f\t Km = %.3f\t gamma = %.3f\t\n", nHill, alpha, K, V, Km, gam);
    fprintf(g_order,"#t_order\tR\n");
    
    int iter_term = ceil(0.1*tf/tau_c); // TERMALIZACIÓN: iteración del total de simulacion a partir de la cual calculamos R
    
    // VARIABLES PARÁMETRO DE ORDEN
    double numt=0, denomt=0;
    double M=0, M_sumt=0, M2_sumt=0; // En M_sumt almacenamos la suma de M para el promedio temporal. Íden con M2_sumt y la suma de cuadrados.
    double x_sumt[N_AGENTS], x2_sumt[N_AGENTS]; // En x_sumt almacenamos la suma de x para el promedio temporal. Íden con x2_sumt y la suma de cuadrados.
    for(int j=0; j<N_AGENTS; j++){
        x_sumt[j] = 0;
        x2_sumt[j] = 0;
    }
    
    //------------- CÁLCULO DE LA SOLUCIÓN NUMÉRICA ------------- START
    
    //INTEGRACIÓN DE LA ED POR INTERVALOS DE TIEMPO DE LONGITUD tau
    int steps_t=1; // Número de puntos usados para calcular el promedio temporal (ontador de iteraciones global)
    
    for (int cont_iter=0; cont_iter<iter_max; cont_iter++)
    {
        // Actualizar array de tiempo en la iteracion cont_iter
        for (int i=0; i<=2*n_tau; i++){
            t[i]=cont_iter*tau_c+(i*h-tau_c);
        }
        
        // Ejecutar RK2 para cada tramo para calcular x_cells en [cont_iter*tau, (cont_iter+1)*tau]
        RK2(t, x_cells, tau, tau_c, tau_max, nHill, K, alpha, Km, V, gam, eps, nHill_c, K_c);
        
        // Pasar concentraciones a fichero
        for(int i=n_tau; i<=2*n_tau; i++){
            fprintf(g,"%.3f",t[i]);
            for(int j=0; j<N_AGENTS; j++){
                //int index = round(i-tau[j]/h);
                fprintf(g,"\t%.3f",x_cells[i][j]);
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
                    printf("t = %.3f, R = %.3f\n", t[i], numt/denomt);
                }
                
                M = 0;
                denomt = 0;
                
                steps_t++;
            }
        }
        
        // Actualizar array de concentraciones
        for(int i=0; i<=n_tau; i++){
            for(int j=0; j<N_AGENTS; j++){
                x_cells[i][j] = x_cells[i+n_tau][j];
            }
        }
    }
    
    //------------- CÁLCULO DE LA SOLUCIÓN NUMÉRICA ------------- END
    
    printf("All good.\n");
    
    fclose(g);
    fclose(g_order);

    return 0;
}


