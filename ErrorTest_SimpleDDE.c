//
//  main.c
//  RK2 (TEST)
//
//  Created by Pilar Gil
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <motivation.h>

// PARÁMETROS GLOBALES:
double tf = 60; // tomamos como tiempo inicial t0 = 0 (a partir de 40 peta)
double h = 0.01;

/** FUNCIÓN A RESOLVER:
    f(t, y(t), y(t-tau))
    *X: Array que contiene la evolucion temporal de la funcion incognita
    i: iteración de RK2 (sirve para acceder a la solución en el tiempo desfasado tau
**/
double f_raw(double *X, int i, double y, double t, double tau)
{
    int index_delayed = t/h;
    return -X[index_delayed];
}

double f_midpoint(double *X, int i, double y, double t, double tau)
{
    double index_delayed = t/h;
    printf("index: %f\n", index_delayed);
    int n1 = floor(index_delayed);
    int n2 = ceil(index_delayed);
    return (-X[n1]-X[n2])/2;
}

double f_euler(double *X, double y, double t)
{
    int index_delayed = t/h;
    return -X[index_delayed]; // ECUACIÓN SENCILLA
}

/** RK
    Algoritmo Runge-Kutta de orden 2 para sistemas de una dimensión
    *f:     Puntero a la funcion vectorial
    t0:     Tiempo inicial
    tf:      Tiempo final
    x0:    Condición inicial
    N:     Dimension del sistema (numero de ecuaciones y de incognitas); dimension de los arrays x0, t, x, y
    X:     Array que contiene la evolucion temporal de la funcion incognita
**/
void RK2_heun(double tf, double tau, double t[], double X[])
{
    int n = tf/h;

    double k1;
    double k2;

    int i;

    for (i=tau/h-1; i<n+tau/h+1; i++) // Usar el punto i-esimo para calcular el (i+1)-esimo
    {
        t[i] = (i-tau/h)*h;
        
        k1 = f_raw(X,i,X[i],t[i],tau);
        k2 = f_raw(X,i,X[i]+h*k1,t[i]+h,tau);
        
        X[i+1] = X[i] + (k1+k2)*h/2.0;
        
        printf("X[i]: %.3f\n", X[i]);
        printf("\n");
    }
}

void RK2_midpoint_raw(double tf, double tau, double t[], double X[])
{
    int n = tf/h;

    double k1;
    double k2;

    int i;

    for (i=tau/h-1; i<n+tau/h+1; i++) // Usar el punto i-esimo para calcular el (i+1)-esimo
    {
        t[i] = (i-tau/h)*h;
        
        k1 = f_raw(X,i,X[i],t[i],tau);
        k2 = f_raw(X,i,X[i]+h*k1/2,t[i]+h/2,tau);
        
        X[i+1] = X[i] + k2*h;
        
        printf("X[i]: %.3f\n", X[i]);
        printf("\n");
    }
}

void RK2_midpoint_interpol(double tf, double tau, double t[], double X[])
{
    int n = tf/h;

    double k1;
    double k2;

    int i;

    for (i=tau/h-1; i<n+tau/h+1; i++) // Usar el punto i-esimo para calcular el (i+1)-esimo
    {
        t[i] = (i-tau/h)*h;
        
        k1 = f_midpoint(X,i,X[i],t[i],tau);
        k2 = f_midpoint(X,i,X[i]+h*k1/2,t[i]+h/2,tau);
        
        X[i+1] = X[i] + k2*h;
        
        printf("X[i]: %.3f\n", X[i]);
        printf("\n");
    }
}

void EULER(double tf, double tau, double *t, double X[])
{
    int n = tf/h;
    int i;

    for (i=tau/h; i<n+tau/h; i++) // Usar el punto i-esimo para calcular el (i+1)-esimo
    {
        X[i+1] = X[i] + f_euler(X, X[i], t[i])*h;
    }
}

int factorial(int k){
    double fact=1;
    if (k==0)
        return fact;
    else for (int i=1; i<=k; i++)
        fact = fact*i;
    return fact;
}


int main()
{
    double tau = 1.5;
    int n_delay = tau/h; //Numero de puntos almacenados con las condiciones iniciales
    int n_interval = tf/h; //Numero de puntos a calcular sin contar con las condiciones iniciales
    int n_total = n_interval + n_delay; //Numero total de puntos que almacenamos
    
    // CÁLCULO DE LA SOLUCIÓN NUMÉRICA
    double t[n_total+1];
    double X_raw[n_total+1], X_midpoint[n_total+1], X_heun[n_total+1], X_euler[n_total+1];
    printf("%.3d \n",n_total);
    
    // INICIALIZACIÓN DEL SISTEMA: entre -tau y 0, X toma valor constante x0=1
    for (int i=0; i<n_delay; i++){
        X_raw[i]=1.0;
        X_midpoint[i]=1.0;
        X_heun[i]=1.0;
        X_euler[i]=1.0;
        t[i]=-(n_delay-i)*h;
        printf("t[i]: %.3f\n", t[i]);
    }
    
    RK2_heun(tf, tau, t, X_heun);
    RK2_midpoint_raw(tf, tau, t, X_raw);
    RK2_midpoint_interpol(tf, tau, t, X_midpoint);
    EULER(tf, tau, t, X_euler);

    printf("All good.\n");
    
    // CÁLCULO DE LA SOLUCIÓN ANALÍTICA
    double X_an[n_total+1];
    double time, aux_pownum=0, aux_sumdenom=0;
    for (int i=0; i<n_delay; i++){
        X_an[i]=1;
    }
    
    int n = 1;
    for (int i=tau/h; i<n_total+1; i++){
        time = (i-tau/h)*h;
        t[i] = time;
        X_an[i]=1.0;
        if (time > n*tau)
            n++;
        for (int k=1; k<=n; k++){
            
            // TOMAMOS LOGARITMOS SOBRE EL COCIENTE
            aux_pownum = k*log(time-(k-1)*tau);
            for (int j=1; j<k+1; j++){
                aux_sumdenom = aux_sumdenom + log(j);
            }
            X_an[i]+= 1.0*pow(-1,k)*exp(aux_pownum-aux_sumdenom);
            aux_sumdenom = 0;
        }
    }
    
    // ALMACENAMIENTO EN FICHERO
    FILE *g;
    g=fopen("evol_temp.txt","w");
    fprintf(g,"#t\tx_mumHeun\tx_numMidpointRaw\tx_numMidPointInterpol\tx_Euler\tx_analitic\n");
    for(int j=0; j<n_total+1; j++){
        fprintf(g,"%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",t[j],X_heun[j],X_raw[j],X_midpoint[j],X_euler[j],X_an[j]);
    }
    fclose(g);

    return 0;
}

