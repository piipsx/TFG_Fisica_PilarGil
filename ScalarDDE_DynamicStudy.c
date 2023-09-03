//
//  main.c
//  DynamicStudy_ScalarDDE
//
//  Created by Pilar Gil
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// PARÁMETROS GLOBALES:
double tf = 200; // tomamos como tiempo inicial t0 = 0
double h = 0.005;

/** FUNCIÓN A RESOLVER:
    f(t, y(t), y(t-tau))
    *X: Array que contiene la evolucion temporal de la funcion incognita
    i: iteración de RK2 (sirve para acceder a la solución en el tiempo desfasado tau
**/
double f_raw(double *X, int i, double y, double t, double tau, double a, double b)
{
    int index_delayed = t/h;
    return (b*X[index_delayed]+a*y)/2; // SIMPLE SCALAR DDE
}

/** RK
    Algoritmo Runge-Kutta de orden 2 para sistemas de una dimensión
    tf:      Tiempo final
    tau:    Retardo del sistema
    X:     Array que contiene los tiempos de simulación
    X:     Array que contiene la evolucion temporal de la funcion incognita
**/
void RK2_midpoint(double tf, double tau, double t[], double X[], double a, double b)
{
    int n = tf/h;

    double k1;
    double k2;

    int i;

    for (i=tau/h-1; i<n+tau/h+1; i++) // Usar el punto i-esimo para calcular el (i+1)-esimo
    {
        t[i] = (i-tau/h)*h;
        
        k1 = f_raw(X,i,X[i],t[i],tau,a,b);
        k2 = f_raw(X,i,X[i]+h*k1,t[i]+h,tau,a,b);
        
        X[i+1] = X[i] + (k1+k2)*h/2.0;
    }
}

// FUNCIÓN CARACTERÍSTICA DEL SISTEMA: f =  f[0] (parte real) + i*f[1] (parte imaginaria)
void f_root(double f[2], double x[2], double a, double b, double tau){ // COMP 2 de FUNCIÓN f
    f[0] = x[0]-b*exp(-x[0]*tau)*cos(x[1]*tau)-a;
    f[1] = x[1]+b*exp(-x[0]*tau)*sin(x[1]*tau);
}

// JACOBIANO DE LA FUNCIÓN f: matriz 2x2
void Jf_root(double Jf[2][2], double x[2], double a, double b, double tau){ // DERIVADA DE f
    Jf[0][0] = 1+b*tau*exp(-x[0]*tau)*cos(x[1]*tau);
    Jf[0][1] = b*tau*exp(-x[0]*tau)*sin(x[1]*tau);
    Jf[1][0] = -b*tau*exp(-x[0]*tau)*sin(x[1]*tau);
    Jf[1][1] = 1+b*tau*exp(-x[0]*tau)*cos(x[1]*tau);
}

int main()
{
    // PARÁMETROS DE LA ECUACIÓN
    double tau = 1.0;
    double a = -4;
    double b = -4;
    
    int n_delay = tau/h; //Numero de puntos almacenados con las condiciones iniciales
    int n_interval = tf/h; //Numero de puntos a calcular sin contar con las condiciones iniciales
    int n_total = n_interval + n_delay; //Numero total de puntos que almacenamos
    
    //------------- CÁLCULO DE LA SOLUCIÓN NUMÉRICA: START -------------//
    
    double t[n_total+1];
    double X_midpoint[n_total+1];
    
    // INICIALIZACIÓN DEL SISTEMA: entre -tau y 0, X toma valor constante x0=1
    for (int i=0; i<n_delay; i++){
        X_midpoint[i]=5.0;
        t[i]=-(n_delay-i)*h;
    }
    
    RK2_midpoint(tf, tau, t, X_midpoint, a, b);

    printf("All good.\n");
    
    //------------- CÁLCULO DE LA SOLUCIÓN NUMÉRICA: END -------------//
    
    //------------- ESTUDIO DE ESTABILIDAD: START -------------//
    
    float tol = 10e-6; // tolerance
    int itmax = 10e2; // maximum number of iterations
    double reachitmax = 0;
    
    // ALMACENAMIENTO EN FICHERO
    FILE *g;
    g=fopen("eigenvalues.txt","w");
    
    fprintf(g,"%d\t", 17); //necesario para la representación en gnuplot
    for(b=-4; b<=4; b=b+0.5){
        fprintf(g,"%.1f\t", b);
    }
    fprintf(g,"\n");
    
    for(a=-4; a<=4; a=a+0.5){
        fprintf(g,"%.1f\t", a);
        
        for(b=-4; b<=4; b=b+0.5){
            
            double x0[2]={0.1, 0.1}, x1[2]={0};
            double f[2] = {0}, f_test[2] = {0};
            double Jf[2][2] = {{0, 0},{0, 0}}, Jf_inv[2][2] = {{0, 0},{0, 0}};
            double detJf=0;
            
            printf("it\t x0\t f0\n");
            
            for (int i=0; i<itmax; i++){
                
                f_root(f, x0, a, b, tau);
                Jf_root(Jf, x0, a, b, tau);
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
                x1[0] = x0[0] -(Jf_inv[0][0]*f[0] + Jf_inv[0][1]*f[1]);
                x1[1] = x0[1] -(Jf_inv[1][0]*f[0] + Jf_inv[1][1]*f[1]);
                
                printf("%d\t%f + i%f\t%f + i%f\n",i,x0[0],x0[1],f[0],f[1]);
                x0[0] = x1[0]; x0[1] = x1[1];
                
                f_root(f_test, x1, a, b, tau);
                if (sqrt(pow(f[0]-f_test[0],2)+pow(f[1]-f_test[1],2)) < tol){
                    printf("Solution reached for a = %.1f, b = %.1f.\n\n", a, b);
                    break;
                }
                
                if (i == itmax-1)
                    reachitmax = 1;
            }
            
            if (a+b>0)
                printf("\nx=0 is unstable.\n");
            else if (b>=a)
                printf("\nx=0 is asymptotically stable.\n");
            else printf("\nThere exists critical delay below which x=0 is asymptotically stable.\n");
            
            printf("Root from characteristic eq. is: %f + i%f\n\n", x1[0],x1[1]);
            
            //------------- ESTUDIO DE ESTABILIDAD: END -------------//
            
            if (reachitmax == 1){
                fprintf(g,"nan\t");
                reachitmax = 0;
            }
            else fprintf(g,"%.3f\t", x1[0]);
            
        }
        fprintf(g,"\n");
    }
    
    fclose(g);

    return 0;
}

