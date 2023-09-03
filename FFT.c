//
//  main.c
//  FAST FOURIER TRANSFORM
//
//  Created by Pilar Gil
//

/* PROPIEDADES DEL PROGRAMA: Cálculo de la Transformada de Fourier Discreta de una serie de datos reales (RDFT)
 
Complex structure: complex number with real and imaginary components.
cooleyTukeyFFT function: recursive function that performs the Cooley-Tukey FFT algorithm. INPUT: data array, array of Complex numbers to store the FFT results, size of the input data, and current stride as arguments.
computeFFT function: it initializes the complex array, calls the cooleyTukeyFFT function, and then prints the FFT results.

 **Note that this implementation does not assume that the input sequence length is a power of 2 for simplicity. It performs zero-padding to the next power of 2 and initializes the remaining elements with zeros.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    double real;
    double imag;
} Complex;

double gaussian(double x, double amplitude, double mean, double dispersion) {
    double exponent = -0.5 * pow((x - mean) / dispersion, 2);
    double coefficient = 1.0 / (dispersion * sqrt(2 * M_PI));
    return amplitude * coefficient * exp(exponent);
}

Complex complexAdd(Complex a, Complex b) {
    Complex result;
    result.real = a.real + b.real;
    result.imag = a.imag + b.imag;
    return result;
}

Complex complexSubtract(Complex a, Complex b) {
    Complex result;
    result.real = a.real - b.real;
    result.imag = a.imag - b.imag;
    return result;
}

Complex complexMultiply(Complex a, Complex b) {
    Complex result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}

void ditfft2(Complex x[], Complex X[], int n, int s) {
    if (n == 1) {
        X[0] = x[0];  // Trivial size-1 DFT base case
    }
    else {
        Complex* X0 = (Complex*)malloc((n/2) * sizeof(Complex));
        Complex* X1 = (Complex*)malloc((n/2) * sizeof(Complex));

        ditfft2(x, X0, n/2, 2*s);                // DFT of (x0, x2s, x4s, ..., x(N-2)s)
        ditfft2(x+s, X1, n/2, 2*s);              // DFT of (xs, xs+2s, xs+4s, ..., x(N-1)s)

        for (int k = 0; k < n/2; k++) {
            Complex p = X0[k];
            Complex q;
            q.real = cos(-2 * M_PI * k / n);
            q.imag = sin(-2 * M_PI * k / n) * (-1);
            q = complexMultiply(q, X1[k]);

            X[k] = complexAdd(p, q);
            X[k + n/2] = complexSubtract(p, q);
        }

        free(X0);
        free(X1);
    }
}

int main() {
    int n = 2048;
    Complex x[n]; // = {{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}};
    Complex X[n];
    
    /*
    // Generate the function test data
    double amplitude = 1.0;   // Amplitude of the function
    double mean = 30.0;
    double dispersion = 5.0;
    double deltaT = 1.0 / n; // Time interval between sampling (normally in the denominator is written an n
    double fs = 1.0/deltaT;  // Sampling rate in Hz
    
    // Save the sample function to a text file
    FILE *file_samplefunc;
    file_samplefunc=fopen("samplefunc.txt","w");
    
    for (int i = 0; i < n; i++) {
        double t = i * deltaT;
        x[i].imag = 0.0; // all of the input sequences are real-valued
    
        x[i].real = amplitude * sin(2 * M_PI * 100 * t) + 4 * amplitude * sin(2 * M_PI * 20 * t); // HARMONIC FUNCTION
        // x[i].real = amplitude*exp(-(pow((t - mean), 2) / (2 * pow(dispersion, 2)))); // GAUSSIAN FUNCTION
        // CODE FOR SQUARE FUNCTION with a certain frequency
        double freq = 5.0;
        double T = 1/freq;
        if (fmod(t,T) < T/2)
            x[i].real = 1.0;
        else x[i].real = -1.0;
        fprintf(file_samplefunc, "%.4f\t %.3f\n", t, x[i].real);
    }
     */

    ditfft2(x, X, n, 1);
    
    printf("Frequency Spectrum:\n");
    printf("k\t Re{Ak}\t Im{Ak}\t |Ak|^2\n");
    
    // Save the FFT results to a text file
    FILE *file_fft;
    file_fft=fopen("fft.txt","w");
    
     
        for (int i = 0; i < n; i++) {
            double magnitude = sqrt(X[i].real * X[i].real + X[i].imag * X[i].imag);
            
            
            printf("%.4f\t %.2f\t %.2f\t %.2f\n", i/(double)n*fs, X[i].real, X[i].imag, magnitude);
            fprintf(file_fft, "%.4f\t %.2f\t %.2f\t %.2f\n", i/(double)n*fs, X[i].real, X[i].imag, magnitude);
        }
     
        fclose(file_fft);
        fclose(file_samplefunc);

    return 0;
}
