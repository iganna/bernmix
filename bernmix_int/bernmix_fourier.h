//
//  bernmix_fourier.h
//  
//
//  Created by Anna Igolkina on 30/09/2017.
//
//

#ifndef bernmix_fourier_h
#define bernmix_fourier_h

#define _USE_MATH_DEFINES
#define EPS 1e-16

#include <stdlib.h>
#include <math.h>
#include <fftw3.h>


void pmf_bern_mixture(int n, double* P, int* W, double* res);


#endif /* bernmix_fourier_h */
