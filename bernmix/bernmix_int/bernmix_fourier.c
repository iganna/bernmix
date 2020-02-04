#include "bernmix_fourier.h"

/*!
 * This function equates to zero the value less than computation error epsilon
 *
 * @param x a number to be corrected
 * @param eps the epsilon error
 *
 * @return corrected value
 */
static double handler(double x, double eps) {

    if (fabs(x) < eps)
        return 0.0;
    else
        return x;
}

/*!
 * A piecewise-defined function required for computations
 *
 * @param y real part of a complex number
 * @param x imaginary part of a complex number
 *
 * @return a computed value
 */
static double arctan2(double y, double x) {

    if ((x == 0.) && (y == 0.))
        return 0.;
    else if ((x == 0.) && (y > 0.))
        return M_PI / 2;
    else if ((x == 0.) && (y < 0.))
        return -M_PI / 2;
    else if ((x < 0.) && (y >= 0.))
        return atan(y / x) + M_PI;
    else if ((x < 0.) && (y < 0.))
        return atan(y / x) - M_PI;
    else
        return atan(y / x);
}

/*!
 * A real part of z value
 */
static double z_real(int l, double p, int w, int sw) {

    return (1. - p + p * cos((double)w * 2. * M_PI / ((double)sw + 1.) * (double)l));
}

/*!
 * An imaginary part of z value
 */
static double z_imag(int l, double p, int w, int sw) {

    return (p * sin((double)w * 2. * M_PI / ((double)sw + 1.) * (double)l));
}

/*!
 * A modulus of z complex value
 */
static double z_modulus(double real, double imag) {

    return sqrt(real * real + imag * imag);
}

/*!
 * An argument of z complex value
 */
static double z_argument(double real, double imag) {

    return arctan2(imag, real);
}

/*!
 * A piecewise-defined function required for computations
 *
 * @param y real part of a complex number
 * @param x imaginary part of a complex number
 *
 * @param res - reference to two returning values
 */
static void comp_numb(int l, double* P, int* W, int sw, int n, double* res) {

    double slg = 0.;
    double sar = 0.;
    double r_tmp = 0.;
    double i_tmp = 0.;
    int j = 0;

    for (; j < n; j++) {
        r_tmp = z_real(l, P[j], W[j], sw);
        i_tmp = z_imag(l, P[j], W[j], sw);
        slg += log(z_modulus(r_tmp, i_tmp));
        sar += z_argument(r_tmp, i_tmp);
    }

    res[0] = exp(slg) * cos(sar);
    res[1] = exp(slg) * sin(sar);
}

/*!
 * Fourier transform
 */
static void fft(int n, double* in_r, double* in_i, double* out_r) {

    fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);

    int i = 0;
    for (; i < n; i++) {
        in[i][0] = in_r[i];
        in[i][1] = in_i[i];
    }

    fftw_plan plan = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    for (i = 0; i < n; i++) {
        out_r[i] = out[i][0];
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
}

/*!
 * Computation of probability for integer values of the weighted sum
 */
void pmf_bern_mixture(int n, double* P, int* W, double* res)
{
    double* in_r;
    double* in_i;
    double* temp = NULL;
    int border;
    int l, i;
    int sw = 0;

    // Compute the sum of
    for (i = 0; i < n; i++)
    {
        sw += W[i];
    }
    border = (int)ceil((double)sw / 2) + 1;
    
    in_r = (double*)malloc(sizeof(double) * (sw + 1));
    in_i = (double*)malloc(sizeof(double) * (sw + 1));
    in_r[0] = 1.0;
    in_i[0] = 0.0;

    for (l = 1; l < sw + 1; l++) {
        if (l < border) {
            temp = (double*)malloc(sizeof(double) * 2);
            comp_numb(l, P, W, sw, n, res);
            in_r[l] = res[0];
            in_i[l] = res[1];
            free(temp);
        }
        else {
            in_r[l] = in_r[sw + 1 - l];
            in_i[l] = -in_i[sw + 1 - l];
        }
    }

    fft(sw + 1, in_r, in_i, res);

    for (l = 0; l < sw + 1; l++) {
        res[l] = handler(res[l], EPS) / (sw + 1);
    }

    free(in_r);
    free(in_i);
}

