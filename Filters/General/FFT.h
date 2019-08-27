#ifndef FFT_H
#define FFT_H

#define FFT_UNIT_TYPE float

#include <math.h>
#include <complex>
#include <vector>

typedef std::complex<FFT_UNIT_TYPE> ComplexNumber;

void fft_direct(const float* in, const int inCount, int* outCount, ComplexNumber*& outData);
void fft_inverse(const float* in, const int inCount, int* outCount, ComplexNumber*& outData);

void fft_direct(const double* in, const int inCount, int* outCount, ComplexNumber*& outData);
void fft_inverse(const double* in, const int inCount, int* outCount, ComplexNumber*& outData);

void fft_direct(const ComplexNumber* in, const int inCount, int* outCount, ComplexNumber*& outData);
void fft_inverse(
  const ComplexNumber* in, const int inCount, int* outCount, ComplexNumber*& outData);

int fft_frame_size_bits(int dataSize);
int fft_frame_size(int dataSize);
void fft_in_out_perm(int* perm, int k);
void fft_roots(ComplexNumber* roots, int n);
void fft_core(ComplexNumber* __restrict out, const ComplexNumber* __restrict in,
  const ComplexNumber* __restrict roots, const int* __restrict rev, int n);
void fft_post_inverse(ComplexNumber* data, int n);

void make_frame(ComplexNumber* __restrict out, const float* __restrict in, const int inCount);
void make_frame(ComplexNumber* __restrict out, const double* __restrict in, const int inCount);
void complexes_to_doubles(
  double* __restrict out, const ComplexNumber* __restrict in, const int inCount);
double* complexes_to_doubles(const ComplexNumber* in, const int inCount);

double complex_module(const ComplexNumber &in);

std::vector<double> rfftfreq(int windowLength, double sampleSpacing);

#endif /* FFT_H */
