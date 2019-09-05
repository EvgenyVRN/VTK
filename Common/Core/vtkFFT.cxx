#include "vtkFFT.h"

#include "vtkMath.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkFFT);

namespace
{

int fft_frame_size_bits(int dataSize)
{
  int k = 0; // Длина n в битах
  while ((1 << k) < dataSize)
    k++;
  return k;
}

int fft_frame_size(int dataSize)
{
  return 1 << fft_frame_size_bits(dataSize);
}

void make_frame(ComplexNumber* __restrict out, const double* __restrict in, const int inCount)
{
  int n = fft_frame_size(inCount);

  for (int i = 0; i < inCount; i++)
  {
    out[i].real((FFT_UNIT_TYPE)in[i]);
    out[i].imag(0.0);
  }

  if (inCount < n)
    memset(out + inCount, 0, sizeof(ComplexNumber) * (n - inCount));
}

ComplexNumber* prepareComplexArray(const double* in, const int inCount, int* outCount)
{
  int n = fft_frame_size(inCount);
  ComplexNumber* tmp = new ComplexNumber[n];
  make_frame(tmp, in, inCount);
  *outCount = n;
  return tmp;
}

ComplexNumber* prepareComplexArray(const ComplexNumber* in, const int inCount, int* outCount)
{
  int n = fft_frame_size(inCount);

  ComplexNumber* tmp = new ComplexNumber[n];
  for (int i = 0; i < inCount; i++)
    tmp[i] = in[i];

  for (int i = inCount; i < n; i++)
    tmp[i] = { 0, 0 };

  *outCount = n;
  return tmp;
}

void fft_in_out_perm(int* perm, int k)
{
  int high1 = -1, n = 1 << k;

  perm[0] = 0;
  for (int i = 1; i < n; i++)
  {
    if ((i & (i - 1)) ==
      0) // Проверка на степень двойки. Если i ей является, то i-1 будет состоять из кучи единиц.
      high1++;
    perm[i] = perm[i ^ (1 << high1)];  // Переворачиваем остаток
    perm[i] |= (1 << (k - high1 - 1)); // Добавляем старший бит
  }
}

void fft_roots(ComplexNumber* roots, int n)
{
  for (int i = 0; i < n / 2; i++)
  {
    double alpha = 2 * vtkMath::Pi() * i / n;
    roots[i].real((FFT_UNIT_TYPE)cos(alpha));
    roots[i].imag((FFT_UNIT_TYPE)sin(alpha));
  }
}

void fft_core(ComplexNumber* __restrict out, const ComplexNumber* __restrict in,
  const ComplexNumber* __restrict roots, const int* __restrict rev, int n)
{
  for (int i = 0; i < n; i++)
    out[i] = in[rev[i]];

  for (int len = 1; len < n; len <<= 1)
  {
    int rstep = n / (len * 2);
    for (int pdest = 0; pdest < n; pdest += len)
    {
      const ComplexNumber* __restrict r = roots;
      for (int i = 0; i < len; i++, pdest++, r += rstep)
      {
        ComplexNumber* __restrict a = out + pdest;
        ComplexNumber* __restrict b = a + len;

        FFT_UNIT_TYPE real = r->real() * b->real() - r->imag() * b->imag();
        FFT_UNIT_TYPE imag = r->imag() * b->real() + r->real() * b->imag();

        b->real(a->real() - real);
        b->imag(a->imag() - imag);
        a->real(a->real() + real);
        a->imag(a->imag() + imag);
      }
    }
  }
}

void fft(const ComplexNumber* in, const int inCount, ComplexNumber*& outData)
{
  int k = fft_frame_size_bits(inCount);
  int n = 1 << k;

  int* rev = new int[n];
  ComplexNumber* roots = new ComplexNumber[n / 2];

  fft_in_out_perm(rev, k);
  fft_roots(roots, n);

  if (outData == 0x0)
    outData = new ComplexNumber[n];

  fft_core(outData, in, roots, rev, n);

  delete[] roots;
  delete[] rev;
}

void fft_post_inverse(ComplexNumber* data, int n)
{
  FFT_UNIT_TYPE tmpCoef = 1.0 / n;
  for (int i = 0; i < n; i++)
    data[i] *= tmpCoef;

  for (int i = 1; i < ((n - 1) / 2 + 1); i++)
    std::swap(data[i], data[n - i]);
}
} // end namespace

void vtkFFT::fft_direct(const double* in, const int inCount, int* outCount, ComplexNumber*& outData)
{
  ComplexNumber* tmp = prepareComplexArray(in, inCount, outCount);
  fft(tmp, *outCount, outData);
  delete[] tmp;
}

void vtkFFT::fft_inverse(
  const ComplexNumber* in, const int inCount, int* outCount, ComplexNumber*& outData)
{
  ComplexNumber* tmp = prepareComplexArray(in, inCount, outCount);
  fft(tmp, *outCount, outData);
  fft_post_inverse(outData, fft_frame_size(inCount));
  delete[] tmp;
}

void vtkFFT::complexes_to_doubles(
  double* __restrict out, const ComplexNumber* __restrict in, const int inCount)
{
  for (int i = 0; i < inCount; i++)
  {
    out[i] = in[i].real();
  }
}

double vtkFFT::complex_module(const ComplexNumber& in)
{
  return sqrt(in.real() * in.real() + in.imag() * in.imag());
}

std::vector<double> vtkFFT::rfftfreq(int windowLength, double sampleSpacing)
{
  std::vector<double> res;
  double val = 1.0 / (windowLength * sampleSpacing);
  int N = (windowLength + 1) / 2;
  for (int i = 0; i < N; i++)
    res.push_back(i * val);

  return res;
}

void vtkFFT::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
