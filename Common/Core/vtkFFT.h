/**
 * @class vtkFFT
 * @brief perform Discrete Fourier Transforms
 *
 * vtkFFT provides methods to perform Discrete Fourier Transforms.
 * These include providing ...
 */

#ifndef vtkFFT_h
#define vtkFFT_h

#include "vtkObject.h"

#include <complex>
#include <vector>


using FFT_UNIT_TYPE = float;
using ComplexNumber = std::complex<FFT_UNIT_TYPE>;

class VTKCOMMONCORE_EXPORT vtkFFT : public vtkObject
{
public:
  static vtkFFT *New();
  vtkTypeMacro(vtkFFT,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  /**
   * Compute the one-dimensional discrete Fourier Transform for real input
   */
  static void fft_direct(const double* in, const int inCount, int* outCount, ComplexNumber*& outData);

  /**
   * Compute the inverse of DFT
   */
  static void fft_inverse(
    const ComplexNumber* in, const int inCount, int* outCount, ComplexNumber*& outData);

  /**
   * Convert an array of complex numbers to array of doubles
   */
  static void complexes_to_doubles(double* __restrict out, const ComplexNumber* __restrict in, const int inCount);

  /**
   * Return the module of complex number
   */
  static double complex_module(const ComplexNumber &in);

  /**
   * Return the Discrete Fourier Transform sample frequencies
   */
  static std::vector<double> rfftfreq(int windowLength, double sampleSpacing);

protected:
  vtkFFT() {}
  ~vtkFFT() override {}
private:
  vtkFFT(const vtkFFT&) = delete;
  void operator=(const vtkFFT&) = delete;
};

#endif
