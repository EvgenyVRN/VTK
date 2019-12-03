/*=========================================================================

  Program:   Visualization Toolkit
  Module:    UnitTestMath.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include <vtkFFT.h>

#include <vtkMathUtilities.h>

#include <algorithm>
#include <memory>


static bool FuzzyCompare(const ComplexNumber& result, const ComplexNumber& test, double epsilon)
{
  return (vtkFFT::complex_module(result - test) < epsilon*epsilon);
}


static int Test_fft_direct(); // TODO: not implemented yet
static int Test_fft_inverse();
static int Test_complexes_to_doubles();
static int Test_complex_module();
static int Test_rfftfreq();
static int Test_fft_direct_inverse();

int UnitTestFFT(int, char* [])
{
  int status = 0;

  status += Test_fft_direct();
  status += Test_fft_inverse();
  status += Test_complexes_to_doubles();
  status += Test_complex_module();
  status += Test_rfftfreq();
  status += Test_fft_direct_inverse();

  if (status != 0)
  {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

int Test_fft_direct()
{
// compare with numpy
//  int status = 0;
//  std::cout << "Test_fft_direct..";

//  static constexpr auto inputSize = 50;
//  static constexpr auto spectrumSize = 25;

//  // test data
//  auto i = 0;
//  std::array<ComplexNumber, spectrumSize> test;
//  // result of np.sin(np.arange(50))
//  test[i++] = ComplexNumber(1.6325e-01, +0.0000e+00);
//  test[i++] = ComplexNumber(1.6381e-01, +3.6392e-02);
//  test[i++] = ComplexNumber(1.6560e-01, +7.6177e-02);
//  test[i++] = ComplexNumber(1.6903e-01, +1.2400e-01);
//  test[i++] = ComplexNumber(1.7506e-01, +1.8809e-01);
//  test[i++] = ComplexNumber(1.8604e-01, +2.8696e-01);
//  test[i++] = ComplexNumber(2.0931e-01, +4.7599e-01);
//  test[i++] = ComplexNumber(2.8296e-01, +1.0408e+00);
//  test[i++] = ComplexNumber(-3.1623e+00, -2.4749e+01);
//  test[i++] = ComplexNumber(2.4794e-03, -1.0365e+00);
//  test[i++] = ComplexNumber(6.7456e-02, -5.3945e-01);
//  test[i++] = ComplexNumber(8.9422e-02, -3.6513e-01);
//  test[i++] = ComplexNumber(1.0032e-01, -2.7419e-01);
//  test[i++] = ComplexNumber(1.0675e-01, -2.1710e-01);
//  test[i++] = ComplexNumber(1.1093e-01, -1.7709e-01);
//  test[i++] = ComplexNumber(1.1383e-01, -1.4690e-01);
//  test[i++] = ComplexNumber(1.1593e-01, -1.2287e-01);
//  test[i++] = ComplexNumber(1.1749e-01, -1.0293e-01);
//  test[i++] = ComplexNumber(1.1867e-01, -8.5828e-02);
//  test[i++] = ComplexNumber(1.1957e-01, -7.0752e-02);
//  test[i++] = ComplexNumber(1.2026e-01, -5.7147e-02);
//  test[i++] = ComplexNumber(1.2078e-01, -4.4614e-02);
//  test[i++] = ComplexNumber(1.2116e-01, -3.2851e-02);
//  test[i++] = ComplexNumber(1.2142e-01, -2.1622e-02);
//  test[i++] = ComplexNumber(1.2157e-01, -1.0730e-02);
//  test[i++] = ComplexNumber(1.2162e-01, +0.0000e+00);

//  // input data
//  i = 0;
//  std::array<double, inputSize> input;
//  std::generate(input.begin(), input.end(), [&i]() { return std::sin(i++); });

//  // fft
//  int outSize;
//  ComplexNumber* result = nullptr;
//  vtkFFT::fft_direct(input.data(), inputSize, &outSize, result);

//  // compare
//  for(i=0; i<spectrumSize; i++)
//  {

//      if (!FuzzyCompare(result[i], test[i], 1e-03))
//      {
//        std::cout << "Expected " << test[i] << " but got " << result[i] << " difference is "
//                  << result[i] - test[i] << std::endl;
//        status++;
//      }
//  }

//  delete[] result;

//  if (status)
//  {
//    std::cout << "..FAILED" << std::endl;
//  }
//  else
//  {
//    std::cout << ".PASSED" << std::endl;
//  }
//  return status;
  return Test_fft_direct_inverse();
}

int Test_fft_inverse()
{
  return Test_fft_direct_inverse();
}

int Test_complexes_to_doubles()
{
  int status = 0;
  std::cout << "Test_complexes_to_doubles..";

  static constexpr int countOut = 6;
  static constexpr int countIn = 10;

  std::array<ComplexNumber, countIn> complexNumbers;
  auto sign = -1;
  for (auto i = 0; i < countIn; i++)
  {
    sign *= -1;
    complexNumbers[i] = ComplexNumber(i, 2 * i * sign);
  }

  std::array<double, countOut> test1 = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
  std::array<double, countOut> result1;
  vtkFFT::complexes_to_doubles(&result1[0], complexNumbers.data(), countOut);

  for (auto i = 0; i < countOut; i++)
  {
    if (!vtkMathUtilities::FuzzyCompare(
          result1.at(i), test1.at(i), std::numeric_limits<double>::epsilon()))
    {
      std::cout << "Expected " << test1.at(i) << " but got " << result1.at(i) << " difference is "
                << test1.at(i) - result1.at(i) << std::endl;
      status++;
    }
  }

  if (status)
  {
    std::cout << "..FAILED" << std::endl;
  }
  else
  {
    std::cout << ".PASSED" << std::endl;
  }
  return status;
}

int Test_complex_module()
{
  int status = 0;
  std::cout << "Test_complex_module..";

  ComplexNumber complexNumber1(3, 4);
  double module1 = vtkFFT::complex_module(complexNumber1);
  double test1 = 5;
  if (!vtkMathUtilities::FuzzyCompare(module1, test1, std::numeric_limits<double>::epsilon()))
  {
    std::cout << "Expected " << test1 << " but got " << module1 << " difference is "
              << module1 - test1 << std::endl;
    status++;
  }

  if (status)
  {
    std::cout << "..FAILED" << std::endl;
  }
  else
  {
    std::cout << ".PASSED" << std::endl;
  }
  return status;
}

int Test_rfftfreq()
{
  int status = 0;
  std::cout << "Test_rfftfreq..";

  constexpr auto samplingFrequency = 1000;
  constexpr auto windowLength = 1000;
  double sampleSpacing = 1.0 / samplingFrequency;
  std::vector<double> frequencies = vtkFFT::rfftfreq(windowLength, sampleSpacing);

  std::vector<double> test1;
  for (auto i = 0; i < windowLength / 2; i++)
    test1.push_back(i);

  if (!(frequencies.size() == test1.size()))
  {
    std::cout << "Difference size: expected " << test1.size() << " but got " << frequencies.size()
              << std::endl;
    status++;
  }

  for (auto i = 0; i < frequencies.size(); i++)
  {
    const auto& expected = test1[i];
    const auto& real = frequencies[i];

    if (!vtkMathUtilities::FuzzyCompare(real, expected, std::numeric_limits<double>::epsilon()))
    {
      std::cout << "Expected " << expected << " but got " << real << " difference is "
                << expected - real << std::endl;
      status++;
    }
  }

  //  if(!(frequencies == test1))
  //    status++;

  if (status)
  {
    std::cout << "..FAILED" << std::endl;
  }
  else
  {
    std::cout << ".PASSED" << std::endl;
  }
  return status;
}

int Test_fft_direct_inverse()
{
  int status = 0;
  std::cout << "Test_fft_direct_inverse..";

  static constexpr auto countIn = 1000;
  std::array<double, countIn> input;
  auto val = 0;
  std::generate(input.begin(), input.end(), [&val]() { return std::sin(val++); });

  int countOut = -1;
  ComplexNumber* spectrum = nullptr;
  vtkFFT::fft_direct(input.data(), countIn, &countOut, spectrum);

  int countRes = -1;
  ComplexNumber* resultComplex = nullptr;
  vtkFFT::fft_inverse(spectrum, countOut, &countRes, resultComplex);

  double result[countRes];
  vtkFFT::complexes_to_doubles(&result[0], resultComplex, countIn);

  delete[] spectrum;
  delete[] resultComplex;

  for (auto i = 0; i < countIn; i++)
  {
    if (!vtkMathUtilities::FuzzyCompare(input[i], result[i], 1e-06))
    {
      std::cout << "Expected " << input[i] << " but got " << result[i] << " difference is "
                << input[i] - result[i] << std::endl;
      status++;
    }
  }

  if (status)
  {
    std::cout << "..FAILED" << std::endl;
  }
  else
  {
    std::cout << ".PASSED" << std::endl;
  }
  return status;
}
