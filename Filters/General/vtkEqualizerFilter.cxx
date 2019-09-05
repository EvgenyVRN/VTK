#include "vtkEqualizerFilter.h"

#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkFFT.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkMath.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
#include "vtkTable.h"
#include "vtkVector.h"

#include <vtksys/SystemTools.hxx>
using namespace vtksys;

#include <string>
#include <vector>

class vtkEqualizerFilter::vtkInternal
{
public:
  void ClearPoints() { this->Points.clear(); }
  void SetTable(vtkTable* input)
  {
    if (input != this->TableSrc)
    {
      this->TableSrc = input;
      this->OriginalSize = 0;
      this->SpectrumSize = 0;
      this->Spectrums.clear();
    }
  }

  vtkIdType GetHalfSize() const { return (OriginalSize + 1) / 2; }

  vtkIdType GetHalfSpectrumSize() const { return (SpectrumSize + 1) / 2; }

  const std::vector<ComplexNumber>& GetSpectrum(vtkDataArray* array)
  {
    auto it = this->Spectrums.find(array->GetName());
    if (it == this->Spectrums.end())
    {
      const vtkIdType tuplesCount = array->GetNumberOfTuples();
      double values[tuplesCount];
      for (vtkIdType tupleId = 0; tupleId < tuplesCount; ++tupleId)
        values[tupleId] = array->GetTuple1(tupleId);

      int spectrumSize = 0;

      ComplexNumber* spectrum = nullptr;
      vtkFFT::fft_direct(&values[0], tuplesCount, &spectrumSize, spectrum);
      this->SpectrumSize = spectrumSize;

      this->Spectrums[array->GetName()] =
        std::vector<ComplexNumber>(spectrum, spectrum + spectrumSize);

      delete[] spectrum;
    }

    return this->Spectrums.at(array->GetName());
  }

  const std::vector<double>& GetNormalizedSpectrum(vtkDataArray* array)
  {
    auto it = this->NormalizedSpectrums.find(array->GetName());
    if( it == this->NormalizedSpectrums.end() )
    {
      const std::vector<ComplexNumber>& spectrum = this->GetSpectrum(array);
      auto maxModule = std::numeric_limits<double>::min();
      std::vector<double> modules;
      for (vtkIdType spectrumId = 0; spectrumId < this->GetHalfSpectrumSize(); ++spectrumId)
      {
        auto module = vtkFFT::complex_module(spectrum[spectrumId]);
        maxModule = (maxModule > module ? maxModule : module);
        modules.push_back(module);
      }

      std::vector<double> normSpectrum;
      for(auto module: modules)
        normSpectrum.push_back(module/maxModule);

      this->NormalizedSpectrums[array->GetName()] = normSpectrum;
    }

    return this->NormalizedSpectrums.at(array->GetName());
  }

  std::vector<std::pair<int, double>> GetModifiers(int samplingFrequency) const
  {
    std::vector<std::pair<int, double>> result;
    if (this->Points.size() < 2)
      return result;

    for (size_t i = 0; i < this->Points.size() - 1; i++)
    {
      vtkVector2f p1, p2;
      p1 = this->Points.at(i);
      p2 = this->Points.at(i + 1);
      // Cases for interval (freq[i], freq[i+1])

      int pos1 = p1.GetX() * this->SpectrumSize / samplingFrequency;
      int pos2 = p2.GetX() * this->SpectrumSize / samplingFrequency;
      // 1. Out of region
      if (pos2 < 0 || pos1 > this->GetHalfSpectrumSize())
        continue;
      // 2. Left is out of region, right is in region
      if (pos1 < 0 && pos2 <= this->GetHalfSpectrumSize())
      {
        double y =
          vtkInternal::lineYValue(0, vtkVector2f(pos1, p1.GetY()), vtkVector2f(pos2, p2.GetY()));
        pos1 = 0;
        p1.SetY(y);
      }
      // 3. Right is out of region, left is in region
      if (pos1 >= 0 && pos2 > this->GetHalfSpectrumSize())
      {
        double y = vtkInternal::lineYValue(
          this->GetHalfSpectrumSize(), vtkVector2f(pos1, p1.GetY()), vtkVector2f(pos2, p2.GetY()));
        pos2 = this->GetHalfSpectrumSize();
        p2.SetY(y);
      }
      // 4. Interval covers region
      if (pos1 < 0 && pos2 > this->GetHalfSpectrumSize())
      {
        double y1 =
          vtkInternal::lineYValue(0, vtkVector2f(pos1, p1.GetY()), vtkVector2f(pos2, p2.GetY()));
        double y2 = vtkInternal::lineYValue(
          this->GetHalfSpectrumSize(), vtkVector2f(pos1, p1.GetY()), vtkVector2f(pos2, p2.GetY()));
        pos1 = 0;
        p1.SetY(y1);
        pos2 = this->GetHalfSpectrumSize();
        p2.SetY(y2);
      }
      // 5. Otherwise, interval is inside region
      double delta = (p2.GetY() - p1.GetY()) / (pos2 - pos1);

      for (int j = pos1; j < pos2; j++)
      {
        float coeff = (p1.GetY() + delta * (j - pos1));
        // double modifier = pow(10, 0.05 * coeff);
        result.emplace_back(j, coeff);
      }
    }

    return result;
  }

  static double lineYValue(double x, vtkVector2f le1, vtkVector2f le2)
  {
    return ((le2.GetY() - le1.GetY()) * x +
             (-le1.GetX() * (le2.GetY() - le1.GetY()) + le1.GetY() * (le2.GetX() - le1.GetX()))) /
      (le2.GetX() - le1.GetX());
  }

  static std::vector<std::string> splitStringByDelimiter(const std::string& source, char delim)
  {
    std::stringstream ss(source);
    std::string item;
    std::vector<std::string> result;
    while (std::getline(ss, item, delim))
      result.push_back(std::move(item));

    return result;
  }

  // attributes
  std::vector<vtkVector2f> Points;
  vtkIdType OriginalSize = 0;
  vtkIdType SpectrumSize = 0;
  vtkTable* TableSrc = nullptr;
  std::map<std::string, std::vector<ComplexNumber> > Spectrums;
  std::map<std::string, std::vector<double> > NormalizedSpectrums;
};

vtkStandardNewMacro(vtkEqualizerFilter);

vtkEqualizerFilter::vtkEqualizerFilter()
  : vtkTableAlgorithm()
  , SamplingFrequency(1000)
  , AllColumns(false)
  , SpectrumGain(0)
  , Internal(new vtkInternal())
{
  this->SetNumberOfOutputPorts(3);
}

vtkEqualizerFilter::~vtkEqualizerFilter()
{
  delete this->Internal;
}

void vtkEqualizerFilter::SetSamplingFrequency(int samplingFrequency)
{
  this->SamplingFrequency = samplingFrequency;
  this->Modified();
}

int vtkEqualizerFilter::GetSamplingFrequency() const
{
  return this->SamplingFrequency;
}

void vtkEqualizerFilter::SetAllColumns(bool useAllColumns)
{
  this->AllColumns = useAllColumns;
  this->Modified();
}

bool vtkEqualizerFilter::GetAllColumns() const
{
  return this->AllColumns;
}

void vtkEqualizerFilter::SetArray(const std::string& name)
{
  this->Array = name;
  this->Modified();
}

const std::string& vtkEqualizerFilter::GetArray() const
{
  return this->Array;
}

void vtkEqualizerFilter::SetPoints(const std::string& pointsStr)
{
  // TODO: refactoring: replace to common function
  this->Internal->ClearPoints();
  std::vector<std::string> vecPointsStr { vtkInternal::splitStringByDelimiter(pointsStr, ';') };

  for (const auto& point : vecPointsStr)
  {
    std::vector<std::string> pointStr { vtkInternal::splitStringByDelimiter(point, ',') } ;
    if(pointStr.size()>1)
    {
      float x = std::stof(pointStr.at(0));
      float y = std::stof(pointStr.at(1));
      this->Internal->Points.emplace_back(x, y);
    }
  }

  this->Modified();
}

std::string vtkEqualizerFilter::GetPoints() const
{
  std::stringstream ss;
  for (auto point : this->Internal->Points)
    ss << point.GetX() << "," << point.GetY() << ";";

  return ss.str();
}

void vtkEqualizerFilter::SetSpectrumGain(int spectrumGain)
{
  this->SpectrumGain = spectrumGain;
  this->Modified();
}

int vtkEqualizerFilter::GetSpectrumGain() const
{
  return this->SpectrumGain;
}

int vtkEqualizerFilter::RequestData(
  vtkInformation*, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkTable* input = vtkTable::GetData(inputVector[0]);
  this->Internal->SetTable(input);

  vtkInformation* outInfo0 = outputVector->GetInformationObject(0);
  vtkInformation* outInfo1 = outputVector->GetInformationObject(1);
  vtkInformation* outInfo2 = outputVector->GetInformationObject(2);

  if (!outInfo0 || !outInfo1 || !outInfo2)
  {
    vtkWarningMacro(<< "No output info.");
    return 0;
  }

  vtkTable* spectrumTable = vtkTable::GetData(outInfo0);
  vtkTable* resultTable = vtkTable::GetData(outInfo1);
  vtkTable* normalizedSpectrumTable = vtkTable::GetData(outInfo2);

  if (!input || !spectrumTable || !resultTable || !normalizedSpectrumTable)
  {
    vtkWarningMacro(<< "No input or output.");
    return 0;
  }

  this->Internal->OriginalSize = input->GetNumberOfRows(); // orig size
  if (this->AllColumns)
  {
    vtkIdType numColumns = input->GetNumberOfColumns();

    for (vtkIdType col = 0; col < numColumns; col++)
    {
      this->UpdateProgress((double)col / numColumns);

      vtkDataArray* array = vtkArrayDownCast<vtkDataArray>(input->GetColumn(col));
      if (!array)
        continue;
      if (array->GetNumberOfComponents() != 1)
        continue;
      if (array->GetName())
      {
        if (SystemTools::Strucmp(array->GetName(), "time") == 0)
        {
          resultTable->AddColumn(array);
          continue;
        }
        if (strcmp(array->GetName(), "vtkValidPointMask") == 0)
        {
          resultTable->AddColumn(array);
          continue;
        }
      }
      if (array->IsA("vtkIdTypeArray"))
        continue;

      ProcessColumn(array, spectrumTable, resultTable, normalizedSpectrumTable);
    }
  }
  else
  {
    vtkDataArray* array =
      vtkArrayDownCast<vtkDataArray>(input->GetColumnByName(this->Array.c_str()));
    if (!array)
    {
      vtkDebugMacro(<< " !array");
      return 1;
    }

    if (array->GetNumberOfComponents() != 1)
    {
      vtkDebugMacro(<< "Number of components != 1");
      return 1;
    }
    if (array->GetName())
    {
      if (SystemTools::Strucmp(array->GetName(), "time") == 0)
      {
        resultTable->AddColumn(array);
        return 1;
      }
      if (strcmp(array->GetName(), "vtkValidPointMask") == 0)
      {
        resultTable->AddColumn(array);
        return 1;
      }
    }
    if (array->IsA("vtkIdTypeArray"))
    {
      vtkDebugMacro(<< "vtkIdTypeArray");
      return 1;
    }

    ProcessColumn(array, spectrumTable, resultTable, normalizedSpectrumTable);
  }

  return 1;
}

void vtkEqualizerFilter::ProcessColumn(
  vtkDataArray* array, vtkTable* spectrumTable, vtkTable* resultTable, vtkTable* normalizedTable)
{
  // FFT
  std::vector<ComplexNumber> spectrum = this->Internal->GetSpectrum(array);
  if (spectrum.empty())
  {
    vtkErrorMacro(<< "Spectrum is empty: " << array->GetName());
    return;
  }

  std::vector<double> normSpectrum = this->Internal->GetNormalizedSpectrum(array);
  if (normSpectrum.empty())
  {
    vtkErrorMacro(<< "Normalized spectrum is empty: " << array->GetName());
    return;
  }
  // end FFT

  // modify by equalizer
  auto modifiers = this->Internal->GetModifiers(this->SamplingFrequency);
  for(const auto& modifier: modifiers)
  {
    spectrum[modifier.first] *= modifier.second;
    spectrum[this->Internal->SpectrumSize - modifier.first - 1] *= modifier.second;

    normSpectrum[modifier.first] *= modifier.second;
  }


  // fill spectrum table
  std::vector<double> freqArray =
    vtkFFT::rfftfreq(this->Internal->SpectrumSize, 1.0 / this->SamplingFrequency);
  // assert(freqArray.size() == this->Internal->GetHalfSpectrumSize());

  vtkSmartPointer<vtkDoubleArray> freqColumn = vtkSmartPointer<vtkDoubleArray>::New();
  freqColumn->SetNumberOfComponents(1);
  freqColumn->SetNumberOfTuples(this->Internal->GetHalfSpectrumSize());
  freqColumn->SetName("Frequency");
  for (vtkIdType spectrumId = 0; spectrumId < this->Internal->GetHalfSpectrumSize(); ++spectrumId)
    freqColumn->SetTuple1(spectrumId, freqArray.at(spectrumId));

  spectrumTable->AddColumn(freqColumn);
  normalizedTable->AddColumn(freqColumn);

  vtkSmartPointer<vtkDoubleArray> leadArray = vtkSmartPointer<vtkDoubleArray>::New();
  leadArray->SetNumberOfComponents(1);
  leadArray->SetNumberOfTuples(this->Internal->GetHalfSpectrumSize());
  leadArray->SetName(array->GetName());

  vtkSmartPointer<vtkDoubleArray> normalizedArray = vtkSmartPointer<vtkDoubleArray>::New();
  normalizedArray->SetNumberOfComponents(1);
  normalizedArray->SetNumberOfTuples(this->Internal->GetHalfSpectrumSize());
  normalizedArray->SetName(array->GetName());

  for (vtkIdType spectrumId = 0; spectrumId < this->Internal->GetHalfSpectrumSize(); ++spectrumId)
  {
    const ComplexNumber& value = spectrum[spectrumId];
    // нас интересует только спектр амплитуд, поэтому используем complex_module
    // делим на число элементов, чтобы амплитуды были в милливольтах, а не в суммах Фурье
    double modifier = pow(10, 0.05 * this->SpectrumGain);
    double module = vtkFFT::complex_module(value) * modifier / this->Internal->GetHalfSpectrumSize();
    leadArray->SetTuple1(spectrumId, module);

    normalizedArray->SetTuple1(spectrumId, normSpectrum[spectrumId]*modifier);
  }
  spectrumTable->AddColumn(leadArray);
  normalizedTable->AddColumn(normalizedArray);
  // end fill spectrum table

  // fill result table
  int outCount;
  ComplexNumber* num = nullptr;
  vtkFFT::fft_inverse(spectrum.data(), this->Internal->SpectrumSize, &outCount, num);
  double* outputData = new double[this->Internal->OriginalSize];

  if (num)
  {
    vtkFFT::complexes_to_doubles(outputData, num, this->Internal->OriginalSize);
    delete [] num;
  }

  vtkSmartPointer<vtkDoubleArray> rfftArray = vtkSmartPointer<vtkDoubleArray>::New();
  rfftArray->SetNumberOfComponents(1);
  rfftArray->SetNumberOfTuples(this->Internal->OriginalSize);
  rfftArray->SetName(array->GetName());

  for (int spectrumId = 0; spectrumId < this->Internal->OriginalSize; ++spectrumId)
  {
    double val = outputData[spectrumId];
    rfftArray->SetTuple1(spectrumId, val);
  }
  resultTable->AddColumn(rfftArray);
  // end fill result table
  delete [] outputData;
}
