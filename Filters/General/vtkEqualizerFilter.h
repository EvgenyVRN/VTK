/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkEqualizerFilter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#ifndef vtkEqualizerFilter_h
#define vtkEqualizerFilter_h

#include "vtkFiltersGeneralModule.h"
#include "vtkTableAlgorithm.h"

#include <string>

/**
 * @class vtkEqualizerFilter
 * @brief The vtkEqualizerFilter class
 */

class VTKFILTERSGENERAL_EXPORT vtkEqualizerFilter : public vtkTableAlgorithm
{
public:
  static vtkEqualizerFilter* New();
  vtkTypeMacro(vtkEqualizerFilter, vtkTableAlgorithm);

  void SetSamplingFrequency(int samplingFrequency);
  int GetSamplingFrequency() const;

  void SetAllColumns(bool useAllColumns);
  bool GetAllColumns() const;

  void SetArray(const std::string& name);
  const std::string& GetArray() const;

  void SetPoints(const std::string& points);
  std::string GetPoints() const;

  void SetSpectrumGain(int spectrumGain);
  int GetSpectrumGain() const;

protected:
  vtkEqualizerFilter();
  ~vtkEqualizerFilter() override;

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

private:
  void ProcessColumn(
    vtkDataArray* array, vtkTable* spectrumTable, vtkTable* resultTable, vtkTable* normalizedTable);

private:
  vtkEqualizerFilter(const vtkEqualizerFilter&) = delete;
  void operator=(const vtkEqualizerFilter&) = delete;

  int SamplingFrequency;
  bool AllColumns;
  std::string Array;
  int SpectrumGain;

  class vtkInternal;
  vtkInternal* Internal;
};

#endif
