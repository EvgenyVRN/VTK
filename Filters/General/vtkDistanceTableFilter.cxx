#include "vtkDistanceTableFilter.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkTable.h"

vtkStandardNewMacro(vtkDistanceTableFilter);

vtkDistanceTableFilter::vtkDistanceTableFilter()
  : vtkTableAlgorithm()
  , Position(0)
{
}

vtkDistanceTableFilter::~vtkDistanceTableFilter() {}

void vtkDistanceTableFilter::SetFunction(const std::string& function)
{
  this->Function = function;
}

void vtkDistanceTableFilter::SetPosition(int pos)
{
    std::cout << "vtkDistanceTableFilter::SetPosition: " << pos << std::endl;
    this->Position = pos;
    this->Modified();
}

int vtkDistanceTableFilter::GetPosition() const
{
    return this->Position;
}

const std::string& vtkDistanceTableFilter::GetFunction() const
{
  return this->Function;
}

int vtkDistanceTableFilter::RequestData(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkTable* input = vtkTable::GetData(inputVector[0]);

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkTable* output = vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  output->ShallowCopy(input);

  return 1;
}
