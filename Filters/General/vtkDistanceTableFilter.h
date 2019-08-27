#ifndef vtkDistanceTableFilter_h
#define vtkDistanceTableFilter_h

#include "vtkFiltersGeneralModule.h"
#include "vtkTableAlgorithm.h"

#include <string>

class VTKFILTERSGENERAL_EXPORT vtkDistanceTableFilter : public vtkTableAlgorithm
{
public:
  static vtkDistanceTableFilter *New();
  vtkTypeMacro(vtkDistanceTableFilter, vtkTableAlgorithm);
  // void PrintSelf(ostream& os, vtkIndent indent) override;

  void SetFunction(const std::string& function);
  const std::string& GetFunction() const;

  //@{
  /**
   * Define one end of the line (small scalar values).  Default is
   * (0,0,0).
   */
  vtkSetVector3Macro(LowPoint,double);
  vtkGetVectorMacro(LowPoint,double,3);
  //@}

  //@{
  /**
   * Define other end of the line (large scalar values).  Default is
   * (0,0,1).
   */
  vtkSetVector3Macro(HighPoint,double);
  vtkGetVectorMacro(HighPoint,double,3);
  //@}

  //@{
  /**
   * Specify range to map scalars into.  Default is [0, 1].
   */
  vtkSetVector2Macro(ScalarRange,double);
  vtkGetVectorMacro(ScalarRange,double,2);
  //@}

  void SetPosition(int pos);
  int GetPosition() const;

protected:
  vtkDistanceTableFilter();
  ~vtkDistanceTableFilter() override;

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

private:
  vtkDistanceTableFilter(const vtkDistanceTableFilter&) = delete;
  void operator=(const vtkDistanceTableFilter&) = delete;

  std::string Function;
  double LowPoint[3];
  double HighPoint[3];
  double ScalarRange[2];
  int Position;
};

#endif // vtkDistanceTableFilter_h
