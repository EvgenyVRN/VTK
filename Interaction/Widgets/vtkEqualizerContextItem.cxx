#include "vtkEqualizerContextItem.h"

#include "vtkBrush.h"
#include "vtkCommand.h"
#include "vtkContext2D.h"
#include "vtkContextKeyEvent.h"
#include "vtkContextMouseEvent.h"
#include "vtkContextScene.h"
#include "vtkContextTransform.h"
#include "vtkObjectFactory.h"
#include "vtkPen.h"
#include "vtkVector.h"

#include <set>
#include <sstream>
#include <vector>

namespace equalizer
{
struct EqualizerPoint
{
  static constexpr uint radius = 5;
  int freq;
  double coef;
  EqualizerPoint(int f, double c)
  {
    freq = f;
    coef = c;
  }
  EqualizerPoint()
    : freq(-1)
    , coef(0)
  {
  }
  EqualizerPoint(const vtkVector2f& vec)
    : freq(vec.GetX())
    , coef(vec.GetY())
  {
  }
  operator vtkVector2f() const { return vtkVector2f(freq, coef); }
  void operator=(const vtkVector2f& pos)
  {
    this->freq = pos.GetX();
    this->coef = pos.GetY();
  }

  bool operator<(const EqualizerPoint& point) { return this->freq < point.freq; }
};

bool operator<(const EqualizerPoint& lhs, const EqualizerPoint& rhs)
{
  return lhs.freq < rhs.freq;
}

bool isNear(double x, double y, double radius)
{
  return (y <= x + radius) && (y >= x - radius);
}

bool isNear(vtkVector2f x, vtkVector2f y, double radius)
{
  return isNear(x.GetX(), y.GetX(), radius) && isNear(x.GetY(), y.GetY(), radius);
}

double lineYValue(double x, vtkVector2f le1, vtkVector2f le2)
{
  return ((le2.GetY() - le1.GetY()) * x +
           (-le1.GetX() * (le2.GetY() - le1.GetY()) + le1.GetY() * (le2.GetX() - le1.GetX()))) /
    (le2.GetX() - le1.GetX());
}

bool isNearLine(vtkVector2f p, vtkVector2f le1, vtkVector2f le2, double radius)
{
  double val = lineYValue(p.GetX(), le1, le2);
  return abs(int(val - p.GetY())) < radius;
}
}

// using namespace equalizer;

class vtkEqualizerContextItem::vtkInternal
{
public:
  typedef std::vector<equalizer::EqualizerPoint> EqualizerPoints;

  static std::vector<std::string> splitStringByDelimiter(const std::string& source, char delim)
  {
    std::stringstream ss(source);
    std::string item;
    std::vector<std::string> result;
    while (std::getline(ss, item, delim))
      result.push_back(std::move(item));

    return result;
  }

  void addPoint(const equalizer::EqualizerPoint& point)
  {
    this->Points.insert((std::lower_bound(this->Points.begin(), this->Points.end(), point)), point);
  }

  std::string pointsToString(vtkContextTransform* transform)
  {
    std::stringstream ss;
    for (auto point : this->Points)
    {
      auto scenePoint = transform->MapFromScene(point);
      ss << scenePoint.GetX() << "," << scenePoint.GetY() << ";";
    }

    return ss.str();
  }

  void setPoints(const std::string& str, vtkContextTransform* transform)
  {
    this->Points.clear();
    // TODO: refactoring, move parsing string to function
    std::vector<std::string> vecPointsStr{ splitStringByDelimiter(str, ';') };

    std::vector<vtkVector2f> points;
    for (auto point : vecPointsStr)
    {
      std::vector<std::string> pointStr{ splitStringByDelimiter(point, ',') };
      if (pointStr.size() > 1)
      {
        float x = std::stof(pointStr.at(0));
        float y = std::stof(pointStr.at(1));
        points.push_back(vtkVector2f(x, y));
      }
    }

    for (auto point : points)
    {
      auto viewPoint = transform->MapToScene(point);
      this->Points.push_back(viewPoint);
    }
    this->TakenPoint = -1;
  }

  std::pair<int, int> GetScopes() const
  {
    auto left = 0;
    auto right = std::numeric_limits<int>::max();
    if (this->Points.size() < 2)
      return std::pair<int, int>(left, right);

    if (this->TakenPoint == -1)
      return std::pair<int, int>(left, right);

    // the first point
    if (this->TakenPoint == 0)
    {
      const equalizer::EqualizerPoint& pointRight = this->Points.at(this->TakenPoint + 1);
      right = pointRight.freq;
    }
    // the last point
    else if (this->TakenPoint == this->Points.size() - 1)
    {
      const equalizer::EqualizerPoint& pointLeft = this->Points.at(this->TakenPoint - 1);
      left = pointLeft.freq;
    }
    else
    {
      const equalizer::EqualizerPoint& pointLeft = this->Points.at(this->TakenPoint - 1);
      const equalizer::EqualizerPoint& pointRight = this->Points.at(this->TakenPoint + 1);
      left = pointLeft.freq;
      right = pointRight.freq;
    }

    return std::pair<int, int>(left, right);
  }

  void LeftButtonPressEvent(const vtkVector2f& pos)
  {
    // 1.Try to find nearest point
    this->TakenPoint = -1;
    for (size_t i = 0; i < this->Points.size(); ++i)
    {
      const auto& point = this->Points.at(i);
      if (equalizer::isNear(pos, point, equalizer::EqualizerPoint::radius))
      {
        this->TakenPoint = i;
        break;
      }
    }

    // 2.Try to find nearest line
    if (this->TakenPoint == -1)
    {
      auto itPrev = this->Points.cbegin();
      auto itCur = itPrev;

      for (++itCur; itCur != this->Points.cend(); ++itCur)
      {
        const equalizer::EqualizerPoint& curPoint = *itCur;
        const equalizer::EqualizerPoint& prevPoint = *itPrev;
        if (equalizer::isNearLine(pos, prevPoint, curPoint, equalizer::EqualizerPoint::radius))
        {
          this->addPoint(pos);
          break;
        }
        itPrev = itCur;
      }
    }
  }

  bool RightButtonPressEvent(const vtkVector2f& pos)
  {
    for (auto it = this->Points.begin(); it != this->Points.end(); ++it)
    {
      const auto& point = *it;
      if (equalizer::isNear(pos, point, equalizer::EqualizerPoint::radius))
      {
        this->Points.erase(it);
        return true;
      }
    }

    return false;
  }

  // attributes
  EqualizerPoints Points;
  int TakenPoint = -1;
};

vtkStandardNewMacro(vtkEqualizerContextItem);

vtkEqualizerContextItem::vtkEqualizerContextItem()
  : vtkContextItem()
  , MouseState(vtkEqualizerContextItem::NO_BUTTON)
  , Pen(vtkPen::New())
  , Internal(new vtkInternal())
{
  this->Pen->SetColor(0, 0, 0);
  this->Pen->SetWidth(3.0);
  // TODO: add real begin and end points (maybe center)
  this->Internal->addPoint(equalizer::EqualizerPoint(0, 100));
  this->Internal->addPoint(equalizer::EqualizerPoint(1000, 100));
  this->Internal->addPoint(equalizer::EqualizerPoint(500, 50));
}

vtkEqualizerContextItem::~vtkEqualizerContextItem()
{
  this->Pen->Delete();
  delete this->Internal;
}

void vtkEqualizerContextItem::Update()
{
  vtkAbstractContextItem::Update();
}

bool vtkEqualizerContextItem::Paint(vtkContext2D* painter)
{
  vtkDebugMacro(<< "Paint event called.");
  if (!this->Visible)
    return false;

  if (this->Internal->Points.size() < 2)
    return false;

  vtkContextScene* scene = this->GetScene();
  if (!scene)
    return false;

  //  auto width = scene->GetViewWidth();

  painter->ApplyPen(this->Pen);
  painter->GetBrush()->SetColor(0, 0, 0);

  auto itPrev = this->Internal->Points.cbegin();
  auto itCur = itPrev;

  const equalizer::EqualizerPoint& curPoint = *itCur;
  painter->DrawEllipse(curPoint.freq, curPoint.coef, equalizer::EqualizerPoint::radius,
    equalizer::EqualizerPoint::radius);
  for (++itCur; itCur != this->Internal->Points.cend(); ++itCur)
  {
    const equalizer::EqualizerPoint& curPoint = *itCur;
    const equalizer::EqualizerPoint& prevPoint = *itPrev;
    painter->DrawLine(prevPoint.freq, prevPoint.coef, curPoint.freq, curPoint.coef);
    painter->DrawEllipse(curPoint.freq, curPoint.coef, equalizer::EqualizerPoint::radius,
      equalizer::EqualizerPoint::radius);
    itPrev = itCur;
  }

  return true;
}

bool vtkEqualizerContextItem::Hit(const vtkContextMouseEvent& mouse)
{
  return this->Visible;
}

bool vtkEqualizerContextItem::MouseEnterEvent(const vtkContextMouseEvent& mouse)
{
  vtkDebugMacro(<< "MouseEnterEvent: pos = " << mouse.GetPos());
  return true;
}

bool vtkEqualizerContextItem::MouseMoveEvent(const vtkContextMouseEvent& mouse)
{
  if (this->MouseState == LEFT_BUTTON_PRESSED)
  {
    vtkContextScene* scene = this->GetScene();
    if (!scene)
      return false;

    if (this->Internal->TakenPoint != -1)
    {
      equalizer::EqualizerPoint& point = this->Internal->Points.at(this->Internal->TakenPoint);
      auto scope = this->Internal->GetScopes();
      auto pos = mouse.GetPos().GetX();
      if (pos < scope.first)
        pos = scope.first;
      if (pos > scope.second)
        pos = scope.second;

      point = equalizer::EqualizerPoint(pos, mouse.GetPos().GetY());

      this->InvokeEvent(vtkCommand::InteractionEvent);
      this->Modified();
      return true;
    }
  }

  return true;
}

bool vtkEqualizerContextItem::MouseLeaveEvent(const vtkContextMouseEvent& mouse)
{
  vtkDebugMacro(<< "MouseLeaveEvent: pos = " << mouse.GetPos());
  return true;
}

bool vtkEqualizerContextItem::MouseButtonPressEvent(const vtkContextMouseEvent& mouse)
{
  vtkDebugMacro(<< "MouseButtonPressEvent: pos = " << mouse.GetPos());
  // if (mouse.GetModifiers() == vtkContextMouseEvent::SHIFT_MODIFIER)
  if (mouse.GetButton() == vtkContextMouseEvent::LEFT_BUTTON)
  {
    this->MouseState = LEFT_BUTTON_PRESSED;
    this->Internal->LeftButtonPressEvent(mouse.GetPos());
  }
  else if (mouse.GetButton() == vtkContextMouseEvent::RIGHT_BUTTON)
  {
    this->MouseState = RIGHT_BUTTON_PRESSED;
    auto removed = this->Internal->RightButtonPressEvent(mouse.GetPos());
  }

  this->InvokeEvent(vtkCommand::StartInteractionEvent);
  return true;
}

bool vtkEqualizerContextItem::MouseButtonReleaseEvent(const vtkContextMouseEvent& mouse)
{
  this->MouseState = NO_BUTTON;
  vtkDebugMacro(<< "MouseButtonReleaseEvent: pos = " << mouse.GetPos());
  this->InvokeEvent(vtkCommand::EndInteractionEvent);
  this->Modified();
  return true;
}

bool vtkEqualizerContextItem::MouseWheelEvent(const vtkContextMouseEvent& mouse, int delta)
{
  // TODO: add a logic
  return true;
}

bool vtkEqualizerContextItem::KeyPressEvent(const vtkContextKeyEvent& key)
{
  vtkDebugMacro(<< "vtkLineContextItem::KeyPressEvent: key = " << key.GetKeyCode());
  return vtkAbstractContextItem::KeyPressEvent(key);
}

void vtkEqualizerContextItem::SetScene(vtkContextScene* scene)
{
  vtkAbstractContextItem::SetScene(scene);
  if (this->Transform && scene)
    this->Internal->setPoints(this->PointsStr, this->Transform);
  this->Modified();
}

void vtkEqualizerContextItem::SetPoints(const std::string& points)
{
  if (!this->Transform)
  {
    vtkDebugMacro(<< "vtkEqualizerContextItem hasn't vtkContextTransform");
    this->PointsStr = points;
    return;
  }
  this->Internal->setPoints(points, this->Transform);
  this->Modified();
}

std::string vtkEqualizerContextItem::GetPoints() const
{
  if (!this->Transform)
  {
    vtkDebugMacro(<< "vtkEqualizerContextItem hasn't vtkContextTransform");
    return std::string();
  }
  return this->Internal->pointsToString(this->Transform);
}

void vtkEqualizerContextItem::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkContextItem::PrintSelf(os, indent);
}
