#include "vtkEqualizerContextItem.h"

#include "vtkBrush.h"
#include "vtkCommand.h"
#include "vtkContext2D.h"
#include "vtkContextKeyEvent.h"
#include "vtkContextMouseEvent.h"
#include "vtkContextScene.h"
#include "vtkContextTransform.h"
#include "vtkLine.h"
#include "vtkMath.h"
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
  static const int radius{ 4 };
  int freq{ -1 };
  double coef{ 0 };
  EqualizerPoint(int f, double c)
  {
    freq = f;
    coef = c;
  }

  EqualizerPoint(const vtkVector2f& vec)
    : freq(vec.GetX())
    , coef(vec.GetY())
  {
  }
  operator vtkVector2f() const { return vtkVector2f(freq, coef); }
  EqualizerPoint& operator=(const vtkVector2f& pos)
  {
    this->freq = pos.GetX();
    this->coef = pos.GetY();

    return *this;
  }

  bool operator<(const EqualizerPoint& point) { return this->freq < point.freq; }
};

bool operator<(const EqualizerPoint& lhs, const EqualizerPoint& rhs)
{
  return lhs.freq < rhs.freq;
}

bool isNear(vtkVector2f pos1, vtkVector2f pos2, double radius)
{
  auto dist = (pos2.GetX() - pos1.GetX()) * (pos2.GetX() - pos1.GetX()) +
    (pos2.GetY() - pos1.GetY()) * (pos2.GetY() - pos1.GetY());
  return dist < radius*radius;
}

bool isNearLine(vtkVector2f p, vtkVector2f le1, vtkVector2f le2, double radius, double closestPoint[3]=nullptr)
{
  double p1[3] = {le1.GetX(), le1.GetY(), 0.0};
  double p2[3] = {le2.GetX(), le2.GetY(), 0.0};
  double xyz[3] = {p.GetX(), p.GetY(), 0.0};
  double t;

  auto onLine = (vtkLine::DistanceToLine(xyz,p1,p2,t,closestPoint) <= radius*radius);
  bool res = { onLine && t < 1.0 && t > 0.0 };
  return res;
}

// TODO: replace to std::lower_bound
// analog std::lower_bound
// use it because
// in gcc available std::lower_bound only with 3 parameters
// in vc available std::lower_bound only with 4 parameters
template<class ForwardIt, class T>
ForwardIt lowerBound3(ForwardIt first, ForwardIt last, const T& value)
{
  ForwardIt it;
  typename std::iterator_traits<ForwardIt>::difference_type count, step;
  count = std::distance(first, last);

  while (count > 0)
  {
    it = first;
    step = count / 2;
    std::advance(it, step);
    if (*it < value)
    {
      first = ++it;
      count -= step + 1;
    }
    else
      count = step;
  }
  return first;
}
}

// using namespace equalizer;

class vtkEqualizerContextItem::vtkInternal
{
public:
  using EqualizerPoints = std::vector<equalizer::EqualizerPoint>;

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
    this->Points.insert(
      (equalizer::lowerBound3(this->Points.begin(), this->Points.end(), point)), point);
  }

  std::string pointsToString()
  {
    std::stringstream ss;
    for (auto point : this->Points)
      ss << point.freq << "," << point.coef << ";";

    return ss.str();
  }

  void setPoints(const std::string& str)
  {
    this->Points.clear();
    // TODO: refactoring, move parsing string to function
    std::vector<std::string> vecPointsStr{ splitStringByDelimiter(str, ';') };

    for (const auto& point : vecPointsStr)
    {
      std::vector<std::string> pointStr{ splitStringByDelimiter(point, ',') };
      if (pointStr.size() > 1)
      {
        float x = std::stof(pointStr.at(0));
        float y = std::stof(pointStr.at(1));
        this->Points.emplace_back(x, y);
      }
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

    // the first or last point
    if ((this->TakenPoint == 0) || (this->TakenPoint == this->Points.size() - 1))
    {
      const equalizer::EqualizerPoint& point = this->Points.at(this->TakenPoint);
      left = point.freq;
      right = point.freq;
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

  void LeftButtonPressEvent(const vtkVector2f& posScreen, vtkContextTransform* transform)
  {
    // 1.Try to find nearest point
    this->TakenPoint = -1;
    for (size_t i = 0; i < this->Points.size(); ++i)
    {
      auto point = transform->MapToScene(this->Points.at(i));
      if (equalizer::isNear(posScreen, point, equalizer::EqualizerPoint::radius))
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
        auto curPoint = transform->MapToScene(*itCur);
        auto prevPoint = transform->MapToScene(*itPrev);
        double closestPoint[3];
        if (equalizer::isNearLine(
              posScreen, prevPoint, curPoint, equalizer::EqualizerPoint::radius, closestPoint))
        {
          vtkVector2f tmp(closestPoint[0], closestPoint[1]);
          // this->addPoint(transform->MapFromScene(tmp));
          this->Points.insert(itCur, transform->MapFromScene(tmp));
          break;
        }
        itPrev = itCur;
      }
    }
  }

  bool RightButtonPressEvent(const vtkVector2f& posScreen, vtkContextTransform* transform)
  {
    if (this->Points.size() < 3)
      return false;

    for (auto it = this->Points.begin(); it != this->Points.end(); ++it)
    {
      auto point = transform->MapToScene(*it);
      if (equalizer::isNear(posScreen, point, equalizer::EqualizerPoint::radius))
      {
        this->Points.erase(it);
        return true;
      }
    }

    return false;
  }

  bool Hit(const vtkVector2f& pos, vtkContextTransform* transform) const
  {
    auto itPrev = this->Points.cbegin();
    auto itCur = itPrev;

    for (++itCur; itCur != this->Points.cend(); ++itCur)
    {
      auto curPoint = transform->MapToScene(*itCur);
      auto prevPoint = transform->MapToScene(*itPrev);
      if (equalizer::isNearLine(pos, prevPoint, curPoint, equalizer::EqualizerPoint::radius))
        return true;
      itPrev = itCur;
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
  , Brush(vtkBrush::New())
  , Internal(new vtkInternal())
{
  this->Pen->SetColor(0, 0, 0);
  this->Pen->SetWidth(3.0);
  this->Brush->SetColor(0, 0, 0);
  // TODO: add real begin and end points (maybe center)
  this->Internal->addPoint(equalizer::EqualizerPoint(0, 100));
  this->Internal->addPoint(equalizer::EqualizerPoint(1000, 100));
  this->Internal->addPoint(equalizer::EqualizerPoint(500, 50));
}

vtkEqualizerContextItem::~vtkEqualizerContextItem()
{
  this->Pen->Delete();
  this->Brush->Delete();
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

  if (!this->Transform)
    return false;

  //  auto width = scene->GetViewWidth();

  painter->ApplyPen(this->Pen);
  painter->ApplyBrush(this->Brush);

  auto itPrev = this->Internal->Points.cbegin();
  auto itCur = itPrev;

  const equalizer::EqualizerPoint& curPoint = this->Transform->MapToScene(*itCur);
  painter->DrawEllipse(curPoint.freq, curPoint.coef, equalizer::EqualizerPoint::radius,
    equalizer::EqualizerPoint::radius);
  for (++itCur; itCur != this->Internal->Points.cend(); ++itCur)
  {
    auto prevPoint = this->Transform->MapToScene(*itPrev);
    auto curPoint = this->Transform->MapToScene(*itCur);
    painter->DrawLine(prevPoint.GetX(), prevPoint.GetY(), curPoint.GetX(), curPoint.GetY());
    painter->DrawEllipse(curPoint.GetX(), curPoint.GetY(), equalizer::EqualizerPoint::radius,
      equalizer::EqualizerPoint::radius);
    itPrev = itCur;
  }

  return true;
}

bool vtkEqualizerContextItem::Hit(const vtkContextMouseEvent& mouse)
{
  if (!this->Transform)
    return false;

  auto hit = this->Internal->Hit(mouse.GetPos(), this->Transform);
  return this->Visible && hit;
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

    if (!this->Transform)
      return false;

    if (this->Internal->TakenPoint != -1)
    {
      equalizer::EqualizerPoint& point = this->Internal->Points.at(this->Internal->TakenPoint);
      auto scope = this->Internal->GetScopes();
      auto posScene = this->Transform->MapFromScene(mouse.GetPos());
      auto posX = vtkMath::ClampValue<int>(posScene.GetX(), scope.first, scope.second);

      auto posY = posScene.GetY();
      if (posY < 0)
        posY = 0;

      point.freq = posX;
      point.coef = posY;

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
  vtkDebugMacro(<< "MouseButtonPressEvent: pos from scene = "
                << this->Transform->MapFromScene(mouse.GetPos()));

  // if (mouse.GetModifiers() == vtkContextMouseEvent::SHIFT_MODIFIER)
  auto pos = mouse.GetPos();
  if (mouse.GetButton() == vtkContextMouseEvent::LEFT_BUTTON)
  {
    this->MouseState = LEFT_BUTTON_PRESSED;
    this->Internal->LeftButtonPressEvent(pos, this->Transform);
  }
  else if (mouse.GetButton() == vtkContextMouseEvent::RIGHT_BUTTON)
  {
    this->MouseState = RIGHT_BUTTON_PRESSED;
    auto removed = this->Internal->RightButtonPressEvent(pos, this->Transform);
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
  vtkDebugMacro(<< "KeyPressEvent: key = " << key.GetKeyCode());
  return vtkAbstractContextItem::KeyPressEvent(key);
}

void vtkEqualizerContextItem::SetScene(vtkContextScene* scene)
{
  vtkAbstractContextItem::SetScene(scene);
  //  if (this->Transform && scene)
  //    this->Internal->setPoints(this->PointsStr, this->Transform);
  this->Modified();
}

void vtkEqualizerContextItem::SetPoints(const std::string& points)
{
  this->Internal->setPoints(points);
  this->Modified();
}

std::string vtkEqualizerContextItem::GetPoints() const
{
  return this->Internal->pointsToString();
}

void vtkEqualizerContextItem::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkContextItem::PrintSelf(os, indent);
}
