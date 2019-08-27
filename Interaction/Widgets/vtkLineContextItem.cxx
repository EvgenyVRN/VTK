#include "vtkLineContextItem.h"

#include "vtkCommand.h"
#include "vtkContext2D.h"
#include "vtkContextKeyEvent.h"
#include "vtkContextMouseEvent.h"
#include "vtkContextScene.h"
#include "vtkContextTransform.h"
#include "vtkObjectFactory.h"
#include "vtkPen.h"

vtkStandardNewMacro(vtkLineContextItem);

vtkLineContextItem::vtkLineContextItem()
  : vtkContextItem()
  , Position(0)
  , MouseState(vtkLineContextItem::NO_BUTTON)
  , Pen(vtkPen::New())
{
  this->Pen->SetColor(0, 0, 0);
  this->Pen->SetWidth(5.0);
}

vtkLineContextItem::~vtkLineContextItem()
{
  this->Pen->Delete();
}

void vtkLineContextItem::PrintSelf(std::ostream& os, vtkIndent indent)
{
  vtkContextItem::PrintSelf(os, indent);
}

void vtkLineContextItem::Update()
{
  vtkAbstractContextItem::Update();
}

bool vtkLineContextItem::Paint(vtkContext2D* painter)
{
  vtkDebugMacro(<< "Paint event called.");
  if (!this->Visible)
    return false;

  vtkContextScene* scene = this->GetScene();
  if (!scene)
    return false;

  auto height = scene->GetViewHeight();

  painter->ApplyPen(this->Pen);
  painter->DrawLine(this->Position, 0, this->Position, height);
  return true;
}

bool vtkLineContextItem::Hit(const vtkContextMouseEvent& mouse)
{
  if (!this->Visible)
    return false;

  vtkVector2f pos = mouse.GetPos();
  auto tolerance = 10;
  if (fabs(pos.GetX() - this->Position) > tolerance)
    return false;

  return true;
}

bool vtkLineContextItem::MouseEnterEvent(const vtkContextMouseEvent& mouse)
{
  vtkDebugMacro(<< "MouseEnterEvent: pos = " << mouse.GetPos());
  return true;
}

bool vtkLineContextItem::MouseMoveEvent(const vtkContextMouseEvent& mouse)
{
  if (this->MouseState == LEFT_BUTTON_PRESSED)
  {
    vtkContextScene* scene = this->GetScene();
    if (!scene)
      return false;

    auto width = scene->GetViewWidth();
    auto new_pos = mouse.GetPos().GetX();
    if (new_pos < 0)
      new_pos = 0;
    if (new_pos > width)
      new_pos = width;

    this->Position = new_pos;

    this->InvokeEvent(vtkCommand::InteractionEvent);
    this->Modified();
    return true;
  }

  return false;
}

bool vtkLineContextItem::MouseLeaveEvent(const vtkContextMouseEvent& mouse)
{
  vtkDebugMacro(<< "MouseLeaveEvent: pos = " << mouse.GetPos());
  return true;
}

bool vtkLineContextItem::MouseButtonPressEvent(const vtkContextMouseEvent& mouse)
{
  this->MouseState = LEFT_BUTTON_PRESSED;
  vtkDebugMacro(<< "MouseButtonPressEvent: pos = " << mouse.GetPos());
  this->InvokeEvent(vtkCommand::StartInteractionEvent);
  return true;
}

bool vtkLineContextItem::MouseButtonReleaseEvent(const vtkContextMouseEvent& mouse)
{
  this->MouseState = NO_BUTTON;
  vtkDebugMacro(<< "MouseButtonReleaseEvent: pos = " << mouse.GetPos());
  this->InvokeEvent(vtkCommand::EndInteractionEvent);
  return true;
}

bool vtkLineContextItem::MouseWheelEvent(const vtkContextMouseEvent& mouse, int delta)
{
  vtkContextScene* scene = this->GetScene();
  if (!scene)
    return false;

  auto width = scene->GetViewWidth();

  auto new_pos = this->Position + delta;
  if (new_pos >= 0 && new_pos <= width)
  {
    this->Position = new_pos;
    this->Modified();
  }

  return true;
}

bool vtkLineContextItem::KeyPressEvent(const vtkContextKeyEvent& key)
{
  vtkDebugMacro(<< "vtkLineContextItem::KeyPressEvent: key = " << key.GetKeyCode());
  return vtkAbstractContextItem::KeyPressEvent(key);
}

void vtkLineContextItem::SetScene(vtkContextScene* scene)
{
  vtkAbstractContextItem::SetScene(scene);
  this->Position = 0;
  this->Modified();
}

void vtkLineContextItem::SetPosition(int pos)
{
  vtkContextScene* scene = this->GetScene();
  if (scene)
  {
    if (this->Transform)
    {
      vtkVector2f mapToScene = this->Transform->MapToScene(vtkVector2f(pos, 0));
      pos = mapToScene.GetX();
    }
  }
  this->Position = pos;
  this->Modified();
}

int vtkLineContextItem::GetPosition() const
{
  auto pos = this->Position;
  if (this->Transform)
  {
    vtkVector2f mapFromScene = this->Transform->MapFromScene(vtkVector2f(pos, 0));
    pos = mapFromScene.GetX();
  }
  return pos;
}
