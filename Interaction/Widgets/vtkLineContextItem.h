#ifndef vtkLineContextItem_h
#define vtkLineContextItem_h

#include "vtkContextItem.h"
#include "vtkInteractionWidgetsModule.h"
#include "vtkObject.h"

class vtkPen;
class vtkContextTransform;

/**
 * @class vtkLineContextItem
 *
 * @brief The vtkLineContextItem class
 */

class VTKINTERACTIONWIDGETS_EXPORT vtkLineContextItem : public vtkContextItem
{
public:
  static vtkLineContextItem* New();

  vtkTypeMacro(vtkLineContextItem, vtkContextItem);
  void PrintSelf(ostream& os, vtkIndent indent) override;


  void InstantiateHandleRepresentation() {}

  /**
   * Perform any updates to the item that may be necessary before rendering.
   * The scene should take care of calling this on all items before their
   * Paint function is invoked.
   */
  void Update() override;

  /**
   * Paint event for the item, called whenever the item needs to be drawn.
   */
  bool Paint(vtkContext2D* painter) override;

  /**
   * Return true if the supplied x, y coordinate is inside the item.
   */
  bool Hit(const vtkContextMouseEvent& mouse) override;

  /**
   * Mouse enter event.
   * Return true if the item holds the event, false if the event can be
   * propagated to other items.
   */
  bool MouseEnterEvent(const vtkContextMouseEvent& mouse) override;
  bool MouseMoveEvent(const vtkContextMouseEvent& mouse) override;
  bool MouseLeaveEvent(const vtkContextMouseEvent& mouse) override;
  bool MouseButtonPressEvent(const vtkContextMouseEvent& mouse) override;
  bool MouseButtonReleaseEvent(const vtkContextMouseEvent& mouse) override;
  /**
   * Mouse wheel event, positive delta indicates forward movement of the wheel.
   */
  bool MouseWheelEvent(const vtkContextMouseEvent& mouse, int delta) override;
  bool KeyPressEvent(const vtkContextKeyEvent& key) override;

  /**
   * Set the vtkContextScene for the item, always set for an item in a scene.
   */
  void SetScene(vtkContextScene* scene) override;

  void SetPosition(int pos);
  int GetPosition() const;

protected:
  enum MouseStates
  {
    NO_BUTTON = 0,
    LEFT_BUTTON_PRESSED = 1
  };

  vtkLineContextItem();
  ~vtkLineContextItem() override;

  int Position;
  MouseStates MouseState;
  vtkPen* Pen;

private:
  vtkLineContextItem(const vtkLineContextItem&) = delete;
  void operator=(const vtkLineContextItem&) = delete;
};

#endif
