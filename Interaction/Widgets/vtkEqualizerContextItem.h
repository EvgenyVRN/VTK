#ifndef vtkEqualizerContextItem_h
#define vtkEqualizerContextItem_h

#include "vtkContextItem.h"
#include "vtkInteractionWidgetsModule.h"
#include "vtkObject.h"

#include <string>

class vtkPen;
class vtkContextTransform;

/**
 * @class vtkEqualizerContextItem
 *
 * @brief The vtkEqualizerContextItem class
 */

class VTKINTERACTIONWIDGETS_EXPORT vtkEqualizerContextItem : public vtkContextItem
{
public:
  static vtkEqualizerContextItem* New();

  vtkTypeMacro(vtkEqualizerContextItem, vtkContextItem);
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

  void SetPoints(const std::string& points);
  std::string GetPoints() const;

protected:
  enum MouseStates
  {
    NO_BUTTON = 0,
    LEFT_BUTTON_PRESSED = 1,
    RIGHT_BUTTON_PRESSED = 2
  };

  vtkEqualizerContextItem();
  ~vtkEqualizerContextItem() override;

  MouseStates MouseState;
  vtkPen* Pen;

  std::string PointsStr;

private:
  vtkEqualizerContextItem(const vtkEqualizerContextItem&) = delete;
  void operator=(const vtkEqualizerContextItem&) = delete;

  class vtkInternal;
  vtkInternal* Internal;
};

#endif
