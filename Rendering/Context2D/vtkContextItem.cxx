/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkContextItem.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkContextItem.h"

//-----------------------------------------------------------------------------
vtkContextItem::vtkContextItem()
  : Transform(nullptr)
{
  this->Opacity = 1.0;
}

//-----------------------------------------------------------------------------
vtkContextItem::~vtkContextItem() = default;

//-----------------------------------------------------------------------------
void vtkContextItem::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

void vtkContextItem::SetTransform(vtkContextTransform *transform)
{
  this->Transform = transform;
}
