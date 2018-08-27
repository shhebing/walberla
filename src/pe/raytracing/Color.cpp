//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file   Color.cpp
//! \author Lukas Werner
//
//======================================================================================================================

#include "Color.h"

namespace walberla {
namespace pe {
namespace raytracing {

/*!\brief Create a Color object from HSV values.
 * \param hue Hue value in degrees from 0-360
 * \param saturation Saturation value from 0-1
 * \param value Value from 0-1
 */
Color Color::colorFromHSV(real_t hue, real_t saturation, real_t value) {
   // based on Max K. Agoston: Computer Graphics and Geometric Modeling - Implementation and Algorithms
   real_t r, g, b;
   
   if (realIsEqual(hue, 360_r)) {
      hue = 0_r;
   } else {
      hue /= 60_r;
   }
   real_t fract = hue - std::floor(hue);
   
   real_t P = value*(1_r - saturation);
   real_t Q = value*(1_r - saturation*fract);
   real_t T = value*(1_r - saturation*(1_r - fract));
   
   if (0_r <= hue && hue < 1_r) {
      r = value;
      g = T;
      b = P;
   } else if (1_r <= hue && hue < 2_r) {
      r = Q;
      g = value,
      b = P;
   } else if (2_r <= hue && hue < 3_r) {
      r = P;
      g = value;
      b = T;
   } else if (3_r <= hue && hue < 4_r) {
      r = P;
      g = Q;
      b = value;
   } else if (4_r <= hue && hue < 5_r) {
      r = T;
      g = P;
      b = value;
   } else if (5_r <= hue && hue < 6_r) {
      r = value;
      g = P;
      b = Q;
   } else {
      r = g = b = 0_r;
   }
   
   return Color(r, g, b);
}

}
}
}
