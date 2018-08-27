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
//! \file Shading.h
//! \author Lukas Werner
//
//======================================================================================================================

#pragma once

#include <pe/basic.h>
#include <pe/Types.h>
#include <pe/raytracing/Color.h>
#include <pe/raytracing/ShadingParameters.h>
#include <core/mpi/MPIWrapper.h>
#include <core/mpi/MPIManager.h>

namespace walberla {
namespace pe {
namespace raytracing {

inline ShadingParameters defaultBodyTypeDependentShadingParams (const BodyID body);
inline ShadingParameters processRankDependentShadingParams (const BodyID body);
inline ShadingParameters defaultShadingParams (const BodyID body);
inline ShadingParameters blackShadingParams (const BodyID body);
inline ShadingParameters whiteShadingParams (const BodyID body);
inline ShadingParameters lightGreyShadingParams (const BodyID body);
inline ShadingParameters greyShadingParams (const BodyID body);
inline ShadingParameters darkGreyShadingParams (const BodyID body);
inline ShadingParameters redShadingParams (const BodyID body);
inline ShadingParameters greenShadingParams (const BodyID body);
inline ShadingParameters blueShadingParams (const BodyID body);
inline ShadingParameters violetShadingParams (const BodyID body);
inline ShadingParameters yellowShadingParams (const BodyID body);

inline ShadingParameters defaultBodyTypeDependentShadingParams (const BodyID body) {
   auto bodyTypeID = body->getTypeID();
   
   if (bodyTypeID == Plane::getStaticTypeID()) {
      return lightGreyShadingParams(body).makeMatte();
   } else if (bodyTypeID == Sphere::getStaticTypeID()) {
      return redShadingParams(body).makeGlossy();
   } else if (bodyTypeID == Capsule::getStaticTypeID()) {
      return blueShadingParams(body).makeGlossy();
   } else if (bodyTypeID == Box::getStaticTypeID()) {
      return violetShadingParams(body);
   } else if (bodyTypeID == Ellipsoid::getStaticTypeID()) {
      return yellowShadingParams(body).makeGlossy(60);
   } else {
      return defaultShadingParams(body);
   }
}

inline ShadingParameters processRankDependentShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   int numProcesses = mpi::MPIManager::instance()->numProcesses();
   int rank = mpi::MPIManager::instance()->rank();
   
   real_t hue = 360_r * real_t(rank)/real_t(numProcesses);
   Color color = Color::colorFromHSV(hue, 1_r, 0.9_r);
   
   return ShadingParameters(color,
                            color*0.5_r,
                            Color(0,0,0),
                            0_r);
}
   
inline ShadingParameters defaultShadingParams (const BodyID body) {
   return greyShadingParams(body);
}
   
inline ShadingParameters whiteShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(1_r, 1_r, 1_r),
                       Color(0.9_r, 0.9_r, 0.9_r),
                       Color(0_r, 0_r, 0_r),
                       0_r);
   return s;
}
   
inline ShadingParameters blackShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(0_r, 0_r, 0_r),
                       Color(0_r, 0_r, 0_r),
                       Color(0.1_r, 0.1_r, 0.1_r),
                       0_r);
   return s;
}

inline ShadingParameters lightGreyShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(0.82_r, 0.82_r, 0.82_r),
                       Color(0.5_r, 0.5_r, 0.5_r),
                       Color(0_r, 0_r, 0_r),
                       0_r);
   return s;
}

inline ShadingParameters greyShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(0.5_r, 0.5_r, 0.5_r),
                       Color(0.4_r, 0.4_r, 0.4_r),
                       Color(0.1_r, 0.1_r, 0.1_r),
                       0_r);
   return s;
}

inline ShadingParameters darkGreyShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(0.2_r, 0.2_r, 0.2_r),
                       Color(0.06_r, 0.06_r, 0.06_r),
                       Color(0.1_r, 0.1_r, 0.1_r),
                       0_r);
   return s;
}
   
inline ShadingParameters redShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(1_r, 0_r, 0_r),
                       Color(0.5_r, 0_r, 0_r),
                       Color(0.1_r, 0.1_r, 0.1_r),
                       0_r);
   return s;
}

inline ShadingParameters greenShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(0_r, 0.72_r, 0_r),
                       Color(0_r, 0.41_r, 0_r),
                       Color(0.1_r, 0.1_r, 0.1_r),
                       0_r);
   return s;
}

inline ShadingParameters blueShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(0.15_r, 0.44_r, 0.91_r),
                       Color(0_r, 0_r, 0.4_r),
                       Color(0.1_r, 0.1_r, 0.1_r),
                       0_r);
   return s;
}
   
inline ShadingParameters yellowShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(1_r, 0.96_r, 0_r),
                       Color(0.5_r, 0.48_r, 0_r),
                       Color(0_r, 0_r, 0_r),
                       0_r);
   return s;
}

inline ShadingParameters violetShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(0.6_r, 0_r, 0.9_r),
                       Color(0.5_r, 0_r, 0.8_r),
                       Color(0_r, 0_r, 0_r),
                       0_r);
   return s;
}

} //namespace raytracing
} //namespace pe
} //namespace walberla
