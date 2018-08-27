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
//! \file PeIntersectionRatioTest.cpp
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"

#include "pe/rigidbody/Ellipsoid.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/SetBodyTypeIDs.h"
#include "pe/Materials.h"

#include "pe_coupling/geometry/PeIntersectionRatio.h"

namespace pe_intersection_ratio_test
{

///////////
// USING //
///////////

using namespace walberla;

typedef boost::tuple<pe::Sphere, pe::Plane, pe::Ellipsoid> BodyTypeTuple;

/*!\brief TODO
 */
//////////
// MAIN //
//////////
int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   pe::SetBodyTypeIDs<BodyTypeTuple>::execute(); //important to be able to compare static body types in intersection function!

   const real_t epsilon( 1e-5_r );

   walberla::id_t sid = 0;
   walberla::id_t uid = 0;

   Vector3<real_t> rPos( 0_r);
   Vector3<real_t> rotationAngles( 0_r);
   Quaternion<real_t> quat( rotationAngles );
   pe::MaterialID material = pe::Material::find("iron");


   ////////////
   // SPHERE //
   ////////////
   {
      Vector3<real_t> bodyPos(1_r, 0_r, 0_r);
      real_t radius = 1_r;

      pe::Sphere sphere(++sid, ++uid, bodyPos, rPos, quat, radius, material, false, false, false);

      pe::RigidBody & rb = sphere; // otherwise not the pe_coupling/geometry version is matched

      Vector3<real_t> pos1(-0.5_r, 0_r, 0_r);
      Vector3<real_t> dir1(1_r, 0_r, 0_r);
      real_t delta1 = walberla::lbm::intersectionRatio(rb, pos1, dir1, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta1, 0.5_r, "Intersection ratio with sphere wrong!");

      Vector3<real_t> pos2(1_r, 1_r, 1_r);
      Vector3<real_t> dir2(0_r, -1_r, -1_r);
      real_t delta2 = walberla::lbm::intersectionRatio(rb, pos2, dir2, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta2, (std::sqrt(2) - 1_r) / std::sqrt(2), "Intersection ratio with sphere wrong!");
   }

   ///////////
   // PLANE //
   ///////////
   {
      Vector3<real_t> bodyPos(1_r, 0_r, 0_r);
      Vector3<real_t> bodyNormal(0_r, 1_r, 1_r);

      bodyNormal = bodyNormal.getNormalized();

      pe::Plane plane(++sid, ++uid, bodyPos, bodyNormal, bodyPos * bodyNormal, material);

      pe::RigidBody & rb = plane; // otherwise not the pe_coupling/geometry version is matched

      Vector3<real_t> pos1(1_r, 0.5_r, 0.5_r);
      Vector3<real_t> dir1(0_r, -1_r, -1_r);
      real_t delta1 = walberla::lbm::intersectionRatio(rb, pos1, dir1, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta1, 0.5_r, "Intersection ratio with plane wrong!");

      Vector3<real_t> dir2(0_r, 0_r, -2_r);
      real_t delta2 = walberla::lbm::intersectionRatio(rb, pos1, dir2, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta2, 0.5_r, "Intersection ratio with plane wrong!");

      Vector3<real_t> dir3(0_r, -3_r, 0_r);
      real_t delta3 = walberla::lbm::intersectionRatio(rb, pos1, dir3, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta3, 1_r/3_r, "Intersection ratio with plane wrong!");
   }

   ///////////////
   // ELLIPSOID //
   ///////////////
   {
      Vector3<real_t> bodyPos(1_r, 0_r, 0_r);
      Vector3<real_t> semiAxes1(1_r, 1_r, 1_r);

      pe::Ellipsoid ellip1(++sid, ++uid, bodyPos, rPos, quat, semiAxes1, material, false, false, false);

      pe::RigidBody & rb1 = ellip1; // otherwise not the pe_coupling/geometry version is matched

      Vector3<real_t> pos1(-0.5_r, 0_r, 0_r);
      Vector3<real_t> dir1(1_r, 0_r, 0_r);
      real_t delta1 = walberla::lbm::intersectionRatio(rb1, pos1, dir1, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta1, 0.5_r, "Intersection ratio with ellipsoid wrong!");

      Vector3<real_t> pos2(1_r, 1_r, 1_r);
      Vector3<real_t> dir2(0_r, -1_r, -1_r);
      real_t delta2 = walberla::lbm::intersectionRatio(rb1, pos2, dir2, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta2, (std::sqrt(2) - 1_r) / std::sqrt(2), "Intersection ratio with ellipsoid wrong!");

      Vector3<real_t> semiAxes2(2_r, 0.5_r, 2_r);
      pe::Ellipsoid ellip2(++sid, ++uid, bodyPos, rPos, quat, semiAxes2, material, false, false, false);

      pe::RigidBody & rb2 = ellip2; // otherwise not the pe_coupling/geometry version is matched

      Vector3<real_t> pos3(1_r, 1_r, 0_r);
      Vector3<real_t> dir3(0_r, -1_r, 0_r);
      real_t delta3 = walberla::lbm::intersectionRatio(rb2, pos3, dir3, epsilon );
      WALBERLA_CHECK_FLOAT_EQUAL(delta3, 0.5_r, "Intersection ratio with ellipsoid wrong!");

   }

   return 0;

}

} //namespace pe_intersection_ratio_test

int main( int argc, char **argv ){
   pe_intersection_ratio_test::main(argc, argv);
}
