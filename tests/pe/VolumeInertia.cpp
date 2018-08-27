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
//! \file VolumeInertia.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <core/debug/TestSubsystem.h>
#include <core/grid_generator/SCIterator.h>
#include <pe/Materials.h>
#include <pe/rigidbody/Box.h>
#include <pe/rigidbody/Capsule.h>
#include <pe/rigidbody/Ellipsoid.h>
#include <pe/rigidbody/Sphere.h>

using namespace walberla;
using namespace walberla::pe;

template< typename ContainmentT >
void calcNumeric( const ContainmentT & body, const AABB & aabb, const real_t spacing, real_t& outVolume, Vec3& outCOM, Mat3& outInertia )
{
   Vector3<real_t> pointOfReference = aabb.min() + Vector3<real_t>( 0.5_r * spacing );

   uint_t volume = 0;
   math::KahanAccumulator<real_t> centroid[3];
   math::KahanAccumulator<real_t> inertiaTensor[6];
   uint_t numPoints = 0;

   for(grid_generator::SCIterator it( aabb, pointOfReference, spacing ); it != grid_generator::SCIterator(); ++it)
   {
      if(body.containsPoint( *it ))
      {
         //volume
         ++volume;

         //center of mass
         centroid[0] += (*it)[0];
         centroid[1] += (*it)[1];
         centroid[2] += (*it)[2];

         //inertia tensor
         const Vector3<real_t> p = *it;
         const real_t & x = p[0];
         const real_t & y = p[1];
         const real_t & z = p[2];

         inertiaTensor[0] += y*y + z*z;
         inertiaTensor[1] += -x*y;
         inertiaTensor[2] += -x*z;
         inertiaTensor[3] += x*x + z*z;
         inertiaTensor[4] += -y*z;
         inertiaTensor[5] += x*x + y*y;
         ++numPoints;
      }
   }

   auto dV = (spacing * spacing * spacing);
   auto dm = dV * Material::getDensity( body.getMaterial() );
   outVolume  = real_c(volume) * dV;
   outCOM     = Vec3( centroid[0].get(), centroid[1].get(), centroid[2].get() ) / real_c(numPoints);
   outInertia = Mat3( inertiaTensor[0].get(), inertiaTensor[1].get(), inertiaTensor[2].get(),
         inertiaTensor[1].get(), inertiaTensor[3].get(), inertiaTensor[4].get(),
         inertiaTensor[2].get(), inertiaTensor[4].get(), inertiaTensor[5].get() ) * dm;
}

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   MaterialID material = createMaterial( "dummy",
                                         real_t( 3900 ),
                                         0,
                                         0,
                                         0,
                                         real_t( 0.23 ),
                                         real_t( 3.6e11 ),
                                         real_t( 8.65e6 ),
                                         real_t( 1.1e1 ),
                                         0 );

   real_t volume;
   Vec3   COM;
   Mat3   inertia;

   Sphere sp(0, 0, Vec3(0,0,0), Vec3(0,0,0), Quat(), 2.34_r, material, false, true, false);
   calcNumeric(sp, sp.getAABB(), 0.01_r, volume, COM, inertia);
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(sp.getVolume(), volume, 10e-4_r) );
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(Sphere::calcVolume( 2.34_r ), volume, 10e-4_r) );
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(sp.getInertia(), inertia, 10_r) );
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(sp.getMass(), volume * Material::getDensity(material), 10e-4_r) );
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(COM, Vec3(0), 10e-4_r) );

   Box bx(0, 0, Vec3(0,0,0), Vec3(0,0,0), Quat(), Vec3(1.5_r, 2.5_r, 3.5_r), material, false, true, false);
   calcNumeric(bx, bx.getAABB(), 0.01_r, volume, COM, inertia);
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(bx.getVolume(), volume, 10e-4_r) );
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(Box::calcVolume( Vec3(1.5_r, 2.5_r, 3.5_r) ), volume, 10e-4_r) );
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(bx.getInertia(), inertia, 10_r) );
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(bx.getMass(), volume * Material::getDensity(material), 10e-4_r) );
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(COM, Vec3(0), 10e-4_r) );

   Ellipsoid el(0, 0, Vec3(0,0,0), Vec3(0,0,0), Quat(), Vec3(1.5_r, 2.5_r, 3.5_r), material, false, true, false);
   calcNumeric(el, el.getAABB(), 0.01_r, volume, COM, inertia);
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(el.getVolume(), volume, 10e-4_r) );
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(Ellipsoid::calcVolume( Vec3(1.5_r, 2.5_r, 3.5_r) ), volume, 10e-4_r) );
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(el.getInertia(), inertia, 10_r) );
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(el.getMass(), volume * Material::getDensity(material), 10e-4_r) );
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(COM, Vec3(0), 10e-4_r) );

   Capsule cp(0, 0, Vec3(0,0,0), Vec3(0,0,0), Quat(), 1.5_r, 2.5_r, material, false, true, false);
   calcNumeric(cp, cp.getAABB(), 0.01_r, volume, COM, inertia);
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(cp.getVolume(), volume, 10e-4_r) );
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(Capsule::calcVolume( 1.5_r, 2.5_r ), volume, 10e-4_r) );
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(cp.getInertia(), inertia, 10_r) );
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(cp.getMass(), volume * Material::getDensity(material), 10e-4_r) );
   WALBERLA_CHECK( walberla::debug::check_functions_detail::check_float_equal_eps(COM, Vec3(0), 10e-4_r) );

   return EXIT_SUCCESS;
}
