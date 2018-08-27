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
//! \file DataTypesTest.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/Environment.h"
#include "core/math/AABB.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"

#include "domain_decomposition/PeriodicIntersectionVolume.h"

#include "stencil/D3Q27.h"

int main( int argc, char** argv )
{
   using namespace walberla;
   using namespace walberla::domain_decomposition;

   debug::enterTestMode();

   Environment walberlaEnv( argc, argv );
   WALBERLA_UNUSED(walberlaEnv);

   std::array< bool, 3 > periodic{{true, true, true}};
   math::AABB domain(0_r,0_r,0_r,100_r,100_r,100_r);
   math::AABB box1(0_r,0_r,0_r,10_r,10_r,10_r);
   math::AABB box2(0_r,0_r,0_r,10_r,10_r,10_r);

   for (int multiple = 0; multiple < 3; ++multiple)
   {
      for (auto dir = stencil::D3Q27::beginNoCenter(); dir != stencil::D3Q27::end(); ++dir)
      {
         Vector3<real_t> shift(
               real_c(dir.cx()) * (9_r + domain.xSize() * real_c(multiple)),
               real_c(dir.cy()) * (9_r + domain.ySize() * real_c(multiple)),
               real_c(dir.cz()) * (9_r + domain.zSize() * real_c(multiple)));
         box2.setCenter(shift + box1.center());
         real_t vol = periodicIntersectionVolume(periodic, domain, box1, box2);

         switch (dir.direction())
         {
         case stencil::BNE:
         case stencil::BNW:
         case stencil::BSE:
         case stencil::BSW:
         case stencil::TNE:
         case stencil::TNW:
         case stencil::TSE:
         case stencil::TSW:
            WALBERLA_CHECK_FLOAT_EQUAL(vol, 1_r);
            break;
         case stencil::BN:
         case stencil::BW:
         case stencil::BS:
         case stencil::BE:
         case stencil::TN:
         case stencil::TW:
         case stencil::TS:
         case stencil::TE:
         case stencil::NE:
         case stencil::NW:
         case stencil::SE:
         case stencil::SW:
            WALBERLA_CHECK_FLOAT_EQUAL(vol, 10_r);
            break;
         case stencil::B:
         case stencil::T:
         case stencil::N:
         case stencil::W:
         case stencil::S:
         case stencil::E:
            WALBERLA_CHECK_FLOAT_EQUAL(vol, 100_r);
            break;
         default:
            WALBERLA_CHECK(false, "Should not end up here!");
            break;
         }
      }
   }

   return EXIT_SUCCESS;
}
