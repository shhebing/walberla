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
//! \file LoadFromConfig.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"

#include "pe/basic.h"
#include "pe/utility/CreateWorld.h"

#include "core/debug/TestSubsystem.h"

#include <boost/tuple/tuple.hpp>

namespace walberla {
using namespace walberla::pe;

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();

   Environment env(argc, argv);
   //! [Load Config]
   auto cfg = env.config();
   if (cfg == nullptr) WALBERLA_ABORT("No config specified!");
   const Config::BlockHandle configBlock  = cfg->getBlock( "LoadFromConfig" );
   //! [Load Config]

   //! [Config Get Parameter]
   real_t radius = configBlock.getParameter<real_t>("radius", real_c(0.4) );
   //! [Config Get Parameter]
   WALBERLA_UNUSED(radius);

   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();

   // create blocks
   //! [Config BlockForest]
   shared_ptr<BlockForest> forest = createBlockForestFromConfig( configBlock );
   //! [Config BlockForest]
   WALBERLA_CHECK_EQUAL( forest->getXSize(), 3 );
   WALBERLA_CHECK_EQUAL( forest->getYSize(), 4 );
   WALBERLA_CHECK_EQUAL( forest->getZSize(), 5 );
   WALBERLA_CHECK( !forest->isXPeriodic() );
   WALBERLA_CHECK(  forest->isYPeriodic() );
   WALBERLA_CHECK( !forest->isZPeriodic() );
   WALBERLA_CHECK_FLOAT_EQUAL( forest->getDomain().minCorner(), Vec3(-15, -15, 0) );
   WALBERLA_CHECK_FLOAT_EQUAL( forest->getDomain().maxCorner(), Vec3(-3, 8, 34) );

   //! [Config HCSITS]
   BlockDataID blockDataID;
   cr::HCSITS hcsits( globalBodyStorage, forest, blockDataID, blockDataID, blockDataID);
   configure(configBlock, hcsits);
   //! [Config HCSITS]
   WALBERLA_CHECK_EQUAL( hcsits.getRelaxationModel(), cr::HCSITS::RelaxationModel::ApproximateInelasticCoulombContactByDecoupling );
   WALBERLA_CHECK_EQUAL( hcsits.getMaxIterations(), 123 );
   WALBERLA_CHECK_FLOAT_EQUAL( hcsits.getRelaxationParameter(), 0.123_r );
   WALBERLA_CHECK_FLOAT_EQUAL( hcsits.getErrorReductionParameter(), 0.123_r );
   WALBERLA_CHECK_FLOAT_EQUAL( hcsits.getGlobalLinearAcceleration(), Vec3(1,-2,3) );

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}