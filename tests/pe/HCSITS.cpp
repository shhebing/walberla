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
//! \file HCSITS.cpp
//! \brief checks equality of hash grids and simple ccd
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "pe/basic.h"

#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"

#include "core/debug/TestSubsystem.h"

namespace walberla {
using namespace walberla::pe;

typedef boost::tuple<Sphere, Plane> BodyTuple ;

void normalReactionTest(cr::HCSITS& cr, SphereID sp)
{
   contactThreshold = Thresholds<real_t>::contactThreshold();
   // plane at 5,5,5
   // radius 1.1
   sp->setPosition(  Vec3(5,5,6) );
   sp->setLinearVel( Vec3(0,0,0) );
   cr.setErrorReductionParameter( 1.0_r );
   cr.timestep( real_c( 1.0_r ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,6.1_r) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,0.1_r) );

   sp->setPosition(  Vec3(5,5,6) );
   sp->setLinearVel( Vec3(0,0,0) );
   cr.setErrorReductionParameter( 0.5_r );
   cr.timestep( real_c( 1.0_r ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,6.05_r) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,0.05_r) );

   sp->setPosition(  Vec3(5,5,6) );
   sp->setLinearVel( Vec3(0,0,0) );
   cr.setErrorReductionParameter( 0.0_r );
   cr.timestep( real_c( 1.0_r ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,6) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,0) );

   sp->setPosition(  Vec3(5,5,6) );
   sp->setLinearVel( Vec3(0,0,-1) );
   cr.setErrorReductionParameter( 1.0_r );
   cr.timestep( real_c( 1.0_r ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,6.1_r) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,0.1_r) );

   sp->setPosition(  Vec3(5,5,6) );
   sp->setLinearVel( Vec3(0,0,-1) );
   cr.setErrorReductionParameter( 0.5_r );
   cr.timestep( real_c( 1.0_r ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,6.05_r) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,0.05_r) );

   sp->setPosition(  Vec3(5,5,6) );
   sp->setLinearVel( Vec3(0,0,-1) );
   cr.setErrorReductionParameter( 0.0_r );
   cr.timestep( real_c( 1.0_r ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,6) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,0) );

   sp->setPosition(  Vec3(5,5,6.2_r) );
   sp->setLinearVel( Vec3(0,0,-0.2_r) );
   cr.setErrorReductionParameter( 1.0_r );
   cr.timestep( real_c( 1.0_r ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,6.0_r) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,-0.2_r) );

   contactThreshold = 1.0_r;
   sp->setPosition(  Vec3(5,5,6.2_r) );
   sp->setLinearVel( Vec3(0,0,-0.2_r) );
   cr.setErrorReductionParameter( 1.0_r );
   cr.timestep( real_c( 1.0_r ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,6.1_r) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,-0.1_r) );

   cr.timestep( real_c( 1.0_r ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,6.1_r) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,0_r) );
   contactThreshold = Thresholds<real_t>::contactThreshold();
}

void speedLimiterTest(cr::HCSITS& cr, SphereID sp)
{
   cr.setErrorReductionParameter( 1.0_r );

   sp->setPosition(  Vec3(5,5,6) );
   sp->setLinearVel( Vec3(0,0,0) );
   cr.setSpeedLimiter( true, 0.2_r );
   cr.timestep( real_c( 1.0_r ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,6.1_r) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,0.1_r) );

   sp->setPosition(  Vec3(5,5,5.5) );
   sp->setLinearVel( Vec3(0,0,0) );
   cr.setSpeedLimiter( true, 0.2_r );
   cr.timestep( real_c( 1.0_r ) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getPosition() , Vec3(5,5,5.94_r) );
   WALBERLA_CHECK_FLOAT_EQUAL( sp->getLinearVel(), Vec3(0,0,0.44_r) );
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();

   // create blocks
   shared_ptr< StructuredBlockForest > forest = blockforest::createUniformBlockGrid(
            math::AABB(0,0,0,10,10,10),
            uint_c( 1), uint_c( 1), uint_c( 1), // number of blocks in x,y,z direction
            uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
            true,                               // max blocks per process
            true, true, true,                   // full periodicity
            false);

   SetBodyTypeIDs<BodyTuple>::execute();

   auto storageID           = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto hccdID              = forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "HCCD");
   auto fcdID               = forest->addBlockData(fcd::createGenericFCDDataHandling<BodyTuple, fcd::AnalyticCollideFunctor>(), "FCD");
   cr::HCSITS cr(globalBodyStorage, forest->getBlockStoragePointer(), storageID, hccdID, fcdID);
   cr.setMaxIterations( 10 );
   cr.setRelaxationParameter    ( 0.7_r );
   cr.setErrorReductionParameter( 1.0_r );
   cr.setGlobalLinearAcceleration( Vec3(0,0,0) );

   pe::createPlane( *globalBodyStorage, 0, Vec3(0, 0, 1), Vec3(5, 5, 5) );

   SphereID sp = pe::createSphere(
            *globalBodyStorage,
            forest->getBlockStorage(),
            storageID,
            999999999,
            Vec3(5,5,6),
            real_c(1.1));

   logging::Logging::instance()->setStreamLogLevel(logging::Logging::PROGRESS);

   WALBERLA_LOG_PROGRESS("Normal Reaction Test: InelasticFrictionlessContact");
   cr.setRelaxationModel( cr::HardContactSemiImplicitTimesteppingSolvers::InelasticFrictionlessContact );
   normalReactionTest(cr, sp);
   WALBERLA_LOG_PROGRESS( "Normal Reaction Test: ApproximateInelasticCoulombContactByDecoupling");
   cr.setRelaxationModel( cr::HardContactSemiImplicitTimesteppingSolvers::ApproximateInelasticCoulombContactByDecoupling );
   normalReactionTest(cr, sp);
   //    WALBERLA_LOG_PROGRESS( "Normal Reaction Test: ApproximateInelasticCoulombContactByOrthogonalProjections");
   //    cr.setRelaxationModel( cr::HardContactSemiImplicitTimesteppingSolvers::ApproximateInelasticCoulombContactByOrthogonalProjections );
   //    normalReactionTest(cr, sp);
   WALBERLA_LOG_PROGRESS( "Normal Reaction Test: InelasticCoulombContactByDecoupling");
   cr.setRelaxationModel( cr::HardContactSemiImplicitTimesteppingSolvers::InelasticCoulombContactByDecoupling );
   normalReactionTest(cr, sp);
   //    WALBERLA_LOG_PROGRESS( "Normal Reaction Test: InelasticCoulombContactByOrthogonalProjections");
   //    cr.setRelaxationModel( cr::HardContactSemiImplicitTimesteppingSolvers::InelasticCoulombContactByOrthogonalProjections );
   //    normalReactionTest(cr, sp);
   WALBERLA_LOG_PROGRESS( "Normal Reaction Test: InelasticGeneralizedMaximumDissipationContact");
   cr.setRelaxationModel( cr::HardContactSemiImplicitTimesteppingSolvers::InelasticGeneralizedMaximumDissipationContact );
   normalReactionTest(cr, sp);

   WALBERLA_LOG_PROGRESS("SpeedLimiter Test: InelasticFrictionlessContact");
   cr.setRelaxationModel( cr::HardContactSemiImplicitTimesteppingSolvers::InelasticFrictionlessContact );
   speedLimiterTest(cr, sp);

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}