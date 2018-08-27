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
//! \file PeSubCyclingTest.cpp
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"

#include "pe/basic.h"
#include "pe/utility/DestroyBody.h"
#include "pe/cr/DEM.h"

#include <pe_coupling/utility/all.h>

namespace pe_sub_cycling_test
{

///////////
// USING //
///////////

using namespace walberla;

using BodyTypeTuple = boost::tuple<pe::Sphere> ;

/*!\brief test case to check functionality of sub cycling in the pe time step provided by the coupling module
 *
 * During this time step, currently acting forces/torques on a body are kept constant.
 *
 * This test first computes the resulting position offset as well as the linear and rotational velocity after a
 * certain force has been applied to a sphere when using several 'real' pe time steps and re-applying the force 10 times.
 *
 * It then checks the pe time step sub-cycling functionality of the coupling module by creating same-sized spheres
 * at different locations and applying the same force. But this time, 10 sub-cycles are carried out.
 * As a result, the same position offset and velocities must be obtained as in the regular case.
 *
 */
//////////
// MAIN //
//////////
int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   // uncomment to have logging
   //logging::Logging::instance()->setLogLevel(logging::Logging::LogLevel::DETAIL);

   const real_t dx     = 1_r;
   const real_t radius = 5_r;

   ///////////////////////////
   // DATA STRUCTURES SETUP //
   ///////////////////////////

   Vector3<uint_t> blocksPerDirection(uint_t(3), uint_t(1), uint_t(1));
   Vector3<uint_t> cellsPerBlock(uint_t(20), uint_t(20), uint_t(20));
   Vector3<bool> periodicity(true, false, false);

   // create fully periodic domain with refined blocks
   auto blocks = blockforest::createUniformBlockGrid( blocksPerDirection[0], blocksPerDirection[1], blocksPerDirection[2],
                                                      cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2],
                                                      dx,
                                                      0, false, false,
                                                      periodicity[0], periodicity[1], periodicity[2],
                                                      false );


   // pe body storage
   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();
   shared_ptr<pe::BodyStorage> globalBodyStorage = make_shared<pe::BodyStorage>();
   auto bodyStorageID = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "Storage");
   auto sphereMaterialID = pe::createMaterial( "sphereMat", 1_r , 0.3_r, 0.2_r, 0.2_r, 0.24_r, 200_r, 200_r, 0_r, 0_r );

   auto ccdID = blocks->addBlockData(pe::ccd::createHashGridsDataHandling( globalBodyStorage, bodyStorageID ), "CCD");
   auto fcdID = blocks->addBlockData(pe::fcd::createGenericFCDDataHandling<BodyTypeTuple, pe::fcd::AnalyticCollideFunctor>(), "FCD");
   pe::cr::DEM cr(globalBodyStorage, blocks->getBlockStoragePointer(), bodyStorageID, ccdID, fcdID, nullptr);

   // set up synchronization procedure
   const real_t overlap = real_t( 1.5 ) * dx;
   std::function<void(void)> syncCall = std::bind( pe::syncNextNeighbors<BodyTypeTuple>, boost::ref(blocks->getBlockForest()), bodyStorageID, static_cast<WcTimingTree*>(nullptr), overlap, false );


   // sphere positions for test scenarios
   Vector3<real_t> positionInsideBlock(10_r, 10_r, 10_r);
   Vector3<real_t> positionAtBlockBorder(19.9_r, 10_r, 10_r);
   Vector3<real_t> positionAtBlockBorder2(20_r + radius + overlap - 0.1_r, 10_r, 10_r);

   Vector3<real_t> testForce(2_r, 1_r, 0_r);
   Vector3<real_t> torqueOffset = Vector3<real_t>(1_r, 0_r, 0_r);

   uint_t peSubCycles( 10 );
   real_t dtPe( 10_r );
   real_t dtPeSubCycle = dtPe / real_c(peSubCycles);

   pe_coupling::TimeStep timestep(blocks, bodyStorageID, cr, syncCall, dtPe, peSubCycles);

   // evaluate how far the sphere will travel with a specific force applied which is the reference distance for later
   // (depends on the chosen time integrator in the DEM and thus can not generally be computed a priori here)

   Vector3<real_t> expectedPosOffset(0_r);
   Vector3<real_t> expectedLinearVel(0_r);
   Vector3<real_t> expectedAngularVel(0_r);
   {

      const Vector3<real_t>& startPos = positionInsideBlock;

      pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                       startPos, radius, sphereMaterialID, false, true, false);

      for( uint_t t = 0; t < peSubCycles; ++t )
      {
         for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
         {
            for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
            {
               bodyIt->addForceAtPos(testForce, bodyIt->getPosition() + torqueOffset);
            }
         }
         cr.timestep(dtPeSubCycle);
         syncCall();
      }

      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
         {
            expectedPosOffset = bodyIt->getPosition() - startPos;
            expectedLinearVel = bodyIt->getLinearVel();
            expectedAngularVel = bodyIt->getAngularVel();
         }
      }

      mpi::allReduceInplace(expectedPosOffset[0], mpi::SUM);
      mpi::allReduceInplace(expectedPosOffset[1], mpi::SUM);
      mpi::allReduceInplace(expectedPosOffset[2], mpi::SUM);
      mpi::allReduceInplace(expectedLinearVel[0], mpi::SUM);
      mpi::allReduceInplace(expectedLinearVel[1], mpi::SUM);
      mpi::allReduceInplace(expectedLinearVel[2], mpi::SUM);
      mpi::allReduceInplace(expectedAngularVel[0], mpi::SUM);
      mpi::allReduceInplace(expectedAngularVel[1], mpi::SUM);
      mpi::allReduceInplace(expectedAngularVel[2], mpi::SUM);

      WALBERLA_LOG_DEVEL_ON_ROOT(" - expecting position offset: " << expectedPosOffset);
      WALBERLA_LOG_DEVEL_ON_ROOT(" - expecting linear vel: " << expectedLinearVel);
      WALBERLA_LOG_DEVEL_ON_ROOT(" - expecting angular vel: " << expectedAngularVel);

      // clean up
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            bodyIt->markForDeletion();
         }
      }
      syncCall();

   }

   //////////////////
   // Inside block //
   //////////////////
   {
      std::string testIdentifier("Test: sphere inside block");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      const Vector3<real_t>& startPos = positionInsideBlock;

      pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                       startPos, radius, sphereMaterialID, false, true, false);

      syncCall();

      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
         {
            auto pos = bodyIt->getPosition();
            bodyIt->addForceAtPos(testForce, pos+torqueOffset);
         }
      }

      timestep();


      Vector3<real_t> curPosOffset(0_r);
      Vector3<real_t> curLinearVel(0_r);
      Vector3<real_t> curAngularVel(0_r);
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
         {
            curPosOffset = bodyIt->getPosition() - startPos;
            curLinearVel = bodyIt->getLinearVel();
            curAngularVel = bodyIt->getAngularVel();

            WALBERLA_CHECK_FLOAT_EQUAL(curPosOffset[0], expectedPosOffset[0], "Mismatch in posOffset0");
            WALBERLA_CHECK_FLOAT_EQUAL(curPosOffset[1], expectedPosOffset[1], "Mismatch in posOffset1");
            WALBERLA_CHECK_FLOAT_EQUAL(curPosOffset[2], expectedPosOffset[2], "Mismatch in posOffset2");

            WALBERLA_CHECK_FLOAT_EQUAL(curLinearVel[0], expectedLinearVel[0], "Mismatch in linearVel0");
            WALBERLA_CHECK_FLOAT_EQUAL(curLinearVel[1], expectedLinearVel[1], "Mismatch in linearVel1");
            WALBERLA_CHECK_FLOAT_EQUAL(curLinearVel[2], expectedLinearVel[2], "Mismatch in linearVel2");

            WALBERLA_CHECK_FLOAT_EQUAL(curAngularVel[0], expectedAngularVel[0], "Mismatch in angularVel0");
            WALBERLA_CHECK_FLOAT_EQUAL(curAngularVel[1], expectedAngularVel[1], "Mismatch in angularVel1");
            WALBERLA_CHECK_FLOAT_EQUAL(curAngularVel[2], expectedAngularVel[2], "Mismatch in angularVel2");


         }
      }

      // clean up
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            bodyIt->markForDeletion();
         }
      }
      syncCall();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended");

   }

   ///////////////////////
   // At block border 1 //
   ///////////////////////
   {
      std::string testIdentifier("Test: sphere at block border 1");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      const Vector3<real_t>& startPos = positionAtBlockBorder;

      pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                       startPos, radius, sphereMaterialID, false, true, false);

      syncCall();

      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
         {
            auto pos = bodyIt->getPosition();
            bodyIt->addForceAtPos(testForce, pos+torqueOffset);
         }
      }

      timestep();


      Vector3<real_t> curPosOffset(0_r);
      Vector3<real_t> curLinearVel(0_r);
      Vector3<real_t> curAngularVel(0_r);
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
         {
            curPosOffset = bodyIt->getPosition() - startPos;
            curLinearVel = bodyIt->getLinearVel();
            curAngularVel = bodyIt->getAngularVel();

            WALBERLA_CHECK_FLOAT_EQUAL(curPosOffset[0], expectedPosOffset[0], "Mismatch in posOffset0");
            WALBERLA_CHECK_FLOAT_EQUAL(curPosOffset[1], expectedPosOffset[1], "Mismatch in posOffset1");
            WALBERLA_CHECK_FLOAT_EQUAL(curPosOffset[2], expectedPosOffset[2], "Mismatch in posOffset2");

            WALBERLA_CHECK_FLOAT_EQUAL(curLinearVel[0], expectedLinearVel[0], "Mismatch in linearVel0");
            WALBERLA_CHECK_FLOAT_EQUAL(curLinearVel[1], expectedLinearVel[1], "Mismatch in linearVel1");
            WALBERLA_CHECK_FLOAT_EQUAL(curLinearVel[2], expectedLinearVel[2], "Mismatch in linearVel2");

            WALBERLA_CHECK_FLOAT_EQUAL(curAngularVel[0], expectedAngularVel[0], "Mismatch in angularVel0");
            WALBERLA_CHECK_FLOAT_EQUAL(curAngularVel[1], expectedAngularVel[1], "Mismatch in angularVel1");
            WALBERLA_CHECK_FLOAT_EQUAL(curAngularVel[2], expectedAngularVel[2], "Mismatch in angularVel2");

         }
      }

      // clean up
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            bodyIt->markForDeletion();
         }
      }
      syncCall();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended");

   }

   ////////////////////////////
   // At block border 1, mod //
   ////////////////////////////
   {
      std::string testIdentifier("Test: sphere at block border 1 mod");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      const Vector3<real_t>& startPos = positionAtBlockBorder;

      pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                       startPos, radius, sphereMaterialID, false, true, false);

      syncCall();

      // also add on shadow copy, but only half of it to have same total force/torque
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            auto pos = bodyIt->getPosition();
            bodyIt->addForceAtPos(0.5_r*testForce, pos+torqueOffset);
         }
      }

      timestep();


      Vector3<real_t> curPosOffset(0_r);
      Vector3<real_t> curLinearVel(0_r);
      Vector3<real_t> curAngularVel(0_r);
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
         {
            curPosOffset = bodyIt->getPosition() - startPos;
            curLinearVel = bodyIt->getLinearVel();
            curAngularVel = bodyIt->getAngularVel();

            WALBERLA_CHECK_FLOAT_EQUAL(curPosOffset[0], expectedPosOffset[0], "Mismatch in posOffset0");
            WALBERLA_CHECK_FLOAT_EQUAL(curPosOffset[1], expectedPosOffset[1], "Mismatch in posOffset1");
            WALBERLA_CHECK_FLOAT_EQUAL(curPosOffset[2], expectedPosOffset[2], "Mismatch in posOffset2");

            WALBERLA_CHECK_FLOAT_EQUAL(curLinearVel[0], expectedLinearVel[0], "Mismatch in linearVel0");
            WALBERLA_CHECK_FLOAT_EQUAL(curLinearVel[1], expectedLinearVel[1], "Mismatch in linearVel1");
            WALBERLA_CHECK_FLOAT_EQUAL(curLinearVel[2], expectedLinearVel[2], "Mismatch in linearVel2");

            WALBERLA_CHECK_FLOAT_EQUAL(curAngularVel[0], expectedAngularVel[0], "Mismatch in angularVel0");
            WALBERLA_CHECK_FLOAT_EQUAL(curAngularVel[1], expectedAngularVel[1], "Mismatch in angularVel1");
            WALBERLA_CHECK_FLOAT_EQUAL(curAngularVel[2], expectedAngularVel[2], "Mismatch in angularVel2");

         }
      }

      // clean up
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            bodyIt->markForDeletion();
         }
      }
      syncCall();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended");

   }

   ///////////////////////
   // At block border 2 //
   ///////////////////////
   {
      std::string testIdentifier("Test: sphere at block border 2");
      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - started");

      const Vector3<real_t>& startPos = positionAtBlockBorder2;

      pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0,
                       startPos, radius, sphereMaterialID, false, true, false);

      syncCall();

      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
         {
            auto pos = bodyIt->getPosition();
            bodyIt->addForceAtPos(testForce, pos+torqueOffset);
         }
      }

      timestep();


      Vector3<real_t> curPosOffset(0_r);
      Vector3<real_t> curLinearVel(0_r);
      Vector3<real_t> curAngularVel(0_r);
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
         {
            curPosOffset = bodyIt->getPosition() - startPos;
            curLinearVel = bodyIt->getLinearVel();
            curAngularVel = bodyIt->getAngularVel();

            WALBERLA_CHECK_FLOAT_EQUAL(curPosOffset[0], expectedPosOffset[0], "Mismatch in posOffset0");
            WALBERLA_CHECK_FLOAT_EQUAL(curPosOffset[1], expectedPosOffset[1], "Mismatch in posOffset1");
            WALBERLA_CHECK_FLOAT_EQUAL(curPosOffset[2], expectedPosOffset[2], "Mismatch in posOffset2");

            WALBERLA_CHECK_FLOAT_EQUAL(curLinearVel[0], expectedLinearVel[0], "Mismatch in linearVel0");
            WALBERLA_CHECK_FLOAT_EQUAL(curLinearVel[1], expectedLinearVel[1], "Mismatch in linearVel1");
            WALBERLA_CHECK_FLOAT_EQUAL(curLinearVel[2], expectedLinearVel[2], "Mismatch in linearVel2");

            WALBERLA_CHECK_FLOAT_EQUAL(curAngularVel[0], expectedAngularVel[0], "Mismatch in angularVel0");
            WALBERLA_CHECK_FLOAT_EQUAL(curAngularVel[1], expectedAngularVel[1], "Mismatch in angularVel1");
            WALBERLA_CHECK_FLOAT_EQUAL(curAngularVel[2], expectedAngularVel[2], "Mismatch in angularVel2");

         }
      }

      // clean up
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            bodyIt->markForDeletion();
         }
      }
      syncCall();

      WALBERLA_LOG_DEVEL_ON_ROOT(testIdentifier << " - ended");

   }

   return 0;

}

} //namespace pe_sub_cycling_test

int main( int argc, char **argv ){
   pe_sub_cycling_test::main(argc, argv);
}
