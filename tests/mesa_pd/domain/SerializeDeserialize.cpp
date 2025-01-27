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
//! \file DynamicRefinement.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/domain/BlockForestDataHandling.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>

#include <blockforest/Initialization.h>
#include <core/debug/TestSubsystem.h>
#include <core/grid_generator/SCIterator.h>
#include <core/mpi/Reduce.h>
#include <core/logging/Logging.h>

namespace walberla {

using namespace walberla::mesa_pd;

void createDump()
{
   const real_t spacing = real_c(1);
   const real_t radius  = real_c(0.5);

   WALBERLA_LOG_INFO_ON_ROOT("*** MESA_PD ***");
   auto ps = std::make_shared<data::ParticleStorage>(100);
   auto ss = std::make_shared<data::ShapeStorage>();

   auto  smallSphere = ss->create<data::Sphere>( radius );
   ss->shapes[smallSphere]->updateMassAndInertia(real_t(2707));

   WALBERLA_LOG_INFO_ON_ROOT("*** BLOCKFOREST ***");
   // create forest
   auto forest = blockforest::createBlockForest( math::AABB(real_t(0),
                                                            real_t(0),
                                                            real_t(0),
                                                            real_t(6),
                                                            real_t(6),
                                                            real_t(6)),
                                                 Vector3<uint_t>(2,2,2),
                                                 Vector3<bool>(true, true, true) );
   forest->saveToFile("SerializeDeserialize.sbf");

   auto bfDataHandlingID = forest->addBlockData(domain::createBlockForestDataHandling(ps), "BFDataHandling");

   WALBERLA_CHECK_EQUAL(forest->size(), 1, "please run with 8 processes -> 1 process per block");

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - START ***");
   //add one global particle
   {
      auto p = ps->create();
      p->setPosition( Vec3(real_t(1)) );
      data::particle_flags::set( p->getFlagsRef(), data::particle_flags::GLOBAL);
      p->setInteractionRadius( radius );
      p->setShapeID( smallSphere );
   }

   for (auto& iBlk : *forest)
   {
      for (auto pt : grid_generator::SCGrid(iBlk.getAABB(), Vector3<real_t>(spacing, spacing, spacing) * real_c(0.5), spacing))
      {
         WALBERLA_CHECK(iBlk.getAABB().contains(pt));

         auto p                       = ps->create();
         p->setPosition( pt );
         p->setInteractionRadius( radius );
         p->setShapeID( smallSphere );
         p->setOwner( mpi::MPIManager::instance()->rank() );
      }
   }
   int64_t numParticles = int64_c(ps->size());
   mpi::reduceInplace(numParticles, mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - END ***");

   WALBERLA_CHECK_EQUAL( ps->size(), 28 );

   WALBERLA_LOG_DEVEL_ON_ROOT("dumping simulation");
   forest->saveBlockData("SerializeDeserialize.dump", bfDataHandlingID);

   WALBERLA_CHECK_EQUAL( ps->size(), 28 );
}

void checkDump()
{
   const real_t spacing = real_c(1);
   const real_t radius  = real_c(0.5);

   WALBERLA_LOG_INFO_ON_ROOT("*** MESA_PD ***");
   auto ps = std::make_shared<data::ParticleStorage>(100);
   auto ss = std::make_shared<data::ShapeStorage>();

   auto  smallSphere = ss->create<data::Sphere>( radius );
   ss->shapes[smallSphere]->updateMassAndInertia(real_t(2707));

   WALBERLA_LOG_INFO_ON_ROOT("*** BLOCKFOREST ***");
   auto forest = make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), "SerializeDeserialize.sbf", true, false );

   //add one global particle
   {
      auto p = ps->create();
      p->setPosition( Vec3(real_t(1)) );
      data::particle_flags::set( p->getFlagsRef(), data::particle_flags::GLOBAL);
      p->setInteractionRadius( radius );
      p->setShapeID( smallSphere );
   }

   WALBERLA_CHECK_EQUAL(ps->size(), 1);
   forest->loadBlockData("SerializeDeserialize.dump", domain::createBlockForestDataHandling(ps), "BFDataHandling");
   WALBERLA_CHECK_EQUAL(ps->size(), 28);

   WALBERLA_LOG_INFO_ON_ROOT("*** CHECKING ***");
   auto pIt = ++ps->begin();
   for (auto& iBlk : *forest)
   {
      for (auto pt : grid_generator::SCGrid(iBlk.getAABB(), Vector3<real_t>(spacing, spacing, spacing) * real_c(0.5), spacing))
      {
         WALBERLA_CHECK(iBlk.getAABB().contains(pt));
         WALBERLA_CHECK_UNEQUAL(pIt, ps->end()); //still particles left
         WALBERLA_CHECK_FLOAT_EQUAL(pIt->getPosition(), pt);
         ++pIt;
      }
   }
   WALBERLA_CHECK_EQUAL(pIt, ps->end()); //all particles have been checked
}

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();

   WALBERLA_MPI_SECTION()
   {
      walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   }

   WALBERLA_LOG_DEVEL_ON_ROOT("*** DUMPING ***");
   createDump();
   WALBERLA_MPI_SECTION()
   {
      WALBERLA_MPI_BARRIER();
   }
   WALBERLA_LOG_DEVEL_ON_ROOT("*** CHECKING ***");
   checkDump();

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
   return walberla::main( argc, argv );
}
