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
//! \file RBGSTest.cpp
//! \ingroup pde
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Abort.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"
#include "core/mpi/MPIManager.h"

#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/VTKWriter.h"

#include "pde/iterations/CGFixedStencilIteration.h"
#include "pde/iterations/CGIteration.h"

#include "stencil/D2Q5.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/VTKOutput.h"

#include <cmath>

namespace walberla {



typedef GhostLayerField< real_t, 1 > PdeField_T;
using Stencil_T = stencil::D2Q5;
using StencilField_T = pde::CGIteration<Stencil_T>::StencilField_T;



void initU( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & uId )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      if( blocks->atDomainYMaxBorder( *block ) )
      {
         PdeField_T * u = block->getData< PdeField_T >( uId );
         CellInterval xyz = u->xyzSizeWithGhostLayer();
         xyz.yMin() = xyz.yMax();
         for( auto cell = xyz.begin(); cell != xyz.end(); ++cell )
         {
            const Vector3< real_t > p = blocks->getBlockLocalCellCenter( *block, *cell );
            u->get( *cell ) = std::sin( 2_r * math::PI * p[0] ) * std::sinh( 2_r * math::PI * p[1] );
         }
      }
   }
}



void initF( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & fId )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      PdeField_T * f = block->getData< PdeField_T >( fId );
      CellInterval xyz = f->xyzSize();
      for( auto cell = xyz.begin(); cell != xyz.end(); ++cell )
      {
         const Vector3< real_t > p = blocks->getBlockLocalCellCenter( *block, *cell );
         f->get( *cell ) = 4_r * math::PI * math::PI * std::sin( 2_r * math::PI * p[0] ) * std::sinh( 2_r * math::PI * p[1] );
      }
   }
}



void copyWeightsToStencilField( const shared_ptr< StructuredBlockStorage > & blocks, const std::vector<real_t> & weights, const BlockDataID & stencilId )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      StencilField_T * stencil = block->getData< StencilField_T >( stencilId );
      
      WALBERLA_FOR_ALL_CELLS_XYZ(stencil,
         for( auto dir = Stencil_T::begin(); dir != Stencil_T::end(); ++dir )
            stencil->get(x,y,z,dir.toIdx()) = weights[ dir.toIdx() ];
      );
   }
}



template <typename Field_T>
void clearField( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & fieldId )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      block->getData< Field_T >( fieldId )->set( typename Field_T::value_type() );
   }
}



int main( int argc, char** argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   const uint_t processes = uint_c( MPIManager::instance()->numProcesses() );
   if( processes != uint_t(1) && processes != uint_t(4) && processes != uint_t(8) )
      WALBERLA_ABORT( "The number of processes must be equal to 1, 4, or 8!" );

   logging::Logging::printHeaderOnStream();
   WALBERLA_ROOT_SECTION() { logging::Logging::instance()->setLogLevel( logging::Logging::PROGRESS ); }

   bool shortrun = false;
   for( int i = 1; i < argc; ++i )
      if( std::strcmp( argv[i], "--shortrun" ) == 0 ) shortrun = true;

   const uint_t xBlocks = ( processes == uint_t(1) ) ? uint_t(1) : ( ( processes == uint_t(4) ) ? uint_t(2) : uint_t(4) );
   const uint_t yBlocks = ( processes == uint_t(1) ) ? uint_t(1) : uint_t(2);
   const uint_t xCells = ( processes == uint_t(1) ) ? uint_t(200) : ( ( processes == uint_t(4) ) ? uint_t(100) : uint_t(50) );
   const uint_t yCells = ( processes == uint_t(1) ) ? uint_t(100) : uint_t(50);
   const real_t xSize = 2_r;
   const real_t ySize = 1_r;
   const real_t dx = xSize / real_c( xBlocks * xCells + uint_t(1) );
   const real_t dy = ySize / real_c( yBlocks * yCells + uint_t(1) );
   auto blocks = blockforest::createUniformBlockGrid( math::AABB( 0.5_r * dx, 0.5_r * dy, 0_r,
                                                                  xSize - 0.5_r * dx, ySize - 0.5_r * dy, dx ),
                                                      xBlocks, yBlocks, uint_t(1),
                                                      xCells, yCells, uint_t(1),
                                                      true,
                                                      false, false, false );

   BlockDataID uId = field::addToStorage< PdeField_T >( blocks, "u", 0_r, field::zyxf, uint_t(1) );
   BlockDataID rId = field::addToStorage< PdeField_T >( blocks, "r", 0_r, field::zyxf, uint_t(1) );
   BlockDataID dId = field::addToStorage< PdeField_T >( blocks, "d", 0_r, field::zyxf, uint_t(1) );
   BlockDataID zId = field::addToStorage< PdeField_T >( blocks, "z", 0_r, field::zyxf, uint_t(1) );

   initU( blocks, uId );

   BlockDataID fId = field::addToStorage< PdeField_T >( blocks, "f", 0_r, field::zyxf, uint_t(1) );

   initF( blocks, fId );   

   SweepTimeloop timeloop( blocks, uint_t(1) );

   blockforest::communication::UniformBufferedScheme< Stencil_T > synchronizeD( blocks );
   synchronizeD.addPackInfo( make_shared< field::communication::PackInfo< PdeField_T > >( dId ) );

   std::vector< real_t > weights( Stencil_T::Size );
   weights[ Stencil_T::idx[ stencil::C ] ] = 2_r / ( blocks->dx() * blocks->dx() ) + 2_r / ( blocks->dy() * blocks->dy() ) + 4_r * math::PI * math::PI;
   weights[ Stencil_T::idx[ stencil::N ] ] = -1_r / ( blocks->dy() * blocks->dy() );
   weights[ Stencil_T::idx[ stencil::S ] ] = -1_r / ( blocks->dy() * blocks->dy() );
   weights[ Stencil_T::idx[ stencil::E ] ] = -1_r / ( blocks->dx() * blocks->dx() );
   weights[ Stencil_T::idx[ stencil::W ] ] = -1_r / ( blocks->dx() * blocks->dx() );
                                                       
   timeloop.addFuncBeforeTimeStep( pde::CGFixedStencilIteration< Stencil_T >( blocks->getBlockStorage(), uId, rId, dId, zId, fId, weights, shortrun ? uint_t(10) : uint_t(10000),
                                                                              synchronizeD, real_c(1e-6) ), "CG iteration" );

   timeloop.run();
   
   // rerun the test with a stencil field
   
   clearField<PdeField_T>( blocks, uId);
   initU( blocks, uId );
   clearField<PdeField_T>( blocks, rId );
   clearField<PdeField_T>( blocks, dId );
   clearField<PdeField_T>( blocks, zId );
   
   BlockDataID stencilId = field::addToStorage< StencilField_T >( blocks, "w" );
   
   SweepTimeloop timeloop2( blocks, uint_t(1) );
   
   copyWeightsToStencilField( blocks, weights, stencilId );
   
   timeloop2.addFuncBeforeTimeStep( pde::CGIteration< Stencil_T >( blocks->getBlockStorage(), uId, rId, dId, zId, fId, stencilId, shortrun ? uint_t(10) : uint_t(10000),
                                                                              synchronizeD, real_c(1e-6) ), "CG iteration" );
   
   timeloop2.run();

   if( !shortrun )
   {
      vtk::writeDomainDecomposition( blocks );
      field::createVTKOutput< PdeField_T >( uId, *blocks, "solution" )();
   }

   logging::Logging::printFooterOnStream();
   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}