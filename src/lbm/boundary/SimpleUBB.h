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
//! \file SimpleUBB.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/field/DensityAndVelocity.h"
#include "lbm/field/PdfField.h"

#include "boundary/Boundary.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/config/Config.h"
#include "core/debug/Debug.h"
#include "core/math/Vector3.h"

#include "field/FlagUID.h"

#include "stencil/Directions.h"

#include <vector>



namespace walberla {
namespace lbm {



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce = false >
class SimpleUBB : public Boundary<flag_t>
{
   typedef PdfField< LatticeModel_T >        PDFField;
   typedef typename LatticeModel_T::Stencil  Stencil;

public:

   static const bool threadsafe = true;

   static shared_ptr<BoundaryConfiguration> createConfiguration( const Config::BlockHandle& /*config*/ )
      { return make_shared<BoundaryConfiguration>(); }



   SimpleUBB( const BoundaryUID& boundaryUID, const FlagUID& uid, PDFField* const pdfField, const Vector3< real_t > & velocity ) :
      Boundary<flag_t>( boundaryUID ), uid_( uid ), pdfField_( pdfField ), velocity_( velocity ) { WALBERLA_ASSERT_NOT_NULLPTR( pdfField_ ); }

   SimpleUBB( const BoundaryUID& boundaryUID, const FlagUID& uid, PDFField* const pdfField, const real_t x, const real_t y, const real_t z ) :
      Boundary<flag_t>( boundaryUID ), uid_( uid ), pdfField_( pdfField ), velocity_( x, y, z ) { WALBERLA_ASSERT_NOT_NULLPTR( pdfField_ ); }

   void pushFlags( std::vector< FlagUID >& uids ) const { uids.push_back( uid_ ); }

   void beforeBoundaryTreatment() const {}
   void  afterBoundaryTreatment() const {}

   template< typename Buffer_T >
   void packCell( Buffer_T &, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const {}

   template< typename Buffer_T >
   void registerCell( Buffer_T &, const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) {}

   void registerCell( const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t, const BoundaryConfiguration& ) {}
   void registerCells( const flag_t, const CellInterval&, const BoundaryConfiguration& ) const {}
   template< typename CellIterator >
   void registerCells( const flag_t, const CellIterator&, const CellIterator&, const BoundaryConfiguration& ) const {}

   void unregisterCell( const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const {}

#ifndef NDEBUG
   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask )
#else
   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t /*mask*/ )
#endif
   {
      WALBERLA_ASSERT_EQUAL( nx, x + cell_idx_c( stencil::cx[ dir ] ) );
      WALBERLA_ASSERT_EQUAL( ny, y + cell_idx_c( stencil::cy[ dir ] ) );
      WALBERLA_ASSERT_EQUAL( nz, z + cell_idx_c( stencil::cz[ dir ] ) );
      WALBERLA_ASSERT_UNEQUAL( mask & this->mask_, numeric_cast<flag_t>(0) );
      WALBERLA_ASSERT_EQUAL( mask & this->mask_, this->mask_ ); // only true if "this->mask_" only contains one single flag, which is the case for the
                                                                // current implementation of this boundary condition (SimpleUBB)

      if( LatticeModel_T::compressible )
      {
         const auto density  = pdfField_->getDensity(x,y,z);
         const auto velocity = AdaptVelocityToExternalForce ? internal::AdaptVelocityToForce<LatticeModel_T>::get( x, y, z, pdfField_->latticeModel(), velocity_, density ) :
                                                              velocity_;

         pdfField_->get( nx, ny, nz, Stencil::invDirIdx(dir) ) = pdfField_->get( x, y, z, Stencil::idx[dir] ) -
                                                                 ( real_c(6) * density * real_c(LatticeModel_T::w[ Stencil::idx[dir] ]) *
                                                                    ( real_c(stencil::cx[ dir ]) * velocity[0] +
                                                                      real_c(stencil::cy[ dir ]) * velocity[1] +
                                                                      real_c(stencil::cz[ dir ]) * velocity[2] ) );
      }
      else
      {
         const auto velocity = AdaptVelocityToExternalForce ? internal::AdaptVelocityToForce<LatticeModel_T>::get( x, y, z, pdfField_->latticeModel(), velocity_, 1_r ) :
                                                              velocity_;

         pdfField_->get( nx, ny, nz, Stencil::invDirIdx(dir) ) = pdfField_->get( x, y, z, Stencil::idx[dir] ) -
                                                                 ( real_c(6) * real_c(LatticeModel_T::w[ Stencil::idx[dir] ]) *
                                                                    ( real_c(stencil::cx[ dir ]) * velocity[0] +
                                                                      real_c(stencil::cy[ dir ]) * velocity[1] +
                                                                      real_c(stencil::cz[ dir ]) * velocity[2] ) );
      }
   }

private:

   const FlagUID uid_;

   PDFField* const pdfField_;

   const Vector3< real_t > velocity_;

}; // class SimpleUBB



} // namespace lbm
} // namespace walberla
