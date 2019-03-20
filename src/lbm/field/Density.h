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
//! \file Density.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

#include <boost/mpl/bool.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/utility/enable_if.hpp>


// Back-end for calculating macroscopic values
// You should never use these functions directly, always refer to the member functions
// of PdfField or the free functions that can be found in MacroscopicValueCalculation.h

namespace walberla {
namespace lbm {



template< typename LatticeModel_T, class Enable = void >
struct Density;

//////////////////
// Compressible //
//////////////////

template< typename LatticeModel_T >
struct Density< LatticeModel_T, typename boost::enable_if_c< LatticeModel_T::compressible >::type >
{
   static_assert( LatticeModel_T::compressible, "Only works with compressible models!" );

   template< typename FieldPtrOrIterator >
   static inline real_t get( const LatticeModel_T & /*latticeModel*/, const FieldPtrOrIterator & it )
   {
      real_t rho = it[0];
      for( uint_t i = 1; i != LatticeModel_T::Stencil::Size; ++i )
         rho += it[i];
      return rho;

      //real_t rho = 0.0_r;
      //for( auto i = LatticeModel_T::Stencil::begin(); i != LatticeModel_T::Stencil::end(); ++i )
      //   rho += it[ i.toIdx() ];
      //return rho;
   }

   template< typename PdfField_T >
   static inline real_t get( const LatticeModel_T & /*latticeModel*/,
                             const PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
      const real_t & xyz0 = pdf(x,y,z,0);
      real_t rho = xyz0;
      for( uint_t i = 1; i != LatticeModel_T::Stencil::Size; ++i )
         rho += pdf.getF( &xyz0, i );
      return rho;
   }
};

////////////////////
// Incompressible //
////////////////////

template< typename LatticeModel_T >
struct Density< LatticeModel_T, typename boost::enable_if< boost::mpl::not_< boost::mpl::bool_< LatticeModel_T::compressible > > >::type >
{
   static_assert( LatticeModel_T::compressible == false, "Only works with incompressible models!" );

   template< typename FieldPtrOrIterator >
   static inline real_t get( const LatticeModel_T & /*latticeModel*/, const FieldPtrOrIterator & it )
   {
      real_t rho = it[0] + 1.0_r;
      for( uint_t i = 1; i != LatticeModel_T::Stencil::Size; ++i )
         rho += it[i];
      return rho;

      //real_t rho = 1.0_r;
      //for( auto i = LatticeModel_T::Stencil::begin(); i != LatticeModel_T::Stencil::end(); ++i )
      //   rho += it[ i.toIdx() ];
      //return rho;
   }

   template< typename PdfField_T >
   static inline real_t get( const LatticeModel_T & /*latticeModel*/,
                             const PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
      const real_t & xyz0 = pdf(x,y,z,0);
      real_t rho = xyz0 + 1.0_r;
      for( uint_t i = 1; i != LatticeModel_T::Stencil::Size; ++i )
         rho += pdf.getF( &xyz0, i );
      return rho;
   }
};



} // namespace lbm
} // namespace walberla
