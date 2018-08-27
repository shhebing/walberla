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
//! \file D3Q19.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "LatticeModelBase.h"
#include "stencil/D3Q19.h"

#include <boost/type_traits/is_same.hpp>


namespace walberla {
namespace lbm {



template< typename CollisionModel_T, bool Compressible = false, typename ForceModel_T = force_model::None, int EquilibriumAccuracyOrder = 2 >
class D3Q19 : public LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >
{
public:

   typedef typename LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >::CollisionModel  CollisionModel;
   typedef typename LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >::ForceModel      ForceModel;

   typedef stencil::D3Q19 Stencil;
   typedef stencil::D3Q19 CommunicationStencil;

   static const char * NAME;

   static const real_t w_0;
   static const real_t w_1;
   static const real_t w_2;
   static const real_t w[19];    // Stencil::Size !
   static const real_t wInv[19]; // Stencil::Size !

   D3Q19( const CollisionModel_T & cm, const ForceModel_T & fm ) :
      LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >( cm, fm ) {}

   // available only if the force model == force_model::None
   D3Q19( const CollisionModel_T & cm ) :
      LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >( cm, force_model::None() )
   {
      static_assert( (boost::is_same< ForceModel_T, force_model::None >::value), "This constructor is only available if the force model is equal to force_model::None!" );
   }

   virtual ~D3Q19() {}

protected:

   virtual void config( IBlock & /*block*/, StructuredBlockStorage & /*sbs*/ ) {}
};

template< typename CM, bool C, typename FM, int EAO > const char*  D3Q19<CM,C,FM,EAO>::NAME = "D3Q19";


template< typename CM, bool C, typename FM, int EAO > const real_t D3Q19<CM,C,FM,EAO>::w_0 = 1.0_r / real_t( 3.0);
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q19<CM,C,FM,EAO>::w_1 = 1.0_r / 18.0_r;
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q19<CM,C,FM,EAO>::w_2 = 1.0_r / 36.0_r;

// must match with the static array 'dir' in stencil::D3Q19
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q19<CM,C,FM,EAO>::w[19] = { 1.0_r / real_t( 3.0),   // C
                                                                                                 1.0_r / 18.0_r,   // N
                                                                                                 1.0_r / 18.0_r,   // S
                                                                                                 1.0_r / 18.0_r,   // W
                                                                                                 1.0_r / 18.0_r,   // E
                                                                                                 1.0_r / 18.0_r,   // T
                                                                                                 1.0_r / 18.0_r,   // B
                                                                                                 1.0_r / 36.0_r,   // NW
                                                                                                 1.0_r / 36.0_r,   // NE
                                                                                                 1.0_r / 36.0_r,   // SW
                                                                                                 1.0_r / 36.0_r,   // SE
                                                                                                 1.0_r / 36.0_r,   // TN
                                                                                                 1.0_r / 36.0_r,   // TS
                                                                                                 1.0_r / 36.0_r,   // TW
                                                                                                 1.0_r / 36.0_r,   // TE
                                                                                                 1.0_r / 36.0_r,   // BN
                                                                                                 1.0_r / 36.0_r,   // BS
                                                                                                 1.0_r / 36.0_r,   // BW
                                                                                                 1.0_r / 36.0_r }; // BE

// must match with the static array 'dir' in stencil::D3Q19
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q19<CM,C,FM,EAO>::wInv[19] = { real_t( 3.0),   // C
                                                                                                    18.0_r,   // N
                                                                                                    18.0_r,   // S
                                                                                                    18.0_r,   // W
                                                                                                    18.0_r,   // E
                                                                                                    18.0_r,   // T
                                                                                                    18.0_r,   // B
                                                                                                    36.0_r,   // NW
                                                                                                    36.0_r,   // NE
                                                                                                    36.0_r,   // SW
                                                                                                    36.0_r,   // SE
                                                                                                    36.0_r,   // TN
                                                                                                    36.0_r,   // TS
                                                                                                    36.0_r,   // TW
                                                                                                    36.0_r,   // TE
                                                                                                    36.0_r,   // BN
                                                                                                    36.0_r,   // BS
                                                                                                    36.0_r,   // BW
                                                                                                    36.0_r }; // BE



} // namespace lbm
} // namespace walberla
