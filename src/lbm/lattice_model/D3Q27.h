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
//! \file D3Q27.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "LatticeModelBase.h"
#include "stencil/D3Q27.h"

#include <boost/mpl/not.hpp>


namespace walberla {
namespace lbm {



template< typename CollisionModel_T, bool Compressible = false, typename ForceModel_T = force_model::None, int EquilibriumAccuracyOrder = 2 >
class D3Q27 : public LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >
{
public:

   static_assert( (boost::mpl::not_< boost::is_same< CollisionModel_T, collision_model::D3Q19MRT > >::value), "D3Q19MRT only works with D3Q19!" );

   typedef typename LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >::CollisionModel  CollisionModel;
   typedef typename LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >::ForceModel      ForceModel;

   typedef stencil::D3Q27 Stencil;
   typedef stencil::D3Q27 CommunicationStencil;

   static const char * NAME;

   static const real_t w_0;
   static const real_t w_1;
   static const real_t w_2;
   static const real_t w_3;
   static const real_t w[27];    // Stencil::Size !
   static const real_t wInv[27]; // Stencil::Size !

   D3Q27( const CollisionModel_T & cm, const ForceModel_T & fm  ) :
      LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >( cm, fm ) {}

   // available only if the force model == force_model::None
   D3Q27( const CollisionModel_T & cm ) :
      LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >( cm, force_model::None() )
   {
      static_assert( (boost::is_same< ForceModel_T, force_model::None >::value), "This constructor is only available if the force model is equal to force_model::None!" );
   }

   virtual ~D3Q27() {}

protected:

   virtual void config( IBlock & /*block*/, StructuredBlockStorage & /*sbs*/ ) {}
};

template< typename CM, bool C, typename FM, int EAO > const char*  D3Q27<CM,C,FM,EAO>::NAME = "D3Q27";

template< typename CM, bool C, typename FM, int EAO > const real_t D3Q27<CM,C,FM,EAO>::w_0 = 8.0_r / real_t( 27.0);
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q27<CM,C,FM,EAO>::w_1 = 2.0_r / real_t( 27.0);
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q27<CM,C,FM,EAO>::w_2 = 1.0_r / real_t( 54.0);
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q27<CM,C,FM,EAO>::w_3 = 1.0_r / 216.0_r;

// must match with the static array 'dir' in stencil::D3Q27
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q27<CM,C,FM,EAO>::w[27] = { 8.0_r / real_t( 27.0),   // C
                                                                                                 2.0_r / real_t( 27.0),   // N
                                                                                                 2.0_r / real_t( 27.0),   // S
                                                                                                 2.0_r / real_t( 27.0),   // W
                                                                                                 2.0_r / real_t( 27.0),   // E
                                                                                                 2.0_r / real_t( 27.0),   // T
                                                                                                 2.0_r / real_t( 27.0),   // B
                                                                                                 1.0_r / real_t( 54.0),   // NW
                                                                                                 1.0_r / real_t( 54.0),   // NE
                                                                                                 1.0_r / real_t( 54.0),   // SW
                                                                                                 1.0_r / real_t( 54.0),   // SE
                                                                                                 1.0_r / real_t( 54.0),   // TN
                                                                                                 1.0_r / real_t( 54.0),   // TS
                                                                                                 1.0_r / real_t( 54.0),   // TW
                                                                                                 1.0_r / real_t( 54.0),   // TE
                                                                                                 1.0_r / real_t( 54.0),   // BN
                                                                                                 1.0_r / real_t( 54.0),   // BS
                                                                                                 1.0_r / real_t( 54.0),   // BW
                                                                                                 1.0_r / real_t( 54.0),   // BE
                                                                                                 1.0_r / 216.0_r,   // TNE
                                                                                                 1.0_r / 216.0_r,   // TNW
                                                                                                 1.0_r / 216.0_r,   // TSE
                                                                                                 1.0_r / 216.0_r,   // TSW
                                                                                                 1.0_r / 216.0_r,   // BNE
                                                                                                 1.0_r / 216.0_r,   // BNW
                                                                                                 1.0_r / 216.0_r,   // BSE
                                                                                                 1.0_r / 216.0_r }; // BSW

// must match with the static array 'dir' in stencil::D3Q27
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q27<CM,C,FM,EAO>::wInv[27] = { 27.0_r / real_t( 8.0),   // C
                                                                                                    27.0_r / real_t( 2.0),   // N
                                                                                                    27.0_r / real_t( 2.0),   // S
                                                                                                    27.0_r / real_t( 2.0),   // W
                                                                                                    27.0_r / real_t( 2.0),   // E
                                                                                                    27.0_r / real_t( 2.0),   // T
                                                                                                    27.0_r / real_t( 2.0),   // B
                                                                                                    real_t( 54.0),                 // NW
                                                                                                    real_t( 54.0),                 // NE
                                                                                                    real_t( 54.0),                 // SW
                                                                                                    real_t( 54.0),                 // SE
                                                                                                    real_t( 54.0),                 // TN
                                                                                                    real_t( 54.0),                 // TS
                                                                                                    real_t( 54.0),                 // TW
                                                                                                    real_t( 54.0),                 // TE
                                                                                                    real_t( 54.0),                 // BN
                                                                                                    real_t( 54.0),                 // BS
                                                                                                    real_t( 54.0),                 // BW
                                                                                                    real_t( 54.0),                 // BE
                                                                                                    216.0_r,                 // TNE
                                                                                                    216.0_r,                 // TNW
                                                                                                    216.0_r,                 // TSE
                                                                                                    216.0_r,                 // TSW
                                                                                                    216.0_r,                 // BNE
                                                                                                    216.0_r,                 // BNW
                                                                                                    216.0_r,                 // BSE
                                                                                                    216.0_r };               // BSW

} // namespace lbm
} // namespace walberla
