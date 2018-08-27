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
//! \file CellwiseSweep.impl.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Ehsan Fattahi <ehsan.fattahi@fau.de>
//! \author Felix Winterhalter <felix.winterhalter@fau.de>
//
//======================================================================================================================

#include "lbm/field/DensityVelocityCallback.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/EquilibriumDistribution.h"
#include "lbm/lattice_model/LatticeModelBase.h"
#include "lbm/sweeps/StreamPull.h"
#include "lbm/sweeps/SweepBase.h"

#include "field/EvaluationFilter.h"
#include "field/iterators/IteratorMacros.h"

#include <boost/mpl/logical.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/int.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>


namespace walberla {
namespace lbm {

//////////////////////////
// D3Q19 SPECIALIZATION //
//////////////////////////

#define WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_MRT_1 \
   boost::mpl::and_< boost::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::MRT_tag >, \
                     boost::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >, \
                     boost::mpl::not_< boost::mpl::bool_< LatticeModel_T::compressible > >, \
                     boost::mpl::equal_to< boost::mpl::int_< LatticeModel_T::equilibriumAccuracyOrder >, boost::mpl::int_< 2 > >, \
                     boost::is_same< DensityVelocityIn_T, DefaultDensityEquilibriumVelocityCalculation > \
   >

WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_MRT_1 )

WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_MRT_1 )
{
   const auto & collisionModel = src->latticeModel().collisionModel();

   const real_t s1  = collisionModel.s1();
   const real_t s2  = collisionModel.s2();
   const real_t s4  = collisionModel.s4();
   const real_t s6  = collisionModel.s6();
   const real_t s8  = collisionModel.s8();
   const real_t s9  = collisionModel.s9();
   const real_t s10 = collisionModel.s10();
   const real_t s11 = collisionModel.s11();
   const real_t s12 = collisionModel.s12();
   const real_t s13 = collisionModel.s13();
   const real_t s14 = collisionModel.s14();
   const real_t s15 = collisionModel.s15();
   const real_t s16 = collisionModel.s16();
   const real_t s17 = collisionModel.s17();
   const real_t s18 = collisionModel.s18();

   const real_t _1_2  = 1_r / 2_r;
   const real_t _1_3  = 1_r / 3_r;
   const real_t _1_4  = 1_r / 4_r;
   const real_t _1_6  = 1_r / 6_r;
   const real_t _1_8  = 1_r / 8_r;
   const real_t _1_12 = 1_r / 12_r;
   const real_t _1_16 = 1_r / 16_r;
   const real_t _1_18 = 1_r / 18_r;
   const real_t _1_24 = 1_r / 24_r;
   const real_t _1_36 = 1_r / 36_r;
   const real_t _1_48 = 1_r / 48_r;
   const real_t _1_72 = 1_r / 72_r;

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         using namespace stencil;

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_STREAM_COLLIDE_PULL()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_INCOMP()

         const Vector3<real_t> velocity( velX, velY, velZ );
         this->densityVelocityOut( x, y, z, lm, velocity, rho + 1_r );

         const real_t velSqr = velX * velX + velY * velY + velZ * velZ;

         const real_t vel9 = 2_r * velX * velX - velY * velY - velZ * velZ;
         const real_t vel11 = velY * velY - velZ * velZ;
         const real_t vel13 = velX * velY;
         const real_t vel14 = velY * velZ;
         const real_t vel15 = velX * velZ;

         const real_t mStar0  = rho;
         const real_t mStar1  = velSqr + ( 1_r - s1 ) * ( -vC  + vNW + vNE + vSW + vSE + vTN + vTS + vTW + vTE + vBN + vBS + vBW + vBE - velSqr );
         const real_t mStar2  = ( 1_r - s2 ) * ( vC - 2_r * ( vN + vS + vW + vE + vT + vB )
                                                          + vNW + vNE + vSW + vSE + vTN + vTS + vTW + vTE + vBN + vBS + vBW + vBE );
         const real_t mStar3  = velX;
         const real_t mStar4  = ( 1_r - s4 ) * ( 2_r * vW - 2_r * vE - vNW + vNE - vSW + vSE - vTW + vTE - vBW + vBE );
         const real_t mStar5  = velY;
         const real_t mStar6  = ( 1_r - s6 ) * ( -2_r * vN + 2_r * vS + vNW + vNE - vSW - vSE + vTN - vTS + vBN - vBS );
         const real_t mStar7  = velZ;
         const real_t mStar8  = ( 1_r - s8 ) * ( -2_r * vT + 2_r * vB + vTN + vTS + vTW + vTE - vBN - vBS - vBW - vBE );
         const real_t mStar9  = vel9 + ( 1_r - s9 ) * ( -vN - vS + 2_r * vW + 2_r * vE - vT - vB + vNW + vNE + vSW + vSE - 2_r * vTN -
                                                              2_r * vTS + vTW + vTE - 2_r * vBN - 2_r * vBS + vBW + vBE - vel9 );
         const real_t mStar10 = ( 1_r - s10 ) * ( vN + vS - 2_r * vW - 2_r * vE + vT + vB + vNW + vNE + vSW + vSE - 2_r * vTN -
                                                        2_r * vTS + vTW + vTE - 2_r * vBN - 2_r * vBS + vBW + vBE );
         const real_t mStar11 = vel11 + ( 1_r - s11 ) * ( vN  + vS  - vT  - vB  + vNW + vNE + vSW + vSE - vTW - vTE - vBW - vBE - vel11 );
         const real_t mStar12 = ( 1_r - s12 ) * ( -vN - vS  + vT  + vB  + vNW + vNE + vSW + vSE - vTW - vTE - vBW - vBE );
         const real_t mStar13 = vel13 + ( 1_r - s13 ) * ( -vNW + vNE + vSW - vSE - vel13 );
         const real_t mStar14 = vel14 + ( 1_r - s14 ) * (  vTN - vTS - vBN + vBS - vel14 );
         const real_t mStar15 = vel15 + ( 1_r - s15 ) * ( -vTW + vTE + vBW - vBE - vel15 );
         const real_t mStar16 = ( 1_r - s16 ) * ( -vNW + vNE - vSW + vSE + vTW - vTE + vBW - vBE );
         const real_t mStar17 = ( 1_r - s17 ) * ( -vNW - vNE + vSW + vSE + vTN - vTS + vBN - vBS );
         const real_t mStar18 = ( 1_r - s18 ) * ( -vTN - vTS + vTW + vTE + vBN + vBS - vBW - vBE );

         dst->get( x, y, z, Stencil_T::idx[C] )  = _1_3  * mStar0  - _1_2  * mStar1  + _1_6  * mStar2;
         dst->get( x, y, z, Stencil_T::idx[N] )  = _1_18 * mStar0  - _1_18 * mStar2  + _1_6  * mStar5  - _1_6  * mStar6  - _1_24 * mStar9  +
                                                   _1_24 * mStar10 + _1_8  * mStar11 - _1_8  * mStar12;
         dst->get( x, y, z, Stencil_T::idx[S] )  = _1_18 * mStar0  - _1_18 * mStar2  - _1_6  * mStar5  + _1_6  * mStar6  - _1_24 * mStar9  +
                                                   _1_24 * mStar10 + _1_8  * mStar11 - _1_8  * mStar12;
         dst->get( x, y, z, Stencil_T::idx[W] )  = _1_18 * mStar0  - _1_18 * mStar2  - _1_6  * mStar3  + _1_6  * mStar4  + _1_12 * mStar9  -
                                                   _1_12 * mStar10;
         dst->get( x, y, z, Stencil_T::idx[E] )  = _1_18 * mStar0  - _1_18 * mStar2  + _1_6  * mStar3  - _1_6  * mStar4  + _1_12 * mStar9  -
                                                   _1_12 * mStar10;
         dst->get( x, y, z, Stencil_T::idx[T] )  = _1_18 * mStar0  - _1_18 * mStar2  + _1_6  * mStar7  - _1_6  * mStar8  - _1_24 * mStar9  +
                                                   _1_24 * mStar10 - _1_8  * mStar11 + _1_8  * mStar12;
         dst->get( x, y, z, Stencil_T::idx[B] )  = _1_18 * mStar0  - _1_18 * mStar2  - _1_6  * mStar7  + _1_6  * mStar8  - _1_24 * mStar9  +
                                                   _1_24 * mStar10 - _1_8  * mStar11 + _1_8  * mStar12;
         dst->get( x, y, z, Stencil_T::idx[NW] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar3  - _1_24 * mStar4  +
                                                   _1_12 * mStar5  + _1_24 * mStar6  + _1_48 * mStar9  + _1_48 * mStar10 + _1_16 * mStar11 +
                                                   _1_16 * mStar12 - _1_4  * mStar13 - _1_8  * mStar16 - _1_8  * mStar17;
         dst->get( x, y, z, Stencil_T::idx[NE] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar3  + _1_24 * mStar4  +
                                                   _1_12 * mStar5  + _1_24 * mStar6  + _1_48 * mStar9  + _1_48 * mStar10 + _1_16 * mStar11 +
                                                   _1_16 * mStar12 + _1_4  * mStar13 + _1_8  * mStar16 - _1_8  * mStar17;
         dst->get( x, y, z, Stencil_T::idx[SW] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar3  - _1_24 * mStar4  -
                                                   _1_12 * mStar5  - _1_24 * mStar6  + _1_48 * mStar9  + _1_48 * mStar10 + _1_16 * mStar11 +
                                                   _1_16 * mStar12 + _1_4  * mStar13 - _1_8  * mStar16 + _1_8  * mStar17;
         dst->get( x, y, z, Stencil_T::idx[SE] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar3  + _1_24 * mStar4  -
                                                   _1_12 * mStar5  - _1_24 * mStar6  + _1_48 * mStar9  + _1_48 * mStar10 + _1_16 * mStar11 +
                                                   _1_16 * mStar12 - _1_4  * mStar13 + _1_8  * mStar16 + _1_8  * mStar17;
         dst->get( x, y, z, Stencil_T::idx[TN] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar5  + _1_24 * mStar6  +
                                                   _1_12 * mStar7  + _1_24 * mStar8  - _1_24 * mStar9  - _1_24 * mStar10 + _1_4  * mStar14 +
                                                   _1_8  * mStar17 - _1_8  * mStar18;
         dst->get( x, y, z, Stencil_T::idx[TS] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar5  - _1_24 * mStar6  +
                                                   _1_12 * mStar7  + _1_24 * mStar8  - _1_24 * mStar9  - _1_24 * mStar10 - _1_4  * mStar14 -
                                                   _1_8  * mStar17 - _1_8  * mStar18;
         dst->get( x, y, z, Stencil_T::idx[TW] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar3  - _1_24 * mStar4  +
                                                   _1_12 * mStar7  + _1_24 * mStar8  + _1_48 * mStar9  + _1_48 * mStar10 - _1_16 * mStar11 -
                                                   _1_16 * mStar12 - _1_4  * mStar15 + _1_8  * mStar16 + _1_8  * mStar18;
         dst->get( x, y, z, Stencil_T::idx[TE] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar3  + _1_24 * mStar4  +
                                                   _1_12 * mStar7  + _1_24 * mStar8  + _1_48 * mStar9  + _1_48 * mStar10 - _1_16 * mStar11 -
                                                   _1_16 * mStar12 + _1_4  * mStar15 - _1_8  * mStar16 + _1_8  * mStar18;
         dst->get( x, y, z, Stencil_T::idx[BN] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar5  + _1_24 * mStar6  -
                                                   _1_12 * mStar7  - _1_24 * mStar8  - _1_24 * mStar9  - _1_24 * mStar10 - _1_4  * mStar14 +
                                                   _1_8  * mStar17 + _1_8  * mStar18;
         dst->get( x, y, z, Stencil_T::idx[BS] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar5  - _1_24 * mStar6  -
                                                   _1_12 * mStar7  - _1_24 * mStar8  - _1_24 * mStar9  - _1_24 * mStar10 + _1_4  * mStar14 -
                                                   _1_8  * mStar17 + _1_8  * mStar18;
         dst->get( x, y, z, Stencil_T::idx[BW] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar3  - _1_24 * mStar4  -
                                                   _1_12 * mStar7  - _1_24 * mStar8  + _1_48 * mStar9  + _1_48 * mStar10 - _1_16 * mStar11 -
                                                   _1_16 * mStar12 + _1_4  * mStar15 + _1_8  * mStar16 - _1_8  * mStar18;
         dst->get( x, y, z, Stencil_T::idx[BE] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar3  + _1_24 * mStar4  -
                                                   _1_12 * mStar7  - _1_24 * mStar8  + _1_48 * mStar9  + _1_48 * mStar10 - _1_16 * mStar11 -
                                                   _1_16 * mStar12 - _1_4  * mStar15 - _1_8  * mStar16 - _1_8  * mStar18;
         
         if (boost::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value == false)
         {
            const auto commonForceTerms = lm.forceModel().template directionIndependentTerms< LatticeModel_T >( x, y, z, velocity, rho + 1.0_r, collisionModel.omega(), collisionModel.omega_bulk() );
            for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
               dst->get( x, y, z, d.toIdx() ) += lm.forceModel().template forceTerm< LatticeModel_T >( x, y, z, velocity, rho + 1.0_r, commonForceTerms, LatticeModel_T::w[ d.toIdx() ], real_c(d.cx()), real_c(d.cy()), real_c(d.cz()), collisionModel.omega(), collisionModel.omega_bulk() );
         }
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT()

WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_MRT_1 )
{
   const auto & collisionModel = src->latticeModel().collisionModel();

   const real_t s1  = collisionModel.s1();
   const real_t s2  = collisionModel.s2();
   const real_t s4  = collisionModel.s4();
   const real_t s6  = collisionModel.s6();
   const real_t s8  = collisionModel.s8();
   const real_t s9  = collisionModel.s9();
   const real_t s10 = collisionModel.s10();
   const real_t s11 = collisionModel.s11();
   const real_t s12 = collisionModel.s12();
   const real_t s13 = collisionModel.s13();
   const real_t s14 = collisionModel.s14();
   const real_t s15 = collisionModel.s15();
   const real_t s16 = collisionModel.s16();
   const real_t s17 = collisionModel.s17();
   const real_t s18 = collisionModel.s18();

   const real_t _1_2  = 1_r / 2_r;
   const real_t _1_3  = 1_r / 3_r;
   const real_t _1_4  = 1_r / 4_r;
   const real_t _1_6  = 1_r / 6_r;
   const real_t _1_8  = 1_r / 8_r;
   const real_t _1_12 = 1_r / 12_r;
   const real_t _1_16 = 1_r / 16_r;
   const real_t _1_18 = 1_r / 18_r;
   const real_t _1_24 = 1_r / 24_r;
   const real_t _1_36 = 1_r / 36_r;
   const real_t _1_48 = 1_r / 48_r;
   const real_t _1_72 = 1_r / 72_r;

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         using namespace stencil;

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_COLLIDE_GET()

         WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_INCOMP()
         
         const Vector3<real_t> velocity( velX, velY, velZ );
         this->densityVelocityOut( x, y, z, lm, velocity, rho + 1_r );

         const real_t velSqr = velX * velX + velY * velY + velZ * velZ;

         const real_t vel9 = 2_r * velX * velX - velY * velY - velZ * velZ;
         const real_t vel11 = velY * velY - velZ * velZ;
         const real_t vel13 = velX * velY;
         const real_t vel14 = velY * velZ;
         const real_t vel15 = velX * velZ;

         const real_t mStar0  = rho;
         const real_t mStar1  = velSqr + ( 1_r - s1 ) * ( -vC  + vNW + vNE + vSW + vSE + vTN + vTS + vTW + vTE + vBN + vBS + vBW + vBE - velSqr );
         const real_t mStar2  = ( 1_r - s2 ) * ( vC - 2_r * ( vN + vS + vW + vE + vT + vB )
                                                          + vNW + vNE + vSW + vSE + vTN + vTS + vTW + vTE + vBN + vBS + vBW + vBE );
         const real_t mStar3  = velX;
         const real_t mStar4  = ( 1_r - s4 ) * ( 2_r * vW - 2_r * vE - vNW + vNE - vSW + vSE - vTW + vTE - vBW + vBE );
         const real_t mStar5  = velY;
         const real_t mStar6  = ( 1_r - s6 ) * ( -2_r * vN + 2_r * vS + vNW + vNE - vSW - vSE + vTN - vTS + vBN - vBS );
         const real_t mStar7  = velZ;
         const real_t mStar8  = ( 1_r - s8 ) * ( -2_r * vT + 2_r * vB + vTN + vTS + vTW + vTE - vBN - vBS - vBW - vBE );
         const real_t mStar9  = vel9 + ( 1_r - s9 ) * ( -vN - vS + 2_r * vW + 2_r * vE - vT - vB + vNW + vNE + vSW + vSE - 2_r * vTN -
                                                              2_r * vTS + vTW + vTE - 2_r * vBN - 2_r * vBS + vBW + vBE - vel9 );
         const real_t mStar10 = ( 1_r - s10 ) * ( vN + vS - 2_r * vW - 2_r * vE + vT + vB + vNW + vNE + vSW + vSE - 2_r * vTN -
                                                        2_r * vTS + vTW + vTE - 2_r * vBN - 2_r * vBS + vBW + vBE );
         const real_t mStar11 = vel11 + ( 1_r - s11 ) * ( vN  + vS  - vT  - vB  + vNW + vNE + vSW + vSE - vTW - vTE - vBW - vBE - vel11 );
         const real_t mStar12 = ( 1_r - s12 ) * ( -vN - vS  + vT  + vB  + vNW + vNE + vSW + vSE - vTW - vTE - vBW - vBE );
         const real_t mStar13 = vel13 + ( 1_r - s13 ) * ( -vNW + vNE + vSW - vSE - vel13 );
         const real_t mStar14 = vel14 + ( 1_r - s14 ) * (  vTN - vTS - vBN + vBS - vel14 );
         const real_t mStar15 = vel15 + ( 1_r - s15 ) * ( -vTW + vTE + vBW - vBE - vel15 );
         const real_t mStar16 = ( 1_r - s16 ) * ( -vNW + vNE - vSW + vSE + vTW - vTE + vBW - vBE );
         const real_t mStar17 = ( 1_r - s17 ) * ( -vNW - vNE + vSW + vSE + vTN - vTS + vBN - vBS );
         const real_t mStar18 = ( 1_r - s18 ) * ( -vTN - vTS + vTW + vTE + vBN + vBS - vBW - vBE );

         src->get( x, y, z, Stencil_T::idx[C] )  = _1_3  * mStar0  - _1_2  * mStar1  + _1_6  * mStar2;
         src->get( x, y, z, Stencil_T::idx[N] )  = _1_18 * mStar0  - _1_18 * mStar2  + _1_6  * mStar5  - _1_6  * mStar6  - _1_24 * mStar9  +
                                                   _1_24 * mStar10 + _1_8  * mStar11 - _1_8  * mStar12;
         src->get( x, y, z, Stencil_T::idx[S] )  = _1_18 * mStar0  - _1_18 * mStar2  - _1_6  * mStar5  + _1_6  * mStar6  - _1_24 * mStar9  +
                                                   _1_24 * mStar10 + _1_8  * mStar11 - _1_8  * mStar12;
         src->get( x, y, z, Stencil_T::idx[W] )  = _1_18 * mStar0  - _1_18 * mStar2  - _1_6  * mStar3  + _1_6  * mStar4  + _1_12 * mStar9  -
                                                   _1_12 * mStar10;
         src->get( x, y, z, Stencil_T::idx[E] )  = _1_18 * mStar0  - _1_18 * mStar2  + _1_6  * mStar3  - _1_6  * mStar4  + _1_12 * mStar9  -
                                                   _1_12 * mStar10;
         src->get( x, y, z, Stencil_T::idx[T] )  = _1_18 * mStar0  - _1_18 * mStar2  + _1_6  * mStar7  - _1_6  * mStar8  - _1_24 * mStar9  +
                                                   _1_24 * mStar10 - _1_8  * mStar11 + _1_8  * mStar12;
         src->get( x, y, z, Stencil_T::idx[B] )  = _1_18 * mStar0  - _1_18 * mStar2  - _1_6  * mStar7  + _1_6  * mStar8  - _1_24 * mStar9  +
                                                   _1_24 * mStar10 - _1_8  * mStar11 + _1_8  * mStar12;
         src->get( x, y, z, Stencil_T::idx[NW] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar3  - _1_24 * mStar4  +
                                                   _1_12 * mStar5  + _1_24 * mStar6  + _1_48 * mStar9  + _1_48 * mStar10 + _1_16 * mStar11 +
                                                   _1_16 * mStar12 - _1_4  * mStar13 - _1_8  * mStar16 - _1_8  * mStar17;
         src->get( x, y, z, Stencil_T::idx[NE] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar3  + _1_24 * mStar4  +
                                                   _1_12 * mStar5  + _1_24 * mStar6  + _1_48 * mStar9  + _1_48 * mStar10 + _1_16 * mStar11 +
                                                   _1_16 * mStar12 + _1_4  * mStar13 + _1_8  * mStar16 - _1_8  * mStar17;
         src->get( x, y, z, Stencil_T::idx[SW] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar3  - _1_24 * mStar4  -
                                                   _1_12 * mStar5  - _1_24 * mStar6  + _1_48 * mStar9  + _1_48 * mStar10 + _1_16 * mStar11 +
                                                   _1_16 * mStar12 + _1_4  * mStar13 - _1_8  * mStar16 + _1_8  * mStar17;
         src->get( x, y, z, Stencil_T::idx[SE] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar3  + _1_24 * mStar4  -
                                                   _1_12 * mStar5  - _1_24 * mStar6  + _1_48 * mStar9  + _1_48 * mStar10 + _1_16 * mStar11 +
                                                   _1_16 * mStar12 - _1_4  * mStar13 + _1_8  * mStar16 + _1_8  * mStar17;
         src->get( x, y, z, Stencil_T::idx[TN] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar5  + _1_24 * mStar6  +
                                                   _1_12 * mStar7  + _1_24 * mStar8  - _1_24 * mStar9  - _1_24 * mStar10 + _1_4  * mStar14 +
                                                   _1_8  * mStar17 - _1_8  * mStar18;
         src->get( x, y, z, Stencil_T::idx[TS] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar5  - _1_24 * mStar6  +
                                                   _1_12 * mStar7  + _1_24 * mStar8  - _1_24 * mStar9  - _1_24 * mStar10 - _1_4  * mStar14 -
                                                   _1_8  * mStar17 - _1_8  * mStar18;
         src->get( x, y, z, Stencil_T::idx[TW] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar3  - _1_24 * mStar4  +
                                                   _1_12 * mStar7  + _1_24 * mStar8  + _1_48 * mStar9  + _1_48 * mStar10 - _1_16 * mStar11 -
                                                   _1_16 * mStar12 - _1_4  * mStar15 + _1_8  * mStar16 + _1_8  * mStar18;
         src->get( x, y, z, Stencil_T::idx[TE] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar3  + _1_24 * mStar4  +
                                                   _1_12 * mStar7  + _1_24 * mStar8  + _1_48 * mStar9  + _1_48 * mStar10 - _1_16 * mStar11 -
                                                   _1_16 * mStar12 + _1_4  * mStar15 - _1_8  * mStar16 + _1_8  * mStar18;
         src->get( x, y, z, Stencil_T::idx[BN] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar5  + _1_24 * mStar6  -
                                                   _1_12 * mStar7  - _1_24 * mStar8  - _1_24 * mStar9  - _1_24 * mStar10 - _1_4  * mStar14 +
                                                   _1_8  * mStar17 + _1_8  * mStar18;
         src->get( x, y, z, Stencil_T::idx[BS] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar5  - _1_24 * mStar6  -
                                                   _1_12 * mStar7  - _1_24 * mStar8  - _1_24 * mStar9  - _1_24 * mStar10 + _1_4  * mStar14 -
                                                   _1_8  * mStar17 + _1_8  * mStar18;
         src->get( x, y, z, Stencil_T::idx[BW] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  - _1_12 * mStar3  - _1_24 * mStar4  -
                                                   _1_12 * mStar7  - _1_24 * mStar8  + _1_48 * mStar9  + _1_48 * mStar10 - _1_16 * mStar11 -
                                                   _1_16 * mStar12 + _1_4  * mStar15 + _1_8  * mStar16 - _1_8  * mStar18;
         src->get( x, y, z, Stencil_T::idx[BE] ) = _1_36 * mStar0  + _1_24 * mStar1  + _1_72 * mStar2  + _1_12 * mStar3  + _1_24 * mStar4  -
                                                   _1_12 * mStar7  - _1_24 * mStar8  + _1_48 * mStar9  + _1_48 * mStar10 - _1_16 * mStar11 -
                                                   _1_16 * mStar12 - _1_4  * mStar15 - _1_8  * mStar16 - _1_8  * mStar18;
         
         if (boost::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value == false)
         {
            const auto commonForceTerms = lm.forceModel().template directionIndependentTerms< LatticeModel_T >( x, y, z, velocity, rho + 1.0_r, collisionModel.omega(), collisionModel.omega_bulk() );
            for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
               src->get( x, y, z, d.toIdx() ) += lm.forceModel().template forceTerm< LatticeModel_T >( x, y, z, velocity, rho + 1.0_r, commonForceTerms, LatticeModel_T::w[ d.toIdx() ], real_c(d.cx()), real_c(d.cy()), real_c(d.cz()), collisionModel.omega(), collisionModel.omega_bulk() );
         }
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
}
WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT()

#undef WALBERLA_LBM_CELLWISE_SWEEP_SPECIALIZATION_MRT_1



} // namespace lbm
} // namespace walberla
