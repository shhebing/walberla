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
//! \file DefaultCellOperation.impl.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "lbm/field/MacroscopicValueCalculation.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/EquilibriumDistribution.h"
#include "lbm/lattice_model/LatticeModelBase.h"

#include <boost/mpl/logical.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>


namespace walberla {
namespace lbm {



////////////////////////////////////////////////////////////////////////
// Available SRT implementations:                                     //
//                                                                    //
// Generic (D*Q*) version:                                            //
//                                      incompressible | compressible //
//                           no forces:       x               x       //
//                                                                    //
// Optimized D3Q19 implementation:                                    //
//                                      incompressible | compressible //
//                           no forces:       x               x       //
////////////////////////////////////////////////////////////////////////


///////////////////////////////
// Specialization for:       //
// - no additional forces    //
///////////////////////////////

template< typename LatticeModel_T >
class DefaultCellOperation< LatticeModel_T, typename boost::enable_if< boost::mpl::and_< boost::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                                                         collision_model::SRT_tag >,
                                                                                         boost::mpl::bool_< LatticeModel_T::CollisionModel::constant >,
                                                                                         boost::mpl::not_< boost::is_same< typename LatticeModel_T::Stencil,
                                                                                                                  stencil::D3Q19 > >,
                                                                                         boost::is_same< typename LatticeModel_T::ForceModel::tag,
                                                                                                         force_model::None_tag > > >::type >
{
public:

   static_assert( (boost::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( (boost::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value == false),              "There is a specialization for D3Q19!" );
   static_assert( (boost::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   typedef PdfField< LatticeModel_T >        PdfField_T;
   typedef typename LatticeModel_T::Stencil  Stencil;

   DefaultCellOperation() : omega_( 0_r ), latticeModel_( NULL ) {}

   void configure( const LatticeModel_T & latticeModel )
   {
      omega_ = latticeModel.collisionModel().omega();
      latticeModel_ = &latticeModel;
   }

   void operator()( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const
   {
      // stream pull
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
         dst->get( x,y,z,d.toIdx() ) = src->get( x-d.cx(), y-d.cy(), z-d.cz(), d.toIdx() );

      Vector3<real_t> velocity;
      real_t rho = dst->getDensityAndVelocity( velocity, x, y, z );

      // collide
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
      {
         dst->get( x, y, z, d.toIdx() ) = ( 1.0_r - omega_ ) * dst->get( x, y, z, d.toIdx() ) +
                                                          omega_   * EquilibriumDistribution< LatticeModel_T >::get( *d, velocity, rho );
      }
   }

   template< typename FieldPtrOrIterator >
   void operator()( FieldPtrOrIterator & src, FieldPtrOrIterator & dst ) const
   {
      // stream pull
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
         dst[ d.toIdx() ] = src.neighbor( d.inverseDir(), d.toIdx() );

      Vector3<real_t> velocity;
      real_t rho = getDensityAndVelocity( velocity, *latticeModel_, dst );

      // collide
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
      {
         dst[ d.toIdx() ] = ( 1.0_r - omega_ ) * dst[ d.toIdx() ] +
                                            omega_   * EquilibriumDistribution< LatticeModel_T >::get( *d, velocity, rho );
      }
   }
private:

   real_t omega_;
   const LatticeModel_T * latticeModel_;
};



/////////////////////////////////////////////
// Specialization for:                     //
// - compressible                          //
// - additional forces (Guo - constant)    //
/////////////////////////////////////////////

template< typename LatticeModel_T >
class DefaultCellOperation< LatticeModel_T, typename boost::enable_if< boost::mpl::and_< boost::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                                                         collision_model::SRT_tag >,
                                                                                         boost::mpl::bool_< LatticeModel_T::CollisionModel::constant >,
                                                                                         boost::is_same< typename LatticeModel_T::ForceModel::tag,
                                                                                                         force_model::Guo_tag > > >::type >
{                                                                                   
public:

   static_assert( (boost::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( (boost::is_same< typename LatticeModel_T::ForceModel::tag, force_model::Guo_tag >::value),         "Only works with Guo constant force model !" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   typedef PdfField< LatticeModel_T >        PdfField_T;
   typedef typename LatticeModel_T::Stencil  Stencil;

   DefaultCellOperation() : omega_( 0_r ), latticeModel_( NULL ) {}

   void configure( const LatticeModel_T & latticeModel )
   {
      omega_ = latticeModel.collisionModel().omega();
      latticeModel_ = &latticeModel;
      force_ = latticeModel.forceModel().force();
   }

   void operator()( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const
   {
      // stream pull
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
         dst->get( x,y,z,d.toIdx() ) = src->get( x-d.cx(), y-d.cy(), z-d.cz(), d.toIdx() );

      Vector3<real_t> velocity;
      real_t rho = dst->getDensityAndVelocity( velocity, x, y, z );

      // collide
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
      {
         const Vector3<real_t> c( real_c(d.cx()), real_c(d.cy()), real_c(d.cz()) );

         const real_t force_trm = 3.0_r * LatticeModel_T::w[ d.toIdx() ] * ( 1_r - 0.5_r * omega_ ) *
                                  ( ( c - velocity + ( 3_r * ( c * velocity ) * c ) ) * force_ );

         dst->get( x, y, z, d.toIdx() ) = ( 1.0_r - omega_ ) * dst->get( x, y, z, d.toIdx() ) +
                                                          omega_   * EquilibriumDistribution< LatticeModel_T >::get( *d, velocity, rho ) +
                                          force_trm;
      }
   }

   template< typename FieldPtrOrIterator >
   void operator()( FieldPtrOrIterator & src, FieldPtrOrIterator & dst ) const
   {
      // stream pull
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
         dst[ d.toIdx() ] = src.neighbor( d.inverseDir(), d.toIdx() );

      Vector3<real_t> velocity;
      real_t rho = getDensityAndVelocity( velocity, *latticeModel_, dst );

      // collide< LatticeModel_T::CollisionModel::constant >
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
      {
         const Vector3<real_t> c( real_c(d.cx()), real_c(d.cy()), real_c(d.cz()) );

         const real_t force_trm = 3.0_r * LatticeModel_T::w[ d.toIdx() ] * ( 1_r - 0.5_r * omega_ ) *
                                  ( ( c - velocity + ( 3_r * ( c * velocity ) * c ) ) * force_ );

         dst[ d.toIdx() ] = ( 1.0_r - omega_ ) * dst[ d.toIdx() ] +
                                            omega_   * EquilibriumDistribution< LatticeModel_T >::get( *d, velocity, rho ) +
                              force_trm;
      }
   }

private:

   real_t omega_;
   const LatticeModel_T * latticeModel_;
   Vector3<real_t> force_;
};







///////////////////////////////
// Specialization for D3Q19: //
// - incompressible          //
// - no additional forces    //
///////////////////////////////

template< typename LatticeModel_T >
class DefaultCellOperation< LatticeModel_T, typename boost::enable_if< boost::mpl::and_< boost::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                                                         collision_model::SRT_tag >,
                                                                                         boost::mpl::bool_< LatticeModel_T::CollisionModel::constant >,
                                                                                         boost::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >,
                                                                                         boost::mpl::not_< boost::mpl::bool_< LatticeModel_T::compressible > >,
                                                                                         boost::is_same< typename LatticeModel_T::ForceModel::tag,
                                                                                                         force_model::None_tag > > >::type >
{
public:

   static_assert( (boost::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( (boost::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value),                       "Only works with D3Q19!" );
   static_assert( LatticeModel_T::compressible == false,                                                             "Only works with incompressible models!" );
   static_assert( (boost::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   typedef PdfField< LatticeModel_T >        PdfField_T;
   typedef typename LatticeModel_T::Stencil  Stencil;

   DefaultCellOperation() : omega_trm_( 0_r ), omega_w0_( 0_r ), omega_w1_( 0_r ), omega_w2_( 0_r ) {}

   void configure( const LatticeModel_T & latticeModel )
   {
      const real_t omega = latticeModel.collisionModel().omega();
      omega_trm_ = 1_r - omega;
      omega_w0_  = 3_r * ( 1_r / real_t( 3) ) * omega;
      omega_w1_  = 3_r * ( 1_r / 18_r ) * omega;
      omega_w2_  = 3_r * ( 1_r / 36_r ) * omega;
   }

   void operator()( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const;

   template< typename FieldPtrOrIterator >
   void operator()( FieldPtrOrIterator & src, FieldPtrOrIterator & dst ) const;

private:

   real_t omega_trm_;
   real_t omega_w0_;
   real_t omega_w1_;
   real_t omega_w2_;
};

template< typename LatticeModel_T >
void DefaultCellOperation< LatticeModel_T, typename boost::enable_if< boost::mpl::and_< boost::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                                                        collision_model::SRT_tag >,
                                                                                        boost::mpl::bool_< LatticeModel_T::CollisionModel::constant >,
                                                                                        boost::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >,
                                                                                        boost::mpl::not_< boost::mpl::bool_< LatticeModel_T::compressible > >,
                                                                                        boost::is_same< typename LatticeModel_T::ForceModel::tag,
                                                                                                        force_model::None_tag > > >::type
   >::operator()( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const
{
   using namespace stencil;

   const real_t dd_tmp_NE = src->get(x-1, y-1, z  , Stencil::idx[NE]);
   const real_t dd_tmp_N  = src->get(x  , y-1, z  , Stencil::idx[N]);
   const real_t dd_tmp_NW = src->get(x+1, y-1, z  , Stencil::idx[NW]);
   const real_t dd_tmp_W  = src->get(x+1, y  , z  , Stencil::idx[W]);
   const real_t dd_tmp_SW = src->get(x+1, y+1, z  , Stencil::idx[SW]);
   const real_t dd_tmp_S  = src->get(x  , y+1, z  , Stencil::idx[S]);
   const real_t dd_tmp_SE = src->get(x-1, y+1, z  , Stencil::idx[SE]);
   const real_t dd_tmp_E  = src->get(x-1, y  , z  , Stencil::idx[E]);
   const real_t dd_tmp_T  = src->get(x  , y  , z-1, Stencil::idx[T]);
   const real_t dd_tmp_TE = src->get(x-1, y  , z-1, Stencil::idx[TE]);
   const real_t dd_tmp_TN = src->get(x  , y-1, z-1, Stencil::idx[TN]);
   const real_t dd_tmp_TW = src->get(x+1, y  , z-1, Stencil::idx[TW]);
   const real_t dd_tmp_TS = src->get(x  , y+1, z-1, Stencil::idx[TS]);
   const real_t dd_tmp_B  = src->get(x  , y  , z+1, Stencil::idx[B]);
   const real_t dd_tmp_BE = src->get(x-1, y  , z+1, Stencil::idx[BE]);
   const real_t dd_tmp_BN = src->get(x  , y-1, z+1, Stencil::idx[BN]);
   const real_t dd_tmp_BW = src->get(x+1, y  , z+1, Stencil::idx[BW]);
   const real_t dd_tmp_BS = src->get(x  , y+1, z+1, Stencil::idx[BS]);
   const real_t dd_tmp_C  = src->get(x  , y  , z  , Stencil::idx[C]);

   const real_t velX_trm = dd_tmp_E + dd_tmp_NE + dd_tmp_SE + dd_tmp_TE + dd_tmp_BE;
   const real_t velY_trm = dd_tmp_N + dd_tmp_NW + dd_tmp_TN + dd_tmp_BN;
   const real_t velZ_trm = dd_tmp_T + dd_tmp_TS + dd_tmp_TW;

   const real_t rho = dd_tmp_C + dd_tmp_S + dd_tmp_W + dd_tmp_B + dd_tmp_SW + dd_tmp_BS + dd_tmp_BW + velX_trm + velY_trm + velZ_trm;

   const real_t velX = velX_trm - dd_tmp_W  - dd_tmp_NW - dd_tmp_SW - dd_tmp_TW - dd_tmp_BW;
   const real_t velY = velY_trm + dd_tmp_NE - dd_tmp_S  - dd_tmp_SW - dd_tmp_SE - dd_tmp_TS - dd_tmp_BS;
   const real_t velZ = velZ_trm + dd_tmp_TN + dd_tmp_TE - dd_tmp_B  - dd_tmp_BN - dd_tmp_BS - dd_tmp_BW - dd_tmp_BE;

   const real_t velXX = velX * velX;
   const real_t velYY = velY * velY;
   const real_t velZZ = velZ * velZ;

   const real_t dir_indep_trm = ( 1_r / 3_r ) * rho - 0.5_r * ( velXX + velYY + velZZ );

   dst->get(x,y,z,Stencil::idx[C]) = omega_trm_ * dd_tmp_C + omega_w0_ * dir_indep_trm;

   const real_t vel_trm_E_W = dir_indep_trm + 1.5_r * velXX;
   const real_t vel_trm_N_S = dir_indep_trm + 1.5_r * velYY;
   const real_t vel_trm_T_B = dir_indep_trm + 1.5_r * velZZ;

   dst->get(x,y,z,Stencil::idx[E]) = omega_trm_ * dd_tmp_E + omega_w1_ * ( vel_trm_E_W + velX );
   dst->get(x,y,z,Stencil::idx[W]) = omega_trm_ * dd_tmp_W + omega_w1_ * ( vel_trm_E_W - velX );
   dst->get(x,y,z,Stencil::idx[N]) = omega_trm_ * dd_tmp_N + omega_w1_ * ( vel_trm_N_S + velY );
   dst->get(x,y,z,Stencil::idx[S]) = omega_trm_ * dd_tmp_S + omega_w1_ * ( vel_trm_N_S - velY );
   dst->get(x,y,z,Stencil::idx[T]) = omega_trm_ * dd_tmp_T + omega_w1_ * ( vel_trm_T_B + velZ );
   dst->get(x,y,z,Stencil::idx[B]) = omega_trm_ * dd_tmp_B + omega_w1_ * ( vel_trm_T_B - velZ );

   const real_t velXmY = velX - velY;
   const real_t vel_trm_NW_SE = dir_indep_trm + 1.5_r * velXmY * velXmY;

   dst->get(x,y,z,Stencil::idx[NW]) = omega_trm_ * dd_tmp_NW + omega_w2_ * ( vel_trm_NW_SE - velXmY );
   dst->get(x,y,z,Stencil::idx[SE]) = omega_trm_ * dd_tmp_SE + omega_w2_ * ( vel_trm_NW_SE + velXmY );

   const real_t velXpY = velX + velY;
   const real_t vel_trm_NE_SW = dir_indep_trm + 1.5_r * velXpY * velXpY;

   dst->get(x,y,z,Stencil::idx[NE]) = omega_trm_ * dd_tmp_NE + omega_w2_ * ( vel_trm_NE_SW + velXpY );
   dst->get(x,y,z,Stencil::idx[SW]) = omega_trm_ * dd_tmp_SW + omega_w2_ * ( vel_trm_NE_SW - velXpY );

   const real_t velXmZ = velX - velZ;
   const real_t vel_trm_TW_BE = dir_indep_trm + 1.5_r * velXmZ * velXmZ;

   dst->get(x,y,z,Stencil::idx[TW]) = omega_trm_ * dd_tmp_TW + omega_w2_ * ( vel_trm_TW_BE - velXmZ );
   dst->get(x,y,z,Stencil::idx[BE]) = omega_trm_ * dd_tmp_BE + omega_w2_ * ( vel_trm_TW_BE + velXmZ );

   const real_t velXpZ = velX + velZ;
   const real_t vel_trm_TE_BW = dir_indep_trm + 1.5_r * velXpZ * velXpZ;

   dst->get(x,y,z,Stencil::idx[TE]) = omega_trm_ * dd_tmp_TE + omega_w2_ * ( vel_trm_TE_BW + velXpZ );
   dst->get(x,y,z,Stencil::idx[BW]) = omega_trm_ * dd_tmp_BW + omega_w2_ * ( vel_trm_TE_BW - velXpZ );

   const real_t velYmZ = velY - velZ;
   const real_t vel_trm_TS_BN = dir_indep_trm + 1.5_r * velYmZ * velYmZ;

   dst->get(x,y,z,Stencil::idx[TS]) = omega_trm_ * dd_tmp_TS + omega_w2_ * ( vel_trm_TS_BN - velYmZ );
   dst->get(x,y,z,Stencil::idx[BN]) = omega_trm_ * dd_tmp_BN + omega_w2_ * ( vel_trm_TS_BN + velYmZ );

   const real_t velYpZ = velY + velZ;
   const real_t vel_trm_TN_BS = dir_indep_trm + 1.5_r * velYpZ * velYpZ;

   dst->get(x,y,z,Stencil::idx[TN]) = omega_trm_ * dd_tmp_TN + omega_w2_ * ( vel_trm_TN_BS + velYpZ );
   dst->get(x,y,z,Stencil::idx[BS]) = omega_trm_ * dd_tmp_BS + omega_w2_ * ( vel_trm_TN_BS - velYpZ );
}

template< typename LatticeModel_T >
template< typename FieldPtrOrIterator >
void DefaultCellOperation< LatticeModel_T, typename boost::enable_if< boost::mpl::and_< boost::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                                                        collision_model::SRT_tag >,
                                                                                        boost::mpl::bool_< LatticeModel_T::CollisionModel::constant >,
                                                                                        boost::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >,
                                                                                        boost::mpl::not_< boost::mpl::bool_< LatticeModel_T::compressible > >,
                                                                                        boost::is_same< typename LatticeModel_T::ForceModel::tag,
                                                                                                        force_model::None_tag > > >::type
   >::operator()( FieldPtrOrIterator & src, FieldPtrOrIterator & dst ) const
{
   using namespace stencil;

   const real_t dd_tmp_NE = src.neighbor(-1, -1,  0, Stencil::idx[NE]);
   const real_t dd_tmp_N  = src.neighbor( 0, -1,  0, Stencil::idx[N] );
   const real_t dd_tmp_NW = src.neighbor(+1, -1,  0, Stencil::idx[NW]);
   const real_t dd_tmp_W  = src.neighbor(+1,  0,  0, Stencil::idx[W] );
   const real_t dd_tmp_SW = src.neighbor(+1, +1,  0, Stencil::idx[SW]);
   const real_t dd_tmp_S  = src.neighbor( 0, +1,  0, Stencil::idx[S] );
   const real_t dd_tmp_SE = src.neighbor(-1, +1,  0, Stencil::idx[SE]);
   const real_t dd_tmp_E  = src.neighbor(-1,  0,  0, Stencil::idx[E] );
   const real_t dd_tmp_T  = src.neighbor( 0,  0, -1, Stencil::idx[T] );
   const real_t dd_tmp_TE = src.neighbor(-1,  0, -1, Stencil::idx[TE]);
   const real_t dd_tmp_TN = src.neighbor( 0, -1, -1, Stencil::idx[TN]);
   const real_t dd_tmp_TW = src.neighbor(+1,  0, -1, Stencil::idx[TW]);
   const real_t dd_tmp_TS = src.neighbor( 0, +1, -1, Stencil::idx[TS]);
   const real_t dd_tmp_B  = src.neighbor( 0,  0, +1, Stencil::idx[B] );
   const real_t dd_tmp_BE = src.neighbor(-1,  0, +1, Stencil::idx[BE]);
   const real_t dd_tmp_BN = src.neighbor( 0, -1, +1, Stencil::idx[BN]);
   const real_t dd_tmp_BW = src.neighbor(+1,  0, +1, Stencil::idx[BW]);
   const real_t dd_tmp_BS = src.neighbor( 0, +1, +1, Stencil::idx[BS]);
   const real_t dd_tmp_C  = src.neighbor( 0,  0,  0, Stencil::idx[C] );

   const real_t velX_trm = dd_tmp_E + dd_tmp_NE + dd_tmp_SE + dd_tmp_TE + dd_tmp_BE;
   const real_t velY_trm = dd_tmp_N + dd_tmp_NW + dd_tmp_TN + dd_tmp_BN;
   const real_t velZ_trm = dd_tmp_T + dd_tmp_TS + dd_tmp_TW;

   const real_t rho = dd_tmp_C + dd_tmp_S + dd_tmp_W + dd_tmp_B + dd_tmp_SW + dd_tmp_BS + dd_tmp_BW + velX_trm + velY_trm + velZ_trm;

   const real_t velX = velX_trm - dd_tmp_W  - dd_tmp_NW - dd_tmp_SW - dd_tmp_TW - dd_tmp_BW;
   const real_t velY = velY_trm + dd_tmp_NE - dd_tmp_S  - dd_tmp_SW - dd_tmp_SE - dd_tmp_TS - dd_tmp_BS;
   const real_t velZ = velZ_trm + dd_tmp_TN + dd_tmp_TE - dd_tmp_B  - dd_tmp_BN - dd_tmp_BS - dd_tmp_BW - dd_tmp_BE;

   const real_t velXX = velX * velX;
   const real_t velYY = velY * velY;
   const real_t velZZ = velZ * velZ;

   const real_t dir_indep_trm = ( 1_r / 3_r ) * rho - 0.5_r * ( velXX + velYY + velZZ );

   dst[ Stencil::idx[C] ] = omega_trm_ * dd_tmp_C + omega_w0_ * dir_indep_trm;

   const real_t vel_trm_E_W = dir_indep_trm + 1.5_r * velXX;
   const real_t vel_trm_N_S = dir_indep_trm + 1.5_r * velYY;
   const real_t vel_trm_T_B = dir_indep_trm + 1.5_r * velZZ;

   dst[ Stencil::idx[E] ] = omega_trm_ * dd_tmp_E + omega_w1_ * ( vel_trm_E_W + velX );
   dst[ Stencil::idx[W] ] = omega_trm_ * dd_tmp_W + omega_w1_ * ( vel_trm_E_W - velX );
   dst[ Stencil::idx[N] ] = omega_trm_ * dd_tmp_N + omega_w1_ * ( vel_trm_N_S + velY );
   dst[ Stencil::idx[S] ] = omega_trm_ * dd_tmp_S + omega_w1_ * ( vel_trm_N_S - velY );
   dst[ Stencil::idx[T] ] = omega_trm_ * dd_tmp_T + omega_w1_ * ( vel_trm_T_B + velZ );
   dst[ Stencil::idx[B] ] = omega_trm_ * dd_tmp_B + omega_w1_ * ( vel_trm_T_B - velZ );

   const real_t velXmY = velX - velY;
   const real_t vel_trm_NW_SE = dir_indep_trm + 1.5_r * velXmY * velXmY;

   dst[ Stencil::idx[NW] ] = omega_trm_ * dd_tmp_NW + omega_w2_ * ( vel_trm_NW_SE - velXmY );
   dst[ Stencil::idx[SE] ] = omega_trm_ * dd_tmp_SE + omega_w2_ * ( vel_trm_NW_SE + velXmY );

   const real_t velXpY = velX + velY;
   const real_t vel_trm_NE_SW = dir_indep_trm + 1.5_r * velXpY * velXpY;

   dst[ Stencil::idx[NE] ] = omega_trm_ * dd_tmp_NE + omega_w2_ * ( vel_trm_NE_SW + velXpY );
   dst[ Stencil::idx[SW] ] = omega_trm_ * dd_tmp_SW + omega_w2_ * ( vel_trm_NE_SW - velXpY );

   const real_t velXmZ = velX - velZ;
   const real_t vel_trm_TW_BE = dir_indep_trm + 1.5_r * velXmZ * velXmZ;

   dst[ Stencil::idx[TW] ] = omega_trm_ * dd_tmp_TW + omega_w2_ * ( vel_trm_TW_BE - velXmZ );
   dst[ Stencil::idx[BE] ] = omega_trm_ * dd_tmp_BE + omega_w2_ * ( vel_trm_TW_BE + velXmZ );

   const real_t velXpZ = velX + velZ;
   const real_t vel_trm_TE_BW = dir_indep_trm + 1.5_r * velXpZ * velXpZ;

   dst[ Stencil::idx[TE] ] = omega_trm_ * dd_tmp_TE + omega_w2_ * ( vel_trm_TE_BW + velXpZ );
   dst[ Stencil::idx[BW] ] = omega_trm_ * dd_tmp_BW + omega_w2_ * ( vel_trm_TE_BW - velXpZ );

   const real_t velYmZ = velY - velZ;
   const real_t vel_trm_TS_BN = dir_indep_trm + 1.5_r * velYmZ * velYmZ;

   dst[ Stencil::idx[TS] ] = omega_trm_ * dd_tmp_TS + omega_w2_ * ( vel_trm_TS_BN - velYmZ );
   dst[ Stencil::idx[BN] ] = omega_trm_ * dd_tmp_BN + omega_w2_ * ( vel_trm_TS_BN + velYmZ );

   const real_t velYpZ = velY + velZ;
   const real_t vel_trm_TN_BS = dir_indep_trm + 1.5_r * velYpZ * velYpZ;

   dst[ Stencil::idx[TN] ] = omega_trm_ * dd_tmp_TN + omega_w2_ * ( vel_trm_TN_BS + velYpZ );
   dst[ Stencil::idx[BS] ] = omega_trm_ * dd_tmp_BS + omega_w2_ * ( vel_trm_TN_BS - velYpZ );
}



///////////////////////////////
// Specialization for D3Q19: //
// - compressible            //
// - no additional forces    //
///////////////////////////////

template< typename LatticeModel_T >
class DefaultCellOperation< LatticeModel_T, typename boost::enable_if< boost::mpl::and_< boost::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                                                         collision_model::SRT_tag >,
                                                                                         boost::mpl::bool_< LatticeModel_T::CollisionModel::constant >,
                                                                                         boost::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >,
                                                                                         boost::mpl::bool_< LatticeModel_T::compressible >,
                                                                                         boost::is_same< typename LatticeModel_T::ForceModel::tag,
                                                                                                         force_model::None_tag > > >::type >
{
public:

   static_assert( (boost::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( (boost::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value),                       "Only works with D3Q19!" );
   static_assert( LatticeModel_T::compressible,                                                                      "Only works with compressible models!" );
   static_assert( (boost::is_same< typename LatticeModel_T::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   typedef PdfField< LatticeModel_T >        PdfField_T;
   typedef typename LatticeModel_T::Stencil  Stencil;

   DefaultCellOperation() : omega_trm_( 0_r ), omega_w0_( 0_r ), omega_w1_( 0_r ), omega_w2_( 0_r ) {}

   void configure( const LatticeModel_T & latticeModel )
   {
      const real_t omega = latticeModel.collisionModel().omega();
      omega_trm_ = 1_r - omega;
      omega_w0_  = 3_r * ( 1_r / real_t( 3) ) * omega;
      omega_w1_  = 3_r * ( 1_r / 18_r ) * omega;
      omega_w2_  = 3_r * ( 1_r / 36_r ) * omega;
   }

   void operator()( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const;

   template< typename FieldPtrOrIterator >
   void operator()( FieldPtrOrIterator & src, FieldPtrOrIterator & dst ) const;


private:

   real_t omega_trm_;
   real_t omega_w0_;
   real_t omega_w1_;
   real_t omega_w2_;
};

template< typename LatticeModel_T >
void DefaultCellOperation< LatticeModel_T, typename boost::enable_if< boost::mpl::and_< boost::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                                                        collision_model::SRT_tag >,
                                                                                        boost::mpl::bool_< LatticeModel_T::CollisionModel::constant >,
                                                                                        boost::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >,
                                                                                        boost::mpl::bool_< LatticeModel_T::compressible >,
                                                                                        boost::is_same< typename LatticeModel_T::ForceModel::tag,
                                                                                                        force_model::None_tag > > >::type
   >::operator()( PdfField_T * src, PdfField_T * dst, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const
{
   using namespace stencil;

   const real_t dd_tmp_NE = src->get(x-1, y-1, z  , Stencil::idx[NE]);
   const real_t dd_tmp_N  = src->get(x  , y-1, z  , Stencil::idx[N]);
   const real_t dd_tmp_NW = src->get(x+1, y-1, z  , Stencil::idx[NW]);
   const real_t dd_tmp_W  = src->get(x+1, y  , z  , Stencil::idx[W]);
   const real_t dd_tmp_SW = src->get(x+1, y+1, z  , Stencil::idx[SW]);
   const real_t dd_tmp_S  = src->get(x  , y+1, z  , Stencil::idx[S]);
   const real_t dd_tmp_SE = src->get(x-1, y+1, z  , Stencil::idx[SE]);
   const real_t dd_tmp_E  = src->get(x-1, y  , z  , Stencil::idx[E]);
   const real_t dd_tmp_T  = src->get(x  , y  , z-1, Stencil::idx[T]);
   const real_t dd_tmp_TE = src->get(x-1, y  , z-1, Stencil::idx[TE]);
   const real_t dd_tmp_TN = src->get(x  , y-1, z-1, Stencil::idx[TN]);
   const real_t dd_tmp_TW = src->get(x+1, y  , z-1, Stencil::idx[TW]);
   const real_t dd_tmp_TS = src->get(x  , y+1, z-1, Stencil::idx[TS]);
   const real_t dd_tmp_B  = src->get(x  , y  , z+1, Stencil::idx[B]);
   const real_t dd_tmp_BE = src->get(x-1, y  , z+1, Stencil::idx[BE]);
   const real_t dd_tmp_BN = src->get(x  , y-1, z+1, Stencil::idx[BN]);
   const real_t dd_tmp_BW = src->get(x+1, y  , z+1, Stencil::idx[BW]);
   const real_t dd_tmp_BS = src->get(x  , y+1, z+1, Stencil::idx[BS]);
   const real_t dd_tmp_C  = src->get(x  , y  , z  , Stencil::idx[C]);

   const real_t velX_trm = dd_tmp_E + dd_tmp_NE + dd_tmp_SE + dd_tmp_TE + dd_tmp_BE;
   const real_t velY_trm = dd_tmp_N + dd_tmp_NW + dd_tmp_TN + dd_tmp_BN;
   const real_t velZ_trm = dd_tmp_T + dd_tmp_TS + dd_tmp_TW;

   const real_t rho = dd_tmp_C + dd_tmp_S + dd_tmp_W + dd_tmp_B + dd_tmp_SW + dd_tmp_BS + dd_tmp_BW + velX_trm + velY_trm + velZ_trm;
   const real_t invRho = 1.0_r / rho;

   const real_t velX = invRho * ( velX_trm - dd_tmp_W  - dd_tmp_NW - dd_tmp_SW - dd_tmp_TW - dd_tmp_BW );
   const real_t velY = invRho * ( velY_trm + dd_tmp_NE - dd_tmp_S  - dd_tmp_SW - dd_tmp_SE - dd_tmp_TS - dd_tmp_BS );
   const real_t velZ = invRho * ( velZ_trm + dd_tmp_TN + dd_tmp_TE - dd_tmp_B  - dd_tmp_BN - dd_tmp_BS - dd_tmp_BW - dd_tmp_BE );

   const real_t velXX = velX * velX;
   const real_t velYY = velY * velY;
   const real_t velZZ = velZ * velZ;

   const real_t dir_indep_trm = ( 1_r / 3_r ) - 0.5_r * ( velXX + velYY + velZZ );

   dst->get(x,y,z,Stencil::idx[C]) = omega_trm_ * dd_tmp_C + omega_w0_ * rho * dir_indep_trm;

   const real_t omega_w1_rho = omega_w1_ * rho;

   const real_t vel_trm_E_W = dir_indep_trm + 1.5_r * velXX;
   const real_t vel_trm_N_S = dir_indep_trm + 1.5_r * velYY;
   const real_t vel_trm_T_B = dir_indep_trm + 1.5_r * velZZ;

   dst->get(x,y,z,Stencil::idx[E]) = omega_trm_ * dd_tmp_E + omega_w1_rho * ( vel_trm_E_W + velX );
   dst->get(x,y,z,Stencil::idx[W]) = omega_trm_ * dd_tmp_W + omega_w1_rho * ( vel_trm_E_W - velX );
   dst->get(x,y,z,Stencil::idx[N]) = omega_trm_ * dd_tmp_N + omega_w1_rho * ( vel_trm_N_S + velY );
   dst->get(x,y,z,Stencil::idx[S]) = omega_trm_ * dd_tmp_S + omega_w1_rho * ( vel_trm_N_S - velY );
   dst->get(x,y,z,Stencil::idx[T]) = omega_trm_ * dd_tmp_T + omega_w1_rho * ( vel_trm_T_B + velZ );
   dst->get(x,y,z,Stencil::idx[B]) = omega_trm_ * dd_tmp_B + omega_w1_rho * ( vel_trm_T_B - velZ );

   const real_t omega_w2_rho = omega_w2_ * rho;

   const real_t velXmY = velX - velY;
   const real_t vel_trm_NW_SE = dir_indep_trm + 1.5_r * velXmY * velXmY;

   dst->get(x,y,z,Stencil::idx[NW]) = omega_trm_ * dd_tmp_NW + omega_w2_rho * ( vel_trm_NW_SE - velXmY );
   dst->get(x,y,z,Stencil::idx[SE]) = omega_trm_ * dd_tmp_SE + omega_w2_rho * ( vel_trm_NW_SE + velXmY );

   const real_t velXpY = velX + velY;
   const real_t vel_trm_NE_SW = dir_indep_trm + 1.5_r * velXpY * velXpY;

   dst->get(x,y,z,Stencil::idx[NE]) = omega_trm_ * dd_tmp_NE + omega_w2_rho * ( vel_trm_NE_SW + velXpY );
   dst->get(x,y,z,Stencil::idx[SW]) = omega_trm_ * dd_tmp_SW + omega_w2_rho * ( vel_trm_NE_SW - velXpY );

   const real_t velXmZ = velX - velZ;
   const real_t vel_trm_TW_BE = dir_indep_trm + 1.5_r * velXmZ * velXmZ;

   dst->get(x,y,z,Stencil::idx[TW]) = omega_trm_ * dd_tmp_TW + omega_w2_rho * ( vel_trm_TW_BE - velXmZ );
   dst->get(x,y,z,Stencil::idx[BE]) = omega_trm_ * dd_tmp_BE + omega_w2_rho * ( vel_trm_TW_BE + velXmZ );

   const real_t velXpZ = velX + velZ;
   const real_t vel_trm_TE_BW = dir_indep_trm + 1.5_r * velXpZ * velXpZ;

   dst->get(x,y,z,Stencil::idx[TE]) = omega_trm_ * dd_tmp_TE + omega_w2_rho * ( vel_trm_TE_BW + velXpZ );
   dst->get(x,y,z,Stencil::idx[BW]) = omega_trm_ * dd_tmp_BW + omega_w2_rho * ( vel_trm_TE_BW - velXpZ );

   const real_t velYmZ = velY - velZ;
   const real_t vel_trm_TS_BN = dir_indep_trm + 1.5_r * velYmZ * velYmZ;

   dst->get(x,y,z,Stencil::idx[TS]) = omega_trm_ * dd_tmp_TS + omega_w2_rho * ( vel_trm_TS_BN - velYmZ );
   dst->get(x,y,z,Stencil::idx[BN]) = omega_trm_ * dd_tmp_BN + omega_w2_rho * ( vel_trm_TS_BN + velYmZ );

   const real_t velYpZ = velY + velZ;
   const real_t vel_trm_TN_BS = dir_indep_trm + 1.5_r * velYpZ * velYpZ;

   dst->get(x,y,z,Stencil::idx[TN]) = omega_trm_ * dd_tmp_TN + omega_w2_rho * ( vel_trm_TN_BS + velYpZ );
   dst->get(x,y,z,Stencil::idx[BS]) = omega_trm_ * dd_tmp_BS + omega_w2_rho * ( vel_trm_TN_BS - velYpZ );
}

template< typename LatticeModel_T >
template< typename FieldPtrOrIterator >
void DefaultCellOperation< LatticeModel_T, typename boost::enable_if< boost::mpl::and_< boost::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                                                        collision_model::SRT_tag >,
                                                                                        boost::mpl::bool_< LatticeModel_T::CollisionModel::constant >,
                                                                                        boost::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >,
                                                                                        boost::mpl::bool_< LatticeModel_T::compressible >,
                                                                                        boost::is_same< typename LatticeModel_T::ForceModel::tag,
                                                                                                        force_model::None_tag > > >::type
   >::operator()( FieldPtrOrIterator & src, FieldPtrOrIterator & dst ) const
{
   using namespace stencil;

   const real_t dd_tmp_NE = src.neighbor(-1, -1,  0, Stencil::idx[NE]);
   const real_t dd_tmp_N  = src.neighbor( 0, -1,  0, Stencil::idx[N] );
   const real_t dd_tmp_NW = src.neighbor(+1, -1,  0, Stencil::idx[NW]);
   const real_t dd_tmp_W  = src.neighbor(+1,  0,  0, Stencil::idx[W] );
   const real_t dd_tmp_SW = src.neighbor(+1, +1,  0, Stencil::idx[SW]);
   const real_t dd_tmp_S  = src.neighbor( 0, +1,  0, Stencil::idx[S] );
   const real_t dd_tmp_SE = src.neighbor(-1, +1,  0, Stencil::idx[SE]);
   const real_t dd_tmp_E  = src.neighbor(-1,  0,  0, Stencil::idx[E] );
   const real_t dd_tmp_T  = src.neighbor( 0,  0, -1, Stencil::idx[T] );
   const real_t dd_tmp_TE = src.neighbor(-1,  0, -1, Stencil::idx[TE]);
   const real_t dd_tmp_TN = src.neighbor( 0, -1, -1, Stencil::idx[TN]);
   const real_t dd_tmp_TW = src.neighbor(+1,  0, -1, Stencil::idx[TW]);
   const real_t dd_tmp_TS = src.neighbor( 0, +1, -1, Stencil::idx[TS]);
   const real_t dd_tmp_B  = src.neighbor( 0,  0, +1, Stencil::idx[B] );
   const real_t dd_tmp_BE = src.neighbor(-1,  0, +1, Stencil::idx[BE]);
   const real_t dd_tmp_BN = src.neighbor( 0, -1, +1, Stencil::idx[BN]);
   const real_t dd_tmp_BW = src.neighbor(+1,  0, +1, Stencil::idx[BW]);
   const real_t dd_tmp_BS = src.neighbor( 0, +1, +1, Stencil::idx[BS]);
   const real_t dd_tmp_C  = src.neighbor( 0,  0,  0, Stencil::idx[C] );

   const real_t velX_trm = dd_tmp_E + dd_tmp_NE + dd_tmp_SE + dd_tmp_TE + dd_tmp_BE;
   const real_t velY_trm = dd_tmp_N + dd_tmp_NW + dd_tmp_TN + dd_tmp_BN;
   const real_t velZ_trm = dd_tmp_T + dd_tmp_TS + dd_tmp_TW;

   const real_t rho = dd_tmp_C + dd_tmp_S + dd_tmp_W + dd_tmp_B + dd_tmp_SW + dd_tmp_BS + dd_tmp_BW + velX_trm + velY_trm + velZ_trm;
   const real_t invRho = 1.0_r / rho;

   const real_t velX = invRho * ( velX_trm - dd_tmp_W  - dd_tmp_NW - dd_tmp_SW - dd_tmp_TW - dd_tmp_BW );
   const real_t velY = invRho * ( velY_trm + dd_tmp_NE - dd_tmp_S  - dd_tmp_SW - dd_tmp_SE - dd_tmp_TS - dd_tmp_BS );
   const real_t velZ = invRho * ( velZ_trm + dd_tmp_TN + dd_tmp_TE - dd_tmp_B  - dd_tmp_BN - dd_tmp_BS - dd_tmp_BW - dd_tmp_BE );

   const real_t velXX = velX * velX;
   const real_t velYY = velY * velY;
   const real_t velZZ = velZ * velZ;

   const real_t dir_indep_trm = ( 1_r / 3_r ) - 0.5_r * ( velXX + velYY + velZZ );

   dst[ Stencil::idx[C] ] = omega_trm_ * dd_tmp_C + omega_w0_ * rho * dir_indep_trm;

   const real_t omega_w1_rho = omega_w1_ * rho;

   const real_t vel_trm_E_W = dir_indep_trm + 1.5_r * velXX;
   const real_t vel_trm_N_S = dir_indep_trm + 1.5_r * velYY;
   const real_t vel_trm_T_B = dir_indep_trm + 1.5_r * velZZ;

   dst[ Stencil::idx[E] ] = omega_trm_ * dd_tmp_E + omega_w1_rho * ( vel_trm_E_W + velX );
   dst[ Stencil::idx[W] ] = omega_trm_ * dd_tmp_W + omega_w1_rho * ( vel_trm_E_W - velX );
   dst[ Stencil::idx[N] ] = omega_trm_ * dd_tmp_N + omega_w1_rho * ( vel_trm_N_S + velY );
   dst[ Stencil::idx[S] ] = omega_trm_ * dd_tmp_S + omega_w1_rho * ( vel_trm_N_S - velY );
   dst[ Stencil::idx[T] ] = omega_trm_ * dd_tmp_T + omega_w1_rho * ( vel_trm_T_B + velZ );
   dst[ Stencil::idx[B] ] = omega_trm_ * dd_tmp_B + omega_w1_rho * ( vel_trm_T_B - velZ );

   const real_t omega_w2_rho = omega_w2_ * rho;

   const real_t velXmY = velX - velY;
   const real_t vel_trm_NW_SE = dir_indep_trm + 1.5_r * velXmY * velXmY;

   dst[ Stencil::idx[NW] ] = omega_trm_ * dd_tmp_NW + omega_w2_rho * ( vel_trm_NW_SE - velXmY );
   dst[ Stencil::idx[SE] ] = omega_trm_ * dd_tmp_SE + omega_w2_rho * ( vel_trm_NW_SE + velXmY );

   const real_t velXpY = velX + velY;
   const real_t vel_trm_NE_SW = dir_indep_trm + 1.5_r * velXpY * velXpY;

   dst[ Stencil::idx[NE] ] = omega_trm_ * dd_tmp_NE + omega_w2_rho * ( vel_trm_NE_SW + velXpY );
   dst[ Stencil::idx[SW] ] = omega_trm_ * dd_tmp_SW + omega_w2_rho * ( vel_trm_NE_SW - velXpY );

   const real_t velXmZ = velX - velZ;
   const real_t vel_trm_TW_BE = dir_indep_trm + 1.5_r * velXmZ * velXmZ;

   dst[ Stencil::idx[TW] ] = omega_trm_ * dd_tmp_TW + omega_w2_rho * ( vel_trm_TW_BE - velXmZ );
   dst[ Stencil::idx[BE] ] = omega_trm_ * dd_tmp_BE + omega_w2_rho * ( vel_trm_TW_BE + velXmZ );

   const real_t velXpZ = velX + velZ;
   const real_t vel_trm_TE_BW = dir_indep_trm + 1.5_r * velXpZ * velXpZ;

   dst[ Stencil::idx[TE] ] = omega_trm_ * dd_tmp_TE + omega_w2_rho * ( vel_trm_TE_BW + velXpZ );
   dst[ Stencil::idx[BW] ] = omega_trm_ * dd_tmp_BW + omega_w2_rho * ( vel_trm_TE_BW - velXpZ );

   const real_t velYmZ = velY - velZ;
   const real_t vel_trm_TS_BN = dir_indep_trm + 1.5_r * velYmZ * velYmZ;

   dst[ Stencil::idx[TS] ] = omega_trm_ * dd_tmp_TS + omega_w2_rho * ( vel_trm_TS_BN - velYmZ );
   dst[ Stencil::idx[BN] ] = omega_trm_ * dd_tmp_BN + omega_w2_rho * ( vel_trm_TS_BN + velYmZ );

   const real_t velYpZ = velY + velZ;
   const real_t vel_trm_TN_BS = dir_indep_trm + 1.5_r * velYpZ * velYpZ;

   dst[ Stencil::idx[TN] ] = omega_trm_ * dd_tmp_TN + omega_w2_rho * ( vel_trm_TN_BS + velYpZ );
   dst[ Stencil::idx[BS] ] = omega_trm_ * dd_tmp_BS + omega_w2_rho * ( vel_trm_TN_BS - velYpZ );
}












} // namespace lbm
} // namespace walberla
