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
//! \file BoundaryHandling.h
//! \ingroup lbm
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Default LBM Boundary Handling
//
//======================================================================================================================

#pragma once

#include "lbm/boundary/FreeSlip.h"
#include "lbm/boundary/NoSlip.h"
#include "lbm/boundary/SimplePressure.h"
#include "lbm/boundary/SimpleUBB.h"

#include "core/math/Vector3.h"
#include "domain_decomposition/BlockDataID.h"
#include "boundary/BoundaryHandling.h"
#include "lbm/field/PdfField.h"

namespace walberla {
namespace lbm {


//**********************************************************************************************************************
/*!
* \brief Creates a default boundary handling for LBM simulations
*
* \ingroup lbm
*
* This functor is usually used like
\code{.cpp}
typedef lbm::DefaultBoundaryHandlingFactory< LatticeModel_T, FlagField_T > Factory;
BlockDataID bid = Factory::addBoundaryHandlingToStorage( blocks, "boundary handling", flagFieldId, pdfFieldId,
                                                         vel0, vel1, pressure0, pressure1, setOfDomainFlags );
\endcode
*
* It adds no slip, free slip, velocity and pressure boundaries to the boundary handler.
* The registered FlagUIDs can be accessed via the get...() functions.
* The following boundary conditions and flags are available:
*
* Boundary Condition | FlagUID        | BoundaryUID
* ------------------ | -------------- | ---------------------------
* NoSlip             | getNoSlip()    | getNoSlipBoundaryUID()
* FreeSlip           | getFreeSlip()  | getFreeSlipBoundaryUID()
* SimpleUBB          | getVelocity0() | getVelocity0BoundaryUID()
* SimpleUBB          | getVelocity1() | getVelocity1BoundaryUID()
* SimplePressure     | getPressure0() | getPressure0BoundaryUID()
* SimplePressure     | getPressure1() | getPressure1BoundaryUID()
*
*
* \tparam LatticeModel  The lattice model used for the simulation
* \tparam FlagFieldT    Type of the used flag field
*/
//**********************************************************************************************************************
template <typename LatticeModel, typename FlagFieldT >
class DefaultBoundaryHandlingFactory
{
public:
   typedef typename FlagFieldT::flag_t            flag_t;
   typedef typename LatticeModel::Stencil         Stencil;
   typedef NoSlip< LatticeModel, flag_t >         BcNoSlip;
   typedef FreeSlip< LatticeModel, FlagFieldT >   BcFreeSlip;
   typedef SimpleUBB< LatticeModel, flag_t >      BcSimpleUBB;
   typedef SimplePressure< LatticeModel, flag_t > BcSimplePressure;
   typedef Vector3<real_t>                        Velocity;
   typedef PdfField< LatticeModel >               PdfFieldLM;

   typedef walberla::boundary::BoundaryHandling< FlagFieldT, Stencil, BcNoSlip, BcFreeSlip, BcSimpleUBB, BcSimpleUBB, BcSimplePressure, BcSimplePressure > BoundaryHandling;

   static BlockDataID addBoundaryHandlingToStorage( const shared_ptr< StructuredBlockStorage > & bs, const std::string & identifier,
                                                    BlockDataID flagFieldID, BlockDataID pdfFieldID, const Set< FlagUID > & flagUIDSet,
                                                    const Vector3<real_t> & velocity0,
                                                    const Vector3<real_t> & velocity1,
                                                    const real_t pressure0,
                                                    const real_t pressure1 )
   {
      DefaultBoundaryHandlingFactory factory ( flagFieldID, pdfFieldID, flagUIDSet, velocity0, velocity1, pressure0, pressure1 );

      return bs->addStructuredBlockData< BoundaryHandling >( factory, identifier );
   }

   static const walberla::FlagUID & getNoSlip()    { static walberla::FlagUID uid( "NoSlip" );    return uid; }
   static const walberla::FlagUID & getFreeSlip()  { static walberla::FlagUID uid( "FreeSlip" );  return uid; }
   static const walberla::FlagUID & getVelocity0() { static walberla::FlagUID uid( "Velocity0" ); return uid; }
   static const walberla::FlagUID & getVelocity1() { static walberla::FlagUID uid( "Velocity1" ); return uid; }
   static const walberla::FlagUID & getPressure0() { static walberla::FlagUID uid( "Pressure0" ); return uid; }
   static const walberla::FlagUID & getPressure1() { static walberla::FlagUID uid( "Pressure1" ); return uid; }

   static const walberla::BoundaryUID & getNoSlipBoundaryUID()    { static walberla::BoundaryUID uid( "NoSlip" );    return uid; }
   static const walberla::BoundaryUID & getFreeSlipBoundaryUID()  { static walberla::BoundaryUID uid( "FreeSlip" );  return uid; }
   static const walberla::BoundaryUID & getVelocity0BoundaryUID() { static walberla::BoundaryUID uid( "Velocity0" ); return uid; }
   static const walberla::BoundaryUID & getVelocity1BoundaryUID() { static walberla::BoundaryUID uid( "Velocity1" ); return uid; }
   static const walberla::BoundaryUID & getPressure0BoundaryUID() { static walberla::BoundaryUID uid( "Pressure0" ); return uid; }
   static const walberla::BoundaryUID & getPressure1BoundaryUID() { static walberla::BoundaryUID uid( "Pressure1" ); return uid; }


   DefaultBoundaryHandlingFactory( const BlockDataID & flagField, const BlockDataID & pdfField, const Set< FlagUID > & flagUIDSet,
                            const Velocity velocity0, const Velocity velocity1,
                            const real_t   pressure0, const real_t   pressure1 );

   BoundaryHandling * operator()( walberla::IBlock * const block, const walberla::StructuredBlockStorage * const storage ) const;

private:
   BlockDataID flagField_;
   BlockDataID pdfField_;

   const Set< FlagUID > flagUIDSet_;

   Velocity velocity0_, velocity1_;
   real_t   pressure0_, pressure1_;

}; // class DefaultBoundaryHandlingFactory


//**********************************************************************************************************************
/*!
* \ingroup lbm
*
* \param flagField  BlockDataID of the flag field used in the simulation
* \param pdfField   BlockDataID of the PDF field used in the simulation
* \param velocity0  Velocity parameter for SimpleUBB "Velocity0"
* \param velocity1  Velocity parameter for SimpleUBB "Velocity1"
* \param pressure0  Pressure parameter for SimplePressure "Pressure0"
* \param pressure1  Pressure parameter for SimplePressure "Pressure1"
*/
//**********************************************************************************************************************
template <typename LatticeModel, typename FlagFieldT >
DefaultBoundaryHandlingFactory<LatticeModel, FlagFieldT>::DefaultBoundaryHandlingFactory(
                                                   const BlockDataID & flagField, const BlockDataID & pdfField, const Set< FlagUID > & flagUIDSet,
                                                   const Velocity velocity0, const Velocity velocity1,
                                                   const real_t   pressure0, const real_t   pressure1 ) :
   flagField_( flagField ), pdfField_( pdfField ), flagUIDSet_(flagUIDSet), velocity0_( velocity0 ), velocity1_( velocity1 ),
   pressure0_( pressure0 ), pressure1_( pressure1 )
{
}

template <typename LatticeModel, typename FlagFieldT >
typename DefaultBoundaryHandlingFactory<LatticeModel, FlagFieldT>::BoundaryHandling *
DefaultBoundaryHandlingFactory<LatticeModel, FlagFieldT>::operator()( walberla::IBlock * const block,
                                                               const walberla::StructuredBlockStorage * const /*storage*/ ) const
{
   PdfFieldLM * const pdfField  = block->getData< PdfFieldLM >( pdfField_  );
   FlagFieldT * const flagField = block->getData< FlagFieldT >( flagField_ );

   flag_t mask = 0;
   for( auto flag = flagUIDSet_.begin(); flag != flagUIDSet_.end(); ++flag )
      mask = static_cast< flag_t >( mask | flagField->getOrRegisterFlag( *flag ) );

   BoundaryHandling * handling = new BoundaryHandling( "default lbm boundary handling", flagField, mask,
        BcNoSlip        ( getNoSlipBoundaryUID(),    getNoSlip(),    pdfField ),
        BcFreeSlip      ( getFreeSlipBoundaryUID(),  getFreeSlip(),  pdfField, flagField, mask ),
        BcSimpleUBB     ( getVelocity0BoundaryUID(), getVelocity0(), pdfField, velocity0_ ),
        BcSimpleUBB     ( getVelocity1BoundaryUID(), getVelocity1(), pdfField, velocity1_ ),
        BcSimplePressure( getPressure0BoundaryUID(), getPressure0(), pdfField, pressure0_ ),
        BcSimplePressure( getPressure1BoundaryUID(), getPressure1(), pdfField, pressure1_ )
    );

   return handling;
}


} // namespace lbm
} // namespace walberla
