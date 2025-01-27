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
//! \file Squirmer.h
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#pragma once

#include "Sphere.h"

namespace walberla {
namespace pe {

class Squirmer: public Sphere
{

public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit Squirmer( id_t sid, id_t uid, const Vec3& gpos, const Vec3& rpos, const Quat& q,
                    real_t radius, real_t squirmerVelocity, real_t squirmerBeta, MaterialID material,
                    const bool global, const bool communicating, const bool infiniteMass );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   virtual ~Squirmer();
   //@}
   //**********************************************************************************************
   //**********************************************************************************************

public:
   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   inline real_t getSquirmerVelocity() const;
   inline real_t getSquirmerBeta()     const;

   inline const Vec3     velFromBF       ( real_t px, real_t py, real_t pz ) const WALBERLA_OVERRIDE;
   inline const Vec3     velFromBF       ( const Vec3& rpos )                const WALBERLA_OVERRIDE;
   inline const Vec3     velFromWF       ( real_t px, real_t py, real_t pz ) const WALBERLA_OVERRIDE;
   inline const Vec3     velFromWF       ( const Vec3& gpos )                const WALBERLA_OVERRIDE;

   static inline id_t getStaticTypeID() WALBERLA_OVERRIDE;
   //@}
   //**********************************************************************************************

protected:

   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   const Vec3 getSquirmerVelocity( const Vec3& rpos ) const;
   //@}

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   real_t squirmerVelocity_;  //!< Velocity of the squirmer.
   real_t squirmerBeta_;  //!< Dipolar character of the squirmer.
   //@}
   //**********************************************************************************************

private:
   static id_t staticTypeID_;  //< type id of sphere, will be set by SetBodyTypeIDs
   static void setStaticTypeID(id_t typeID) WALBERLA_OVERRIDE {staticTypeID_ = typeID;}

   //** friend declaration
   /// needed to be able to set static type ids with setStaticTypeID
   template <class T, int N>
   friend struct SetBodyTypeIDs;
};

//=================================================================================================
//
//  GET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the squirmer's velocity.
 *
 * \return The velocity of the squirmer.
 */
inline real_t Squirmer::getSquirmerVelocity() const
{
   return squirmerVelocity_;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns the squirmer's dipole characteristic.
 *
 * \return The dipole characteristic of the squirmer.
 */
inline real_t Squirmer::getSquirmerBeta() const
{
   return squirmerBeta_;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns unique type id of this type
 *
 * \return geometry specific type id
 */
inline id_t Squirmer::getStaticTypeID()
{
   return staticTypeID_;
}

} // namespace pe
} // namespace walberla
