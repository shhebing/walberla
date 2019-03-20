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
//! \file Union.h
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <pe/Config.h>
#include <pe/rigidbody/BodyStorage.h>
#include <pe/rigidbody/RigidBody.h>
#include <pe/rigidbody/RigidBodyCastIterator.h>
#include <pe/rigidbody/RigidBodyIterator.h>
#include <pe/Types.h>

#include <core/debug/Debug.h>
#include <core/math/Matrix3.h>
#include <core/math/Vector3.h>
#include <core/DataTypes.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include <boost/tuple/tuple.hpp>

namespace walberla {
namespace pe {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/**
 * \ingroup pe
 * \brief Base class for the sphere geometry.
 *
 * The Sphere class represents the base class for the sphere geometry. It provides
 * the basic functionality of a sphere. For a full description of the sphere geometry,
 * see the Sphere class description.
 */
template <typename BodyTypeTuple>
class Union : public RigidBody
{
public:
   //**Type definitions****************************************************************************
   using BodyTypeTupleT        = BodyTypeTuple;

   //**********************************************************************************************

   using size_type             = BodyStorage::size_type;           //!< Size type of the body storage.
   using iterator              = BodyStorage::iterator;            //!< Iterator over non-const bodies.
   using const_iterator        = BodyStorage::const_iterator;      //!< Iterator over constant bodies.
   template <typename C> //cast type
   using cast_iterator         = BodyStorage::cast_iterator<C>;
   template <typename C> //cast type
   using const_cast_iterator   = BodyStorage::const_cast_iterator<C>;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit Union( id_t sid, id_t uid, const Vec3& gpos, const Vec3& rpos, const Quat& q,
                   const bool global, const bool communicating, const bool infiniteMass );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   virtual ~Union();
   //@}
   //**********************************************************************************************
   //**********************************************************************************************

public:
   //**Get functions*******************************************************************************
   /*!\name BodyStorage functions */
   //@{
   inline bool                 isEmpty() const {return bodies_.isEmpty();}
   inline size_t               size()    const {return bodies_.size();}

   inline iterator       begin   ()       {return bodies_.begin();}
   inline const_iterator begin   () const {return bodies_.begin();}
   inline const_iterator cbegin  () const {return bodies_.cbegin();}
   inline iterator       end     ()       {return bodies_.end();}
   inline const_iterator end     () const {return bodies_.end();}
   inline const_iterator cend    () const {return bodies_.cend();}

   template< typename C >
   inline cast_iterator<C> begin()              {return bodies_.begin<C>();}
   template< typename C >
   inline const_cast_iterator<C> begin() const  {return bodies_.begin<C>();}
   template< typename C >
   inline const_cast_iterator<C> cbegin() const {return bodies_.cbegin<C>();}
   template< typename C >
   inline cast_iterator<C> end()                {return bodies_.end<C>();}
   template< typename C >
   inline const_cast_iterator<C> end() const    {return bodies_.end<C>();}
   template< typename C >
   inline const_cast_iterator<C> cend() const   {return bodies_.cend<C>();}
   //@}
   //**********************************************************************************************

   virtual inline real_t getVolume()         const;

   //**Set functions*******************************************************************************
   /*!\name Set functions */
   //@{
   virtual void setRemote( bool remote ) WALBERLA_OVERRIDE;
   //@}
   //**********************************************************************************************

   virtual inline bool   hasSubBodies()      const WALBERLA_OVERRIDE { return true; }

   //**Simulation functions************************************************************************
   /*!\name Simulation functions */
   //@{
   virtual void update( const Vec3& dp ) WALBERLA_OVERRIDE;  // Translation update of a subordinate rigid body
   virtual void update( const Quat& dq ) WALBERLA_OVERRIDE;  // Rotation update of a subordinate rigid body
   //@}
   //**********************************************************************************************

   //**Signal functions***************************************************************************
   /*!\name Signal functions */
   //@{
   virtual void handleModification() WALBERLA_OVERRIDE;
   virtual void handleTranslation() WALBERLA_OVERRIDE;
   virtual void handleRotation() WALBERLA_OVERRIDE;
   //@}
   //**********************************************************************************************

   //**Rigid body manager functions****************************************************************
   /*!\name Rigid body manager functions */
   //@{
   inline RigidBody& add( std::unique_ptr<RigidBody>&& body );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   static inline id_t getStaticTypeID();
   //@}
   //**********************************************************************************************

   //**Output functions****************************************************************************
   /*!\name Output functions */
   //@{
   virtual void print( std::ostream& os, const char* tab ) const;
   //@}
   //**********************************************************************************************

protected:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   virtual void setPositionImpl       ( real_t px, real_t py, real_t pz )         WALBERLA_OVERRIDE;
   virtual void setOrientationImpl    ( real_t r, real_t i, real_t j, real_t k )  WALBERLA_OVERRIDE;
   virtual void translateImpl         ( real_t dx, real_t dy, real_t dz )         WALBERLA_OVERRIDE;
   virtual void rotateImpl            ( const Quat& dq )                          WALBERLA_OVERRIDE;
   virtual void rotateAroundOriginImpl( const Quat& dq )                          WALBERLA_OVERRIDE;
   virtual void rotateAroundPointImpl ( const Vec3& point, const Quat& dq )       WALBERLA_OVERRIDE;
   virtual bool containsRelPointImpl   ( real_t px, real_t py, real_t pz ) const  WALBERLA_OVERRIDE;
   virtual bool isSurfaceRelPointImpl  ( real_t px, real_t py, real_t pz ) const  WALBERLA_OVERRIDE;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline virtual void calcBoundingBox() WALBERLA_OVERRIDE;  // Calculation of the axis-aligned bounding box
   void calcCenterOfMass(); ///< updates the center of mass (gpos)
   inline         void calcInertia();      // Calculation of the moment of inertia
   //@}
   //**********************************************************************************************

private:
   BodyStorage bodies_;  //!< Rigid bodies contained in the union.
   std::vector<id_t> containedTypeIDs_;

   static id_t staticTypeID_;  //< type id of sphere, will be set by SetBodyTypeIDs
   static void setStaticTypeID(id_t typeID) {staticTypeID_ = typeID;}

   //** friend declaration
   /// needed to be able to set static type ids with setStaticTypeID
   template <class T>
   friend struct SetBodyTypeIDs;
};
//*************************************************************************************************

//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================


//*************************************************************************************************
//*************************************************************************************************
/*!\brief Constructor for the Sphere class.
 *
 * \param sid Unique system-specific ID for the sphere.
 * \param uid User-specific ID for the sphere.
 * \param gpos Global geometric center of the sphere.
 * \param q The orientation of the sphere's body frame in the global world frame.
 * \param global specifies if the sphere should be created in the global storage
 * \param communicating specifies if the sphere should take part in synchronization (syncNextNeighbour, syncShadowOwner)
 * \param infiniteMass specifies if the sphere has infinite mass and will be treated as an obstacle
 */
template <typename BodyTypeTuple>
Union<BodyTypeTuple>::Union( id_t sid, id_t uid, const Vec3& gpos, const Vec3& rpos, const Quat& q,
                             const bool global, const bool communicating, const bool /*infiniteMass*/ )
   : RigidBody( getStaticTypeID(), sid, uid )  // Initialization of the parent class
{
   // Initializing the instantiated union
   gpos_   = gpos;                   // Setting the global center of mass
   rpos_   = rpos;                   // Setting the relative position
   q_      = q;                      // Setting the orientation
   R_      = q_.toRotationMatrix();  // Setting the rotation matrix

   calcCenterOfMass();
   calcInertia();

   setGlobal( global );
   //   setMass( infiniteMass );
   setCommunicating( communicating );
   setFinite( true );
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the Sphere class.
 */
template <typename BodyTypeTuple>
Union<BodyTypeTuple>::~Union()
{
   // Clearing the bodies
   bodies_.clear();

   // Logging the destruction of the sphere
   WALBERLA_LOG_DETAIL( "Destroyed union " << sid_ );
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns the volume of the union.
 *
 * \return The volume of the union.
 */
template <typename BodyTypeTuple>
inline real_t Union<BodyTypeTuple>::getVolume() const
{
   real_t volume(0);
   for( auto bodyIt=bodies_.begin(); bodyIt!=bodies_.end(); ++bodyIt )
      volume += bodyIt->getVolume();
   return volume;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Setting the remote flag of the union and all contained bodies.
 *
 * \param remote \a true to declare the union remote, \a false declare it local.
 * \return void
 *
 * This function sets the remote flag of the union and all contained rigid bodies. Note that
 * this function should not be used explicitly, but is automatically called during the MPI
 * communication to set the remote status of a union within the simulation world. Using
 * this function explicitly may lead to simulation errors during a parallel simulation!
 */
template <typename BodyTypeTuple>
void Union<BodyTypeTuple>::setRemote( bool remote )
{
   remote_ = remote;

   // Setting the remote flag of the contained rigid bodies
   for( auto bIt=bodies_.begin(); bIt!=bodies_.end(); ++bIt )
      bIt->setRemote( remote );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the bounding box of the union.
 *
 * \return void
 *
 * This function updates the axis-aligned bounding box of the union according to the current
 * position and orientation of the contained rigid bodies.
 */
template <typename BodyTypeTuple>
void Union<BodyTypeTuple>::calcBoundingBox()
{
   // Setting the bounding box of an empty union
   if( bodies_.isEmpty() ) {
      aabb_ = math::AABB(
                 gpos_[0] - 0.01_r,
            gpos_[1] - 0.01_r,
            gpos_[2] - 0.01_r,
            gpos_[0] + 0.01_r,
            gpos_[1] + 0.01_r,
            gpos_[2] + 0.01_r
            );
   }

   // Using the bounding box of the first contained bodies as initial bounding box
   // and merging it with the bounding boxes of all other bodies
   else {
      aabb_ = bodies_.begin()->getAABB();
      for( auto b=bodies_.begin()+1; b!=bodies_.end(); ++b )
         aabb_.merge( b->getAABB() );
   }

   WALBERLA_ASSERT( aabb_.checkInvariant() , "Invalid bounding box detected" );
   WALBERLA_ASSERT( aabb_.contains( gpos_ ), "Invalid bounding box detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the total mass and center of mass.
 *
 * \return void
 */
template <typename BodyTypeTuple>
void Union<BodyTypeTuple>::calcCenterOfMass()
{
   // Checking the state of the union
   WALBERLA_ASSERT( checkInvariants(), "Invalid union state detected" );

   // Initializing the total mass and the inverse mass
   mass_    = 0_r;
   invMass_ = 0_r;

   // Don't calculate the center of mass of an empty union
   if( bodies_.isEmpty() ) return;

   // Calculating the center of mass of a single body
   if( bodies_.size() == 1 )
   {
      const BodyID body( bodies_.begin().getBodyID() );
      gpos_ = body->getPosition();
      mass_ = body->getMass();
      if( !isFixed() && mass_ > 0_r )
         invMass_ = 1_r / mass_;
   }

   // Calculating the center of mass of a union containing several bodies
   else
   {
      // Resetting the center of mass
      gpos_.reset();

      // Calculating the center of mass of a finite union
      if( finite_ )
      {
         // Accumulating the mass of all contained rigid bodies
         for( auto b=bodies_.begin(); b!=bodies_.end(); ++b ) {
            WALBERLA_ASSERT( b->isFinite(), "Invalid infinite body in finite union detected" );
            mass_ += b->getMass();
            gpos_ += b->getPosition() * b->getMass();
         }

         // Calculating the center of mass for unions with non-zero mass
         if( mass_ > 0_r ) {
            if( !isFixed() ) invMass_ = 1_r / mass_;
            gpos_ /= mass_;
         }

         // Calculating the center of mass for unions with a mass of zero
         else {
            size_t counter( 0 );
            gpos_.reset();

            for( auto b=bodies_.begin(); b!=bodies_.end(); ++b ) {
               gpos_ += b->getPosition();
               ++counter;
            }

            gpos_ /= counter;
         }
      }

      // Calculating the center of mass of an infinite union
      else {
         size_t counter( 0 );

         for( auto b=bodies_.begin(); b!=bodies_.end(); ++b ) {
            if( b->isFinite() ) continue;
            gpos_ += b->getPosition();
            ++counter;
         }

         gpos_ /= counter;
      }
   }

   // Checking the state of the union
   WALBERLA_ASSERT( checkInvariants(), "Invalid union state detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the moment of inertia in reference to the body frame of the union.
 *
 * \return void
 */
template <typename BodyTypeTuple>
void Union<BodyTypeTuple>::calcInertia()
{
   // Checking the state of the union
   WALBERLA_ASSERT( checkInvariants(), "Invalid union state detected" );

   // Initializing the body moment of inertia and the inverse moment of inertia
   I_    = 0_r;
   Iinv_ = 0_r;

   // Don't calculate the moment of inertia of an infinite or empty union
   if( !isFinite() || bodies_.isEmpty() || floatIsEqual(mass_, 0_r) ) return;

   // Calculating the global moment of inertia
   real_t mass;
   Vec3   pos;

   for( auto b=bodies_.begin(); b!=bodies_.end(); ++b )
   {
      mass = b->getMass();
      pos  = b->getPosition() - gpos_;

      I_ += b->getInertia();
      I_[0] += mass * ( pos[1]*pos[1] + pos[2]*pos[2] );
      I_[1] -= mass * pos[0] * pos[1];
      I_[2] -= mass * pos[0] * pos[2];
      I_[3] -= mass * pos[0] * pos[1];
      I_[4] += mass * ( pos[0]*pos[0] + pos[2]*pos[2] );
      I_[5] -= mass * pos[1] * pos[2];
      I_[6] -= mass * pos[0] * pos[2];
      I_[7] -= mass * pos[1] * pos[2];
      I_[8] += mass * ( pos[0]*pos[0] + pos[1]*pos[1] );
   }

   // Rotating the moment of inertia from the global frame of reference to the body frame of reference
   I_ = R_ * I_ * R_;

   // Calculating the inverse of the body moment of inertia
   if( !isFixed() ) Iinv_ = I_.getInverse();

   // Checking the state of the union
   WALBERLA_ASSERT( checkInvariants(), "Invalid union state detected" );
}
//*************************************************************************************************




//=================================================================================================
//
//  SET FUNCTIONS
//
//=================================================================================================


//*************************************************************************************************
/*!\brief Setting the global position of the union.
 *
 * \param gpos The global position.
 * \return void
 * \exception std::logic_error Invalid translation of a global union inside an exclusive section.
 *
 * Setting the global position of the union's center of mass. If the union contains an infinite
 * rigid body, the function shifts the union to reposition its anchor point according to the
 * given global coordinate.
 */
template <typename BodyTypeTuple>
void Union<BodyTypeTuple>::setPositionImpl( real_t px, real_t py, real_t pz )
{
   Vec3 gpos = Vec3( px, py, pz );

   // Calculating the position change
   const Vec3 dp( gpos - gpos_ );

   // Setting the global position
   gpos_ = gpos;

   // Updating the contained bodies
   for( auto& b : bodies_ )
      b.update( dp );

   Union<BodyTypeTuple>::calcBoundingBox();    // Setting the axis-aligned bounding box
   wake();               // Waking the union from sleep mode
   signalTranslation();  // Signaling the position change to the superordinate body
}
//*************************************************************************************************



//*************************************************************************************************
/*!\brief Setting the global orientation of the union.
 *
 * \param r The value for the real part.
 * \param i The value for the first imaginary part.
 * \param j The value for the second imaginary part.
 * \param k The value for the third imaginary part.
 * \return void
 * \exception std::logic_error Invalid rotation of a global union inside an exclusive section.
 *
 * Setting the orientation/rotation of the entire union. This function causes all contained
 * primitives to rotate around the center of mass of the union (if the union is finite) or around
 * the anchor point of the union (if the union is infinite). The orientation of the rigid bodies
 * within the union in reference to the body frame of the union is not changed.
 */
template <typename BodyTypeTuple>
void Union<BodyTypeTuple>::setOrientationImpl( real_t r, real_t i, real_t j, real_t k )
{
   const Quat q ( r, i, j, k );
   const Quat dq( q * q_.getInverse() );

   q_ = q;
   R_ = q_.toRotationMatrix();

   // Updating the contained bodies
   for( auto& b : bodies_ )
      b.update( dq );

   Union<BodyTypeTuple>::calcBoundingBox();  // Setting the axis-aligned bounding box
   wake();             // Waking the union from sleep mode
   signalRotation();   // Signaling the change of orientation to the superordinate body
}
//*************************************************************************************************



//=================================================================================================
//
//  SIMULATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Translation update of a subordinate union.
 *
 * \param dp Change in the global position of the superordinate union.
 * \return void
 *
 * This update function is triggered by the superordinate body in case of a translational
 * movement. This movement involves a change in the global position and the axis-aligned
 * bounding box.
 */
template <typename BodyTypeTuple>
void Union<BodyTypeTuple>::update( const Vec3& dp )
{
   // Checking the state of the union
   WALBERLA_ASSERT( checkInvariants(), "Invalid union state detected" );
   WALBERLA_ASSERT( hasSuperBody(), "Invalid superordinate body detected" );

   // Updating the global position
   gpos_ += dp;

   // Updating the contained bodies
   for( auto& b : bodies_ )
      b.update( dp );

   // Setting the axis-aligned bounding box
   Union<BodyTypeTuple>::calcBoundingBox();

   // Checking the state of the union
   WALBERLA_ASSERT( checkInvariants(), "Invalid union state detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Rotation update of a subordinate union.
 *
 * \param dq Change in the orientation of the superordinate union.
 * \return void
 *
 * This update function is triggered by the superordinate body in case of a rotational movement.
 * This movement involves a change in the global position, the orientation/rotation and the
 * axis-aligned bounding box of the box.
 */
template <typename BodyTypeTuple>
void Union<BodyTypeTuple>::update( const Quat& dq )
{
   // Checking the state of the union
   WALBERLA_ASSERT( checkInvariants(), "Invalid union state detected" );
   WALBERLA_ASSERT( hasSuperBody(), "Invalid superordinate body detected" );

   // Calculating the new global position
   gpos_ = sb_->getPosition() + ( sb_->getRotation() * rpos_ );

   // Calculating the new orientation and rotation
   q_ = dq * q_;
   R_ = q_.toRotationMatrix();

   // Updating the contained bodies
   for( auto& b : bodies_ )
      b.update( dq );

   // Setting the axis-aligned bounding box
   Union<BodyTypeTuple>::calcBoundingBox();

   // Checking the state of the union
   WALBERLA_ASSERT( checkInvariants(), "Invalid union state detected" );
}
//*************************************************************************************************




//=================================================================================================
//
//  RIGID BODY MANAGER FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Adding a rigid body to the union.
 *
 * \param body The rigid body to be added to the union.
 * \return void
 * \exception std::logic_error Invalid adding into a global union inside an exclusive section.
 * \exception std::logic_error Global flags of body and union do not match.
 *
 * This function adds another rigid body to the union (as for instance a sphere, a box, a
 * capsule, or another union). It updates all union properties that change due to the new
 * rigid body: the center of mass, the relative position of all contained rigid bodies, the
 * attached sections of all contained links, the moment of inertia, and the axis-aligned
 * bounding box.\n
 * The union takes full responsibility for the newly added body, including the necessary
 * memory management. After adding the body to the union, the body is considered part of the
 * union. All functions called on the union (as for instance all kinds of set functions,
 * translation or rotation functions) additionally affect the rigid body. However, the body
 * can still individually be positioned, oriented, translated, rotated, or made (in-)visible
 * within the union.\n\n
 *
 *
 * \section union_add_infinite Adding infinite rigid bodies
 *
 * Adding an infinite rigid body (as for instance a plane) to the union makes the union an
 * infinite rigid body. This additionally resets the linear and angular velocity of the union
 * and fixes its global position. Note that removing the last infinite body from an union will
 * not restore previous settings such as the velocities and will not unfix the union!\n\n
 *
 *
 * \section union_add_global Global bodies/unions
 *
 * Adding a global rigid body (i.e. a body defined inside the pe_GLOBAL_SECTION) to the union
 * requires the union to be also global. Adding a non-global rigid body to a union requires
 * the union to be also non-global. The attempt to add a global rigid body to a non-global
 * union or a non-global body to a global union results in a \a std::logic_error exception.
 *
 *
 * \section union_add_rules Additional rules
 *
 * The following rules apply for the mobility and visibility of the resulting union:
 *  - If either the union or the added rigid body is fixed, the new compound will also be fixed.
 *    For instance, adding a fixed rigid body to an union will fix the union, and adding a rigid
 *    body to a fixed union will fix the body.
 *  - Neither the (in-)visibility of the added rigid body nor (in-)visibility of the union will
 *    change due to the add operation. For instance adding a visible rigid body to an invisible
 *    union will not change the visibility of the body. Neither is the visibility of the union
 *    changed. In order to change the visiblity, the setVisible() function can be either called
 *    individually for the rigid body (to exclusively make the body (in-)visible) or the entire
 *    union (to make the entire union (in-)visible.
 */
template <typename BodyTypeTuple>
RigidBody& Union<BodyTypeTuple>::add( std::unique_ptr<RigidBody>&& body )
{
   // Checking for "self-assignment"
   if( body.get() == BodyID( this ) ) return *this;

   // Checking the global flags of the body and the union in MPI parallel simulations
   if( body->isGlobal() ^ global_ )
      throw std::logic_error( "Global flags of body and union do not match" );

   // Registering the rigid body
   auto& bd = bodies_.add( std::move(body) );

   Vec3 oldCenterOfMass = getPosition();
   Vec3 oldImpulse      = getLinearVel() * getMass();

   Vec3 bodyCenterOfMass = bd.getPosition();
   Vec3 bodyImpulse      = bd.getLinearVel() * bd.getMass();

   bd.setSB(this); //having a superbody will forward all getVel calls to superbody!!!

   // Updating the axis-aligned bounding box
   if( bodies_.size() == 1 )
      aabb_ = bd.getAABB();
   else
      aabb_.merge( bd.getAABB() );

   // Setting the union's total mass and center of mass
   calcCenterOfMass();

   // Setting the contained primitives' relative position in reference to the center of mass
   for( auto& b : bodies_ )
      b.calcRelPosition();

   // Setting the moment of inertia
   calcInertia();

   //moving impuls to union
   setLinearVel( Vec3(0,0,0) );
   setAngularVel( Vec3(0,0,0) );
   addImpulseAtPos( oldImpulse, oldCenterOfMass);
   addImpulseAtPos( bodyImpulse, bodyCenterOfMass);
   bd.setLinearVel( Vec3(0,0,0) );
   bd.setAngularVel( Vec3(0,0,0) );

   // Signaling the internal modification to the superordinate body
   signalModification();

   return bodies_.back();
}
//*************************************************************************************************




//=================================================================================================
//
//  TRANSLATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Translation of the center of mass of the union by the displacement vector \a dp.
 *
 * \param dp The displacement vector.
 * \return void
 * \exception std::logic_error Invalid translation of a global union inside an exclusive section.
 *
 * Changing the global position of the entire union's center of mass. All contained rigid bodies
 * are moved by the same displacement.
 */
template <typename BodyTypeTuple>
void Union<BodyTypeTuple>::translateImpl( real_t dx, real_t dy, real_t dz )
{
   Vec3 dp(dx, dy, dz);

   // Changing the global position/reference point
   gpos_ += dp;

   // Updating the contained bodies
   for( auto& b : bodies_ )
      b.update( dp );

   Union<BodyTypeTuple>::calcBoundingBox();    // Setting the axis-aligned bounding box
   wake();               // Waking the union from sleep mode
   signalTranslation();  // Signaling the position change to the superordinate body
}
//*************************************************************************************************




//=================================================================================================
//
//  ROTATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Rotation of the union by the quaternion \a dq.
 *
 * \param dq The quaternion for the rotation.
 * \return void
 * \exception std::logic_error Invalid rotation of a global union inside an exclusive section.
 *
 * Changing the orientation/rotation of the entire union. This function causes all contained
 * rigid bodies to rotate around the center of mass of the union (if the union is finite) or
 * around the anchor point of the union (if the union is infinite). The orientation of the
 * bodies within the union in reference to the body frame of the union is not changed.
 */
template <typename BodyTypeTuple>
void Union<BodyTypeTuple>::rotateImpl( const Quat& dq )
{
   if (isFixed())
      WALBERLA_ABORT("Trying to rotate a fixed body: " << *this);

   q_ = dq * q_;                // Updating the orientation of the union
   R_ = q_.toRotationMatrix();  // Updating the rotation of the union

   // Updating the contained bodies
   for( auto& b : bodies_ )
      b.update( dq );

   Union<BodyTypeTuple>::calcBoundingBox();  // Setting the axis-aligned bounding box
   wake();             // Waking the union from sleep mode
   signalRotation();   // Signaling the change of orientation to the superordinate body
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Rotation of the union around the origin of the global world frame.
 *
 * \param xangle Rotation around the x-axis (radian measure).
 * \param yangle Rotation around the y-axis (radian measure).
 * \param zangle Rotation around the z-axis (radian measure).
 * \return void
 * \exception std::logic_error Invalid rotation of a global union inside an exclusive section.
 *
 * This function rotates the entire union around the origin of the global world frame and
 * changes both the global position and the orientation/rotation of the union. Additionally,
 * all contained rigid bodies change their position and orientation accordingly. The rotations
 * are applied in the order x, y, z. The orientation of the bodies within the union in
 * reference to the body frame of the union is not changed.
 */
template <typename BodyTypeTuple>
void Union<BodyTypeTuple>::rotateAroundOriginImpl( const Quat& dq )
{
   q_    = dq * q_;                // Updating the orientation of the union
   R_    = q_.toRotationMatrix();  // Updating the rotation of the union
   gpos_ = dq.rotate( gpos_ );     // Updating the global position of the union

   // Updating the contained bodies
   for( auto& b : bodies_ )
      b.update( dq );

   Union<BodyTypeTuple>::calcBoundingBox();    // Setting the axis-aligned bounding box
   wake();               // Waking the union from sleep mode
   signalTranslation();  // Signaling the position change to the superordinate body
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Rotation of the union around a specific global coordinate.
 *
 * \param point The global center of the rotation.
 * \param axis The global rotation axis.
 * \param angle The rotation angle (radian measure).
 * \return void
 * \exception std::logic_error Invalid rotation of a global union inside an exclusive section.
 *
 * This function rotates the entire union around the given global coordiante \a point and
 * changes both the global position and the orientation/rotation of the union. Additionally,
 * all contained rigid bodies change their position and orientation accordingly. The orientation
 * of the bodies within the union in reference to the body frame of the union is not changed.
 */
template <typename BodyTypeTuple>
void Union<BodyTypeTuple>::rotateAroundPointImpl( const Vec3& point, const Quat &dq )
{
   const Vec3 dp( gpos_ - point );

   q_    = dq * q_;                  // Updating the orientation of the union
   R_    = q_.toRotationMatrix();    // Updating the rotation of the union
   gpos_ = point + dq.rotate( dp );  // Updating the global position of the union

   // Updating the contained bodies
   for( auto& b : bodies_ )
      b.update( dq );

   Union<BodyTypeTuple>::calcBoundingBox();    // Setting the axis-aligned bounding box
   wake();               // Waking the union from sleep mode
   signalTranslation();  // Signaling the position change to the superordinate body
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies inside the union.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies inside the sphere, \a false if not.
 */
template <typename BodyTypeTuple>
bool Union<BodyTypeTuple>::containsRelPointImpl( real_t px, real_t py, real_t pz ) const
{
   const Vec3 gpos( pointFromBFtoWF( px, py, pz ) );
   for( auto& b : bodies_ )
      if( b.containsPoint( gpos ) ) return true;
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies on the surface of the union.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies on the surface of the sphere, \a false if not.
 *
 * The tolerance level of the check is pe::surfaceThreshold.
 */
template <typename BodyTypeTuple>
bool Union<BodyTypeTuple>::isSurfaceRelPointImpl( real_t px, real_t py, real_t pz ) const
{
   bool surface( false );
   const Vec3 gpos( pointFromBFtoWF( px, py, pz ) );

   for( auto& b : bodies_ )
   {
      if( b.containsPoint( gpos ) ) return false;
      else if( b.isSurfacePoint( gpos ) ) surface = true;
   }

   return surface;
}
//*************************************************************************************************




//=================================================================================================
//
//  FUNCTIONS FOR INTERNAL CHANGES IN COMPOUND GEOMETRIES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Signals an internal modification of a contained subordiante rigid body.
 *
 * \return void
 *
 * In case one of the contained rigid bodies is interally modified, this function is called to
 * recalculate the changed properties of the union.
 */
template <typename BodyTypeTuple>
void Union<BodyTypeTuple>::handleModification()
{
   // Setting the finite flag
   finite_ = true;
   for( auto b=bodies_.begin(); b!=bodies_.end(); ++b ) {
      if( !b->isFinite() ) {
         finite_ = false;
         break;
      }
   }

   // Setting the finiteness of the union
   if( !finite_ ) setFinite( false );

   // Setting the union's total mass and center of mass
   calcCenterOfMass();

   // Setting the contained primitives' relative position in reference to the center of mass
   for( auto& b : bodies_ )
      b.calcRelPosition();

   // Setting the moment of inertia
   calcInertia();

   // Updating the axis-aligned bounding box
   Union<BodyTypeTuple>::calcBoundingBox();

   // Signaling the internal modification to the superordinate body
   signalModification();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Signals a position change of a contained subordiante rigid body.
 *
 * \return void
 *
 * In case one of the contained rigid bodies changes its position, this function is called
 * to recalculate the changed properties of the union.
 */
template <typename BodyTypeTuple>
void Union<BodyTypeTuple>::handleTranslation()
{
   // Setting the union's total mass and center of mass
   calcCenterOfMass();

   // Setting the contained bodies' relative position in reference to the center of mass
   for( auto& b : bodies_ )
      b.calcRelPosition();

   // Setting the moment of inertia
   calcInertia();

   Union<BodyTypeTuple>::calcBoundingBox();    // Setting the axis-aligned bounding box
   wake();               // Waking the union from sleep mode
   signalTranslation();  // Signaling the position change to the superordinate body
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Signals an orientation change of a contained subordiante rigid body.
 *
 * \return void
 *
 * In case one of the contained rigid bodies changes its orientation, this function is called
 * to recalculate the changed properties of the union.
 */
template <typename BodyTypeTuple>
void Union<BodyTypeTuple>::handleRotation()
{
   // Setting the moment of inertia
   calcInertia();

   Union<BodyTypeTuple>::calcBoundingBox();  // Setting the axis-aligned bounding box
   wake();             // Waking the union from sleep mode
   signalRotation();   // Signaling the change of orientation to the superordinate body
}
//*************************************************************************************************




//=================================================================================================
//
//  OUTPUT FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Output of the current state of an union.
 *
 * \param os Reference to the output stream.
 * \param tab Indentation in front of every line of the union output.
 * \return void
 */
template <typename BodyTypeTuple>
void Union<BodyTypeTuple>::print( std::ostream& os, const char* tab ) const
{
   using std::setw;

   os << tab << " Union " << uid_ << " with " << bodies_.size();
   if( bodies_.size() == 1 ) os << " rigid body\n";
   else os << " rigid bodies\n";

   os << tab << "   Fixed: " << isFixed() << " , sleeping: " << !isAwake() << "\n";

   os << tab << "   System ID         = " << sid_ << "\n"
      << tab << "   Total mass        = ";
   if( isFinite() ) os << mass_ << "\n";
   else os << "*infinite*\n";

   os << tab << "   Global position   = " << gpos_ << "\n"
      << tab << "   Relative position = " << rpos_ << "\n"
      << tab << "   Linear velocity   = " << v_ << "\n"
      << tab << "   Angular velocity  = " << w_ << "\n";

   os << tab << "   Bounding box      = " << aabb_ << "\n"
      << tab << "   Quaternion        = " << q_ << "\n"
      << tab << "   Rotation matrix   = ( " << setw(9) << R_[0] << " , " << setw(9) << R_[1] << " , " << setw(9) << R_[2] << " )\n"
      << tab << "                       ( " << setw(9) << R_[3] << " , " << setw(9) << R_[4] << " , " << setw(9) << R_[5] << " )\n"
      << tab << "                       ( " << setw(9) << R_[6] << " , " << setw(9) << R_[7] << " , " << setw(9) << R_[8] << " )\n";

   if( isFinite() ) {
      os << std::setiosflags(std::ios::right)
         << tab << "   Moment of inertia = ( " << setw(9) << I_[0] << " , " << setw(9) << I_[1] << " , " << setw(9) << I_[2] << " )\n"
         << tab << "                       ( " << setw(9) << I_[3] << " , " << setw(9) << I_[4] << " , " << setw(9) << I_[5] << " )\n"
         << tab << "                       ( " << setw(9) << I_[6] << " , " << setw(9) << I_[7] << " , " << setw(9) << I_[8] << " )\n"
         << std::resetiosflags(std::ios::right);
   }

   // Extending the indentation for the nested bodies and links
   std::string longtab( tab );
   longtab.append( "  " );

   // Printing all contained bodies
   for( auto bIt=bodies_.begin(); bIt!=bodies_.end(); ++bIt ) {
      os << "\n";
      bIt->print( os, longtab.c_str() );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Global output operator for unions.
 *
 * \param os Reference to the output stream.
 * \param u Reference to a constant union object.
 * \return Reference to the output stream.
 */
template <typename BodyTypeTuple>
std::ostream& operator<<( std::ostream& os, const Union<BodyTypeTuple>& u )
{
   os << "--" << "UNION PARAMETERS"
      << "--------------------------------------------------------------\n";
   u.print( os, "" );
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for union handles.
 *
 * \param os Reference to the output stream.
 * \param u Constant union handle.
 * \return Reference to the output stream.
 */
template <typename BodyTypeTuple>
std::ostream& operator<<( std::ostream& os, Union<BodyTypeTuple> const * u )
{
   os << "--" << "UNION PARAMETERS"
      << "--------------------------------------------------------------\n";
   u->print( os, "" );
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************


//=================================================================================================
//
//  GET FUNCTIONS
//
//=================================================================================================


//*************************************************************************************************
/*!\brief Returns unique type id of this type
 *
 * \return geometry specific type id
 */
template <typename BodyTypeTuple>
inline id_t Union<BodyTypeTuple>::getStaticTypeID()
{
   return staticTypeID_;
}

template <typename BodyTypeTuple>
id_t Union<BodyTypeTuple>::staticTypeID_ = std::numeric_limits<id_t>::max();

} // namespace pe
}  // namespace walberla
