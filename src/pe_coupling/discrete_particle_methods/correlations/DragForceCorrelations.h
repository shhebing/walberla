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
//! \file DragForceCorrelations.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/Constants.h"

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {


/*!\brief Various correlation functions for the drag force.
 *
 * These functions calculate the drag force for fluid-particle interactions based on different empirical correlations
 * from literature.
 * Always be aware that those empirical formulas were obtained by different setups and are thus not generally applicable!
 *
 * To be compatible with the interface of InteractionForceEvaluator, all functions use the following signature:
 * Vector3<real_t> ( const Vector3<real_t> & fluidVel, const Vector3<real_t> & particleVel, real_t solidVolumeFraction,
 *                   real_t diameter, real_t fluidDynamicViscosity, real_t fluidDensity )
 *
 */


// helper functions -> often used in more advanced drag correlations

// equation to calculate the drag coefficient on isolated spherical particle
// Schiller, L., Naumann, A., 1935. A drag coefficient correlation. Vdi Zeitung 77, 318-320.
real_t dragCoeffSchillerNaumann( real_t reynoldsNumber )
{
   WALBERLA_ASSERT_GREATER_EQUAL( reynoldsNumber, 0_r );

   return ( reynoldsNumber < 1000_r ) ? 24_r * ( 1_r + 0.15_r * std::pow(reynoldsNumber, 0.687_r ) ) / reynoldsNumber
                                            : 0.44_r;
}

// Coefficient from Stokes' law for drag, only valid for Stokes regime (low Reynolds numbers)
// = 3 * PI * mu * D * fluidVolumeFraction
real_t dragCoeffStokes ( real_t fluidVolumeFraction, real_t diameter, real_t fluidDynamicViscosity )
{
   return 3_r * math::M_PI * diameter * fluidDynamicViscosity * fluidVolumeFraction;
}

// threshold value for absolute relative velocity
// if it is below this value, a drag force of 0 is set, to avoid instabilities stemming from divisions by this small value
const real_t thresholdAbsoluteVelocityDifference = 1e-10_r;


//////////////////////
//                  //
//   CORRELATIONS   //
//                  //
//////////////////////

// Stokes drag law
Vector3<real_t> dragForceStokes( const Vector3<real_t> & fluidVel, const Vector3<real_t> & particleVel,
                                 real_t solidVolumeFraction, real_t diameter, real_t fluidDynamicViscosity, real_t /*fluidDensity*/ )
{
   WALBERLA_ASSERT_GREATER_EQUAL( solidVolumeFraction, 0_r );
   WALBERLA_ASSERT_LESS_EQUAL( solidVolumeFraction, 1_r );

   Vector3<real_t> velDiff = fluidVel - particleVel;
   real_t absVelDiff = velDiff.length();

   if( absVelDiff < thresholdAbsoluteVelocityDifference ) return Vector3<real_t>(0_r);

   real_t fluidVolumeFraction = 1_r - solidVolumeFraction;

   return dragCoeffStokes( fluidVolumeFraction, diameter, fluidDynamicViscosity ) * velDiff;
}


// S. Ergun, Fluid flow through packed columns. Chemical Engineering Progress 48 (1952), 89-94.
// Y. C. Wen, Y.H. Yu, Mechanics of fluidization. Chemical Engineering Progress Symposium Series 62 (1966), 100-111.
// see also Beetstra, van der Hoef, Kuipers, "Drag Force of Intermediate Reynolds Number Flow Past Mono- and Bidisperse Arrays of Spheres" (2007)
Vector3<real_t>  dragForceErgunWenYu( const Vector3<real_t> & fluidVel, const Vector3<real_t> & particleVel,
                                      real_t solidVolumeFraction, real_t diameter, real_t fluidDynamicViscosity, real_t fluidDensity )
{
   WALBERLA_ASSERT_GREATER_EQUAL( solidVolumeFraction, 0_r );
   WALBERLA_ASSERT_LESS_EQUAL( solidVolumeFraction, 1_r );

   Vector3<real_t> velDiff = fluidVel - particleVel;
   real_t absVelDiff = velDiff.length();

   if( absVelDiff < thresholdAbsoluteVelocityDifference ) return Vector3<real_t>(0_r);

   real_t fluidVolumeFraction = 1_r - solidVolumeFraction;

   if( fluidVolumeFraction < 0.8_r )
   {
      // Ergun relation
      real_t reynoldsNumber = fluidVolumeFraction * fluidDensity * absVelDiff * diameter / fluidDynamicViscosity;
      real_t fDrag = 150_r * solidVolumeFraction / ( 18_r * fluidVolumeFraction * fluidVolumeFraction ) +
                     1.75_r / ( 18_r * fluidVolumeFraction * fluidVolumeFraction ) * reynoldsNumber;
      return fDrag * dragCoeffStokes( fluidVolumeFraction, diameter, fluidDynamicViscosity ) * velDiff;
   } else
   {
      // Wen & Yu correlation
      real_t reynoldsNumber = fluidVolumeFraction * fluidDensity * absVelDiff * diameter / fluidDynamicViscosity;
      real_t fDrag = dragCoeffSchillerNaumann( reynoldsNumber ) * reynoldsNumber / 24_r * std::pow( fluidVolumeFraction, -3.7_r );
      return fDrag * dragCoeffStokes( fluidVolumeFraction, diameter, fluidDynamicViscosity ) * velDiff;
   }
}

// drag correlation proposed by Tang et al. - "A New Drag Correlation from Fully Resolved Simulations of Flow Past
// Monodisperse Static Arrays of Spheres", AiChE, 2014
Vector3<real_t> dragForceTang( const Vector3<real_t> & fluidVel, const Vector3<real_t> & particleVel,
                               real_t solidVolumeFraction, real_t diameter, real_t fluidDynamicViscosity, real_t fluidDensity )
{
   WALBERLA_ASSERT_GREATER_EQUAL( solidVolumeFraction, 0_r );
   WALBERLA_ASSERT_LESS_EQUAL( solidVolumeFraction, 1_r );

   Vector3<real_t> velDiff = fluidVel - particleVel;
   real_t absVelDiff = velDiff.length();

   if( absVelDiff < thresholdAbsoluteVelocityDifference ) return Vector3<real_t>(0_r);

   real_t fluidVolumeFraction = 1_r - solidVolumeFraction;
   real_t fluidVolumeFractionP2 = fluidVolumeFraction * fluidVolumeFraction;
   real_t inv_fluidVolumeFractionP4 = 1_r / (fluidVolumeFractionP2 * fluidVolumeFractionP2);
   real_t reynoldsNumber = fluidVolumeFraction * fluidDensity * absVelDiff * diameter / fluidDynamicViscosity;

   // Eq.21 from the paper
   real_t fDrag = 10_r * solidVolumeFraction / fluidVolumeFractionP2 + fluidVolumeFractionP2 * ( 1_r + 1.5_r * std::sqrt(solidVolumeFraction) )
                + ( 0.11_r * solidVolumeFraction * ( 1_r + solidVolumeFraction ) - 0.00456_r * inv_fluidVolumeFractionP4
                + ( 0.169_r * fluidVolumeFraction + 0.0644_r * inv_fluidVolumeFractionP4 ) * std::pow( reynoldsNumber, -0.343_r ) ) * reynoldsNumber;

   return fDrag * dragCoeffStokes( fluidVolumeFraction, diameter, fluidDynamicViscosity ) * velDiff;

}


// drag correlation based on findings from Felice (1994)
// used e.g. in Kafui et al (2002)
Vector3<real_t> dragForceFelice( const Vector3<real_t> & fluidVel, const Vector3<real_t> & particleVel,
                                 real_t solidVolumeFraction, real_t diameter, real_t fluidDynamicViscosity, real_t fluidDensity )
{
   WALBERLA_ASSERT_GREATER_EQUAL( solidVolumeFraction, 0_r );
   WALBERLA_ASSERT_LESS_EQUAL( solidVolumeFraction, 1_r );

   Vector3<real_t> velDiff = fluidVel - particleVel;
   real_t absVelDiff = velDiff.length();

   if( absVelDiff < thresholdAbsoluteVelocityDifference ) return Vector3<real_t>(0_r);

   real_t fluidVolumeFraction = 1_r - solidVolumeFraction;

   real_t reynoldsNumber = fluidVolumeFraction * fluidDensity * absVelDiff * diameter / fluidDynamicViscosity;

   real_t temp1 = ( 0.63_r + 4.8_r / std::sqrt( reynoldsNumber ) );
   real_t dragCoeff = temp1 * temp1;

   real_t temp2 = 1.5_r - std::log10( reynoldsNumber );
   real_t chi = 3.7_r - std::pow( 0.65_r, (- 0.5_r * temp2 * temp2 ) );

   return 0.125_r * dragCoeff * fluidDensity * math::M_PI * diameter * diameter * absVelDiff *
          std::pow( fluidVolumeFraction, 2_r - chi) * velDiff;

}

// drag correlation based on findings from Tenneti, Garg, Subramaniam (2011)
// used e.g. in Finn, Li, Apte - Particle based modelling and simulation of natural sand dynamics in the wave bottom boundary layer (2016)
// could be generalized also for non-spherical particles, see Finn et al (2016)
Vector3<real_t> dragForceTenneti( const Vector3<real_t> & fluidVel, const Vector3<real_t> & particleVel,
                                  real_t solidVolumeFraction, real_t diameter, real_t fluidDynamicViscosity, real_t fluidDensity )
{
   WALBERLA_ASSERT_GREATER_EQUAL( solidVolumeFraction, 0_r );
   WALBERLA_ASSERT_LESS_EQUAL( solidVolumeFraction, 1_r );

   Vector3<real_t> velDiff = fluidVel - particleVel;
   const real_t absVelDiff = velDiff.length();

   if( absVelDiff < thresholdAbsoluteVelocityDifference ) return Vector3<real_t>(0_r);

   const real_t fluidVolumeFraction = 1_r - solidVolumeFraction;

   const real_t reynoldsNumber = fluidVolumeFraction * fluidDensity * absVelDiff * diameter / fluidDynamicViscosity;

   const real_t fvfCubed = fluidVolumeFraction * fluidVolumeFraction * fluidVolumeFraction;
   const real_t A = 5.81_r * solidVolumeFraction / fvfCubed + 0.48_r * std::cbrt( solidVolumeFraction ) / ( fvfCubed * fluidVolumeFraction );

   const real_t svfCubed = solidVolumeFraction * solidVolumeFraction * solidVolumeFraction;
   const real_t B = svfCubed * reynoldsNumber * ( 0.95_r + 0.61_r * svfCubed / ( fluidVolumeFraction * fluidVolumeFraction ) );

   // version from Finn et al.
   const real_t CdRe0Sphere = 1_r + 0.15_r *  std::pow( reynoldsNumber, 0.687_r );

   const real_t CdRe = fluidVolumeFraction * ( CdRe0Sphere / fvfCubed + A + B );

   return 3_r * math::M_PI * diameter * fluidDynamicViscosity * fluidVolumeFraction * CdRe * velDiff;

}


Vector3<real_t> noDragForce( const Vector3<real_t> & /*fluidVel*/, const Vector3<real_t> & /*particleVel*/,
                             real_t /*solidVolumeFraction*/, real_t /*diameter*/, real_t /*fluidDynamicViscosity*/, real_t /*fluidDensity*/ )
{
   return Vector3<real_t>(0_r);
}

} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
