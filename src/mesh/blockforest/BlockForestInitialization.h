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
//! \file Initialization.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/loadbalancing/StaticCurve.h"

#include "core/math/AABB.h"
#include "core/math/IntegerFactorization.h"
#include "core/math/Primes.h"

#include <cmath>
#include <map>

namespace walberla {
namespace mesh {

class ComplexGeometryBlockforestCreator
{
public:
   ComplexGeometryBlockforestCreator( const AABB & aabb, const blockforest::SetupBlockForest::RootBlockExclusionFunction & rootBlockExclusionFunction );

   void setWorkloadMemorySUIDAssignmentFunction( const blockforest::SetupBlockForest::WorkloadMemorySUIDAssignmentFunction & workloadMemorySUIDAssignmentFunction ) { workloadMemorySUIDAssignmentFunction_ = workloadMemorySUIDAssignmentFunction; }
   void setTargetProcessAssignmentFunction( const blockforest::SetupBlockForest::TargetProcessAssignmentFunction & targetProcessAssignmentFunction ) { targetProcessAssignmentFunction_ = targetProcessAssignmentFunction; }
   void setProcessMemoryLimit( const real_t processMemoryLimit ) { processMemoryLimit_ = processMemoryLimit; }
   void setMaxBlockSkewness( const real_t maxBlockSkewness ) { maxBlockSkewness_ = maxBlockSkewness; }

   void setPeriodicity( const Vector3< bool > & periodicity ) { periodicity_ = periodicity; }

   shared_ptr<SetupBlockForest>      createSetupBlockForest     ( const uint_t targetNumRootBlocks, const uint_t numProcesses = uint_c( MPIManager::instance()->numProcesses() ) ) const;
   shared_ptr<SetupBlockForest>      createSetupBlockForest     ( const Vector3<real_t> & blockSize, const uint_t numProcesses = uint_c( MPIManager::instance()->numProcesses() ) ) const;

   shared_ptr<BlockForest>           createBlockForest          ( const uint_t targetNumRootBlocks ) const;
   shared_ptr<BlockForest>           createBlockForest          ( const Vector3<real_t> & blockSize ) const;

private:

   uint_t findNumBlocks( const Vector3<uint_t> & numRootBlocks3D ) const;

   AABB aabb_;
   uint_t maxIterations_;
   real_t acceptableRelativeError_;
   real_t maxBlockSkewness_; // maximum allowed factor between shortest and longest block edge
   real_t processMemoryLimit_;
   Vector3< bool > periodicity_;

   blockforest::SetupBlockForest::RootBlockExclusionFunction           rootBlockExclusionFunction_;
   blockforest::SetupBlockForest::WorkloadMemorySUIDAssignmentFunction workloadMemorySUIDAssignmentFunction_;
   blockforest::SetupBlockForest::TargetProcessAssignmentFunction      targetProcessAssignmentFunction_;
};


class ComplexGeometryStructuredBlockforestCreator
{
public:
   ComplexGeometryStructuredBlockforestCreator( const AABB & aabb, const Vector3<real_t> & cellSize, const blockforest::SetupBlockForest::RootBlockExclusionFunction & rootBlockExclusionFunction );

   void setWorkloadMemorySUIDAssignmentFunction( const blockforest::SetupBlockForest::WorkloadMemorySUIDAssignmentFunction & workloadMemorySUIDAssignmentFunction ) { workloadMemorySUIDAssignmentFunction_ = workloadMemorySUIDAssignmentFunction; }
   void setTargetProcessAssignmentFunction( const blockforest::SetupBlockForest::TargetProcessAssignmentFunction & targetProcessAssignmentFunction ) { targetProcessAssignmentFunction_ = targetProcessAssignmentFunction; }
   void setRefinementSelectionFunction( const blockforest::SetupBlockForest::RefinementSelectionFunction & refinementSelectionFunction ) { refinementSelectionFunction_ = refinementSelectionFunction; }
   void setProcessMemoryLimit( const real_t processMemoryLimit ) { processMemoryLimit_ = processMemoryLimit; }

   void setPeriodicity( const Vector3< bool > & periodicity ) { periodicity_ = periodicity; }

   shared_ptr<SetupBlockForest>      createSetupBlockForest     ( const uint_t targetNumRootBlocks, const uint_t numProcesses = uint_c( MPIManager::instance()->numProcesses() ) ) const;
   shared_ptr<SetupBlockForest>      createSetupBlockForest     ( const Vector3<uint_t> & blockSize, const uint_t numProcesses = uint_c( MPIManager::instance()->numProcesses() ) ) const;

   shared_ptr<StructuredBlockForest> createStructuredBlockForest( const uint_t targetNumRootBlocks ) const;
   shared_ptr<StructuredBlockForest> createStructuredBlockForest( const Vector3<uint_t> & blockSize ) const;



private:

   uint_t findNumBlocks( const Vector3<uint_t> & blockSize ) const;

   Vector3<uint_t> domainSizeCells() const { return  Vector3<uint_t>( uint_c( std::ceil( aabb_.xSize() / cellSize_[0] ) ), uint_c( std::ceil( aabb_.ySize() / cellSize_[1] ) ), uint_c( std::ceil( aabb_.zSize() / cellSize_[2] ) ) ); }

   AABB aabb_;
   Vector3<real_t> cellSize_;
   uint_t maxIterations_;
   real_t acceptableRelativeError_;
   real_t processMemoryLimit_;
   Vector3< bool > periodicity_;

   blockforest::SetupBlockForest::RootBlockExclusionFunction           rootBlockExclusionFunction_;
   blockforest::SetupBlockForest::WorkloadMemorySUIDAssignmentFunction workloadMemorySUIDAssignmentFunction_;
   blockforest::SetupBlockForest::TargetProcessAssignmentFunction      targetProcessAssignmentFunction_;
   blockforest::SetupBlockForest::RefinementSelectionFunction          refinementSelectionFunction_;
};


} // namespace mesh
} // namespace walberla