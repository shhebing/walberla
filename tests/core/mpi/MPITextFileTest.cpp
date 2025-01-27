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
//! \file MPITextFileText.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "core/Abort.h"
#include "core/DataTypes.h"

#include "core/debug/TestSubsystem.h"

#include "core/mpi/MPIManager.h"
#include "core/mpi/MPITextFile.h"

#include "core/Filesystem.h"

#include <vector>
#include <sstream>
#include <fstream>


void testSameSizeFile( const std::string & filename, const size_t chunkSize )
{
   using namespace walberla;

   WALBERLA_CHECK_GREATER( chunkSize, 0 );

   const int rank = MPIManager::instance()->rank();
   std::ostringstream oss;
   oss << rank;

   std::string chunk( chunkSize, char( 'A' + static_cast<char>( rank % 26 ) ) );
   chunk[ chunk.size() - size_t(1) ] = '\n';
   
   mpi::writeMPITextFile( filename, chunk );

   std::ifstream ifs( filename.c_str() );
   ifs.seekg( numeric_cast< std::ifstream::off_type >( uint_c(rank) * chunkSize ), std::ios_base::beg );
   std::vector<char> buffer( chunkSize );
   ifs.read( &(buffer[0]), numeric_cast< std::streamsize >( chunkSize ) );
   ifs.close();
   std::string referenceChunk( buffer.begin(), buffer.end() );

   WALBERLA_CHECK_EQUAL( chunk, referenceChunk );

   WALBERLA_MPI_BARRIER();
   WALBERLA_ROOT_SECTION()
   {
      if( filesystem::exists( filename ) )
         filesystem::remove( filename );
   }
   WALBERLA_MPI_BARRIER();
}



void testDifferentSizeFile( const std::string & filename, const size_t minChunkSize )
{
   using namespace walberla;

   WALBERLA_CHECK_GREATER( minChunkSize, 0 );

   const int rank = MPIManager::instance()->rank();
   std::ostringstream oss;
   oss << rank;

   const size_t chunkSize = minChunkSize * uint_c(rank + 1);

   std::string chunk(chunkSize , char( 'A' + static_cast<char>( rank % 26 ) ) );
   chunk[ chunk.size() - size_t(1) ] = '\n';

   mpi::writeMPITextFile( filename, chunk );

   std::ifstream ifs( filename.c_str() );
   ifs.seekg( numeric_cast< std::ifstream::off_type >( uint_c( ( rank * rank + rank ) / 2 ) * minChunkSize ), std::ios_base::beg );
   std::vector<char> buffer( chunkSize );
   ifs.read( &(buffer[0]), numeric_cast< std::streamsize >( chunkSize ) );
   ifs.close();
   std::string referenceChunk( buffer.begin(), buffer.end() );

   WALBERLA_CHECK_EQUAL( chunk, referenceChunk );

   WALBERLA_MPI_BARRIER();
   WALBERLA_ROOT_SECTION()
   {
      if( filesystem::exists( filename ) )
         filesystem::remove( filename );
   }
   WALBERLA_MPI_BARRIER();
}



int main( int argc, char * argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->useWorldComm();

   std::vector<std::string> args( argv, argv + argc );

   size_t      chunkSize;
   std::string filename;
   try
   {
      chunkSize = walberla::string_to_num<size_t>( args.at(2) );
      filename  = args.at( 1 );
   }
   catch( ... )
   {
      WALBERLA_ABORT_NO_DEBUG_INFO( "Usage:\n" << args[0] << " FILENAME CHUNK_SIZE" );
   }

   testSameSizeFile( filename, chunkSize );
   testDifferentSizeFile( filename, chunkSize );
}
