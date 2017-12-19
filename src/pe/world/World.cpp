#include "World.h"

#include "pe/cr/HCSITS.h"

namespace walberla
{
namespace pe
{



World::World(shared_ptr<BlockForest> forest, shared_ptr<BodyStorage> globalBodyStorage,
      BlockDataID storageID, BlockDataID ccdID, BlockDataID fcdID, boost::function<void(void)> syncCall,
      boost::function<void(void)> syncCallWithoutTT, Config::BlockHandle peConfig, shared_ptr<WcTimingTree> timingTree) :
      blockForest_(forest), globalBodyStorage_(globalBodyStorage), timingTree_(timingTree), storageID_(storageID), ccdID_(ccdID), fcdID_(fcdID), syncCall_(syncCall), syncCallWithoutTT_(syncCallWithoutTT)
{

  shortestDistanceInBlock_ = 0;
  longestDistanceFromParticleCenterToEdge_ = 0;

  real_t blockXSize = forest->getRootBlockXSize();
  real_t blockYSize = forest->getRootBlockYSize();
  real_t blockZSize = forest->getRootBlockZSize();

  if (blockXSize <= blockYSize && blockXSize <= blockZSize)
  {
    shortestDistanceInBlock_ = blockXSize;
  }
  else if (blockYSize <= blockZSize)
  {
    shortestDistanceInBlock_ = blockYSize;
  }
  else
  {
    shortestDistanceInBlock_ = blockZSize;
  }

  serializeFilename_ = peConfig.getParameter<std::string>("SuspendFilename", "");

  std::string collisionResolverString = peConfig.getParameter<std::string>("CollisionResolver", "DEM");

  Config::BlockHandle collisionConf = peConfig.getBlock("CollisonResolverConfig");

  if (collisionResolverString == "DEM")
  {
    shared_ptr<cr::DEM> cr = make_shared<cr::DEM>(globalBodyStorage_, forest, storageID_, ccdID_, fcdID_, timingTree_.get());
    collisionResolver_ = cr;
  }
  else if (collisionResolverString == "HCSITS")
  {
    shared_ptr<cr::HCSITS> cr = make_shared<cr::HCSITS>(globalBodyStorage_, forest, storageID_, ccdID_, fcdID_, timingTree_.get());
    configure(collisionConf, *cr);
    collisionResolver_ = cr;
  }
  else
  {
    WALBERLA_ABORT("Solver not implemented in World" << collisionResolverString);
  }

  Vec3 globalLinearAcceleration = peConfig.getParameter<Vec3>("globalLinearAcceleration", Vec3(0, 0, 0));
  WALBERLA_LOG_INFO_ON_ROOT("globalLinearAcceleration: " << globalLinearAcceleration);
  collisionResolver_->setGlobalLinearAcceleration(globalLinearAcceleration);
}


void World::syncAfterParticleCreation()
{
  mpi::allReduceInplace(longestDistanceFromParticleCenterToEdge_, mpi::MAX);
  uint_t maximumNumberOfSpannedBlocks = uint_c(ceil(longestDistanceFromParticleCenterToEdge_ / shortestDistanceInBlock_));
  if (maximumNumberOfSpannedBlocks < 1)
  {
    maximumNumberOfSpannedBlocks = 1;
  }
  for (uint_t i = 0; i < maximumNumberOfSpannedBlocks; ++i)
  {
    syncCallWithoutTT_();
  }
}


PlaneID World::createPlane(id_t uid, Vec3 normal, const Vec3 &gpos, MaterialID material)
{
  return pe::createPlane(*globalBodyStorage_, uid, normal, gpos, material);
}


SphereID World::createSphere(id_t uid, const Vec3 &gpos, real_t radius, MaterialID material, bool global, bool communicating, bool infiniteMass)
{
  if (radius > longestDistanceFromParticleCenterToEdge_)
  {
    longestDistanceFromParticleCenterToEdge_ = radius;
  }
  return pe::createSphere(*globalBodyStorage_, *blockForest_, storageID_, uid, gpos, radius, material, global, communicating, infiniteMass);
}


BoxID World::createBox(id_t uid, const Vec3& gpos, const Vec3& lengths, MaterialID material, bool global, bool communicating, bool infiniteMass)
{
  for(uint_t i = 0; i < 3; i++) {
    if ((lengths[i]/2.0) > longestDistanceFromParticleCenterToEdge_)
    {
      longestDistanceFromParticleCenterToEdge_ = real_c(lengths[i]/2.0);
    }
  }

  return pe::createBox(*globalBodyStorage_, *blockForest_, storageID_, uid, gpos, lengths, material, global, communicating, infiniteMass);
}

CapsuleID World::createCapsule(id_t uid, const Vec3& gpos, const real_t radius, const real_t length, MaterialID material, bool global, bool communicating, bool infiniteMass)
{

    if((length/2.0 + radius) > longestDistanceFromParticleCenterToEdge_)
    {
      longestDistanceFromParticleCenterToEdge_ = real_c(length/2.0 + radius);
    }
    return pe::createCapsule(*globalBodyStorage_, *blockForest_, storageID_, uid, gpos, radius, length, material, global, communicating, infiniteMass);
}




} // namespace pe
} // namespace walberla
