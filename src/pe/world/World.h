#pragma once

#include "blockforest/BlockForest.h"
#include "core/DataTypes.h"
#include "core/config/Config.h"
#include "domain_decomposition/BlockDataID.h"
#include "pe/basic.h"

#include <fstream>

namespace walberla
{
namespace pe
{

class World
{
public:
  World(shared_ptr<BlockForest> forest, shared_ptr<BodyStorage> globalBodyStorage,
        BlockDataID storageID, BlockDataID ccdID, BlockDataID fcdID, boost::function<void(void)> syncCall,
        boost::function<void(void)> syncCallWithoutTT, Config::BlockHandle peConfig,  shared_ptr<WcTimingTree> timingTree = nullptr);

  void syncAfterParticleCreation();
  void sync() { syncCall_(); }
  void timestep(const real_t dt = 1) { collisionResolver_->timestep(dt); }
  void simulationStep(const real_t dt = 1)
  {
    timestep(dt);
    sync();
  }

  void serialize() {
    if(serializeFilename_ == "") return;
    blockForest_->saveBlockData(serializeFilename_, storageID_);
  };

  void serialize(std::string filename) {
    if(filename == "") return;
    blockForest_->saveBlockData(filename, storageID_);
  }

  const BlockDataID &getStorageID() const { return storageID_; }
  const BlockDataID &getFCDID() const { return fcdID_; }
  shared_ptr<cr::ICR> getCollisionResolver() { return collisionResolver_; }
  const Vec3 &getGlobalLinearAcceleration() const { return collisionResolver_->getGlobalLinearAcceleration(); };
  void setGlobalLinearAcceleration(const Vec3 &newAcceleration) { collisionResolver_->setGlobalLinearAcceleration(newAcceleration); };

  PlaneID createPlane(id_t uid, Vec3 normal, const Vec3 &gpos, MaterialID material = Material::find("iron"));
  SphereID createSphere(id_t uid, const Vec3 &gpos, real_t radius, MaterialID material = Material::find("iron"), bool global = false, bool communicating = true, bool infiniteMass = false);
  BoxID createBox(id_t uid, const Vec3& gpos, const Vec3& lengths, MaterialID material = Material::find("iron"), bool global = false, bool communicating = true, bool infiniteMass = false);
  CapsuleID createCapsule(id_t uid, const Vec3& gpos, const real_t radius, const real_t length, MaterialID material = Material::find("iron"), bool global = false, bool communicating = true, bool infiniteMass = false);

  template<typename functionType, typename... Rest>
  BodyID createBody(functionType creationFunction, const real_t long_axis_length, Rest... args) {
    if(long_axis_length/2 > longestDistanceFromParticleCenterToEdge_) {
      longestDistanceFromParticleCenterToEdge_ = long_axis_length/2;
    }

    return creationFunction(*globalBodyStorage_, *blockForest_, storageID_, args...);
  }


private:
  const shared_ptr<BlockForest> blockForest_;
  const shared_ptr<BodyStorage> globalBodyStorage_;
  const shared_ptr<WcTimingTree> timingTree_;

  const BlockDataID storageID_;
  const BlockDataID ccdID_;
  const BlockDataID fcdID_;

  shared_ptr<cr::ICR> collisionResolver_;

  real_t shortestDistanceInBlock_;
  real_t longestDistanceFromParticleCenterToEdge_;

  std::string serializeFilename_;

  const boost::function<void(void)> syncCall_;
  const boost::function<void(void)> syncCallWithoutTT_;
};

template <typename BodyTuple>
World createWorldFromConfig(shared_ptr<BlockForest> forest, Config::BlockHandle peConf, shared_ptr<WcTimingTree> timingTree) {

  SetBodyTypeIDs<BodyTuple>::execute();

  auto globalBodyStorage = make_shared<BodyStorage>();

  std::string resumeString = peConf.getParameter<std::string>("ResumeFilename", "");

  BlockDataID storageID;

  // Support suspend resume
  if(resumeString != "") {
    std::ifstream f(resumeString.c_str());
    if(f.good()) {
      f.close();
      storageID = forest->loadBlockData(resumeString.c_str(), createStorageDataHandling<BodyTuple>());
    }
    else {
      f.close();
      storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
    }
  }
  else {
    storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
  }

  //auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
  auto ccdID = forest->addBlockData(ccd::createHashGridsDataHandling(globalBodyStorage, storageID), "CCD");
  auto fcdID = forest->addBlockData(fcd::createSimpleFCDDataHandling<BodyTuple>(), "FCD");

   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt) { ccd::ICCD* ccd = blockIt->getData< ccd::ICCD >( ccdID ); ccd->reloadBodies(); }

  boost::function<void(void)> syncCall;
  boost::function<void(void)> syncCallWithoutTT;

  bool largeParticleSupport = peConf.getParameter<bool>("LargeParticleSupport", false);
  if (!largeParticleSupport)
  {
    syncCall = boost::bind(pe::syncNextNeighbors<BodyTuple>, boost::ref(*forest), storageID, timingTree.get(), real_c(0.0), false);
  }
  else
  {
    syncCall = boost::bind(pe::syncShadowOwners<BodyTuple>, boost::ref(*forest), storageID, timingTree.get(), real_c(0.0), false);
  }

  if (!largeParticleSupport)
  {
    syncCallWithoutTT = boost::bind(pe::syncNextNeighbors<BodyTuple>, boost::ref(*forest), storageID, static_cast<WcTimingTree *>(NULL), real_c(0.0), false);
  }
  else
  {
    syncCallWithoutTT = boost::bind(pe::syncShadowOwners<BodyTuple>, boost::ref(*forest), storageID, static_cast<WcTimingTree *>(NULL), real_c(0.0), false);
  }

  return World(forest, globalBodyStorage, storageID, ccdID, fcdID, syncCall, syncCallWithoutTT, peConf, timingTree);
}

} // namespace pe
} // namespace walberla
