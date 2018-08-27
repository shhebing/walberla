#include <pe/basic.h>
#include "pe/utility/BodyCast.h"

#include "pe/Materials.h"

#include "pe/rigidbody/Box.h"
#include "pe/rigidbody/Capsule.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Union.h"
#include "pe/rigidbody/Ellipsoid.h"

#include "pe/rigidbody/SetBodyTypeIDs.h"
#include "pe/Types.h"

#include "core/debug/TestSubsystem.h"
#include "core/DataTypes.h"
#include <core/math/Random.h>
#include "core/math/Vector3.h"

#include <pe/raytracing/Ray.h>
#include <pe/raytracing/Intersects.h>
#include <pe/raytracing/Raytracer.h>
#include <pe/raytracing/Color.h>
#include <pe/raytracing/ShadingFunctions.h>

#include <pe/ccd/HashGrids.h>
#include "pe/rigidbody/BodyStorage.h"
#include <core/timing/TimingTree.h>

#include <pe/utility/GetBody.h>

#include <sstream>
#include <tuple>

namespace walberla {
using namespace walberla::pe;
using namespace walberla::pe::raytracing;

typedef boost::tuple<Box, Plane, Sphere, Capsule, Ellipsoid> BodyTuple ;

void SphereIntersectsTest()
{
   MaterialID iron = Material::find("iron");
   Sphere sp1(123, 1, Vec3(3,3,3), Vec3(0,0,0), Quat(), 2, iron, false, true, false);
   real_t t;
   Vec3 n;
   
   // ray through the center
   Ray ray1(Vec3(3,-5,3), Vec3(0,1,0));
   WALBERLA_LOG_INFO("RAY -> SPHERE");
   
   WALBERLA_CHECK(intersects(&sp1, ray1, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL(t, 6_r);
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], 0_r);
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], -1_r);
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], 0_r);
   
   // ray tangential
   Ray ray2(Vec3(3,-5,3), Vec3(0,7.5,real_t(std::sqrt(15_r)/2_r)).getNormalized());
   WALBERLA_CHECK(intersects(&sp1, ray2, t, n));
   
   // sphere behind ray origin
   Sphere sp2(123, 1, Vec3(3,-8,3), Vec3(0,0,0), Quat(), 2, iron, false, true, false);
   WALBERLA_CHECK(!intersects(&sp2, ray1, t, n));
   
   // sphere around ray origin
   Sphere sp3(123, 1, Vec3(3,-5,3), Vec3(0,0,0), Quat(), 2, iron, false, true, false);
   WALBERLA_CHECK(intersects(&sp3, ray1, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL(t, 2_r);
}

void PlaneIntersectsTest() {
   MaterialID iron = Material::find("iron");
   // plane with center 3,3,3 and parallel to y-z plane
   Plane pl1(1, 1, Vec3(3, 3, 3), Vec3(1, 0, 0), 1.0_r, iron);
   
   Ray ray1(Vec3(-5,3,3), Vec3(1,0,0));
   real_t t;
   Vec3 n;
   
   WALBERLA_LOG_INFO("RAY -> PLANE");
   WALBERLA_CHECK(intersects(&pl1, ray1, t, n), "ray through center did not hit");
   WALBERLA_CHECK_FLOAT_EQUAL(t, 8_r, "distance between ray and plane is incorrect");
   
   Ray ray2(Vec3(-5,3,3), Vec3(1,0,-1).getNormalized());
   WALBERLA_CHECK(intersects(&pl1, ray2, t, n), "ray towards random point on plane didn't hit");
   WALBERLA_CHECK_FLOAT_EQUAL(t, real_t(sqrt(128_r)), "distance between ray and plane is incorrect");
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], -1_r, "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], 0_r, "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], 0_r, "incorrect normal calculated");
   
   Plane pl1neg(1, 1, Vec3(3, 3, 3), Vec3(-1, 0, 0), 1.0_r, iron);
   WALBERLA_CHECK(intersects(&pl1neg, ray2, t, n), "ray towards random point on plane didn't hit");
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], -1_r, "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], 0_r, "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], 0_r, "incorrect normal calculated");
   
   Ray ray3(Vec3(-5,3,3), Vec3(-1,0,0).getNormalized());
   Plane pl5(1, 1, Vec3(-7, 3, 3), Vec3(1, 0, 0), 1.0_r, iron);
   WALBERLA_CHECK(intersects(&pl5, ray3, t, n), "ray towards random point on plane didn't hit");
   WALBERLA_CHECK_FLOAT_EQUAL(t, 2_r, "distance between ray and plane is incorrect");
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], 1_r, "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], 0_r, "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], 0_r, "incorrect normal calculated");
   
   // plane with center 3,3,3 and parallel to x-z plane
   Plane pl2(1, 1, Vec3(3, 3, 3), Vec3(0, 1, 0), 1.0_r, iron);
   WALBERLA_CHECK(!intersects(&pl2, ray1, t, n), "ray parallel to plane shouldnt hit");
   
   // plane with center -10,3,3 and parallel to y-z plane
   Plane pl4(1, 1, Vec3(-10, 3, 3), Vec3(1, 0, 0), 1.0_r, iron);
   WALBERLA_CHECK(!intersects(&pl4, ray1, t, n), "ray hit plane behind origin");
   
   Plane pl6(1, 1, Vec3(3, 3, 0), Vec3(-1, 0, 0), 1.0_r, iron);
   Ray ray4(Vec3(0,0,5), Vec3(1, 0, -1).getNormalized());
   WALBERLA_CHECK(intersects(&pl6, ray4, t, n), "ray didnt hit");
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], -1_r, "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], 0_r, "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], 0_r, "incorrect normal calculated");
}

void BoxIntersectsTest() {
   WALBERLA_LOG_INFO("RAY -> BOX");
   
   MaterialID iron = Material::find("iron");
   real_t t;
   Vec3 n;
   
   Box box1(127, 5, Vec3(0, -15, 0), Vec3(0, 0, 0), Quat(), Vec3(10, 10, 10), iron, false, true, false);
   Ray ray1(Vec3(3,-5,3), Vec3(0,1,0));
   WALBERLA_CHECK(!intersects(&box1, ray1, t, n));
   
   Box box2(128, 5, Vec3(0, -2, 0), Vec3(0, 0, 0), Quat(), Vec3(10, 10, 10), iron, false, true, false);
   WALBERLA_CHECK(intersects(&box2, ray1, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(t, 8_r, 1e-7_r);
   
   Box box3(128, 5, Vec3(0, 5, 0), Vec3(0, 0, 0), Quat(), Vec3(10, 10, 10), iron, false, true, false);
   WALBERLA_CHECK(intersects(&box3, ray1, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL(t, 5_r);
   
   Ray ray6(Vec3(-8,5,0), Vec3(1,0,0));
   WALBERLA_CHECK(intersects(&box3, ray6, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL(t, 3_r);
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], -1_r, "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], 0_r, "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], 0_r, "incorrect normal calculated");
   
   Ray ray7(Vec3(8,5,0), Vec3(-1,0,0));
   WALBERLA_CHECK(intersects(&box3, ray7, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL(t, 3_r);
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], 1_r, "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], 0_r, "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], 0_r, "incorrect normal calculated");
   
   // ray origin within box
   Ray ray2(Vec3(-2,0,0), Vec3(1,0,1).getNormalized());
   WALBERLA_CHECK(intersects(&box3, ray2, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(t, 7.0710_r, 1e-4_r);
   
   Ray ray3(Vec3(3,-5,3), Vec3(2, -1.5_r, 0.5_r).getNormalized());
   Box box4(128, 5, Vec3(0, 8, 0), Vec3(0, 0, 0), Quat(), Vec3(10, 10, 10), iron, false, true, false);
   WALBERLA_CHECK(!intersects(&box4, ray3, t, n));
   
   Ray ray4(Vec3(3,-5,3), Vec3(-2, 3, 0.5_r).getNormalized());
   WALBERLA_CHECK(intersects(&box4, ray4, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(t, 9.7068_r, 1e-4_r);
   
   Box box5(128, 5, Vec3(4, 0, 0), Vec3(0, 0, 0), Quat(), Vec3(4, 4, 4), iron, false, true, false);
   box5.rotate(0,0,math::M_PI/4);
   Ray ray5(Vec3(0,1.5,0), Vec3(1,0,0));
   WALBERLA_CHECK(intersects(&box5, ray5, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(t, 2.67157_r, 1e-4_r);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(n[0], -0.707107_r, 1e-5_r, "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(n[1], 0.707107_r, 1e-5_r, "incorrect normal calculated");
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], 0_r, "incorrect normal calculated");
}

void AABBIntersectsTest() {
   WALBERLA_LOG_INFO("RAY -> AABB");
   
   Ray ray1(Vec3(-5,5,5), Vec3(1,0,0));
   real_t t;
   
   AABB aabb(0,0,0,
             10,10,10);
   
   WALBERLA_CHECK(intersects(aabb, ray1, t));
   WALBERLA_CHECK_FLOAT_EQUAL(t, 5_r);
   
   WALBERLA_CHECK(intersects(aabb, ray1, t, 1.0));
   WALBERLA_CHECK_FLOAT_EQUAL(t, 4_r);
   
   Ray ray2(Vec3(-5,5,10.5), Vec3(1,0,0)); // ray shooting over aabb, but within padding passed to intersects
   WALBERLA_CHECK(intersects(aabb, ray1, t, 1.0));
   WALBERLA_CHECK_FLOAT_EQUAL(t, 4_r);
}

void CapsuleIntersectsTest() {
   MaterialID iron = Material::find("iron");
   real_t t;
   Vec3 n;
   
   Capsule cp1(0, 0, Vec3(2,3,3), Vec3(0,0,0), Quat(), 2_r, 2_r, iron, false, true, false);
   
   // ray through the center
   Ray ray1(Vec3(3,-5,3), Vec3(0,1,0));
   WALBERLA_LOG_INFO("RAY -> CAPSULE");
   
   WALBERLA_CHECK(intersects(&cp1, ray1, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL(t, 6_r);
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], 0_r);
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], -1_r);
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], 0_r);
   
   Ray ray2(Vec3(-5,3,3), Vec3(1,0,0));
   WALBERLA_CHECK(intersects(&cp1, ray2, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL(t, 4_r);
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], -1_r);
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], 0_r);
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], 0_r);
}

void EllipsoidTest() {
   WALBERLA_LOG_INFO("RAY -> ELLIPSOID");
   
   MaterialID iron = Material::find("iron");
   real_t t;
   Vec3 n;
   
   Ellipsoid el1(0,0, Vec3(2,3,3), Vec3(0,0,0), Quat(), Vec3(2,3,1), iron, false, true, false);
   
   Ray ray1(Vec3(-2,3,3), Vec3(1,0,0).getNormalized());
   WALBERLA_CHECK(intersects(&el1, ray1, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL(t, 2_r);
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], -1_r);
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], 0_r);
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], 0_r);

   Ray ray2(Vec3(2,3,0), Vec3(0,0,-1).getNormalized());
   WALBERLA_CHECK(!intersects(&el1, ray2, t, n));
   
   Ray ray3(Vec3(2,3,5), Vec3(0,0,-1).getNormalized());
   WALBERLA_CHECK(intersects(&el1, ray3, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL(t, 1_r);
   WALBERLA_CHECK_FLOAT_EQUAL(n[0], 0_r);
   WALBERLA_CHECK_FLOAT_EQUAL(n[1], 0_r);
   WALBERLA_CHECK_FLOAT_EQUAL(n[2], 1_r);

   Ray ray4(Vec3(-2,2_r,2_r), Vec3(1,0_r,0.5_r).getNormalized());
   WALBERLA_CHECK(intersects(&el1, ray4, t, n));
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(t, 2.36809_r, 1e-5_r);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(n[0], -0.78193_r, 1e-5_r);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(n[1], -0.62324_r, 1e-5_r);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(n[2], 0.012265_r, 1e-5_r);
}

ShadingParameters customBodyToShadingParams(const BodyID body) {
   if (body->getID() == 10) {
      return greenShadingParams(body).makeGlossy(30);
   } else if (body->getID() == 7) {
      return greenShadingParams(body).makeGlossy(10);
   } else if (body->getID() == 9) {
      return darkGreyShadingParams(body).makeGlossy(50);
   } else if (body->getID() == 3) {
      return redShadingParams(body).makeGlossy(200);
   } else {
      return defaultBodyTypeDependentShadingParams(body);
   }
}

void RaytracerTest(Raytracer::Algorithm raytracingAlgorithm = Raytracer::RAYTRACE_HASHGRIDS, walberla::uint8_t antiAliasFactor = 1) {
   WALBERLA_LOG_INFO("Raytracer");
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   shared_ptr<BlockForest> forest = createBlockForest(AABB(0,0,0,10,10,10), Vector3<uint_t>(1,1,1), Vector3<bool>(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID = forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "CCD");
   
   Lighting lighting(Vec3(0, 5, 8), // 8, 5, 9.5 gut für ebenen, 0,5,8
                     Color(1, 1, 1), //diffuse
                     Color(1, 1, 1), //specular
                     Color(0.4_r, 0.4_r, 0.4_r)); //ambient
   Raytracer raytracer(forest, storageID, globalBodyStorage, ccdID,
                       size_t(640), size_t(480),
                       49.13_r, antiAliasFactor,
                       Vec3(-5,5,5), Vec3(-1,5,5), Vec3(0,0,1), //-5,5,5; -1,5,5
                       lighting,
                       Color(0.2_r, 0.2_r, 0.2_r),
                       customBodyToShadingParams);
   
   MaterialID iron = Material::find("iron");
   
   //PlaneID xNegPlane = createPlane(*globalBodyStorage, 0, Vec3(-1,0,0), Vec3(5,0,0), iron);
   // xNegPlane obstructs only the top left sphere and intersects some objects
   
   //PlaneID xNegPlaneClose = createPlane(*globalBodyStorage, 0, Vec3(-1,0,0), Vec3(1,0,0), iron);
   
   // Test Scene v1 - Spheres, (rotated) boxes, confining walls, tilted plane in right bottom back corner
   createPlane(*globalBodyStorage, 0, Vec3(0,-1,0), Vec3(0,10,0), iron); // left wall
   createPlane(*globalBodyStorage, 0, Vec3(0,1,0), Vec3(0,0,0), iron); // right wall
   createPlane(*globalBodyStorage, 0, Vec3(0,0,1), Vec3(0,0,0), iron); // floor
   createPlane(*globalBodyStorage, 0, Vec3(0,0,-1), Vec3(0,0,10), iron); // ceiling
   createPlane(*globalBodyStorage, 0, Vec3(-1,0,0), Vec3(10,0,0), iron); // back wall
   createPlane(*globalBodyStorage, 0, Vec3(1,0,0), Vec3(0,0,0), iron); // front wall, should not get rendered
   
   createPlane(*globalBodyStorage, 0, Vec3(-1,1,1), Vec3(8,2,2), iron); // tilted plane in right bottom back corner
   
   createSphere(*globalBodyStorage, *forest, storageID, 2, Vec3(6,9.5_r,9.5_r), 0.5_r);
   createSphere(*globalBodyStorage, *forest, storageID, 3, Vec3(4,5.5_r,5), 1_r);
   createSphere(*globalBodyStorage, *forest, storageID, 6, Vec3(3,8.5_r,5), 1_r);
   BoxID box = createBox(*globalBodyStorage, *forest, storageID, 7, Vec3(5,6.5_r,5), Vec3(2,4,3));
   if (box != nullptr) box->rotate(0,math::M_PI/4,math::M_PI/4);
   createBox(*globalBodyStorage, *forest, storageID, 8, Vec3(5,1,8), Vec3(2,2,2));
   // Test scene v1 end
   
   // Test scene v2 additions start
   createBox(*globalBodyStorage, *forest, storageID, 9, Vec3(9,9,5), Vec3(1,1,10));
   createCapsule(*globalBodyStorage, *forest, storageID, 10, Vec3(3, 9, 1), 0.5_r, 7_r, iron);
   CapsuleID capsule = createCapsule(*globalBodyStorage, *forest, storageID, 11, Vec3(7, 3.5_r, 7.5_r), 1_r, 2_r, iron);
   if (capsule != nullptr) capsule->rotate(0,math::M_PI/3,math::M_PI/4-math::M_PI/8);
   // Test scene v2 end
   
   // Test scene v3 additions start
   EllipsoidID ellipsoid = createEllipsoid(*globalBodyStorage, *forest, storageID, 12, Vec3(6,2,2.5_r), Vec3(3,2,1.2_r));
   ellipsoid->rotate(0, math::M_PI/6_r, 0);
   // Test scene v3 end
   
   //raytracer.setTBufferOutputDirectory("tbuffer");
   //raytracer.setTBufferOutputEnabled(true);
   raytracer.setImageOutputEnabled(true);
   //raytracer.setLocalImageOutputEnabled(true);
   
   raytracer.setRaytracingAlgorithm(raytracingAlgorithm);
   raytracer.generateImage<BodyTuple>(0);
}

ShadingParameters customSpheresBodyToShadingParams(const BodyID body) {
   if (body->getTypeID() == Plane::getStaticTypeID()) {
      return greyShadingParams(body);
   }
   
   switch (body->getID()) {
      case 0:
         return blueShadingParams(body).makeGlossy(1);
      case 1:
         return blueShadingParams(body).makeGlossy(10);
      case 2:
         return blueShadingParams(body).makeGlossy(30);
      case 3:
         return blueShadingParams(body).makeGlossy(80);
      case 4:
         return whiteShadingParams(body);
      case 5:
         return lightGreyShadingParams(body);
      case 6:
         return greyShadingParams(body);
      case 7:
         return darkGreyShadingParams(body);
      case 8:
         return blackShadingParams(body).makeGlossy(100);
      case 9:
         return redShadingParams(body);
      case 10:
         return blueShadingParams(body);
      case 11:
         return violetShadingParams(body);
      case 12:
         return greenShadingParams(body);
      case 13:
         return greenShadingParams(body).makeGlossy(30);
      case 14:
         return blueShadingParams(body).makeGlossy(1000);
      default:
         return lightGreyShadingParams(body);
   }
}

void RaytracerSpheresTestScene(Raytracer::Algorithm raytracingAlgorithm = Raytracer::RAYTRACE_HASHGRIDS, walberla::uint8_t antiAliasFactor = 1) {
   WALBERLA_LOG_INFO("Raytracer Spheres Scene");
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   shared_ptr<BlockForest> forest = createBlockForest(AABB(0,0,0,10,10,10), Vector3<uint_t>(1,1,1), Vector3<bool>(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID = forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "CCD");
   
   Lighting lighting(Vec3(0, 5, 8), // 8, 5, 9.5 gut für ebenen, 0,5,8
                     Color(1, 1, 1), //diffuse
                     Color(1, 1, 1), //specular
                     Color(0.4_r, 0.4_r, 0.4_r)); //ambient
   Raytracer raytracer(forest, storageID, globalBodyStorage, ccdID,
                       size_t(640), size_t(480),
                       49.13_r, antiAliasFactor,
                       Vec3(-5,5,5), Vec3(-1,5,5), Vec3(0,0,1), //-5,5,5; -1,5,5
                       lighting,
                       Color(0.2_r,0.2_r,0.2_r),
                       customSpheresBodyToShadingParams);
   
   MaterialID iron = Material::find("iron");
   
   createPlane(*globalBodyStorage, 0, Vec3(0,-1,0), Vec3(0,10,0), iron); // left wall
   createPlane(*globalBodyStorage, 0, Vec3(0,1,0), Vec3(0,0,0), iron); // right wall
   createPlane(*globalBodyStorage, 0, Vec3(0,0,1), Vec3(0,0,0), iron); // floor
   createPlane(*globalBodyStorage, 0, Vec3(0,0,-1), Vec3(0,0,10), iron); // ceiling
   createPlane(*globalBodyStorage, 0, Vec3(-1,0,0), Vec3(10,0,0), iron); // back wall
   createPlane(*globalBodyStorage, 0, Vec3(1,0,0), Vec3(0,0,0), iron); // front wall, should not get rendered
   
   walberla::id_t id=0;
   for (int j=0; j<4; j++) {
      for (int i=0; i<4; i++) {
         createSphere(*globalBodyStorage, *forest, storageID, id, Vec3(6,real_c(i+1)*2_r,real_c(j+1)*2_r), 0.9_r);
         id++;
      }
   }
   
   raytracer.setImageOutputEnabled(true);
   
   raytracer.setRaytracingAlgorithm(raytracingAlgorithm);
   raytracer.generateImage<BodyTuple>(1);
}

ShadingParameters customHashGridsBodyToShadingParams(const BodyID body) {
   if (body->getTypeID() == Sphere::getStaticTypeID()) {
      return yellowShadingParams(body);
   }
   
   return defaultBodyTypeDependentShadingParams(body);
}


void HashGridsTest(Raytracer::Algorithm raytracingAlgorithm, walberla::uint8_t antiAliasFactor,
                   size_t boxes, size_t capsules, size_t spheres, size_t numberOfViews = 1,
                   real_t boxLenMin = 0.1_r, real_t boxLenMax = 0.2_r, bool boxRotation = false,
                   real_t capRadiusMin = 0.1_r, real_t capRadiusMax = 0.2_r, real_t capLenMin = 0.1_r, real_t capLenMax = 0.3_r,
                   real_t sphereRadiusMin = 0.1_r, real_t sphereRadiusMax = 0.3_r) {
   WALBERLA_LOG_INFO("Generating " << boxes << " boxes, " << capsules << " capsules and " << spheres << " spheres");
   
   using namespace walberla::pe::ccd;
   WcTimingTree tt;
   tt.start("Setup");
   
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   shared_ptr<BlockForest> forest = createBlockForest(AABB(0,0,0,4,4,4), Vector3<uint_t>(2,3,1), Vector3<bool>(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID = forest->addBlockData(ccd::createHashGridsDataHandling(globalBodyStorage, storageID), "CCD");
   
   const AABB& forestAABB = forest->getDomain();
   
   bool removeUnproblematic = false;
   std::vector<walberla::id_t> problematicBodyIDs = {165, 5, 31}; //{50, 44, 66, 155, 170, 51};
   std::vector<walberla::id_t> bodySIDs;
   
   // generate bodies for test
   std::vector<BodyID> bodies;
   for (size_t i = 0; i < boxes; i++) {
      real_t len = math::realRandom(boxLenMin, boxLenMax); //0.2 0.5
      real_t x = math::realRandom(forestAABB.xMin(), forestAABB.xMax());
      real_t y = math::realRandom(forestAABB.yMin(), forestAABB.yMax());
      real_t z = math::realRandom(forestAABB.zMin(), forestAABB.zMax());
      walberla::id_t id = walberla::id_t(i);
      BoxID box_ = createBox(*globalBodyStorage, *forest, storageID, id, Vec3(x, y, z), Vec3(len, len, len));
      WALBERLA_CHECK(box_ != nullptr);
      if (boxRotation) {
         box_->rotate(0, math::realRandom(0_r, 1_r)*math::M_PI, math::realRandom(0_r, 1_r)*math::M_PI);
      }
      bodies.push_back(box_);
      bodySIDs.push_back(box_->getSystemID());
   }
   
   for (size_t i = 0; i < capsules; i++) {
      real_t len = math::realRandom(capLenMin, capLenMax); // 0.2 0.5
      real_t radius = math::realRandom(capRadiusMin, capRadiusMax);
      real_t x = math::realRandom(forestAABB.xMin(), forestAABB.xMax());
      real_t y = math::realRandom(forestAABB.yMin(), forestAABB.yMax());
      real_t z = math::realRandom(forestAABB.zMin(), forestAABB.zMax());
      walberla::id_t id = walberla::id_t(boxes+i);
      CapsuleID capsule = createCapsule(*globalBodyStorage, *forest, storageID, id, Vec3(x, y, z), radius, len);
      WALBERLA_CHECK(capsule != nullptr);
      capsule->rotate(0, math::realRandom(0_r, 1_r)*math::M_PI, math::realRandom(0_r, 1_r)*math::M_PI);
      bodies.push_back(capsule);
      bodySIDs.push_back(capsule->getSystemID());
   }
   
   for (size_t i = 0; i < spheres; i++) {
      real_t radius = math::realRandom(sphereRadiusMin, sphereRadiusMax); // 0.2 0.3
      real_t x = math::realRandom(forestAABB.xMin(), forestAABB.xMax());
      real_t y = math::realRandom(forestAABB.yMin(), forestAABB.yMax());
      real_t z = math::realRandom(forestAABB.zMin(), forestAABB.zMax());
      walberla::id_t id = walberla::id_t(boxes+capsules+i);
      SphereID sphere = createSphere(*globalBodyStorage, *forest, storageID, id, Vec3(x, y, z), radius);
      WALBERLA_CHECK(sphere != nullptr);
      bodies.push_back(sphere);
      bodySIDs.push_back(sphere->getSystemID());
   }
   
   for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt) {
      ccd::HashGrids* hashgrids = blockIt->getData<ccd::HashGrids>(ccdID);
      hashgrids->update();
      for (auto bodyIt = LocalBodyIterator::begin(*blockIt, storageID); bodyIt != LocalBodyIterator::end(); ++bodyIt) {
         if (removeUnproblematic && std::find(problematicBodyIDs.begin(), problematicBodyIDs.end(), bodyIt->getID()) == problematicBodyIDs.end()) {
            bodyIt->setPosition(-100, -100, -100);
         }
      }
   }
   
   MaterialID iron = Material::find("iron");
   createPlane(*globalBodyStorage, 0, Vec3(0,-1,0), Vec3(0,forestAABB.yMax(),0), iron); // left wall
   createPlane(*globalBodyStorage, 0, Vec3(0,1,0), Vec3(0,forestAABB.yMin(),0), iron); // right wall
   createPlane(*globalBodyStorage, 0, Vec3(0,0,1), Vec3(0,0,forestAABB.zMin()), iron); // floor
   createPlane(*globalBodyStorage, 0, Vec3(0,0,-1), Vec3(0,0,forestAABB.zMax()), iron); // ceiling
   createPlane(*globalBodyStorage, 0, Vec3(-1,0,0), Vec3(forestAABB.xMax(),0,0), iron); // back wall
   createPlane(*globalBodyStorage, 0, Vec3(1,0,0), Vec3(forestAABB.xMin(),0,0), iron); // front wall, should not get rendered
   
   
   std::vector<std::tuple<Vec3, Vec3, Vec3>> viewVectors;
   
   // y up, in negative z direction
   viewVectors.emplace_back(Vec3(2, 2.1_r, 7),
                                     Vec3(2.1_r, 2, 4),
                                     Vec3(0,1,0));
   // y up, in positive z direction
   viewVectors.emplace_back(Vec3(2, 2, -3),
                                     Vec3(2, 2.1_r, 0.1_r),
                                     Vec3(0,1,0));
   // x up, in positive z direction
   viewVectors.emplace_back(Vec3(2, 2, -3),
                                     Vec3(2, 2.1_r, 0.1_r),
                                     Vec3(1,0,0));
   // y and x up, in positive z direction
   viewVectors.emplace_back(Vec3(2, 2, -3),
                                     Vec3(2, 2.1_r, 0.1_r),
                                     Vec3(1,1,0));
   // y and x up, in negative z direction
   viewVectors.emplace_back(Vec3(2, 2, 6.5),
                                     Vec3(2.1_r, 2.1_r, 4),
                                     Vec3(0.5_r,1,0));
   // z up, in positive x direction
   viewVectors.emplace_back(Vec3(-3, 2, 1.9_r),
                                     Vec3(0, 2.1_r, 2),
                                     Vec3(0,0,1));
   // z up, in negative x direction
   viewVectors.emplace_back(Vec3(7, 2, 1.9_r),
                                     Vec3(4, 2.1_r, 2),
                                     Vec3(0,0,1));
   // z and y up, in negative x direction
   viewVectors.emplace_back(Vec3(7, 2, 1.9_r),
                                     Vec3(4, 2.1_r, 2),
                                     Vec3(0,1,1));
   // z and x up, in negative y direction
   viewVectors.emplace_back(Vec3(2, 6, 1.9_r),
                                     Vec3(2.3_r, 4, 2),
                                     Vec3(1,0,1));
   // z up, in positive y direction
   viewVectors.emplace_back(Vec3(2, -3.6_r, 1.9_r),
                                     Vec3(2.3_r, 0, 2.1_r),
                                     Vec3(0,0,1));
   
   Lighting lighting0(Vec3(forestAABB.xSize()/2_r+1, forestAABB.ySize()/2_r,
                           2_r*forestAABB.zMax()+2), // 8, 5, 9.5 gut für ebenen, 0,5,8
                      Color(1, 1, 1), //diffuse
                      Color(1, 1, 1), //specular
                      Color(0.4_r, 0.4_r, 0.4_r)); //ambient
   tt.stop("Setup");

   size_t i = 0;
   for (auto& vector: viewVectors) {
      if (i == numberOfViews) {
         break;
      }
      
      Raytracer raytracer(forest, storageID, globalBodyStorage, ccdID,
                           size_t(640), size_t(480),
                           49.13_r, antiAliasFactor,
                           std::get<0>(vector),
                           std::get<1>(vector),
                           std::get<2>(vector),
                           lighting0,
                           Color(0.2_r,0.2_r,0.2_r),
                           customHashGridsBodyToShadingParams);
      raytracer.setImageOutputEnabled(true);
      raytracer.setFilenameTimestepWidth(12);
      WALBERLA_LOG_INFO("output #" << i << " to: " << (boxes*100000000 + capsules*10000 + spheres) << " in " << raytracer.getImageOutputDirectory());
      raytracer.setRaytracingAlgorithm(raytracingAlgorithm);
      raytracer.generateImage<BodyTuple>(boxes*100000000 + capsules*10000 + spheres, &tt);
      i++;
   }
   
   auto temp = tt.getReduced();
   WALBERLA_ROOT_SECTION() {
      std::cout << temp;
   }
}

ShadingParameters customArtifactsBodyToShadingParams(const BodyID body) {
   if (body->getTypeID() == Box::getStaticTypeID()) {
      return lightGreyShadingParams(body);
   }
   return defaultShadingParams(body);
}

void raytraceArtifactsForest(Raytracer::Algorithm raytracingAlgorithm, walberla::uint8_t antiAliasFactor,
                             const shared_ptr<BlockStorage> forest, const BlockDataID storageID,
                             const shared_ptr<BodyStorage> globalBodyStorage,
                             const BlockDataID ccdID,
                             const Vec3& cameraPosition, const Vec3& lookAtPoint, const Vec3& upVector,
                             size_t numberOfBoxes, walberla::uint8_t timestepWidth) {
   WcTimingTree tt;

   Lighting lighting(cameraPosition,
                     Color(1, 1, 1), //diffuse
                     Color(1, 1, 1), //specular
                     Color(0.4_r, 0.4_r, 0.4_r)); //ambient
   
   Raytracer raytracer(forest, storageID, globalBodyStorage, ccdID,
                       size_t(640), size_t(480),
                       49.13_r, antiAliasFactor,
                       cameraPosition,
                       lookAtPoint,
                       upVector,
                       lighting,
                       Color(0.2_r,0.2_r,0.2_r),
                       customArtifactsBodyToShadingParams);
   raytracer.setImageOutputEnabled(true);
   raytracer.setFilenameTimestepWidth(timestepWidth);
   raytracer.setRaytracingAlgorithm(raytracingAlgorithm);
   WALBERLA_LOG_INFO("output to: " << numberOfBoxes << " in " << raytracer.getImageOutputDirectory());
   raytracer.generateImage<BodyTuple>(numberOfBoxes, &tt);
   
   auto temp = tt.getReduced();
   WALBERLA_ROOT_SECTION() {
      std::cout << temp;
   }
}

void HashGridsArtifactsTest(Raytracer::Algorithm raytracingAlgorithm, walberla::uint8_t antiAliasFactor,
                            size_t boxes, real_t boxLenMin = 0.1_r, real_t boxLenMax = 0.2_r) {
   WALBERLA_LOG_INFO_ON_ROOT("HashGrids Artifacts Test - In negative Z direction");
   
   WALBERLA_LOG_INFO(" Generating " << boxes << " boxes");
   
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   shared_ptr<BlockForest> forest = createBlockForest(AABB(0,0,0,4,4,4), Vector3<uint_t>(1,1,1), Vector3<bool>(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID = forest->addBlockData(ccd::createHashGridsDataHandling(globalBodyStorage, storageID), "CCD");
   
   const AABB& forestAABB = forest->getDomain();
   
   // generate bodies for test
   for (size_t i = 0; i < boxes; i++) {
      real_t len = math::realRandom(boxLenMin, boxLenMax); //0.2 0.5
      real_t x_min = math::realRandom(forestAABB.xMin()+len/2_r, forestAABB.xMax());
      real_t y_min = math::realRandom(forestAABB.yMin()+len/2_r, forestAABB.yMax());
      real_t z_min = math::realRandom(forestAABB.zMin()+len/2_r, forestAABB.zMax());
      if (i%5 == 0) {
         x_min = forestAABB.xMax() - math::realRandom(len/2_r, len);
      } else if (i%5 == 1){
         x_min = forestAABB.xMin() + math::realRandom(0_r, len/2_r);
      } else if (i%5 == 2){
         y_min = forestAABB.yMax() - math::realRandom(len/2_r, len);
      } else if (i%5 == 3){
         y_min = forestAABB.yMin() + math::realRandom(0_r, len/2_r);
      } else if (i%5 == 4){
         z_min = forestAABB.zMin() + math::realRandom(0_r, len/2_r);
      }
      walberla::id_t id = walberla::id_t(i);
      BoxID box_ = createBox(*globalBodyStorage, *forest, storageID, id, Vec3(x_min, y_min, z_min), Vec3(len, len, len));
      WALBERLA_CHECK(box_ != nullptr);
   }
   
   raytraceArtifactsForest(raytracingAlgorithm, antiAliasFactor,
                           forest, storageID, globalBodyStorage, ccdID,
                           Vec3(2, 2, 7), Vec3(2, 2, 4), Vec3(0,1,0),
                           boxes, 3);
}

void HashGridsFromNegativeArtifactsTest(Raytracer::Algorithm raytracingAlgorithm, walberla::uint8_t antiAliasFactor,
                                        size_t boxes, real_t boxLenMin = 0.1_r, real_t boxLenMax = 0.2_r) {
   WALBERLA_LOG_INFO_ON_ROOT("HashGrids Artifacts Test - In positive Z direction");
   
   WALBERLA_LOG_INFO_ON_ROOT(" Generating " << boxes << " boxes");
   
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   shared_ptr<BlockForest> forest = createBlockForest(AABB(0,0,0,4,4,4), Vector3<uint_t>(1,1,1), Vector3<bool>(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID = forest->addBlockData(ccd::createHashGridsDataHandling(globalBodyStorage, storageID), "CCD");
   
   const AABB& forestAABB = forest->getDomain();

   // generate bodies for test
   std::vector<BodyID> bodies;
   for (size_t i = 0; i < boxes; i++) {
      real_t len = math::realRandom(boxLenMin, boxLenMax); //0.2 0.5
      real_t x_min = math::realRandom(forestAABB.xMin()+len/2_r, forestAABB.xMax());
      real_t y_min = math::realRandom(forestAABB.yMin()+len/2_r, forestAABB.yMax());
      real_t z_min = math::realRandom(forestAABB.zMin()+len/2_r, forestAABB.zMax());
      
      if (i%5 == 0) {
         x_min = forestAABB.xMax() - math::realRandom(len/2_r, len);
      } else if (i%5 == 1){
         x_min = forestAABB.xMin() + math::realRandom(0_r, len/2_r);
      } else if (i%5 == 2){
         y_min = forestAABB.yMax() - math::realRandom(len/2_r, len);
      } else if (i%5 == 3){
         y_min = forestAABB.yMin() + math::realRandom(0_r, len/2_r);
      } else if (i%5 == 4){
         z_min = forestAABB.zMax() - math::realRandom(len/2_r, len);
      }
      
      //real_t z_min = len+0.1;
      walberla::id_t id = walberla::id_t(i);
      BoxID box_ = createBox(*globalBodyStorage, *forest, storageID, id, Vec3(x_min, y_min, z_min), Vec3(len, len, len));
      WALBERLA_CHECK(box_ != nullptr);
   }
   
   raytraceArtifactsForest(raytracingAlgorithm, antiAliasFactor,
                           forest, storageID, globalBodyStorage, ccdID,
                           Vec3(2, 2, -3), Vec3(2, 2, 0), Vec3(0,1,0),
                           boxes, 4);
}

void HashGridsFromNegativeXArtifactsTest(Raytracer::Algorithm raytracingAlgorithm, walberla::uint8_t antiAliasFactor,
                                         size_t boxes, real_t boxLenMin = 0.1_r, real_t boxLenMax = 0.2_r) {
   WALBERLA_LOG_INFO_ON_ROOT("HashGrids Artifacts Test - In positive X direction");
   WALBERLA_LOG_INFO_ON_ROOT(" Generating " << boxes << " boxes");
   
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   shared_ptr<BlockForest> forest = createBlockForest(AABB(0,0,0,4,4,4), Vector3<uint_t>(1,1,1), Vector3<bool>(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID = forest->addBlockData(ccd::createHashGridsDataHandling(globalBodyStorage, storageID), "CCD");
   
   const AABB& forestAABB = forest->getDomain();
   
   // generate bodies for test
   for (size_t i = 0; i < boxes; i++) {
      real_t len = math::realRandom(boxLenMin, boxLenMax); //0.2 0.5
      real_t x_min = math::realRandom(forestAABB.xMin()+len/2_r, forestAABB.xMax());
      real_t y_min = math::realRandom(forestAABB.yMin()+len/2_r, forestAABB.yMax());
      real_t z_min = math::realRandom(forestAABB.zMin()+len/2_r, forestAABB.zMax());
      
      if (i%5 == 0) {
         z_min = forestAABB.zMax() - math::realRandom(len/2_r, len);
      } else if (i%5 == 1){
         z_min = forestAABB.zMin() + math::realRandom(0_r, len/2_r);
      } else if (i%5 == 2){
         y_min = forestAABB.yMax() - math::realRandom(len/2_r, len);
      } else if (i%5 == 3){
         y_min = forestAABB.yMin() + math::realRandom(0_r, len/2_r);
      } else if (i%5 == 4){
         x_min = forestAABB.xMax() - math::realRandom(len/2_r, len);
      }
      
      //real_t z_min = len+0.1;
      walberla::id_t id = walberla::id_t(i);
      BoxID box_ = createBox(*globalBodyStorage, *forest, storageID, id, Vec3(x_min, y_min, z_min), Vec3(len, len, len));
      WALBERLA_CHECK(box_ != nullptr);
   }
   
   raytraceArtifactsForest(raytracingAlgorithm, antiAliasFactor,
                           forest, storageID, globalBodyStorage, ccdID,
                           Vec3(-3, 2, 2), Vec3(0, 2, 2), Vec3(0,0,1),
                           boxes, 6);
}


Vec3 minCornerToGpos(const Vec3& minCorner, real_t lengths) {
   return minCorner + Vec3(lengths/2, lengths/2, lengths/2);
}

void HashGridsTestScene(Raytracer::Algorithm raytracingAlgorithm = Raytracer::RAYTRACE_HASHGRIDS, walberla::uint8_t antiAliasFactor = 1) {
   WALBERLA_LOG_INFO_ON_ROOT("HashGrids Test Scene");
   
   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();
   shared_ptr<BlockForest> forest = createBlockForest(AABB(0,0,0,8,8,8), Vector3<uint_t>(1,1,1), Vector3<bool>(false, false, false));
   auto storageID = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto ccdID = forest->addBlockData(ccd::createHashGridsDataHandling(globalBodyStorage, storageID), "CCD");
   
   const AABB& forestAABB = forest->getDomain();
   
   std::vector<BodyID> bodies;
   
   // create bodies
   size_t id = 0;
   real_t len = 0.6_r;
   
   real_t x_min = 0, y_min = 0;
   len = 1.2_r;
   real_t gap = 0.4_r;
   
   // cubes on z = 0 plane
   for (int i = 0; ; ++i) {
      x_min = forestAABB.xMin() + real_c(i)*(gap+len);
      if (x_min > forestAABB.max(0)) {
         break;
      }
      for (int j = 0; ; ++j) {
         y_min = forestAABB.yMin() + real_c(j)*(gap+len);
         if (y_min > forestAABB.max(1)) {
            break;
         }
         
         bodies.push_back(createBox(*globalBodyStorage, *forest, storageID, ++id, minCornerToGpos(Vec3(x_min, y_min, 0), len), Vec3(len, len, len)));
      }
   }
   
   //cubes on z = max plane
   for (int i = 0; ; ++i) {
      x_min = forestAABB.xMin() + real_c(i)*(gap+len);
      if (x_min > forestAABB.max(0)) {
         break;
      }
      for (int j = 0; ; ++j) {
         y_min = forestAABB.yMin() + real_c(j)*(gap+len);
         if (y_min > forestAABB.max(1)) {
            break;
         }
         
         bodies.push_back(createBox(*globalBodyStorage, *forest, storageID, ++id, minCornerToGpos(Vec3(x_min, y_min, forestAABB.zMax()-len), len), Vec3(len, len, len)));
      }
   }
   
   std::vector<std::tuple<Vec3, Vec3, Vec3>> viewVectors;
   
   // in negative x direction -> cubes to the right
   viewVectors.emplace_back(Vec3(15,4,4),
                                         Vec3(8,4,4),
                                         Vec3(0,1,0));
   // in negative x direction and negative z direction, up vector in y direction -> cubes from the right tilted
   viewVectors.emplace_back(Vec3(12,4,8),
                                         Vec3(6,4,2),
                                         Vec3(0,1,0));
   // in negative x direction and negative z direction, up vector in negative y direction
   viewVectors.emplace_back(Vec3(12,4,8),
                                         Vec3(6,4,2),
                                         Vec3(0,-1,0));
   // in positive x direction
   viewVectors.emplace_back(Vec3(-7,4,4),
                                         Vec3(0,4,4),
                                         Vec3(0,1,0));
   // in negative x direction
   viewVectors.emplace_back(Vec3(4,4,15),
                                         Vec3(4,4,8),
                                         Vec3(0,1,0));
   
   WcTimingTree tt;
   
   Lighting lighting(Vec3(1,2,15),
                     Color(1, 1, 1), //diffuse
                     Color(1, 1, 1), //specular
                     Color(0.4_r, 0.4_r, 0.4_r)); //ambient
   
   int i = 0;
   for (auto& vector: viewVectors) {
      Raytracer raytracer(forest, storageID, globalBodyStorage, ccdID,
                          size_t(640), size_t(480),
                          49.13_r, antiAliasFactor,
                          std::get<0>(vector),
                          std::get<1>(vector),
                          std::get<2>(vector),
                          lighting,
                          Color(0.2_r,0.2_r,0.2_r));
      
      raytracer.setRaytracingAlgorithm(raytracingAlgorithm);
      raytracer.setImageOutputEnabled(true);
      raytracer.setFilenameTimestepWidth(1);
      WALBERLA_LOG_INFO("output to: " << i << " in " << raytracer.getImageOutputDirectory());
      raytracer.setRaytracingAlgorithm(raytracingAlgorithm);
      raytracer.generateImage<BodyTuple>(size_t(i), &tt);
      
      auto temp = tt.getReduced();
      WALBERLA_ROOT_SECTION() {
         std::cout << temp;
      }
      
      i++;
   }
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   SetBodyTypeIDs<BodyTuple>::execute();
   math::seedRandomGenerator( static_cast<unsigned int>(1337 * mpi::MPIManager::instance()->worldRank()) );
   
   SphereIntersectsTest();
   PlaneIntersectsTest();
   BoxIntersectsTest();
   AABBIntersectsTest();
   CapsuleIntersectsTest();
   EllipsoidTest();

   const Raytracer::Algorithm algorithm = Raytracer::RAYTRACE_COMPARE_BOTH_STRICTLY;
   const walberla::uint8_t antiAliasFactor = 1;
   RaytracerTest(algorithm, antiAliasFactor);
   RaytracerSpheresTestScene(algorithm, antiAliasFactor);
   HashGridsTestScene(algorithm, antiAliasFactor);
   
   if (argc >= 2 && strcmp(argv[1], "--longrun") == 0) {
      HashGridsTest(algorithm, antiAliasFactor,
                    50, 30, 130,
                    10);
      HashGridsTest(algorithm, antiAliasFactor,
                    60, 60, 3,
                    1,
                    0.1_r, 0.3_r, true,
                    0.1_r, 0.2_r, 0.1_r, 0.2_r,
                    0.5_r, 0.6_r);
      HashGridsArtifactsTest(algorithm, antiAliasFactor, 750, 0.2_r, 0.3_r);
      HashGridsFromNegativeArtifactsTest(algorithm, antiAliasFactor, 750, 0.2_r, 0.3_r);
      HashGridsFromNegativeXArtifactsTest(algorithm, antiAliasFactor, 750, 0.2_r, 0.3_r);
   }
   
   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
