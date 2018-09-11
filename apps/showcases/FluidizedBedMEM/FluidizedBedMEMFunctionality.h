#pragma once

namespace FBfunc {

    using namespace walberla;
    using walberla::uint_t;

    typedef walberla::uint8_t flag_t;

////////////////
// Parameters //
////////////////

    struct SetupFB {

        // geometry
        uint_t xlength;                //[m]
        uint_t ylength;                //[m]
        uint_t zlength;               //[m]
        real_t offsetBot = real_c(0);               //[m]
        real_t offsetTop = real_c(0);               //[m]

        real_t diameterA;     //[m]
        real_t diameterB;     //[m]

        // fluidized bed properties
        real_t kinViscosity;         //[m2/s]

        real_t densityRatioA;
        real_t densityRatioB;

        real_t velocity;             //[m/s]
        real_t gravity;              //[m/s2]

        uint_t nrParticlesA;
        uint_t nrParticlesB;

        //conversion
        real_t dt_SI;                   //[s]
        real_t dx_SI;                   //[m]
        real_t densityFluid_SI;         //[kg/m3]

        // simulation parameters
        real_t extForce;

        uint_t evaluationFreq;
        uint_t velCheckFreq;
        real_t stopVel;
        uint_t peSubCycles;
        real_t lbmSubCycles;

        uint_t checkpointFreq;

        // post processing

        std::vector<CellInterval> probesCell;
        std::vector<AABB> probesAABB;

    };

    class VelocityFunctor_T{

    public:

        VelocityFunctor_T(const real_t xCenter, const real_t zCenter, const real_t spotDiameter , const real_t  velMax, const real_t initialRamping, const real_t rampingTime, const bool spotInflow):
                xCenter_(xCenter), zCenter_(zCenter), spotRadius_(spotDiameter * 0.5), velMax_(velMax)
                ,initialRamping_(initialRamping), rampingTime_(rampingTime), spotInflow_(spotInflow)
        {}

        void operator()( const  real_t t  )
        {
            (void)t;
        }

        Vector3< real_t > operator()( const Vector3< real_t > & x, const real_t t ){

            Vector3<real_t> velocity;

            if (spotInflow_){
                velocity = Vector3<real_t>(0.0);
                real_t distance = sqrt((xCenter_ - x[0])*(xCenter_ - x[0]) + (zCenter_ - x[2])*(zCenter_ - x[2]));

                if (distance < spotRadius_){

                    velocity[1] = (1.0 - ( distance * distance )/( spotRadius_ * spotRadius_ ) ) * velMax_;

                    if (t < rampingTime_){
                        velocity[1] *= (initialRamping_ + (t / rampingTime_)*(1.0 - initialRamping_));
                    }
                }

            } else {

                velocity[1] = velMax_;

                if (t < rampingTime_){
                    velocity[1] *= (initialRamping_ + (t / rampingTime_)*(1.0 - initialRamping_));
                }

            }

            return velocity;
        }


    private:

        const real_t xCenter_;
        const real_t zCenter_;
        const real_t spotRadius_;
        const real_t velMax_;
        const real_t initialRamping_;
        const real_t rampingTime_;
        const bool spotInflow_;
    };


///////////
// FLAGS //
///////////

    const FlagUID Fluid_Flag("fluid");
    const FlagUID UBB_Flag("velocity bounce back");
    const FlagUID NoSlip_Flag("no slip");
    const FlagUID Outlet_Flag("outlet");
    const FlagUID MO_Flag("moving obstacle");
    const FlagUID FormerMO_Flag("former moving obstacle");


// ObstacleLocationCheck
    class ObstacleLocationCheck {
    public:

        ObstacleLocationCheck(const shared_ptr<StructuredBlockStorage> &blocks,
                              const BlockDataID &bodyStorageID,
                              const AABB &aabb, const real_t additionalSpace)
                : blocks_(blocks), bodyStorageID_(bodyStorageID),
                  aabb_(aabb.xMin() - additionalSpace, aabb.yMin() - additionalSpace, aabb.zMin() - additionalSpace,
                        aabb.xMax() + additionalSpace, aabb.yMax() + additionalSpace, aabb.zMax() + additionalSpace) {}

        void operator()() {
            for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt) {
                for (auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt) {
                    if (bodyIt->isFixed() || !bodyIt->isFinite())
                        continue;

                    WALBERLA_CHECK(aabb_.contains(bodyIt->getPosition()[0], bodyIt->getPosition()[1], bodyIt->getPosition()[2]));
                }
            }
        }

    private:
        shared_ptr<StructuredBlockStorage> blocks_;
        const BlockDataID bodyStorageID_;
        const AABB aabb_;

    };

// CALCULATE PRESSURE DROP
    class PressureDropper {
    public:

        PressureDropper(const shared_ptr<StructuredBlockStorage> &blocks,
                        const BlockDataID &bodyStorageID,
                        SetupFB &config, const bool preliminaryTimeloop ) :
                blocks_(blocks), bodyStorageID_(bodyStorageID),
                checkFrequency_((config.evaluationFreq > 0) ? config.evaluationFreq : uint_c(1)), executionCounter_(uint_c(0)),
                preFac_((config.densityFluid_SI * config.dx_SI * config.dx_SI / (config.dt_SI * config.dt_SI)) / (real_c(config.xlength * config.zlength))),
                preliminaryTimeloop_(preliminaryTimeloop), currentPressureDrop_(real_t(0)){

            std::ofstream output;
            output.open("PressureDrop.txt", std::ofstream::out | std::ofstream::trunc);
            output.close();
        }

        void operator()() {
            ++executionCounter_;
            if (((executionCounter_ - uint_c(1)) % checkFrequency_ != 0) && !preliminaryTimeloop_)
                return;

            pe::Vec3 force = pe::Vec3(0.0, 0.0, 0.0);

            for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt) {
                for (auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt) {
                    force += bodyIt->getForce();
                }
            }

            WALBERLA_MPI_SECTION() {
                mpi::reduceInplace(force[0], mpi::SUM);
                mpi::reduceInplace(force[1], mpi::SUM);
                mpi::reduceInplace(force[2], mpi::SUM);
            }

            WALBERLA_ROOT_SECTION() {

                currentPressureDrop_ = force.length() * preFac_;

                if (!preliminaryTimeloop_) {
                    std::ofstream output;
                    output.open("PressureDrop.txt", std::ofstream::out | std::ofstream::app);
                    output << executionCounter_ << "\t" << currentPressureDrop_ << "\t[Pa]" << std::endl;
                    output.close();
                }
            }

        }

        real_t getPressureDrop () const
        {
            return currentPressureDrop_;
        }

    private:
        shared_ptr<StructuredBlockStorage> blocks_;
        const BlockDataID bodyStorageID_;
        const uint_t checkFrequency_;
        uint_t executionCounter_;
        const real_t preFac_;
        const bool preliminaryTimeloop_;

        real_t currentPressureDrop_;

    };

// ADD EXTERNAL FORCES
    class ExternalForce {

    public:
        ExternalForce(const shared_ptr<StructuredBlockStorage> &blocks,
                       const BlockDataID &bodyStorageID,
                       SetupFB &config) :
                blocks_(blocks), bodyStorageID_(bodyStorageID),
                preFac_(config.gravity * (config.densityRatioA - real_c(1.0)) * (real_c(4.0) / real_c(3.0)) * math::PI) {
        }

        void operator()() {

            real_t radius;

            for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt) {
                for (auto bodyIt = walberla::pe::LocalBodyIterator::begin<walberla::pe::Sphere>(*blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end<walberla::pe::Sphere>(); ++bodyIt) {
                    radius = bodyIt->getRadius();
                    bodyIt->addForce(0.0, -preFac_ * (radius * radius * radius), 0.0);
                }
            }
        }

    private:
        shared_ptr<StructuredBlockStorage> blocks_;
        const BlockDataID bodyStorageID_;
        const real_t preFac_;
    };

// CALCULATE Granular Temperature
    class GranularTemperature {

    public:
        GranularTemperature(const shared_ptr<StructuredBlockStorage> &blocks,
                            const BlockDataID &bodyStorageID,
                            SetupFB &config) :
                blocks_(blocks), bodyStorageID_(bodyStorageID),
                nrParticles_(real_c(config.nrParticlesA)),
                preFac_(config.dx_SI * config.dx_SI / (config.dt_SI * config.dt_SI)),
                checkFrequency_((config.evaluationFreq > 0) ? config.evaluationFreq : uint_c(1)), executionCounter_(0)
        {
            std::ofstream output;
            output.open("GranularTemperatureSQRT.txt", std::ofstream::out | std::ofstream::trunc);
            output.close();
        }

        void operator()() {
            ++executionCounter_;
            if ((executionCounter_ - uint_c(1)) % checkFrequency_ != 0)
                return;

            real_t velSumX = real_c(0);
            real_t velSumY = real_c(0);
            real_t velSumZ = real_c(0);

            for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt) {
                for (auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt) {
                    velSumX += bodyIt->getLinearVel()[0];
                    velSumY += bodyIt->getLinearVel()[1];
                    velSumZ += bodyIt->getLinearVel()[2];
                }
            }

            WALBERLA_MPI_SECTION() {
                mpi::allReduceInplace(velSumX, mpi::SUM);
                mpi::allReduceInplace(velSumY, mpi::SUM);
                mpi::allReduceInplace(velSumZ, mpi::SUM);
            }

            real_t velMeanX = velSumX / nrParticles_;
            real_t velMeanY = velSumY / nrParticles_;
            real_t velMeanZ = velSumZ / nrParticles_;

            real_t granTempX = real_c(0);
            real_t granTempY = real_c(0);
            real_t granTempZ = real_c(0);

            for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt) {
                for (auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt) {

                    granTempX += (bodyIt->getLinearVel()[0] - velMeanX) * (bodyIt->getLinearVel()[0] - velMeanX);
                    granTempY += (bodyIt->getLinearVel()[1] - velMeanY) * (bodyIt->getLinearVel()[1] - velMeanY);
                    granTempZ += (bodyIt->getLinearVel()[2] - velMeanZ) * (bodyIt->getLinearVel()[2] - velMeanZ);
                }
            }

            WALBERLA_MPI_SECTION() {
                mpi::reduceInplace(granTempX, mpi::SUM);
                mpi::reduceInplace(granTempY, mpi::SUM);
                mpi::reduceInplace(granTempZ, mpi::SUM);
            }

            WALBERLA_ROOT_SECTION() {
                granTempX = granTempX * preFac_;
                granTempY = granTempY * preFac_;
                granTempZ = granTempZ * preFac_;
                std::ofstream output;
                output.open("GranularTemperatureSQRT.txt", std::ofstream::out | std::ofstream::app);
                output << executionCounter_ << "\t" << sqrt(granTempX) << "\t" << sqrt(granTempY) << "\t" << sqrt(granTempZ) << "\t[m/s]" << std::endl;
                output.close();
            }

        }

    private:
        shared_ptr<StructuredBlockStorage> blocks_;
        const BlockDataID bodyStorageID_;
        const real_t nrParticles_;
        const real_t preFac_;
        const uint_t checkFrequency_;
        uint_t executionCounter_;

    };

    // CALCULATE Granular Temperature at local volume
    class GranularTemperatureLocally {

    public:
        GranularTemperatureLocally(const shared_ptr<StructuredBlockStorage> &blocks,
                                   const BlockDataID &bodyStorageID,
                                   SetupFB &config) :
                blocks_(blocks), bodyStorageID_(bodyStorageID),
                nrParticles_(real_c(config.nrParticlesA)),
                preFac_(config.dx_SI * config.dx_SI / (config.dt_SI * config.dt_SI)),
                probes_(config.probesAABB),
                checkFrequency_((config.evaluationFreq > 0) ? config.evaluationFreq : uint_c(1)), executionCounter_(0)
        {
            std::ofstream output;
            output.open("GranularTemperatureLocalSQRT.txt", std::ofstream::out | std::ofstream::trunc);
            output.close();
        }

        void operator()() {
            ++executionCounter_;
            if ((executionCounter_ - uint_c(1)) % checkFrequency_ != 0)
                return;

            uint_t probeCounter = uint_c(1);

            real_t velSumX;
            real_t velSumY;
            real_t velSumZ;
            real_t velMeanX;
            real_t velMeanY;
            real_t velMeanZ;
            real_t granTempX;
            real_t granTempY;
            real_t granTempZ;

            for (std::vector<AABB>::const_iterator it = probes_.begin(); it != probes_.end(); ++it) {

                velSumX = real_c(0);
                velSumY = real_c(0);
                velSumZ = real_c(0);

                for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt) {

                    for (auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt) {
                        if ( it->contains(bodyIt->getPosition()) ) {
                            velSumX += bodyIt->getLinearVel()[0];
                            velSumY += bodyIt->getLinearVel()[1];
                            velSumZ += bodyIt->getLinearVel()[2];
                        }
                    }
                }

                WALBERLA_MPI_SECTION() {
                    mpi::allReduceInplace(velSumX, mpi::SUM);
                    mpi::allReduceInplace(velSumY, mpi::SUM);
                    mpi::allReduceInplace(velSumZ, mpi::SUM);
                }

                velMeanX = velSumX / nrParticles_;
                velMeanY = velSumY / nrParticles_;
                velMeanZ = velSumZ / nrParticles_;

                granTempX = real_c(0);
                granTempY = real_c(0);
                granTempZ = real_c(0);

                for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt) {
                    for (auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt) {
                        if ( it->contains(bodyIt->getPosition()) ) {

                            granTempX += (bodyIt->getLinearVel()[0] - velMeanX) * (bodyIt->getLinearVel()[0] - velMeanX);
                            granTempY += (bodyIt->getLinearVel()[1] - velMeanY) * (bodyIt->getLinearVel()[1] - velMeanY);
                            granTempZ += (bodyIt->getLinearVel()[2] - velMeanZ) * (bodyIt->getLinearVel()[2] - velMeanZ);
                        }
                    }
                }

                WALBERLA_MPI_SECTION() {
                    mpi::reduceInplace(granTempX, mpi::SUM);
                    mpi::reduceInplace(granTempY, mpi::SUM);
                    mpi::reduceInplace(granTempZ, mpi::SUM);
                }

                WALBERLA_ROOT_SECTION() {
                    granTempX = granTempX * preFac_;
                    granTempY = granTempY * preFac_;
                    granTempZ = granTempZ * preFac_;
                    std::ofstream output;
                    output.open("GranularTemperatureLocalSQRT" + std::to_string(probeCounter) + ".txt", std::ofstream::out | std::ofstream::app);
                    output << executionCounter_ << "\t" << sqrt(granTempX) << "\t" << sqrt(granTempY) << "\t" << sqrt(granTempZ) << "\t[m/s]" << std::endl;

                    output.close();
                }

                ++probeCounter;

            }
        }

    private:
        shared_ptr<StructuredBlockStorage> blocks_;
        const BlockDataID bodyStorageID_;
        const real_t nrParticles_;
        const real_t preFac_;
        const std::vector<AABB> probes_;
        const uint_t checkFrequency_;
        uint_t executionCounter_;

    };

// CALCULATE CENTER OF MASS
    class CenterOfMass {

    public:
        CenterOfMass(const shared_ptr<StructuredBlockStorage> &blocks,
                     const BlockDataID &bodyStorageID,
                     SetupFB &config, const bool calculateErgun) :
                blocks_(blocks), bodyStorageID_(bodyStorageID),
                dx_SI_(config.dx_SI), offsetBot_(real_c(config.offsetBot)), diamater_(config.diameterA),
                totalParticleVolume_(real_c(config.nrParticlesA) * 1.0 / 6.0 * math::PI * config.diameterA * config.diameterA * config.diameterA),
                crossSectionArea_( real_c( config.xlength * config.zlength ) ),
                dynViscosity_(config.kinViscosity),
                densityFluid_(real_c(1.0)),
                inflowVel_(config.velocity),
                pascalUnitConversion_(config.densityFluid_SI * config.dx_SI * config.dx_SI / ( config.dt_SI * config.dt_SI )),
                checkFrequency_((config.evaluationFreq > 0) ? config.evaluationFreq : uint_c(1)), executionCounter_(uint_c(0)),
                calculateErgun_(calculateErgun) {

            std::ofstream output1;
            output1.open("CenterOfMass.txt", std::ofstream::out | std::ofstream::trunc);
            output1.close();
            if (calculateErgun) {
                std::ofstream output2;
                output2.open("ErgunPressure.txt", std::ofstream::out | std::ofstream::trunc);
                output2.close();
            }

        }

        void operator()() {

            ++executionCounter_;
            if ((executionCounter_ - uint_c(1)) % checkFrequency_ != 0)
                return;

            real_t posSum = real_c(0);
            real_t massSum = real_c(0);
            real_t massCenter = real_c(0);
            real_t voidFraction;
            real_t pressureDropErgun;

            for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt) {
                for (auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt) {
                    massSum += bodyIt->getMass();
                    posSum += bodyIt->getPosition()[1] * bodyIt->getMass();
                }
            }

            WALBERLA_MPI_SECTION() {
                mpi::reduceInplace(posSum, mpi::SUM);
                mpi::reduceInplace(massSum, mpi::SUM);
            }

            WALBERLA_ROOT_SECTION() {

                massCenter = (posSum / massSum - offsetBot_);

                std::ofstream output1;
                output1.open("CenterOfMass.txt", std::ofstream::out | std::ofstream::app);
                output1 << executionCounter_ << "\t" << massCenter * dx_SI_ * 1000 << "\t[mm]" << std::endl;
                output1.close();

                if (calculateErgun_) {
                    std::ofstream output2;
                    voidFraction = 1.0 - totalParticleVolume_ / (crossSectionArea_ * massCenter * 2.0);
                    pressureDropErgun = ((150 * (1 - voidFraction) * (1 - voidFraction) * dynViscosity_ * inflowVel_) /
                                         (voidFraction * voidFraction * voidFraction * diamater_ * diamater_)
                                         + (1.75 * (1 - voidFraction) * densityFluid_ * inflowVel_ * inflowVel_) / (voidFraction * voidFraction * diamater_))
                                        * (massCenter * 2.0);

                    output2.open("ErgunPressure.txt", std::ofstream::out | std::ofstream::app);
                    output2 << executionCounter_ << "\t" << "VoidFraction: " << voidFraction << "\t\tErgunPressureDrop: " << pressureDropErgun * pascalUnitConversion_ << "\t[Pa]" << std::endl;

                    output2.close();
                }
            }
        }


    private:
        shared_ptr<StructuredBlockStorage> blocks_;
        const BlockDataID bodyStorageID_;

        const real_t dx_SI_;
        const real_t offsetBot_;
        const real_t diamater_;
        const real_t totalParticleVolume_;
        const real_t crossSectionArea_;
        const real_t dynViscosity_;
        const real_t densityFluid_;
        const real_t inflowVel_;
        const real_t pascalUnitConversion_;
        const uint_t checkFrequency_;
        uint_t executionCounter_;
        const bool calculateErgun_;
    };

// CALCULATE PARTICLE FLUX
    class ParticleFlux {

    public:

        ParticleFlux(const shared_ptr<StructuredBlockStorage> &blocks,
                       const BlockDataID &bodyStorageID,
                       SetupFB &config) :
                blocks_(blocks), bodyStorageID_(bodyStorageID),
                checkFrequency_((config.evaluationFreq > 0) ? config.evaluationFreq : uint_c(1)), executionCounter_(uint_c(0)),
                preFac_((real_c(4.0) / real_c(3.0) * math::PI * config.densityRatioA * config.densityFluid_SI * config.dx_SI / config.dt_SI) /
                        (real_c(config.xlength * config.zlength))) {
            std::ofstream output;
            output.open("MassFluxY.txt", std::ofstream::out | std::ofstream::trunc);
            output.close();
            output.open("MassFluxAbs.txt", std::ofstream::out | std::ofstream::trunc);
            output.close();
        }

        void operator()() {
            ++executionCounter_;
            if ((executionCounter_ - uint_c(1)) % checkFrequency_ != 0)
                return;

            real_t velSumY = real_c(0);
            real_t velSumAbs = real_c(0);

            for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt) {
                for (auto bodyIt = walberla::pe::LocalBodyIterator::begin<walberla::pe::Sphere>(*blockIt, bodyStorageID_); bodyIt != walberla::pe::LocalBodyIterator::end<walberla::pe::Sphere>(); ++bodyIt) {

                    velSumY += bodyIt->getLinearVel()[1] * bodyIt->getRadius() * bodyIt->getRadius();
                    velSumAbs += fabs(bodyIt->getLinearVel()[1]) * bodyIt->getRadius() * bodyIt->getRadius();
                }
            }

            // Equation inspired from:
            // Interface-resolved direct numerical simulation of the erosion of a sediment bed sheared by laminar channel flow
            // Kidanemariam and Uhlmann, 2014

            WALBERLA_MPI_SECTION() {
                mpi::reduceInplace(velSumY, mpi::SUM);
                mpi::reduceInplace(velSumAbs, mpi::SUM);
            }
            WALBERLA_ROOT_SECTION() {
                velSumY = velSumY * preFac_;
                velSumAbs = velSumAbs * preFac_;

                std::ofstream output;
                output.open("MassFluxY.txt", std::ofstream::out | std::ofstream::app);
                output << executionCounter_ << "\t" << velSumY << "\t[kg/sm2]" << std::endl;
                output.close();
                output.open("MassFluxAbs.txt", std::ofstream::out | std::ofstream::app);
                output << executionCounter_ << "\t" << velSumAbs << "\t[kg/sm2]" << std::endl;
                output.close();
            }

        }

    private:
        shared_ptr<StructuredBlockStorage> blocks_;
        const BlockDataID bodyStorageID_;
        const uint_t checkFrequency_;
        uint_t executionCounter_;
        const real_t preFac_;
    };
/*
// CALCULATE AVERAGE NUMBER OF COLLISIONS
    class CollisionCounter {

    public:
        CollisionCounter(const shared_ptr<FBfunc::TimeStep> peStep, SetupFB &config) :
                peStep_(peStep), checkFrequency_((config.evaluationFreq > 0) ? config.evaluationFreq : uint_c(1)), executionCounter_(uint_c(0)) {
            std::ofstream output;
            output.open("Collisions.txt", std::ofstream::out | std::ofstream::trunc);
            output.close();
        }

        void operator()() {

            ++executionCounter_;
            if ((executionCounter_ - uint_c(1)) % checkFrequency_ != 0)
                return;

            real_t avgCollisions = real_c(0);
            avgCollisions = peStep_->getAverageContacts();

            WALBERLA_MPI_SECTION() {
                mpi::reduceInplace(avgCollisions, mpi::SUM);
            }
            WALBERLA_ROOT_SECTION() {
                std::ofstream output;
                output.open("Collisions.txt", std::ofstream::out | std::ofstream::app);
                output << executionCounter_ << "\t" << avgCollisions << std::endl;
                output.close();
            }
        }

    private:
        const shared_ptr<FBfunc::TimeStep> peStep_;
        const uint_t checkFrequency_;
        uint_t executionCounter_;
    };
*/
// CALCULATE PRESSURE WITH DENSITY
    template<typename LatticeModel_T, typename FlagField_T>
    class PressureByDensity {

    public:

        typedef lbm::PdfField<LatticeModel_T> PdfField;

        PressureByDensity(const shared_ptr<StructuredBlockStorage> &blocks, const ConstBlockDataID &pdfFieldId, const ConstBlockDataID &flagFieldId,
                          const CellInterval &globalLowerPlane, const CellInterval &globalUpperPlane, const FlagUID &flagToCheck, SetupFB & config) :

                blocks_(blocks), pdfFieldId_(pdfFieldId), flagFieldId_(flagFieldId),  globalLowerPlane_(globalLowerPlane), globalUpperPlane_(globalUpperPlane),
                flagToCheck_(flagToCheck),
                checkFrequency_((config.evaluationFreq > 0) ? config.evaluationFreq : uint_c(1)), executionCounter_(uint_c(0)),
                unitConversion_(config.densityFluid_SI * config.dx_SI * config.dx_SI / (config.dt_SI * config.dt_SI)),
                numCellsInPlane_(real_c(globalLowerPlane.xSize() * globalLowerPlane.zSize())) {

            std::ofstream output;
            output.open("PressuredropThroughDensity.txt", std::ofstream::out | std::ofstream::trunc);
            output.close();
        }

        void operator()() {
            ++executionCounter_;
            if ((executionCounter_ - uint_c(1)) % checkFrequency_ != 0)
                return;

            real_t lowerDensitySum = real_c(0);
            real_t upperDensitySum = real_c(0);

            for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt) {

                const FlagField_T *flagField = blockIt->getData<FlagField_T>(flagFieldId_);
                const PdfField *pdfField = blockIt->getData<PdfField>(pdfFieldId_);
                const flag_t fluid = flagField->getFlag(flagToCheck_);

                CellInterval localLowerPlane = CellInterval();
                CellInterval localUpperPlane = CellInterval();

                blocks_->transformGlobalToBlockLocalCellInterval(localLowerPlane, *blockIt, globalLowerPlane_);
                blocks_->transformGlobalToBlockLocalCellInterval(localUpperPlane, *blockIt, globalUpperPlane_);

                localLowerPlane.intersect(pdfField->xyzSize());
                localUpperPlane.intersect(pdfField->xyzSize());

                for (cell_idx_t z = localLowerPlane.zMin(); z <= localLowerPlane.zMax(); ++z) {
                    for (cell_idx_t y = localLowerPlane.yMin(); y <= localLowerPlane.yMax(); ++y) {
                        for (cell_idx_t x = localLowerPlane.xMin(); x <= localLowerPlane.xMax(); ++x) {
                            if (flagField->isFlagSet(x, y, z, fluid)) {
                                lowerDensitySum += pdfField->getDensity(x, y, z);
                            }
                        }
                    }
                }
                for (cell_idx_t z = localUpperPlane.zMin(); z <= localUpperPlane.zMax(); ++z) {
                    for (cell_idx_t y = localUpperPlane.yMin(); y <= localUpperPlane.yMax(); ++y) {
                        for (cell_idx_t x = localUpperPlane.xMin(); x <= localUpperPlane.xMax(); ++x) {
                            if (flagField->isFlagSet(x, y, z, fluid)) {
                                upperDensitySum += pdfField->getDensity(x, y, z);
                            }
                        }
                    }
                }
            }

            WALBERLA_MPI_SECTION() {
                mpi::reduceInplace(lowerDensitySum, mpi::SUM);
                mpi::reduceInplace(upperDensitySum, mpi::SUM);
            }
            WALBERLA_ROOT_SECTION() {
                real_t pressureDifference = ((lowerDensitySum - upperDensitySum) / numCellsInPlane_ / 3) * unitConversion_;
                std::ofstream output;
                output.open("PressuredropThroughDensity.txt", std::ofstream::out | std::ofstream::app);
                output << executionCounter_ << "\t" << pressureDifference << "\t[Pa]" << std::endl;
                output.close();
            }
        }

    private:

        shared_ptr<StructuredBlockStorage> blocks_;
        const ConstBlockDataID pdfFieldId_;
        const ConstBlockDataID flagFieldId_;
        const CellInterval globalLowerPlane_;
        const CellInterval globalUpperPlane_;
        const FlagUID flagToCheck_;
        const uint_t checkFrequency_;
        uint_t executionCounter_;
        const real_t unitConversion_;
        const real_t numCellsInPlane_;
    };

// OBSTACLE LINEAR VELOCITY CHECK
    class ObstacleLinearVelocityCheck {

    public:

        ObstacleLinearVelocityCheck(const shared_ptr<StructuredBlockStorage> &blocks, const BlockDataID &bodyStorageID,
                                    SetupFB &config) :
                blocks_(blocks), bodyStorageID_(bodyStorageID), uMax_(config.stopVel * config.stopVel), radius_(0.5 * config.diameterA),
                executionCounter_(uint_c(0)), checkFrequency_((config.velCheckFreq > 0) ? config.velCheckFreq : uint_c(1)) {
            std::ofstream output;
            output.open("MaxBodyVelocity.txt", std::ofstream::out | std::ofstream::trunc);
            output.close();
        }

        void operator()() {
            ++executionCounter_;
            if ((executionCounter_ - uint_c(1)) % checkFrequency_ != 0)
                return;

            real_t maxTransVel = real_c(0.0);
            real_t maxAngularVel = real_c(0.0);
            for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt) {
                for (auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt) {
                    if (bodyIt->isFixed() || !bodyIt->isFinite())
                        continue;

                    const auto &u = bodyIt->getLinearVel();
                    const auto &w = bodyIt->getAngularVel();
                    WALBERLA_CHECK((u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) <= uMax_, "Maximum body velocity is " << u.length() << std::endl
                                                                                                                   << "x " << u[0] << " y " << u[1] << " z " << u[2] << std::endl);
                    if (maxTransVel < u.sqrLength()) {
                        maxTransVel = u.sqrLength();
                    }
                    if (maxAngularVel < w.sqrLength()) {
                        maxAngularVel = w.sqrLength();
                    }
                }
            }

            WALBERLA_MPI_SECTION() {
                mpi::reduceInplace(maxTransVel, mpi::MAX);
                mpi::reduceInplace(maxAngularVel, mpi::MAX);
            }
            WALBERLA_ROOT_SECTION() {
                std::ofstream output;
                output.open("MaxBodyVelocity.txt", std::ofstream::out | std::ofstream::app);
                output << executionCounter_ << "\t" << "Max translational velocity: " << sqrt(maxTransVel) << "\t"
                       << "Max tangential velocity: " << sqrt(maxAngularVel) * radius_ << std::endl;
                output.close();
            }

        }

    private:

        shared_ptr<StructuredBlockStorage> blocks_;
        const BlockDataID bodyStorageID_;
        const real_t uMax_;
        const real_t radius_;
        uint_t executionCounter_;
        const uint_t checkFrequency_;
    };

    struct TimingPoolLogger {
        TimingPoolLogger(WcTimingPool &timingPool, const SweepTimeloop &timeloop, const uint_t interval)
                : timingPool_(timingPool), timeloop_(timeloop), interval_(interval) {
        }

        void operator()() {
            if (timeloop_.getCurrentTimeStep() > uint_t(10) && timeloop_.getCurrentTimeStep() % interval_ == 0) {
                timingPool_.logResultOnRoot();
            }
        }

    private:
        WcTimingPool &timingPool_;
        const SweepTimeloop &timeloop_;
        uint_t interval_;
    };

    // CALCULATE AREA FRACTION
    template<typename FlagField_T>
    class GasAreaFractionCalculator {
    public:

        GasAreaFractionCalculator(const shared_ptr<StructuredBlockStorage> &blocks, const ConstBlockDataID &flagFieldId,
                                  FBfunc::SetupFB &config) :
                blocks_(blocks), flagFieldId_(flagFieldId), areas_(config.probesCell),
                checkFrequency_((config.evaluationFreq > 0) ? config.evaluationFreq : uint_c(1)), executionCounter_(uint_c(0)) {

            for (std::vector<CellInterval>::const_iterator it = areas_.begin(); it != areas_.end(); ++it) {
                it ->  zMin() =  blocks->getDomainCellBB().zMin();
                it ->  zMax() =  blocks->getDomainCellBB().zMax();
            }

        }

        void operator()() {
            ++executionCounter_;
            if ((executionCounter_ - uint_c(1)) % checkFrequency_ != 0)
                return;

            uint_t probeCounter = uint_c(1);
            CellInterval globalFieldCellBB;

            for (std::vector<CellInterval>::const_iterator it = areas_.begin(); it != areas_.end(); ++it) {

                std::vector<uint_t> area = std::vector<uint_t>(it->xSize() * it->ySize(), uint_c(0));

                for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt) {

                    const FlagField_T *flagField = blockIt->getData<FlagField_T>(flagFieldId_);
                    const flag_t solid = flagField->getFlag(FBfunc::MO_Flag);

                    blocks_->transformBlockLocalToGlobalCellInterval(globalFieldCellBB, *blockIt, flagField->xyzSize());
                    globalFieldCellBB.intersect(*it);

                    for (cell_idx_t x = it->xMin(); x <= it->xMax(); ++x) {
                        for (cell_idx_t y = it->yMin(); y <= it->yMax(); ++y) {
                            for (cell_idx_t z = it->zMin(); z <= it->zMax(); ++z) {
                                if (globalFieldCellBB.contains(x, y, z)) {
                                    if (flagField->isFlagSet(x, y, z, solid)) {
                                        area[(x - it->xMin()) + ((y - it->yMin()) * it->xSize())] += uint_c(1);
                                    }
                                }
                            }
                        }
                    }
                }

                WALBERLA_MPI_SECTION() {
                    mpi::reduceInplace(area, mpi::SUM);
                }

                WALBERLA_ROOT_SECTION() {
                    uint_t linesWithoutSolids = uint_c(0);
                    for (std::vector<uint_t>::const_iterator lines = area.begin(); lines !=  area.end(); ++lines){
                        if ( *lines == uint_c(0) ){
                            ++linesWithoutSolids;
                        }
                    }
                    std::ofstream output;
                    output.open("GasAreaFraction" + std::to_string(probeCounter) + ".txt", std::ofstream::out | std::ofstream::app);
                    output << executionCounter_ << "\t" << linesWithoutSolids / area.size() << std::endl;
                    output.close();
                }

                ++probeCounter;
            }
        }

    private:

        shared_ptr<StructuredBlockStorage> blocks_;
        const ConstBlockDataID flagFieldId_;
        const std::vector<CellInterval> areas_;
        const uint_t checkFrequency_;
        uint_t executionCounter_;

    }; // Area Calculator

/////////////////////////////////
// OUTSIDE TIME LOOP FUNCTIONS
/////////////////////////////////

// SET VELOCITY FORCE AND TORQUE OF SPHERES TO ZERO
    void restSphere(const shared_ptr<StructuredBlockStorage> &blocks, const BlockDataID &bodyStorageID) {
        for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) {
            for (auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end(); ++bodyIt) {
                bodyIt->setLinearVel(0.0, 0.0, 0.0);
                bodyIt->resetForceAndTorque();
            }
        }
    }

// GET MAXIMUM VEL OF PARTICLES
    real_t getMaxVel(const shared_ptr<StructuredBlockStorage> &blocks, const BlockDataID &bodyStorageID) {

        real_t maxVel = real_c(0.0);

        for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) {
            for (auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt) {
                real_t sphereVel = bodyIt->getLinearVel().sqrLength();
                if (maxVel < sphereVel) {
                    maxVel = sphereVel;
                }
            }
        }

        WALBERLA_MPI_SECTION() {
            mpi::allReduceInplace(maxVel, mpi::MAX);
        }

        return sqrt(maxVel);
    }

// GET AVERAGE VEL OF PARTICLES
    real_t getAvgVel(const shared_ptr<StructuredBlockStorage> &blocks,
                     const BlockDataID &bodyStorageID, const uint_t nrParticles) {

        real_t avgVel = real_c(0.0);

        for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) {
            for (auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt) {
                avgVel += bodyIt->getLinearVel().length();
            }
        }


        WALBERLA_MPI_SECTION() {
            mpi::allReduceInplace(avgVel, mpi::SUM);
        }

        return avgVel / real_c(nrParticles);
    }

    // GET SUMMED FORCE Y DIRECTION
    real_t getForceY(const shared_ptr<StructuredBlockStorage> &blocks,
                     const BlockDataID &bodyStorageID) {

        real_t forceSum = real_c(0.0);

        for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) {
            for (auto bodyIt = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt) {
                forceSum += bodyIt->getForce()[1];
            }
        }

        WALBERLA_MPI_SECTION() {
            mpi::allReduceInplace(forceSum, mpi::SUM);
        }

        return forceSum;
    }


    // GET AVERAGE FLUID VELOCITY
    template<typename LatticeModel_T, typename FlagField_T, typename BodyField_T>
    class AverageFluidVelocity {
    public:

        typedef lbm::PdfField<LatticeModel_T> PdfField;

        AverageFluidVelocity(const shared_ptr<StructuredBlockStorage> &blocks, const BlockDataID &bodyStorageID, const ConstBlockDataID &bodyFieldId,
                             const ConstBlockDataID &pdfFieldId, const ConstBlockDataID &flagFieldId, const Set<FlagUID> &cellsToCheck) :
                blocks_(blocks), bodyStorageID_(bodyStorageID), bodyFieldId_(bodyFieldId), pdfFieldId_(pdfFieldId), flagFieldId_(flagFieldId), cellsToCheck_(cellsToCheck)
        {}

        real_t getAverageFluidLbmVelocity() {

            Vector3 <real_t> velSum   = Vector3 <real_t> (0.0,0.0,0.0);
            uint_t           numCells = uint_c(0);

            for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt) {

                const FlagField_T *flagField = blockIt->getData<FlagField_T>(flagFieldId_);
                const PdfField    *pdfField  = blockIt->getData<PdfField>(pdfFieldId_);
                //const BodyField_T *bodyField = blockIt->getData<BodyField_T>(bodyFieldId_);

                const CellInterval &cellBB = pdfField->xyzSize();
                WALBERLA_ASSERT_EQUAL(flagField->xyzSize(), cellBB);

                typename FlagField_T::flag_t mask = 0;
                for (auto flag = cellsToCheck_.begin(); flag != cellsToCheck_.end(); ++flag)
                    mask = static_cast< typename FlagField_T::flag_t >( mask | flagField->getFlag(*flag));

                for (cell_idx_t z = cellBB.zMin(); z <= cellBB.zMax(); ++z) {
                    for (cell_idx_t y = cellBB.yMin(); y <= cellBB.yMax(); ++y) {
                        for (cell_idx_t x = cellBB.xMin(); x <= cellBB.xMax(); ++x) {
                            if (flagField->isPartOfMaskSet(x, y, z, mask)) {

                                velSum += pdfField->getVelocity(x, y, z);
                                numCells += 1;

                            }
                        }
                    }
                }
            }

            WALBERLA_MPI_SECTION() {
                mpi::allReduceInplace(velSum[0], mpi::SUM);
                mpi::allReduceInplace(velSum[1], mpi::SUM);
                mpi::allReduceInplace(velSum[2], mpi::SUM);
                mpi::allReduceInplace(numCells , mpi::SUM);
            }

            return (velSum.length()/real_c(numCells));
        }

    private:
        shared_ptr<StructuredBlockStorage> blocks_;
        const BlockDataID bodyStorageID_;
        const ConstBlockDataID bodyFieldId_;
        const ConstBlockDataID pdfFieldId_;
        const ConstBlockDataID flagFieldId_;
        const Set<FlagUID> cellsToCheck_;

    }; //

    template<typename LatticeModel_T, typename FlagField_T, typename BodyField_T>
    class VelocityCheckMEM {
    public:

        typedef lbm::PdfField<LatticeModel_T> PdfField;

        VelocityCheckMEM(const shared_ptr<StructuredBlockStorage> &blocks, const BlockDataID &bodyStorageID, const ConstBlockDataID &bodyFieldId,
                         const ConstBlockDataID &pdfFieldId, const ConstBlockDataID &flagFieldId, const Set<FlagUID> &cellsToCheck,
                         FBfunc::SetupFB & config) :
                blocks_(blocks), bodyStorageID_(bodyStorageID), bodyFieldId_(bodyFieldId), pdfFieldId_(pdfFieldId), flagFieldId_(flagFieldId), cellsToCheck_(cellsToCheck),
                uMax_(config.stopVel * config.stopVel), checkFrequency_((config.velCheckFreq > 0) ? config.velCheckFreq : uint_c(1)), dx_SI_(config.dx_SI), executionCounter_(uint_c(0)) {

            std::ofstream output;
            output.open("MaximalLBMvelocity.txt", std::ofstream::out | std::ofstream::trunc);
            output.close();
        }

        void operator()(const IBlock *const block) {
            ++executionCounter_;
            if ((executionCounter_ - uint_c(1)) % checkFrequency_ != 0)
                return;

            const FlagField_T *flagField = block->getData<FlagField_T>(flagFieldId_);
            const PdfField *pdfField = block->getData<PdfField>(pdfFieldId_);
            const BodyField_T *bodyField = block->getData<BodyField_T>(bodyFieldId_);


            const CellInterval &cellBB = pdfField->xyzSize();
            WALBERLA_ASSERT_EQUAL(flagField->xyzSize(), cellBB);

            real_t currentMax = real_c(0.0);

            typename FlagField_T::flag_t mask = 0;
            for (auto flag = cellsToCheck_.begin(); flag != cellsToCheck_.end(); ++flag)
                mask = static_cast< typename FlagField_T::flag_t >( mask | flagField->getFlag(*flag));

            for (cell_idx_t z = cellBB.zMin(); z <= cellBB.zMax(); ++z) {
                for (cell_idx_t y = cellBB.yMin(); y <= cellBB.yMax(); ++y) {
                    for (cell_idx_t x = cellBB.xMin(); x <= cellBB.xMax(); ++x) {
                        if (flagField->isPartOfMaskSet(x, y, z, mask)) {
                            if (pdfField->getVelocity(x, y, z).sqrLength() > uMax_) {

                                Cell globalPos;
                                blocks_->transformBlockLocalToGlobalCell(globalPos, *block, Cell(x, y, z));
                                std::ofstream output;
                                output.open("BrokenArea" + std::to_string(MPIManager::instance()->rank()) + ".txt");

                                output << "Position of maximum value is x: " << globalPos.x() * dx_SI_ * 1000 << "   y: " << globalPos.y() * dx_SI_ * 1000 << "   z: "
                                       << globalPos.z() * dx_SI_ * 1000 << " [mm]" << std::endl;
                                output << "Broken velocity is " << pdfField->getVelocity(x, y, z).length() << std::endl << std::endl;

                                for (auto dir = stencil::D3Q27::beginNoCenter(); dir != stencil::D3Q27::end(); ++dir) {

                                    if (flagField->isPartOfMaskSet(x + dir.cx(), y + dir.cy(), z + dir.cz(), mask)) {
                                        output << "Cell in direction\t" << dir.cx() << "\t" << dir.cy() << "\t" << dir.cz() << "\thas velocity "
                                               << pdfField->getVelocity(x + dir.cx(), y + dir.cy(), z + dir.cz()).length() << std::endl;
                                    } else if (flagField->isPartOfMaskSet(x + dir.cx(), y + dir.cy(), z + dir.cz(),
                                                                          static_cast< typename FlagField_T::flag_t >( mask | flagField->getFlag(FBfunc::MO_Flag)))) {
                                        Vector3<real_t> gpos;
                                        blocks_->getBlockLocalCellCenter(*block, Cell(x + dir.cx(), y + dir.cy(), z + dir.cz()), gpos);
                                        blocks_->transformBlockLocalToGlobal(gpos, *block);
                                        output << "Cell in direction\t" << dir.cx() << "\t" << dir.cy() << "\t" << dir.cz() << "\tis solid with Body ID: "
                                               << (*bodyField)(x + dir.cx(), y + dir.cy(), z + dir.cz())->getID()
                                               << " and velocity: " << ((*bodyField)(x + dir.cx(), y + dir.cy(), z + dir.cz())->velFromWF(gpos)).length() << std::endl;
                                    } else {
                                        output << "Cell in direction\t" << dir.cx() << "\t" << dir.cy() << "\t" << dir.cz() << "\t is neither fluid nor particle" << std::endl;
                                    }
                                }
                                output.close();

                                WALBERLA_CHECK(false, "Position of maximum value is x: " << globalPos.x() * dx_SI_ * 1000 << "   y: " << globalPos.y() * dx_SI_ * 1000
                                                                                         << "   z: " << globalPos.z() * dx_SI_ * 1000 << " [mm]" << std::endl << std::endl);
                            }

                            if (currentMax < pdfField->getVelocity(x, y, z).sqrLength()) {
                                currentMax = pdfField->getVelocity(x, y, z).sqrLength();
                            }
                        }
                    }
                }
            }
            WALBERLA_MPI_SECTION() {
                mpi::reduceInplace(currentMax, mpi::MAX);
            }
            WALBERLA_ROOT_SECTION() {
                //std::cout << "Maximal LBM Velocity: " << sqrt(currentMax) << std::endl;
                std::ofstream output;
                output.open( "MaximalLBMvelocity.txt", std::ofstream::out | std::ofstream::app);
                output << executionCounter_ << "\t" << sqrt(currentMax) << std::endl;
                output.close();
            }
        }

    private:
        shared_ptr<StructuredBlockStorage> blocks_;
        const BlockDataID bodyStorageID_;
        const ConstBlockDataID bodyFieldId_;
        const ConstBlockDataID pdfFieldId_;
        const ConstBlockDataID flagFieldId_;
        const Set<FlagUID> cellsToCheck_;
        const real_t uMax_;
        const uint_t checkFrequency_;
        const real_t dx_SI_;
        uint_t executionCounter_;

    }; // VelocityCheckMEM

    // CALCULATE SOLID FRACTION
    template<typename FlagField_T>
    class FractionCalculator {
    public:

        FractionCalculator(const shared_ptr<StructuredBlockStorage> &blocks, const ConstBlockDataID &flagFieldId,
                              FBfunc::SetupFB &config) :
                blocks_(blocks), flagFieldId_(flagFieldId), probes_(config.probesCell),
                checkFrequency_((config.evaluationFreq > 0) ? config.evaluationFreq : uint_c(1)), executionCounter_(uint_c(0)) {
        }

        void operator()() {
            ++executionCounter_;
            if ((executionCounter_ - uint_c(1)) % checkFrequency_ != 0)
                return;

            uint_t probeCounter = uint_c(1);
            CellInterval localBox = CellInterval();

            for (std::vector<CellInterval>::const_iterator it = probes_.begin(); it != probes_.end(); ++it) {

                real_t boxVolume = real_c(it->numCells());
                real_t numberMOcells = real_t(0);

                for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt) {

                    const FlagField_T *flagField = blockIt->getData<FlagField_T>(flagFieldId_);
                    const flag_t solid           = flagField->getFlag(FBfunc::MO_Flag);

                    blocks_->transformGlobalToBlockLocalCellInterval(localBox, *blockIt, *it);
                    localBox.intersect(flagField->xyzSize());

                    for (cell_idx_t z = localBox.zMin(); z <= localBox.zMax(); ++z) {
                        for (cell_idx_t y = localBox.yMin(); y <= localBox.yMax(); ++y) {
                            for (cell_idx_t x = localBox.xMin(); x <= localBox.xMax(); ++x) {
                                if (flagField->isFlagSet(x, y, z, solid)) {
                                    ++numberMOcells;
                                }
                            }
                        }
                    }
                }
                WALBERLA_MPI_SECTION() {
                    mpi::reduceInplace(numberMOcells, mpi::SUM);
                }
                WALBERLA_ROOT_SECTION() {
                    std::ofstream output;
                    output.open("SolidFraction" + std::to_string(probeCounter) + ".txt", std::ofstream::out | std::ofstream::app);
                    output << executionCounter_ << "\t" << numberMOcells / boxVolume << std::endl;
                    output.close();
                }

                ++probeCounter;
            }
        }

    private:

        shared_ptr<StructuredBlockStorage> blocks_;
        const ConstBlockDataID flagFieldId_;
        const std::vector<CellInterval> probes_;
        const uint_t checkFrequency_;
        uint_t executionCounter_;

    }; // Fraction Calculator

// CHECK POINT CREATER

    class CheckpointCreater {

    public:
        CheckpointCreater( BlockStorage & blockStorage,
                           const BlockDataID & pdfFieldID,
                           const BlockDataID & bodyStorageID,
                           FBfunc::SetupFB &config)
                : blockStorage_(blockStorage), pdfFieldID_(pdfFieldID), bodyStorageID_(bodyStorageID),
                  checkFrequency_((config.checkpointFreq > 0) ? config.checkpointFreq : uint_c(1)), executionCounter_(uint_c(0))
        {
        }

        void operator()() {
            ++executionCounter_;
            if ((executionCounter_ ) % checkFrequency_ != 0)
                return;

            WALBERLA_LOG_RESULT_ON_ROOT("Writing data to file");

            blockStorage_.saveBlockData("checkpoint_bodies.dat", bodyStorageID_);
            blockStorage_.saveBlockData("checkpoint_pdf.dat", pdfFieldID_);

        }
    private:
        BlockStorage & blockStorage_;
        const BlockDataID pdfFieldID_;
        const BlockDataID bodyStorageID_;

        uint_t checkFrequency_;
        uint_t executionCounter_;
    };


} //namespace

