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
//! \file FluidizedBedMEM.cpp
//! \author Dominik Schuster <dominik.schuster@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/all.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/SharedFunctor.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"
#include "core/math/Random.h"
#include "core/timing/RemainingTimeLogger.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"
#include "field/StabilityChecker.h"

#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SplitPureSweep.h"
#include "lbm/lattice_model/SmagorinskyLES.h"
#include "lbm/PerformanceLogger.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "pe/basic.h"
#include "pe/cr/ICR.h"
#include "pe/Types.h"
#include "pe/vtk/BodyVtkOutput.h"
#include "pe/vtk/SphereVtkOutput.h"

#include "pe_coupling/mapping/all.h"
#include "pe_coupling/momentum_exchange_method/all.h"
#include "pe_coupling/utility/all.h"

#include <random>

#include "vtk/all.h"
#include "field/vtk/all.h"
#include "lbm/vtk/all.h"

#include "FluidizedBedMEMFunctionality.h"

namespace fluidizedBed {

///////////
// USING //
///////////

    using namespace walberla;
    using walberla::uint_t;

    const uint_t FieldGhostLayers = 1;

    typedef GhostLayerField<real_t, FieldGhostLayers> ScalarField_T;
    typedef GhostLayerField<pe::BodyID, FieldGhostLayers> BodyField_T;

    typedef walberla::uint8_t flag_t;
    typedef FlagField<flag_t> FlagField_T;

    typedef boost::tuple<pe::Sphere, pe::Plane> BodyTypeTuple;


/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////

// class MyBoundaryHandling is responsible for creating the boundary handling and setting up the geometry at the outer domain boundaries

    template <typename BoundaryHandling_T, typename FlagField_T, typename PdfField_T, typename NoSlip_T, typename SimplePressure_T, typename UBB_T, typename MO_T>
    class MyBoundaryHandling {
    public:

        MyBoundaryHandling(const BlockDataID &flagField, const BlockDataID &pdfField, const BlockDataID &bodyField,
                //const shared_ptr<Vector3<real_t>> velocity,
                           const bool xPeriodic, const bool zPeriodic,
                           const bool spotInflow, const real_t spotDiameter, const real_t xCenter, const real_t  zCenter,
                           const real_t velocityMax, const real_t initialRamping, const real_t rampingTime,
                           const shared_ptr< lbm::TimeTracker > & timeTracker):
                flagField_(flagField), pdfField_(pdfField), bodyField_(bodyField)
                ,xPeriodic_(xPeriodic), zPeriodic_(zPeriodic), spotInflow_(spotInflow), spotDiameter_(spotDiameter)
                ,xCenter_(xCenter), zCenter_(zCenter), velocityMax_(velocityMax), initialRamping_(initialRamping)
                ,rampingTime_(rampingTime), timeTracker_(timeTracker)
                ,velocityFunctor_(FBfunc::VelocityFunctor_T (xCenter, zCenter, spotDiameter, velocityMax, initialRamping, rampingTime, spotInflow))
        {
        }

        BoundaryHandling_T *operator()(IBlock *const block, const StructuredBlockStorage *const storage) const;

    private:

        const BlockDataID flagField_;
        const BlockDataID pdfField_;
        const BlockDataID bodyField_;
        const bool xPeriodic_;
        const bool zPeriodic_;
        const bool spotInflow_;
        const real_t spotDiameter_;
        const real_t xCenter_;
        const real_t zCenter_;
        const real_t velocityMax_;
        const real_t initialRamping_;
        const real_t rampingTime_;
        const shared_ptr< lbm::TimeTracker > timeTracker_;
        FBfunc::VelocityFunctor_T velocityFunctor_;

    }; // class MyBoundaryHandling


    template <typename BoundaryHandling_T, typename FlagField_T, typename PdfField_T, typename NoSlip_T, typename SimplePressure_T, typename UBB_T, typename MO_T>
    BoundaryHandling_T *MyBoundaryHandling<BoundaryHandling_T, FlagField_T, PdfField_T, NoSlip_T, SimplePressure_T, UBB_T, MO_T>::operator()(IBlock *const block, const StructuredBlockStorage *const storage) const {
        WALBERLA_ASSERT_NOT_NULLPTR(block);
        WALBERLA_ASSERT_NOT_NULLPTR(storage);

        FlagField_T *flagField = block->getData<FlagField_T>(flagField_);
        PdfField_T *pdfField = block->getData<PdfField_T>(pdfField_);
        BodyField_T *bodyField = block->getData<BodyField_T>(bodyField_);

        const auto fluid = flagField->flagExists(FBfunc::Fluid_Flag) ? flagField->getFlag(FBfunc::Fluid_Flag) : flagField->registerFlag(FBfunc::Fluid_Flag);

        BoundaryHandling_T *handling = new BoundaryHandling_T("cf boundary handling", flagField, fluid,
                                                              boost::tuples::make_tuple(NoSlip_T("NoSlip", FBfunc::NoSlip_Flag, pdfField),
                                                                                        SimplePressure_T("Outlet", FBfunc::Outlet_Flag, pdfField, real_c(1.0)),
                                                                                        UBB_T("UBB", FBfunc::UBB_Flag, pdfField, timeTracker_, uint_t(0), velocityFunctor_, block->getAABB() ),
                                                                                        MO_T("MO", FBfunc::MO_Flag, pdfField, flagField, bodyField, fluid, *storage, *block)));

        const auto ubb = flagField->getFlag(FBfunc::UBB_Flag);
        const auto noSlip = flagField->getFlag(FBfunc::NoSlip_Flag);
        const auto outlet = flagField->getFlag(FBfunc::Outlet_Flag);

        CellInterval domainBB = storage->getDomainCellBB();

        domainBB.xMin() -= cell_idx_c(FieldGhostLayers);
        domainBB.xMax() += cell_idx_c(FieldGhostLayers);

        if (!xPeriodic_) {
// LEFT
            CellInterval left(domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMin(), domainBB.yMax(), domainBB.zMax());
            storage->transformGlobalToBlockLocalCellInterval(left, *block);
            handling->forceBoundary(noSlip, left);

// RIGHT
            CellInterval right(domainBB.xMax(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax());
            storage->transformGlobalToBlockLocalCellInterval(right, *block);
            handling->forceBoundary(noSlip, right);
        }

        domainBB.zMin() -= cell_idx_c(FieldGhostLayers);
        domainBB.zMax() += cell_idx_c(FieldGhostLayers);

        if (!zPeriodic_) {
// FRONT
            CellInterval front(domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMin());
            storage->transformGlobalToBlockLocalCellInterval(front, *block);
            handling->forceBoundary(noSlip, front);

// BACK
            CellInterval back(domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax());
            storage->transformGlobalToBlockLocalCellInterval(back, *block);
            handling->forceBoundary(noSlip, back);
        }


        domainBB.yMin() -= cell_idx_c(FieldGhostLayers);
        domainBB.yMax() += cell_idx_c(FieldGhostLayers);

// BOTTOM

        if (spotInflow_){

            CellInterval bottom(domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMin(), domainBB.zMax());
            storage->transformGlobalToBlockLocalCellInterval(bottom, *block);
            handling->forceBoundary(noSlip, bottom);

            for (cell_idx_t xIterator = domainBB.xMin(); xIterator <= domainBB.xMax(); ++xIterator ){
                for (cell_idx_t zIterator = domainBB.zMin(); zIterator <= domainBB.zMax(); ++zIterator ){
                    if ( sqrt ( (real_c(xIterator) - xCenter_)*(real_c(xIterator) - xCenter_) + (real_c(zIterator) - zCenter_)*(real_c(zIterator) - zCenter_) ) <= spotDiameter_ * 0.5 ){
                        CellInterval inflow (xIterator, domainBB.yMin(), zIterator, xIterator, domainBB.yMin(), zIterator);
                        storage->transformGlobalToBlockLocalCellInterval( inflow , *block);
                        handling->forceBoundary(ubb, inflow);
                    }
                }

            }


        } else {

            CellInterval bottom(domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMin(), domainBB.zMax());
            storage->transformGlobalToBlockLocalCellInterval(bottom, *block);
            handling->forceBoundary(ubb, bottom);

        }



// TOP
        CellInterval top(domainBB.xMin(), domainBB.yMax(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax());
        storage->transformGlobalToBlockLocalCellInterval(top, *block);
        handling->forceBoundary(outlet, top);

        handling->fillWithDomain(FieldGhostLayers);

        return handling;
    }

    template <typename CollisionModel_T>
    struct LatticeModelCreator
    {

    };

    template<>
    struct LatticeModelCreator< lbm::collision_model::TRT >
    {
        static lbm::D3Q19<lbm::collision_model::TRT, false> makeLatticeModel(BlockDataID /*omegaFieldID*/, real_t omega)
        {
            return lbm::D3Q19<lbm::collision_model::TRT, false> (lbm::collision_model::TRT::constructWithMagicNumber(omega));
        }
    };


    template<>
    struct LatticeModelCreator< lbm::collision_model::SRT >
    {
        static lbm::D3Q19<lbm::collision_model::SRT, false> makeLatticeModel(BlockDataID /*omegaFieldID*/, real_t omega)
        {
            return lbm::D3Q19<lbm::collision_model::SRT, false> (lbm::collision_model::SRT(omega));
        }
    };

    template<>
    struct LatticeModelCreator< lbm::collision_model::SRTField<ScalarField_T> >
    {

        static lbm::D3Q19<lbm::collision_model::SRTField<ScalarField_T>, false> makeLatticeModel(BlockDataID omegaFieldID, real_t /*omega*/)
        {
            return lbm::D3Q19<lbm::collision_model::SRTField<ScalarField_T>, false> (lbm::collision_model::SRTField<ScalarField_T>(omegaFieldID));
        }
    };


    template <typename  collisionModel_T>
    int runSimulation ( Environment & env) {

        uint_t processes = MPIManager::instance()->numProcesses();

        FBfunc::SetupFB config;

        auto model_parameters = env.config()->getBlock("Model");

        const bool turb = model_parameters.getParameter<bool>("turbulence_model", false);

//////////////
// TYPEDEFS //
//////////////

        typedef lbm::D3Q19<collisionModel_T, false> LatticeModel_T;
        typedef typename LatticeModel_T::Stencil Stencil_T;
        typedef lbm::PdfField<LatticeModel_T> PdfField_T;

// boundary handling

        typedef lbm::NoSlip<LatticeModel_T, flag_t> NoSlip_T;
        typedef lbm::SimplePressure<LatticeModel_T, flag_t> SimplePressure_T;
        typedef lbm::DynamicUBB<LatticeModel_T, flag_t, FBfunc::VelocityFunctor_T> UBB_T;
        typedef pe_coupling::CurvedLinear<LatticeModel_T, FlagField_T> MO_T;

        typedef boost::tuples::tuple<NoSlip_T, SimplePressure_T, UBB_T, MO_T> BoundaryConditions_T;
        typedef BoundaryHandling<FlagField_T, Stencil_T, BoundaryConditions_T> BoundaryHandling_T;

        typedef MyBoundaryHandling<BoundaryHandling_T, FlagField_T, PdfField_T, NoSlip_T, SimplePressure_T, UBB_T, MO_T> MyBoundaryHandling_T;


///////////////////////////
// SIMULATION PROPERTIES //
///////////////////////////

        // GEOMETRY PARAMETERS
        auto geo_parameters = env.config()->getBlock("Geometry");

        const real_t depth_SI = geo_parameters.getParameter<real_t>("depth", real_c(0.02));                         //[m]
        const real_t width_SI = geo_parameters.getParameter<real_t>("width", real_c(0.10));                         //[m]
        const real_t height_SI = geo_parameters.getParameter<real_t>("height", real_c(0.40));                         //[m]
        const real_t offsetTop_SI = geo_parameters.getParameter<real_t>("offset_top", real_c(0.00));                         //[m]
        const real_t offsetBot_SI = geo_parameters.getParameter<real_t>("offset_bot", real_c(0.00));                         //[m]

        if (!geo_parameters.isDefined("depth")) {
            WALBERLA_ABORT("You need to specify \"depth\" in the configuration file");
        }
        if (!geo_parameters.isDefined("width")) {
            WALBERLA_ABORT("You need to specify \"width\" in the configuration file");
        }
        if (!geo_parameters.isDefined("height")) {
            WALBERLA_ABORT("You need to specify \"height\" in the configuration file");
        }

        const real_t crossSectionArea_SI = width_SI * depth_SI;

        const bool xPeriodic = geo_parameters.getParameter<bool>("width_periodic", false);
        const bool zPeriodic = geo_parameters.getParameter<bool>("depth_periodic", false);

        const bool twoSpecies = geo_parameters.getParameter<bool>("two_species", false);
        //const bool createSpeciesWithRandomGenerator = geo_parameters.getParameter<bool>("create_species_with_random_generator", false);

        const real_t diameterA_SI = geo_parameters.getParameter<real_t>("diameter_A", real_c(0.00215));   //[m]
        const real_t diameterB_SI = geo_parameters.getParameter<real_t>("diameter_B", real_c(0.00215));   //[m]
        const real_t radiusStandardDeviationA = 0.5 * geo_parameters.getParameter<real_t>("diameter_standard_distribution_A", real_c(0.00000));
        const real_t radiusStandardDeviationB = 0.5 * geo_parameters.getParameter<real_t>("diameter_standard_distribution_B", real_c(0.00000));
        config.dx_SI = geo_parameters.getParameter<real_t>("cell_size", real_c(0.0002));
        config.diameterA = diameterA_SI / config.dx_SI;
        config.diameterB = diameterB_SI / config.dx_SI;
        const real_t radiusA = real_c(config.diameterA * 0.5);
        const real_t radiusB = real_c(config.diameterB * 0.5);

        const uint_t xBlocks = geo_parameters.getParameter<uint_t>("blocks_width", uint_c(5));
        const uint_t yBlocks = geo_parameters.getParameter<uint_t>("blocks_height", uint_c(20));
        const uint_t zBlocks = geo_parameters.getParameter<uint_t>("blocks_depth", uint_c(1));

        if (processes != (xBlocks * yBlocks * zBlocks)) {
            WALBERLA_ABORT("number of processes must be equal to number of blocks ( here: " << (xBlocks * yBlocks * zBlocks) << " )");
        }

        config.xlength = uint_c(std::round(width_SI / config.dx_SI));
        config.ylength = uint_c(std::round((height_SI + offsetTop_SI + offsetBot_SI) / config.dx_SI));
        config.zlength = uint_c(std::round(depth_SI / config.dx_SI));
        config.offsetTop = offsetBot_SI / config.dx_SI;
        config.offsetBot = offsetBot_SI / config.dx_SI;

        const bool spotInflow = geo_parameters.getParameter<bool>("spot_inflow", false);
        real_t spotDiameter = geo_parameters.getParameter<real_t>("spot_diameter", real_c(0.01));
        spotDiameter = spotDiameter / config.dx_SI;

        // FLUIDIZED BED PARAMETERS
        auto fb_parameters = env.config()->getBlock("FluidizedBed");

        config.nrParticlesA = fb_parameters.getParameter<uint_t>("#particles_A", uint_c(0));
        config.nrParticlesB = fb_parameters.getParameter<uint_t>("#particles_B", uint_c(0));
        const real_t densitySolidA_SI = fb_parameters.getParameter<real_t>("density_particle_A",
                                                                           real_c(22.5));         //[kg/m3] // alumina katalyst
        const real_t densitySolidB_SI = fb_parameters.getParameter<real_t>("density_particle_B", real_c(22.5));         //[kg/m3]
        config.densityFluid_SI = fb_parameters.getParameter<real_t>("density_fluid", real_c(1.2));          //[kg/m3]      //water
        const real_t dynViscosity_SI = fb_parameters.getParameter<real_t>("dynamic_viscosity",
                                                                          real_c(1.84e-5));       //[Pa*s]  // viscosity of water
        const real_t kinViscosity_SI = dynViscosity_SI / config.densityFluid_SI;                                       //[m2/s]
        const real_t gravity_SI = fb_parameters.getParameter<real_t>("gravity", real_c(9.81));

        config.densityRatioA = densitySolidA_SI / config.densityFluid_SI;
        config.densityRatioB = densitySolidB_SI / config.densityFluid_SI;

        real_t velocity_SI;
        real_t particleReynoldsNumber_SI;
        if (fb_parameters.isDefined("velocity")) {
            velocity_SI = fb_parameters.getParameter<real_t>("velocity", real_c(0.35));
            particleReynoldsNumber_SI = velocity_SI * diameterA_SI / kinViscosity_SI;
        } else if (fb_parameters.isDefined("particle_RE")) {
            particleReynoldsNumber_SI = fb_parameters.getParameter<real_t>("particle_RE", real_c(48));
            velocity_SI = particleReynoldsNumber_SI * kinViscosity_SI / diameterA_SI;
        } else if (fb_parameters.isDefined("volume_flow_rate")) {
            real_t volumeFlowRate_SI = fb_parameters.getParameter<real_t>("volume_flow_rate", real_c(2.5));
            velocity_SI = real_c(volumeFlowRate_SI / crossSectionArea_SI / 3600);
            particleReynoldsNumber_SI = velocity_SI * diameterA_SI / kinViscosity_SI;
        } else {
            WALBERLA_ABORT("You need to specify either \"volume_flow_rate\" or \"velocity\" or \"particle_RE\" in the configuration file");
        }

        // SIMULATION PARAMETERS
        auto sim_parameters = env.config()->getBlock("Simulation");

        const bool useLubrication = sim_parameters.getParameter<bool>("use_lubrication", false);

        const bool createPackedBed = sim_parameters.getParameter<bool>("create_packed_bed", false);
        const bool checkMaximumParticleVelocity = sim_parameters.getParameter<bool>("check_maximum_particle_velocity", true);
        const real_t settleVelocity = sim_parameters.getParameter<real_t>("settle_velocity", real_c(0.002));
        const uint_t settleCheckFrequency = sim_parameters.getParameter<uint_t>("settle_check_frequency", uint_t(100));
        const real_t peAccelerator = sim_parameters.getParameter<real_t>("pe_accelerator", real_c(50));

        const bool runPreliminaryTimeloop = sim_parameters.getParameter<bool>("run_preliminary_timeloop", false);
        const real_t eps = sim_parameters.getParameter<real_t>("eps", real_c(0.05));

        const bool initializeFromCheckPointFile = sim_parameters.getParameter<bool>("initialize_from_checkpoint", false);
        const bool enableCheckpointing = sim_parameters.getParameter<bool>("enable_checkpointing", false);
        const std::string checkpointPathBodies = sim_parameters.getParameter<std::string>("path_to_checkpoint_bodies", "checkpoint_bodies.dat");
        const std::string checkpointPathPdf = sim_parameters.getParameter<std::string>("path_to_checkpoint_pdf", "checkpoint_pdf.dat");

        config.evaluationFreq = sim_parameters.getParameter<uint_t>("evaluation_frequency", uint_c(100));
        config.velCheckFreq = sim_parameters.getParameter<uint_t>("velocity_check_frequency", uint_c(10));
        config.checkpointFreq = sim_parameters.getParameter<uint_t>("checkpoint_frequency", uint_c(100000)); //has to be large
        config.stopVel = sim_parameters.getParameter<real_t>("abort_LBM_velocity", real_c(0.50));
        const bool useForceAveraging = sim_parameters.getParameter<bool>("use_force_averaging", true);
        config.lbmSubCycles = useForceAveraging ? real_t(2) : real_t(1);

        if (sim_parameters.isDefined("time_step")) {
            config.dt_SI = sim_parameters.getParameter<real_t>("time_step", real_c(1e-6));
            config.velocity = velocity_SI * config.dt_SI / config.dx_SI;
        } else if (sim_parameters.isDefined("LBM_velocity")) {
            config.velocity = sim_parameters.getParameter("LBM_velocity", real_c(0.06));
            config.dt_SI = config.velocity * config.dx_SI / velocity_SI;
        } else if (sim_parameters.isDefined("omega")) {
            real_t omegaTemp = sim_parameters.getParameter<real_t>("omega", real_c(1.90));
            config.kinViscosity = lbm::collision_model::viscosityFromOmega(omegaTemp);
            config.dt_SI = config.kinViscosity / kinViscosity_SI * config.dx_SI * config.dx_SI;
            config.velocity = velocity_SI * config.dt_SI / config.dx_SI;
        } else if (sim_parameters.isDefined("tau")) {
            real_t tauTemp = sim_parameters.getParameter<real_t>("tau", real_c(0.52));
            config.kinViscosity = lbm::collision_model::viscosityFromOmega(real_c(1.0) / tauTemp);
            config.dt_SI = config.kinViscosity / kinViscosity_SI * config.dx_SI * config.dx_SI;
            config.velocity = velocity_SI * config.dt_SI / config.dx_SI;
        } else {
            WALBERLA_ABORT("You need to specify either \"time_step\" or \"LBM_velocity\" or \"omega\" or \"tau\" in the configuration file");
        }


        const real_t smagorinskyConstant = sim_parameters.getParameter<real_t>("smagorinsky_constant", real_c(0.1));

        uint_t timesteps;
        uint_t timestepsRamping;

        if (sim_parameters.isDefined("simulation_time")) {
            real_t simTime = sim_parameters.getParameter<real_t>("simulation_time", real_c(1.0));
            timesteps = uint_c(simTime / (config.lbmSubCycles * config.dt_SI));
        } else if (sim_parameters.isDefined("#time_steps")) {
            timesteps = sim_parameters.getParameter<uint_t>("#time_steps", uint_c(100000));
        } else {
            WALBERLA_ABORT("You need to specify either \"simulation_time\" or \"#time_steps\"  in the configuration file");
        }

        if (sim_parameters.isDefined("ramping_time")) {
            real_t rampingTime = sim_parameters.getParameter<real_t>("ramping_time", real_c(0.0));
            timestepsRamping = uint_c(rampingTime / (config.lbmSubCycles * config.dt_SI));
        } else {
            timestepsRamping = sim_parameters.getParameter<uint_t>("#ramping_steps", uint_c(0));
        }

        const real_t initialRampingVelocityFactor = sim_parameters.getParameter<real_t>("initial_ramping_scaling", real_c(0.50));

        const real_t initialSphereVel =
                geo_parameters.getParameter<real_t>("initial_sphere_random_velocity", real_c(0.0)) * (config.dt_SI / config.dx_SI);
        const real_t initialSphereDistance =
                geo_parameters.getParameter<real_t>("initial_sphere_distance", real_c(diameterA_SI * 1.5)) / config.dx_SI;


        config.kinViscosity = kinViscosity_SI * config.dt_SI / config.dx_SI / config.dx_SI;
        config.gravity = gravity_SI * config.dt_SI * config.dt_SI / config.dx_SI;

        const real_t FroudeNumber_SI = real_c(sqrt(velocity_SI * velocity_SI / (gravity_SI * diameterA_SI) / (config.densityRatioA - 1)));
        const real_t ArrheniusNumber_SI =
                (densitySolidA_SI - config.densityFluid_SI) * config.densityFluid_SI * diameterA_SI * diameterA_SI * diameterA_SI * gravity_SI /
                (dynViscosity_SI * dynViscosity_SI);

        const real_t particleReynoldsNumber_LBM = config.velocity * config.diameterA / config.kinViscosity;
        const real_t FroudeNumber_LBM = real_c(
                sqrt(config.velocity * config.velocity / (config.gravity * config.diameterA) / (config.densityRatioA - 1)));
        const real_t galileiNumber =
                config.gravity * config.diameterA * config.diameterA * config.diameterA / (config.kinViscosity * config.kinViscosity);

        const real_t UmfFroudenumber = real_c(33.7 * (sqrt(1 + 3.6e-5 * ArrheniusNumber_SI) - 1) / sqrt(ArrheniusNumber_SI));
        const real_t REpmf = real_c(sqrt(33.7 * 33.7 + 0.0494 * ArrheniusNumber_SI) - 33.7);            // Wen You 1966
        const real_t gravityPressure = real_c(config.nrParticlesA) * gravity_SI * (densitySolidA_SI - config.densityFluid_SI)
                                       * (1.0 / 6.0 * math::PI * diameterA_SI * diameterA_SI * diameterA_SI) / crossSectionArea_SI;

        const real_t omega = lbm::collision_model::omegaFromViscosity(config.kinViscosity);

        // MATERIAL PARAMETERS
        auto mat_parameters = env.config()->getBlock("Material");

        const bool useDEM = mat_parameters.getParameter<bool>("use_DEM", true);

        const uint_t HCSITSMaxIterations = mat_parameters.getParameter<uint_t>("HCSITS_max_iterations", uint_c(10));
        const real_t HCSITSRelaxationParameter = mat_parameters.getParameter<real_t>("HCSITS_relaxation_parameter", real_c(0.7));

        const real_t volumeA = real_c(4.0 / 3.0 * math::PI * radiusA * radiusA * radiusA);
        const real_t volumeB = real_c(4.0 / 3.0 * math::PI * radiusB * radiusB * radiusB);
        const real_t MijA = real_c(0.5 * config.densityRatioA * volumeA);
        const real_t MijB = real_c(0.5 * config.densityRatioB * volumeB);

        real_t MijWall = MijA;
        if (twoSpecies) {
            MijWall = 0.5 * (MijA + MijB);
        }

        const real_t cdsA = mat_parameters.getParameter<real_t>("static_friction_particles_A", real_c(0.5));           // static friction
        const real_t cdsB = mat_parameters.getParameter<real_t>("static_friction_particles_B", real_c(0.5));           // static friction
        const real_t cdsWall = mat_parameters.getParameter<real_t>("static_friction_particles_wall", real_c(0.5));        // static friction
        const real_t cdfA = mat_parameters.getParameter<real_t>("dynamic_friction_particles_A", real_c(0.5));          // dynamic friction
        const real_t cdfB = mat_parameters.getParameter<real_t>("dynamic_friction_particles_B", real_c(0.5));          // dynamic friction
        const real_t cdfWall = mat_parameters.getParameter<real_t>("dynamic_friction_wall", real_c(0.5));                 // dynamic friction

        const real_t poisson = real_c(0.2);                  // not used in DEM
        const real_t youngsModulus = 100000;                 // not used in DEM
        //const real_t restitution = real_c(0.8);              // not used in DEM

        // found in: Interface-resolved direct numerical simulation of the erosion of a sediment bed sheared by laminar channel flow
        //           Kidanemariam and Uhlmann
        //           DOI: 10.1016/j.ijmultiphaseflow.2014.08.008

        //const real_t normalizedstiffness = real_c(13000);
        //const real_t stiffness = normalizedstiffness * ((densityRatio-1)*gravity*volume/diameter);

        // found in: Comparison between finite volume and lattice-Boltzmann method simulations of gas-fluidised beds: bed expansion and particle–fluid interaction force
        // R. Third, Y. Chen, C. R. Müller

        const real_t restitutionParticlesA = mat_parameters.getParameter<real_t>("restitution_particles_A", real_c(0.2));
        const real_t restitutionParticlesB = mat_parameters.getParameter<real_t>("restitution_particles_B", real_c(0.2));
        const real_t restitutionWall = mat_parameters.getParameter<real_t>("restitution_wall", real_c(0.2));
        const real_t collisionTimeParticlesA = mat_parameters.getParameter<real_t>("collision_time_particles_A", real_c(10));
        const real_t collisionTimeParticlesB = mat_parameters.getParameter<real_t>("collision_time_particles_B", real_c(10));
        const real_t collisionTimeWall = mat_parameters.getParameter<real_t>("collision_time_wall", real_c(10));
        const real_t tangentialDampingFactorParticlesA = mat_parameters.getParameter<real_t>("tangential_damping_factor_particles_A",
                                                                                             real_c(1.0));
        const real_t tangentialDampingFactorParticlesB = mat_parameters.getParameter<real_t>("tangential_damping_factor_particles_B",
                                                                                             real_c(1.0));
        const real_t tangentialDampingFactorWall = mat_parameters.getParameter<real_t>("tangential_damping_factor_wall", real_c(1.0));
        const real_t timestepsPerCollision = mat_parameters.getParameter<real_t>("timesteps_per_collision", real_c(100));

        const real_t minCollisionTime = std::min(collisionTimeParticlesA, std::min(collisionTimeParticlesB, collisionTimeWall));

        config.peSubCycles = uint_c(std::round(timestepsPerCollision / minCollisionTime));

        // put eq(14) in eq(16) and solve for k_n

        // PARTICLES A
        const real_t stiffnessNparticlesA = (MijA * math::PI * math::PI) /
                                            (collisionTimeParticlesA * collisionTimeParticlesA *
                                             (1.0 - (log(restitutionParticlesA) * log(restitutionParticlesA))
                                                    / (math::PI * math::PI + log(restitutionParticlesA) * log(restitutionParticlesA))));
        const real_t dampingNparticlesA = -2.0 * sqrt(MijA * stiffnessNparticlesA) * log(restitutionParticlesA) /
                                          sqrt(math::PI * math::PI + log(restitutionParticlesA) * log(restitutionParticlesA));
        const real_t dampingTparticlesA = dampingNparticlesA * tangentialDampingFactorParticlesA;

        // PARTICLES B
        const real_t stiffnessNparticlesB = (MijB * math::PI * math::PI) /
                                            (collisionTimeParticlesB * collisionTimeParticlesB *
                                             (1.0 - (log(restitutionParticlesB) * log(restitutionParticlesB))
                                                    / (math::PI * math::PI + log(restitutionParticlesB) * log(restitutionParticlesB))));
        const real_t dampingNparticlesB = -2.0 * sqrt(MijB * stiffnessNparticlesB) * log(restitutionParticlesB) /
                                          sqrt(math::PI * math::PI + log(restitutionParticlesB) * log(restitutionParticlesB));
        const real_t dampingTparticlesB = dampingNparticlesB * tangentialDampingFactorParticlesB;

        // WALL
        const real_t stiffnessNwall = (MijWall * math::PI * math::PI) /
                                      (collisionTimeWall * collisionTimeWall * (1.0 - (log(restitutionWall) * log(restitutionWall))
                                                                                      / (math::PI * math::PI +
                                                                                         log(restitutionWall) * log(restitutionWall))));
        const real_t dampingNwall = -2.0 * sqrt(MijWall * stiffnessNwall) * log(restitutionWall) /
                                    sqrt(math::PI * math::PI + log(restitutionWall) * log(restitutionWall));
        const real_t dampingTwall = dampingNwall * tangentialDampingFactorWall;

        // cdn=cdt

        // EVALUATION PARAMETERS
        auto eval_parameters = env.config()->getBlock("Evaluation");

        const bool calculateCenterOfMass = eval_parameters.getParameter<bool>("mass_center", false);
        const bool calculatePressureDrop = eval_parameters.getParameter<bool>("pressure_drop", false);
        const bool calculateErgunPressure = eval_parameters.getParameter<bool>("ergun_pressure_drop", false);
        const bool calculatePressureDifference = eval_parameters.getParameter<bool>("pressure_difference", false);
        const bool calculateParticleFlux = eval_parameters.getParameter<bool>("particle_flux", false);
        const bool calculateGranularTempGlobal = eval_parameters.getParameter<bool>("granular_temperature_global", false);
        const bool calculateGranularTempLocal = eval_parameters.getParameter<bool>("granular_temperature_local", false);
        const bool calculateSolidVolumeFraction = eval_parameters.getParameter<bool>("solid_volume_fraction", false);
        //const bool calculateCollisions = eval_parameters.getParameter<bool>("collision_count", false);

        const bool box1 = eval_parameters.getParameter<bool>("box1", false);

        real_t box1xSize = eval_parameters.getParameter<real_t>("box1_Size_x", real_t(0.02));
        box1xSize = box1xSize / config.dx_SI;
        real_t box1ySize = eval_parameters.getParameter<real_t>("box1_Size_y", real_t(0.02));
        box1ySize = box1ySize / config.dx_SI;
        real_t box1zSize = eval_parameters.getParameter<real_t>("box1_Size_z", real_t(0.02));
        box1zSize = box1zSize / config.dx_SI;

        real_t box1xCenter = eval_parameters.getParameter<real_t>("box1_Center_x", real_t(0.05));
        box1xCenter = box1xCenter / config.dx_SI;
        real_t box1yCenter = eval_parameters.getParameter<real_t>("box1_Center_y", real_t(0.05));
        box1yCenter = box1yCenter / config.dx_SI;
        real_t box1zCenter = eval_parameters.getParameter<real_t>("box1_Center_z", real_t(0.01));
        box1zCenter = box1zCenter / config.dx_SI;

        const bool box2 = eval_parameters.getParameter<bool>("box2", false);

        real_t box2xSize = eval_parameters.getParameter<real_t>("box2_Size_x", real_t(0.02));
        box2xSize = box2xSize / config.dx_SI;
        real_t box2ySize = eval_parameters.getParameter<real_t>("box2_Size_y", real_t(0.02));
        box2ySize = box2ySize / config.dx_SI;
        real_t box2zSize = eval_parameters.getParameter<real_t>("box2_Size_z", real_t(0.02));
        box2zSize = box2zSize / config.dx_SI;

        real_t box2xCenter = eval_parameters.getParameter<real_t>("box2_Center_x", real_t(0.05));
        box2xCenter = box2xCenter / config.dx_SI;
        real_t box2yCenter = eval_parameters.getParameter<real_t>("box2_Center_y", real_t(0.10));
        box2yCenter = box2yCenter / config.dx_SI;
        real_t box2zCenter = eval_parameters.getParameter<real_t>("box2_Center_z", real_t(0.01));
        box2zCenter = box2zCenter / config.dx_SI;

        const bool box3 = eval_parameters.getParameter<bool>("box3", false);

        real_t box3xSize = eval_parameters.getParameter<real_t>("box3_Size_x", real_t(0.02));
        box3xSize = box3xSize / config.dx_SI;
        real_t box3ySize = eval_parameters.getParameter<real_t>("box3_Size_y", real_t(0.02));
        box3ySize = box3ySize / config.dx_SI;
        real_t box3zSize = eval_parameters.getParameter<real_t>("box3_Size_z", real_t(0.02));
        box3zSize = box3zSize / config.dx_SI;

        real_t box3xCenter = eval_parameters.getParameter<real_t>("box3_Center_x", real_t(0.05));
        box3xCenter = box3xCenter / config.dx_SI;
        real_t box3yCenter = eval_parameters.getParameter<real_t>("box3_Center_y", real_t(0.15));
        box3yCenter = box3yCenter / config.dx_SI;
        real_t box3zCenter = eval_parameters.getParameter<real_t>("box3_Center_z", real_t(0.01));
        box3zCenter = box3zCenter / config.dx_SI;

        // VTK PARAMETERS
        auto vtk_parameters = env.config()->getBlock("VTK");

        const bool writePackedBed = vtk_parameters.getParameter<bool>("write_packed_bed", false);
        const uint_t frequencyPackedBed = vtk_parameters.getParameter<uint_t>("frequency_packed_bed", uint_c(500));

        const bool writeSpheres = vtk_parameters.getParameter<bool>("write_spheres", false);
        const uint_t frequencySpheres = vtk_parameters.getParameter<uint_t>("frequency_spheres", uint_c(200));

        const bool writeFlagField = vtk_parameters.getParameter<bool>("write_flag_field", false);
        const uint_t frequencyFlagField = vtk_parameters.getParameter<uint_t>("frequency_flag_field", uint_c(200000));

        const bool writeFluidField = vtk_parameters.getParameter<bool>("write_fluid_field", false);
        const uint_t frequencyFluidField = vtk_parameters.getParameter<uint_t>("frequency_fluid_field", uint_c(50000));
        const real_t samplingResolutionFluidField = vtk_parameters.getParameter<real_t>("sampling_resolution_fluid_field", real_c(2.0));

        const bool writeOmega = vtk_parameters.getParameter<bool>("write_omega", false);
        const uint_t frequencyOmega = vtk_parameters.getParameter<uint_t>("frequency_omega", uint_c(50000));
        const real_t samplingResolutionOmega = vtk_parameters.getParameter<real_t>("sampling_resolution_omega", real_c(2.0));


        ///////////////////////////
        // BLOCK STRUCTURE SETUP //
        ///////////////////////////

        const uint_t xCells = config.xlength / xBlocks;
        const uint_t yCells = config.ylength / yBlocks;
        const uint_t zCells = config.zlength / zBlocks;

        // ASSERTS
        WALBERLA_CHECK_EQUAL(config.xlength, xCells * xBlocks);
        WALBERLA_CHECK_EQUAL(config.ylength, yCells * yBlocks);
        WALBERLA_CHECK_EQUAL(config.zlength, zCells * zBlocks);

        const real_t dx = real_t(1);
        const real_t overlap = real_c(1.5) * dx;

        auto blocks = blockforest::createUniformBlockGrid(xBlocks, yBlocks, zBlocks, xCells, yCells, zCells, dx, true,
                                                          xPeriodic, false, zPeriodic);

        // LOCAL VOLUME FOR SOLID FRACTION CALCULATION AND GRANULAR TEMPERATURE LOCALLY
        const AABB probe1(
                (blocks->getDomain().xMin() < (box1xCenter - box1xSize / 2.0)) ? (box1xCenter - box1xSize / 2.0) : blocks->getDomain().xMin(),
                (blocks->getDomain().yMin() < (box1yCenter - box1ySize / 2.0)) ? (box1yCenter - box1ySize / 2.0) : blocks->getDomain().yMin(),
                (blocks->getDomain().zMin() < (box1zCenter - box1zSize / 2.0)) ? (box1zCenter - box1zSize / 2.0) : blocks->getDomain().zMin(),
                (blocks->getDomain().xMax() > (box1xCenter + box1xSize / 2.0)) ? (box1xCenter + box1xSize / 2.0) : blocks->getDomain().xMax(),
                (blocks->getDomain().yMax() > (box1yCenter + box1ySize / 2.0)) ? (box1yCenter + box1ySize / 2.0) : blocks->getDomain().yMax(),
                (blocks->getDomain().zMax() > (box1zCenter + box1zSize / 2.0)) ? (box1zCenter + box1zSize / 2.0) : blocks->getDomain().zMax());

        const AABB probe2(
                (blocks->getDomain().xMin() < (box2xCenter - box2xSize / 2.0)) ? (box2xCenter - box2xSize / 2.0) : blocks->getDomain().xMin(),
                (blocks->getDomain().yMin() < (box2yCenter - box2ySize / 2.0)) ? (box2yCenter - box2ySize / 2.0) : blocks->getDomain().yMin(),
                (blocks->getDomain().zMin() < (box2zCenter - box2zSize / 2.0)) ? (box2zCenter - box2zSize / 2.0) : blocks->getDomain().zMin(),
                (blocks->getDomain().xMax() > (box2xCenter + box2xSize / 2.0)) ? (box2xCenter + box2xSize / 2.0) : blocks->getDomain().xMax(),
                (blocks->getDomain().yMax() > (box2yCenter + box2ySize / 2.0)) ? (box2yCenter + box2ySize / 2.0) : blocks->getDomain().yMax(),
                (blocks->getDomain().zMax() > (box2zCenter + box2zSize / 2.0)) ? (box2zCenter + box2zSize / 2.0) : blocks->getDomain().zMax());

        const AABB probe3(
                (blocks->getDomain().xMin() < (box3xCenter - box3xSize / 2.0)) ? (box3xCenter - box3xSize / 2.0) : blocks->getDomain().xMin(),
                (blocks->getDomain().yMin() < (box3yCenter - box3ySize / 2.0)) ? (box3yCenter - box3ySize / 2.0) : blocks->getDomain().yMin(),
                (blocks->getDomain().zMin() < (box3zCenter - box3zSize / 2.0)) ? (box3zCenter - box3zSize / 2.0) : blocks->getDomain().zMin(),
                (blocks->getDomain().xMax() > (box3xCenter + box3xSize / 2.0)) ? (box3xCenter + box3xSize / 2.0) : blocks->getDomain().xMax(),
                (blocks->getDomain().yMax() > (box3yCenter + box3ySize / 2.0)) ? (box3yCenter + box3ySize / 2.0) : blocks->getDomain().yMax(),
                (blocks->getDomain().zMax() > (box3zCenter + box3zSize / 2.0)) ? (box3zCenter + box3zSize / 2.0) : blocks->getDomain().zMax());

        // Transformation of AABB to CellInterval for evaluation functions
        if (box1) {
            config.probesAABB.push_back(probe1);
            config.probesCell.push_back(blocks->getCellBBFromAABB(probe1));
        }
        if (box2) {
            config.probesAABB.push_back(probe2);
            config.probesCell.push_back(blocks->getCellBBFromAABB(probe2));
        }
        if (box3) {
            config.probesAABB.push_back(probe3);
            config.probesCell.push_back(blocks->getCellBBFromAABB(probe3));
        }

        const real_t yLowerPlane = config.offsetBot <= 0.0 ? 0.5 : config.offsetBot - 0.5;
        const real_t yUpperPlane = config.offsetTop <= 0.0 ? real_c(config.ylength) - 0.5 : real_c(config.ylength) - config.offsetTop + 0.5;

        const AABB lowerPlaneAABB(blocks->getDomain().xMin(), yLowerPlane, blocks->getDomain().zMin(),
                                  blocks->getDomain().xMax(), yLowerPlane, blocks->getDomain().zMax());
        const AABB upperPlaneAABB(blocks->getDomain().xMin(), yUpperPlane, blocks->getDomain().zMin(),
                                  blocks->getDomain().xMax(), yUpperPlane, blocks->getDomain().zMax());

        const CellInterval lowerCellPlane = blocks->getCellBBFromAABB(lowerPlaneAABB);
        const CellInterval upperCellPlane = blocks->getCellBBFromAABB(upperPlaneAABB);

        /////////////////
        // PE COUPLING //
        /////////////////

        shared_ptr<pe::BodyStorage> globalBodyStorage = make_shared<pe::BodyStorage>();
        pe::SetBodyTypeIDs<BodyTypeTuple>::execute();

        auto bodyStorageID = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "pe Body Storage");

        const auto materialA = pe::createMaterial("myMatA", config.densityRatioA,
                                                  restitutionParticlesA, cdsA, cdfA, poisson, youngsModulus, stiffnessNparticlesA,
                                                  dampingNparticlesA, dampingTparticlesA);
        const auto materialB = pe::createMaterial("myMatB", config.densityRatioB,
                                                  restitutionParticlesB, cdsB, cdfB, poisson, youngsModulus, stiffnessNparticlesB,
                                                  dampingNparticlesA, dampingTparticlesB);
        const auto materialWall = pe::createMaterial("myMatWall", 0.5 * (config.densityRatioA + config.densityRatioB),
                                                     restitutionWall, cdsWall, cdfWall, poisson, youngsModulus, stiffnessNwall, dampingNwall,
                                                     dampingTwall);

        // fixed geometry
        if (!zPeriodic) {
            pe::createPlane(*globalBodyStorage, 0, Vector3<real_t>(0, 0, 1), Vector3<real_t>(0, 0, 0),
                            materialWall);                                                           // front
            pe::createPlane(*globalBodyStorage, 0, Vector3<real_t>(0, 0, -1), Vector3<real_t>(0, 0, real_c(config.zlength)),
                            materialWall);                                     // back
        }

        if (!xPeriodic) {
            pe::createPlane(*globalBodyStorage, 0, Vector3<real_t>(1, 0, 0), Vector3<real_t>(0, 0, 0),
                            materialWall);                                                           // left
            pe::createPlane(*globalBodyStorage, 0, Vector3<real_t>(-1, 0, 0), Vector3<real_t>(real_c(config.xlength), 0, 0),
                            materialWall);                                     // right
        }

        pe::createPlane(*globalBodyStorage, 0, Vector3<real_t>(0, 1, 0), Vector3<real_t>(0, config.offsetBot, 0),
                        materialWall);                            // bottom
        pe::createPlane(*globalBodyStorage, 0, Vector3<real_t>(0, -1, 0), Vector3<real_t>(0, real_c(config.ylength) - config.offsetTop, 0),
                        materialWall);                            // top

        vtk::writeDomainDecomposition(blocks, "domain_decomposition");

        ///////////////////
        // PE SIMULATION //
        ///////////////////

        //create empty pdf field
        BlockDataID pdfFieldID;

        // create omega  field for turbulence model
        BlockDataID omegaFieldID = field::addToStorage<ScalarField_T>(blocks, "omega field", omega, field::fzyx, FieldGhostLayers);

        // create appropriate lattice model
        LatticeModel_T latticeModel = LatticeModelCreator<collisionModel_T>::makeLatticeModel(omegaFieldID, omega);

        uint_t numberCreatedParticlesA = uint_c(0);
        uint_t numberCreatedParticlesB = uint_c(0);

        auto ccdID = blocks->addBlockData(pe::ccd::createHashGridsDataHandling(globalBodyStorage, bodyStorageID), "CCD");

        // set up collision response, here DEM solver
        auto fcdID = blocks->addBlockData(pe::fcd::createGenericFCDDataHandling<BodyTypeTuple, pe::fcd::AnalyticCollideFunctor>(), "FCD");

        std::unique_ptr<pe::cr::ICR> cr;
        if (useDEM) {
            cr = std::make_unique<pe::cr::DEM>(globalBodyStorage, blocks->getBlockStoragePointer(), bodyStorageID, ccdID, fcdID);
            WALBERLA_LOG_INFO_ON_ROOT("Using DEM!");

        } else {
            cr = std::make_unique<pe::cr::HCSITS>(globalBodyStorage, blocks->getBlockStoragePointer(), bodyStorageID, ccdID, fcdID);
            pe::cr::HCSITS *hcsits = static_cast<pe::cr::HCSITS *>(cr.get());
            hcsits->setMaxIterations(HCSITSMaxIterations);
            hcsits->setRelaxationParameter(HCSITSRelaxationParameter);
            WALBERLA_LOG_INFO_ON_ROOT("Using HCSITS!");
        }

        // set up synchronization procedure
        std::function<void(void)> syncCall = std::bind(pe::syncNextNeighbors<BodyTypeTuple>, std::ref(blocks->getBlockForest()), bodyStorageID,
                                                       static_cast<WcTimingTree *>(NULL), overlap, false);

        if (initializeFromCheckPointFile) {

            WALBERLA_LOG_INFO_ON_ROOT("Initializing data from checkpoint file!");

            WALBERLA_LOG_RESULT_ON_ROOT("Reading pdf from checkpoint file");

            shared_ptr<lbm::internal::PdfFieldHandling<LatticeModel_T> > pdfDataHandling =
                    make_shared<lbm::internal::PdfFieldHandling<LatticeModel_T> >(blocks, latticeModel, false, Vector3<real_t>(0), real_t(1),
                                                                                  FieldGhostLayers, field::fzyx);
            // add pdf field
            pdfFieldID = (blocks->getBlockStorage()).loadBlockData(checkpointPathPdf, pdfDataHandling, "pdf field (fzyx)");

            WALBERLA_LOG_RESULT_ON_ROOT("Reading body from checkpoint file");

            bodyStorageID = blocks->loadBlockData(checkpointPathBodies, pe::createStorageDataHandling<BodyTypeTuple>());

        } else {

            std::mt19937 rnd(MPIManager::instance()->rank());

            std::default_random_engine gen(MPIManager::instance()->rank());
            std::normal_distribution<real_t> distributionA(radiusA, radiusStandardDeviationA);
            std::normal_distribution<real_t> distributionB(radiusB, radiusStandardDeviationB);

            pe::SphereID sphere;

            real_t ratioAB = real_c(config.nrParticlesA) / real_c(config.nrParticlesB);

            real_t xInterval = initialSphereDistance;
            real_t yInterval = initialSphereDistance;
            real_t zInterval = initialSphereDistance;

            const uint_t rows = uint_c(real_c(config.xlength) / xInterval);
            const uint_t columns = uint_c(real_c(config.zlength) / zInterval);

            xInterval = real_c(config.xlength) / real_c(rows);
            zInterval = real_c(config.zlength) / real_c(columns);

            if (twoSpecies) {
                for (unsigned int i = 0; i < (config.nrParticlesA + config.nrParticlesB); ++i) {
                    if (real_c(numberCreatedParticlesA + 1) < ratioAB * real_c(numberCreatedParticlesB + 1)) {

                        sphere = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, i + 1,
                                                  Vector3<real_t>(0.5 * xInterval + xInterval * real_c((i % (rows * columns)) / columns),
                                                                  config.offsetBot + 0.5 * yInterval + yInterval * real_c(i / (rows * columns)),
                                                                  0.5 * zInterval + zInterval * real_c((i % columns))),
                                                  distributionA(gen), materialA);

                        if (sphere != NULL) {

                            ++numberCreatedParticlesA;

                            if (initialSphereVel > 0.0) {
                                sphere->setLinearVel(math::realRandom(-initialSphereVel, initialSphereVel, rnd),
                                                     math::realRandom(-initialSphereVel, initialSphereVel, rnd),
                                                     math::realRandom(-initialSphereVel, initialSphereVel, rnd));
                            }
                        }
                    } else {

                        sphere = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, i + 1,
                                                  Vector3<real_t>(0.5 * xInterval + xInterval * real_c((i % (rows * columns)) / columns),
                                                                  config.offsetBot + 0.5 * yInterval + yInterval * real_c(i / (rows * columns)),
                                                                  0.5 * zInterval + zInterval * real_c((i % columns))),
                                                  distributionB(gen), materialB);

                        if (sphere != NULL) {

                            ++numberCreatedParticlesB;

                            if (initialSphereVel > 0.0) {
                                sphere->setLinearVel(math::realRandom(-initialSphereVel, initialSphereVel, rnd),
                                                     math::realRandom(-initialSphereVel, initialSphereVel, rnd),
                                                     math::realRandom(-initialSphereVel, initialSphereVel, rnd));
                            }
                        }

                    }
                }
            } else {
                for (unsigned int i = 0; i < (config.nrParticlesA); ++i) {

                    sphere = pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, i + 1,
                                              Vector3<real_t>(0.5 * xInterval + xInterval * real_c((i % (rows * columns)) / columns),
                                                              config.offsetBot + 0.5 * yInterval + yInterval * real_c(i / (rows * columns)),
                                                              0.5 * zInterval + zInterval * real_c((i % columns))),
                                              distributionA(gen), materialA);

                    if (sphere != NULL) {

                        ++numberCreatedParticlesA;

                        if (initialSphereVel > 0.0) {
                            sphere->setLinearVel(math::realRandom(-initialSphereVel, initialSphereVel, rnd),
                                                 math::realRandom(-initialSphereVel, initialSphereVel, rnd),
                                                 math::realRandom(-initialSphereVel, initialSphereVel, rnd));
                        }
                    }
                }
            }

            WALBERLA_MPI_SECTION() {
                mpi::allReduceInplace(numberCreatedParticlesA, mpi::SUM);
                mpi::allReduceInplace(numberCreatedParticlesB, mpi::SUM);
            }

            WALBERLA_ROOT_SECTION() {
                std::cout << "Number of created particles A is:     " << numberCreatedParticlesA << std::endl;
                std::cout << "Number of created particles B is:     " << numberCreatedParticlesB << std::endl << std::endl;
            }

            ///////////////////
            // PACKED BED ////
            //////////////////

            if (createPackedBed) {

                // spheres_init
                auto sphereVtkOutputInit = make_shared<pe::SphereVtkOutput>(bodyStorageID, blocks->getBlockStorage());
                auto sphereVtkWriterInit = vtk::createVTKOutput_PointData(sphereVtkOutputInit, "spheres_init");

                sphereVtkWriterInit->write();

                cr->setGlobalLinearAcceleration(pe::Vec3(0.0, -(config.gravity), 0.0));

                real_t avgVel = real_c(100.0);
                real_t maxVel = real_c(100.0);
                real_t currentSettleVelocity = real_c(100.0);
                uint_t itPE = uint_c(0);
                real_t velocityConversion = config.dx_SI / config.dt_SI;

                while ((currentSettleVelocity > settleVelocity) || (itPE < 4999)) {

                    cr->timestep(peAccelerator / real_c(config.peSubCycles));
                    syncCall();
                    ++itPE;

                    if (writePackedBed && itPE % frequencyPackedBed == 0) {
                        sphereVtkWriterInit->write();
                    }

                    if ((itPE) % settleCheckFrequency == 0) {

                        avgVel =
                                FBfunc::getAvgVel(blocks, bodyStorageID, numberCreatedParticlesA + numberCreatedParticlesB) * velocityConversion;
                        maxVel = FBfunc::getMaxVel(blocks, bodyStorageID) * velocityConversion;

                        if (checkMaximumParticleVelocity) {
                            currentSettleVelocity = maxVel;
                        } else {
                            currentSettleVelocity = avgVel;
                        }

                        WALBERLA_ROOT_SECTION() {
                            std::cout << "Time step " << itPE << std::endl;
                            std::cout << "Average Velocity:   " << avgVel << " m/s" << std::endl;
                            std::cout << "Maximum Velocity:   " << maxVel << " m/s" << std::endl << std::endl;
                        }
                    }
                }


                cr->setGlobalLinearAcceleration(pe::Vec3(0.0, 0.0, 0.0));
                FBfunc::restSphere(blocks, bodyStorageID);

                avgVel = FBfunc::getAvgVel(blocks, bodyStorageID, (numberCreatedParticlesA + numberCreatedParticlesB)) * velocityConversion;
                avgVel = FBfunc::getAvgVel(blocks, bodyStorageID, (numberCreatedParticlesA + numberCreatedParticlesB)) * velocityConversion;
                maxVel = FBfunc::getMaxVel(blocks, bodyStorageID) * velocityConversion;

                WALBERLA_ROOT_SECTION() {
                    std::cout << "||||||||||||| END OF PE SIM ||||||||||||||" << std::endl;
                    std::cout << "Average Velocity:   " << avgVel << " m/s" << std::endl;
                    std::cout << "Maximum Velocity:   " << maxVel << " m/s" << std::endl << std::endl;
                }
            }


            ////////////////////
            // PACKED BED END //
            ////////////////////


            WALBERLA_ROOT_SECTION() {
                std::cout << "\nSetup Information\n" << std::endl;
                std::cout << "Physical time step is:                " << config.dt_SI << " s" << std::endl;
                std::cout << "Physical cell size is:                " << config.dx_SI << " m" << std::endl;
                std::cout << "Particle Reynoldsnumber(LBM) is:      " << particleReynoldsNumber_LBM << std::endl;
                std::cout << "Froudenumber(LBM) is:                 " << FroudeNumber_LBM << std::endl;
                std::cout << "Minimum Fluidization Re number is:    " << REpmf << std::endl;
                std::cout << "LBM velocity is:                      " << config.velocity << std::endl;
                std::cout << "LBM viscosity is:                     " << config.kinViscosity << std::endl;
                std::cout << "Omega is:                             " << omega << std::endl << std::endl;

                std::ofstream output;
                output.open("Setup.txt");
                output << "Domain size in X direction:           " << config.xlength << std::endl;
                output << "Domain size in Y direction:           " << config.ylength << std::endl;
                output << "Domain size in Z direction:           " << config.zlength << std::endl;
                output << "Number of processes:                  " << processes << std::endl;
                output << "Physical time step is:                " << config.dt_SI << " s" << std::endl;
                output << "Physical cell size is:                " << config.dx_SI << " m" << std::endl;
                output << "Number of particles A is:             " << config.nrParticlesA << std::endl;
                if (twoSpecies) {
                    output << "Number of particles B is:             " << config.nrParticlesA << std::endl;
                }
                if (!initializeFromCheckPointFile) {
                    output << "Number of created particles A is:     " << numberCreatedParticlesA << std::endl;
                    if (twoSpecies) {
                        output << "Number of created particles B is:     " << numberCreatedParticlesB << std::endl;
                    }

                }
                output << "Cells per diameter :                  " << config.diameterA << std::endl;
                output << "Fluid density is:                     " << config.densityFluid_SI << " kg/m3" << std::endl;
                output << "Fluid viscosity is:                   " << dynViscosity_SI << std::endl;
                output << "Particle density is:                  " << densitySolidA_SI << " kg/m3" << std::endl;
                output << "Inflow velocity is:                   " << velocity_SI << " m/s" << std::endl;
                output << "Normal Stiffness A is:                " << stiffnessNparticlesA << std::endl;
                if (twoSpecies) {
                    output << "Normal Stiffness B is:                " << stiffnessNparticlesB << std::endl;
                }
                output << "Normal Stiffness Wall is:             " << stiffnessNwall << std::endl;
                output << "Normal Damping A is:                  " << dampingNparticlesA << std::endl;
                if (twoSpecies) {
                    output << "Normal Damping B is:                    " << dampingNparticlesB << std::endl;
                }
                output << "Normal Damping Wall is:               " << dampingNwall << std::endl;
                output << "Tangential Damping A is:              " << dampingTparticlesA << std::endl;
                if (twoSpecies) {
                    output << "Tangential Damping B is:                " << dampingTparticlesB << std::endl;
                }
                output << "Tangential Damping Wall is:           " << dampingTwall << std::endl;
                output << "Particle Reynoldsnumber(SI) is:       " << particleReynoldsNumber_SI << std::endl;
                output << "Particle Reynoldsnumber(LBM) is:      " << particleReynoldsNumber_LBM << std::endl;
                output << "Froudenumber(SI) is:                  " << FroudeNumber_SI << std::endl;
                output << "Froudenumber(LBM) is:                 " << FroudeNumber_LBM << std::endl;
                output << "Arrhenius number is:                  " << ArrheniusNumber_SI << std::endl;
                output << "Galilei number is:                    " << galileiNumber << std::endl;
                output << "Minimum Fluidization Re number is:    " << REpmf << std::endl;
                output << "Minimum Fluidization Fr number is:    " << UmfFroudenumber << std::endl;
                output << "LBM velocity is:                      " << config.velocity << std::endl;
                output << "LBM gravity is:                       " << config.gravity << std::endl;
                output << "LBM viscosity is:                     " << config.kinViscosity << std::endl;
                output << "Pressure due to gravity is:           " << gravityPressure << " Pa" << std::endl;
                output << "Omega is:                             " << omega << std::endl << std::endl;
                output.close();
            }


            // add pdf field
            pdfFieldID = lbm::addPdfFieldToStorage(blocks, "pdf field (fzyx)", latticeModel,
                                                   Vector3<real_t>(real_t(0),
                                                                   spotInflow ? real_t(0.0) : initialRampingVelocityFactor * config.velocity,
                                                                   real_t(0)), real_t(1),
                                                   uint_t(1), field::fzyx);
        }

        // Communication scheme
        std::function<void()> commFunction;

        blockforest::communication::UniformBufferedScheme<Stencil_T> scheme(blocks);
        scheme.addPackInfo(make_shared<lbm::PdfFieldPackInfo<LatticeModel_T> >(pdfFieldID));
        commFunction = scheme;

        ////////////////////////
        // ADD DATA TO BLOCKS //
        ////////////////////////

        // add flag field
        BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>(blocks, "flag field");

        // add body field
        BlockDataID bodyFieldID = field::addToStorage<BodyField_T>(blocks, "body field", NULL, field::fzyx);

        std::function<void(void)> syncCall2 = std::bind(pe::syncNextNeighbors<BodyTypeTuple>, std::ref(blocks->getBlockForest()), bodyStorageID,
                                                        static_cast<WcTimingTree *>(NULL), overlap, false);
        syncCall2();

        for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) {
            pe::ccd::ICCD* ccd = blockIt->getData<pe::ccd::ICCD>(ccdID);
            ccd->reloadBodies();
        }

        int numberSpheres = 0;

        for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) {
            for (auto bodyIt = pe::ShadowBodyIterator::begin(*blockIt, bodyStorageID); bodyIt != pe::ShadowBodyIterator::end(); ++bodyIt) {
                if (bodyIt->isFixed() || !bodyIt->isFinite())
                    continue;

                ++numberSpheres;
            }
        }

        // Object for keeping track of time
        shared_ptr<lbm::TimeTracker> timeTrack = make_shared<lbm::TimeTracker>();

        if(initializeFromCheckPointFile){
            timestepsRamping = uint_t(0);
        }

        BlockDataID boundaryHandlingID = blocks->addStructuredBlockData<BoundaryHandling_T>(
                MyBoundaryHandling_T(flagFieldID, pdfFieldID, bodyFieldID, xPeriodic, zPeriodic, spotInflow, spotDiameter, real_c(config.xlength)*0.5,
                                     real_c(config.zlength)*0.5, config.velocity, initialRampingVelocityFactor, real_c(timestepsRamping), timeTrack), "boundary handling");

        // map pe bodies into the LBM simulation
        // moving bodies are handled by the momentum exchange method
        pe_coupling::mapMovingBodies<BoundaryHandling_T>(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, FBfunc::MO_Flag, pe_coupling::selectRegularBodies);

        /////////////////////
        // PRE TIME LOOP   //
        /////////////////////

        SweepTimeloop timeloopPRE(blocks->getBlockStorage(), 1000000);

        // add LBM communication function and boundary handling sweep
        timeloopPRE.add() << BeforeFunction(commFunction, "LBM Communication") << Sweep(BoundaryHandling_T::getBlockSweep(boundaryHandlingID), "Boundary Handling");
        // LBM stream collide sweep
        timeloopPRE.add() << Sweep(makeSharedSweep(lbm::makeCellwiseSweep<LatticeModel_T, FlagField_T>(pdfFieldID, flagFieldID, FBfunc::Fluid_Flag)), "LBM DEFAULT");

        // Velocity Check
        timeloopPRE.add() << Sweep(FBfunc::VelocityCheckMEM<LatticeModel_T, FlagField_T, BodyField_T>(blocks, bodyStorageID, bodyFieldID, pdfFieldID, flagFieldID, FBfunc::Fluid_Flag, config),
                                   "LBM Velocity Check");

        shared_ptr<FBfunc::PressureDropper> PressureDropper = make_shared<FBfunc::PressureDropper>(blocks, bodyStorageID, config, true);
        timeloopPRE.addFuncAfterTimeStep( SharedFunctor<FBfunc::PressureDropper>(PressureDropper), "Calculate Pressure Drop");

        // Reset Forces
        timeloopPRE.addFuncAfterTimeStep(pe_coupling::ForceTorqueOnBodiesResetter(blocks, bodyStorageID), "Reset Body Forces");

        WcTimingPool timeloopTimingPRE;

        ////////////////////
        // MAIN TIME LOOP //
        ////////////////////

        SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);

        ////////////////
        // VTK OUTPUT //
        ////////////////

        if (writeSpheres) {
            auto sphereVtkOutput = make_shared<pe::SphereVtkOutput>(bodyStorageID, blocks->getBlockStorage());
            auto sphereVTK = vtk::createVTKOutput_PointData(sphereVtkOutput, "spheres", frequencySpheres);
            timeloop.addFuncAfterTimeStep(vtk::writeFiles(sphereVTK), "VTK (sphere data)");
        }

        if (writeOmega) {
            auto omegaVTK = vtk::createVTKOutput_BlockData(blocks, "omega", frequencyOmega, 0);          // 0 => no fieldghostlayer
            omegaVTK->addCellDataWriter(make_shared<field::VTKWriter<ScalarField_T, float >>(omegaFieldID, "omega_field"));
            omegaVTK->setSamplingResolution(samplingResolutionOmega);
            timeloop.addFuncAfterTimeStep(vtk::writeFiles(omegaVTK), "VTK (omega field data)");
        }

        if (writeFlagField) {
            // flag field (ghost layers are also written)
            auto flagFieldVTK = vtk::createVTKOutput_BlockData(blocks, "flag_field", frequencyFlagField, FieldGhostLayers);
            flagFieldVTK->addCellDataWriter(make_shared<field::VTKWriter<FlagField_T> >(flagFieldID, "FlagField"));
            timeloop.addFuncAfterTimeStep(vtk::writeFiles(flagFieldVTK), "VTK (flag field data)");
        }

        if (writeFluidField) {
            // pdf field (ghost layers cannot be written because re-sampling/coarsening is applied)
            auto pdfFieldVTK = vtk::createVTKOutput_BlockData(blocks, "fluid_field", frequencyFluidField, 0, false);      // 0 => no fieldghostlayer
            pdfFieldVTK->setSamplingResolution(samplingResolutionFluidField);

            blockforest::communication::UniformBufferedScheme<stencil::D3Q27> pdfGhostLayerSync(blocks);
            pdfGhostLayerSync.addPackInfo(make_shared<field::communication::PackInfo<PdfField_T> >(pdfFieldID));
            pdfFieldVTK->addBeforeFunction(pdfGhostLayerSync);

            field::FlagFieldCellFilter<FlagField_T> fluidFilter(flagFieldID);
            fluidFilter.addFlag(FBfunc::Fluid_Flag);
            pdfFieldVTK->addCellInclusionFilter(fluidFilter);

            pdfFieldVTK->addCellDataWriter(make_shared<lbm::VelocityVTKWriter<LatticeModel_T, float> >(pdfFieldID, "VelocityFromPDF"));
            pdfFieldVTK->addCellDataWriter(make_shared<lbm::DensityVTKWriter<LatticeModel_T, float> >(pdfFieldID, "DensityFromPDF"));

            timeloop.addFuncAfterTimeStep(vtk::writeFiles(pdfFieldVTK), "VTK (fluid field data)");
        }

        ////////////////
        // VTK END    //
        ////////////////

        timeloop.addFuncAfterTimeStep(
                makeSharedFunctor( field::makeStabilityChecker< lbm::PdfField< LatticeModel_T >, FlagField_T >(
                        blocks, pdfFieldID, flagFieldID, FBfunc::Fluid_Flag, uint_t(1), false, true ) ), "LBM stability check" );

        // sweep for updating the pe body mapping into the LBM simulation
        timeloop.add() << Sweep(
                pe_coupling::BodyMapping<BoundaryHandling_T>(blocks, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID,
                                                             FBfunc::MO_Flag, FBfunc::FormerMO_Flag), "Body Mapping");

        // sweep for restoring PDFs in cells previously occupied by pe bodies
        typedef pe_coupling::EquilibriumReconstructor<LatticeModel_T, BoundaryHandling_T> Reconstructor_T;
        Reconstructor_T reconstructor(blocks, boundaryHandlingID, pdfFieldID, bodyFieldID);
        timeloop.add() << Sweep(pe_coupling::PDFReconstruction<LatticeModel_T, BoundaryHandling_T, Reconstructor_T>
                                        (blocks, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID, reconstructor,
                                         FBfunc::FormerMO_Flag, FBfunc::Fluid_Flag), "PDF Restore");

        if (turb) {
            // create flag field filter for turbulence model
            field::FlagFieldEvaluationFilter<FlagField_T> filter = field::FlagFieldEvaluationFilter<FlagField_T>(flagFieldID, FBfunc::Fluid_Flag);

            // Turbulence model
            timeloop.addFuncAfterTimeStep(lbm::SmagorinskyLES<LatticeModel_T, field::FlagFieldEvaluationFilter<FlagField_T> >
                                                  (blocks, filter, pdfFieldID, omegaFieldID, config.kinViscosity, smagorinskyConstant), "Turbulence Model");
        }


        if (enableCheckpointing) {
            timeloop.addFuncAfterTimeStep(FBfunc::CheckpointCreater(blocks->getBlockStorage(), pdfFieldID, bodyStorageID, config), "Checkpointing");
        }

        shared_ptr<pe_coupling::BodiesForceTorqueContainer> bodiesFTContainer1 = make_shared<pe_coupling::BodiesForceTorqueContainer>(blocks, bodyStorageID);
        std::function<void(void)> storeForceTorqueInCont1 = std::bind(&pe_coupling::BodiesForceTorqueContainer::store, bodiesFTContainer1);
        shared_ptr<pe_coupling::BodiesForceTorqueContainer> bodiesFTContainer2 = make_shared<pe_coupling::BodiesForceTorqueContainer>(blocks, bodyStorageID);
        std::function<void(void)> setForceTorqueOnBodiesFromCont2 = std::bind(&pe_coupling::BodiesForceTorqueContainer::setOnBodies, bodiesFTContainer2);
        shared_ptr<pe_coupling::ForceTorqueOnBodiesScaler> forceScaler = make_shared<pe_coupling::ForceTorqueOnBodiesScaler>(blocks, bodyStorageID, real_t(1));
        std::function<void(void)> setForceScalingFactorToHalf = std::bind(&pe_coupling::ForceTorqueOnBodiesScaler::resetScalingFactor,forceScaler,real_t(0.5));

        if( useForceAveraging ) {
            bodiesFTContainer2->store();
        }

        // add LBM communication function and boundary handling sweep
        timeloop.add() << BeforeFunction(commFunction, "LBM Communication") << Sweep(BoundaryHandling_T::getBlockSweep(boundaryHandlingID), "Boundary Handling");
        // LBM stream collide sweep
        timeloop.add() << Sweep(makeSharedSweep(lbm::makeCellwiseSweep<LatticeModel_T, FlagField_T>(pdfFieldID, flagFieldID, FBfunc::Fluid_Flag)), "LBM DEFAULT");
        //timeloop.add() << Sweep( lbm::SplitPureSweep< LatticeModel_T >( pdfFieldID ), "LBM SPLIT PURE" );

        // Velocity Check
        timeloop.add() << Sweep(FBfunc::VelocityCheckMEM<LatticeModel_T, FlagField_T, BodyField_T>
                                        (blocks, bodyStorageID, bodyFieldID, pdfFieldID, flagFieldID, FBfunc::Fluid_Flag, config), "LBM Velocity Check");

        if( useForceAveraging ) {

            // store force/torque from hydrodynamic interactions in container1
            timeloop.addFuncAfterTimeStep(storeForceTorqueInCont1, "Force Storing");

            // set force/torque from previous time step (in container2)
            timeloop.addFuncAfterTimeStep(setForceTorqueOnBodiesFromCont2, "Force setting");

            // average the force/torque by scaling it with factor 1/2 (except in first timestep, there it is 1, which it is initially)
            timeloop.addFuncAfterTimeStep( pe_coupling::ForceTorqueOnBodiesScaler(blocks, bodyStorageID, real_t(0.5)),  "Force averaging");
            timeloop.addFuncAfterTimeStep( setForceScalingFactorToHalf, "Force scaling adjustment" );

            // swap containers
            timeloop.addFuncAfterTimeStep( pe_coupling::BodyContainerSwapper( bodiesFTContainer1, bodiesFTContainer2 ), "Swap FT container" );

        }

        if (useLubrication) {
            timeloop.addFuncAfterTimeStep(pe_coupling::LubricationCorrection(blocks, globalBodyStorage, bodyStorageID, config.kinViscosity), "Lubrication Force");
        }

        if (calculatePressureDrop) {
            timeloop.addFuncAfterTimeStep(FBfunc::PressureDropper(blocks, bodyStorageID, config, false), "Calculate Pressure Drop");
        }

        if (calculateParticleFlux) {
            timeloop.addFuncAfterTimeStep(FBfunc::ParticleFlux(blocks, bodyStorageID, config), "Calculate Particle Flux");
        }

        if (calculateGranularTempGlobal) {
            timeloop.addFuncAfterTimeStep(FBfunc::GranularTemperature(blocks, bodyStorageID, config), "Calculate Granular Temperature");
        }

        if (calculateGranularTempLocal) {
            timeloop.addFuncAfterTimeStep(FBfunc::GranularTemperatureLocally(blocks, bodyStorageID, config), "Calculate Granular Temperature Locally");
        }

        if (calculatePressureDifference) {
            timeloop.addFuncAfterTimeStep(FBfunc::PressureByDensity<LatticeModel_T, FlagField_T>(blocks, pdfFieldID, flagFieldID, lowerCellPlane, upperCellPlane, FBfunc::Fluid_Flag, config),
                                          "Calculate Pressure Drop with Density");
        }

        timeloop.addFuncAfterTimeStep(FBfunc::ExternalForce(blocks, bodyStorageID, config), "Add External Forces");

        if (calculateSolidVolumeFraction) {
            timeloop.addFuncAfterTimeStep(FBfunc::FractionCalculator<FlagField_T>(blocks, flagFieldID, config), "Calculate Solid Volume Fraction");
        }

        if (calculateCenterOfMass) {
            timeloop.addFuncAfterTimeStep(FBfunc::CenterOfMass(blocks, bodyStorageID, config, calculateErgunPressure), "Calculate Center of Mass");
        }

        // advance pe rigid body simulation
        timeloop.addFuncAfterTimeStep( pe_coupling::TimeStep( blocks, bodyStorageID, *cr, syncCall, config.lbmSubCycles, config.peSubCycles ), "pe Time Step" );

        // Obstacle test
        timeloop.addFuncAfterTimeStep(FBfunc::ObstacleLocationCheck(blocks, bodyStorageID, blocks->getDomain(),
                                                                    (xPeriodic || zPeriodic) ? radiusA + real_c(2) * overlap : real_c(2) * overlap), "Obstacle Location Check");

        timeloop.addFuncAfterTimeStep(FBfunc::ObstacleLinearVelocityCheck(blocks, bodyStorageID, config), "Obstacle Velocity Check");

        WcTimingPool timeloopTiming;

        timeloop.addFuncAfterTimeStep(RemainingTimeLogger(timeloop.getNrOfTimeSteps(), 60), "Remaining Time Logger");

        timeloop.addFuncAfterTimeStep(FBfunc::TimingPoolLogger(timeloopTiming, timeloop, 1000), "Timing Pool Logger");

        timeloop.addFuncAfterTimeStep(lbm::PerformanceLogger<FlagField_T>(blocks, flagFieldID, FBfunc::Fluid_Flag, uint_c(1000)), "Performance Logger");

        timeloop.addFuncAfterTimeStep(SharedFunctor<lbm::TimeTracker>(timeTrack), "Time Tracker");


        //////////////////////////
        // EXECUTE SIMULATIONS  //
        //////////////////////////

        if (runPreliminaryTimeloop && !initializeFromCheckPointFile) {
            WALBERLA_ROOT_SECTION()
            {
                std::cout << std::endl << "||||||||||||| PRELIMINARY SIM ||||||||||||||" << std::endl;
            }

            real_t relativePressureDifference = real_t(1000);
            real_t previousPressure = real_t(0);

            // perform a single simulation step
            timeloopPRE.singleStep(timeloopTimingPRE);
            previousPressure = PressureDropper->getPressureDrop();

            while (relativePressureDifference > eps) {

                // perform a single simulation step
                timeloopPRE.singleStep(timeloopTimingPRE);

                if (timeloopPRE.getCurrentTimeStep() % uint_t(20) == uint_t(0)) {

                    relativePressureDifference = fabs(PressureDropper->getPressureDrop() - previousPressure) / previousPressure;

                    WALBERLA_LOG_INFO_ON_ROOT(
                            "Current pressure drop over whole bed: " << PressureDropper->getPressureDrop()
                                                                     << " [Pa] after "
                                                                     << timeloopPRE.getCurrentTimeStep()
                                                                     << " timesteps!");
                    WALBERLA_LOG_INFO_ON_ROOT(
                            "Relativ change: " << relativePressureDifference );

                    previousPressure = PressureDropper->getPressureDrop();

                    WALBERLA_MPI_SECTION() {
                        mpi::broadcastObject(relativePressureDifference);
                    }
                }
            }

            WALBERLA_LOG_INFO_ON_ROOT(
                    "Pressure drop has converged to at: " << PressureDropper->getPressureDrop() << "[Pa]");
        }

        WALBERLA_ROOT_SECTION() {
            std::cout << std::endl << "||||||||||||| MAIN SIM ||||||||||||||" << std::endl;
        }

        timeloop.run(timeloopTiming);
        timeloopTiming.logResultOnRoot();

        return 0;
    }

//////////
// MAIN //
//////////

    int main(int argc, char **argv) {

        using namespace walberla;

        Environment env(argc, argv);

        auto model_parameters = env.config()->getBlock("Model");

        if (argc < 2) {
            WALBERLA_ROOT_SECTION() {
                std::cout << "Usage: " << argv[0] << " path-to-configuration-file " << std::endl;
            }
            return EXIT_SUCCESS;
        }


        if (model_parameters.getParameter<bool>("turbulence_model", false)){

            return runSimulation<lbm::collision_model::SRTField<ScalarField_T>>(env);

        } else {

            return runSimulation<lbm::collision_model::TRT>(env);
        }
    }

}

int main( int argc, char ** argv )
{
    return fluidizedBed::main( argc, argv );
}

