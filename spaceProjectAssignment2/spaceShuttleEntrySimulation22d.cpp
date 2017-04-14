/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>

#include <spaceProjectAssignment2/inputOutputDirectories.h>

namespace tudat
{


///                                                                 ////
/// Barebones class for aerodynamic guidance of the space shuttle   ////
///                                                                 ////
class SpaceShuttleAerodynamicGuidance: public aerodynamics::AerodynamicGuidance
{

public:
     //! Constructor
    SpaceShuttleAerodynamicGuidance( const simulation_setup::NamedBodyMap& bodyMap, const std::string vehicleName )

    {
        shuttleflightConditions_ =
             bodyMap.at( vehicleName )->getFlightConditions( );
        coefficientInterface_
             = bodyMap.at( vehicleName )->getAerodynamicCoefficientInterface( );
        bodyMap_ = bodyMap;

    }
private:

    boost::shared_ptr< aerodynamics::FlightConditions > shuttleflightConditions_;
    //! The aerodynamic angles are to be computed here
    void updateGuidance( const double time )
    {
        // Define the mach number, altitude, angle of attack and sideslip angle for the all assignment
        double Mach = shuttleflightConditions_->getCurrentMachNumber( );
        double Altitude = shuttleflightConditions_->getCurrentAltitude( );
        currentAngleOfSideslip_ = 0.0;

        if(  Mach > 12.0  || Altitude > 150000.0)
        {
            currentAngleOfAttack_ = 40.0  / 180.0 * mathematical_constants::PI;
        }
        else if( Mach < 6.0 )
        {
            currentAngleOfAttack_ = 10.0 / 180.0 * mathematical_constants::PI;
        }
        else
        {
            currentAngleOfAttack_ = ( 10.0 + (Mach - 6.0) * 5.0 ) / 180.0 * mathematical_constants::PI;

        }


         std::vector< double > currentAerodynamicParametersInput_;
         currentAerodynamicParametersInput_.push_back( currentAngleOfAttack_ );
         currentAerodynamicParametersInput_.push_back( shuttleflightConditions_->getCurrentMachNumber( ) );

         Eigen::Vector3d currentAerodynamicParameters = coefficientInterface_->getCurrentForceCoefficients( );


         double currentLiftCoefficient = currentAerodynamicParameters[2];

         double currentLattitude = shuttleflightConditions_->getAerodynamicAngleCalculator( )
                 ->getAerodynamicAngle( reference_frames::latitude_angle );

         double currentLongitude = shuttleflightConditions_->getAerodynamicAngleCalculator( )
                 ->getAerodynamicAngle( reference_frames::longitude_angle );

         double currentFlightPathAngle = shuttleflightConditions_->getAerodynamicAngleCalculator( )
                 ->getAerodynamicAngle( reference_frames::flight_path_angle );

         double currentHeadingAngle = shuttleflightConditions_->getAerodynamicAngleCalculator( )
                 ->getAerodynamicAngle( reference_frames::heading_angle );


         double currentDensity = shuttleflightConditions_->getCurrentDensity( );
         double currentAirspeed = shuttleflightConditions_->getCurrentAirspeed( );
         double omegaE = (2.0 * mathematical_constants::PI)/(physical_constants::SIDEREAL_DAY);

         double Area = 2690.0 * 0.3048 * 0.3048;

         double GravitationalConstantEarth = bodyMap_["Earth"]->getGravityFieldModel()->getGravitationalParameter();


         Eigen::Vector6d States = bodyMap_["STS"]->getState();
         double RadialDistance = States.segment(0,3).norm();
         double MassVehicle = 165000 * 0.45;


         // Create position vectors
         Eigen::Vector3d positionOfBodySubjectToAcceleration;
         positionOfBodySubjectToAcceleration << RadialDistance, 0.0, 0.0;

         const Eigen::Vector3d positionOfBodyExertingAcceleration = Eigen::Vector3d::Zero( );

         // Determine acceleration vectors.
         const Eigen::Vector3d gravitationalAcceleration =
                 tudat::gravitation::computeGravitationalAcceleration(
                         positionOfBodySubjectToAcceleration, GravitationalConstantEarth,
                         positionOfBodyExertingAcceleration );

         double g = gravitationalAcceleration.norm();
         // Determine the lift L= Clv^2 rho S/2
         double Lift = (pow(currentAirspeed,2)*currentDensity*currentLiftCoefficient*Area)/2.0;

         double flightPathAngleChangePart2 = -(g-(pow(currentAirspeed,2)/RadialDistance))*
                 cos(currentFlightPathAngle) + 2.0*omegaE*currentAirspeed *cos(currentLattitude)*sin(currentHeadingAngle);

         double flightPathAngleChangeIdeal = -1*(MassVehicle/Lift)*(flightPathAngleChangePart2);

         // Ideally change in flight path angle is 0 but if not minimize absolut value of flight path angle rate


            if (flightPathAngleChangeIdeal <= 1.0 && flightPathAngleChangeIdeal >= -1.0 )
            {
                currentBankAngle_ = acos(flightPathAngleChangeIdeal);
            }
            else if(flightPathAngleChangePart2 > 0.0)
            {
              currentBankAngle_ =   mathematical_constants::PI;
            }


         else if(flightPathAngleChangePart2 < 0.0)
         {
            if (flightPathAngleChangeIdeal <= 1.0 && flightPathAngleChangeIdeal >= -1.0 )
            {
                currentBankAngle_ = acos(flightPathAngleChangeIdeal);
            }
            else
            {
              currentBankAngle_ =  0.0;
            }

         }


     }


private:

    boost::shared_ptr< aerodynamics::FlightConditions > shuttleFlightConditions_;
    boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > coefficientInterface_;
    simulation_setup::NamedBodyMap bodyMap_;


};


}

//! Execute propagation of orbits of STS during entry.
int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    using namespace tudat::ephemerides;
    using namespace tudat::interpolators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::spice_interface;
    using namespace tudat::simulation_setup;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::propagators;
    using namespace tudat::aerodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::input_output;
    using namespace tudat;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );

    // Set simulation start and end epoch.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = 2.0*24.0*3600.0;


    // Step size of 10 s, changeit if needed
    const double fixedStepSize = 0.1;


//    // Define propagation termination conditions (stop after 2 days).
//    boost::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
//            boost::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define simulation body settings.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth" } );
    bodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< simulation_setup::ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< AtmosphereSettings >( nrlmsise00 );

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

    // Create vehicle objects
    bodyMap[ "STS" ] = boost::make_shared< simulation_setup::Body >( );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Create vehicle  coefficients
    std::map< int, std::string > forceCoefficientFiles;
    forceCoefficientFiles[ 0 ] =
             tudat_applications::getStsInputPath( ) + "STS_CD.dat";
    forceCoefficientFiles[ 2 ] =
            tudat_applications::getStsInputPath( ) + "STS_CL.dat";

    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            simulation_setup::readTabulatedAerodynamicCoefficientsFromFiles(
                forceCoefficientFiles, 2690.0 * 0.3048 * 0.3048,
                boost::assign::list_of( aerodynamics::angle_of_attack_dependent )(  aerodynamics::mach_number_dependent ),
                true, true );

    bodyMap[ "STS" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "STS" ) );
    bodyMap[ "STS" ]->setConstantBodyMass( 165000 * 0.45 );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > AccelerationsOfSts;
    AccelerationsOfSts[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );
    AccelerationsOfSts[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );

    accelerationMap[  "STS" ] = AccelerationsOfSts;
    bodiesToPropagate.push_back( "STS" );
    centralBodies.push_back( "Earth" );


  //  basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(bodyMap, AccelerationMap, bodiesToPropagate, centralBodies );

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );
    boost::shared_ptr<aerodynamics::AerodynamicGuidance> aerodynamicGuidanceSettings =
            boost::make_shared< SpaceShuttleAerodynamicGuidance >(
        bodyMap, "STS" );
    setGuidanceAnglesFunctions( aerodynamicGuidanceSettings, bodyMap.at( "STS" ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             DEFINE INITIAL STATE                   ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set spherical elements for STS.
    Eigen::Vector6d stsSphericalEntryState;
    stsSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
            spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
    stsSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) = 0.3;
    stsSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex ) = 1.2;
    stsSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 7.45E3;
    stsSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
            -1.2 * mathematical_constants::PI / 180.0;
    stsSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex ) = 0.6;

    // Convert sts state from spherical elements to Cartesian elements.
    Eigen::Vector6d systemInitialState = transformStateToGlobalFrame(
                convertSphericalOrbitalToCartesianState(
                                stsSphericalEntryState ),
                simulationStartEpoch, bodyMap.at( "Earth" )->getRotationalEphemeris( ) );


    // Define list of dependent variables to save.
        std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;

        dependentVariablesList.push_back(boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(aerodynamic, "STS", "Earth", 0 ) );

        // Saving the latitude and longitude
        //Latitude Angle
        dependentVariablesList.push_back(boost::make_shared< BodyAerodynamicAngleVariableSaveSettings>("STS", reference_frames::latitude_angle ));

        //longitude Angle
        dependentVariablesList.push_back(boost::make_shared< BodyAerodynamicAngleVariableSaveSettings>("STS",reference_frames::longitude_angle));
        //Others
        dependentVariablesList.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >( mach_number_dependent_variable, "STS" ) );
        dependentVariablesList.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >( airspeed_dependent_variable, "STS" ) );
        dependentVariablesList.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >( local_density_dependent_variable, "STS" ) );
        dependentVariablesList.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >(
                        altitude_dependent_variable, "STS", "Earth" ) );

       /* dependentVariablesList.push_back(
                    boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        aerodynamic, "STS", "Earth", 1 ) ); */

        dependentVariablesList.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >(
                        aerodynamic_force_coefficients_dependent_variable, "STS" ) );

        // Saving the angles
        dependentVariablesList.push_back(
                    boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                        "STS", reference_frames::AerodynamicsReferenceFrameAngles::bank_angle ) );
        dependentVariablesList.push_back(
                    boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                        "STS", reference_frames::AerodynamicsReferenceFrameAngles::flight_path_angle ) );
        dependentVariablesList.push_back(
                    boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                        "STS", reference_frames::AerodynamicsReferenceFrameAngles::angle_of_attack ) );

        // Angle of SideSlip
        dependentVariablesList.push_back(boost::make_shared< BodyAerodynamicAngleVariableSaveSettings>("STS", reference_frames::AerodynamicsReferenceFrameAngles::angle_of_sideslip));

        // Create object with list of dependent variables
        boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
                boost::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

        // Define termination conditions
        boost::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable1 =
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable, "STS", "Earth" );

        boost::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable2 =
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    mach_number_dependent_variable, "STS", "Earth" );

        boost::shared_ptr< PropagationTerminationSettings > terminationSettings1 =
                boost::make_shared< PropagationDependentVariableTerminationSettings >(
                    terminationDependentVariable1, 10.0E3, true );

        boost::shared_ptr< PropagationTerminationSettings > terminationSettings2 =
                boost::make_shared< PropagationDependentVariableTerminationSettings >(
                    terminationDependentVariable2, 3.0, true );

        boost::shared_ptr< PropagationTerminationSettings > terminationSettings3 =
                boost::make_shared< PropagationTimeTerminationSettings >(
                    simulationEndEpoch);

        //PropagationTimeTerminationSettings( const double terminationTime )
            std::vector< boost::shared_ptr< PropagationTerminationSettings > > terminationvector;
            terminationvector.push_back(terminationSettings1);
            terminationvector.push_back(terminationSettings2);
            terminationvector.push_back(terminationSettings3);

            boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
                    boost::make_shared< PropagationHybridTerminationSettings >(
                        terminationvector, true );
            // Create propagation settings.
               boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                       boost::make_shared< TranslationalStatePropagatorSettings< double > >
                       ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                         terminationSettings, cowell, dependentVariablesToSave );
               boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                       boost::make_shared< IntegratorSettings< > >
                       ( rungeKutta4, simulationStartEpoch, fixedStepSize );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////                            PROPAGATE ORBIT                              ///////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



     // Create simulation object and propagate dynamics.
         SingleArcDynamicsSimulator< > dynamicsSimulator(
                     bodyMap, integratorSettings, propagatorSettings, true, false, false );


   //   std::string outputDirectory = tudat_applications::getOutputPath( );
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////                                SAVE OUTPUT                              ///////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Write Asterix propagation history to file.
    input_output::writeDataMapToTextFile( dynamicsSimulator.getEquationsOfMotionNumericalSolution( ),
                                          "Results22d.dat",
                                          tudat_applications::getStsOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write Asterix propagation history to file.
    input_output::writeDataMapToTextFile( dynamicsSimulator.getDependentVariableHistory( ),
                                          "Results_Dependent_Variables22d.dat",
                                          tudat_applications::getStsOutputPath( ),
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );



    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
