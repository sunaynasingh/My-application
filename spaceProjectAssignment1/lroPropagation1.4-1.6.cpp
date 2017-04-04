/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

namespace tudat_applications
{

//! Get path for output directory.
static inline std::string getOutputPath( )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    return ( filePath_.substr( 0, filePath_.length( ) -
                               std::string( "lroPropagation.cpp" ).length( ) ) ) + "PropagationResults";
}

}

int main( )
{

    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::gravitation;
    using namespace tudat::numerical_integrators;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    // Set simulation time settings.
    const double simulationStartEpoch = physical_constants::JULIAN_YEAR;
    const double simulationEndEpoch = physical_constants::JULIAN_YEAR + 7.0 * tudat::physical_constants::JULIAN_DAY;


    // Define propagation termination conditions (stop after 2 weeks).
    boost::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            boost::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch );

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Jupiter" );

    // Create body objects.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }

    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Create spacecraft object.
       double vehicleMass = 1200.0;
       bodyMap[ "LRO" ] = boost::make_shared< simulation_setup::Body >( );
       bodyMap[ "LRO" ]->setConstantBodyMass( vehicleMass );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );



    for( unsigned int assignmentQuestion = 4; assignmentQuestion <= 12; assignmentQuestion++ )
    {
        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;


        // Define list of dependent variables to save.
        std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariableslist;
        // Define propagation settings.
         std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfLRO;



    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            CREATE ACCELERATIONS          ////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



        // Create radiation pressure settings
        double referenceAreaRadiation = 8.0;
        double radiationPressureCoefficient = 1.2;
        std::vector< std::string > occultingBodies;
        occultingBodies.push_back( "Moon" );
        boost::shared_ptr< RadiationPressureInterfaceSettings > LRORadiationPressureSettings =
                boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                    "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

        // Create and set radiation pressure settings
        bodyMap[ "LRO" ]->setRadiationPressureInterface(
                    "Sun", createRadiationPressureInterface(
                        LRORadiationPressureSettings, "LRO", bodyMap ) );
//        // Define propagation settings.
//        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfLRO;


        if ( assignmentQuestion == 4 )
        {

            accelerationsOfLRO[ "Moon" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 16, 16 ) );
            dependentVariableslist.push_back(
                        boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(basic_astrodynamics::spherical_harmonic_gravity, "LRO", "Moon", 1 ) );
            accelerationsOfLRO[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >(
                                                        basic_astrodynamics::central_gravity ) );
            dependentVariableslist.push_back(
                        boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(basic_astrodynamics::third_body_central_gravity, "LRO", "Mars", 1 ) );
            accelerationsOfLRO[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );
            dependentVariableslist.push_back(
                        boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(basic_astrodynamics::third_body_central_gravity, "LRO", "Earth", 1 ) );

            accelerationsOfLRO[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >(
                                                           basic_astrodynamics::central_gravity ) );
            dependentVariableslist.push_back(
                        boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(basic_astrodynamics::third_body_central_gravity, "LRO", "Jupiter", 1 ) );

            accelerationsOfLRO[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::central_gravity ) );
            dependentVariableslist.push_back(
                        boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(basic_astrodynamics::third_body_central_gravity, "LRO", "Sun", 1 ) );
            accelerationsOfLRO[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::cannon_ball_radiation_pressure ) );
            dependentVariableslist.push_back(
                                  boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(basic_astrodynamics::cannon_ball_radiation_pressure, "LRO", "Sun", 1) );

        }

       else if( assignmentQuestion == 5)    //without radiaiton pressure 1
        {
            accelerationsOfLRO[ "Moon" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 16, 16 ) );

            accelerationsOfLRO[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >(
                                                        basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >(
                                                           basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::central_gravity ) );

            }

        else if( assignmentQuestion == 6)    //without sun case 2
        {
            accelerationsOfLRO[ "Moon" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 16, 16 ) );

            accelerationsOfLRO[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >(
                                                        basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >(
                                                           basic_astrodynamics::central_gravity ) );
            //  accelerationsOfLRO[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
            //                                                   basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(basic_astrodynamics::cannon_ball_radiation_pressure ) );
        }

        else if( assignmentQuestion == 7 )    //without mars case
        {
            accelerationsOfLRO[ "Moon" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 16, 16 ) );


            //accelerationsOfLRO[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >(
            //                                                basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >(
                                                           basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::cannon_ball_radiation_pressure ));
        }
        else if( assignmentQuestion == 8)    //without Jupiter case 2
        {
            accelerationsOfLRO[ "Moon" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 16, 16 ) );

            accelerationsOfLRO[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >(
                                                        basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );
            //accelerationsOfLRO[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >(
            //                                                 basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::cannon_ball_radiation_pressure ));
        }
        else if( assignmentQuestion == 9 )    //without Earth case 2
        {
            accelerationsOfLRO[ "Moon" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 16, 16 ) );


            accelerationsOfLRO[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >(
                                                        basic_astrodynamics::central_gravity ) );
            //accelerationsOfLRO[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
            //                                                  basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >(
                                                           basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::cannon_ball_radiation_pressure ));
        }
        else if( assignmentQuestion == 10 )    //without spherical case 2
        {

            accelerationsOfLRO[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >(
                                                        basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >(
                                                        basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >(
                                                             basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::central_gravity ) );
            accelerationsOfLRO[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::cannon_ball_radiation_pressure ));
        }

      if (assignmentQuestion == 11 ||assignmentQuestion == 12)
        {
          // Define acceleration model settings
          accelerationsOfLRO[ "Moon" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 16, 16 ) );

          accelerationsOfLRO[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >(
                                                      basic_astrodynamics::central_gravity ) );
          accelerationsOfLRO[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::central_gravity ) );
          accelerationsOfLRO[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );
          accelerationsOfLRO[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                             basic_astrodynamics::central_gravity ) );
          accelerationsOfLRO[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(basic_astrodynamics::cannon_ball_radiation_pressure ) );
          // Define thrust settings
                      double thrustMagnitude = 1.0;
                      double specificImpulse = 200.0;
                      boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings;
                      if (assignmentQuestion == 11)
                      {
                              // VeloThrust
                          thrustDirectionGuidanceSettings = boost::make_shared < ThrustDirectionFromStateGuidanceSettings >(
                                      "Moon", true, false);
                      }
                              else if(assignmentQuestion == 12)
                      {
                                      // PosiThrust
                          thrustDirectionGuidanceSettings = boost::make_shared < ThrustDirectionFromStateGuidanceSettings >(
                                      "Moon", false, true);
                      }
                      dependentVariableslist.push_back(
                                             boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                                                 basic_astrodynamics::thrust_acceleration, "LRO", "LRO", 0 ) );

                      boost::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings =
                              boost::make_shared< ConstantThrustEngineSettings >(
                                  thrustMagnitude, specificImpulse);

                      // Define acceleration model settings
                      accelerationsOfLRO[ "LRO" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                                 thrustDirectionGuidanceSettings, thrustMagnitudeSettings) );
                  }


        accelerationMap[  "LRO" ] = accelerationsOfLRO;
        bodiesToPropagate.push_back( "LRO" );
        centralBodies.push_back( "Moon" );

        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        ///                                                           ///
        ///  MODIFY CENTRAL BODY SETTINGS FOR SUB-QUESTION IF NEEDED  ///
        ///                                                           ///
        ///


        // Set initial Kepler elements for LRO.
        Eigen::Vector6d initialKeplerElements;
        initialKeplerElements[ semiMajorAxisIndex ] = 1900.0E3;
        initialKeplerElements[ eccentricityIndex ] = 0.05;
        initialKeplerElements[ inclinationIndex ] = 89.7 * mathematical_constants::PI / 180.0;

        ///                                                           ///
        ///  MODIFY INITIAL STATE ACCORDING TO STUDENT NUMBER         ///
        ///                                                           ///
        ///                                                           ///
        initialKeplerElements[ argumentOfPeriapsisIndex ] = 1.885;
        initialKeplerElements[ longitudeOfAscendingNodeIndex ] = 5.655;
        initialKeplerElements[ trueAnomalyIndex ] = 5.655;


        Eigen::Vector6d lroInitialCartesianState =  convertKeplerianToCartesianElements(
                    initialKeplerElements, bodyMap[ "Moon" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
        if( centralBodies.at( 0 ) == "SSB" )
        {
            lroInitialCartesianState += bodyMap[ "Moon" ]->getStateInBaseFrameFromEphemeris( simulationStartEpoch );
        }

        ///                                                         ///
        ///  DEFINE PROPAGATION AND INTEGRATION SETTINGS HERE       ///
        ///                                                         ///
        ///
        ///                                                         ///

        // Define settings for propagation of translational dynamics.

           boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                   boost::make_shared< TranslationalStatePropagatorSettings < double > >
                   ( centralBodies, accelerationModelMap, bodiesToPropagate, lroInitialCartesianState, terminationSettings);

          // Create object with list of dependent variables
                          boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
                                  boost::make_shared< DependentVariableSaveSettings >( dependentVariableslist );


        // Create list of propagation settings.
        std::vector<boost::shared_ptr< PropagatorSettings< double > > > propagatorSettingsVector;
        propagatorSettingsVector.push_back(translationalPropagatorSettings);

if (assignmentQuestion == 11  || assignmentQuestion == 12 )
{

// Create mass rate models
                    boost::shared_ptr< MassRateModelSettings > massRateModelSettings =
                            boost::make_shared< FromThrustMassModelSettings>(true) ;
                    std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
                    massRateModels[ "LRO" ] = createMassRateModel(
                                "LRO", massRateModelSettings, bodyMap, accelerationModelMap );

                    // Create settings for propagating the mass of the vehicle.
                    std::vector< std::string > bodiesWithMassToPropagate;
                    bodiesWithMassToPropagate.push_back( "LRO" );

                    Eigen::VectorXd initialBodyMasses = Eigen::VectorXd( 1 );
                    initialBodyMasses( 0 ) = vehicleMass;

                    boost::shared_ptr< PropagatorSettings< double > > massPropagatorSettings =
                            boost::make_shared< MassPropagatorSettings< double > >(
                                bodiesWithMassToPropagate, massRateModels, initialBodyMasses, terminationSettings );

                    propagatorSettingsVector.push_back(massPropagatorSettings);
            }
// Create propagation settings for concurrent dynamics
      boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
              boost::make_shared< MultiTypePropagatorSettings< double > >(
                  propagatorSettingsVector, terminationSettings, dependentVariablesToSave );


        const double fixedStepSize = 30.0;

//        boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
//             boost::make_shared< TranslationalStatePropagatorSettings< double > >
//                ( centralBodies, accelerationModelMap, bodiesToPropagate, lroInitialCartesianState, simulationEndEpoch, cowell,
//                  boost::make_shared< DependentVariableSaveSettings >( dependentVariableslist ));

        boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                boost::make_shared< IntegratorSettings< > >
                ( rungeKutta4, simulationStartEpoch, fixedStepSize );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true, false, false );
        std::map< double, Eigen::VectorXd > cartesianIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > keplerianIntegrationResult;
        // std::map< double, Eigen::VectorXd > dependentVariables = dynamicsSimulator.getDependentVariableHistory();
        // Compute map of Kepler elements
        Eigen::Vector6d currentCartesianState;
        for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = cartesianIntegrationResult.begin( );
             stateIterator != cartesianIntegrationResult.end( ); stateIterator++ )
        {
            // Retrieve current Cartesian state (convert to Moon-centered frame if needed)
            currentCartesianState = stateIterator->second;


            if( centralBodies.at( 0 ) == "SSB" )
            {
                currentCartesianState -=
                        bodyMap[ "Moon" ]->getStateInBaseFrameFromEphemeris( stateIterator->first );
            }

            keplerianIntegrationResult[ stateIterator->first ] =
                    convertCartesianToKeplerianElements(
                        currentCartesianState, bodyMap[ "Moon" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
        }


        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////        PROVIDE OUTPUT TO FILE                        //////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // Write Asterix propagation history to file.
        input_output::writeDataMapToTextFile( cartesianIntegrationResult,
                                              "lroOrbit_" +
                                              boost::lexical_cast< std::string >( assignmentQuestion ) + ".dat",
                                              tudat_applications::getOutputPath( ),
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );

        // Write Asterix propagation history to file.
        input_output::writeDataMapToTextFile( keplerianIntegrationResult,
                                              "lroOrbitKeplerian_" +
                                              boost::lexical_cast< std::string >( assignmentQuestion ) + ".dat",
                                              tudat_applications::getOutputPath( ),
                                              "",
                                              std::numeric_limits< double >::digits10,
                                              std::numeric_limits< double >::digits10,
                                              "," );


        // Save barycentric Moon state (only for 1.4)
        if( assignmentQuestion == 4 )
        {

            // Write dependent LRO propagation history to file.
            input_output::writeDataMapToTextFile( dynamicsSimulator.getDependentVariableHistory(),
                                                  "dependentvariableshistory.dat",
                                                  tudat_applications::getOutputPath( ),
                                                  "",
                                                  std::numeric_limits< double >::digits10,
                                                  std::numeric_limits< double >::digits10,
                                                  "," );
      }
    }

}
