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
#include<math.h>

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
                               std::string( "lroPropagationPart2.cpp" ).length( ) ) ) + "PropagationResults";
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
    bodiesToCreate.push_back( "Saturn" );

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

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );


    for( unsigned int assignmentQuestion = 1; assignmentQuestion <= 5; assignmentQuestion++ )
    {
        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;


        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            CREATE ACCELERATIONS          ////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Define propagation settings.
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfLRO;


        accelerationsOfLRO[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >(
                                                    basic_astrodynamics::central_gravity ) );
        accelerationsOfLRO[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
        accelerationsOfLRO[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::central_gravity ) );
        accelerationsOfLRO[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );
        accelerationsOfLRO[ "Saturn" ].push_back( boost::make_shared< AccelerationSettings >(
                                                      basic_astrodynamics::central_gravity ) );
        //accelerationsOfLRO[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >(
        //                                            basic_astrodynamics::central_gravity ) );

        ////////////////////////////for part 2.4/////////////////////////////////////////////

        if (assignmentQuestion == 4)
        {

            accelerationsOfLRO[ "Moon" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 12, 12 ) );

            accelerationsOfLRO[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(basic_astrodynamics::cannon_ball_radiation_pressure ) );
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        accelerationMap[  "LRO" ] = accelerationsOfLRO;
        bodiesToPropagate.push_back( "LRO" );
        centralBodies.push_back( "Moon" );

        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

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


        Eigen::Vector6d Part2InitialCartesianState =  convertKeplerianToCartesianElements(
                    initialKeplerElements, bodyMap[ "Moon" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
        if( centralBodies.at( 0 ) == "SSB" )
        {
            Part2InitialCartesianState += bodyMap[ "Moon" ]->getStateInBaseFrameFromEphemeris( simulationStartEpoch );
        }

        ///                                                         ///
        ///  DEFINE PROPAGATION AND INTEGRATION SETTINGS HERE       ///
        ///                                                         ///
        ///
        ///



        //        if ( assignmentQuestion == 1 )
        //        {
        //                                                             ///
        //for (int fixedStepSize=1; fixedStepSize<=1000;fixedStepSize=fixedStepSize*10)
        //{
        //        // Define settings for propagation of translational dynamics.

        //           boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
        //                   boost::make_shared< TranslationalStatePropagatorSettings < double > >
        //                   ( centralBodies, accelerationModelMap, bodiesToPropagate, Part2InitialCartesianState, terminationSettings,cowell);

        //           boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
        //                   boost::make_shared< TranslationalStatePropagatorSettings< double > >
        //                   ( centralBodies, accelerationModelMap, bodiesToPropagate, Part2InitialCartesianState, simulationEndEpoch,cowell );

        //        // Create list of propagation settings.
        //        std::vector<boost::shared_ptr< PropagatorSettings< double > > > propagatorSettingsVector;
        //        propagatorSettingsVector.push_back(translationalPropagatorSettings);

        //        boost::shared_ptr< IntegratorSettings< > > integratorSettings =
        //                boost::make_shared< IntegratorSettings< > >
        //                ( rungeKutta4, simulationStartEpoch, fixedStepSize );

        //        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        //        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //        // Create simulation object and propagate dynamics.
        //        SingleArcDynamicsSimulator< > dynamicsSimulator(
        //                    bodyMap, integratorSettings, propagatorSettings, true, false, false );
        //        std::map< double, Eigen::VectorXd > cartesianIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        //        std::map< double, Eigen::VectorXd > keplerianIntegrationResult;
        //        // Compute map of Kepler elements
        //        Eigen::Vector6d currentCartesianState;
        //        for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = cartesianIntegrationResult.begin( );
        //             stateIterator != cartesianIntegrationResult.end( ); stateIterator++ )
        //        {
        //            // Retrieve current Cartesian state (convert to Moon-centered frame if needed)
        //            currentCartesianState = stateIterator->second;


        //            if( centralBodies.at( 0 ) == "SSB" )
        //            {
        //                currentCartesianState -=
        //                        bodyMap[ "Moon" ]->getStateInBaseFrameFromEphemeris( stateIterator->first );
        //            }

        //            keplerianIntegrationResult[ stateIterator->first ] =
        //                    convertCartesianToKeplerianElements(
        //                        currentCartesianState, bodyMap[ "Moon" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
        //        }


        //        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //        ///////////////////////        PROVIDE OUTPUT TO FILE                        //////////////////////////////////////////
        //        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        //        // Write Asterix propagation history to file.
        //        input_output::writeDataMapToTextFile( cartesianIntegrationResult,
        //                                              "lroOrbit_Cowell_" +
        //                                              boost::lexical_cast< std::string >( assignmentQuestion ) + ".dat"+ boost::lexical_cast< std::string >( fixedStepSize )+".dat",
        //                                              tudat_applications::getOutputPath( ),
        //                                              "",
        //                                              std::numeric_limits< double >::digits10,
        //                                              std::numeric_limits< double >::digits10,
        //                                              "," );

        //        // Write Asterix propagation history to file.
        //        input_output::writeDataMapToTextFile( keplerianIntegrationResult,
        //                                              "lroOrbitKeplerian_Cowell_" +
        //                                              boost::lexical_cast< std::string >( assignmentQuestion ) + ".dat"+ boost::lexical_cast< std::string >( fixedStepSize  )+".dat",
        //                                              tudat_applications::getOutputPath( ),
        //                                              "",
        //                                              std::numeric_limits< double >::digits10,
        //                                              std::numeric_limits< double >::digits10,
        //                                              "," );

        //    }


        //        for (int fixedStepSize=1; fixedStepSize<=1000;fixedStepSize=fixedStepSize*10)
        //        {

        //           // Define settings for propagation of translational dynamics.

        //////                   boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
        //////                           boost::make_shared< TranslationalStatePropagatorSettings < double > >
        //////                           ( centralBodies, accelerationModelMap, bodiesToPropagate, Part2InitialCartesianState, terminationSettings, encke);

        //                   boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
        //                           boost::make_shared< TranslationalStatePropagatorSettings< double > >
        //                           ( centralBodies, accelerationModelMap, bodiesToPropagate, Part2InitialCartesianState, simulationEndEpoch, encke );

        //////                // Create list of propagation settings.
        //////                std::vector<boost::shared_ptr< PropagatorSettings< double > > > propagatorSettingsVector;
        //////                propagatorSettingsVector.push_back(translationalPropagatorSettings);

        //                boost::shared_ptr< IntegratorSettings< > > integratorSettings =
        //                        boost::make_shared< IntegratorSettings< > >
        //                        ( rungeKutta4, simulationStartEpoch, fixedStepSize );

        //                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        //                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //                // Create simulation object and propagate dynamics.
        //                SingleArcDynamicsSimulator< > dynamicsSimulator(
        //                            bodyMap, integratorSettings, propagatorSettings, true, false, false );
        //                std::map< double, Eigen::VectorXd > cartesianIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        //                std::map< double, Eigen::VectorXd > keplerianIntegrationResult;
        //                // Compute map of Kepler elements
        //                Eigen::Vector6d currentCartesianState;
        //                for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = cartesianIntegrationResult.begin( );
        //                     stateIterator != cartesianIntegrationResult.end( ); stateIterator++ )
        //                {
        //                    // Retrieve current Cartesian state (convert to Moon-centered frame if needed)
        //                    currentCartesianState = stateIterator->second;


        //                    if( centralBodies.at( 0 ) == "SSB" )
        //                    {
        //                        currentCartesianState -=
        //                                bodyMap[ "Moon" ]->getStateInBaseFrameFromEphemeris( stateIterator->first );
        //                    }

        //                    keplerianIntegrationResult[ stateIterator->first ] =
        //                            convertCartesianToKeplerianElements(
        //                                currentCartesianState, bodyMap[ "Moon" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
        //                }


        //                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                ///////////////////////        PROVIDE OUTPUT TO FILE                        //////////////////////////////////////////
        //                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        //                // Write Asterix propagation history to file.
        //                input_output::writeDataMapToTextFile( cartesianIntegrationResult,
        //                                                      "lroOrbit_Encke_" +
        //                                                      boost::lexical_cast< std::string >( assignmentQuestion ) + ".dat"+ boost::lexical_cast< std::string >( fixedStepSize )+".dat",
        //                                                      tudat_applications::getOutputPath( ),
        //                                                      "",
        //                                                      std::numeric_limits< double >::digits10,
        //                                                      std::numeric_limits< double >::digits10,
        //                                                      "," );

        //            }

//        if (assignmentQuestion == 2 || assignmentQuestion == 4)
//        {
//            for (long double tolerance=pow(10,-14); tolerance<=pow(10,-10);tolerance=tolerance*10)
//            {

//                // Define settings for propagation of translational dynamics.

//                ////                   boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
//                ////                           boost::make_shared< TranslationalStatePropagatorSettings < double > >
//                ////                           ( centralBodies, accelerationModelMap, bodiesToPropagate, Part2InitialCartesianState, terminationSettings, encke);

//                boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
//                        boost::make_shared< TranslationalStatePropagatorSettings< double > >
//                        ( centralBodies, accelerationModelMap, bodiesToPropagate, Part2InitialCartesianState, simulationEndEpoch, encke );

//                ////                // Create list of propagation settings.
//                ////                std::vector<boost::shared_ptr< PropagatorSettings< double > > > propagatorSettingsVector;
//                ////                propagatorSettingsVector.push_back(translationalPropagatorSettings);

//                boost::shared_ptr< IntegratorSettings< > > integratorSettings =
//                        boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
//                        ( rungeKuttaVariableStepSize, simulationStartEpoch, 10.0, RungeKuttaCoefficients::rungeKuttaFehlberg78,
//                          0.05, 1000.0, tolerance, tolerance);

//                ////////////////////////////////////////////////////////////////////////////////////
//                ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
//                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//                // Create simulation object and propagate dynamics.
//                SingleArcDynamicsSimulator< > dynamicsSimulator(
//                            bodyMap, integratorSettings, propagatorSettings, true, false, false );
//                std::map< double, Eigen::VectorXd > cartesianIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

//                // Compute map of Kepler elements
//                Eigen::Vector6d currentCartesianState;
//                for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = cartesianIntegrationResult.begin( );
//                     stateIterator != cartesianIntegrationResult.end( ); stateIterator++ )
//                {
//                    // Retrieve current Cartesian state (convert to Moon-centered frame if needed)
//                    currentCartesianState = stateIterator->second;

//                }


//                //                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                //                ///////////////////////        PROVIDE OUTPUT TO FILE                        //////////////////////////////////////////
//                //                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//                // Write Asterix propagation history to file.
//                input_output::writeDataMapToTextFile( cartesianIntegrationResult,
//                                                      "lroOrbit_Encke_RK78" +
//                                                      boost::lexical_cast< std::string >( assignmentQuestion ) + boost::lexical_cast< std::string >( tolerance )+".dat",
//                                                      tudat_applications::getOutputPath( ),
//                                                      "",
//                                                      std::numeric_limits< double >::digits10,
//                                                      std::numeric_limits< double >::digits10,
//                                                      "," );


//            }


                if (assignmentQuestion == 2 || assignmentQuestion == 4)
             {
                                for (long double tolerance=pow(10,-14); tolerance<=pow(10,-10);tolerance=tolerance*10)
                                {

                                    // Define settings for propagation of translational dynamics.

                        ////                   boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                        ////                           boost::make_shared< TranslationalStatePropagatorSettings < double > >
                        ////                           ( centralBodies, accelerationModelMap, bodiesToPropagate, Part2InitialCartesianState, terminationSettings, encke);

                                           boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                                                   boost::make_shared< TranslationalStatePropagatorSettings< double > >
                                                   ( centralBodies, accelerationModelMap, bodiesToPropagate, Part2InitialCartesianState, simulationEndEpoch, cowell );

                        ////                // Create list of propagation settings.
                        ////                std::vector<boost::shared_ptr< PropagatorSettings< double > > > propagatorSettingsVector;
                        ////                propagatorSettingsVector.push_back(translationalPropagatorSettings);

                                           boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                                         boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                                         ( rungeKuttaVariableStepSize, simulationStartEpoch, 10.0, RungeKuttaCoefficients::rungeKuttaFehlberg78,
                                           0.05, 1000.0, tolerance, tolerance);

            ////////////////////////////////////////////////////////////////////////////////////
                                        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
                                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                                        // Create simulation object and propagate dynamics.
                                        SingleArcDynamicsSimulator< > dynamicsSimulator(
                                                    bodyMap, integratorSettings, propagatorSettings, true, false, false );
                                        std::map< double, Eigen::VectorXd > cartesianIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

                                        // Compute map of Kepler elements
                                        Eigen::Vector6d currentCartesianState;
                                        for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = cartesianIntegrationResult.begin( );
                                             stateIterator != cartesianIntegrationResult.end( ); stateIterator++ )
                                        {
                                            // Retrieve current Cartesian state (convert to Moon-centered frame if needed)
                                            currentCartesianState = stateIterator->second;

            }



                        //                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        //                ///////////////////////        PROVIDE OUTPUT TO FILE                        //////////////////////////////////////////
                        //                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


                                        // Write Asterix propagation history to file.
                                        input_output::writeDataMapToTextFile( cartesianIntegrationResult,
                                                                              "lroOrbit_Cowell_RK78" +
                                                                              boost::lexical_cast< std::string >( assignmentQuestion ) + boost::lexical_cast< std::string >( tolerance )+".dat",
                                                                              tudat_applications::getOutputPath( ),
                                                                              "",
                                                                              std::numeric_limits< double >::digits10,
                                                                              std::numeric_limits< double >::digits10,
                                                                              "," );



                }
        }
    }
}
