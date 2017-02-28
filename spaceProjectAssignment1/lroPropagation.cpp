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

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Venus" );

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
    bodyMap[ "LRO" ] = boost::make_shared< simulation_setup::Body >( );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );



    for( unsigned int assignmentQuestion = 1; assignmentQuestion <= 3; assignmentQuestion++ )
    {
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            CREATE ACCELERATIONS          ////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        // Define propagation settings.
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfLRO;

        ///                                         ///
        ///  DEFINE YOUR ACCELERATION TYPES HERE    ///
        ///                                         ///
        ///                                         ///

        accelerationMap[  "LRO" ] = accelerationsOfLRO;
        bodiesToPropagate.push_back( "LRO" );

        ///                                                           ///
        ///  MODIFY CENTRAL BODY SETTINGS FOR SUB-QUESTION IF NEEDED  ///
        ///                                                           ///
        ///                                                           ///
        if( assignmentQuestion == 1 )
        {
            centralBodies.push_back( "SSB" );
        }
        else
        {
            std::cerr<<"Error, central bodies not specified for this subquestion."<<std::endl;
        }

        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        // Set initial Kepler elements for LRO.
        Eigen::Vector6d initialKeplerElements;
        initialKeplerElements[ semiMajorAxisIndex ] = 1900.0E3;
        initialKeplerElements[ eccentricityIndex ] = 0.05;
        initialKeplerElements[ inclinationIndex ] =
                89.7 * mathematical_constants::PI / 180.0;

        ///                                                           ///
        ///  MODIFY INITIAL STATE ACCORDING TO STUDENT NUMBER         ///
        ///                                                           ///
        ///                                                           ///
        initialKeplerElements[ argumentOfPeriapsisIndex ] = 30.0;
        initialKeplerElements[ longitudeOfAscendingNodeIndex ] = 30.0;
        initialKeplerElements[ trueAnomalyIndex ] = 30.0;


        Eigen::Vector6d lroInitialCartesianState =  convertKeplerianToCartesianElements(
                    initialKeplerElements, bodyMap[ "Moon" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
        if( centralBodies.at( 0 ) == "SSB" )
        {
            lroInitialCartesianState += bodyMap[ "Moon" ]->getStateInBaseFrameFromEphemeris( simulationStartEpoch );
        }

        ///                                                         ///
        ///  DEFINE PROPAGATION AND INTEGRATION SETTINGS HERE       ///
        ///                                                         ///
        ///                                                         ///
        boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings;
        boost::shared_ptr< IntegratorSettings< > > integratorSettings;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true, false, false );
        std::map< double, Eigen::VectorXd > cartesianIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > keplerianIntegrationResult;

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


        // Save barycentric Moon state (only for 1.1)
        if( assignmentQuestion == 1  )
        {
            std::map< double, Eigen::VectorXd > moonBarycentricStates;
            for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = cartesianIntegrationResult.begin( );
                 stateIterator != cartesianIntegrationResult.end( ); stateIterator++ )
            {
                moonBarycentricStates[ stateIterator->first ] = bodyMap.at( "Moon" )->getStateInBaseFrameFromEphemeris(
                            stateIterator->first );
            }

            // Write Asterix propagation history to file.
            input_output::writeDataMapToTextFile( cartesianIntegrationResult,
                                                  "moonBarycentricStates.dat",
                                                  tudat_applications::getOutputPath( ),
                                                  "",
                                                  std::numeric_limits< double >::digits10,
                                                  std::numeric_limits< double >::digits10,
                                                  "," );
        }
    }

}


