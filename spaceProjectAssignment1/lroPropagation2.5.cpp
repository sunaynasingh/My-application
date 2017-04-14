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
                             std::string( "lroPropagation_part2_5.cpp" ).length( ) ) ) + "PropagationResultsPart3";
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
    bodiesToCreate.push_back( "Saturn" );
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
    bodyMap[ "LRO" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "LRO" ]->setConstantBodyMass( 1200.0 );



    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );



    for( unsigned int assignmentQuestion = 1; assignmentQuestion <= 1; assignmentQuestion++ )
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
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfMoon;
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfSun;
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfMars;
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfJupiter;

        ///                                         ///
        ///  DEFINE YOUR ACCELERATION TYPES HERE    ///
        ///                                         ///
        ///                                         ///

        accelerationsOfLRO["Moon"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfLRO["Earth"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfLRO["Sun"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfLRO["Mars"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfLRO["Jupiter"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfLRO["Saturn"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationMap[  "LRO" ] = accelerationsOfLRO;
        bodiesToPropagate.push_back( "LRO" );
        centralBodies.push_back( "Moon" ); // propagate LRO w.r.t. Moon


        accelerationsOfEarth["Moon"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfEarth["Sun"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfEarth["Mars"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfEarth["Jupiter"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfEarth["Saturn"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationMap[  "Earth" ] = accelerationsOfEarth;
        bodiesToPropagate.push_back( "Earth" );
        centralBodies.push_back( "Sun" ); // propagate Earth w.r.t. Sun


        accelerationsOfMoon["Earth"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfMoon["Sun"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfMoon["Mars"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfMoon["Jupiter"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfMoon["Saturn"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationMap[  "Moon" ] = accelerationsOfMoon;
        bodiesToPropagate.push_back( "Moon" );
        centralBodies.push_back( "Earth" ); // propagate Moon w.r.t. Earth


        accelerationsOfSun["Moon"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfSun["Earth"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfSun["Mars"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfSun["Jupiter"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfSun["Saturn"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationMap[  "Sun" ] = accelerationsOfSun;
        bodiesToPropagate.push_back( "Sun" );
        centralBodies.push_back( "SSB" ); // propagate sun w.r.t. barycenter


        accelerationsOfMars["Moon"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfMars["Earth"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfMars["Sun"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfMars["Jupiter"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfMars["Saturn"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationMap[  "Mars" ] = accelerationsOfMars;
        bodiesToPropagate.push_back( "Mars" );
        centralBodies.push_back( "Sun" ); // propagate Mars w.r.t. Sun


        accelerationsOfJupiter["Moon"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfJupiter["Earth"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfJupiter["Sun"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfJupiter["Mars"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationsOfJupiter["Saturn"].push_back( boost::make_shared<AccelerationSettings>(basic_astrodynamics::central_gravity));
        accelerationMap[  "Jupiter" ] = accelerationsOfJupiter;
        bodiesToPropagate.push_back( "Jupiter" );
        centralBodies.push_back( "Sun" ); // propagate Jupiter w.r.t. Sun

unsigned int numberOfNumericalBodies = bodiesToPropagate.size( );
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
        initialKeplerElements[ argumentOfPeriapsisIndex ] = 36.0*3.0/180.0*tudat::mathematical_constants::PI;
        initialKeplerElements[ longitudeOfAscendingNodeIndex ] = 36.0*9.0/180.0*tudat::mathematical_constants::PI;
        initialKeplerElements[ trueAnomalyIndex ] = 36.0*9.0/180.0*tudat::mathematical_constants::PI;


        Eigen::VectorXd initialCartesianState =  Eigen::VectorXd::Zero(36);
        initialCartesianState.segment(0,6) = convertKeplerianToCartesianElements(
                    initialKeplerElements, bodyMap[ "Moon" ]->getGravityFieldModel( )->getGravitationalParameter( ) );
        if( centralBodies.at( 0 ) == "SSB" )
        {
            initialCartesianState += bodyMap[ "Moon" ]->getStateInBaseFrameFromEphemeris( simulationStartEpoch );
        }

        Eigen::Vector6d systemInitialStateJupiter = getInitialStatesOfBodies(
                   { "Jupiter" }, { "Sun" }, bodyMap, simulationStartEpoch );
        Eigen::Vector6d systemInitialStateMars = getInitialStatesOfBodies(
                    { "Mars" }, { "Sun" }, bodyMap, simulationStartEpoch );
        Eigen::Vector6d systemInitialStateMoon = getInitialStatesOfBodies(
                    {"Moon"}, {"Earth"}, bodyMap, simulationStartEpoch );
        Eigen::Vector6d systemInitialStateEarth = getInitialStatesOfBodies(
                    {"Earth"}, {"Sun"}, bodyMap, simulationStartEpoch );
        Eigen::Vector6d systemInitialStateSun = getInitialStatesOfBodies(
                    {"Sun"}, {"SSB"}, bodyMap, simulationStartEpoch );


        initialCartesianState.segment(6,6) = systemInitialStateEarth;
        initialCartesianState.segment(12,6) = systemInitialStateMoon;
        initialCartesianState.segment(18,6) =systemInitialStateSun;
        initialCartesianState.segment(24,6) = systemInitialStateMars;
        initialCartesianState.segment(30,6) =systemInitialStateJupiter;

        ///                                                         ///
        ///  DEFINE PROPAGATION AND INTEGRATION SETTINGS HERE       ///
        ///                                                         ///
        ///                                                         ///



       //  boost::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings;
       //  boost::shared_ptr< IntegratorSettings< > > integratorSettings;

       // std::cout<<"assignment" >>assignmentQuestion<<std::endl;
        TranslationalPropagatorType propagatortype = cowell;

       // propagation settings modified

        boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                       boost::make_shared< TranslationalStatePropagatorSettings< double > >
                       ( centralBodies, accelerationModelMap, bodiesToPropagate, initialCartesianState, simulationEndEpoch, propagatortype);

        // Apply the range Kutta 4 with several intervals

        double  fixedStepSize = 10.0;


        // Select the propagtor
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



    }

}
