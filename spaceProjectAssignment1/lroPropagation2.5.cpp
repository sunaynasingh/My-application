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
                               std::string( "lroPropagation2.cpp" ).length( ) ) ) + "PropagationResults";
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
    using namespace tudat::basic_astrodynamics;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "naif0009.tls");

    // Set simulation time settings.
    const double simulationStartEpoch = physical_constants::JULIAN_YEAR;
    const double simulationEndEpoch = physical_constants::JULIAN_YEAR + 7.0 * tudat::physical_constants::JULIAN_DAY;


    // Define propagation termination conditions (stop after 2 weeks).
    boost::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            boost::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch );

    // Define bodies in simulation.
    unsigned int totalNumberOfBodies = 6;
    std::vector< std::string > bodyNames;
    bodyNames.resize( totalNumberOfBodies );
    bodyNames[ 0 ] = "Moon";
    bodyNames[ 1 ] = "Earth";
    bodyNames[ 2 ] = "Mars";
    bodyNames[ 3 ] = "Jupiter";
    bodyNames[ 4 ] = "Sun";

    // Create bodies needed in simulation
                  std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
                          getDefaultBodySettings( bodyNames );
                  NamedBodyMap bodyMap = createBodies( bodySettings );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

       // Create spacecraft object.
    double vehicleMass = 1200.0;
    bodyMap[ "LRO" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "LRO" ]->setConstantBodyMass( vehicleMass );

     // setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    SelectedAccelerationMap accelerationMap;
    for( unsigned int i = 0; i < bodyNames.size( ); i++ )
                      {
                          std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > currentAccelerations;
     for( unsigned int j = 0; j < bodyNames.size( ); j++ )
     {
         // Create central gravity acceleration between each 2 bodies.
         if( i != j )
         {
             currentAccelerations[ bodyNames.at( j ) ].push_back(
                         boost::make_shared< AccelerationSettings >( central_gravity ) );\
         }
     }
     accelerationMap[ bodyNames.at( i ) ] = currentAccelerations;
 }
    unsigned int numberOfBodies = bodyNames.size( );


        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            CREATE ACCELERATIONS          ////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        // Define propagation settings.
//        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfLRO;



        // Define central bodies to use in propagation.
        std::vector< std::string > centralBodies;
        centralBodies.resize( numberOfBodies );
//        SelectedAccelerationMap accelerationMap;

            for( unsigned int i = 0; i < numberOfBodies; i++ )
            {
                // Set Earth as central body for Moon
                if( i == 0 )
                {
                    centralBodies[ i ] = "Earth";
                }
                // Set barycenter as central 'body' for Sun
                else if( i == 4 )
                {
                    centralBodies[ i ] = "SSB";
                }
                // Set Sun as central body for all planets
                else
                {
                    centralBodies[ i ] = "Sun";
                }
            }
            std::vector< std::string > bodiesToPropagate_a = bodyNames;

            // Get initial state vector as input to integration.
                                  Eigen::VectorXd PlanetInitialState = getInitialStatesOfBodies(
                                              bodiesToPropagate_a, centralBodies, bodyMap, simulationStartEpoch );
            bodyNames.resize(7);
            bodyNames[6] = "LRO";

       // Define list of bodies to propagate
        std::vector< std::string > bodiesToPropagate = bodyNames;
         unsigned int numberOfNumericalBodies = bodiesToPropagate.size( );
         centralBodies.resize( numberOfNumericalBodies );
                               centralBodies[6] = "Moon";
         std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > currentAccelerations2;
          for( unsigned int j = 0; j < 5; j++ )
          {
         // Create central gravity acceleration between each 2 bodies.
          currentAccelerations2[ bodyNames.at( j ) ].push_back(boost::make_shared< AccelerationSettings >( central_gravity ) );\
           }
          accelerationMap[ bodyNames.at( 6 ) ] = currentAccelerations2;
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


        ///                                                         ///
        ///  DEFINE PROPAGATION AND INTEGRATION SETTINGS HERE       ///
        ///                                                         ///
        ///
        ///
        Eigen::VectorXd systemInitialState(PlanetInitialState.rows() + Part2InitialCartesianState.rows()); // used to be  lroInitialCartesianStateSSB.rows()
                         systemInitialState << PlanetInitialState, Part2InitialCartesianState;  // used to be lroInitialCartesianStateSSB

                  double const fixedStepSize =10;
      // Define settings for propagation of translational dynamics.

                boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                           boost::make_shared< TranslationalStatePropagatorSettings < double > >
                           ( centralBodies, accelerationModelMap, bodiesToPropagate, Part2InitialCartesianState, terminationSettings,cowell);

                boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                           boost::make_shared< TranslationalStatePropagatorSettings< double > >
                           ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch,cowell );

                // Create list of propagation settings.
                std::vector<boost::shared_ptr< PropagatorSettings< double > > > propagatorSettingsVector;
                propagatorSettingsVector.push_back(translationalPropagatorSettings);

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

                // Compute map of Kepler elements
                std::vector< std::map< double, Eigen::VectorXd > > currentCartesianState;
                currentCartesianState.resize( numberOfNumericalBodies );
                for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = cartesianIntegrationResult.begin( );
                     stateIterator != cartesianIntegrationResult.end( ); stateIterator++ )

                    // Retrieve current Cartesian state (convert to Moon-centered frame if needed)

                {
                     for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
                      {
                        currentCartesianState[ i ][ stateIterator->first ] = stateIterator->second.segment( i * 6, 6 );
                       }
                 }


                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                ///////////////////////        PROVIDE OUTPUT TO FILE                        //////////////////////////////////////////
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
{
                // Write Asterix propagation history to file.
                input_output::writeDataMapToTextFile( currentCartesianState[ i ],
                                                      "lroOrbit_Soln_" +
                                                       bodyNames.at(i) + ".dat",
                                                      tudat_applications::getOutputPath( ),
                                                      "",
                                                      std::numeric_limits< double >::digits10,
                                                      std::numeric_limits< double >::digits10,
                                                      "," );

           }
         }


