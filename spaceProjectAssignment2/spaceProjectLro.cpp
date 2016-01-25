/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      150101    D. Dirkx          File created
 *
 *    References
 *
 *    Notes
 *
 */

#include <fstream>
#include <limits>
#include <string>
#include <utility>
#include <iostream>

#include <boost/assign/list_of.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include <Eigen/Core>

#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h>
#include <Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h>
#include <Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h>
#include <Tudat/Astrodynamics/Gravitation/centralGravityModel.h>
#include <Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h>
#include <Tudat/Astrodynamics/Gravitation/thirdBodyPerturbation.h>
#include <Tudat/Astrodynamics/StateDerivativeModels/cartesianStateDerivativeModel.h>
#include <Tudat/Astrodynamics/StateDerivativeModels/compositeStateDerivativeModel.h>
#include <Tudat/External/SpiceInterface/spiceEphemeris.h>
#include <Tudat/InputOutput/basicInputOutput.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
#include <Tudat/Mathematics/Interpolators/lagrangeInterpolator.h>

#include "LroPropagation/body.h"

//! Function to return the path of the current file.
/*!
 *  Function to return the path of the current file.
 *  \return Path of current file.
 */
static inline std::string getSpaceProjectApplicationPath( )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path in the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    return filePath_.substr( 0, filePath_.length( ) -
                                std::string( "spaceProjectLro.cpp" ).length( ) );
}

//! Function to read a gravity field file from the Planetary Data System (PDS).
/*!
 *  Function to read a gravity field file from the Planetary Data System (PDS),
 *  returns the gravitational parameters, reference radius and sine and cosine
 *  spherical harmomic coefficients.
 *  \param fileName Name of PDS gravity field file to be loaded.
 *  \param maximumDegree Maximum degree of gravity field to be loaded.
 *  \param maximumOrder Maximum order of gravity field to be loaded.
 *  \param gravitationalParameter Gravitational parameter of body (reference).
 *  \param referenceRadius Reference radius of body (reference).
 *  \return Spherical harmonics coefficients.
 */
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > readPdsGravityFieldFile(
        const std::string& fileName, const int maximumDegree, const int maximumOrder,
        double& gravitationalParameter, double& referenceRadius )
{
    // Attempt to open gravity file.
    std::fstream stream( fileName.c_str( ), std::ios::in );
    if( stream.fail( ) )
    {
        boost::throw_exception(
                    std::runtime_error( "Pds gravity field data file could not be opened." ) );
    }

    // Declare variables for reading file.
    std::vector< std::string > vectorOfIndividualStrings;
    vectorOfIndividualStrings.resize( 4 );
    std::string line;

    // Get first line of file.
    std::getline( stream, line );

    // Get reference radius and gravitational parameter from first line of file.
    boost::algorithm::trim( line );
    boost::algorithm::split( vectorOfIndividualStrings,
                             line,
                             boost::algorithm::is_any_of( ", " ),
                             boost::algorithm::token_compress_on );
    gravitationalParameter = boost::lexical_cast< double >( vectorOfIndividualStrings[ 1 ] );
    referenceRadius = boost::lexical_cast< double >( vectorOfIndividualStrings[ 0 ] );

    // Declare variables for reading in cosine and sine coefficients.
    int currentDegree = 0, currentOrder = 0;
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd( maximumDegree + 1, maximumOrder + 1 );
    cosineCoefficients.setZero( );
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd( maximumDegree + 1, maximumOrder + 1 );
    sineCoefficients.setZero( );

    // Read coefficients up to required maximum degree and order.
    while ( !stream.fail( ) && !stream.eof( ) &&
            ( currentDegree <= maximumDegree || currentOrder <= maximumOrder )  )
    {
        // Read current line
        std::getline( stream, line );

        // Trim input string (removes all leading and trailing whitespaces).
        boost::algorithm::trim( line );

        // Split string into multiple strings, each containing one element from a line from the
        // data file.
        boost::algorithm::split( vectorOfIndividualStrings,
                                 line,
                                 boost::algorithm::is_any_of( ", " ),
                                 boost::algorithm::token_compress_on );

        // Check current line for consistency
        if( vectorOfIndividualStrings.size( ) < 4 )
        {
            std::cerr<<"Error when reading pds gravity field file, number of fields is "
                    <<vectorOfIndividualStrings.size( )<<std::endl;
        }
        else
        {
            // Read current degree and orde from line.
            currentDegree = boost::lexical_cast< int >( vectorOfIndividualStrings[ 0 ] );
            currentOrder = boost::lexical_cast< int >( vectorOfIndividualStrings[ 1 ] );

            // Set cosine and sine coefficients for current degree and order.
            if( currentDegree <= maximumDegree && currentOrder <= maximumOrder )
            {
                cosineCoefficients( currentDegree, currentOrder ) =
                        boost::lexical_cast< double >( vectorOfIndividualStrings[ 2 ] );
                sineCoefficients( currentDegree, currentOrder ) =
                        boost::lexical_cast< double >( vectorOfIndividualStrings[ 3 ] );
            }
        }
    }

    // Set cosine coefficient at (0,0) to 1.
    cosineCoefficients( 0, 0 ) = 1.0;

    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > coefficients;
    coefficients = std::make_pair( cosineCoefficients, sineCoefficients );

    return coefficients;
}

//! Function to rotate a vector from a body-fixed to an inertial frame.
/*!
 *  Function to rotate a vector from a body-fixed to an inertial frame, using curret
 *  rotational state of a given body.
 *  \param vectorInBodyFixedFrame Original vector, given in body-fixed frame.
 *  \param bodyWithFixedFrame Body object from which current rotational state is taken.
 *  \return vectorInBodyFixedFrame Rotated to inertial frame.
 */
Eigen::Vector3d rotateVectorFromBodyFixedToInertialFrame(
        const Eigen::Vector3d& vectorInBodyFixedFrame,
        const boost::shared_ptr< tudat_application::Body > bodyWithFixedFrame )
{
    return bodyWithFixedFrame->getCurrentRotationToGlobalFrame( ) * vectorInBodyFixedFrame;
}

//! Function to get position of body A in a frame corotating w.r.t.
//! (but not necesarilly centerd on) body B.
/*!
 *  \param bodyWithPosition Body for which its current position is rotated to a corotating frame.
 *  \param bodyWithFixedFrame Body to the corotating frame of which the current position of
 *  bodyWithPosition is to be transformed
 *  \return Current position of bodyWithPosition in a frame corotating with bodyWithFixedFrame.
 */
Eigen::Vector3d getPositionVectorWrtBodyFixedFrame(
        const boost::shared_ptr< tudat_application::Body > bodyWithPosition,
        const boost::shared_ptr< tudat_application::Body > bodyWithFixedFrame )
{
    return bodyWithFixedFrame->getCurrentRotationToLocalFrame( ) *
            bodyWithPosition->getCurrentPosition( );
}

//! Function to create a central gravitational acceleration
/*!
 *  Function to create a central gravitational acceleration, i.e. as exerted by a radially
 *  homogeneous body.
 *  \param bodyExertingAcceleration Object for body with gravity field exerting acceleration.
 *  \param bodyUndergoingAcceleration Object for body with undergoing acceleration.
 *  \return Point mass gravitational acceleration from bodyExertingAcceleration on
 *  bodyUndergoingAcceleration
 */
boost::shared_ptr< tudat::gravitation::CentralGravitationalAccelerationModel3d >
createCentralGravityAcceleration(
        const boost::shared_ptr< tudat_application::Body > bodyExertingAcceleration,
        const boost::shared_ptr< tudat_application::Body > bodyUndergoingAcceleration )
{
    // Retrieve the gravitational parameter from the bodyExertingAcceleration,
    // give warning if unsuccesful.
    double bodyGravitationalParameter = 0.0;
    try
    {
        bodyGravitationalParameter =
                bodyExertingAcceleration->getGravityFieldModel( )->getGravitationalParameter( );
    }
    catch( std::runtime_error )
    {
        std::cerr<<"Error when making central gravity acceleration, "<<
                   "did not find a gravity field model"<<std::endl;
    }

    // Create functions retrieve current state of bodyu undergoing and body exerting acceleration.
    boost::function< Eigen::Vector3d( ) > positionFunctionOfBodyExertingAcceleration =
            boost::bind( &tudat_application::Body::getCurrentPosition, bodyExertingAcceleration );
    boost::function< Eigen::Vector3d( ) > positionFunctionOfBodyUndergoingAcceleration =
            boost::bind( &tudat_application::Body::getCurrentPosition, bodyUndergoingAcceleration );

    // Create and return gravitational acceleration
    return boost::make_shared< tudat::gravitation::CentralGravitationalAccelerationModel3d >(
                positionFunctionOfBodyUndergoingAcceleration, bodyGravitationalParameter,
                positionFunctionOfBodyExertingAcceleration );
}

//! Function to create a third body central gravitational acceleration
/*!
 *  Function to create a third body central gravitational acceleration, i.e. the gravitational
 *  acceleration due to a (so-called third) body on a body w.r.t. to a central body/
 *  \param bodyExertingAcceleration Object for body with gravity field exerting acceleration.
 *  \param bodyUndergoingAcceleration Object for body with undergoing acceleration.
 *  \param centralBody Central body w.r.t. which acceleration on bodyUndergoingAcceleration
 *  due to bodyExertingAcceleration is calculated.
 *  \return Third body central body gravitational acceleration from bodyExertingAcceleration on
 *  bodyUndergoingAcceleration, w.r.t. centralBody.
 */
boost::shared_ptr< tudat::gravitation::ThirdBodyCentralGravityAcceleration >
createThirdBodyCentralGravityAcceleration(
        const boost::shared_ptr< tudat_application::Body > bodyExertingAcceleration,
        const boost::shared_ptr< tudat_application::Body > bodyUndergoingAcceleration,
        const boost::shared_ptr< tudat_application::Body > centralBody)
{
    // Create respective central gravitational acceleration and create third body acceleration.
    return boost::make_shared< tudat::gravitation::ThirdBodyCentralGravityAcceleration >(
                createCentralGravityAcceleration( bodyExertingAcceleration,
                                                  bodyUndergoingAcceleration ),
                createCentralGravityAcceleration( bodyExertingAcceleration,
                                                  centralBody ) );
}

//! Function to create a spherical harmonic gravitational acceleration
/*!
  *  Function to create a spherical harmonic gravitational acceleration, i.e. with
  *  gravity field expanded as a spherical harmonic series.
  *  \param bodyExertingAcceleration Object for body with spherical harmonic gravity
  *  field exerting acceleration.
  *  \param bodyUndergoingAcceleration Object for body with undergoing acceleration.
  *  \param maximumDegree Degree to which to expand spherical harmonic gravity field.
  *  \param maximumOrder Order to which to expand spherical harmonic gravity field.
  *  \return Spherical harmonic gravitational acceleration from bodyExertingAcceleration on
  *  bodyUndergoingAcceleration.
  */
boost::shared_ptr< tudat::basic_astrodynamics::AccelerationModel3d >
createSphericalHarmonicGravityAcceleration(
        const boost::shared_ptr< tudat_application::Body > bodyExertingAcceleration,
        const boost::shared_ptr< tudat_application::Body > bodyUndergoingAcceleration,
        const int maximumDegree, const int maximumOrder )
{
    // Declare characteristics of spherical harmonic gravity field.
    double bodyGravitationalParameter = 0.0;
    double referenceRadius = 0.0;
    Eigen::MatrixXd cosineCoefficients, sineCoefficients;

    // Check if body has spherical harmonic gravity field
    if( boost::dynamic_pointer_cast< tudat::gravitation::SphericalHarmonicsGravityField >(
                bodyExertingAcceleration->getGravityFieldModel( ) ) )
    {

        // Get characteristics of spherical harmonic gravity field.
        boost::shared_ptr< tudat::gravitation::SphericalHarmonicsGravityField >
                sphericalHarmonicField = boost::dynamic_pointer_cast<
                tudat::gravitation::SphericalHarmonicsGravityField >(
                    bodyExertingAcceleration->getGravityFieldModel( ) ) ;

        bodyGravitationalParameter = sphericalHarmonicField->getGravitationalParameter( );
        referenceRadius = sphericalHarmonicField->getReferenceRadius( );
        cosineCoefficients = sphericalHarmonicField->getCosineCoefficients( ).
                block( 0, 0, maximumDegree + 1, maximumOrder + 1 );
        sineCoefficients = sphericalHarmonicField->getSineCoefficients( ).
                block( 0, 0, maximumDegree + 1, maximumOrder + 1 );
    }
    else
    {
        std::cerr<<"Error when making sh gravity acceleration, did not find an sh "<<
                   "gravity field model"<<std::endl;
    }

    // Create functions to retrieve current state of bodies undergoing and exerting acceleration.
    boost::function< Eigen::Vector3d( ) > positionFunctionOfBodyExertingAcceleration =
            boost::bind( &getPositionVectorWrtBodyFixedFrame, bodyExertingAcceleration,
                         bodyExertingAcceleration );
    boost::function< Eigen::Vector3d( ) > positionFunctionOfBodyUndergoingAcceleration =
            boost::bind( &getPositionVectorWrtBodyFixedFrame, bodyUndergoingAcceleration,
                         bodyExertingAcceleration );

    // Create and return gravitational acceleration
    return boost::make_shared<
            tudat::gravitation::SphericalHarmonicsGravitationalAccelerationModelXd >(
                positionFunctionOfBodyUndergoingAcceleration, bodyGravitationalParameter,
                referenceRadius, cosineCoefficients, sineCoefficients,
                positionFunctionOfBodyExertingAcceleration );
}

//! Update a list of bodies to the current time and propagated state.
/*!
 *  Update a list of bodies to the current time and propagated state, by setting the current
 *  independent variable and state values, which are provided directly as input, and
 *  the dependent variables.
 *  \param bodyMap Named list of body objects.
 *  \param currentTime Current time, to be used for calculation of dependent variables.
 *  \param currentState To be set as state of body with name: nameOfPropagatedBody.
 *  \param nameOfPropagatedBody Name of body for which the currentState input value is the
 *  current state.
 *  \param propagatedFrameToEphemerisFrameFunction Function which calculates a vector
 *  (as a function of time) to translate the currentState to the frame in which the inertial state
 *  of the bodies is defined. That is, for a numerical integration in frame A, and an inertial
 *  frame B, this function should return the origin of A in frame B as a function of time.
 */
void updateBodies(
        const std::map< std::string, boost::shared_ptr< tudat_application::Body > >& bodyMap,
        const double currentTime, const tudat::basic_mathematics::Vector6d& currentState,
        const std::string nameOfPropagatedBody,
        const boost::function< tudat::basic_mathematics::Vector6d( const double ) >
        propagatedFrameToEphemerisFrameFunction =
        boost::lambda::constant( tudat::basic_mathematics::Vector6d::Zero( ) ) )
{
    // Iterate over all bodies in bodyMap
    for( std::map< std::string, boost::shared_ptr< tudat_application::Body > >::const_iterator
         bodyIterator = bodyMap.begin( ); bodyIterator != bodyMap.end( ); bodyIterator++ )
    {
        // Check if name of current body is name of propagated body; if so set its state from the
        // currentState passed to this function (and the time for any additional time-dependent
        // environment models).
        if( bodyIterator->first == nameOfPropagatedBody )
        {
            bodyIterator->second->setCurrentTimeAndState(
                        currentTime, currentState + propagatedFrameToEphemerisFrameFunction(
                            currentTime ) );
        }
        // Otherwise set the current state of the body from the environment models in the body
        else
        {
            bodyIterator->second->updateStateFromEphemeris( currentTime );
        }
    }
}

//! Function to create a state derivative model object from a list of acceleration models.
/*!
 *  Function to create a state derivative model object from a list of acceleration models.
 *  \param accelerationModels Vector of acceleration model objects.
 *  \param bodyMap Named list of body objects.
 *  \param nameOfPropagatedBody Name of body for the propagation of which a state derivative
 *  model is to be created.
 *  \param propagatedFrameToEphemerisFrameFunction Function which calculates a vector
 *  (as a function of time) to translate the currentState to the frame in which the inertial state
 *  of the bodies is defined. That is, for a numerical integration in frame A, and an inertial
 *  frame B, this function should return the origin of A in frame B as a function of time.
 *  \return State derivative model using accelerationModels to calculate state derivative
 *  of nameOfPropagatedBody.
 */
boost::shared_ptr< tudat::state_derivative_models::StateDerivativeModel
< double, tudat::basic_mathematics::Vector6d > >
createStateDerivativeModel(
        const std::vector< boost::shared_ptr< tudat::basic_astrodynamics::AccelerationModel3d > >&
        accelerationModels,
        const std::map< std::string, boost::shared_ptr< tudat_application::Body > >& bodyMap,
        const std::string& nameOfPropagatedBody,
        const boost::function< tudat::basic_mathematics::Vector6d( const double ) >
        propagatedFrameToEphemerisFrameFunction =
        boost::lambda::constant( tudat::basic_mathematics::Vector6d::Zero( ) ) )
{
    // Create function which will update the environmemnt to the current state and time during
    // the numerical integration
    boost::function< void( const double, const tudat::basic_mathematics::Vector6d& ) >
            updateFunction =
            boost::bind( &updateBodies, bodyMap, _1, _2, nameOfPropagatedBody,
                         propagatedFrameToEphemerisFrameFunction );

    // Create and return state derivative model
    return boost::make_shared< tudat::state_derivative_models::CartesianStateDerivativeModel
            < double, tudat::basic_mathematics::Vector6d > >(
                accelerationModels, updateFunction );
}

//! Function to create a state derivative model object from a list of acceleration models,
//! one of which is a spherical harmonic acceleration.
/*!
 *  Function to create a state derivative model object from a list of acceleration models,
 *  one of which is a spherical harmonic acceleration.
 *  \param accelerationModels Vector of acceleration model objects.
 *  \param bodyMap Named list of body objects.
 *  \param nameOfPropagatedBody Name of body for the propagation of which a state
 *  derivative model is to be created.
 *  \param bodyWithSphericalHarmonicField Name of the body for which the gravity field is
 *  used for the calculation of the spherical harmonic acceleration.
 *  \param propagatedFrameToEphemerisFrameFunction Function which calculates a vector
 *  (as a function of time) to translate the currentState to the frame in which the inertial state
 *  of the bodies is defined. That is, for a numerical integration in frame A, and an inertial
 *  frame B, this function should return the origin of A in frame B as a function of time.
 *  \return State derivative model using accelerationModels to calculate state derivative of
 *  nameOfPropagatedBody.
 */
boost::shared_ptr< tudat::state_derivative_models::StateDerivativeModel
< double, tudat::basic_mathematics::Vector6d > >
createStateDerivativeModel(
        const std::vector< boost::shared_ptr< tudat::basic_astrodynamics::AccelerationModel3d > >&
        accelerationModels,
        const std::map< std::string, boost::shared_ptr< tudat_application::Body > >& bodyMap,
        const std::string& nameOfPropagatedBody,
        const std::string bodyWithSphericalHarmonicField,
        const boost::function< tudat::basic_mathematics::Vector6d( const double ) >
        propagatedFrameToEphemerisFrameFunction =
        boost::lambda::constant( tudat::basic_mathematics::Vector6d::Zero( ) ) )
{
    // Create function which will update the environmemnt to the current state and time
    // during the numerical integration.
    boost::function< void( const double, const tudat::basic_mathematics::Vector6d& ) >
            updateFunction = boost::bind(
                &updateBodies, bodyMap, _1, _2, nameOfPropagatedBody,
                propagatedFrameToEphemerisFrameFunction );

    // Define reference frame transformation typedef.
    using tudat::basic_astrodynamics::AccelerationModel3d;
    typedef boost::function< Eigen::Vector3d ( const Eigen::Vector3d& ) >
            ReferenceFrameTransformationFunction;

    // Declare list of acceleration models with rotations.
    std::vector< std::pair< boost::shared_ptr< AccelerationModel3d >,
            std::vector< ReferenceFrameTransformationFunction > > >
            accelerationModelsWithRotations;

    // Define dummy for central and thrid body accelerations.
    std::vector< ReferenceFrameTransformationFunction > dummyRotations;
    dummyRotations.push_back(
                &tudat::state_derivative_models::transformNothing< Eigen::Vector3d > );

    // Iterate over full list of acceleration models and pair with associated rotation.
    for( unsigned int i = 0; i < accelerationModels.size( ); i++ )
    {
        // If not a spherical harmonic acceleration, pair with dummy rotation.
        if( boost::dynamic_pointer_cast
                < tudat::gravitation::SphericalHarmonicsGravitationalAccelerationModelXd >(
                    accelerationModels.at( i ) ) == NULL )
        {
            accelerationModelsWithRotations.push_back(
                        std::make_pair( accelerationModels.at( i ), dummyRotations ) );
        }
        // Else, pair acceleration with the transformation from body-fixed to inertial.
        else
        {
            std::vector< ReferenceFrameTransformationFunction > rotationFromBodyFixedFrame;
            rotationFromBodyFixedFrame.push_back(
                        boost::bind( &rotateVectorFromBodyFixedToInertialFrame, _1,
                                     bodyMap.at( bodyWithSphericalHarmonicField ) ) );
            accelerationModelsWithRotations.push_back(
                        std::make_pair( accelerationModels.at( i ), rotationFromBodyFixedFrame ) );

        }
    }

    // Create and return state derivative model
    return boost::make_shared< tudat::state_derivative_models::CartesianStateDerivativeModel
            < double, tudat::basic_mathematics::Vector6d > >(
                accelerationModelsWithRotations, updateFunction );
}

//! Function to create map of time as key and state as value, where the states describe a
//! Keplerian orbit about a body A, expressed in a frame with origin B (where A=B by default).
/*!
 *  Function to create map (table_ of time as key and state as value, where the states describe a
 *  fully Keplerian orbit.
 *  \param startTime Time from which to start the generation of the table.
 *  \param endTime Time from which to end the generation of the table.
 *  \param timeResolution Size of time step to take between succesive entries of the table.
 *  \param initialKeplerElements Kepler elements at startTime.
 *  \param centralBodyGravitationalParameter Gravitational parameter of the body A
 *  about which the Keplerian orbit is described.
 *  \param centralBodyStateFunction Function to determine the state of body A w.r.t. to the origin
 *  of the reference frame in which the states in the table are to be expressed.
 *  \return Table with states of a Keplerian orbit about a body A, expressed in a frame with
 *  origin B.
 */
std::map< double, tudat::basic_mathematics::Vector6d > createKeplerOrbitStateHistory(
        const double startTime, const double endTime, const double timeResolution,
        const tudat::basic_mathematics::Vector6d& initialKeplerElements,
        const double centralBodyGravitationalParameter,
        const boost::function< tudat::basic_mathematics::Vector6d( const double ) >
        centralBodyStateFunction )
{
    std::map< double, tudat::basic_mathematics::Vector6d > stateHistory;

    tudat::basic_mathematics::Vector6d currentKeplerElements;

    // Iterate each time step until end time reached.
    double currentTime = startTime;
    while( currentTime <= endTime )
    {
        // Calculate current Kepler elements
        currentKeplerElements = tudat::orbital_element_conversions::
                propagateKeplerOrbit(
                    initialKeplerElements, currentTime - startTime,
                    centralBodyGravitationalParameter );

        // Convert to Cartesian elements and set in state history.
        stateHistory[ currentTime ] = tudat::orbital_element_conversions::
                convertKeplerianToCartesianElements(
                    currentKeplerElements, centralBodyGravitationalParameter ) +
                centralBodyStateFunction( currentTime );

        // Update current time.
        currentTime += timeResolution;
    }

    // Return state history.
    return stateHistory;
}

//! Function to create map of time as key and state as value, where the states are obtained
//! from an existing Ephemeris object.
/*!
 *  Function to create map of time as key and state as value, where the states are obtained
 *  from an existing Ephemeris object. This function may be useful for writing the state
 *  history of a body to a file.
 *  \param bodyEphemeris Ephemeris from which the states are to be computed
 *  \param startTime Time from which to start the generation of the table.
 *  \param endTime Time from which to end the generation of the table.
 *  \param timeResolution Size of time step to take between succesive entries of the table.
 *  \return Table with states taken from bodyEphemeris at given time interval.
 */
std::map< double, tudat::basic_mathematics::Vector6d > createBodyStateHistoryFromEphemeris(
        const boost::shared_ptr< tudat::ephemerides::Ephemeris > bodyEphemeris,
        const double startTime, const double endTime, const double timeResolution )
{
    std::map< double, tudat::basic_mathematics::Vector6d > stateHistory;

    // Iterate each time step until end time reached.
    double currentTime = startTime;
    while( currentTime <= endTime )
    {
        // Get state and set in state history.
        stateHistory[ currentTime ] = bodyEphemeris->getCartesianStateFromEphemeris(
                    currentTime, tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000 );

        // Update current time.
        currentTime += timeResolution;
    }

    return stateHistory;
}

//! Execute propagation of LRO about the Moon
int main(int argc, char *argv[])
{
    using namespace tudat_application;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::ephemerides;
    using namespace tudat::gravitation;
    using namespace tudat::input_output;
    using namespace tudat::interpolators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::state_derivative_models;
    using namespace tudat::unit_conversions;
    using namespace tudat::spice_interface;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::mathematical_constants;


    // Load Spice kernels.
    loadSpiceKernelInTudat( getTudatRootPath( ) + "External/SpiceInterface/Kernels/pck00009.tpc" );
    loadSpiceKernelInTudat( getTudatRootPath( ) +
                            "External/SpiceInterface/Kernels/de-403-masses.tpc" );
    loadSpiceKernelInTudat( getTudatRootPath( ) + "External/SpiceInterface/Kernels/de421.bsp" );
    loadSpiceKernelInTudat( getTudatRootPath( ) + "External/SpiceInterface/Kernels/naif0009.tls" );


    // Create solar system bodies with ephemeris and gravitational parameters from Spice
    std::map< std::string, boost::shared_ptr< Body > > bodyMap;
    bodyMap[ "Earth" ] = boost::make_shared< Body >( );
    bodyMap[ "Earth" ]->setEphemeris( boost::make_shared< SpiceEphemeris >(
                                          "Earth", "SSB", false, false, false, "J2000" ) );
    bodyMap[ "Earth" ]->setGravityFieldModel ( boost::make_shared< GravityFieldModel >(
                                                   getBodyGravitationalParameter( "Earth" ) ) );
    bodyMap[ "Moon" ] = boost::make_shared< Body >( );
    bodyMap[ "Moon" ]->setEphemeris( boost::make_shared< SpiceEphemeris >(
                                         "Moon", "SSB", false, false, false, "J2000" ) );
    bodyMap[ "Moon" ]->setGravityFieldModel ( boost::make_shared< GravityFieldModel >(
                                                  getBodyGravitationalParameter( "Moon" ) ) );
    bodyMap[ "Sun" ] = boost::make_shared< Body >( );
    bodyMap[ "Sun" ]->setEphemeris( boost::make_shared< SpiceEphemeris >(
                                        "Sun", "SSB", false, false, false, "J2000" ) );
    bodyMap[ "Sun" ]->setGravityFieldModel ( boost::make_shared< GravityFieldModel >(
                                                 getBodyGravitationalParameter( "Sun" ) ) );

    // Create body for LRO with 1000.0 kg mass
    bodyMap[ "LRO" ] = boost::make_shared< Body >( Vector6d::Zero( ), 0.0, 1000.0 );

    // Set initial Kepler elements for LRO.
    tudat::basic_mathematics::Vector6d initialKeplerElements;
    initialKeplerElements[ semiMajorAxisIndex ] = 1900.0E3;
    initialKeplerElements[ eccentricityIndex ] = 0.05;
    initialKeplerElements[ inclinationIndex ] =
            89.7 * PI / 180.0;

    ////////////////////////////////////////////////////////////////////////////////////////
    ///////  CHANGE THESE VALUES ACCORDING TO YOUR STUDENT NUMBER (SEE ASSIGNMENT)  ////////
    ////////////////////////////////////////////////////////////////////////////////////////
    initialKeplerElements[ argumentOfPeriapsisIndex ] = 0.0;
    initialKeplerElements[ longitudeOfAscendingNodeIndex ] = 0.0;
    initialKeplerElements[ trueAnomalyIndex ] = 0.0;

    // Create Keplerian orbit state history for LRO
    std::map< double, tudat::basic_mathematics::Vector6d > lroStateHistory =
            createKeplerOrbitStateHistory(
                -1000.0, 1.0E6 + 1000, 30.0, initialKeplerElements,
                bodyMap[ "Moon" ]->getGravityFieldModel( )->getGravitationalParameter( ),
                boost::bind( &Ephemeris::getCartesianStateFromEphemeris,
                             bodyMap[ "Moon" ]->getEphemeris( ), _1, JULIAN_DAY_ON_J2000 ) );

    // Create orbit interpolator for LRO and set as Ephemeris
    boost::shared_ptr< LagrangeInterpolator< double, Vector6d > > lagrangeInterpolator =
            boost::make_shared< LagrangeInterpolator< double, Vector6d > >( lroStateHistory, 8 );
    boost::shared_ptr< TabulatedCartesianEphemeris > lroEphemeris =
            boost::make_shared< TabulatedCartesianEphemeris >( lagrangeInterpolator );
    bodyMap[ "LRO" ]->setEphemeris( lroEphemeris );
}
