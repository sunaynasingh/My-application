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
 *      YYMMDD    Author            Comment
 *      120824    P. Musegaas       First creation of code.
 *      120914    P. Musegaas       Update after discussion on course.
 *      121218    K. Kumar          Added option to set output directory (defaults to
 *                                  root-directory).
 *
 *    References
 *
 *    Notes
 *      Part (f) of the assignment: performing a 4D grid search.
 *
 */

#include <cmath>
#include <fstream>
#include <iostream>

#include <Eigen/Core>

#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include "Tudat/Astrodynamics/MissionSegments/lambertRoutines.h"
#include "Tudat/Astrodynamics/MissionSegments/escapeAndCapture.h"

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

#include "SampleReturnMission/Ephemeris/getEphemeris.h"

#include "SampleReturnMission/basicInputOutput.h"
#include "SampleReturnMission/missionDefinitions.h"

// Declare that we want to use the namespace belonging to this course and the tudat namespace.
// NEVER do this in a header file!
using namespace tudat_course::sample_return_mission;
using namespace tudat;

int main( )
{
    // Set directory where output files will be stored. By default, this is your project
    // root-directory.
    const std::string outputDirectory = basic_input_output::getApplicationRootPath( ) + "/";

    // Set bounds for the departure dates from Earth to the target planet, its step-size and number
    // of steps. Notice how an integer is cast into a double.
    const double lowerBoundDepartureDate = 8400.0;
    const double upperBoundDepartureDate = 9860.0;
    const unsigned int numberOfStepsDepartureDate = 21;
    const double stepSizeDepartureDate = ( upperBoundDepartureDate - lowerBoundDepartureDate ) /
            static_cast< double >( numberOfStepsDepartureDate - 1 );

    // Set bounds for the transfer time from Earth to the target planet, its step-size and number
    // of steps.
    const double lowerBoundTransferTime = 100.0;
    const double upperBoundTransferTime = 400.0;
    const unsigned int numberOfStepsTransferTime = 21;
    const double stepSizeTransferTime = ( upperBoundTransferTime - lowerBoundTransferTime ) /
            static_cast< double >( numberOfStepsTransferTime - 1 );

    // Set bounds for the stay time at the target planet, its step-size and number of steps.
    const double lowerBoundStayTime = 200.0;
    const double upperBoundStayTime = 500.0;
    const unsigned int numberOfStepsStayTime = 21;
    const double stepSizeStayTime = ( upperBoundStayTime - lowerBoundStayTime ) /
            static_cast< double >( numberOfStepsStayTime - 1 );

    // Set bounds for the return time from the target planet to Earth, its step-size and number
    // of steps.
    const double lowerBoundReturnTime = 100.0;
    const double upperBoundReturnTime = 400.0;
    const unsigned int numberOfStepsReturnTime = 21;
    const double stepSizeReturnTime = ( upperBoundReturnTime - lowerBoundReturnTime ) /
            static_cast< double >( numberOfStepsReturnTime - 1 );

    // Calculate the total number of samples.
    const int totalNumberOfSamples = numberOfStepsDepartureDate * numberOfStepsTransferTime *
            numberOfStepsStayTime * numberOfStepsReturnTime;

    // Initialize the departure dates vector, the transfer times vector, the stay times vector and
    // the return times vector.
    Eigen::VectorXd departureDates( totalNumberOfSamples ),
            transferTimes( totalNumberOfSamples ),
            stayTimes( totalNumberOfSamples ),
            returnTimes( totalNumberOfSamples );

    // Fill the vectors with values.
    unsigned int totalCounter = 0;
    for ( unsigned int counter1 = 0; counter1 < numberOfStepsDepartureDate; counter1++ )
    {
        for ( unsigned int counter2 = 0; counter2 < numberOfStepsTransferTime; counter2++ )
        {
            for ( unsigned int counter3 = 0; counter3 < numberOfStepsStayTime; counter3++ )
            {
                for ( unsigned int counter4 = 0; counter4 < numberOfStepsReturnTime; counter4++ )
                {
                    departureDates( totalCounter ) =
                            lowerBoundDepartureDate + counter1 * stepSizeDepartureDate;
                    transferTimes( totalCounter ) =
                            lowerBoundTransferTime + counter2 * stepSizeTransferTime;
                    stayTimes( totalCounter ) = lowerBoundStayTime + counter3 * stepSizeStayTime;
                    returnTimes( totalCounter ) =
                            lowerBoundReturnTime + counter4 * stepSizeReturnTime;
                    totalCounter++;
                }
            }
        }
    }

    // Initialize a matrix that will keep track of the deltaVs.
    Eigen::VectorXd deltaVs( totalNumberOfSamples );

    // Calculate the trajectory for all the dates.
    for ( int counter = 0; counter < totalNumberOfSamples; counter++ )
    {
        // Store all the relevant dates from the previously filled vectors.
        const double earthDepartureDateInMJD2000 = departureDates( counter );
        const double transferTimeInDays = transferTimes( counter );
        const double targetArrivalDateInMJD2000 = earthDepartureDateInMJD2000 + transferTimeInDays;
        const double targetDepartureDateInMJD2000 =
                targetArrivalDateInMJD2000 + stayTimes( counter );
        const double returnTimeInDays = returnTimes( counter );
        const double earthArrivalDateInMJD2000 = targetDepartureDateInMJD2000 + returnTimeInDays;

        // Obtain Keplerian elements of the planets.
        const basic_mathematics::Vector6d earthDepartureKeplerianElements =
                ephemeris::getKeplerianElementsOfPlanet( mission_definitions::ephemerisEarth,
                                                         earthDepartureDateInMJD2000 );
        const basic_mathematics::Vector6d targetArrivalKeplerianElements =
                ephemeris::getKeplerianElementsOfPlanet( mission_definitions::ephemerisTarget,
                                                         targetArrivalDateInMJD2000 );
        const basic_mathematics::Vector6d targetDepartureKeplerianElements =
                ephemeris::getKeplerianElementsOfPlanet( mission_definitions::ephemerisTarget,
                                                         targetDepartureDateInMJD2000 );
        const basic_mathematics::Vector6d earthArrivalKeplerianElements =
                ephemeris::getKeplerianElementsOfPlanet( mission_definitions::ephemerisEarth,
                                                         earthArrivalDateInMJD2000 );

        // Convert Keplerian elements into Cartesian position and velocity of the planets.
        const basic_mathematics::Vector6d earthDepartureCartesianElements =
                orbital_element_conversions::convertKeplerianToCartesianElements(
                    earthDepartureKeplerianElements,
                    mission_definitions::gravitationalParameterSun );
        const basic_mathematics::Vector6d targetArrivalCartesianElements =
                orbital_element_conversions::convertKeplerianToCartesianElements(
                    targetArrivalKeplerianElements,
                    mission_definitions::gravitationalParameterSun );
        const basic_mathematics::Vector6d targetDepartureCartesianElements =
                orbital_element_conversions::convertKeplerianToCartesianElements(
                    targetDepartureKeplerianElements,
                    mission_definitions::gravitationalParameterSun );
        const basic_mathematics::Vector6d earthArrivalCartesianElements =
                orbital_element_conversions::convertKeplerianToCartesianElements(
                    earthArrivalKeplerianElements,
                    mission_definitions::gravitationalParameterSun );

        // Split the cartesian elements into position and velocity.
        const Eigen::Vector3d earthDeparturePosition =
                earthDepartureCartesianElements.segment( 0, 3 );
        const Eigen::Vector3d earthDepartureVelocity =
                earthDepartureCartesianElements.segment( 3, 3 );
        const Eigen::Vector3d targetArrivalPosition =
                targetArrivalCartesianElements.segment( 0, 3 );
        const Eigen::Vector3d targetArrivalVelocity =
                targetArrivalCartesianElements.segment( 3, 3 );
        const Eigen::Vector3d targetDeparturePosition =
                targetDepartureCartesianElements.segment( 0, 3 );
        const Eigen::Vector3d targetDepartureVelocity =
                targetDepartureCartesianElements.segment( 3, 3 );
        const Eigen::Vector3d earthArrivalPosition =
                earthArrivalCartesianElements.segment( 0, 3 );
        const Eigen::Vector3d earthArrivalVelocity =
                earthArrivalCartesianElements.segment( 3, 3 );

        // Convert transfer time and return time to seconds.
        const double transferTime =
                unit_conversions::convertJulianDaysToSeconds( transferTimeInDays );
        const double returnTime =
                unit_conversions::convertJulianDaysToSeconds( returnTimeInDays );

        // Initialize the vectors that will contain the spacecraft velocities.
        Eigen::Vector3d velocityAtEarthDeparture, velocityAtTargetArrival,
                velocityAtTargetDeparture, velocityAtEarthArrival;

        // Compute these velocities using the Lambert targeter.
        mission_segments::solveLambertProblemIzzo(
                    earthDeparturePosition, targetArrivalPosition, transferTime,
                    mission_definitions::gravitationalParameterSun,
                    velocityAtEarthDeparture, velocityAtTargetArrival );
        mission_segments::solveLambertProblemIzzo(
                    targetDeparturePosition, earthArrivalPosition, returnTime,
                    mission_definitions::gravitationalParameterSun,
                    velocityAtTargetDeparture, velocityAtEarthArrival );

        // Compute the velocity increments at Earth and the target.
        const Eigen::Vector3d velocityIncrementAtEarthDeparture =
                velocityAtEarthDeparture - earthDepartureVelocity;
        const Eigen::Vector3d velocityIncrementAtTargetArrival =
                velocityAtTargetArrival - targetArrivalVelocity;
        const Eigen::Vector3d velocityIncrementAtTargetDeparture =
                velocityAtTargetDeparture - targetDepartureVelocity;
        const Eigen::Vector3d velocityIncrementAtEarthArrival =
                velocityAtEarthArrival - earthArrivalVelocity;

        // Compute the excess velocities at Earth and the target.
        const double excessVelocityAtEarthDeparture = velocityIncrementAtEarthDeparture.norm( );
        const double excessVelocityAtTargetArrival = velocityIncrementAtTargetArrival.norm( );
        const double excessVelocityAtTargetDeparture = velocityIncrementAtTargetDeparture.norm( );
        const double excessVelocityAtEarthArrival = velocityIncrementAtEarthArrival.norm( );

        // Compute the deltaV for the escape and capture maneuvers.
        const double deltaVatEarthDeparture = mission_segments::computeEscapeOrCaptureDeltaV(
                    mission_definitions::gravitationalParameterEarth,
                    mission_definitions::parkingOrbitSemiMajorAxis,
                    mission_definitions::parkingOrbitEccentricity,
                    excessVelocityAtEarthDeparture );
        const double deltaVAtTargetArrival = mission_segments::computeEscapeOrCaptureDeltaV(
                    mission_definitions::gravitationalParameterTarget,
                    mission_definitions::captureOrbitSemiMajorAxis,
                    mission_definitions::captureOrbitEccentricity,
                    excessVelocityAtTargetArrival );
        const double deltaVAtTargetDeparture = mission_segments::computeEscapeOrCaptureDeltaV(
                    mission_definitions::gravitationalParameterTarget,
                    mission_definitions::captureOrbitSemiMajorAxis,
                    mission_definitions::captureOrbitEccentricity,
                    excessVelocityAtTargetDeparture );
        const double deltaVatEarthArrival = mission_segments::computeEscapeOrCaptureDeltaV(
                    mission_definitions::gravitationalParameterEarth,
                    mission_definitions::earthReturnOrbitSemiMajorAxis,
                    mission_definitions::earthReturnOrbitEccentricity,
                    excessVelocityAtEarthArrival );

        // Compute the total deltaV and store it in the deltaV vector.
        deltaVs( counter ) = deltaVatEarthDeparture + deltaVAtTargetArrival +
                deltaVAtTargetDeparture + deltaVatEarthArrival;
    }

    // Write results to one comma separated file. First column is departure dates, second column is
    // the transfer times, third column is the stay times, fourth colum is the return times and
    // finally the fifth column is the deltaV values.

    // Combine all the vectors in one matrix.
    Eigen::MatrixXd results( totalNumberOfSamples, 5 );
    results << departureDates, transferTimes, stayTimes, returnTimes, deltaVs;

    // Specify the output format, such that Matlab can read it with csvread.
    Eigen::IOFormat csvFormat( 15, 0, ", ", "\n" );

    // Set absolute path to file containing data.
    const std::string outputFileAbsolutePath = outputDirectory + "ResultsPartF.csv";

    // Export the matrix.
    std::ofstream exportFile( outputFileAbsolutePath.c_str( ) );
    exportFile << results.format( csvFormat );
    exportFile.close( );

    return 0;
}
