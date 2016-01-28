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
 *      Parts (d) and (e) of the assignment: planning a 2D and 3D mission to the target planet.
 *
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include <Eigen/Core>

#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include "Tudat/Astrodynamics/MissionSegments/lambertRoutines.h"
#include "Tudat/Astrodynamics/MissionSegments/escapeAndCapture.h"

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

#include "spaceProjectAssignment1/Ephemeris/getEphemeris.h"

#include "spaceProjectAssignment1/basicInputOutput.h"
#include "spaceProjectAssignment1/missionDefinitions.h"

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
    const unsigned int numberOfStepsDepartureDate = 1461;
    const double stepSizeDepartureDate = ( upperBoundDepartureDate - lowerBoundDepartureDate ) /
            static_cast< double >( numberOfStepsDepartureDate - 1 );

    // Set bounds for the transfer time from Earth to the target planet, its step-size and number
    // of steps.
    const double lowerBoundTransferTime = 100.0;
    const double upperBoundTransferTime = 400.0;
    const unsigned int numberOfStepsTransferTime = 301;
    const double stepSizeTransferTime = ( upperBoundTransferTime - lowerBoundTransferTime ) /
            static_cast< double >( numberOfStepsTransferTime - 1 );

    // Fill the departure dates vector and the transfer times vector.
    // Note that you can use the variable "counter" in both for-loops, since it is in local scope.
    Eigen::VectorXd departureDates( numberOfStepsDepartureDate ),
            transferTimes( numberOfStepsTransferTime );
    for ( unsigned int counter = 0; counter < numberOfStepsDepartureDate; counter++ )
    {
        departureDates( counter ) = lowerBoundDepartureDate + counter * stepSizeDepartureDate;
    }
    for ( unsigned int counter = 0; counter < numberOfStepsTransferTime; counter++ )
    {
        transferTimes( counter ) = lowerBoundTransferTime + counter * stepSizeTransferTime;
    }

    // Initialize a matrix that will keep track of the deltaVs.
    Eigen::MatrixXd deltaVs( numberOfStepsTransferTime, numberOfStepsDepartureDate );

    // Start the grid search.
    // Loop through all the transfer times.
    // Note in this case that because there are nested loops, the counter variables must be named
    // differently!
    for ( unsigned int counter1 = 0; counter1 < numberOfStepsTransferTime; counter1++ )
    {
        // Set the transfer time.
        const double transferTimeInDays = transferTimes( counter1 );

        // Loop through all the departure dates.
        for( unsigned int counter2 = 0; counter2 < numberOfStepsDepartureDate; counter2++ )
        {
            // Set the departure date.
            const double departureDateInMJD2000 = departureDates( counter2 );

            // Calculate the arrival date.
            const double arrivalDateInMJD2000 = departureDateInMJD2000 + transferTimeInDays;

            // Obtain Keplerian elements of the planets.
            const basic_mathematics::Vector6d earthKeplerianElements =
                    ephemeris::getKeplerianElementsOfPlanet(
                        mission_definitions::ephemerisEarth, departureDateInMJD2000 );
            const basic_mathematics::Vector6d targetKeplerianElements
                    = ephemeris::getKeplerianElementsOfPlanet(
                        mission_definitions::ephemerisTarget, arrivalDateInMJD2000 );

            // Convert Keplerian elements into Cartesian position and velocity of the planets.
            const basic_mathematics::Vector6d earthCartesianElements =
                    orbital_element_conversions::convertKeplerianToCartesianElements(
                        earthKeplerianElements, mission_definitions::gravitationalParameterSun );
            const basic_mathematics::Vector6d targetCartesianElements =
                    orbital_element_conversions::convertKeplerianToCartesianElements(
                        targetKeplerianElements, mission_definitions::gravitationalParameterSun );

            // Split the cartesian elements into position and velocity.
            const Eigen::Vector3d earthPosition = earthCartesianElements.segment( 0, 3 );
            const Eigen::Vector3d earthVelocity = earthCartesianElements.segment( 3, 3 );
            const Eigen::Vector3d targetPosition = targetCartesianElements.segment( 0, 3 );
            const Eigen::Vector3d targetVelocity = targetCartesianElements.segment( 3, 3 );

            // Convert transfer time to seconds.
            const double transferTime = tudat::
                    unit_conversions::convertJulianDaysToSeconds( transferTimeInDays );

            // Initialize the vectors that will contain the spacecraft velocities.
            Eigen::Vector3d velocityAtEarth, velocityAtTarget;

            // Compute these velocities using the Lambert targeter.
            mission_segments::solveLambertProblemIzzo(
                        earthPosition, targetPosition, transferTime,
                        mission_definitions::gravitationalParameterSun,
                        velocityAtEarth, velocityAtTarget );

            // Compute the velocity increments at Earth and the target.
            const Eigen::Vector3d velocityIncrementAtEarth = velocityAtEarth - earthVelocity;
            const Eigen::Vector3d velocityIncrementAtTarget = velocityAtTarget - targetVelocity;

            // Compute the excess velocities at Earth and the target.
            const double excessVelocityAtEarth = velocityIncrementAtEarth.norm( );
            const double excessVelocityAtTarget = velocityIncrementAtTarget.norm( );

            // Compute the deltaV for the escape and capture maneuvers.
            const double deltaVatEarth = mission_segments::computeEscapeOrCaptureDeltaV(
                        mission_definitions::gravitationalParameterEarth,
                        mission_definitions::parkingOrbitSemiMajorAxis,
                        mission_definitions::parkingOrbitEccentricity, excessVelocityAtEarth );
            const double deltaVAtTarget = mission_segments::computeEscapeOrCaptureDeltaV(
                        mission_definitions::gravitationalParameterTarget,
                        mission_definitions::captureOrbitSemiMajorAxis,
                        mission_definitions::captureOrbitEccentricity, excessVelocityAtTarget );

            // Compute the total deltaV and store it in the matrix.
            deltaVs( counter1, counter2 ) = deltaVatEarth + deltaVAtTarget;
        }
    }

    // Write the results to comma separated files. These are ready to be read in by MATLAB using
    // the csvread() function and the matrix can be plotted directly as a mesh.

    // Set output format for matrix output.
    Eigen::IOFormat csvFormat( 15, 0, ", ", "\n" );

    // Set absolute path to file containing deltaVs.
    const std::string deltaVOutputFileAbsolutePath = outputDirectory + "DeltaVsPartD.csv";

    // Export the deltaV matrix.
    std::ofstream exportFile1( deltaVOutputFileAbsolutePath.c_str( ) );
    exportFile1 << deltaVs.format( csvFormat );
    exportFile1.close( );

    // Set absolute path to file containing departure dates.
    const std::string departureDateOutputFileAbsolutePath
            = outputDirectory + "DepartureDatesPartD.csv";

    // Export the departure dates vector.
    std::ofstream exportFile2( departureDateOutputFileAbsolutePath.c_str( ) );
    exportFile2.precision( 15 );
    exportFile2 << departureDates;
    exportFile2.close( );

    // Set absolute path to file containing transfer times.
    const std::string transferTimesOutputFileAbsolutePath
            = outputDirectory + "TransferTimesPartD.csv";

    // Export the transfer times vector.
    std::ofstream exportFile3( transferTimesOutputFileAbsolutePath.c_str( ) );
    exportFile3.precision( 15 );
    exportFile3 << transferTimes;
    exportFile3.close( );

    return 0;
}
