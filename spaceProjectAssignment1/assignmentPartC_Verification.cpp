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
 *      120824    P. Musegaas       File created.
 *      120914    P. Musegaas       Update after discussion on course.
 *
 *    References
 *
 *    Notes
 *      Part (c) of the assignment: Verifying the Lambert targeter and the escape and capture
 *      maneuver.
 *
 */

#include <cmath>
#include <iostream>

#include <Eigen/Core>

#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include "Tudat/Astrodynamics/MissionSegments/lambertRoutines.h"
#include "Tudat/Astrodynamics/MissionSegments/escapeAndCapture.h"

#include "SampleReturnMission/missionDefinitions.h"

// Declare that we want to use the namespace belonging to this course and the tudat namespace.
// NEVER do this in a header file!
using namespace tudat_course::sample_return_mission;
using namespace tudat;

int main( )
{
    /////////      Calculate the Hohmann transfer delta V.    ///////////

    // Calculate the orbital velocities of Earth and the target planet.
    const double orbitalVelocityEarth = std::sqrt( mission_definitions::gravitationalParameterSun /
                                                   mission_definitions::orbitalRadiusEarth );
    const double orbitalVelocityTarget
            = std::sqrt( mission_definitions::gravitationalParameterSun /
                         mission_definitions::orbitalRadiusTarget );

    // Calculate the required heliocentric velocity at Earth and the target for Hohmann.
    // Note the .0's in the divisions. These are very important, because otherwise the numbers
    // are treated as integers. For example, 3 / 2 = 1! Note that 2. is equal to 2.0.
    const double semiMajorAxisTransferOrbit = ( mission_definitions::orbitalRadiusEarth +
                                                mission_definitions::orbitalRadiusTarget ) / 2.0;
    const double hohmannVelocityAtEarth =
            std::sqrt( mission_definitions::gravitationalParameterSun *
                       ( 2.0 / mission_definitions::orbitalRadiusEarth -
                         1.0 / semiMajorAxisTransferOrbit ) );
    const double hohmannVelocityAtTarget =
            std::sqrt( mission_definitions::gravitationalParameterSun *
                       ( 2.0 / mission_definitions::orbitalRadiusTarget -
                         1.0 / semiMajorAxisTransferOrbit ) );

    // Calculate the excess velocities.
    const double hohmannExcessEarth = std::fabs( hohmannVelocityAtEarth - orbitalVelocityEarth );
    const double hohmannExcessTarget = std::fabs(
                hohmannVelocityAtTarget - orbitalVelocityTarget );
    const double hohmannExcessSummed = hohmannExcessEarth + hohmannExcessTarget;

    // Print the result.
    std::cout << "The Hohmann excess velocity at Earth is: " << hohmannExcessEarth << " m/s."
              << std::endl;
    std::cout << "The Hohmann excess velocity at the target is: " << hohmannExcessTarget << " m/s."
              << std::endl;
    std::cout << "The Hohmann excess velocities summed yield: " << hohmannExcessSummed << " m/s."
              << std::endl << std::endl;




    ////////////       Calculate Hohmann transfer using the Lambert targeter.    ///////////

    // Compute the Hohmann time of flight for a mission to the target in seconds.
    const double hohmannTimeOfFlight =
            tudat::mathematical_constants::PI
            * std::sqrt( semiMajorAxisTransferOrbit *
                         semiMajorAxisTransferOrbit *
                         semiMajorAxisTransferOrbit /
                         mission_definitions::gravitationalParameterSun );

    // Fill the initial position and final position with that of circular coplanar positions of the
    // Earth and the target. Note that an offset is generated in the Earth's position of
    // 0.001 meters, such that the Lambert targeter can determine the plane of the orbit.
    const Eigen::Vector3d earthPosition( mission_definitions::orbitalRadiusEarth, 0.001, 0. );
    const Eigen::Vector3d targetPosition(
                -1.0 * mission_definitions::orbitalRadiusTarget, 0., 0. );

    // Initialize variables that will store the initial and final velocities.
    Eigen::Vector3d velocityAtEarth, velocityAtTarget;

    // Call the Lambert targeter. Izzo's method is used, because it is a stable and fast Lambert
    // targeter.
    mission_segments::solveLambertProblemIzzo( earthPosition, targetPosition,
                                               hohmannTimeOfFlight,
                                               mission_definitions::gravitationalParameterSun,
                                               velocityAtEarth, velocityAtTarget );

    // Fill vectors containing the velocities of the Earth and the target.
    const Eigen::Vector3d earthVelocity( 0.0, orbitalVelocityEarth, 0.0 );
    const Eigen::Vector3d targetVelocity( 0.0, -1.0 * orbitalVelocityTarget, 0.0 );

    // Compute the 3D excess velocities required.
    const Eigen::Vector3d lambertExcessEarth = velocityAtEarth - earthVelocity;
    const Eigen::Vector3d lambertExcessTarget = velocityAtTarget - targetVelocity;

    // Compute the sum.
    const double lambertExcessSummed = lambertExcessEarth.norm( ) + lambertExcessTarget.norm( );

    // Print the result.
    std::cout << "The excess velocity at Earth calculated with Tudat is: "
              << lambertExcessEarth.norm( ) << " m/s." << std::endl;
    std::cout << "The excess velocity at the target calculated with Tudat is: "
              << lambertExcessTarget.norm( ) << " m/s." << std::endl;
    std::cout << "The excess velocities calculated with Tudat yield summed: "
              << lambertExcessSummed << " m/s." << std::endl << std::endl;




    ////////////      'Manually' calculate the escape and capture maneuver costs.    ///////////

    // Compute the velocity in the parking orbit.
    const double velocityEarthParkingOrbit =
            std::sqrt( mission_definitions::gravitationalParameterEarth /
                       mission_definitions::parkingOrbitSemiMajorAxis );

    // Compute the escape velocity.
    const double escapeVelocityEarthParking =
            std::sqrt( 2.0 * mission_definitions::gravitationalParameterEarth /
                       mission_definitions::parkingOrbitSemiMajorAxis );

    // Compute the deltaV for escaping.
    const double escapeEarthDeltaV = escapeVelocityEarthParking - velocityEarthParkingOrbit;

    // Print the result.
    std::cout << "To escape Earth orbit from the parking orbit, a deltaV of "
              << escapeEarthDeltaV << " m/s is required." << std::endl;

    // Compute excess velocity required.
    const double excessVelocityEarth = std::abs( hohmannVelocityAtEarth - orbitalVelocityEarth );

    // Compute deltaV to achieve that.
    const double velocityForTransferToTarget =
            std::sqrt( excessVelocityEarth * excessVelocityEarth +
                       escapeVelocityEarthParking * escapeVelocityEarthParking );
    const double deltaVforTransferToTarget = velocityForTransferToTarget -
            velocityEarthParkingOrbit;

    // Print the result.
    std::cout << "To escape Earth orbit to target transfer orbit from the parking orbit, a deltaV"
              << " of " << deltaVforTransferToTarget << " is required." << std::endl;

    // Calculate the pericenter radius of the capture orbit.
    const double pericenterRadiusTargetCaptureOrbit =
            mission_definitions::captureOrbitSemiMajorAxis *
            ( 1.0 - mission_definitions::captureOrbitEccentricity );

    // Calculate the velocity at pericenter.
    const double pericenterVelocityTargetCaptureOrbit =
            std::sqrt( mission_definitions::gravitationalParameterTarget /
                       mission_definitions::captureOrbitSemiMajorAxis *
                       ( 1.0 + mission_definitions::captureOrbitEccentricity ) /
                       ( 1.0 - mission_definitions::captureOrbitEccentricity ) );

    // Calculate the escape velocity from the target capture orbit.
    const double escapeVelocityTargetCapture =
            std::sqrt( 2.0 * mission_definitions::gravitationalParameterTarget /
                       pericenterRadiusTargetCaptureOrbit );

    // Calculate escape deltaV from target capture orbit.
    const double escapeTargetDeltaV =
            escapeVelocityTargetCapture - pericenterVelocityTargetCaptureOrbit;

    // Print the result.
    std::cout << "To escape the target from the capture orbit, a deltaV of "
              << escapeTargetDeltaV << " m/s is required." << std::endl;

    // Compute excess velocity required.
    const double excessVelocityTarget = std::fabs(
                hohmannVelocityAtTarget - orbitalVelocityTarget );

    // Compute deltaV to achieve that.
    const double velocityForTransferFromEarth =
            std::sqrt( excessVelocityTarget * excessVelocityTarget +
                       escapeVelocityTargetCapture * escapeVelocityTargetCapture );
    const double deltaVforTransferFromEarth =
            velocityForTransferFromEarth - pericenterVelocityTargetCaptureOrbit;

    // Print the result.
    std::cout << "To arrive in the target transfer orbit from the Earth, a deltaV"
              << " of " << deltaVforTransferFromEarth << " m/s is required." << std::endl
              << std::endl;




    ////////////       Calculate the escape and capture maneuvers using Tudat.     ///////////
    const double earthEscapeDeltaVusingTudat = mission_segments::computeEscapeOrCaptureDeltaV(
                mission_definitions::gravitationalParameterEarth,
                mission_definitions::parkingOrbitSemiMajorAxis,
                mission_definitions::parkingOrbitEccentricity, 0.0 );
    const double earthTransferDeltaVusingTudat = mission_segments::computeEscapeOrCaptureDeltaV(
                mission_definitions::gravitationalParameterEarth,
                mission_definitions::parkingOrbitSemiMajorAxis,
                mission_definitions::parkingOrbitEccentricity, excessVelocityEarth );
    const double targetEscapeDeltaVusingTudat = mission_segments::computeEscapeOrCaptureDeltaV(
                mission_definitions::gravitationalParameterTarget,
                mission_definitions::captureOrbitSemiMajorAxis,
                mission_definitions::captureOrbitEccentricity, 0.0 );
    const double targetTransferDeltaVusingTudat = mission_segments::computeEscapeOrCaptureDeltaV(
                mission_definitions::gravitationalParameterTarget,
                mission_definitions::captureOrbitSemiMajorAxis,
                mission_definitions::captureOrbitEccentricity, excessVelocityTarget );

    // Print the results.
    std::cout << "Escape Earth orbit calculated using Tudat yields a deltaV of "
              << earthEscapeDeltaVusingTudat << " m/s." << std::endl;
    std::cout << "Escape Earth orbit to the target transfer calculated using Tudat yields a "
              << "deltaV of " << earthTransferDeltaVusingTudat << " m/s." << std::endl;
    std::cout << "Escape the target calculated using Tudat yields a deltaV of "
              << targetEscapeDeltaVusingTudat << " m/s." << std::endl;
    std::cout << "Capture in target orbit from Earth transfer calculated using Tudat yields a "
              << "deltaV of " << targetTransferDeltaVusingTudat << " m/s." << std::endl
              << std::endl;



    // Total deltaV from Earth parking orbit to the target capture orbit.
    const double totalDeltaVHand = deltaVforTransferToTarget + deltaVforTransferFromEarth;
    const double totalDeltaVTudat = earthTransferDeltaVusingTudat + targetTransferDeltaVusingTudat;

    std::cout << "Total deltaV for the mission from Earth by hand is " << totalDeltaVHand
              << " m/s." << std::endl;
    std::cout << "Total deltaV for the mission from Earth using Tudat is " << totalDeltaVTudat
              << " m/s." << std::endl;

    return 0;
}
