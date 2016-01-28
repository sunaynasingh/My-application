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
 *      120914    P. Musegaas       Adaptation of own code for course. Temporary file.
 *      130225    K. Kumar          Removed tudat namespace; updated tudat references.
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the Major Planets,
 *          http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf, last accessed: 24 February, 2011.
 *
 *    Notes
 *      This file is to be replaced once there is a better version of ephemeris available in Tudat.
 *
 */

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h>

#include <Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

#include "spaceProjectAssignment1/Ephemeris/approximatePlanets.h"

namespace tudat_course
{
namespace sample_return_mission
{
namespace ephemeris
{

//! Get Keplerian Elements.
Eigen::VectorXd EphemerisApproximatePlanetPositionsInner::
        getKeplerianElements( const double timeInMJD2000 )
{
    // Declare some namespaces that we want to use.
    using namespace tudat::unit_conversions;
    using namespace tudat::orbital_element_conversions;

    // Initialize the keplerian state vector.
    tudat::basic_mathematics::Vector6d keplerianState;

    // Calculate the number of centuries after MJD 2000.
    const double centuriesAfterMJD2000 = timeInMJD2000 / 36525.0;

    // Calculate the semi major axis, eccentricity, inclination and longitude of ascending node.
    // Note that some variables have to be converted to SI units.
    keplerianState( semiMajorAxisIndex ) = convertAstronomicalUnitsToMeters(
            semiMajorAxis + centuriesAfterMJD2000 * rateOfChangeSemiMajorAxis );
    keplerianState( eccentricityIndex ) = eccentricity + centuriesAfterMJD2000 *
            rateOfChangeEccentricity;
    keplerianState( inclinationIndex ) = convertDegreesToRadians(
            inclination + centuriesAfterMJD2000 * rateOfChangeInclination );
    keplerianState( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians(
                longitudeOfAscendingNode + centuriesAfterMJD2000 *
                rateOfChangeLongitudeOfAscendingNode );

    // Calculate the longitude of perihelion and mean longitude at the given date. Also converted
    // to SI units.
    const double longitudeOfPerihelionAtGivenDate = convertDegreesToRadians(
                longitudeOfPerihelion + centuriesAfterMJD2000 * rateOfChangeLongitudeOfPerihelion );
    const double meanLongitudeAtGivenDate = convertDegreesToRadians(
                meanLongitude + centuriesAfterMJD2000 * rateOfChangeMeanLongitude );

    // Calculate the argument of periapsis.
    keplerianState( argumentOfPeriapsisIndex ) = longitudeOfPerihelionAtGivenDate -
            keplerianState( longitudeOfAscendingNodeIndex );

    // Calculate mean anomaly at the given date. Needs to be converted to the 0 - 2*pi spectrum.
    const double meanAnomaly = tudat::basic_mathematics::computeModulo(
                meanLongitudeAtGivenDate - longitudeOfPerihelionAtGivenDate,
                2.0 * tudat::mathematical_constants::PI );

    // Calculate the eccentric anomaly.                                                                );
    const double eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly(
                keplerianState( eccentricityIndex ), meanAnomaly );

    // Calculate the true anomaly.
    keplerianState( trueAnomalyIndex ) = convertEccentricAnomalyToTrueAnomaly(
                eccentricAnomaly, keplerianState( eccentricityIndex ) );

    // Return the keplerian elements.
    return keplerianState;
}

//! Get Keplerian Elements.
Eigen::VectorXd EphemerisApproximatePlanetPositionsOuter::
        getKeplerianElements( const double timeInMJD2000 )
{
    // Declare some namespaces that we want to use.
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::unit_conversions;
    using tudat::mathematical_constants::PI;

    // Initialize the keplerian state vector.
    tudat::basic_mathematics::Vector6d keplerianState;

    // Calculate the number of centuries after MJD 2000.
    const double centuriesAfterMJD2000 = timeInMJD2000 / 36525.0;

    // Calculate the semi major axis, eccentricity, inclination and longitude of ascending node.
    // Note that some variables have to be converted to SI units.
    keplerianState( semiMajorAxisIndex ) = convertAstronomicalUnitsToMeters(
            semiMajorAxis + centuriesAfterMJD2000 * rateOfChangeSemiMajorAxis );
    keplerianState( eccentricityIndex ) = eccentricity + centuriesAfterMJD2000 *
            rateOfChangeEccentricity;
    keplerianState( inclinationIndex ) = convertDegreesToRadians(
            inclination + centuriesAfterMJD2000 * rateOfChangeInclination );
    keplerianState( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians(
                longitudeOfAscendingNode + centuriesAfterMJD2000 *
                rateOfChangeLongitudeOfAscendingNode );

    // Calculate the longitude of perihelion and mean longitude at the given date. Also converted
    // to SI units.
    const double longitudeOfPerihelionAtGivenDate = convertDegreesToRadians(
                longitudeOfPerihelion + centuriesAfterMJD2000 * rateOfChangeLongitudeOfPerihelion );
    const double meanLongitudeAtGivenDate = convertDegreesToRadians(
                meanLongitude + centuriesAfterMJD2000 * rateOfChangeMeanLongitude );

    // Calculate the argument of periapsis.
    keplerianState( argumentOfPeriapsisIndex ) = longitudeOfPerihelionAtGivenDate -
            keplerianState( longitudeOfAscendingNodeIndex );

    // Calculate mean anomaly at the given date. Needs to be converted to the 0 - 2*pi spectrum.
    const double meanAnomaly = tudat::basic_mathematics::computeModulo(
                meanLongitudeAtGivenDate - longitudeOfPerihelionAtGivenDate +
                PI / 180.0 * bParameter * centuriesAfterMJD2000 * centuriesAfterMJD2000 +
                PI / 180.0 * cParameter
                * std::cos( fParameter * centuriesAfterMJD2000 * PI / 180.0 ) +
                PI / 180.0 * sParameter *
                std::sin( fParameter * centuriesAfterMJD2000 * PI / 180.0 ),
                2.0 * PI );

    // Calculate the eccentric anomaly.
    const double eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly(
                keplerianState( eccentricityIndex ), meanAnomaly );

    // Calculate the true anomaly.
    keplerianState( trueAnomalyIndex ) = convertEccentricAnomalyToTrueAnomaly(
                eccentricAnomaly, keplerianState( eccentricityIndex ) );

    // Return the keplerian elements.
    return keplerianState;
}

//! Get Keplerian Elements.
Eigen::VectorXd EphemerisApproximatePlanetPositionsPluto::
        getKeplerianElements( const double timeInMJD2000 )
{
    // Declare some namespaces that we want to use.
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::unit_conversions;
    using tudat::mathematical_constants::PI;

    // Initialize the keplerian state vector.
    tudat::basic_mathematics::Vector6d keplerianState;

    // Calculate the number of centuries after MJD 2000.
    const double centuriesAfterMJD2000 = timeInMJD2000 / 36525.0;

    // Calculate the semi major axis, eccentricity, inclination and longitude of ascending node.
    // Note that some variables have to be converted to SI units.
    keplerianState( semiMajorAxisIndex ) = convertAstronomicalUnitsToMeters(
            semiMajorAxis + centuriesAfterMJD2000 * rateOfChangeSemiMajorAxis );
    keplerianState( eccentricityIndex ) = eccentricity + centuriesAfterMJD2000 *
            rateOfChangeEccentricity;
    keplerianState( inclinationIndex ) = convertDegreesToRadians(
            inclination + centuriesAfterMJD2000 * rateOfChangeInclination );
    keplerianState( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians(
                longitudeOfAscendingNode + centuriesAfterMJD2000 *
                rateOfChangeLongitudeOfAscendingNode );

    // Calculate the longitude of perihelion and mean longitude at the given date. Also converted
    // to SI units.
    const double longitudeOfPerihelionAtGivenDate = convertDegreesToRadians(
                longitudeOfPerihelion + centuriesAfterMJD2000 * rateOfChangeLongitudeOfPerihelion );
    const double meanLongitudeAtGivenDate = convertDegreesToRadians(
                meanLongitude + centuriesAfterMJD2000 * rateOfChangeMeanLongitude );

    // Calculate the argument of periapsis.
    keplerianState( argumentOfPeriapsisIndex ) = longitudeOfPerihelionAtGivenDate -
            keplerianState( longitudeOfAscendingNodeIndex );

    // Calculate mean anomaly at the given date. Needs to be converted to the 0 - 2*pi spectrum.
    const double meanAnomaly = tudat::basic_mathematics::computeModulo(
                meanLongitudeAtGivenDate - longitudeOfPerihelionAtGivenDate +
                PI / 180.0 * bParameter * centuriesAfterMJD2000 * centuriesAfterMJD2000,
                2.0 * PI );

    // Calculate the eccentric anomaly.
    const double eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly(
                keplerianState( eccentricityIndex ), meanAnomaly );

    // Calculate the true anomaly.
    keplerianState( trueAnomalyIndex ) = convertEccentricAnomalyToTrueAnomaly(
                eccentricAnomaly, keplerianState( eccentricityIndex ) );

    // Return the keplerian elements.
    return keplerianState;
}

} // namespace ephemeris
} // namespace sample_return_mission
} // namespace tudat_course
