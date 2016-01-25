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
 *      120914    P. Musegaas       File created.
 *      130225    K. Kumar          Removed tudat namespace and updated include guard name.
 *
 *    References
 *      Wertz, J., 2001, Orbit & Constellation Design & Management, Microcosm Press, New York, USA.
 *
 *    Notes
 *      This file contains the mission definitions for the mission.
 *
 */

#ifndef SAMPLE_RETURN_MISSION_MISSION_DEFINITIONS_H
#define SAMPLE_RETURN_MISSION_MISSION_DEFINITIONS_H

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include "Ephemeris/getEphemeris.h"

namespace tudat_course
{
namespace sample_return_mission
{
namespace mission_definitions
{

// Declare constants for calculations. Values taken from (Wertz, 2001) in SI units.
// Note that 'const' prevents the parameters from being changeable. Always use it if a
// parameter should not be changed: it is much safer and equally fast!
const double gravitationalParameterSun = 1.327178e20;
const double gravitationalParameterEarth = 3.98600441e14;
const double gravitationalParameterTarget = 4.2832e13;
const double radiusEarth = 6378136.0;

// Declare mission parameters.
const double parkingOrbitSemiMajorAxis = radiusEarth + 2.0e5;
const double parkingOrbitEccentricity = 0.0;
const double captureOrbitSemiMajorAxis = 4.5e6;
const double captureOrbitEccentricity = 0.2;
const double earthReturnOrbitSemiMajorAxis = radiusEarth + 4.0e5;
const double earthReturnOrbitEccentricity = 0.0;

// Set parameters for ephemeris. Changing between 2D and 3D ephemeris is done here.
const ephemeris::PlanetIndices ephemerisEarth = ephemeris::Earth2D;
const ephemeris::PlanetIndices ephemerisTarget = ephemeris::Mars2D;

// Declare additional constants for part (a).
const double orbitalRadiusEarth = tudat::unit_conversions::
        convertAstronomicalUnitsToMeters( 1.00000011 );
const double orbitalRadiusTarget = tudat::unit_conversions::
        convertAstronomicalUnitsToMeters( 1.52366231 );

} // namespace mission_definitions
} // namespace sample_return_mission
} // namespace tudat_course

#endif // SAMPLE_RETURN_MISSION_MISSION_DEFINITIONS_H
