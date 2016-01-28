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
 *      120824    P. Musegaas       Creation of the code.
 *      120914    P. Musegaas       Update after discussion on course.
 *      130225    K. Kumar          Removed tudat namespace.
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the Major Planets,
 *          http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf, last accessed: 24 February, 2011.
 *
 *    Notes
 *      This file is to be replaced once there is a better version of ephemeris available in Tudat.
 *
 */

#include <stdexcept>

#include "spaceProjectAssignment1/Ephemeris/approximatePlanets.h"
#include "spaceProjectAssignment1/Ephemeris/approximatePlanetsEcliptic.h"
#include "spaceProjectAssignment1/Ephemeris/getEphemeris.h"

namespace tudat_course
{
namespace sample_return_mission
{
namespace ephemeris
{

//! Get the Keplerian elements of the specified planet at the specified time in MJD2000.
Eigen::VectorXd getKeplerianElementsOfPlanet( int Planet, double timeInMJD2000)
{
    Eigen::VectorXd keplerianElements( 6 );
    switch( Planet )
    {
        case Mercury3D:
        {
            APP_Mercury ephemerisObject;
            keplerianElements = ephemerisObject.getKeplerianElements( timeInMJD2000 );
            break;
        }
        case Venus3D:
        {
            APP_Venus ephemerisObject;
            keplerianElements = ephemerisObject.getKeplerianElements( timeInMJD2000 );
            break;
        }
        case Earth3D:
        {
            APP_Earth ephemerisObject;
            keplerianElements = ephemerisObject.getKeplerianElements( timeInMJD2000 );
            break;
        }
        case Mars3D:
        {
            APP_Mars ephemerisObject;
            keplerianElements = ephemerisObject.getKeplerianElements( timeInMJD2000 );
            break;
        }
        case Jupiter3D:
        {
            APP_Jupiter ephemerisObject;
            keplerianElements = ephemerisObject.getKeplerianElements( timeInMJD2000 );
            break;
        }
        case Saturn3D:
        {
            APP_Saturn ephemerisObject;
            keplerianElements = ephemerisObject.getKeplerianElements( timeInMJD2000 );
            break;
        }
        case Uranus3D:
        {
            APP_Uranus ephemerisObject;
            keplerianElements = ephemerisObject.getKeplerianElements( timeInMJD2000 );
            break;
        }
        case Neptune3D:
        {
            APP_Neptune ephemerisObject;
            keplerianElements = ephemerisObject.getKeplerianElements( timeInMJD2000 );
            break;
        }
        case Pluto3D:
        {
            APP_Pluto ephemerisObject;
            keplerianElements = ephemerisObject.getKeplerianElements( timeInMJD2000 );
            break;
        }
        case Mercury2D:
        {
            APP_Mercury_Ecliptic ephemerisObject;
            keplerianElements = ephemerisObject.getKeplerianElements( timeInMJD2000 );
            break;
        }
        case Venus2D:
        {
            APP_Venus_Ecliptic ephemerisObject;
            keplerianElements = ephemerisObject.getKeplerianElements( timeInMJD2000 );
            break;
        }
        case Earth2D:
        {
            APP_Earth_Ecliptic ephemerisObject;
            keplerianElements = ephemerisObject.getKeplerianElements( timeInMJD2000 );
            break;
        }
        case Mars2D:
        {
            APP_Mars_Ecliptic ephemerisObject;
            keplerianElements = ephemerisObject.getKeplerianElements( timeInMJD2000 );
            break;
        }
        case Jupiter2D:
        {
            APP_Jupiter_Ecliptic ephemerisObject;
            keplerianElements = ephemerisObject.getKeplerianElements( timeInMJD2000 );
            break;
        }
        case Saturn2D:
        {
            APP_Saturn_Ecliptic ephemerisObject;
            keplerianElements = ephemerisObject.getKeplerianElements( timeInMJD2000 );
            break;
        }
        case Uranus2D:
        {
            APP_Uranus_Ecliptic ephemerisObject;
            keplerianElements = ephemerisObject.getKeplerianElements( timeInMJD2000 );
            break;
        }
        case Neptune2D:
        {
            APP_Neptune_Ecliptic ephemerisObject;
            keplerianElements = ephemerisObject.getKeplerianElements( timeInMJD2000 );
            break;
        }
        case Pluto2D:
        {
            APP_Pluto_Ecliptic ephemerisObject;
            keplerianElements = ephemerisObject.getKeplerianElements( timeInMJD2000 );
            break;
        }
        default:
        {
            std::runtime_error( "An invalid argument for the planet was given." );
        }
    }
    return keplerianElements;
}
} // namespace ephemeris
} // namespace sample_return_mission
} // namespace tudat_course
