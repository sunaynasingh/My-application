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
 *      130225    K. Kumar          Removed tudat namespace and updated include guard name.
 *
 *    References
 *      Standish, E.M. Keplerian Elements for Approximate Positions of the Major Planets,
 *          http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf, last accessed: 24 February, 2011.
 *
 *    Notes
 *      This file is to be replaced once there is a better version of ephemeris available in Tudat.
 *
 */

#ifndef SAMPLE_RETURN_MISSION_GET_EPHEMERIS_H
#define SAMPLE_RETURN_MISSION_GET_EPHEMERIS_H

#include <Eigen/Core>

namespace tudat_course
{
namespace sample_return_mission
{
namespace ephemeris
{

//! Ephemeris enumeration
/*!
 * An enumeration to define the planet of which the ephemeris is desired. To be replaced once a
 * better version of ephemeris is available.
 * 1 = Mercury in 3D,
 * 2 = Venus in 3D,
 * 3 = Earth in 3D,
 * 4 = Mars in 3D,
 * 5 = Jupiter in 3D,
 * 6 = Saturn in 3D,
 * 7 = Uranus in 3D,
 * 8 = Neptune in 3D,
 * 9 = Pluto in 3D,
 * 11 = Mercury in 2D,
 * 12 = Venus in 2D,
 * 13 = Earth in 2D,
 * 14 = Mars in 2D,
 * 15 = Jupiter in 2D,
 * 16 = Saturn in 2D,
 * 17 = Uranus in 2D,
 * 18 = Neptune in 2D,
 * 19 = Pluto in 2D.
 */
enum PlanetIndices
{
    Mercury3D = 1, Venus3D, Earth3D, Mars3D, Jupiter3D, Saturn3D, Uranus3D, Neptune3D, Pluto3D,
    Mercury2D = 11, Venus2D, Earth2D, Mars2D, Jupiter2D, Saturn2D, Uranus2D, Neptune2D, Pluto2D
};

//! Get the Keplerian elements of the specified planet at the specified time in MJD2000.
/*!
 * Calculates the Keplerian elements of Mars at a specified time in MJD2000. Obtained using
 * approximate planet position equations described in (Standish, 2011).
 * \param planet Index of the planet of which the ephemeris is desired.
 * \param timeInMJD2000 Time at which the ephemeris is desired measured in MJD2000.
 * \return Keplerian elements of the planet orbit at the desired date.
 */
Eigen::VectorXd getKeplerianElementsOfPlanet( int Planet, double timeInMJD2000 );

} // namespace ephemeris
} // namespace sample_return_mission
} // namespace tudat_course

#endif // SAMPLE_RETURN_MISSION_GET_EPHEMERIS_H
