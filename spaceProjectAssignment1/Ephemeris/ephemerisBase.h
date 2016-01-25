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

#ifndef SAMPLE_RETURN_MISSION_EPHEMERIS_BASE_H
#define SAMPLE_RETURN_MISSION_EPHEMERIS_BASE_H

#include <Eigen/Core>

namespace tudat_course
{
namespace sample_return_mission
{
namespace ephemeris
{

//! Ephemeris object.
/*!
 * Object containing ephemeris data that can be extracted to yield positions of bodies in the
 * heliocentric ecliptic reference frame.
 */
class EphemerisObject
{
public:

    //! Virtual destructor.
    /*!
     * Virtual destructor.
     */
    virtual ~EphemerisObject( ) { }

    //! Get Keplerian Elements.
    /*!
     * Obtain the Keplerian Elements of the object at a given date.
     */
    virtual Eigen::VectorXd getKeplerianElements( const double timeInMJD2000 ) = 0;

protected:

private:
};

} // namespace ephemeris
} // namespace sample_return_mission
} // namespace tudat_course

#endif // SAMPLE_RETURN_MISSION_EPHEMERIS_BASE_H
