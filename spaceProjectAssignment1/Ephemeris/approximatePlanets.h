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

#ifndef SAMPLE_RETURN_MISSION_APPROXIMATE_PLANETS_H
#define SAMPLE_RETURN_MISSION_APPROXIMATE_PLANETS_H

#include <Eigen/Core>

#include "spaceProjectAssignment1/Ephemeris/ephemerisBase.h"

namespace tudat_course
{
namespace sample_return_mission
{
namespace ephemeris
{

//! Ephemeris object for approximate planet positions for the inner planets.
/*!
 * Object containing ephemeris data for approximate planet positions for the inner planets.
 */
class EphemerisApproximatePlanetPositionsInner : public EphemerisObject
{
public:
    //! Constructor.
    /*!
     * Constructor that sets all the constants.
     */
    EphemerisApproximatePlanetPositionsInner( const double semiMajorAxis_,
                                              const double eccentricity_,
                                              const double inclination_,
                                              const double meanLongitude_,
                                              const double longitudeOfPerihelion_,
                                              const double longitudeOfAscendingNode_,
                                              const double rateOfChangeSemiMajorAxis_,
                                              const double rateOfChangeEccentricity_,
                                              const double rateOfChangeInclination_,
                                              const double rateOfChangeMeanLongitude_,
                                              const double rateOfChangeLongitudeOfPerihelion_,
                                              const double rateOfChangeLongitudeOfAscendingNode_ )
        : semiMajorAxis( semiMajorAxis_ ),
          eccentricity( eccentricity_ ),
          inclination( inclination_ ),
          meanLongitude( meanLongitude_ ),
          longitudeOfPerihelion( longitudeOfPerihelion_ ),
          longitudeOfAscendingNode( longitudeOfAscendingNode_ ),
          rateOfChangeSemiMajorAxis( rateOfChangeSemiMajorAxis_ ),
          rateOfChangeEccentricity( rateOfChangeEccentricity_ ),
          rateOfChangeInclination( rateOfChangeInclination_ ),
          rateOfChangeMeanLongitude( rateOfChangeMeanLongitude_ ),
          rateOfChangeLongitudeOfPerihelion( rateOfChangeLongitudeOfPerihelion_ ),
          rateOfChangeLongitudeOfAscendingNode( rateOfChangeLongitudeOfAscendingNode_ )
    { }

    //! Get Keplerian Elements.
    /*!
     * Obtain the Keplerian Elements of the object at a given date.
     */
    Eigen::VectorXd getKeplerianElements( const double timeInMJD2000 );

protected:

    // Variables from the approximate planet positions database for inner planets.
    const double semiMajorAxis;
    const double eccentricity;
    const double inclination;
    const double meanLongitude;
    const double longitudeOfPerihelion;
    const double longitudeOfAscendingNode;
    const double rateOfChangeSemiMajorAxis;
    const double rateOfChangeEccentricity;
    const double rateOfChangeInclination;
    const double rateOfChangeMeanLongitude;
    const double rateOfChangeLongitudeOfPerihelion;
    const double rateOfChangeLongitudeOfAscendingNode;

private:
};

//! Ephemeris object for approximate planet positions for the outer planets.
/*!
 * Object containing ephemeris data for approximate planet positions for the outer planets.
 */
class EphemerisApproximatePlanetPositionsOuter : public EphemerisObject
{
public:
    //! Constructor.
    /*!
     * Constructor that sets all the constants.
     */
    EphemerisApproximatePlanetPositionsOuter( const double semiMajorAxis_,
                                              const double eccentricity_,
                                              const double inclination_,
                                              const double meanLongitude_,
                                              const double longitudeOfPerihelion_,
                                              const double longitudeOfAscendingNode_,
                                              const double rateOfChangeSemiMajorAxis_,
                                              const double rateOfChangeEccentricity_,
                                              const double rateOfChangeInclination_,
                                              const double rateOfChangeMeanLongitude_,
                                              const double rateOfChangeLongitudeOfPerihelion_,
                                              const double rateOfChangeLongitudeOfAscendingNode_,
                                              const double bParameter_,
                                              const double cParameter_,
                                              const double sParameter_,
                                              const double fParameter_ )
        : semiMajorAxis( semiMajorAxis_ ),
          eccentricity( eccentricity_ ),
          inclination( inclination_ ),
          meanLongitude( meanLongitude_ ),
          longitudeOfPerihelion( longitudeOfPerihelion_ ),
          longitudeOfAscendingNode( longitudeOfAscendingNode_ ),
          rateOfChangeSemiMajorAxis( rateOfChangeSemiMajorAxis_ ),
          rateOfChangeEccentricity( rateOfChangeEccentricity_ ),
          rateOfChangeInclination( rateOfChangeInclination_ ),
          rateOfChangeMeanLongitude( rateOfChangeMeanLongitude_ ),
          rateOfChangeLongitudeOfPerihelion( rateOfChangeLongitudeOfPerihelion_ ),
          rateOfChangeLongitudeOfAscendingNode( rateOfChangeLongitudeOfAscendingNode_ ),
          bParameter( bParameter_ ),
          cParameter( cParameter_ ),
          sParameter( sParameter_ ),
          fParameter( fParameter_ )
    { }

    //! Get Keplerian Elements.
    /*!
     * Obtain the Keplerian Elements of the object at a given date.
     */
    Eigen::VectorXd getKeplerianElements( const double timeInMJD2000 );

protected:

    // Variables from the approximate planet positions database for inner planets.
    const double semiMajorAxis;
    const double eccentricity;
    const double inclination;
    const double meanLongitude;
    const double longitudeOfPerihelion;
    const double longitudeOfAscendingNode;
    const double rateOfChangeSemiMajorAxis;
    const double rateOfChangeEccentricity;
    const double rateOfChangeInclination;
    const double rateOfChangeMeanLongitude;
    const double rateOfChangeLongitudeOfPerihelion;
    const double rateOfChangeLongitudeOfAscendingNode;
    const double bParameter;
    const double cParameter;
    const double sParameter;
    const double fParameter;

private:
};

//! Ephemeris object for approximate planet positions for Pluto.
/*!
 * Object containing ephemeris data for approximate planet positions for Pluto.
 */
class EphemerisApproximatePlanetPositionsPluto : public EphemerisObject
{
public:
    //! Constructor.
    /*!
     * Constructor that sets all the constants.
     */
    EphemerisApproximatePlanetPositionsPluto( const double semiMajorAxis_,
                                              const double eccentricity_,
                                              const double inclination_,
                                              const double meanLongitude_,
                                              const double longitudeOfPerihelion_,
                                              const double longitudeOfAscendingNode_,
                                              const double rateOfChangeSemiMajorAxis_,
                                              const double rateOfChangeEccentricity_,
                                              const double rateOfChangeInclination_,
                                              const double rateOfChangeMeanLongitude_,
                                              const double rateOfChangeLongitudeOfPerihelion_,
                                              const double rateOfChangeLongitudeOfAscendingNode_,
                                              const double bParameter_ )
        : semiMajorAxis( semiMajorAxis_ ),
          eccentricity( eccentricity_ ),
          inclination( inclination_ ),
          meanLongitude( meanLongitude_ ),
          longitudeOfPerihelion( longitudeOfPerihelion_ ),
          longitudeOfAscendingNode( longitudeOfAscendingNode_ ),
          rateOfChangeSemiMajorAxis( rateOfChangeSemiMajorAxis_ ),
          rateOfChangeEccentricity( rateOfChangeEccentricity_ ),
          rateOfChangeInclination( rateOfChangeInclination_ ),
          rateOfChangeMeanLongitude( rateOfChangeMeanLongitude_ ),
          rateOfChangeLongitudeOfPerihelion( rateOfChangeLongitudeOfPerihelion_ ),
          rateOfChangeLongitudeOfAscendingNode( rateOfChangeLongitudeOfAscendingNode_ ),
          bParameter( bParameter_ )
    { }

    //! Get Keplerian Elements.
    /*!
     * Obtain the Keplerian Elements of the object at a given date.
     */
    Eigen::VectorXd getKeplerianElements( const double timeInMJD2000 );

protected:

    // Variables from the approximate planet positions database for inner planets.
    const double semiMajorAxis;
    const double eccentricity;
    const double inclination;
    const double meanLongitude;
    const double longitudeOfPerihelion;
    const double longitudeOfAscendingNode;
    const double rateOfChangeSemiMajorAxis;
    const double rateOfChangeEccentricity;
    const double rateOfChangeInclination;
    const double rateOfChangeMeanLongitude;
    const double rateOfChangeLongitudeOfPerihelion;
    const double rateOfChangeLongitudeOfAscendingNode;
    const double bParameter;

private:
};

// From here on initialize the implementations of the different ephemerides.

class APP_Mercury : public EphemerisApproximatePlanetPositionsInner
{
public:

    APP_Mercury( )
        : EphemerisApproximatePlanetPositionsInner(
              0.38709843, 0.20563661, 7.00559432, 252.25166724, 77.45771895, 48.33961819,
              0.00000000, 0.00002123, -0.00590158, 149472.67486623, 0.15940013, -0.12214182 )
    { }

protected:

private:
};

class APP_Venus : public EphemerisApproximatePlanetPositionsInner
{
public:

    APP_Venus( )
        : EphemerisApproximatePlanetPositionsInner(
              0.72332102, 0.00676399, 3.39777545, 181.97970850, 131.76755713, 76.67261496,
              -0.00000026, -0.00005107, 0.00043494, 58517.81560260, 0.05679648, -0.27274174 )
    { }
protected:

private:

};

class APP_Earth : public EphemerisApproximatePlanetPositionsInner
{
public:

    APP_Earth( )
        : EphemerisApproximatePlanetPositionsInner(
              1.00000018, 0.01673163, -0.00054346, 100.46691572, 102.93005885, -5.11260389,
              -0.00000003, -0.00003661, -0.01337178, 35999.37306329, 0.31795260, -0.24123856 )
    { }

protected:

private:
};

class APP_Mars : public EphemerisApproximatePlanetPositionsInner
{
public:

    APP_Mars( )
        : EphemerisApproximatePlanetPositionsInner(
              1.52371243, 0.09336511, 1.85181869, -4.56813164, -23.91744784, 49.71320984,
              0.00000097, 0.00009149, -0.00724757, 19140.29934243, 0.45223625, -0.26852431 )
    { }

protected:

private:
};

class APP_Jupiter : public EphemerisApproximatePlanetPositionsOuter
{
public:

    APP_Jupiter( )
        : EphemerisApproximatePlanetPositionsOuter(
              5.20248019, 0.04853590, 1.29861416, 34.33479152, 14.27495244, 100.29282654,
              -0.00002864, 0.00018026, -0.00322699, 3034.90371757, 0.18199196, 0.13024619,
              -0.00012452, 0.06064060, -0.35635438, 38.35125000 )
    { }

protected:

private:
};

class APP_Saturn : public EphemerisApproximatePlanetPositionsOuter
{
public:

    APP_Saturn( )
        : EphemerisApproximatePlanetPositionsOuter(
              9.54149883, 0.05550825, 2.49424102, 50.07571329, 92.86136063, 113.63998702,
              -0.00003065, -0.00032044, 0.00451969, 1222.11494724, 0.54179478, -0.25015002,
              0.00025899, -0.13434469, 0.87320147, 38.35125000 )
    { }

protected:

private:
};

class APP_Uranus : public EphemerisApproximatePlanetPositionsOuter
{
public:

    APP_Uranus( )
        : EphemerisApproximatePlanetPositionsOuter(
              19.18797948, 0.04685740, 0.77298127, 314.20276625, 172.43404441, 73.96250215,
              -0.00020455, -0.00001550, -0.00180155, 428.49512595, 0.09266985, 0.05739699,
              0.00058331, -0.97731848, 0.17689245, 7.67025000 )
    { }

protected:

private:
};

class APP_Neptune : public EphemerisApproximatePlanetPositionsOuter
{
public:

    APP_Neptune( )
        : EphemerisApproximatePlanetPositionsOuter(
              30.06952752, 0.00895439, 1.77005520, 304.22289287, 46.68158724, 131.78635853,
              0.00006447, 0.00000818, 0.00022400, 218.46515314, 0.01009938, -0.00606302,
              -0.00041348, 0.68346318, -0.10162547, 7.67025000 )
    { }

protected:

private:
};

class APP_Pluto : public EphemerisApproximatePlanetPositionsPluto
{
public:

    APP_Pluto( )
        : EphemerisApproximatePlanetPositionsPluto(
              39.48686035, 0.24885238, 17.14104260, 238.96535011, 224.09702598, 110.30167986,
              0.00449751, 0.00006016, 0.00000501, 145.18042903, -0.00968827, -0.00809981,
              -0.01262724 )
    { }

protected:

private:
};

} // namespace ephemeris
} // namespace sample_return_mission
} // namespace tudat_course

#endif // SAMPLE_RETURN_MISSION_APPROXIMATE_PLANETS_H
