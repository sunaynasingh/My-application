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
 *      121030    K. Kumar          File created.
 *      130225    K. Kumar          Updated include-guard and namespace names; updated Vector6d
 *                                  references to use Tudat definition.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef SATELLITE_PROPAGATOR_EXAMPLES_BODY_H
#define SATELLITE_PROPAGATOR_EXAMPLES_BODY_H

#include <map>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
#include <Tudat/Astrodynamics/Ephemerides/ephemeris.h>
#include <Tudat/Astrodynamics/Gravitation/gravityFieldModel.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>
#include <Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h>

namespace tudat_application
{

//! Test body class.
/*!
 *  This class serves as an example of how a container can be constructed that stores state and
 *  time information, which can be used in conjunction with acceleration models, the Cartesian state
 *  derivative model, and the composite state derivative model. It should be noted that this class
 *  should can be used "as is", but need to be may be modified taking intp consideration
 *  the application at hand. Classes such as this are application-specific, hence unavailable
 *  in the Tudat libraries.
 */
class Body
{
public:

    //! Constructor for a body
    /*!
     * Constructor for a body, sets current time, state, rotation and mass values
     * (all with default parameters). The input state is used internally to
     * set the current position (taken as a segment of the input state given by the indices
     * (0, 3)) and the current velocity (taken as a segment of the input state given by the indices
     * (3, 3).
     * \param state Current state of body at initialization (default = zeroes).
     * \param time Current time of body at initialization (default = zeroes).
     * \param bodyMass Current mass of body at initialization (default = zeroes).
     * \param currentRotationToGlobalFrame Current rotation of from body-fixed to inertial frames
     *  at initialization (default = identity)
     */
    Body( const tudat::basic_mathematics::Vector6d& state =
            tudat::basic_mathematics::Vector6d::Zero( ),
          const double time = 0.0, const double bodyMass = 0.0,
          const Eigen::Quaterniond currentRotationToGlobalFrame =
            Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) )
        : currentState( state ),
          currentPosition( state.segment( 0, 3 ) ),
          currentVelocity( state.segment( 3, 3 ) ),
          currentTime( time ),
          currentRotationToGlobalFrame_( currentRotationToGlobalFrame ),
          bodyMass_( bodyMass )
    { }

    //! Set current time and state.
    /*!
     *  Sets the current time, position and current velocity of the body based on the input
     *  arguments. The current position is taken as a segment of the input state given by the
     *  indices (0, 3)), and the current velocity is taken as a segment of the input state given by
     *  the indices (3, 3).
     *  Note: any updates of dependent variables which depend on time should be made here.
     *  \param time Current time of body (from which to calculate any dependent variables in
     *  future code modifications).
     *  \param state Current state of body.
     */
    void setCurrentTimeAndState( const double time,
                                 const tudat::basic_mathematics::Vector6d& state )
    {
        currentTime = time;
        currentState = state;
        currentPosition = state.segment( 0, 3 );
        currentVelocity = state.segment( 3, 3 );
    }

    //! Update body to current time
    /*!
     *  Update body to current time, calculating the current state from the ephemeris_ member
     *  variable.
     *  \param time Current time of body (from which to calculate the state, as well as any dependent
     *  variables infuture code modifications).
     */
    void updateStateFromEphemeris( const double time )
    {
        setCurrentTimeAndState(
                    time, bodyEphemeris_->getCartesianStateFromEphemeris(
                        time, tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000 ) );
    }

    //! Get current state.
    /*!
     * Returns the internally stored current state vector.
     * \return Current state.
     */
    tudat::basic_mathematics::Vector6d getCurrentState( ) { return currentState; }

    //! Get current position.
    /*!
     * Returns the internally stored current position vector.
     * \return Current position.
     */
    Eigen::Vector3d getCurrentPosition( ) { return currentPosition; }

    //! Get current velocity.
    /*!
     * Returns the internally stored current velocity vector.
     * \return Current velocity.
     */
    Eigen::Vector3d getCurrentVelocity( ) { return currentVelocity; }

    //! Get current time.
    /*!
     * Returns the internally stored current time.
     * \return Current time.
     */
    double getCurrentTime( ) { return currentTime; }

    //! Function to set the ephemeris of the body.
    /*!
     *  Function to set the ephemeris of the body, which is used to represent the (a priori)
     *  state history of the body.
     *  \param New ephemeris of the body.
     */
    void setEphemeris( const boost::shared_ptr< tudat::ephemerides::Ephemeris > bodyEphemeris )
    {
        bodyEphemeris_ = bodyEphemeris;
    }

    //! Function to set the gravity field of the body.
    /*!
     *  Function to set the gravity field of the body; input is also used to (re)set the mass
     *  of the body.
     *  \param New gravity field of the body.
     */
    void setGravityFieldModel(
            const boost::shared_ptr< tudat::gravitation::GravityFieldModel > gravityFieldModel )
    {
        gravityFieldModel_ = gravityFieldModel;
        bodyMass_ = gravityFieldModel_->getGravitationalParameter( );
    }

    //! Function to get the gravity field model of the body.
    /*!
     *  Function to get the gravity field model of the body.
     *  \return Gravity field model of the body.
     */
    boost::shared_ptr< tudat::gravitation::GravityFieldModel > getGravityFieldModel( )
    {
        return gravityFieldModel_;
    }

    //! Function to get the ephemeris of the body.
    /*!
     *  Function to get the ephemeris of the body.
     *  \return Ephemeris of the body.
     */
    boost::shared_ptr< tudat::ephemerides::Ephemeris > getEphemeris( )
    {
        return bodyEphemeris_;
    }

    //! Get current rotation from body-fixed to inertial frame.
    /*!
     *  Get current rotation from body-fixed to inertial frame.
     *  NOTE: This rotation is currently always the identity rotation, in the future, when the
     *  Body is endowed with a RotationalEphemeris, updates of the currentRotationToGlobalFrame_
     *  by the setCurrentTimeAndState function will allow the rotation to take on different values.
     *  \return Current rotation from body-fixed to inertial frame
     */
    Eigen::Quaterniond getCurrentRotationToGlobalFrame( )
    {
        return currentRotationToGlobalFrame_;
    }

    //! Get current rotation from inertial to body-fixed frame.
    /*!
     *  Get current rotation from inertial to body-fixed frame.
     *  NOTE: This rotation is currently always the identity rotation, in the future, when the
     *  Body is endowed with a RotationalEphemeris, updates of the currentRotationToGlobalFrame_
     *  by the setCurrentTimeAndState function will allow the rotation to take on different values.
     *  \return Current rotation from inertial to body-fixed frame
     */
    Eigen::Quaterniond getCurrentRotationToLocalFrame( )
    {
        return currentRotationToGlobalFrame_.inverse( );
    }

protected:

private:

    //! Current state.
    tudat::basic_mathematics::Vector6d currentState;

    //! Current position.
    Eigen::Vector3d currentPosition;

    //! Current position.
    Eigen::Vector3d currentVelocity;

    //! Current time.
    double currentTime;

    //! Current rotation from body-fixed to inertial frame.
    Eigen::Quaterniond currentRotationToGlobalFrame_;

    //! Mass of body (default set to zero, calculated from GravityFieldModel when it is set).
    double bodyMass_;

    //! Ephemeris of body.
    boost::shared_ptr< tudat::ephemerides::Ephemeris > bodyEphemeris_;

    //! Gravity field model of body.
    boost::shared_ptr< tudat::gravitation::GravityFieldModel > gravityFieldModel_;    

};

} // namespace tudat_application

#endif // SATELLITE_PROPAGATOR_EXAMPLES_BODY_H
