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
 *      121218    K. Kumar          Code created.
 *      130225    K. Kumar          Removed tudat namespace and updated include guard name.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef SAMPLE_RETURN_MISSION_BASIC_INPUT_OUTPUT_H
#define SAMPLE_RETURN_MISSION_BASIC_INPUT_OUTPUT_H

#include <string>

namespace tudat_course
{
namespace sample_return_mission
{
namespace basic_input_output
{

//! Get root-path for Tudat Core library.
/*!
 * Returns root-path corresponding with root-directory of Tudat Core library as a string with
 * trailing slash included.
 * \return Tudat Core root-path.
 */
static std::string getApplicationRootPath( )
{
#ifdef SAMPLE_RETURN_MISSION_CUSTOM_ROOT_PATH
    return std::string( SAMPLE_RETURN_MISSION_CUSTOM_ROOT_PATH );
#else
    // Declare file path string assigned to filePath.
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    return filePath_.substr( 0, filePath_.length( ) -
                                std::string( "basicInputOutput.h" ).length( ) );
#endif
}

} // namespace basic_input_output
} // namespace sample_return_mission
} // namespace tudat_course

#endif // SAMPLE_RETURN_MISSION_BASIC_INPUT_OUTPUT_H
