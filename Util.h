//--------------------------------------------------------------------------
// This file is part of makeDHexaMesh.
//
// makeDHexaMesh is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// makeDHexaMesh is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with makeDHexaMesh. If not, see <http://www.gnu.org/licenses/>.
//--------------------------------------------------------------------------
#ifndef DBLDEF_UTIL
#define DBLDEF_UTIL

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

#include "CommonParameters.h"

namespace Util{

// Calculate length between inputted two coordinates
double calcLength( const CommonParameters::XY& coord0, const CommonParameters::XY& coord1 );

// Calculate length between inputted two coordinates
double calcLength( const CommonParameters::XYZ& coord0, const CommonParameters::XYZ& coord1 );

// Average two coordinates
CommonParameters::XYZ averageTwoCoords( const CommonParameters::XYZ& coord1, const CommonParameters::XYZ& coord2 );

// Check whether the location of the two points are same
bool isSameLocation( const CommonParameters::XYZ& coord1, const CommonParameters::XYZ& coord2 );

}

#endif
