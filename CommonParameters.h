﻿//--------------------------------------------------------------------------
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
#ifndef DBLDEF_COMMON_PARAMETERS
#define DBLDEF_COMMON_PARAMETERS

namespace CommonParameters{

struct XY{
	double X;
	double Y;
};

struct XYZ{
	double X;
	double Y;
	double Z;
};

// Circular constant
const static double PI = 3.14159265359;

// Factor converting values from radians to degrees
const static double RAD2DEG = 180.0 / PI;

// Factor converting values from degrees to radians
const static double DEG2RAD = PI / 180.0;

static char programName[]="makeDHexaMesh";

static char version[]="v1.2";

}

#endif
