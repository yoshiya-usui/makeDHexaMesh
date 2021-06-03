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
#include "math.h"
#include "Util.h"

// Calculate length between inputted two coordinates
double Util::calcLength( const CommonParameters::XY& coord0, const CommonParameters::XY& coord1 ){

	return hypot( coord0.X - coord1.X,  coord0.Y - coord1.Y );

}

// Calculate length between inputted two coordinates
double Util::calcLength( const CommonParameters::XYZ& coord0, const CommonParameters::XYZ& coord1 ){
	
	double length(0.0);
	length += pow( coord0.X - coord1.X, 2 );
	length += pow( coord0.Y - coord1.Y, 2 );
	length += pow( coord0.Z - coord1.Z, 2 );

	return sqrt( length );

}

// Average two coordinates
CommonParameters::XYZ Util::averageTwoCoords( const CommonParameters::XYZ& coord1, const CommonParameters::XYZ& coord2 ){

	const CommonParameters::XYZ coordAvg = {
		0.5 * ( coord1.X + coord2.X ),
		0.5 * ( coord1.Y + coord2.Y ),
		0.5 * ( coord1.Z + coord2.Z )
	};

	return coordAvg;

}

// Check whether the location of the two points are same
bool Util::isSameLocation( const CommonParameters::XYZ& coord1, const CommonParameters::XYZ& coord2 ){

	const double eps = 1.0e-6;
	if( fabs(coord1.X - coord2.X) < eps && fabs(coord1.Y - coord2.Y) < eps && fabs(coord1.Z - coord2.Z) < eps ){
		return true;
	}
	return false;

}
