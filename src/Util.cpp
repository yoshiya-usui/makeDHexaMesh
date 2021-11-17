//--------------------------------------------------------------------------
// MIT License
//
// Copyright (c) 2021 Yoshiya Usui
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
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
