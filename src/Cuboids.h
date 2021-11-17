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
#ifndef DBLDEF_CUBOIDS
#define DBLDEF_CUBOIDS

#include "CommonParameters.h"
#include <fstream>
#include <iostream>
#include <string>

// Class of the cuboids deciding maximum edge length
class Cuboids{

public:
	// Constructer
	Cuboids();

	// Destructer
	~Cuboids();

	// Read parameters
	void readParameters( const std::string& fileName  );

	// Calculate maximum horizontal length of the specified point
	double calcMaximumEdgeLengthHorizontal( const CommonParameters::XYZ& coord ) const;

	// Calculate maximum horizontal length of the specified point
	double calcMaximumEdgeLengthVertical( const CommonParameters::XYZ& coord ) const;

private:
	enum TypeID{
		EARTH = 0,
		AIR = 1
	};

	// Copy constructer
	Cuboids(const Cuboids& rhs);

	// Assignment operator
	Cuboids& operator=(const Cuboids& rhs);

	// Center coord
	CommonParameters::XYZ m_centerCoord;

	// Rotation angle of elliptical fine zone in X-Y plane
	double m_rotationAngle;

	// Total number of ellipsoids
	int m_numCuboids;

	// Length along x-axis
	double* m_lengthX;

	// Length along y-axis
	double* m_lengthY;

	// Length along z-axis
	double* m_lengthZ;

	// Maximum horizontal edge length within ellipsoids
	double* m_edgeLengthHorizontal;

	// Maximum vertical edge length within ellipsoids
	double* m_edgeLengthVertical;

	// Get total number of cuboids
	int getNumOfCuboids() const;

	// Check whether specified point is located in the cuboid
	bool locateInCuboids( const CommonParameters::XYZ& coord, const int iCuboid ) const;

};

#endif
