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
