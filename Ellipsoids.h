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
#ifndef DBLDEF_ELLIPSOIDS
#define DBLDEF_ELLIPSOIDS

#include "CommonParameters.h"
#include <fstream>
#include <iostream>
#include <string>

// Class of the ellipsoids deciding maximum edge length
class Ellipsoids{

public:
	// Constructer
	Ellipsoids();

	// Destructer
	~Ellipsoids();

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
	Ellipsoids(const Ellipsoids& rhs);

	// Assignment operator
	Ellipsoids& operator=(const Ellipsoids& rhs);

	// Center coord
	CommonParameters::XYZ m_centerCoord;

	// Rotation angle of elliptical fine zone in X-Y plane
	double m_rotationAngle;

	// Total number of ellipsoids
	int m_numEllipsoids;

	// Length of long axis
	double* m_radius;

	// Maximum horizontal edge length within ellipsoids
	double* m_edgeLengthHorizontal;

	// Maximum vertical edge length within ellipsoids
	double* m_edgeLengthVertical;

	// Oblateness of horizontal direction
	double* m_oblatenessHorizontal;

	// Oblateness of vertical direction
	double* m_oblateness[2];

	// Get total number of ellipsoids
	int getNumOfEllipsoids() const;

	// Check whether specified point is located in the Ellipsoid
	bool locateInEllipsoid( const CommonParameters::XYZ& coord, const int iEllipsoid ) const;

	// Calculate maximum horizontal length between two spheres
	double calcEdgeLengthTransitionHorizontal( const CommonParameters::XYZ& coord, const int iEllipsoid ) const;

	// Calculate maximum vertical length between two spheres
	double calcEdgeLengthTransitionVertical( const CommonParameters::XYZ& coord, const int iEllipsoid ) const;

	// Calculate factor for the maximum edge length between two spheres
	double calcFactorForEdgeLengthTransitionRegion( const CommonParameters::XYZ& coord, const int iEllipsoid ) const;

	// Calculate length on ellipsoid
	double calculateLengthOnEllipsoid( const double angleH, const double angleV, const int iEllipsoid, const int iType ) const;

};

#endif
