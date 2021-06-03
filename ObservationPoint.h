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
#ifndef DBLDEF_OBSERVATION_POINT
#define DBLDEF_OBSERVATION_POINT

#include <fstream>
#include <iostream>
#include "CommonParameters.h"

// Class of observation point
class ObservationPoint{

public:
	// Default Constructer
	ObservationPoint();

	// Destructer
	~ObservationPoint();

	// Read data of observation point from input file
	void readObservationPointData( std::ifstream& ifs );

	// Calculate maximum horizontal length
	double calcMaximumHorizontalLength( const CommonParameters::XYZ& coord ) const;

	// Calculate maximum vertical length
	double calcMaximumVerticalLength( const CommonParameters::XYZ& coord ) const;

	// Calculate factor for intermediate region
	double calcFactorForIntermediateRegion( const CommonParameters::XYZ& coord, const int iSphere ) const;

	//// Judge the distance between the observation station and the specified point is lower than threshold value
	//bool judgeDistanceObsSiteAndPointLowerThanThreshold( const CommonParameters::XYZ& coord, const double threshold ) const;

private:
	enum TypeOfRegion{
		SPHERE = 0,
		CUBOID,
	};

	// Copy constructer
	ObservationPoint(const ObservationPoint& rhs);

	// Assignment operator
	ObservationPoint& operator=(ObservationPoint& rhs);

	// Coordinate of the point
	CommonParameters::XYZ m_pointCoord;

	// Type of region
	int m_typeOfRegion;

	// Total number of spheres
	int m_numRegions;

	// Length for local mesh control
	double* m_length;

	// Maximum horizontal edge length within sphere around the point
	double* m_maxEdgeLengthHorizontal;

	// Maximum vertical edge length within sphere around the point
	double* m_maxEdgeLengthVertical;

	// Oblateness of regions
	double* m_oblateness;

	bool locateInRegion( const CommonParameters::XYZ& coord, const int iRegion ) const;

};

#endif
