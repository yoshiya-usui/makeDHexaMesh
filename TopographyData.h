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
#ifndef DBLDEF_TOPOGRAPHY_DATA
#define DBLDEF_TOPOGRAPHY_DATA

#include "TopographyData.h"
#include "CommonParameters.h"
#include <vector>
#include <string>

// Class of topography data
class TopographyData{

public:
	// Default constructer
	TopographyData();

	// Constructer
	TopographyData( const std::string& inFileName, const int maxNumPoint, const double maxDist, const double eps );

	// Destructer
	~TopographyData();

	// Read topography data
	void readTopographyData( const double xMin, const double xMax, const double yMin, const double yMax );

	// Interpolate Z coordinate
	double interpolateZCoord( const CommonParameters::XY& coord ) const;

	// Calculated average Z coordinate in a rectangular area
	double calcAverageZCoord( const double xMin, const double xMax, const double yMin, const double yMax );

private:
	// Copy constructer
	TopographyData(const TopographyData& rhs);

	// Assignment operator
	TopographyData& operator=(TopographyData& rhs);

	// Coordinates
	std::vector<CommonParameters::XYZ> m_coords; 

	// File name of topographyt data
	std::string m_fileNameOfTopographytData;

	// Maximum number of points used for interpolating altitudes
	int m_maxNumOfPointsForInterpolating;

	// Maximum distance used for interpolating altitudes
	double m_maxDistanceForInterpolating;

	// Distance used to avoid too small denominator in inverse distance weighting
	double m_distanceUsedToAvoidTooSmallDenominator;

	void outputVTK( const std::string& fname ) const;

#ifdef _TOPO_FUNC
	// Calculate z coordinate by a function
	double calcZByFunction( const CommonParameters::XY& coord ) const;
#endif

};

#endif
