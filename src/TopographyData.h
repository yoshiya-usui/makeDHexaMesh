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
#ifndef DBLDEF_TOPOGRAPHY_DATA
#define DBLDEF_TOPOGRAPHY_DATA

#include "TopographyData.h"
#include "CommonParameters.h"
#include "Lake.h"

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

	// Read lake data
	void readLakeData();

	// Interpolate Z coordinate
	double interpolateZCoord( const CommonParameters::XY& coord ) const;

	// Calculated average Z coordinate in a rectangular area
	double calcAverageZCoord( const double xMin, const double xMax, const double yMin, const double yMax );

	// Calculated average Z coordinate in a rectangular area for lake
	double calcAverageZCoordForLake(const double xMin, const double xMax, const double yMin, const double yMax, const int lakeIndex);

	// Read lake data
	void readLakeData(std::ifstream& ifs);

	// Check if the point is located in a lake
	bool isPointLocatedInALake(const double x, const double y, int& lakeIndex) const;

	// Check if lake is included in the model
	bool isLakeIncluded() const;

	// Get lake surface altitude
	double getLakeSurfaceAltitude(const int lakeIndex) const;

	// Get lake resistivity
	double getLakeResistivity(const int lakeIndex) const;

	// Select lake areas
	void selectLakeAreas(const double xMin, const double xMax, const double yMin, const double yMax);

	// Get number of selected lake areas
	int getNumOfSelectLakeAreas() const;

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

	// List of lake objects
	std::vector<Lake> m_lakes;

	void outputVTK( const std::string& fname ) const;

#ifdef _TOPO_FUNC
	// Calculate z coordinate by a function
	double calcZByFunction( const CommonParameters::XY& coord ) const;
#endif

};
#endif
