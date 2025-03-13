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
#ifndef DBLDEF_LAKE
#define DBLDEF_LAKE

#include <vector>
#include <fstream>
#include "CommonParameters.h"

// Class of lake
class Lake{

public:

	// Default Constructer
	Lake();

	// Copy constructer
	Lake(const Lake& rhs);

	// Destructer
	~Lake();

	// Read lake data
	void readLakeData(std::ifstream& ifs);

	// Read topography data
	void readLakeBottomDepths();

	// Check if the point is included in an area
	bool isPointIncludedInAnArea(const double x, const double y) const;

	// Get altitude of lake surface
	double getAltitudeOfLakeSurface() const;

	// Get lake resistivity
	double getLakeResistivity() const;

	// Calculated average lake depth in a rectangular area
	double calcAverageLakeDepth(const double xMin, const double xMax, const double yMin, const double yMax, int& numData) const;

	// Calculated average Z coordinate in a rectangular area
	double calcAverageZCoord(const double xMin, const double xMax, const double yMin, const double yMax) const;

	// Select and add lake area
	void selectLakeArea(const double xMin, const double xMax, const double yMin, const double yMax);

	// Get threshould of average lake depth
	double getThreshouldOfAverageLakeDepth() const;

	// Get threshould of number of data point in each element face
	int getThreshouldOfNumOfDataPointInEachElementFace() const;

	// Get number of selected lake areas
	int getNumOfSelectLakeAreas() const;

private:

	struct Area {
		double XStart;
		double XEnd;
		double YStart;
		double YEnd;
	};

	// File name of lake bottom depths
	std::string m_fileNameOfLakeBottomDepths;

	// Altitude of lake surface
	double m_altitudeOfLakeSurface;

	// Resistivity of lake
	double m_lakeResistivity;

	// Threshould of average lake depth
	double m_threshouldOfAverageLakeDepth;

	// Threshould of number of data point in each element face
	int m_threshouldOfNumOfDataPointInEachElementFace;

	// Coordinates
	std::vector<CommonParameters::XYZ> m_coords;

	// List of lake areas
	std::vector<Area> m_listOfLakeAreas;

	// Assignment operator
	Lake& operator=(Lake& rhs);

	// Add area
	void addArea(const double xMin, const double xMax, const double yMin, const double yMax);

};

#endif
