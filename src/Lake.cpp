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
#include "Lake.h"
#include "math.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <assert.h> 

// Default constructer
Lake::Lake():
	m_altitudeOfLakeSurface(0.0),
	m_lakeResistivity(0.0),
	m_threshouldOfAverageLakeDepth(0.1),
	m_threshouldOfNumOfDataPointInEachElementFace(3)
{
}

// Copy constructer
Lake::Lake(const Lake& rhs) {

	m_fileNameOfLakeBottomDepths = rhs.m_fileNameOfLakeBottomDepths;
	m_altitudeOfLakeSurface = rhs.m_altitudeOfLakeSurface;
	m_lakeResistivity = rhs.m_lakeResistivity;
	m_threshouldOfAverageLakeDepth = rhs.m_threshouldOfAverageLakeDepth;
	m_threshouldOfNumOfDataPointInEachElementFace = rhs.m_threshouldOfNumOfDataPointInEachElementFace;
	m_coords = rhs.m_coords;
	m_listOfLakeAreas = rhs.m_listOfLakeAreas;

}

// Destructer
Lake::~Lake(){
}

// Read lake data
void Lake::readLakeData(std::ifstream& ifs) {

	ifs >> m_altitudeOfLakeSurface;
	std::cout << "Altitude of lake surface (km): " << m_altitudeOfLakeSurface << std::endl;

	ifs >> m_lakeResistivity;
	std::cout << "Lake resistivity (Ohm-m): " << m_lakeResistivity << std::endl;

	ifs >> m_threshouldOfAverageLakeDepth;
	std::cout << "Threshould of average lake depth (km): " << m_threshouldOfAverageLakeDepth << std::endl;
	if (m_threshouldOfAverageLakeDepth <= 0.0)
	{
		std::cerr << " Error : threshould of average lake depth should be positive." << std::endl;
		exit(1);
	}
	ifs >> m_threshouldOfNumOfDataPointInEachElementFace;
	std::cout << "Threshould of number of data point in each element face: " << m_threshouldOfNumOfDataPointInEachElementFace << std::endl;
	if (m_threshouldOfNumOfDataPointInEachElementFace <= 0)
	{
		std::cerr << " Error : threshould of number of data point in each element face should be positive." << std::endl;
		exit(1);
	}

	ifs >> m_fileNameOfLakeBottomDepths;

}

// Read topography data
void Lake::readLakeBottomDepths() {

	// Open input file
	std::ifstream ifs(m_fileNameOfLakeBottomDepths.c_str());
	if (!ifs) {
		std::cerr << "Cannot open file " << m_fileNameOfLakeBottomDepths << std::endl;
		exit(1);
	}
	std::cout << "Read lake bottom depths from " << m_fileNameOfLakeBottomDepths << std::endl;

	std::string sbuf;
	int icount(0);
	while (std::getline(ifs, sbuf)) {
		CommonParameters::XYZ coord = { 0.0, 0.0, 0.0 };
		std::istringstream iss(sbuf);
		iss >> coord.X >> coord.Y >> coord.Z;
		if (coord.Z <= 0.0)
		{
			std::cerr << " Error : Lake bottom depth should be positive." << std::endl;
			exit(1);
		}
		m_coords.push_back(coord);
		++icount;
	}
	ifs.close();

	std::cout << "Total number of data is " << icount << std::endl;

}

// Check if the point is included in an area
bool Lake::isPointIncludedInAnArea(const double x, const double y) const {

	for (std::vector<Area>::const_iterator itr = m_listOfLakeAreas.begin(); itr != m_listOfLakeAreas.end(); ++itr) {
		if (x >= itr->XStart && x <= itr->XEnd && y >= itr->YStart && y <= itr->YEnd)
		{
			return true;
		}
	}
	return false;

}

// Get altitude of lake surface
double Lake::getAltitudeOfLakeSurface() const {
	return m_altitudeOfLakeSurface;
}

// Get lake resistivity
double Lake::getLakeResistivity() const {
	return m_lakeResistivity;
}

// Calculated average lake depth in a rectangular area
double Lake::calcAverageLakeDepth(const double xMin, const double xMax, const double yMin, const double yMax, int& numData) const
{
	numData = 0;
	double depth(0.0);
	const std::vector<CommonParameters::XYZ>::const_iterator itrEnd = m_coords.end();
	for (std::vector<CommonParameters::XYZ>::const_iterator itr = m_coords.begin(); itr != itrEnd; ++itr) {
		if (itr->X >= xMin && itr->X <= xMax && itr->Y >= yMin && itr->Y <= yMax) {
			++numData;
			depth += itr->Z;
		}
	}

	if (numData <= 0) {
		return -9999.0;
	}
	else
	{
		return depth / static_cast<double>(numData);
	}

}

// Calculated average Z coordinate in a rectangular area
double Lake::calcAverageZCoord(const double xMin, const double xMax, const double yMin, const double yMax) const {
	int numData(0);
	return calcAverageLakeDepth(xMin, xMax, yMin, yMax, numData) + getAltitudeOfLakeSurface();
}

// Select and add lake area
void Lake::selectLakeArea(const double xMin, const double xMax, const double yMin, const double yMax) {

	int numData(0);
	const double averageDepth = calcAverageLakeDepth(xMin, xMax, yMin, yMax, numData);
	if (averageDepth >= getThreshouldOfAverageLakeDepth() && numData >= getThreshouldOfNumOfDataPointInEachElementFace())
	{
		addArea(xMin, xMax, yMin, yMax);
	}

}

// Get threshould of average lake depth
double Lake::getThreshouldOfAverageLakeDepth() const {
	return m_threshouldOfAverageLakeDepth;
}

// Get threshould of number of data point in each element face
int Lake::getThreshouldOfNumOfDataPointInEachElementFace() const{
	return m_threshouldOfNumOfDataPointInEachElementFace;
}

// Get number of selected lake areas
int Lake::getNumOfSelectLakeAreas() const {
	return m_listOfLakeAreas.size();
}

// Add area
void Lake::addArea(const double xMin, const double xMax, const double yMin, const double yMax){

	Area obj;
	obj.XStart = xMin;
	obj.XEnd = xMax;
	obj.YStart = yMin;
	obj.YEnd = yMax;
	m_listOfLakeAreas.push_back(obj);

}

