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
#include "Cuboids.h"
#include "math.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <assert.h> 
#include <stdlib.h>

// Default constructer
Cuboids::Cuboids():
	m_rotationAngle(0.0),
	m_numCuboids(0),
	m_lengthX(NULL),
	m_lengthY(NULL),
	m_lengthZ(NULL),
	m_edgeLengthHorizontal(NULL),
	m_edgeLengthVertical(NULL)
{
	m_centerCoord.X = 0.0;
	m_centerCoord.Y = 0.0;
}

// Destructer
Cuboids::~Cuboids(){

	if( m_lengthX != NULL ){
		delete [] m_lengthX;
		m_lengthX = NULL;
	}

	if( m_lengthY != NULL ){
		delete [] m_lengthY;
		m_lengthY = NULL;
	}

	if( m_lengthZ != NULL ){
		delete [] m_lengthZ;
		m_lengthZ = NULL;
	}

	if( m_edgeLengthHorizontal != NULL ){
		delete [] m_edgeLengthHorizontal;
		m_edgeLengthHorizontal = NULL;
	}

	if( m_edgeLengthVertical != NULL ){
		delete [] m_edgeLengthVertical;
		m_edgeLengthVertical = NULL;
	}

}

// Read control parameters
void Cuboids::readParameters(std::ifstream& ifs) {

	std::cout << "Read data of ellipsoids used for specifing maximum edge length." << std::endl;

	ifs >> m_centerCoord.X >> m_centerCoord.Y >> m_centerCoord.Z;
	ifs >> m_rotationAngle;
	ifs >> m_numCuboids;
	if (m_numCuboids < 0) {
		std::cerr << "Error : Total number of cuboids (" << m_numCuboids << ") is less than 0 !!" << std::endl;
		exit(1);
	}
	if (m_numCuboids == 0) {
		return;
	}
	m_lengthX = new double[m_numCuboids];
	m_lengthY = new double[m_numCuboids];
	m_lengthZ = new double[m_numCuboids];
	m_edgeLengthHorizontal = new double[m_numCuboids];
	m_edgeLengthVertical = new double[m_numCuboids];
	for (int i = 0; i < m_numCuboids; ++i) {
#ifdef _VERTICAL_LENGTH_CTRL
		ifs >> m_lengthX[i] >> m_lengthY[i] >> m_lengthZ[i]
			>> m_edgeLengthHorizontal[i] >> m_edgeLengthVertical[i];
#else
		ifs >> m_lengthX[i] >> m_lengthY[i] >> m_lengthZ[i]
			>> m_edgeLengthHorizontal[i];
		m_edgeLengthVertical[i] = m_edgeLengthHorizontal[i];
#endif
	}
	for (int i = 1; i < m_numCuboids; ++i) {
		if (m_lengthX[i] < m_lengthX[i - 1]) {
			std::cerr << "Length along x-axis of the region " << i << " is smaller than that of the previous region." << std::endl;
			exit(1);
		}
		if (m_lengthY[i] < m_lengthY[i - 1]) {
			std::cerr << "Length along y-axis of the region " << i << " is smaller than that of the previous region." << std::endl;
			exit(1);
		}
		if (m_lengthZ[i] < m_lengthZ[i - 1]) {
			std::cerr << "Length along z-axis of the region " << i << " is smaller than that of the previous region." << std::endl;
			exit(1);
		}
		if (m_edgeLengthHorizontal[i] < m_edgeLengthHorizontal[i - 1]) {
			std::cerr << "Horizontal edge length of the region " << i << " is smaller than that of the previous region." << std::endl;
			exit(1);
		}
#ifdef _VERTICAL_LENGTH_CTRL
		if (m_edgeLengthVertical[i] < m_edgeLengthVertical[i - 1]) {
			std::cerr << "Vertical edge length of the region " << i << " is smaller than that of the previous region." << std::endl;
			exit(1);
		}
#endif
	}

	if( m_numCuboids > 0 ){
		std::cout << "Center coordinte [km] : " << m_centerCoord.X << " " << m_centerCoord.Y  << " " << m_centerCoord.Z  << std::endl;
		std::cout << "Rotation angle [deg] : " << m_rotationAngle << std::endl;
		std::cout << "Total number of cuboids : " << m_numCuboids << std::endl;
		std::cout.precision(6);
		for( int i = 0; i < m_numCuboids; ++i ){
			std::cout << " X-Length [km] : " << m_lengthX[i]
					<< " Y-Length [km] : " << m_lengthY[i]
					<< " Z-Length [km] : " << m_lengthZ[i]
					<< ", Edge length horizontal [km] : " << m_edgeLengthHorizontal[i]
#ifdef _VERTICAL_LENGTH_CTRL
					<< ", Edge length vertical [km] : " << m_edgeLengthVertical[i]
#endif
					<< std::endl;
		}
	}

	// Degrees => Radians
	m_rotationAngle *= CommonParameters::DEG2RAD;

}

// Calculate maximum horizontal length of the specified point
double Cuboids::calcMaximumEdgeLengthHorizontal( const CommonParameters::XYZ& coord ) const{

	double minVal = 1.0e20;

	if( m_numCuboids < 1 ){
		return minVal;
	}

	for( int iCuboids = 0; iCuboids < m_numCuboids; ++iCuboids ){
		if( locateInCuboids( coord, iCuboids ) ){
			return std::min( m_edgeLengthHorizontal[iCuboids], minVal );
			break;
		}
	}

	return minVal;

}

// Calculate maximum vertical length of the specified point
double Cuboids::calcMaximumEdgeLengthVertical( const CommonParameters::XYZ& coord ) const{

	double minVal = 1.0e20;

	if( m_numCuboids < 1 ){
		return minVal;
	}

	for( int iCuboids = 0; iCuboids < m_numCuboids; ++iCuboids ){
		if( locateInCuboids( coord, iCuboids ) ){
			return std::min( m_edgeLengthVertical[iCuboids], minVal );
			break;
		}
	}

	return minVal;

}

// Check whether specified point is located in the cuboid
bool Cuboids::locateInCuboids( const CommonParameters::XYZ& coord, const int iCuboid ) const{
	
	if( iCuboid < 0 || iCuboid >= m_numCuboids ){
		std::cerr << "Wrong cuboid ID:  " << iCuboid << std::endl;
		exit(1);
	}

	const double vecXOrg = coord.X - m_centerCoord.X;
	const double vecYOrg = coord.Y - m_centerCoord.Y;
	// Coordinate transform
	const double vecX = vecXOrg * cos( - m_rotationAngle ) - vecYOrg * sin( - m_rotationAngle );
	const double vecY = vecXOrg * sin( - m_rotationAngle ) + vecYOrg * cos( - m_rotationAngle );
	const double vecZ = coord.Z - m_centerCoord.Z;

	if( fabs(vecX) <= 0.5 * m_lengthX[iCuboid] && fabs(vecY) <= 0.5 * m_lengthY[iCuboid] && fabs(vecZ) <= 0.5 * m_lengthZ[iCuboid] ){
		return true;
	}

	return false;

}

// Get total number of cuboids
int Cuboids::getNumOfCuboids() const{

	return m_numCuboids;

}
