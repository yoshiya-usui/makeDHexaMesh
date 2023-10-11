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
#include "Ellipsoids.h"
#include "math.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <assert.h> 
#include <stdlib.h>

// Default constructer
Ellipsoids::Ellipsoids():
	m_rotationAngle(0.0),
	m_numEllipsoids(0),
	m_radius(NULL),
	m_edgeLengthHorizontal(NULL),
	m_edgeLengthVertical(NULL),
	m_oblatenessHorizontal(NULL)
{

	m_centerCoord.X = 0.0;
	m_centerCoord.Y = 0.0;

	m_oblateness[0] = NULL;
	m_oblateness[1] = NULL;
}

// Destructer
Ellipsoids::~Ellipsoids(){

	if( m_radius != NULL ){
		delete [] m_radius;
		m_radius = NULL;
	}

	if( m_edgeLengthHorizontal != NULL ){
		delete [] m_edgeLengthHorizontal;
		m_edgeLengthHorizontal = NULL;
	}

	if( m_edgeLengthVertical != NULL ){
		delete [] m_edgeLengthVertical;
		m_edgeLengthVertical = NULL;
	}

	if( m_oblatenessHorizontal != NULL ){
		delete [] m_oblatenessHorizontal;
		m_oblatenessHorizontal = NULL;
	}

	if( m_oblateness[0] != NULL ){
		delete [] m_oblateness[0];
		m_oblateness[0] = NULL;
	}

	if( m_oblateness[1] != NULL ){
		delete [] m_oblateness[1];
		m_oblateness[1] = NULL;
	}

}

// Read control parameters
void Ellipsoids::readParameters(std::ifstream& ifs) {

	std::cout << "Read data of ellipsoids used for specifing maximum edge length." << std::endl;

	ifs >> m_centerCoord.X >> m_centerCoord.Y >> m_centerCoord.Z;
	ifs >> m_rotationAngle;
	ifs >> m_numEllipsoids;

	if (m_numEllipsoids < 0) {
		std::cerr << "Error : Total number of ellipsoids (" << m_numEllipsoids << ") is less than 0 !!" << std::endl;
		exit(1);
	}
	if (m_numEllipsoids == 0) {
		return;
	}

	m_radius = new double[m_numEllipsoids];
	m_edgeLengthHorizontal = new double[m_numEllipsoids];
	m_edgeLengthVertical = new double[m_numEllipsoids];
	m_oblatenessHorizontal = new double[m_numEllipsoids];
	m_oblateness[Ellipsoids::EARTH] = new double[m_numEllipsoids];
	m_oblateness[Ellipsoids::AIR] = new double[m_numEllipsoids];
	for (int i = 0; i < m_numEllipsoids; ++i) {
#ifdef _VERTICAL_LENGTH_CTRL
		ifs >> m_radius[i] >> m_edgeLengthHorizontal[i] >> m_edgeLengthVertical[i]
			>> m_oblatenessHorizontal[i] >> m_oblateness[Ellipsoids::EARTH][i] >> m_oblateness[Ellipsoids::AIR][i];
#else
		ifs >> m_radius[i] >> m_edgeLengthHorizontal[i]
			>> m_oblatenessHorizontal[i] >> m_oblateness[Ellipsoids::EARTH][i] >> m_oblateness[Ellipsoids::AIR][i];
		m_edgeLengthVertical[i] = m_edgeLengthHorizontal[i];
#endif		
	}
	for (int i = 1; i < m_numEllipsoids; ++i) {
		if (m_radius[i] < m_radius[i - 1]) {
			std::cerr << "Radius of the region " << i << " is smaller than that of the previous region." << std::endl;
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
		if (m_oblatenessHorizontal[i] < 0 || m_oblatenessHorizontal[i] > 1) {
			std::cerr << "Oblateness of horizontal ellipsoid must be smaller than 1 and larger than 0." << std::endl;
			exit(1);
		}
		if (m_oblateness[Ellipsoids::EARTH][i] >= 1.0) {
			std::cerr << "Oblateness must be smaller than 1." << std::endl;
			exit(1);
		}
		if (m_oblateness[Ellipsoids::AIR][i] >= 1.0) {
			std::cerr << "Oblateness must be smaller than 1." << std::endl;
			exit(1);
		}
		if (m_radius[i] * (1.0 - m_oblateness[Ellipsoids::EARTH][i]) < m_radius[i - 1] * (1.0 - m_oblateness[Ellipsoids::EARTH][i - 1])) {
			std::cerr << "Depth of sphere " << i << " is shallower than that of the previous sphere in the earth." << std::endl;
			exit(1);
		}
		if (m_radius[i] * (1.0 - m_oblateness[Ellipsoids::AIR][i]) < m_radius[i - 1] * (1.0 - m_oblateness[Ellipsoids::AIR][i - 1])) {
			std::cerr << "Depth of sphere " << i << " is shallower than that of the previous sphere in the air." << std::endl;
			exit(1);
		}
	}

	if( m_numEllipsoids > 0 ){
		std::cout << "Center coordinte [km] : " << m_centerCoord.X << " " << m_centerCoord.Y  << " " << m_centerCoord.Z  << std::endl;
		std::cout << "Rotation angle [deg] : " << m_rotationAngle << std::endl;
		std::cout << "Total number of ellipsoids : " << m_numEllipsoids << std::endl;
		std::cout.precision(6);
		for( int i = 0; i < m_numEllipsoids; ++i ){
			std::cout << " Radius [km] : " << m_radius[i]
					<< ", Edge length horizontal [km] : " << m_edgeLengthHorizontal[i]
#ifdef _VERTICAL_LENGTH_CTRL
					<< ", Edge length vertical [km] : " << m_edgeLengthVertical[i]
#endif
					<< ", Oblateness horizontal : " << m_oblatenessHorizontal[i]
					<< ", Oblateness lower : " << m_oblateness[Ellipsoids::EARTH][i]
					<< ", Oblateness upper : " << m_oblateness[Ellipsoids::AIR][i]
					<< std::endl;
		}
	}

	// Degrees => Radians
	m_rotationAngle *= CommonParameters::DEG2RAD;

}

// Calculate maximum horizontal length of the specified point
double Ellipsoids::calcMaximumEdgeLengthHorizontal( const CommonParameters::XYZ& coord ) const{

	double minVal = 1.0e20;

	if( m_numEllipsoids < 1 ){
		return minVal;
	}

	if( locateInEllipsoid( coord, 0 ) ){
		return std::min( m_edgeLengthHorizontal[0], minVal );
	}

	for( int iSphere = 1; iSphere < m_numEllipsoids; ++iSphere ){
		if( locateInEllipsoid( coord, iSphere ) ){
			return std::min( calcEdgeLengthTransitionHorizontal( coord, iSphere ), minVal );
			break;
		}
	}

	return minVal;

}

// Calculate maximum vertical length of the specified point
double Ellipsoids::calcMaximumEdgeLengthVertical( const CommonParameters::XYZ& coord ) const{

	double minVal = 1.0e20;

	if( m_numEllipsoids < 1 ){
		return minVal;
	}

	if( locateInEllipsoid( coord, 0 ) ){
		return std::min( m_edgeLengthVertical[0], minVal );
	}

	for( int iSphere = 1; iSphere < m_numEllipsoids; ++iSphere ){
		if( locateInEllipsoid( coord, iSphere ) ){
			return std::min( calcEdgeLengthTransitionVertical( coord, iSphere ), minVal );
			break;
		}
	}

	return minVal;

}

// Check whether specified point is located in the Ellipsoid
bool Ellipsoids::locateInEllipsoid( const CommonParameters::XYZ& coord, const int iEllipsoid ) const{
	
	if( iEllipsoid < 0 || iEllipsoid >= m_numEllipsoids ){
		std::cerr << "Wrong sphere ID:  " << iEllipsoid << std::endl;
		exit(1);
	}

	const double vecXOrg = coord.X - m_centerCoord.X;
	const double vecYOrg = coord.Y - m_centerCoord.Y;
	// Coordinate transform
	const double vecX = vecXOrg * cos( - m_rotationAngle ) - vecYOrg * sin( - m_rotationAngle );
	const double vecY = vecXOrg * sin( - m_rotationAngle ) + vecYOrg * cos( - m_rotationAngle );
	const double vecZ = coord.Z - m_centerCoord.Z;

	const int iType = vecZ < 0 ? Ellipsoids::AIR : Ellipsoids::EARTH;

	const double longAxisLength = m_radius[iEllipsoid];
	const double shortAxisLength = longAxisLength * ( 1.0 - m_oblatenessHorizontal[iEllipsoid] );
	const double depth = longAxisLength * ( 1.0 - m_oblateness[iType][iEllipsoid] );

	double val = pow( vecX / longAxisLength, 2 )
			   + pow( vecY / shortAxisLength, 2 )
			   + pow( vecZ / depth, 2 );

	if( val <= 1.0 ){
		return true;
	}

	return false;

}

// Calculate maximum horizontal length between two spheres
double Ellipsoids::calcEdgeLengthTransitionHorizontal( const CommonParameters::XYZ& coord, const int iEllipsoid ) const{

	const double factor = calcFactorForEdgeLengthTransitionRegion( coord, iEllipsoid );
	return factor * ( m_edgeLengthHorizontal[iEllipsoid] - m_edgeLengthHorizontal[iEllipsoid-1] ) + m_edgeLengthHorizontal[iEllipsoid-1];

}

// Calculate maximum vertical length between two spheres
double Ellipsoids::calcEdgeLengthTransitionVertical( const CommonParameters::XYZ& coord, const int iEllipsoid ) const{

	const double factor = calcFactorForEdgeLengthTransitionRegion( coord, iEllipsoid );
	return factor * ( m_edgeLengthVertical[iEllipsoid] - m_edgeLengthVertical[iEllipsoid-1] ) + m_edgeLengthVertical[iEllipsoid-1];

}

// Calculate factor for the maximum edge length between two spheres
double Ellipsoids::calcFactorForEdgeLengthTransitionRegion( const CommonParameters::XYZ& coord, const int iEllipsoid ) const{

	if( iEllipsoid < 1 || iEllipsoid >= m_numEllipsoids ){
		std::cerr << "Wrong sphere ID:  " << iEllipsoid << std::endl;
		exit(1);
	}

	const double vecXOrg = coord.X - m_centerCoord.X;
	const double vecYOrg = coord.Y - m_centerCoord.Y;
	// Coordinate transform
	const double vecX = vecXOrg * cos( - m_rotationAngle ) - vecYOrg * sin( - m_rotationAngle );
	const double vecY = vecXOrg * sin( - m_rotationAngle ) + vecYOrg * cos( - m_rotationAngle );
	const double vecZ = coord.Z - m_centerCoord.Z;

	const int iType = vecZ < 0 ? Ellipsoids::AIR : Ellipsoids::EARTH;

	const double angleHorizontal = atan2( vecY, vecX );
	const double lengthHorizontal = hypot( vecY, vecX );
	const double angleVertical = atan2( vecZ, lengthHorizontal );
	double length = hypot( lengthHorizontal, vecZ );

	const double length0 = calculateLengthOnEllipsoid( angleHorizontal, angleVertical, iEllipsoid-1, iType );
	const double length1 = calculateLengthOnEllipsoid( angleHorizontal, angleVertical, iEllipsoid, iType );

	if( length < length0 ){
		length = length0;
	}
	if( length > length1 ){
		length = length1;
	}

	return ( length - length0 ) / ( length1 - length0 );

}

// Calculate length on ellipsoid
double Ellipsoids::calculateLengthOnEllipsoid( const double angleH, const double angleV, const int iEllipsoid, const int iType ) const{

	if( iEllipsoid < 0 || iEllipsoid >= m_numEllipsoids ){
		std::cerr << "Wrong sphere ID:  " << iEllipsoid << std::endl;
		exit(1);
	}

	if( angleH < -CommonParameters::PI || angleH > CommonParameters::PI ){
		std::cerr << "Horizontal angle is improper : " << angleH << std::endl;
		exit(1);
	}

	if( angleV < -CommonParameters::PI || angleV > CommonParameters::PI ){
		std::cerr << "Vertical angle is improper : " << angleV << std::endl;
		exit(1);
	}

	const double longAxisLength = m_radius[iEllipsoid];
	const double shortAxisLength = longAxisLength * ( 1.0 - m_oblatenessHorizontal[iEllipsoid] );
	const double verticalLength = longAxisLength * ( 1.0 - m_oblateness[iType][iEllipsoid] );

	const double eps = 1.0e-9;
	double lengthH(-1.0);
	if( fabs(angleH - CommonParameters::PI*0.5) < eps || fabs(angleH + CommonParameters::PI*0.5) < eps ){
		lengthH = shortAxisLength;
	}
	else{
		const double constValH = longAxisLength * shortAxisLength / hypot( shortAxisLength, longAxisLength * tan(angleH) );
		lengthH = hypot( constValH, constValH * tan(angleH) );
	}
	
	if( fabs(angleV - CommonParameters::PI*0.5) < eps || fabs(angleV + CommonParameters::PI*0.5) < eps ){
		return verticalLength;
	}
	else{
		const double constValV = lengthH * verticalLength / hypot( verticalLength, lengthH * tan(angleV) );
		return hypot( constValV, constValV * tan(angleV) );
	}

}

// Get total number of ellipsoids
int Ellipsoids::getNumOfEllipsoids() const{

	return m_numEllipsoids;

}
