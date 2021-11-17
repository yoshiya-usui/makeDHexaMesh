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
#include "TopographyData.h"
#include "math.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>

// Default constructer
TopographyData::TopographyData():
	m_maxNumOfPointsForInterpolating(3),
	m_maxDistanceForInterpolating(100.0),
	m_distanceUsedToAvoidTooSmallDenominator(1.0e-6)
{
}

// Constructer
TopographyData::TopographyData( const std::string& inFileName, const int maxNumPoint, const double maxDist, const double eps ):
	m_fileNameOfTopographytData(inFileName),
	m_maxNumOfPointsForInterpolating(maxNumPoint),
	m_maxDistanceForInterpolating(maxDist),
	m_distanceUsedToAvoidTooSmallDenominator(eps)
{
}

// Destructer
TopographyData::~TopographyData(){
}

// Read topography data
void TopographyData::readTopographyData( const double xMin, const double xMax, const double yMin, const double yMax ){

	// Open input file -----
	std::ifstream ifs( m_fileNameOfTopographytData.c_str() );
	if( !ifs ) {
		std::cerr << "Cannot open file " << m_fileNameOfTopographytData << std::endl;
		exit(1);
	}
	std::cout << "Read data from " << m_fileNameOfTopographytData << std::endl;
	//----------------------

	std::string sbuf;
	int icount(0);
    while( std::getline(ifs, sbuf) ){
		CommonParameters::XYZ coord = { 0.0, 0.0, 0.0 };
		std::istringstream iss( sbuf );
		iss >> coord.X >> coord.Y >> coord.Z;
		if( coord.X < xMin - m_maxDistanceForInterpolating ||
			coord.X > xMax + m_maxDistanceForInterpolating ||
			coord.Y < yMin - m_maxDistanceForInterpolating ||
			coord.Y > yMax + m_maxDistanceForInterpolating ){
			continue;
		}
		m_coords.push_back( coord );
		++icount;
    }
	ifs.close();

	std::cout << "Total number of data is " << icount << std::endl;

}

#ifdef _TOPO_FUNC
// Calculate z coordinate by a function
double TopographyData::calcZByFunction( const CommonParameters::XY& coord ) const{

	const double xOrg = coord.X;
	const double yOrg = coord.Y;

	const bool includeSea(true);
	if( includeSea ){
//		const double zSeaBot = m_seaDepth;
//#if 0
//		// Schwalenberg & Edwards 2004
//		const double pi = acos(-1.0);
//		const double deg = xOrg * 360.0;
//		const double rad = deg * pi / 180.0;
//		const double depth = 0.1 * cos(rad) + zSeaBot;
//#else
//		// Gaussian sea mount
//		const double sigma = 0.25;
//		const double maxHeight = 0.2;
//		const double factor = maxHeight * sqrt(2.0*CommonParameters::PI) * sigma;
//		const double dist = hypot(xOrg, yOrg);
//		const double height = factor / ( sqrt(2.0*CommonParameters::PI) * sigma ) * exp( - 0.5 * pow(dist/sigma,2) );
//		const double depth = zSeaBot - height;
//#endif

		const double maxDepth = 3.0;
		double seaDepth = -1.0;
		const double slopeH =  10.0;
		const double coastX =  55.0;
		const double coastY = -55.0;
		const double radius = hypot( coord.X - 45.0, coord.Y + 45.0 );
		bool locateInDeepSea(false);
		if( coord.X >= coastX - slopeH && coord.Y <= coastY + slopeH && radius <= slopeH ){
			seaDepth = 0.5 * maxDepth - 0.5 * maxDepth * cos( CommonParameters::PI * ( slopeH - radius ) / slopeH );
		}else if( coord.Y >= -45.0 && coord.X <= coastX && coord.X >= coastX - slopeH ){
			seaDepth = 0.5 * maxDepth - 0.5 * maxDepth * cos( - CommonParameters::PI * ( coord.X - coastX ) / slopeH );
		}else if( coord.X <=  45.0 && coord.Y >= coastY && coord.Y <= coastY + slopeH ){
			seaDepth = 0.5 * maxDepth - 0.5 * maxDepth * cos( CommonParameters::PI * ( coord.Y - coastY ) / slopeH );
		}else if( coord.X <= coastX - slopeH && coord.Y >= coastY + slopeH  ){
			seaDepth = maxDepth;
			locateInDeepSea = true;
		}
		if(locateInDeepSea){
			const double waveLengthSmall = 3.0;
			const double amplitudeSmall = 0.1;
			const double radianX = 2.0 * CommonParameters::PI * coord.X / waveLengthSmall;
			const double radianY = 2.0 * CommonParameters::PI * coord.Y / waveLengthSmall;
			seaDepth -= amplitudeSmall * sin(radianX) * sin(radianY);
		}
		if(locateInDeepSea){
			const double waveLengthSmall = 10.0;
			const double amplitudeSmall = 0.25;
			const double radianX = 2.0 * CommonParameters::PI * coord.X / waveLengthSmall;
			const double radianY = 2.0 * CommonParameters::PI * coord.Y / waveLengthSmall;
			seaDepth -= amplitudeSmall * sin(radianX) * sin(radianY);
		}
		if(locateInDeepSea){
			const double waveLengthSmall = 30.0;
			const double amplitudeSmall = 0.5;
			const double radianX = 2.0 * CommonParameters::PI * coord.X / waveLengthSmall;
			const double radianY = 2.0 * CommonParameters::PI * coord.Y / waveLengthSmall;
			seaDepth -= amplitudeSmall * sin(radianX) * sin(radianY);
		}
		return seaDepth;
	}else{
		double height(0.0);
//#if 1
//		const double length = 0.45 * 0.5;
//		if( xOrg > 1.0 || xOrg < -1.0 || yOrg > 1.0 || yOrg < -1.0 ){
//			return zOrg;
//		}
//		if( xOrg < length && xOrg > -length && yOrg < length && yOrg > -length ){
//			height = 0.45;
//		}else if( fabs(xOrg) >= fabs(yOrg) ){
//			if( xOrg > 0.0 ){
//				height = 0.45 + (xOrg - length) / (1.0 - length) * (0.0 - 0.45);
//			}else{
//				height = 0.0 + (xOrg + 1.0) / (-length + 1.0) * (0.45 - 0.0);
//			}
//		}
//		else{
//			if( yOrg > 0.0 ){
//				height = 0.45 + (yOrg - length) / (1.0 - length) * (0.0 - 0.45);
//			}else{
//				height = 0.0 + (yOrg + 1.0) / (-length + 1.0) * (0.45 - 0.0);
//			}
//		}
//#else
//		const double sigma = 20.0;
//		const double maxHeight = 10.0;
//		const double factor = maxHeight * sqrt(2.0*CommonParameters::PI)*sigma;
//		const double dist = hypot(xOrg, yOrg);
//		const double height = factor / ( sqrt( 2.0 * CommonParameters::PI ) * sigma ) * exp( - 0.5 * pow(dist/sigma,2) );
//#endif
		return -height;
	}

}
#endif

// Interpolate Z coordinate by inverse distance weighting
double TopographyData::interpolateZCoord( const CommonParameters::XY& coord ) const{

#ifdef _TOPO_FUNC
	return calcZByFunction(coord);
#endif

	std::vector< std::pair<double,int> > stackDistances;
	for( int i = 0; i < m_maxNumOfPointsForInterpolating; ++i ){
		stackDistances.push_back( std::make_pair( 1.0e12, -1 ) );// Initialize 
	}

	int iDataMinDist(-1);
	double minDist(1.0e12);

	int iData(0);
	const std::vector<CommonParameters::XYZ>::const_iterator itrEnd = m_coords.end();
	for( std::vector<CommonParameters::XYZ>::const_iterator itr = m_coords.begin(); itr != itrEnd; ++itr, ++iData ){
		const double distance = hypot( itr->X - coord.X, itr->Y - coord.Y );
		if( distance < minDist ){
			iDataMinDist = iData;
			minDist = distance;
		}
		if( distance > m_maxDistanceForInterpolating ){
			continue;
		}
		if( distance < stackDistances.back().first ){
			stackDistances.back().first = distance;
			stackDistances.back().second = iData;
			sort( stackDistances.begin(), stackDistances.end() );
		}
	}

	if( stackDistances.front().second < 0 ){
		std::cerr << " Warning : No topography data were found near the point (X,Y) = ( " << coord.X << ", " << coord.Y << " )" << std::endl;
		if( iDataMinDist < 0 ){
			std::cerr << " Error : Point serial of the nearest point is negative !!" << std::endl;
			exit(1);
		}
		std::cerr << "           Thus, the value of the nearest point (distance = " << minDist << ", value = " << m_coords[iDataMinDist].Z << ") is given." << std::endl;
		return m_coords[iDataMinDist].Z;
	}

	double zCoord = 0.0;
	double weightSum(0.0);
	for( std::vector< std::pair<double,int> >::iterator itr = stackDistances.begin(); itr != stackDistances.end(); ++itr ){
		if( itr->second < 0 ){
			continue;
		}
		const double weight = 1.0 / ( m_distanceUsedToAvoidTooSmallDenominator + itr->first );
		const double z = m_coords[itr->second].Z;
		zCoord += weight * z;
		weightSum += weight;
	}
	zCoord /= weightSum;
	return zCoord;

}

// Calculated average Z coordinate in a rectangular area
double TopographyData::calcAverageZCoord( const double xMin, const double xMax, const double yMin, const double yMax ){

#ifdef _TOPO_FUNC
	double z(0.0);
	const int numDiv = 5;
	const double dx = ( xMax - xMin ) / static_cast<double>(numDiv);
	const double dy = ( yMax - yMin ) / static_cast<double>(numDiv);
	int icount(0);
	for( int ix = 0; ix <= numDiv; ++ix ){
		const double x = xMin + dx * static_cast<double>(ix);
		for( int iy = 0; iy <= numDiv; ++iy ){
			const double y = yMin + dy * static_cast<double>(iy);
			const CommonParameters::XY coord = {x, y};
			z += calcZByFunction(coord);
			++icount;
		}
	}
	z /= static_cast<double>(icount);
	return z;
#endif

	int numData(0);
	double depth(0.0);
	const std::vector<CommonParameters::XYZ>::const_iterator itrEnd = m_coords.end();
	for( std::vector<CommonParameters::XYZ>::const_iterator itr = m_coords.begin(); itr != itrEnd; ++itr ){
		if( itr->X >= xMin && itr->X <= xMax && itr->Y >= yMin && itr->Y <= yMax ){
			++numData;
			depth += itr->Z;
		}
	}

	if( numData == 0 ){
		std::cerr << " Warning : Topography data cannot be found in the area (xMin,xMax,yMin,yMax) = ( "
			<< xMin << "," << xMax << "," << yMin << "," << yMax << " )" << std::endl;
		const CommonParameters::XY coord = { 0.5*(xMin+xMax), 0.5*(yMin+yMax) }; 
		return interpolateZCoord(coord);
	}

	return depth / static_cast<double>(numData);

}

void TopographyData::outputVTK( const std::string& fname ) const{

	std::ofstream ofsVTK( fname.c_str() );
	if( !ofsVTK ) {
		std::cerr << "Cannot open file " << fname << std::endl;
		exit(1);
	}

	ofsVTK << "# vtk DataFile Version 2.0" << std::endl;
	ofsVTK << "Altitude" << std::endl;
	ofsVTK << "ASCII" << std::endl;
	ofsVTK << "DATASET UNSTRUCTURED_GRID" << std::endl;

	// output data to vtk file -----
	const int numPoints = static_cast<int>( m_coords.size() );

	ofsVTK.precision(9);
	ofsVTK << std::fixed;

	ofsVTK << "POINTS " << numPoints << " double" << std::endl;
	ofsVTK.precision(6);
	for( std::vector<CommonParameters::XYZ>::const_iterator itr = m_coords.begin(); itr != m_coords.end(); ++itr ){
		ofsVTK << std::setw(15) << std::scientific << itr->X;
		ofsVTK << std::setw(15) << std::scientific << itr->Y;
		ofsVTK << std::setw(15) << std::scientific << itr->Z << std::endl;
	}

	ofsVTK << "CELLS " << numPoints << " " << numPoints*2 << std::endl;
	for( int i = 0; i < numPoints; ++i ){
		ofsVTK << std::setw(10) << 1;
		ofsVTK << std::setw(10) << i << std::endl;
	}

	ofsVTK << "CELL_TYPES " << numPoints << std::endl;
	for( int i = 0; i < numPoints; ++i ){
		ofsVTK << std::setw(10) << 1;
	}

	ofsVTK << "POINT_DATA " << numPoints << std::endl;
	ofsVTK << "SCALARS Height(km) float" << std::endl;
	ofsVTK << "LOOKUP_TABLE default" << std::endl;
	ofsVTK.precision(6);
	for( std::vector<CommonParameters::XYZ>::const_iterator itr = m_coords.begin(); itr != m_coords.end(); ++itr ){
		ofsVTK << std::setw(15) << std::scientific << itr->Z << std::endl;
	}
	//------------------------------

	ofsVTK.close();
}
