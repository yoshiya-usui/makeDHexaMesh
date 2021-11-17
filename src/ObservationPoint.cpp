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
#include "ObservationPoint.h"
#include "math.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <assert.h> 

// Default constructer
ObservationPoint::ObservationPoint():
	m_typeOfRegion(CUBOID),
	m_numRegions(0),
	m_length(NULL),
	m_maxEdgeLengthHorizontal(NULL),
	m_maxEdgeLengthVertical(NULL),
	m_oblateness(NULL)
{
	m_pointCoord.X = 0.0;
	m_pointCoord.Y = 0.0;
}

// Destructer
ObservationPoint::~ObservationPoint(){

	if( m_length != NULL ){
		delete [] m_length;
		m_length = NULL;
	}
	
	if( m_maxEdgeLengthHorizontal != NULL ){
		delete [] m_maxEdgeLengthHorizontal;
		m_maxEdgeLengthHorizontal = NULL;
	}
	
	if( m_maxEdgeLengthVertical != NULL ){
		delete [] m_maxEdgeLengthVertical;
		m_maxEdgeLengthVertical = NULL;
	}

	if( m_oblateness != NULL ){
		delete [] m_oblateness;
		m_oblateness = NULL;
	}

}

// Read data of observation point from input file
void ObservationPoint::readObservationPointData( std::ifstream& ifs ){

	ifs >> m_pointCoord.X >> m_pointCoord.Y >> m_pointCoord.Z;

	std::cout << "Coordinate : " << m_pointCoord.X << " " << m_pointCoord.Y << " " << m_pointCoord.Z << std::endl;

	ifs >> m_numRegions;
	if( m_numRegions <= 0 ){
		std::cerr << "Error : Total number of regions is less than 1. : " << m_numRegions << std::endl;
		exit(1);
	}
	std::cout << "Number of regions : " << m_numRegions << std::endl;

	m_length = new double[m_numRegions];
	m_maxEdgeLengthHorizontal = new double[m_numRegions];
	m_maxEdgeLengthVertical = new double[m_numRegions];
	m_oblateness = new double[m_numRegions];

	for( int i = 0; i < m_numRegions; ++i ){
#ifdef _VERTICAL_LENGTH_CTRL
		ifs >> m_length[i] >> m_maxEdgeLengthHorizontal[i] >> m_maxEdgeLengthVertical[i] >> m_oblateness[i];
#else
		ifs >> m_length[i] >> m_maxEdgeLengthHorizontal[i] >> m_oblateness[i];
		m_maxEdgeLengthVertical[i] = m_maxEdgeLengthHorizontal[i] ;
#endif
		if( i > 0 && ( m_length[i] < m_length[i-1] ) ){
			std::cerr << "Error : Inner length ( " << m_length[i-1] << " ) is greater than outer length ( " << m_length[i] << " ) !! " << std::endl;
			exit(1);
		}
		if( i > 0 && ( m_maxEdgeLengthHorizontal[i] < m_maxEdgeLengthHorizontal[i-1] ) ){
			std::cerr << "Error : Inner edge length ( " << m_maxEdgeLengthHorizontal[i-1] 
				<< " ) is greater than outer edge length ( " << m_maxEdgeLengthHorizontal[i] << " ) !! " << std::endl;
			exit(1);
		}
		if( i > 0 && ( m_maxEdgeLengthVertical[i] < m_maxEdgeLengthVertical[i-1] ) ){
			std::cerr << "Error : Inner edge length ( " << m_maxEdgeLengthVertical[i-1] 
				<< " ) is greater than outer edge length ( " << m_maxEdgeLengthVertical[i] << " ) !! " << std::endl;
			exit(1);
		}
		if( m_oblateness[i] >= 1.0 ){
			std::cerr << "Error : Oblateness must be smaller than 1." << std::endl;
			exit(1);
		}
	}

	for( int i = 0; i < m_numRegions; ++i ){
		std::cout << " length [km] : " << m_length[i]
				<< ", Edge length horizontal [km] : " << m_maxEdgeLengthHorizontal[i] 
#ifdef _VERTICAL_LENGTH_CTRL
				<< ", Edge length vertical [km] : " << m_maxEdgeLengthVertical[i] 
#endif
				<< ", Oblateness : " << m_oblateness[i] << std::endl;
	}
	
}

// Calculate maximum horizontal length
double ObservationPoint::calcMaximumHorizontalLength( const CommonParameters::XYZ& coord ) const{

	const double minVal = 1.0e20;
	
	if( m_typeOfRegion == SPHERE ){
		if( locateInRegion( coord, 0 ) ){
			return m_maxEdgeLengthHorizontal[0];
		}
		assert( m_numRegions > 0 );
		if( !locateInRegion( coord, m_numRegions - 1 ) ){
			return minVal;
		}
		assert( m_numRegions > 1 );
		int iRegion = m_numRegions - 2;
		for( ; iRegion >= 1; --iRegion ){
			if( !locateInRegion( coord, iRegion ) ){
				break;
			}
		}
		const double factor = calcFactorForIntermediateRegion( coord, iRegion );
		return m_maxEdgeLengthHorizontal[iRegion] + ( m_maxEdgeLengthHorizontal[iRegion+1] - m_maxEdgeLengthHorizontal[iRegion] ) * factor;
	}else if( m_typeOfRegion == CUBOID ){
		if( m_numRegions < 1 ){
			return minVal;
		}
		for( int iCuboids = 0; iCuboids < m_numRegions; ++iCuboids ){
			if( locateInRegion( coord, iCuboids ) ){
				return std::min( m_maxEdgeLengthHorizontal[iCuboids], minVal );
				break;
			}
		}
		return minVal;
	}

	return minVal;
}

// Calculate maximum vertical length
double ObservationPoint::calcMaximumVerticalLength( const CommonParameters::XYZ& coord ) const{

	const double minVal = 1.0e20;
	
	if( m_typeOfRegion == SPHERE ){
		if( locateInRegion( coord, 0 ) ){
			return m_maxEdgeLengthVertical[0];
		}
		assert( m_numRegions > 0 );
		if( !locateInRegion( coord, m_numRegions - 1 ) ){
			return minVal;
		}
		assert( m_numRegions > 1 );
		int iRegion = m_numRegions - 2;
		for( ; iRegion >= 1; --iRegion ){
			if( !locateInRegion( coord, iRegion ) ){
				break;
			}
		}
		const double factor = calcFactorForIntermediateRegion( coord, iRegion );
		return m_maxEdgeLengthVertical[iRegion] + ( m_maxEdgeLengthVertical[iRegion+1] - m_maxEdgeLengthVertical[iRegion] ) * factor;
	}else if( m_typeOfRegion == CUBOID ){
		if( m_numRegions < 1 ){
			return minVal;
		}
		for( int iCuboids = 0; iCuboids < m_numRegions; ++iCuboids ){
			if( locateInRegion( coord, iCuboids ) ){
				return std::min( m_maxEdgeLengthVertical[iCuboids], minVal );
				break;
			}
		}
		return minVal;
	}

	return minVal;
}

// Calculate factor for intermediate region
double ObservationPoint::calcFactorForIntermediateRegion( const CommonParameters::XYZ& coord, const int iRegion ) const{

	const double vecX = coord.X - m_pointCoord.X;
	const double vecY = coord.Y - m_pointCoord.Y;
	const double vecZ = coord.Z - m_pointCoord.Z;

	double val = pow( vecX / m_length[iRegion], 2 ) + pow( vecY / m_length[iRegion], 2 ) + pow( vecZ / ( m_length[iRegion] * ( 1.0 - m_oblateness[iRegion] ) ), 2 );
	val = sqrt(val);

	double val1 = sqrt( 3.0 * pow( m_length[iRegion+1] / m_length[iRegion], 2 ) );

	return ( val - 1.0 ) / ( val1 - 1.0 );

}

//// Judge the distance between the observation station and the specified point is lower than threshold value
//bool ObservationPoint::judgeDistanceObsSiteAndPointLowerThanThreshold( const CommonParameters::XYZ& coord, const double threshold ) const{
//
//	const double vecX = coord.X - m_pointCoord.X;
//	const double vecY = coord.Y - m_pointCoord.Y;
//	const double vecZ = coord.Z - m_pointCoord.Z;
//
//	double val = pow( vecX / threshold, 2 ) + pow( vecY / threshold, 2 ) + pow( vecZ / threshold, 2 );
//
//	if( val < 1.0 ){
//		return true;
//	}
//
//	return false;
//
//}

bool ObservationPoint::locateInRegion( const CommonParameters::XYZ& coord, const int iRegion ) const{

	assert( iRegion >= 0 || iRegion < m_numRegions );

	const double vecX = coord.X - m_pointCoord.X;
	const double vecY = coord.Y - m_pointCoord.Y;
	const double vecZ = coord.Z - m_pointCoord.Z;
	if( m_typeOfRegion == SPHERE ){
		const double depth = m_length[iRegion] * ( 1.0 - m_oblateness[iRegion] );
		double val = pow( vecX / m_length[iRegion], 2 )
				   + pow( vecY / m_length[iRegion], 2 )
				   + pow( vecZ / depth, 2 );
		if( val <= 1.0 ){
			return true;
		}
		return false;
	}else if( m_typeOfRegion == CUBOID ){
		if( fabs(vecX) <= 0.5 * m_length[iRegion] &&
			fabs(vecY) <= 0.5 * m_length[iRegion] &&
			fabs(vecZ) <= 0.5 * m_length[iRegion] * ( 1.0 - m_oblateness[iRegion] ) ){
			return true;
		}
		return false;
	}

	return false;

}
