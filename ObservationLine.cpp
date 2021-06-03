﻿//--------------------------------------------------------------------------
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
#include "ObservationLine.h"
#include "math.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>

// Default constructer
ObservationLine::ObservationLine():
	m_numLayers(0),
	m_radius(NULL),
	m_maxEdgeLengthHorizontal(NULL),
	m_maxEdgeLengthVertical(NULL)
{
}

// Destructer
ObservationLine::~ObservationLine(){
	if( m_radius != NULL ){
		delete [] m_radius;
		m_radius = NULL;
	}
	
	if( m_maxEdgeLengthHorizontal != NULL ){
		delete [] m_maxEdgeLengthHorizontal;
		m_maxEdgeLengthHorizontal = NULL;
	}

	if( m_maxEdgeLengthVertical != NULL ){
		delete [] m_maxEdgeLengthVertical;
		m_maxEdgeLengthVertical = NULL;
	}

}

// Calculate maximum horizontal length
double ObservationLine::calcMaximumHorizontalLength( const CommonParameters::XYZ& coord ) const{

	double edgeLeng = 1.0e+20;
	for( int iPoint = 0; iPoint < 2; ++iPoint ){// Loop of the end points
		const double vecX = coord.X - m_endPoints[iPoint].X;
		const double vecY = coord.Y - m_endPoints[iPoint].Y;
		const double vecZ = coord.Z - m_endPoints[iPoint].Z;

		const double radius = sqrt( pow( vecX, 2.0 ) + pow( vecY, 2.0 ) + pow( vecZ, 2.0 ) );

		if( radius <= m_radius[0] ){
			edgeLeng = std::min( m_maxEdgeLengthHorizontal[0], edgeLeng );
			break;
		}

		for( int i = m_numLayers - 1; i >= 1; --i ){
			if( radius <= m_radius[i] ){
				 edgeLeng = m_maxEdgeLengthHorizontal[i-1] + ( radius - m_radius[i-1] ) / ( m_radius[i] - m_radius[i-1] ) * ( m_maxEdgeLengthHorizontal[i] - m_maxEdgeLengthHorizontal[i-1] );
			}else{
				break;
			}
		}
	}

	const double lineLength = sqrt( 
		pow(m_endPoints[1].X - m_endPoints[0].X, 2) + 
		pow(m_endPoints[1].Y - m_endPoints[0].Y, 2) + 
		pow(m_endPoints[1].Z - m_endPoints[0].Z, 2) );
	const double innerProduct =
		(coord.X - m_endPoints[0].X) * (m_endPoints[1].X - m_endPoints[0].X) + 
		(coord.Y - m_endPoints[0].Y) * (m_endPoints[1].Y - m_endPoints[0].Y) + 
		(coord.Z - m_endPoints[0].Z) * (m_endPoints[1].Z - m_endPoints[0].Z);
	const double innerProductDiv = innerProduct / lineLength;
	if( innerProductDiv >= 0.0 && innerProductDiv <= lineLength ){// Within the two end points
		const CommonParameters::XYZ vectorAlongLine = {
			innerProductDiv / lineLength * (m_endPoints[1].X - m_endPoints[0].X),
			innerProductDiv / lineLength * (m_endPoints[1].Y - m_endPoints[0].Y),
			innerProductDiv / lineLength * (m_endPoints[1].Z - m_endPoints[0].Z) };
		const CommonParameters::XYZ vectorPerpendicularToLine = {
			coord.X - m_endPoints[0].X - vectorAlongLine.X,
			coord.Y - m_endPoints[0].Y - vectorAlongLine.Y,
			coord.Z - m_endPoints[0].Z - vectorAlongLine.Z };
		const double distanceFromLine = sqrt(pow(vectorPerpendicularToLine.X, 2) + pow(vectorPerpendicularToLine.Y, 2) + pow(vectorPerpendicularToLine.Z, 2) );
		if( distanceFromLine <= m_radius[0] ){
			return std::min( m_maxEdgeLengthHorizontal[0], edgeLeng );
		}

		for( int i = m_numLayers - 1; i >= 1; --i ){
			if( distanceFromLine <= m_radius[i] ){
				edgeLeng = m_maxEdgeLengthHorizontal[i-1] + ( distanceFromLine - m_radius[i-1] ) / ( m_radius[i] - m_radius[i-1] ) * ( m_maxEdgeLengthHorizontal[i] - m_maxEdgeLengthHorizontal[i-1] );
			}else{
				break;
			}
		}
	}
	return edgeLeng;

}

// Calculate maximum vertical length
double ObservationLine::calcMaximumVerticalLength( const CommonParameters::XYZ& coord ) const{

	double edgeLeng = 1.0e+20;
	for( int iPoint = 0; iPoint < 2; ++iPoint ){// Loop of the end points
		const double vecX = coord.X - m_endPoints[iPoint].X;
		const double vecY = coord.Y - m_endPoints[iPoint].Y;
		const double vecZ = coord.Z - m_endPoints[iPoint].Z;

		const double radius = sqrt( pow( vecX, 2.0 ) + pow( vecY, 2.0 ) + pow( vecZ, 2.0 ) );

		if( radius <= m_radius[0] ){
			edgeLeng = std::min( m_maxEdgeLengthVertical[0], edgeLeng );
			break;
		}

		for( int i = m_numLayers - 1; i >= 1; --i ){
			if( radius <= m_radius[i] ){
				 edgeLeng = m_maxEdgeLengthVertical[i-1] + ( radius - m_radius[i-1] ) / ( m_radius[i] - m_radius[i-1] ) * ( m_maxEdgeLengthVertical[i] - m_maxEdgeLengthVertical[i-1] );
			}else{
				break;
			}
		}
	}

	const double lineLength = sqrt( 
		pow(m_endPoints[1].X - m_endPoints[0].X, 2) + 
		pow(m_endPoints[1].Y - m_endPoints[0].Y, 2) + 
		pow(m_endPoints[1].Z - m_endPoints[0].Z, 2) );
	const double innerProduct =
		(coord.X - m_endPoints[0].X) * (m_endPoints[1].X - m_endPoints[0].X) + 
		(coord.Y - m_endPoints[0].Y) * (m_endPoints[1].Y - m_endPoints[0].Y) + 
		(coord.Z - m_endPoints[0].Z) * (m_endPoints[1].Z - m_endPoints[0].Z);
	const double innerProductDiv = innerProduct / lineLength;
	if( innerProductDiv >= 0.0 && innerProductDiv <= lineLength ){// Within the two end points
		const CommonParameters::XYZ vectorAlongLine = {
			innerProductDiv / lineLength * (m_endPoints[1].X - m_endPoints[0].X),
			innerProductDiv / lineLength * (m_endPoints[1].Y - m_endPoints[0].Y),
			innerProductDiv / lineLength * (m_endPoints[1].Z - m_endPoints[0].Z) };
		const CommonParameters::XYZ vectorPerpendicularToLine = {
			coord.X - m_endPoints[0].X - vectorAlongLine.X,
			coord.Y - m_endPoints[0].Y - vectorAlongLine.Y,
			coord.Z - m_endPoints[0].Z - vectorAlongLine.Z };
		const double distanceFromLine = sqrt(pow(vectorPerpendicularToLine.X, 2) + pow(vectorPerpendicularToLine.Y, 2) + pow(vectorPerpendicularToLine.Z, 2) );
		if( distanceFromLine <= m_radius[0] ){
			return std::min( m_maxEdgeLengthVertical[0], edgeLeng );
		}

		for( int i = m_numLayers - 1; i >= 1; --i ){
			if( distanceFromLine <= m_radius[i] ){
				edgeLeng = m_maxEdgeLengthVertical[i-1] + ( distanceFromLine - m_radius[i-1] ) / ( m_radius[i] - m_radius[i-1] ) * ( m_maxEdgeLengthVertical[i] - m_maxEdgeLengthVertical[i-1] );
			}else{
				break;
			}
		}
	}
	return edgeLeng;

}

// Read data of observation points consisting line from input file
void ObservationLine::readObservationLineData( std::ifstream& ifs ){

	for( int i = 0; i < 2; ++i ){
		ifs >> m_endPoints[i].X >> m_endPoints[i].Y >> m_endPoints[i].Z;
	}

	ifs >> m_numLayers;
	if( m_numLayers <= 0 ){
		std::cerr << "Error : Total number of layers is less than 1 : " << m_numLayers << std::endl;
		exit(1);
	}

	m_radius = new double[m_numLayers];
	m_maxEdgeLengthHorizontal = new double[m_numLayers];
	m_maxEdgeLengthVertical = new double[m_numLayers];

	for( int i = 0; i < m_numLayers; ++i ){
#ifdef _VERTICAL_LENGTH_CTRL
		ifs >> m_radius[i] >> m_maxEdgeLengthHorizontal[i] >> m_maxEdgeLengthVertical[i];
#else
		ifs >> m_radius[i] >> m_maxEdgeLengthHorizontal[i];
		m_maxEdgeLengthVertical[i] = m_maxEdgeLengthHorizontal[i];
#endif
		if( i > 0 && ( m_radius[i] < m_radius[i-1] ) ){
			std::cerr << "Error : Inner radius ( " << m_radius[i-1] << " ) is greater than outer radius ( " << m_radius[i] << " ) !! " << std::endl;
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
	}

	std::cout << " First coordinate  : " << m_endPoints[0].X << " " << m_endPoints[0].Y << " " << m_endPoints[0].Z << std::endl;
	std::cout << " Second coordinate : " << m_endPoints[1].X << " " << m_endPoints[1].Y << " " << m_endPoints[1].Z << std::endl;
	std::cout << " Number of layers  : " << m_numLayers << std::endl;
	for( int i = 0; i < m_numLayers; ++i ){
		std::cout << " Radius [km] : " << m_radius[i]
				<< ", Edge length horizontal [km] : " << m_maxEdgeLengthHorizontal[i]
#ifdef _VERTICAL_LENGTH_CTRL
				<< ", Edge length vertical [km] : " << m_maxEdgeLengthVertical[i]
#endif
				<< std::endl;
	}

}

//// Judge the distance between the observation station and the specified point is lower than threshold value
//bool ObservationLine::judgeDistanceObsSiteAndPointLowerThanThreshold( const CommonParameters::XYZ& coord, const double threshold ) const{
//
//	for( int iPoint = 0; iPoint < 2; ++iPoint ){// Loop of the end points
//		const double vecX = coord.X - m_endPoints[iPoint].X;
//		const double vecY = coord.Y - m_endPoints[iPoint].Y;
//		const double vecZ = coord.Z - m_endPoints[iPoint].Z;
//		const double radius = sqrt( pow( vecX, 2.0 ) + pow( vecY, 2.0 ) + pow( vecZ, 2.0 ) );
//		if( radius < threshold){
//			return true;
//		}
//	}
//
//	const double lineLength = sqrt( 
//		pow(m_endPoints[1].X - m_endPoints[0].X, 2) + 
//		pow(m_endPoints[1].Y - m_endPoints[0].Y, 2) + 
//		pow(m_endPoints[1].Z - m_endPoints[0].Z, 2) );
//	const double innerProduct =
//		(coord.X - m_endPoints[0].X) * (m_endPoints[1].X - m_endPoints[0].X) + 
//		(coord.Y - m_endPoints[0].Y) * (m_endPoints[1].Y - m_endPoints[0].Y) + 
//		(coord.Z - m_endPoints[0].Z) * (m_endPoints[1].Z - m_endPoints[0].Z);
//	const double innerProductDiv = innerProduct / lineLength;
//	if( innerProductDiv >= 0.0 && innerProductDiv <= lineLength ){// Within the two end points
//		const CommonParameters::XYZ vectorAlongLine = {
//			innerProductDiv / lineLength * (m_endPoints[1].X - m_endPoints[0].X),
//			innerProductDiv / lineLength * (m_endPoints[1].Y - m_endPoints[0].Y),
//			innerProductDiv / lineLength * (m_endPoints[1].Z - m_endPoints[0].Z) };
//		const CommonParameters::XYZ vectorPerpendicularToLine = {
//			coord.X - m_endPoints[0].X - vectorAlongLine.X,
//			coord.Y - m_endPoints[0].Y - vectorAlongLine.Y,
//			coord.Z - m_endPoints[0].Z - vectorAlongLine.Z };
//		const double distanceFromLine = sqrt(pow(vectorPerpendicularToLine.X, 2) + pow(vectorPerpendicularToLine.Y, 2) + pow(vectorPerpendicularToLine.Z, 2) );
//		if( distanceFromLine < threshold ){
//			return true;
//		}
//	}
//
//	return false;
//
//}
//
