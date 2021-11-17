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
#include "ObservingSiteList.h"
#include <stdlib.h>
#include <string>  
#include <fstream>
#include <iostream>
#include <assert.h>

// Default constructer
ObservingSiteList::ObservingSiteList():
	m_numObservationPoint(0),
	m_numObservationLine(0),
	m_obsPoint(NULL),
	m_obsLine(NULL)
{
}

// Destructer
ObservingSiteList::~ObservingSiteList(){

	if( m_obsPoint != NULL ){
		delete [] m_obsPoint;
		m_obsPoint = NULL;
	}

	if( m_obsLine != NULL ){
		delete [] m_obsLine;
		m_obsLine = NULL;
	}

}

// Read data of observing site list from input file
void ObservingSiteList::readObservingSiteData(){

	std::ifstream ifs( "obs_site.dat", std::ios::in );

	if( ifs.fail() ){
		std::cerr << "Error : File open error : obs_site.dat !!" << std::endl;
		exit(1);
	}

	std::cout << "Read data of observation site list." << std::endl;

	ifs >> m_numObservationPoint;
	if( m_numObservationPoint > 0 ){
		m_obsPoint = new ObservationPoint[m_numObservationPoint];
	}

	for( int i = 0; i < m_numObservationPoint; ++i ){
		m_obsPoint[i].readObservationPointData( ifs );
    }

	ifs >> m_numObservationLine;
	if( m_numObservationLine > 0 ){
		m_obsLine = new ObservationLine[m_numObservationLine];
	}

	for( int i = 0; i < m_numObservationLine; ++i ){
		m_obsLine[i].readObservationLineData( ifs );
    }

	ifs.close();

}

// Calculate maximum horizontal length of specified coordinate
double ObservingSiteList::calcMaximumLengthHorizontal( const CommonParameters::XYZ& coord ) const{

	double val = 1.0e20;

	for( int i = 0; i < m_numObservationPoint; ++i ){
		const double dbuf = m_obsPoint[i].calcMaximumHorizontalLength( coord );
		if( dbuf < val ){
			val = dbuf;
		}
    }

	for( int i = 0; i < m_numObservationLine; ++i ){
		const double dbuf = m_obsLine[i].calcMaximumHorizontalLength( coord );
		if( dbuf < val ){
			val = dbuf;
		}
    }

	return val;

}

// Calculate maximum vertical length of specified coordinate
double ObservingSiteList::calcMaximumLengthVertical( const CommonParameters::XYZ& coord ) const{

	double val = 1.0e20;

	for( int i = 0; i < m_numObservationPoint; ++i ){
		const double dbuf = m_obsPoint[i].calcMaximumVerticalLength( coord );
		if( dbuf < val ){
			val = dbuf;
		}
    }

	for( int i = 0; i < m_numObservationLine; ++i ){
		const double dbuf = m_obsLine[i].calcMaximumVerticalLength( coord );
		if( dbuf < val ){
			val = dbuf;
		}
    }

	return val;

}

//// Judge the distance between the observation stations and the specified points is lower than threshold value
//bool ObservingSiteList::judgeDistanceObsSitesAndPointLowerThanThreshold( const CommonParameters::XYZ& coord, const double threshold ) const{
//
//	for( int i = 0; i < m_numObservationPoint; ++i ){
//		if( m_obsPoint[i].judgeDistanceObsSiteAndPointLowerThanThreshold( coord, threshold ) ){
//			return true;
//		}
//    }
//
//	for( int i = 0; i < m_numObservationLine; ++i ){
//		if( m_obsLine[i].judgeDistanceObsSiteAndPointLowerThanThreshold( coord, threshold ) ){
//			return true;
//		}
//    }
//
//	return false;
//
//}
