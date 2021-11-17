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
#ifndef DBLDEF_OBSERVATION_POINT
#define DBLDEF_OBSERVATION_POINT

#include <fstream>
#include <iostream>
#include "CommonParameters.h"

// Class of observation point
class ObservationPoint{

public:
	// Default Constructer
	ObservationPoint();

	// Destructer
	~ObservationPoint();

	// Read data of observation point from input file
	void readObservationPointData( std::ifstream& ifs );

	// Calculate maximum horizontal length
	double calcMaximumHorizontalLength( const CommonParameters::XYZ& coord ) const;

	// Calculate maximum vertical length
	double calcMaximumVerticalLength( const CommonParameters::XYZ& coord ) const;

	// Calculate factor for intermediate region
	double calcFactorForIntermediateRegion( const CommonParameters::XYZ& coord, const int iSphere ) const;

	//// Judge the distance between the observation station and the specified point is lower than threshold value
	//bool judgeDistanceObsSiteAndPointLowerThanThreshold( const CommonParameters::XYZ& coord, const double threshold ) const;

private:
	enum TypeOfRegion{
		SPHERE = 0,
		CUBOID,
	};

	// Copy constructer
	ObservationPoint(const ObservationPoint& rhs);

	// Assignment operator
	ObservationPoint& operator=(ObservationPoint& rhs);

	// Coordinate of the point
	CommonParameters::XYZ m_pointCoord;

	// Type of region
	int m_typeOfRegion;

	// Total number of spheres
	int m_numRegions;

	// Length for local mesh control
	double* m_length;

	// Maximum horizontal edge length within sphere around the point
	double* m_maxEdgeLengthHorizontal;

	// Maximum vertical edge length within sphere around the point
	double* m_maxEdgeLengthVertical;

	// Oblateness of regions
	double* m_oblateness;

	bool locateInRegion( const CommonParameters::XYZ& coord, const int iRegion ) const;

};

#endif
