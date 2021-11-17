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
#ifndef DBLDEF_OBSERVING_SITE_LIST
#define DBLDEF_OBSERVING_SITE_LIST

#include "ObservationPoint.h"
#include "ObservationLine.h"
#include "CommonParameters.h"

// Class of the list of the observing site
class ObservingSiteList{

public:
	// Default constructer
	ObservingSiteList();

	// Destructer
	~ObservingSiteList();

	// Read data of observing site list from input file
	void readObservingSiteData();

	// Calculate maximum horizontal length of specified coordinate
	double calcMaximumLengthHorizontal( const CommonParameters::XYZ& coord ) const;

	// Calculate maximum vertical length of specified coordinate
	double calcMaximumLengthVertical( const CommonParameters::XYZ& coord ) const;

	//// Judge the distance between the observation stations and the specified point is lower than threshold value
	//bool judgeDistanceObsSitesAndPointLowerThanThreshold( const CommonParameters::XYZ& coord, const double threshold ) const;

private:

	// Copy constructer
	ObservingSiteList(const ObservingSiteList& rhs);

	// Assignment operator
	ObservingSiteList& operator=(const ObservingSiteList& rhs);

	// Total number of observation point
	int m_numObservationPoint;

	// Total number of observation line
	int m_numObservationLine;

	// Array of observation point
	ObservationPoint* m_obsPoint;

	// Array of observation line
	ObservationLine* m_obsLine;

};

#endif
