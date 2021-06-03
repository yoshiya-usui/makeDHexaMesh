//--------------------------------------------------------------------------
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
