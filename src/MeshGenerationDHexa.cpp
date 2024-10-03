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
#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include "Util.h"
#include "CommonParameters.h"
#include "MeshGenerationDHexa.h"

double MeshGenerationDHexa::distanceConversion = 1000;

// Constructer
MeshGenerationDHexa::MeshGenerationDHexa() :
	m_numXInit(0),
	m_numYInit(0),
	m_numZInit(0),
	m_CoordinatesXInit(NULL),
	m_CoordinatesYInit(NULL),
	m_CoordinatesZInit(NULL),
#ifdef _LAYERS
	m_numLayers(-1),
	//m_numResistivityBlocks(0),
	m_elemGroupingZ(NULL),
	m_numElemGroupX(NULL),
	m_numElemGroupY(NULL),
	m_elemGroupingX(NULL),
	m_elemGroupingY(NULL),
	m_numResisivityBlockAccumulated(NULL),
#endif
	m_initialResistivity(100.0),
	m_airResistivity(1.0e10),
	m_seaResistivity(0.25),
	m_seaDepth(-9999.999),
	m_numResisivityAnomalies(0),
	m_anomalyData(NULL),
	m_edgeIndexOfSeaSurface(-1),
	m_edgeIndexOfEarthSurface(-1),
	m_levelLimitParameterCellPartitioning(100),
	m_hasDivisionNumberRead(false),
	m_hasXCoordinatesRead(false),
	m_hasYCoordinatesRead(false),
	m_hasZCoordinatesRead(false),
	m_hasInitialResistivityRead(false),
	m_hasAirResistivityRead(false),
	m_hasSeaResistivityRead(false),
	m_hasSeaDepthRead(false),
	m_hasAnomaliesRead(false),
#ifdef _LAYERS
	m_hasLayersRead(false),
	m_hasLayerRead(NULL),
#endif
	m_includeSea(false),
	m_incorporateTopo(false),
	m_isUsedThresholdsLimitDistanceBetweenLevelBoundaryAndObsSites(false),
	m_partitioningType(HORIZONTAL_ONLY),
	m_minimumSeaDepth(0.1),
	m_topographyData(NULL),
	m_is2DStructureAssumedForMesh(false)
{
	// Abscissas of two point Gauss quadrature
	m_abscissas2Point[0] = -1.0 / sqrt(3.0);
	m_abscissas2Point[1] = 1.0 / sqrt(3.0);

	// Weights of two point Gauss quadrature
	m_weights2Point[0] = 1.0;
	m_weights2Point[1] = 1.0;

	// Calculate integral points and weights of two point Gauss quadrature
	int ip(0);
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < 2; ++k) {
				m_integralPointXi[ip] = m_abscissas2Point[i];
				m_integralPointEta[ip] = m_abscissas2Point[j];
				m_integralPointZeta[ip] = m_abscissas2Point[k];
				m_weights[ip] = m_weights2Point[i] * m_weights2Point[j] * m_weights2Point[k];
				++ip;
			}
		}
	}

	// Array of reference coord xi values for each node
	m_xiAtNode[0] = -1.0;
	m_xiAtNode[1] = 1.0;
	m_xiAtNode[2] = 1.0;
	m_xiAtNode[3] = -1.0;
	m_xiAtNode[4] = -1.0;
	m_xiAtNode[5] = 1.0;
	m_xiAtNode[6] = 1.0;
	m_xiAtNode[7] = -1.0;

	// Array of reference coord eta values for each node
	m_etaAtNode[0] = -1.0;
	m_etaAtNode[1] = -1.0;
	m_etaAtNode[2] = 1.0;
	m_etaAtNode[3] = 1.0;
	m_etaAtNode[4] = -1.0;
	m_etaAtNode[5] = -1.0;
	m_etaAtNode[6] = 1.0;
	m_etaAtNode[7] = 1.0;

	// Array of reference coord zeta values for each node
	m_zetaAtNode[0] = -1.0;
	m_zetaAtNode[1] = -1.0;
	m_zetaAtNode[2] = -1.0;
	m_zetaAtNode[3] = -1.0;
	m_zetaAtNode[4] = 1.0;
	m_zetaAtNode[5] = 1.0;
	m_zetaAtNode[6] = 1.0;
	m_zetaAtNode[7] = 1.0;

}

// Destructer
MeshGenerationDHexa::~MeshGenerationDHexa(){

	if( m_CoordinatesXInit != NULL ){
		delete[] m_CoordinatesXInit;
	}

	if( m_CoordinatesYInit != NULL ){
		delete[] m_CoordinatesYInit;
	}

	if( m_CoordinatesZInit != NULL ){
		delete[] m_CoordinatesZInit;
	}

#ifdef _LAYERS
	if( m_elemGroupingZ != NULL ){
		delete[] m_elemGroupingZ;
	}
	
	if( m_numElemGroupX != NULL ){
		delete[] m_numElemGroupX;
	}

	if( m_numElemGroupY != NULL ){
		delete[] m_numElemGroupY;
	}

	if( m_elemGroupingX != NULL ){
		for( int i = 0 ; i < m_numLayers; ++i ){
			if( m_elemGroupingX[i] != NULL ){
				delete[] m_elemGroupingX[i];
			}
		}
		delete[] m_elemGroupingX;
	}

	if( m_elemGroupingY != NULL ){
		for( int i = 0 ; i < m_numLayers; ++i ){
			if( m_elemGroupingY[i] != NULL ){
				delete[] m_elemGroupingY[i];
			}
		}
		delete[] m_elemGroupingY;
	}

	if( m_numResisivityBlockAccumulated != NULL ){
		delete[] m_numResisivityBlockAccumulated;
	}

	if( m_hasLayerRead != NULL ){
		delete[] m_hasLayerRead;
	}
#endif
	
	if( m_anomalyData != NULL ){
		delete[] m_anomalyData;
	}


	if( m_topographyData != NULL ){
		delete m_topographyData;
	}

	for (std::vector<Ellipsoids*>::iterator itr = m_ellipsoids.begin(); itr != m_ellipsoids.end(); ++itr) {
		delete* itr;
	}

	for (std::vector<Cuboids*>::iterator itr = m_cuboids.begin(); itr != m_cuboids.end(); ++itr) {
		delete* itr;
	}
	
}

void MeshGenerationDHexa::readInputFile(){

	std::cout << "Read Input file." << std::endl;	

	const std::string fileName = "meshgen.inp";

	// Read parameters from input file
	std::ifstream inputFile( fileName.c_str(), std::ios::in );
	if( inputFile.fail() ){
		std::cerr << "File open error : " << fileName << " !!" << std::endl;
		exit(1);
	}

	while(!inputFile.eof()){
		std::string line;
		inputFile >> line;
		if( line == "DIVISION_NUMBERS" ){
			//readDivisionNumbers(inputFile);
			if( m_hasDivisionNumberRead == true ){
				std::cerr << "Division numbers has already been read." << std::endl;	
				exit(1);
			}
			else{
				std::cout << "Read division numbers." << std::endl;	
			}
			inputFile >> m_numXInit;
			inputFile >> m_numYInit;
			inputFile >> m_numZInit;
			std::cout << m_numXInit << " " << m_numYInit << " " << m_numZInit << std::endl;
			m_hasDivisionNumberRead = true;
		}
		else if( line == "X_COORDINATES" ){
			//readXCoordinates(inputFile);
			if( m_hasXCoordinatesRead == true ){
				std::cerr << "X coordinates have already been read." << std::endl;	
				exit(1);
			}
			else{
				std::cout << "Read X coordinates." << std::endl;
			}
			m_CoordinatesXInit = new double[ m_numXInit + 1 ];// X coorinates
			for( int i = 0; i < m_numXInit + 1; ++i ){
				inputFile >> m_CoordinatesXInit[i];
				std::cout << m_CoordinatesXInit[i] << " ";
			}
			for (int i = 1; i < m_numXInit + 1; ++i) {
				if( fabs(m_CoordinatesXInit[i] - m_CoordinatesXInit[i - 1]) < CommonParameters::EPS ){
					std::cerr << "Two edge locations (" << m_CoordinatesXInit[i - 1] <<  " km and " << m_CoordinatesXInit[i] <<  " km) are too close!" << std::endl;
					exit(1);
				}
			}
			std::cout << std::endl;
			m_hasXCoordinatesRead = true;
		}
		else if( line == "Y_COORDINATES" ){
			//readYCoordinates(inputFile);
			if( m_hasYCoordinatesRead == true ){
				std::cerr << "Y coordinates have already been read." << std::endl;	
				exit(1);
			}
			else{
				std::cout << "Read Y coordinates." << std::endl;
			}
			m_CoordinatesYInit = new double[ m_numYInit + 1 ];// Y coorinates
			for( int i = 0; i < m_numYInit + 1; ++i ){
				inputFile >> m_CoordinatesYInit[i];
				std::cout << m_CoordinatesYInit[i] << " " ;
			}
			for (int i = 1; i < m_numYInit + 1; ++i) {
				if (fabs(m_CoordinatesYInit[i] - m_CoordinatesYInit[i - 1]) < CommonParameters::EPS) {
					std::cerr << "Two edge locations (" << m_CoordinatesYInit[i - 1] << " km and " << m_CoordinatesYInit[i] << " km) are too close!" << std::endl;
					exit(1);
				}
			}
			std::cout << std::endl;
			m_hasYCoordinatesRead = true;
		}
		else if( line == "Z_COORDINATES" ){
			//readZCoordinates(inputFile);
			if( m_hasZCoordinatesRead == true ){
				std::cerr << "Z coordinates have already been read." << std::endl;	
				exit(1);
			}
			else{
				std::cout << "Read Z coordinates." << std::endl;
			}
			m_CoordinatesZInit = new double[ m_numZInit + 1 ];// Z coorinates
			for( int i = 0; i < m_numZInit + 1; ++i ){
				inputFile >> m_CoordinatesZInit[i];
				std::cout << m_CoordinatesZInit[i] << " ";
			}
			for (int i = 1; i < m_numZInit + 1; ++i) {
				if (fabs(m_CoordinatesZInit[i] - m_CoordinatesZInit[i - 1]) < CommonParameters::EPS) {
					std::cerr << "Two edge locations (" << m_CoordinatesZInit[i - 1] << " km and " << m_CoordinatesZInit[i] << " km) are too close!" << std::endl;
					exit(1);
				}
			}
			std::cout << std::endl;
			m_hasZCoordinatesRead = true;
		}
#ifdef _LAYERS
		else if( line == "LAYERS" ){
			readLayers(inputFile);
		}
		else if( line == "LAYER" ){
			int layerIDStart(0);
			int layerIDEnd(0);
			inputFile >> layerIDStart;
			inputFile >> layerIDEnd;
			if( layerIDStart == 1  ){
				std::cerr << "Resistivity block division cannot be specified for layer one, in which resistivity value is uniform regardless of the input." << std::endl;
				exit(1);
			}else if( layerIDStart < 2 ){
				std::cerr << "IDs of layers must be larger than 1 . : StartID = " << layerIDStart << ", EndID = " << layerIDEnd << std::endl;
				exit(1);
			}else if( layerIDEnd > m_numLayers ){
				std::cerr << "IDs of layers must less than or equal to total number of layer. StartID = " << layerIDStart << ", EndID = " << layerIDEnd << std::endl;
				exit(1);
			}else{
				readLayer( inputFile, layerIDStart, layerIDEnd );
			}
		}
#endif
		else if( line == "INITIAL_RESISTIVITY" ){
			//readInitialResistivity(inputFile);
			if( m_hasInitialResistivityRead == true ){
				std::cerr << "Initial resistivity has already been read." << std::endl;	
				exit(1);
			}
			else{
				std::cout << "Read initial resistivity." << std::endl;	
			}
			inputFile >> m_initialResistivity;
			if( m_initialResistivity < 0 ){
				std::cerr << "Initial resistivity value is less than zero !" << std::endl;
				exit(1);			
			}
			std::cout << m_initialResistivity << std::endl;
			m_hasInitialResistivityRead = true;
		}
		else if( line == "AIR_RESISTIVITY" ){
			if( m_hasAirResistivityRead == true ){
				std::cerr << "Air resistivity has already been read." << std::endl;	
				exit(1);
			}
			else{
				std::cout << "Read resistivity of the air." << std::endl;	
			}
			inputFile >> m_airResistivity;
			if( m_airResistivity < 0 ){
				std::cerr << "Resistivity value of the air is less than zero !" << std::endl;
				exit(1);			
			}
			std::cout << m_airResistivity << std::endl;
			m_hasAirResistivityRead = true;
		}
		else if( line == "SEA_RESISTIVITY" ){
			if( m_hasSeaResistivityRead == true ){
				std::cerr << "Sea resistivity has already been read." << std::endl;	
				exit(1);
			}
			else{
				std::cout << "Read resistivity of the sea." << std::endl;	
			}
			inputFile >> m_seaResistivity;
			if( m_seaResistivity < 0 ){
				std::cerr << "Resistivity value of the sea is less than zero !" << std::endl;
				exit(1);			
			}
			std::cout << m_seaResistivity << std::endl;
			m_hasSeaResistivityRead = true;
		}
		else if( line == "SEA_DEPTH" ){
			if( m_hasSeaDepthRead == true ){
				std::cerr << "Initial sea depth has already been read." << std::endl;	
				exit(1);
			}
			else{
				std::cout << "Read initial sea depth." << std::endl;	
			}
			inputFile >> m_seaDepth;
			if( m_seaDepth < 0 ){
				std::cerr << "Initial sea depth is less than zero !" << std::endl;
				exit(1);			
			}
			m_includeSea = true;
			std::cout << m_seaDepth << std::endl;
			m_hasSeaDepthRead = true;
		}
		else if( line == "TOPO" ){
			m_incorporateTopo = true;
			std::string inFileName;
			int maxNumPoint(0);
			double maxDist(0.0);
			double eps(0.0);
			inputFile >> inFileName;
			inputFile >> maxNumPoint;
			inputFile >> maxDist;
			inputFile >> eps;
			m_topographyData = new TopographyData( inFileName, maxNumPoint, maxDist, eps );
			std::cout << "Topography/bathymetry is incorporated." << std::endl;
			std::cout << "Maximum number of points used for interpolating of z-coordinates: " << maxNumPoint << std::endl;	
			std::cout << "Search radius for interpolating [km]: " << maxDist << std::endl;	
			std::cout << "Distance used to avoid too small denominator in inverse distance weighting [km]: " << eps << std::endl;	
		}
		else if( line == "ANOMALIES" ){
			if( m_hasAnomaliesRead == true ){
				std::cerr << "Data of resistivity anomalies has already been read." << std::endl;	
				exit(1);
			}
			else{
				std::cout << "Read data of resistivity anomalies." << std::endl;	
			}
			inputFile >> m_numResisivityAnomalies;
			std::cout << m_numResisivityAnomalies << std::endl;
			m_anomalyData = new struct DataOfAnomaly[m_numResisivityAnomalies];
			for( int i = 0; i < m_numResisivityAnomalies; ++i ){
				inputFile >> m_anomalyData[i].XStart; 
				inputFile >> m_anomalyData[i].XEnd; 
				inputFile >> m_anomalyData[i].YStart;
				inputFile >> m_anomalyData[i].YEnd;
				inputFile >> m_anomalyData[i].ZStart;
				inputFile >> m_anomalyData[i].ZEnd;
				inputFile >> m_anomalyData[i].resistivityValue; 
				if( m_anomalyData[i].resistivityValue < 0 ){
					std::cerr << "Resistivity value of anomaly " << i << " is less than zero !" << std::endl;
					exit(1);			
				}
				int ibuf(0);
				inputFile >> ibuf;
				m_anomalyData[i].fixResistivityValue = ( ibuf == 1 ? true : false );
				std::cout << m_anomalyData[i].XStart << " " << m_anomalyData[i].XEnd << " "
						  << m_anomalyData[i].YStart << " " << m_anomalyData[i].YEnd << " "
						  << m_anomalyData[i].ZStart << " " << m_anomalyData[i].ZEnd << " "
						  << m_anomalyData[i].resistivityValue;

				if( m_anomalyData[i].fixResistivityValue ){
					std::cout << " Fix"  << std::endl;
				}else{
					std::cout << " Free"  << std::endl;
				}
			}
			m_hasAnomaliesRead = true;
		}
		else if( line == "LEVEL_LIMIT_PARAM_CELL" ){
			int ibuf(0);
			inputFile >> ibuf;
			m_levelLimitParameterCellPartitioning = ibuf;
		}
		else if( line == "TYPE" ){
			int ibuf(0);
			inputFile >> ibuf;
			m_partitioningType = ibuf;
		}
		else if( line == "THRE_SEA_DEPTH" ){
			inputFile >> m_minimumSeaDepth;
			if( m_minimumSeaDepth <= 0.0 ){
				std::cerr << "The minimum sea depth must be positive !" << std::endl;
				exit(1);			
			}
		}
		else if( line == "2D_STRUCTURE" ){
			m_is2DStructureAssumedForMesh = true;
			std::cout << "2-D structure is assumed for the mesh"  << std::endl;
		}
		else if (line == "CUBOIDS") {
			Cuboids* obj = new Cuboids();
			obj->readParameters(inputFile);
			m_cuboids.push_back(obj);
		}
		else if (line == "ELLIPSOIDS") {
		Ellipsoids* obj = new Ellipsoids();
			obj->readParameters(inputFile);
			m_ellipsoids.push_back(obj);
		}
		else if( line == "END" ){
			std::cout << "End of the data." << std::endl;
			break;
		}
	}
	inputFile.close();

	if (m_includeSea) {
		if (m_seaDepth < CommonParameters::EPS ) {
			std::cerr << "Initial sea depth (" << m_seaDepth << " km) is too shallow !" << std::endl;
			exit(1);
		}
	}

	switch (getPartitioningType()){
		case NO_PARTITIONING:
			std::cout << "No partitioning." << std::endl;	
			break;
		case HORIZONTAL_ONLY:
			std::cout << "Partitioning type 1." << std::endl;	
			break;
		case FULL_PARTITIONING:
			std::cout << "Partitioning type 2." << std::endl;	
			break;
		default:
			std::cerr << "Wrong type of partitioning ! : " << m_partitioningType << std::endl;
			exit(1);			
			break;
	}
	if(m_is2DStructureAssumedForMesh && m_levelLimitParameterCellPartitioning > 0 ){
		std::cout << "Warning: Because 2-D structure is assumed, the upper limit of the level parameter cell partitioning is forced to be zero." << std::endl;	
		m_levelLimitParameterCellPartitioning = 0;
	}
	std::cout << "Upper limit of the level parameter cell partitioning : " << m_levelLimitParameterCellPartitioning << std::endl;	

	if( m_incorporateTopo ){
		std::cout << "Threshold of sea depth : " << m_minimumSeaDepth << " (km)" << std::endl;	
	}

	if( getPartitioningType() != NO_PARTITIONING ){
		m_observingSiteList.readObservingSiteData();
	}

}

// Caluculate initial mesh data
void MeshGenerationDHexa::calcInitialMeshData(){

	std::cout << "Generating initial mesh data." << std::endl;

	const int numNodeX = m_numXInit + 1;
	const int numNodeY = m_numYInit + 1;
	const int numNodeZ = m_numZInit + 1;
	const int numNodeTotalInit = numNodeX * numNodeY * numNodeZ; 
	m_nodeCoordinates.reserve(numNodeTotalInit);
	for( int iz = 0 ; iz < numNodeZ; ++iz ){
		for( int iy = 0 ; iy < numNodeY; ++iy ){
			for( int ix = 0 ; ix < numNodeX; ++ix ){
				// Location of nodes
				CommonParameters::XYZ coord = { m_CoordinatesXInit[ix], m_CoordinatesYInit[iy], m_CoordinatesZInit[iz] };
				m_nodeCoordinates.push_back(coord);
			}
		}
	}

	// Search edge indexes of the sea surface and the Earth's surface
	searchEdgeIndexesOfSeaAndEarthSurface();

	// Caluculate neighbor elements and locations of center points.
	const int numElemTotalInit = m_numXInit * m_numYInit * m_numZInit; 
	//m_nodesOfElements.reserve(numElemTotalInit);
	//m_elementType.reserve(numElemTotalInit);
	//m_parameterCellOfElements.reserve(numElemTotalInit);
	//m_neighborElements.reserve(numElemTotalInit);
	//m_childElements.reserve(numElemTotalInit);
	//m_treeLevelOfElements.reserve(numElemTotalInit);
	m_elemInfo.reserve(numElemTotalInit);
	int*** xyzToElemID = NULL;// Three dimensional array containing element ID
	xyzToElemID = new int**[ m_numXInit ];
	for( int ix = 0 ; ix < m_numXInit; ++ix ){
		xyzToElemID[ix] = new int*[ m_numYInit ];
		for( int iy = 0 ; iy < m_numYInit; ++iy ){
			xyzToElemID[ix][iy] = new int[ m_numZInit ];
			for( int iz = 0 ; iz < m_numZInit; ++iz ){
				xyzToElemID[ix][iy][iz] = -1;// Initialization
			}
		}
	}
	const int numElemOfXYPlane = m_numXInit * m_numYInit;
	const int numNodeOfXYPlane = numNodeX * numNodeY;
	int iElem(0);
	int paramCellID = 1;
	m_parameterCellToResistivity.insert( std::make_pair(0, m_airResistivity) );
	m_parameterCellToFixFlag.insert( std::make_pair(0, 1) );
	if( m_includeSea ){
		m_parameterCellToResistivity.insert( std::make_pair(1, m_seaResistivity) );
		m_parameterCellToFixFlag.insert( std::make_pair(1, 1) );
		paramCellID = 2;
	}

#ifdef _LAYERS
	calcNumResisivityBlockAccumulated();
#endif

	for( int iz = 0 ; iz < m_numZInit; ++iz ){
		for( int iy = 0 ; iy < m_numYInit; ++iy ){
			for( int ix = 0 ; ix < m_numXInit; ++ix ){
				ElementInfo info;
				// IDs of nodes belonging to each element
				const int node0 = iz * numNodeOfXYPlane + iy * numNodeX + ix;
				const int node1 = node0 + 1;
				const int node2 = node0 + numNodeX + 1;
				const int node3 = node0 + numNodeX;
				info.nodes[0] = node0;
				info.nodes[1] = node1;
				info.nodes[2] = node2;
				info.nodes[3] = node3;
				info.nodes[4] = node0 + numNodeOfXYPlane;
				info.nodes[5] = node1 + numNodeOfXYPlane;
				info.nodes[6] = node2 + numNodeOfXYPlane;
				info.nodes[7] = node3 + numNodeOfXYPlane;

				// IDs of neighbor elements
				int neibElem[6];
				neibElem[0] = iElem - 1;
				neibElem[1] = iElem + 1;
				neibElem[2] = iElem - m_numXInit;
				neibElem[3] = iElem + m_numXInit;
				neibElem[4] = iElem - numElemOfXYPlane;
				neibElem[5] = iElem + numElemOfXYPlane;
				// IDs of neighbor elements are setted to be -1 at boundaries.
				if( ix == 0 ) {
					neibElem[0] = -1;
				}
				if( ix == m_numXInit - 1 ) {
					neibElem[1] = -1;
				}
				if( iy == 0 ) {
					neibElem[2] = -1;
				}
				if( iy == m_numYInit - 1 ) {
					neibElem[3] = -1;
				}
				if( iz == 0 ) {
					neibElem[4] = -1;
				}
				if( iz == m_numZInit - 1 ) {
					neibElem[5] = -1;
				}
				for( int iFace = 0; iFace < 6; ++iFace ){
					for( int iNeib = 0; iNeib < 4; ++iNeib ){
						info.neib[iFace][iNeib] = -1;
					}
					info.neib[iFace][0] = neibElem[iFace];
				}

				info.isActive = true;
				info.parent = -1;
				for( int i = 0; i < 8; ++i ){
					info.childs[i] = -1;
				}

				// Tree level of elements
				info.level = 0;
				for( int i = 0; i < 6; ++i ){
					info.levelNeib[i] = 0;
				}
				
				// Type of element
				int elementType(LAND);
				if( m_includeSea ){
					if( iz < m_edgeIndexOfSeaSurface ){
						elementType = AIR;
					}else if( iz < m_edgeIndexOfEarthSurface ){
						elementType = SEA;
					}
				}else{
					if( iz < m_edgeIndexOfEarthSurface ){
						elementType = AIR;
					}
				}
				info.type = elementType;

				// Parameter cell
				if( m_includeSea ){
					switch (elementType){
						case AIR:
							info.parameterCell = 0;
							break;
						case SEA:
							info.parameterCell = 1;
							break;
						default:
							if( m_is2DStructureAssumedForMesh ){
								paramCellID = calcResisivityBlockIDOfInitialMeshFor2DStructure(ix, iy, iz);
							}
#ifdef _LAYERS
							paramCellID = calcResisivityBlockID(ix, iy, iz);
#endif
							info.parameterCell = paramCellID;
							m_parameterCellToResistivity.insert( std::make_pair(paramCellID, m_initialResistivity) );
							m_parameterCellToFixFlag.insert( std::make_pair(paramCellID, 0) );
							++paramCellID;
							break;
					}
				}else{
					switch (elementType){
						case AIR:
							info.parameterCell = 0;
							break;
						default:
							if( m_is2DStructureAssumedForMesh ){
								paramCellID = calcResisivityBlockIDOfInitialMeshFor2DStructure(ix, iy, iz);
							}
#ifdef _LAYERS
							paramCellID = calcResisivityBlockID(ix, iy, iz);
#endif
							info.parameterCell = paramCellID;
							m_parameterCellToResistivity.insert( std::make_pair(paramCellID, m_initialResistivity) );
							m_parameterCellToFixFlag.insert( std::make_pair(paramCellID, 0) );
							++paramCellID;
							break;
					}
				}

				m_elemInfo.push_back(info);

				// Three dimensional array containing element ID
				xyzToElemID[ix][iy][iz] = iElem;

				++iElem;
			}
		}
	}
	
	//--------------------------------------------------
	// Elements and nodes belonging to boundary planes
	//--------------------------------------------------
	// Y-Z Plane ( Minus Side )
	for( int iz = 0; iz < m_numZInit; ++iz ){
		for( int iy = 0; iy < m_numYInit; ++iy ){
			//const ElementAndFace elemAndFace = { xyzToElemID[0][iy][iz], 0 };
			//m_elementAndFaceBoundaryPlanes[YZMinus].push_back(elemAndFace);
			m_elementOnBoundaryPlanes[YZMinus].insert(xyzToElemID[0][iy][iz]);
		}
	}

	// Y-Z Plane ( Plus Side )
	for( int iz = 0; iz < m_numZInit; ++iz ){
		for( int iy = 0; iy < m_numYInit; ++iy ){
			//const ElementAndFace elemAndFace = { xyzToElemID[m_numXInit-1][iy][iz], 1 };
			//m_elementAndFaceBoundaryPlanes[YZPlus].push_back(elemAndFace);
			m_elementOnBoundaryPlanes[YZPlus].insert(xyzToElemID[m_numXInit-1][iy][iz]);
		}
	}

	// Z-X Plane ( Minus Side )
	for( int iz = 0; iz < m_numZInit; ++iz ){
		for( int ix = 0; ix < m_numXInit; ++ix ){
			//const ElementAndFace elemAndFace = { xyzToElemID[ix][0][iz], 2 };
			//m_elementAndFaceBoundaryPlanes[ZXMinus].push_back(elemAndFace);
			m_elementOnBoundaryPlanes[ZXMinus].insert(xyzToElemID[ix][0][iz]);
		}
	}

	// Z-X Plane ( Plus Side )
	for( int iz = 0; iz < m_numZInit; ++iz ){
		for( int ix = 0; ix < m_numXInit; ++ix ){
			//const ElementAndFace elemAndFace = { xyzToElemID[ix][m_numYInit-1][iz], 3 };
			//m_elementAndFaceBoundaryPlanes[ZXPlus].push_back(elemAndFace);
			m_elementOnBoundaryPlanes[ZXPlus].insert(xyzToElemID[ix][m_numYInit-1][iz]);
		}
	}

	// X-Y Plane ( Minus Side ) => Top Boundary
	for( int iy = 0; iy < m_numYInit; ++iy ){
		for( int ix = 0; ix < m_numXInit; ++ix ){
			//const ElementAndFace elemAndFace = { xyzToElemID[ix][iy][0], 4 };
			//m_elementAndFaceBoundaryPlanes[XYMinus].push_back(elemAndFace);
			m_elementOnBoundaryPlanes[XYMinus].insert(xyzToElemID[ix][iy][0]);
		}
	}

	// X-Y Plane ( Plus Side ) => Bottom Boundary
	for( int iy = 0; iy < m_numYInit; ++iy ){
		for( int ix = 0; ix < m_numXInit; ++ix ){
			//const ElementAndFace elemAndFace = { xyzToElemID[ix][iy][m_numZInit-1], 5 };
			//m_elementAndFaceBoundaryPlanes[XYPlus].push_back(elemAndFace);
			m_elementOnBoundaryPlanes[XYPlus].insert(xyzToElemID[ix][iy][m_numZInit-1]);
		}
	}

	// Earth surface
	for( int iy = 0; iy < m_numYInit; ++iy ){
		for( int ix = 0; ix < m_numXInit; ++ix ){
			//const ElementAndFace elemAndFace = { xyzToElemID[ix][iy][m_edgeIndexOfEarthSurface], 4 };
			//m_elementAndFaceEarthSurface.push_back(elemAndFace);
			m_elementOfEarthSurface.insert(xyzToElemID[ix][iy][m_edgeIndexOfEarthSurface]);
		}
	}
	m_elementOfEarthSurfaceWithInactiveElements = m_elementOfEarthSurface;

	// Sea surface
	if( m_includeSea ){
		for( int iy = 0; iy < m_numYInit; ++iy ){
			for( int ix = 0; ix < m_numXInit; ++ix ){
				m_elementOfSeaSurface.insert(xyzToElemID[ix][iy][m_edgeIndexOfSeaSurface]);
			}
		}
		m_elementOfSeaSurfaceWithInactiveElements = m_elementOfSeaSurface;
	}

	for( int ix = 0 ; ix < m_numXInit; ++ix ){
		for( int iy = 0 ; iy < m_numYInit; ++iy ){
			if( xyzToElemID[ix][iy] != NULL ){
				delete[] xyzToElemID[ix][iy];
			}
		}
		if( xyzToElemID[ix] != NULL ){
			delete[] xyzToElemID[ix];
		}
	}
	delete[] xyzToElemID;

}

// Caluculate resistivity distribution for initial mesh
void MeshGenerationDHexa::calcResisivityDistributionForInitialMesh(){

	std::cout << "Caluculate resistivity distribution for initial mesh." << std::endl;

	const double eps = 1.0e-6;
	for( std::vector<ElementInfo>::const_iterator itr = m_elemInfo.begin(); itr != m_elemInfo.end(); ++itr ){
		if( !itr->isActive ){
			// Skip inactive element
			continue;
		}
		const int nodeMinusSide = itr->nodes[0];
		const int nodePlusSide = itr->nodes[6];
		const double blkX1 = m_nodeCoordinates[nodeMinusSide].X;
		const double blkY1 = m_nodeCoordinates[nodeMinusSide].Y;
		const double blkZ1 = m_nodeCoordinates[nodeMinusSide].Z;
		const double blkX2 = m_nodeCoordinates[nodePlusSide].X;
		const double blkY2 = m_nodeCoordinates[nodePlusSide].Y;
		const double blkZ2 = m_nodeCoordinates[nodePlusSide].Z;
		for( int iano = 0; iano < m_numResisivityAnomalies; ++iano ){
			const double x1 = m_anomalyData[iano].XStart;
			const double x2 = m_anomalyData[iano].XEnd;
			const double y1 = m_anomalyData[iano].YStart;
			const double y2 = m_anomalyData[iano].YEnd;
			const double z1 = m_anomalyData[iano].ZStart;
			const double z2 = m_anomalyData[iano].ZEnd;
			if( x1 <= blkX1 + eps && x2 >= blkX2 - eps &&
				y1 <= blkY1 + eps && y2 >= blkY2 - eps &&
				z1 <= blkZ1 + eps && z2 >= blkZ2 - eps ){
				const int paramCell = itr->parameterCell;
				m_parameterCellToResistivity[paramCell] = m_anomalyData[iano].resistivityValue;
				m_parameterCellToFixFlag[paramCell] = m_anomalyData[iano].fixResistivityValue ? 1 : 0;
			}
		}
	}

}

// Reconstruct resistivity distribution to assign different parameter cell to each subsurface element
void MeshGenerationDHexa::reconstructResisivityDistribution( const std::set<int>& elemsSeaToLand, const std::map<int, double>& elemsLandToSea ){

	assert( elemsSeaToLand.empty() || elemsLandToSea.empty() );

	std::cout << "Reconstruct resistivity distribution." << std::endl;

	std::map<int, double> parameterCellToResistivityNew;
	std::map<int, int> parameterCellToFixFlagNew;

	const int levelIimit = m_levelLimitParameterCellPartitioning;
	parameterCellToResistivityNew.insert( std::make_pair(0, m_airResistivity) );
	parameterCellToFixFlagNew.insert( std::make_pair(0, 1) );
	int paramCellIDNew = 1;
	if( m_includeSea ){
		parameterCellToResistivityNew.insert( std::make_pair(paramCellIDNew, m_seaResistivity) );
		parameterCellToFixFlagNew.insert( std::make_pair(paramCellIDNew, 1) );
		++paramCellIDNew;
	}else{
		for( std::map<int, double>::const_iterator itr = elemsLandToSea.begin(); itr != elemsLandToSea.end(); ++itr ){
			ElementInfo& info = m_elemInfo[itr->first];
			if( info.isActive ){
				info.parameterCell = paramCellIDNew;
				parameterCellToResistivityNew.insert( std::make_pair(paramCellIDNew, itr->second) );
				parameterCellToFixFlagNew.insert( std::make_pair(paramCellIDNew, 1) );
				++paramCellIDNew;
			}
		}
	}
	if( m_is2DStructureAssumedForMesh ){
		int elemIndex(0);
		std::vector<int> parameterCellsOrg;
		for( std::vector<ElementInfo>::iterator itr = m_elemInfo.begin(); itr != m_elemInfo.end(); ++itr, ++elemIndex ){
			if( itr->type == AIR || itr->type == SEA ){
				continue;
			}
			if( itr->level > 0 ){
				continue;
			}
			if( !itr->isActive && itr->level != 0 ){
				continue;
			}
			if( allChildrenAreConvertedToSea(elemIndex, elemsLandToSea) ){
				continue;
			}
			const int paramCellIDOrg = itr->parameterCell;
			parameterCellsOrg.push_back(paramCellIDOrg);
		}
		std::sort(parameterCellsOrg.begin(), parameterCellsOrg.end());
		parameterCellsOrg.erase( std::unique( parameterCellsOrg.begin(), parameterCellsOrg.end() ), parameterCellsOrg.end() );
		for( std::vector<int>::const_iterator itrParamCells = parameterCellsOrg.begin(); itrParamCells != parameterCellsOrg.end(); ++itrParamCells ){
			int elemIndex = 0;
			for( std::vector<ElementInfo>::iterator itr = m_elemInfo.begin(); itr != m_elemInfo.end(); ++itr, ++elemIndex ){
				if( itr->type == AIR || itr->type == SEA ){
					continue;
				}
				if( itr->level > 0 ){
					continue;
				}
				if( !itr->isActive && itr->level != 0 ){
					continue;
				}
				if( allChildrenAreConvertedToSea(elemIndex, elemsLandToSea) ){
					continue;
				}
				const int paramCellIDOrg = itr->parameterCell;
				if( paramCellIDOrg != *itrParamCells ){
					continue;
				}
				const double resistivity = m_parameterCellToResistivity[paramCellIDOrg];
				const int fixFlag = m_parameterCellToFixFlag[paramCellIDOrg];
				itr->parameterCell = paramCellIDNew;
				giveSameParamCellIDToChildren(elemIndex, paramCellIDNew, elemsLandToSea);
				parameterCellToResistivityNew.insert( std::make_pair(paramCellIDNew, resistivity) );
				parameterCellToFixFlagNew.insert( std::make_pair(paramCellIDNew, fixFlag) );
			}
			++paramCellIDNew;
		}
	}else{
		int elemIndex(0);
		for( std::vector<ElementInfo>::iterator itr = m_elemInfo.begin(); itr != m_elemInfo.end(); ++itr, ++elemIndex ){
			if( itr->type == AIR || itr->type == SEA ){
				continue;
			}
			if( itr->level > levelIimit ){
				continue;
			}
			if( !itr->isActive && itr->level != levelIimit ){
				continue;
			}
			if( allChildrenAreConvertedToSea(elemIndex, elemsLandToSea) ){
				continue;
			}
			const int paramCellIDOrg = itr->parameterCell;
			const double resistivity = m_parameterCellToResistivity[paramCellIDOrg];
			const int fixFlag = m_parameterCellToFixFlag[paramCellIDOrg];
			itr->parameterCell = paramCellIDNew;
			giveSameParamCellIDToChildren(elemIndex, paramCellIDNew, elemsLandToSea);
			parameterCellToResistivityNew.insert( std::make_pair(paramCellIDNew, resistivity) );
			parameterCellToFixFlagNew.insert( std::make_pair(paramCellIDNew, fixFlag) );
			++paramCellIDNew;
		}
	}
	int elemIndex = 0;
	for( std::vector<ElementInfo>::const_iterator itr = m_elemInfo.begin(); itr != m_elemInfo.end(); ++itr, ++elemIndex ){
		if( itr->type != SEA ){
			continue;
		}
		if( itr->level > levelIimit ){
			continue;
		}
		if( !itr->isActive && itr->level != levelIimit ){
			continue;
		}
		if( changeSeaToLand(elemIndex, elemsSeaToLand, paramCellIDNew) ){
			parameterCellToResistivityNew.insert( std::make_pair(paramCellIDNew, m_initialResistivity) );
			parameterCellToFixFlagNew.insert( std::make_pair(paramCellIDNew, 0) );
			++paramCellIDNew;
		}
	}

	m_parameterCellToResistivity.swap(parameterCellToResistivityNew);
	m_parameterCellToFixFlag.swap(parameterCellToFixFlagNew);

}

// Reconstruct the array of elements just below the Earth surface
void MeshGenerationDHexa::reconstructElementOfEarthSurface(){

	m_elementOfEarthSurface.clear();
	
	int elemIndex(0);
	for( std::vector<ElementInfo>::const_iterator itr = m_elemInfo.begin(); itr != m_elemInfo.end(); ++itr, ++elemIndex ){
		if( !itr->isActive ){
			continue;
		}
		if( itr->type != LAND ){
			continue;
		}
		const int elemIndexUpper = itr->neib[4][0];
		if( elemIndexUpper < 0 ){
			continue;
		}
		const int typeUpper = m_elemInfo[elemIndexUpper].type ;
		if( typeUpper != LAND ){
			m_elementOfEarthSurface.insert(elemIndex);
		}
	}

	m_elementOfEarthSurfaceWithInactiveElements = m_elementOfEarthSurface;

}

#ifdef _LAYERS
// Caluculate accumulated numbers of resisivity blocks 
void MeshGenerationDHexa::calcNumResisivityBlockAccumulated(){
	
	m_numResisivityBlockAccumulated = new int[m_numLayers+1];

	m_numResisivityBlockAccumulated[0] = 0;
	int icount(0);
	for( int iLayer = 0; iLayer < m_numLayers; ++iLayer ){
		icount += m_numElemGroupX[iLayer] * m_numElemGroupY[iLayer];
		m_numResisivityBlockAccumulated[iLayer + 1] = icount;
	}

}

// Caluculate ID of resisivity block
int MeshGenerationDHexa::calcResisivityBlockID( const int ix, const int iy, const int iz ){

	if( iz < m_elemGroupingZ[1] ){ // corresponds to the air
		return 0;
	}

	int iLayer(1);
	for( ; iLayer < m_numLayers + 1; ++iLayer ){
		if( iz >= m_elemGroupingZ[iLayer] && iz < m_elemGroupingZ[iLayer+1] ){
			break;
		}
	}

	int iGx(0);	
	for( ; iGx < m_numElemGroupX[iLayer]; ++iGx ){
		if( ix >= m_elemGroupingX[iLayer][iGx] && ix < m_elemGroupingX[iLayer][iGx+1] ){
			break;
		}
	}

	int iGy(0);	
	for( ; iGy < m_numElemGroupY[iLayer]; ++iGy ){
		if( iy >= m_elemGroupingY[iLayer][iGy] && iy < m_elemGroupingY[iLayer][iGy+1] ){
			break;
		}
	}

	int	blkID = m_numResisivityBlockAccumulated[iLayer] + m_numElemGroupX[iLayer]*iGy + iGx;

	return blkID;

}
#endif

// Caluculate ID of resisivity block of the initial mesh for 2D structure
int MeshGenerationDHexa::calcResisivityBlockIDOfInitialMeshFor2DStructure( const int ix, const int iy, const int iz ) const{

	const int izInEarth = iz - m_edgeIndexOfEarthSurface;;
	int blkID = -1;
	if( m_includeSea ){
		blkID = m_numYInit * izInEarth + iy + 2;
	}else{
		blkID = m_numYInit * izInEarth + iy + 1;
	}

	return blkID;

}

// Output mesh data and resistivity model data
void MeshGenerationDHexa::outputMeshData() const{

	std::cout << "Outputing mesh data." << std::endl;

	// Post processing
	FILE *fp;

	//------------------------
	//--- Output mesh data ---
	//------------------------
	if( (fp = fopen("mesh.dat", "w")) == NULL ) {
		printf("File open error : model_0.dat !! \n");
		exit(1);
	}

	fprintf(fp, "%s\n", "DHEXA" );
	const int numNodeTotal = static_cast<int>( m_nodeCoordinates.size() );
	fprintf(fp, "%10d\n", numNodeTotal );
	for( int iNode = 0; iNode < numNodeTotal; ++iNode ){
		fprintf(fp, "%10d%20f%20f%20f\n",
			iNode,
			m_nodeCoordinates[iNode].X * distanceConversion,			
			m_nodeCoordinates[iNode].Y * distanceConversion,
			m_nodeCoordinates[iNode].Z * distanceConversion );
	}

	//int numElemActive(0);
	//const int numElemTotal = static_cast<int>( m_elemInfo.size() );
	//for( int iElem = 0; iElem < numElemTotal; ++iElem ){	
	//	if( !isActive(iElem) ){
	//		// Skip inactive element
	//		continue;
	//	}
	//	++numElemActive;
	//}
	const int numElemTotal = static_cast<int>( m_elemInfo.size() );
	fprintf(fp, "%10d\n", numElemTotal );
	//int iElemActive(0);
	for( int iElem = 0; iElem < numElemTotal; ++iElem ){	
		assert(isActive(iElem));
		fprintf(fp, "%10d%10d%10d%10d%10d%10d%10d%10d%10d\n",
			iElem,
			m_elemInfo[iElem].nodes[0],
			m_elemInfo[iElem].nodes[1],
			m_elemInfo[iElem].nodes[2],
			m_elemInfo[iElem].nodes[3],
			m_elemInfo[iElem].nodes[4],
			m_elemInfo[iElem].nodes[5],
			m_elemInfo[iElem].nodes[6],
			m_elemInfo[iElem].nodes[7] );
		for( int iNeib = 0; iNeib < 6; ++iNeib ){
			if( m_elemInfo[iElem].neib[iNeib][0] < 0 ){
				// Determined by the first face 
				fprintf(fp, "%s%10d\n", "          ", -1);
			}else{
				int num = 1;
				if( m_elemInfo[iElem].levelNeib[iNeib] > m_elemInfo[iElem].level ){
					if( getPartitioningType() == FULL_PARTITIONING ){
						num = 4;
					}else if( getPartitioningType() == HORIZONTAL_ONLY ){
						if( iNeib == 4 || iNeib == 5 ){
							num = 4;
						}else{
							num = 2;
						}
					}
				}
				fprintf(fp, "%s%10d", "          ", num);
				for( int i = 0; i < num; ++i ){
					fprintf(fp, "%10d", m_elemInfo[iElem].neib[iNeib][i] );
				}
				fprintf(fp, "\n");
			}
		}
		//++iElemActive;
	}

	//--------------------------------------------------
	// Elements and nodes belonging to boundary planes
	//--------------------------------------------------
	for( int iboun = 0; iboun < 6; ++iboun ){
		const std::set<int>& elems = m_elementOnBoundaryPlanes[iboun];
		const int numElems = static_cast<int>( elems.size() );
		fprintf(fp, "%10d\n", numElems );
		int faceID = -1;
		switch (iboun){
			case YZMinus:
				faceID = 0;
				break;
			case YZPlus:
				faceID = 1;
				break;
			case ZXMinus:
				faceID = 2;
				break;
			case ZXPlus:
				faceID = 3;
				break;
			case XYMinus:
				faceID = 4;
				break;
			case XYPlus:
				faceID = 5;
				break;
			default:
				break;
		}
		for( std::set<int>::const_iterator itr = elems.begin(); itr != elems.end(); ++itr ){
			fprintf(fp, "%10d%10d\n", *itr, faceID);
		}
	}

	// Earth surface
	const std::set<int>& elems = m_elementOfEarthSurface;
	const int numElems = static_cast<int>( elems.size() );
	fprintf(fp, "%10d\n", numElems );
	for( std::set<int>::const_iterator itr = elems.begin(); itr != elems.end(); ++itr ){
			fprintf(fp, "%10d%10d\n", *itr, 4);
	}

	fclose(fp);

}

// Output VTK file
void MeshGenerationDHexa::outputVTK( const std::string& fileName ) const{

	std::cout << "Output mesh data to a VTK file : " << fileName << std::endl;

	// VTK file
	std::ofstream vtkFile( fileName.c_str(), std::ios::out );

	vtkFile << "# vtk DataFile Version 2.0" << std::endl;
	vtkFile << "MeshData" << std::endl;
	vtkFile << "ASCII" << std::endl;
	vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
	const int numNodeTotal = static_cast<int>( m_nodeCoordinates.size() );
	vtkFile << "POINTS " << numNodeTotal << " float" << std::endl;
	for( int iNode = 0; iNode < numNodeTotal; ++iNode ){
		vtkFile << m_nodeCoordinates[iNode].X * distanceConversion << " "
			    << m_nodeCoordinates[iNode].Y * distanceConversion << " "
				<< m_nodeCoordinates[iNode].Z * distanceConversion << std::endl;
	}

	//int numElemActive(0);
	//const int numElemTotal = static_cast<int>( m_elemInfo.size() );
	//for( int iElem = 0; iElem < numElemTotal; ++iElem ){	
	//	if( !isActive(iElem) ){
	//		// Skip inactive element
	//		continue;
	//	}
	//	++numElemActive;
	//}
	const int numElemTotal = static_cast<int>( m_elemInfo.size() );
	vtkFile << "CELLS " << numElemTotal << " " << numElemTotal * 9 << std::endl;
	int iElemActive(0);
	for( int iElem = 0; iElem < numElemTotal; ++iElem ){	
		assert(isActive(iElem));
		vtkFile << 8 << " "
				<< m_elemInfo[iElem].nodes[0] << " " 
				<< m_elemInfo[iElem].nodes[1] << " " 
				<< m_elemInfo[iElem].nodes[2] << " " 
				<< m_elemInfo[iElem].nodes[3] << " " 
				<< m_elemInfo[iElem].nodes[4] << " " 
				<< m_elemInfo[iElem].nodes[5] << " " 
				<< m_elemInfo[iElem].nodes[6] << " " 
				<< m_elemInfo[iElem].nodes[7] << std::endl;
		++iElemActive;
	}

	vtkFile << "CELL_TYPES " << numElemTotal << std::endl;
	for( int iElem = 0 ; iElem < numElemTotal; ++iElem ){
		vtkFile << "12" << std::endl;
	}

	vtkFile << "CELL_DATA " << numElemTotal << std::endl;
	vtkFile << "SCALARS BlockID int" <<  std::endl;
	vtkFile << "LOOKUP_TABLE default" <<  std::endl;
	for( int iElem = 0; iElem < numElemTotal; ++iElem ){	
		vtkFile << m_elemInfo[iElem].parameterCell << std::endl;
	}

	vtkFile << "SCALARS Resistivity[Ohm-m] float" <<  std::endl;
	vtkFile << "LOOKUP_TABLE default" <<  std::endl;
	for( int iElem = 0; iElem < numElemTotal; ++iElem ){	
		const int paramCell = m_elemInfo[iElem].parameterCell;
		std::map<int, double>::const_iterator itr = m_parameterCellToResistivity.find(paramCell);
		if( itr == m_parameterCellToResistivity.end() ){
			std::cerr << "Parameter cell " << paramCell << " is not found in m_parameterCellToResistivity." << std::endl;
			exit(1);
		}
		vtkFile << itr->second << std::endl;
	}

	vtkFile << "SCALARS ElemID int" <<  std::endl;
	vtkFile << "LOOKUP_TABLE default" <<  std::endl;
	for( int iElem = 0 ; iElem < numElemTotal; ++iElem ){
		vtkFile << iElem << std::endl;
	}

	vtkFile << "SCALARS Level int" <<  std::endl;
	vtkFile << "LOOKUP_TABLE default" <<  std::endl;
	for( int iElem = 0 ; iElem < numElemTotal; ++iElem ){
		const int level = m_elemInfo[iElem].level;
		vtkFile << level << std::endl;
	}

	vtkFile << "SCALARS Volume float" << std::endl;
	vtkFile << "LOOKUP_TABLE default" << std::endl;
	for (int iElem = 0; iElem < numElemTotal; ++iElem) {
		const double volume = calculateVolume(iElem);
		vtkFile << volume << std::endl;
	}

	vtkFile << "POINT_DATA " << numNodeTotal << std::endl;
	vtkFile << "SCALARS NodeID int" <<  std::endl;
	vtkFile << "LOOKUP_TABLE default" <<  std::endl;
	for( int iNode = 0 ; iNode < numNodeTotal; ++iNode ){
		vtkFile << iNode << std::endl;
	}

	vtkFile.close();

}

// Output resistivity data
void MeshGenerationDHexa::outputResistivityData() const{

	std::cout << "Outputing resistivity data." << std::endl;

	// Post processing
	FILE *fp;
	if( (fp = fopen("resistivity_block_iter0.dat", "w")) == NULL ) {
		printf("File open error : resistivity_block_iter0.dat !! \n");
		exit(1);
	}

	const int numElemTotal = static_cast<int>( m_elemInfo.size() );
	const int numParamCell = static_cast<int>( m_parameterCellToResistivity.size() );
	fprintf(fp, "%10d%10d\n",numElemTotal, numParamCell );
	int iElem(0);
	for( std::vector<ElementInfo>::const_iterator itr = m_elemInfo.begin(); itr != m_elemInfo.end(); ++itr, ++iElem ){
		fprintf(fp, "%10d%10d\n", iElem, itr->parameterCell );
	}

	//fprintf(fp, "%10d%5s%15e%15e%15e%15e%10d\n", 0, "     ", m_airResistivity, 1.0e-20, 1.0e+20, 1.0, 1 );
	//std::map<int, double>::const_iterator itrP = m_parameterCellToResistivity.begin();
	//++itrP;
	//int iBlk(1);
	//for( ;itrP != m_parameterCellToResistivity.end(); ++itrP, ++iBlk ){
	//	assert(itrP->first == iBlk);
	//	fprintf(fp, "%10d%5s%15e%15e%15e%15e%10d\n", iBlk, "     ", itrP->second, 1.0e-20, 1.0e+20, 1.0, 0 );
	//}
	int iBlk(0);
	for( std::map<int, double>::const_iterator itrP = m_parameterCellToResistivity.begin(); 
		itrP != m_parameterCellToResistivity.end(); ++itrP, ++iBlk ){
		assert(itrP->first == iBlk);
		std::map<int, int>::const_iterator itrFlag = m_parameterCellToFixFlag.find(iBlk);
		assert(itrFlag != m_parameterCellToFixFlag.end() );
		const int ifix = itrFlag->second;
		fprintf(fp, "%10d%5s%15e%15e%15e%15e%10d\n", iBlk, "     ", itrP->second, 1.0e-20, 1.0e+20, 1.0, ifix );
	}

	fclose(fp);

}

// Remove inactive elements and merge nodes
void MeshGenerationDHexa::removeInactiveElementsAndMergeNodes(){

	std::cout << "Remove inactive elements and merge nodes" << std::endl;

	//std::map<int,int> nodeIDOrgToMerged;
	//std::vector<CommonParameters::XYZ> nodeCoordinatesNew;
	//int nodeIDOrg(0);
	//int nodeIDNew(0);
	//for( std::vector<CommonParameters::XYZ>::const_iterator itr = m_nodeCoordinates.begin(); itr != m_nodeCoordinates.end(); ++itr, ++nodeIDOrg ){
	//	bool newCoord(true);
	//	int counter(0);
	//	for( std::vector<CommonParameters::XYZ>::iterator itrNew = nodeCoordinatesNew.begin(); itrNew != nodeCoordinatesNew.end(); ++itrNew, ++counter ){
	//		if( Util::isSameLocation( *itr, *itrNew ) ){
	//			newCoord = false;
	//			nodeIDOrgToMerged.insert( std::make_pair(nodeIDOrg, counter) );
	//			break;
	//		}
	//	}
	//	if( newCoord ){
	//		nodeCoordinatesNew.push_back(*itr);
	//		nodeIDOrgToMerged.insert( std::make_pair(nodeIDOrg, nodeIDNew) );
	//		++nodeIDNew;
	//	}
	//}
	//m_nodeCoordinates.swap(nodeCoordinatesNew);

	//std::vector<ElementInfo> elemInfoNew;
	//std::map<int,int> elemIDTotalToActive;

	//int elemIDActive(0);
	//int elemIDTotal(0);
	//for( std::vector<ElementInfo>::const_iterator itr = m_elemInfo.begin(); itr != m_elemInfo.end(); ++itr, ++elemIDTotal ){
	//	if( !itr->isActive ){
	//		// Skip inactive element
	//		continue;
	//	}
	//	elemIDTotalToActive.insert( std::make_pair(elemIDTotal, elemIDActive) );
	//	elemInfoNew.push_back(*itr);
	//	++elemIDActive;
	//}

	//for( std::vector<ElementInfo>::iterator itr = elemInfoNew.begin(); itr != elemInfoNew.end(); ++itr ){
	//	for( int iNode = 0; iNode < 8; ++iNode ){
	//		const int nodeIDOrg = itr->nodes[iNode];
	//		itr->nodes[iNode] = nodeIDOrgToMerged[nodeIDOrg];
	//	}
	//	for( int iNeib = 0; iNeib < 6; ++iNeib ){
	//		for( int i = 0; i < 4; ++i ){
	//			const int elemNeib = itr->neib[iNeib][i];
	//			if( elemNeib < 0 ){
	//				continue;
	//			}
	//			itr->neib[iNeib][i] = elemIDTotalToActive[elemNeib];
	//		}
	//	}
	//}

	//m_elemInfo.swap(elemInfoNew);

	//for( int iboun = 0; iboun < 6; ++iboun ){
	//	const std::set<int>& elems = m_elementOnBoundaryPlanes[iboun];
	//	std::set<int> elemsNew;
	//	for( std::set<int>::const_iterator itr = elems.begin(); itr != elems.end(); ++itr ){
	//		assert(*itr >= 0);
	//		elemsNew.insert( elemIDTotalToActive[*itr] );
	//	}
	//	m_elementOnBoundaryPlanes[iboun].swap(elemsNew);
	//}

	//// Earth surface
	//const std::set<int>& elems = m_elementOfEarthSurface;
	//std::set<int> elemsNew;
	//for( std::set<int>::const_iterator itr = elems.begin(); itr != elems.end(); ++itr ){
	//	assert(*itr >= 0);
	//	elemsNew.insert( elemIDTotalToActive[*itr] );
	//}
	//m_elementOfEarthSurface.swap(elemsNew);

#ifdef _MERGE
	mergeNodes();
#endif
	removeInactiveElements();

	std::cout << "Number of nodes : " << m_nodeCoordinates.size() << ", Number of elements : " << m_elemInfo.size() << std::endl;

}

// Merge nodes
void MeshGenerationDHexa::mergeNodes(){

	std::cout << "Merge nodes." << std::endl;

	const int numNodeTotalOrg = static_cast<int>( m_nodeCoordinates.size() );

	// Divided to 8 regions and a boundary region
	//std::vector<CommonParameters::XYZ> nodeCoordinatesOrg[8];
	//for( int i = 0; i < 8; ++i ){
	//	nodeCoordinatesOrg[i].reserve(numNodeTotalOrg/8);
	//}
	//std::vector<CommonParameters::XYZ> nodeCoordinatesOrgBoundary;
	//nodeCoordinatesOrgBoundary.reserve(numNodeTotalOrg/8);
	std::map<int, CommonParameters::XYZ> nodeCoordinatesOrg[8];
	std::map<int, CommonParameters::XYZ> nodeCoordinatesOrgBoundary;
	const double eps = 1.0e-6;
	int nodeIDOrg(0);
	for( std::vector<CommonParameters::XYZ>::const_iterator itr = m_nodeCoordinates.begin(); itr != m_nodeCoordinates.end(); ++itr, ++nodeIDOrg ){
		int iRegion(-1);
		if( fabs(itr->X) < eps || fabs(itr->Y) < eps || fabs(itr->Z) < eps ){
			// Boundary reegion
			//nodeCoordinatesOrgBoundary.push_back(*itr);
			nodeCoordinatesOrgBoundary.insert(std::make_pair(nodeIDOrg, *itr));
			continue;
		}
		if( itr->X < 0.0 && itr->Y < 0.0 && itr->Z < 0.0 ){
			iRegion = 0;
		}else if( itr->X > 0.0 && itr->Y < 0.0 && itr->Z < 0.0 ){
			iRegion = 1;
		}else if( itr->X > 0.0 && itr->Y > 0.0 && itr->Z < 0.0 ){
			iRegion = 2;
		}else if( itr->X < 0.0 && itr->Y > 0.0 && itr->Z < 0.0 ){
			iRegion = 3;
		}else if( itr->X < 0.0 && itr->Y < 0.0 && itr->Z > 0.0 ){
			iRegion = 4;
		}else if( itr->X > 0.0 && itr->Y < 0.0 && itr->Z > 0.0 ){
			iRegion = 5;
		}else if( itr->X > 0.0 && itr->Y > 0.0 && itr->Z > 0.0 ){
			iRegion = 6;
		}else if( itr->X < 0.0 && itr->Y > 0.0 && itr->Z > 0.0 ){
			iRegion = 7;
		}
		assert(iRegion >= 0);
		//nodeCoordinatesOrg[iRegion].push_back(*itr);
		nodeCoordinatesOrg[iRegion].insert(std::make_pair(nodeIDOrg, *itr));
	}

	std::map<int,int> nodeIDOrgToNew;
	std::vector<CommonParameters::XYZ> nodeCoordinatesNew[8];
	int nodeIDNew(0);
	int counterOffset(0);
	for( int iRegion = 0; iRegion < 8; ++iRegion ){// Loop of regions except boundary region
		const std::map<int, CommonParameters::XYZ>& nodesOrg = nodeCoordinatesOrg[iRegion];
		std::vector<CommonParameters::XYZ>& nodesNew = nodeCoordinatesNew[iRegion];
		nodesNew.reserve(nodesOrg.size());
		for( std::map<int, CommonParameters::XYZ>::const_iterator itrOrg = nodesOrg.begin(); itrOrg != nodesOrg.end(); ++itrOrg ){
			bool newCoord(true);
			int counter = counterOffset;
			for( std::vector<CommonParameters::XYZ>::iterator itrNew = nodesNew.begin(); itrNew != nodesNew.end(); ++itrNew, ++counter ){
				if( Util::isSameLocation(itrOrg->second, *itrNew) ){
					newCoord = false;
					nodeIDOrgToNew.insert( std::make_pair(itrOrg->first, counter) );
					break;
				}
			}
			if( newCoord ){
				nodesNew.push_back(itrOrg->second);
				nodeIDOrgToNew.insert( std::make_pair(itrOrg->first, nodeIDNew) );
				++nodeIDNew;
			}
		}
		counterOffset = nodeIDNew;
	}

	std::vector<CommonParameters::XYZ> nodeCoordinatesNewMerged;
	nodeCoordinatesNewMerged.reserve(numNodeTotalOrg);
	for( int iRegion = 0; iRegion < 8; ++iRegion ){
		nodeCoordinatesNewMerged.insert(nodeCoordinatesNewMerged.end(), nodeCoordinatesNew[iRegion].begin(), nodeCoordinatesNew[iRegion].end());
	}

	for( std::map<int, CommonParameters::XYZ>::const_iterator itrOrg = nodeCoordinatesOrgBoundary.begin(); 
		itrOrg != nodeCoordinatesOrgBoundary.end(); ++itrOrg ){
		bool newCoord(true);
		int counter = 0;
		for( std::vector<CommonParameters::XYZ>::iterator itrNew = nodeCoordinatesNewMerged.begin();
			itrNew != nodeCoordinatesNewMerged.end(); ++itrNew, ++counter ){
			if( Util::isSameLocation(itrOrg->second, *itrNew) ){
				newCoord = false;
				nodeIDOrgToNew.insert( std::make_pair(itrOrg->first, counter) );
				break;
			}
		}
		if( newCoord ){
			nodeCoordinatesNewMerged.push_back(itrOrg->second);
			nodeIDOrgToNew.insert( std::make_pair(itrOrg->first, nodeIDNew) );
			++nodeIDNew;
		}
	}

	m_nodeCoordinates.swap(nodeCoordinatesNewMerged);

	for( std::vector<ElementInfo>::iterator itr = m_elemInfo.begin(); itr != m_elemInfo.end(); ++itr ){
		for( int iNode = 0; iNode < 8; ++iNode ){
			const int nodeIDOrg = itr->nodes[iNode];
			itr->nodes[iNode] = nodeIDOrgToNew[nodeIDOrg];
		}
	}

}

// Remove inactive elements
void MeshGenerationDHexa::removeInactiveElements(){
	
	std::cout << "Remove inactive elements." << std::endl;

	const int numElemTotal = static_cast<int>(m_elemInfo.size());
	std::vector<ElementInfo> elemInfoNew;
	elemInfoNew.reserve(numElemTotal);
	std::map<int,int> elemIDTotalToActive;

	int elemIDActive(0);
	int elemIDTotal(0);
	for( std::vector<ElementInfo>::const_iterator itr = m_elemInfo.begin(); itr != m_elemInfo.end(); ++itr, ++elemIDTotal ){
		if( !itr->isActive ){
			// Skip inactive element
			continue;
		}
		elemIDTotalToActive.insert( std::make_pair(elemIDTotal, elemIDActive) );
		elemInfoNew.push_back(*itr);
		++elemIDActive;
	}

	for( std::vector<ElementInfo>::iterator itr = elemInfoNew.begin(); itr != elemInfoNew.end(); ++itr ){
		for( int iNeib = 0; iNeib < 6; ++iNeib ){
			for( int i = 0; i < 4; ++i ){
				const int elemNeib = itr->neib[iNeib][i];
				if( elemNeib < 0 ){
					continue;
				}
				itr->neib[iNeib][i] = elemIDTotalToActive[elemNeib];
			}
		}
	}

	m_elemInfo.swap(elemInfoNew);

	for( int iboun = 0; iboun < 6; ++iboun ){
		const std::set<int>& elems = m_elementOnBoundaryPlanes[iboun];
		std::set<int> elemsNew;
		for( std::set<int>::const_iterator itr = elems.begin(); itr != elems.end(); ++itr ){
			assert(*itr >= 0);
			elemsNew.insert( elemIDTotalToActive[*itr] );
		}
		m_elementOnBoundaryPlanes[iboun].swap(elemsNew);
	}

	{// Earth surface
		const std::set<int>& elems = m_elementOfEarthSurface;
		std::set<int> elemsNew;
		for( std::set<int>::const_iterator itr = elems.begin(); itr != elems.end(); ++itr ){
			assert(*itr >= 0);
			elemsNew.insert( elemIDTotalToActive[*itr] );
		}
		m_elementOfEarthSurface.swap(elemsNew);	
	}

	// Sea surface
	if( m_includeSea ){
		const std::set<int>& elems = m_elementOfSeaSurface;
		std::set<int> elemsNew;
		for( std::set<int>::const_iterator itr = elems.begin(); itr != elems.end(); ++itr ){
			assert(*itr >= 0);
			elemsNew.insert( elemIDTotalToActive[*itr] );
		}
		m_elementOfSeaSurface.swap(elemsNew);
	}

}

// Type of partitioning
int MeshGenerationDHexa::getPartitioningType() const{
	return m_partitioningType;
}

// Check whether each parameter cell contains at least one active element
void MeshGenerationDHexa::checkWhetherEachParameterCellContainsAtLeastOneActiveElement() const{

	const int numParamCells = static_cast<int>( m_parameterCellToResistivity.size() );
	int* paramCellToNumElems = new int [numParamCells];
	for( int i = 0; i < numParamCells; ++i ){
		paramCellToNumElems[i] = 0;
	}

	for( std::vector<ElementInfo>::const_iterator itr = m_elemInfo.begin(); itr != m_elemInfo.end(); ++itr ){
		if( !itr->isActive ){
			continue;
		}
		if( itr->parameterCell < 0 ){
			std::cerr << "Parameter cell index is negative : " << itr->parameterCell  << std::endl;
			exit(1);
		}
		if( itr->parameterCell >= numParamCells ){
			std::cerr << "Parameter cell index is larger than or equal to the number of parameter cells : " << itr->parameterCell  << std::endl;
			exit(1);
		}
		++paramCellToNumElems[itr->parameterCell];
	}

	for( int i = 0; i < numParamCells; ++i ){
		if( paramCellToNumElems[i] == 0 ){
			std::cerr << "Parameter cell " << i << " has no element !!" << std::endl;
			exit(1);
		}
#ifdef _DEBUG_WRITE
		std::cout << i << " " << paramCellToNumElems[i] << std::endl;
#endif
	}

	delete [] paramCellToNumElems;

}

// Include topography
void MeshGenerationDHexa::includeTopography( std::set<int>& elemsSeaToLand ){

	if(!m_incorporateTopo){
		std::map<int, double> dummy;
		reconstructResisivityDistribution(elemsSeaToLand, dummy);
		return;
	}

	std::cout << "Topography/bathymetry is incorporated." << std::endl;

	const double xMin = m_CoordinatesXInit[0];
	const double xMax = m_CoordinatesXInit[m_numXInit];
	const double yMin = m_CoordinatesYInit[0];
	const double yMax = m_CoordinatesYInit[m_numYInit];
	m_topographyData->readTopographyData( xMin, xMax, yMin, yMax );

	const std::vector<CommonParameters::XYZ> nodeCoordinatesOrg = m_nodeCoordinates;

	std::multimap<int, int> nodeToElem;
	int iElem(0);
	int maxLevel(0);
	for( std::vector<ElementInfo>::const_iterator itr = m_elemInfo.begin(); itr != m_elemInfo.end(); ++itr, ++iElem ){
		for( int i = 0; i < 8; ++i ){
			const int nodeID = itr->nodes[i];
			nodeToElem.insert( std::make_pair(nodeID, iElem) );
		}
		if( itr->level > maxLevel ){
			maxLevel = itr->level;
		}
	}

	for( int iLevel = 0; iLevel <= maxLevel; ++iLevel ){
		includeTopographyAux( iLevel, maxLevel, nodeCoordinatesOrg, nodeToElem );
	}

	std::map<int, double> elemsLandToSea;
	std::map<int, double> elemEarthSurfToSeaDepth;

	if( m_includeSea ){
		selectElementsToBeChangedToLandFromSea( elemsSeaToLand );
	}else{
		selectElementsToBeChangedToSeaFromLand( elemsLandToSea, elemEarthSurfToSeaDepth );
	}

	reconstructResisivityDistribution(elemsSeaToLand, elemsLandToSea);

	if( !elemsSeaToLand.empty() ){
		reconstructElementOfEarthSurface();

		const std::vector<CommonParameters::XYZ> nodeCoordinatesOrg2 = m_nodeCoordinates;
		
		std::cout << "Topography is incorporated for land area." << std::endl;
		for( int iLevel = 0; iLevel <= maxLevel; ++iLevel ){
			includeTopographyAux2( iLevel, maxLevel, nodeCoordinatesOrg2, nodeToElem );
		}
	}

	outputEarthSurfaceDepthByVtk(elemEarthSurfToSeaDepth);

}

//// Include topography for land area changed from the sea
//void MeshGenerationDHexa::includeTopographyForLandAreaChangedFromSea( std::set<int>& surfElemsLand, std::set<int>& surfElemSea ){
//
//	if( !m_incorporateTopo || !m_includeSea ){
//		return;
//	}
//
//	std::cout << "Topography is incorporated for land area changed from the sea." << std::endl;
//
//	const std::vector<CommonParameters::XYZ> nodeCoordinatesOrg = m_nodeCoordinates;
//
//	std::multimap<int, int> nodeToElem;
//	int iElem(0);
//	int maxLevel(0);
//	for( std::vector<ElementInfo>::const_iterator itr = m_elemInfo.begin(); itr != m_elemInfo.end(); ++itr, ++iElem ){
//		for( int i = 0; i < 8; ++i ){
//			const int nodeID = itr->nodes[i];
//			nodeToElem.insert( std::make_pair(nodeID, iElem) );
//		}
//		if( itr->level > maxLevel ){
//			maxLevel = itr->level;
//		}
//	}
//
//	for( int iLevel = 0; iLevel <= maxLevel; ++iLevel ){
//		includeTopographyForLandAreaChangedFromSeaAux( iLevel, nodeCoordinatesOrg, surfElemsLand, surfElemSea, nodeToElem );
//	}
//}

// Get flag specifing whether inputted element is active 
bool MeshGenerationDHexa::isActive( const int elemID ) const{
	
	return m_elemInfo[elemID].isActive;

}

// Add children to candidates of the elements to be changed to the land
void MeshGenerationDHexa::addChildrenToLandCandidates( const int elemIndex, std::set<int>& elemsSeaToLand ){

	if( elemIndex < 0 ){
		return;
	}
	ElementInfo& elemInfo = m_elemInfo[elemIndex];
	if( elemInfo.type == SEA && elemsSeaToLand.find(elemIndex) == elemsSeaToLand.end() ){
		elemsSeaToLand.insert( elemIndex );
	}
	for( int i = 0; i < 8; ++i ){
		if( elemInfo.childs[i] < 0 ){
			continue;
		}
		addChildrenToLandCandidates(elemInfo.childs[i], elemsSeaToLand);
	}
	
}

// Get whether all children are converted to the sea
bool MeshGenerationDHexa::allChildrenAreConvertedToSea( const int elemIndex, const std::map<int, double>& elemsLandToSea ) const{

	if( elemsLandToSea.find(elemIndex) == elemsLandToSea.end() ){
		// This element is not converted to the sea
		return false;
	}
	const ElementInfo& elemInfo = m_elemInfo[elemIndex];
	for( int i = 0; i < 8; ++i ){
		if( elemInfo.childs[i] < 0 ){
			continue;
		}
		if( !allChildrenAreConvertedToSea(elemInfo.childs[i], elemsLandToSea) ){
			return false;
		}
	}
	return true;

}

// Search edge indexes of the sea surface and the Earth's surface
void MeshGenerationDHexa::searchEdgeIndexesOfSeaAndEarthSurface(){

	const double eps = 1.0e-6;
	if( m_includeSea ){
		bool foundSeaSurface(false);
		bool foundEarthSurface(false);
		for( int iz = 0 ; iz < m_numZInit + 1; ++iz ){
			const double z = m_CoordinatesZInit[iz];
			if( fabs(z) < eps ){
				m_edgeIndexOfSeaSurface = iz;
				foundSeaSurface = true;
			}
			if( fabs(z - m_seaDepth) < eps ){
				m_edgeIndexOfEarthSurface = iz;
				foundEarthSurface = true;
			}
		}
		if( !foundSeaSurface ){
			std::cerr << "Edge corresponding the sea surface (z = 0 km) is not found." << std::endl;
			exit(1);
		}
		if( !foundEarthSurface ){
			std::cerr << "Edge corresponding the sea floor (z = " << m_seaDepth << " km) is not found." << std::endl;
			exit(1);
		}
	}else{
		bool found(false);
		for( int iz = 0 ; iz < m_numZInit + 1; ++iz ){
			const double z = m_CoordinatesZInit[iz];
			if( fabs(z) < eps ){
				m_edgeIndexOfEarthSurface = iz;
				found = true;
			}
		}
		if( !found ){
			std::cerr << "Edge corresponding the sea floor (z = 0 km) is not found." << std::endl;
			exit(1);
		}
	}

}

// Calculate the maximum horizontal edge length at each node;
void MeshGenerationDHexa::calcMaxEdgeLengthHorizontalAtNode( std::multimap<int, std::pair<int, double> >& nodeToElemsAndEdgeLenth ) const{

	const int nodeToEdge[8][2] = {
		{ 0, 4 }, // 0
		{ 0, 6 }, // 1
		{ 1, 6 }, // 2
		{ 1, 4 }, // 3
		{ 2, 5 }, // 4
		{ 2, 7 }, // 5
		{ 3, 7 }, // 6
		{ 3, 5 }, // 7
	};

	const int numElemTotal = static_cast<int>( m_elemInfo.size() );
	int iElemActive(0);
	for( int iElem = 0; iElem < numElemTotal; ++iElem ){	
		if( !isActive(iElem) ){
			// Skip inactive element
			continue;
		}
		CommonParameters::XYZ coords[8];
		for( int i = 0; i < 8; ++i ){
			const int nodeID = m_elemInfo[iElem].nodes[i];
			coords[i].X = m_nodeCoordinates[nodeID].X;
			coords[i].Y = m_nodeCoordinates[nodeID].Y;
			coords[i].Z = m_nodeCoordinates[nodeID].Z;
		}
		const double edgeLength[12] = {
			Util::calcLength( coords[0], coords[1] ),
			Util::calcLength( coords[2], coords[3] ),
			Util::calcLength( coords[4], coords[5] ),
			Util::calcLength( coords[6], coords[7] ),
			Util::calcLength( coords[0], coords[3] ),
			Util::calcLength( coords[4], coords[7] ),
			Util::calcLength( coords[1], coords[2] ),
			Util::calcLength( coords[5], coords[6] ),
			Util::calcLength( coords[0], coords[4] ),
			Util::calcLength( coords[1], coords[5] ),
			Util::calcLength( coords[3], coords[7] ),
			Util::calcLength( coords[2], coords[6] )
		};
		for( int iNode = 0; iNode < 8; ++iNode ){
			const int nodeID = m_elemInfo[iElem].nodes[iNode];
			double lengthMax = 0.0;
			for( int i = 0; i < 2; ++i ){
				const double length = edgeLength[ nodeToEdge[iNode][i] ];
				if( length > lengthMax ){
					lengthMax = length;
				}
			}
			std::pair<int, double> elemAndEdge = std::make_pair(iElem, lengthMax);
			nodeToElemsAndEdgeLenth.insert( std::make_pair( nodeID, elemAndEdge ) );
		}
	}

}

// Calculate the maximum vertical edge length at each node;
void MeshGenerationDHexa::calcMaxEdgeLengthVerticalAtNode( std::multimap<int, std::pair<int, double> >& nodeToElemsAndEdgeLenth ) const{

	const int nodeToEdge[8] = {
		8, // 0
		9, // 1
		11,// 2
		10,// 3
		8, // 4
		9, // 5
		11,// 6
		10,// 7
	};

	const int numElemTotal = static_cast<int>( m_elemInfo.size() );
	int iElemActive(0);
	for( int iElem = 0; iElem < numElemTotal; ++iElem ){	
		if( !isActive(iElem) ){
			// Skip inactive element
			continue;
		}
		CommonParameters::XYZ coords[8];
		for( int i = 0; i < 8; ++i ){
			const int nodeID = m_elemInfo[iElem].nodes[i];
			coords[i].X = m_nodeCoordinates[nodeID].X;
			coords[i].Y = m_nodeCoordinates[nodeID].Y;
			coords[i].Z = m_nodeCoordinates[nodeID].Z;
		}
		const double edgeLength[12] = {
			Util::calcLength( coords[0], coords[1] ),
			Util::calcLength( coords[2], coords[3] ),
			Util::calcLength( coords[4], coords[5] ),
			Util::calcLength( coords[6], coords[7] ),
			Util::calcLength( coords[0], coords[3] ),
			Util::calcLength( coords[4], coords[7] ),
			Util::calcLength( coords[1], coords[2] ),
			Util::calcLength( coords[5], coords[6] ),
			Util::calcLength( coords[0], coords[4] ),
			Util::calcLength( coords[1], coords[5] ),
			Util::calcLength( coords[3], coords[7] ),
			Util::calcLength( coords[2], coords[6] )
		};
		for( int iNode = 0; iNode < 8; ++iNode ){
			const int nodeID = m_elemInfo[iElem].nodes[iNode];
			double lengthMax = edgeLength[ nodeToEdge[iNode] ];
			std::pair<int, double> elemAndEdge = std::make_pair(iElem, lengthMax);
			nodeToElemsAndEdgeLenth.insert( std::make_pair( nodeID, elemAndEdge ) );
		}
	}


}

// Partition mesh
void MeshGenerationDHexa::partitionMesh(){

	if( getPartitioningType() == NO_PARTITIONING ){
		return;
	}

	const int maxIteration = 100;
	int iter(0);
	while (true){
		if( iter > maxIteration ){
			std::cerr << "Rearch max iteration number : " << maxIteration << std::endl;
			exit(1);
		}
		++iter;
		// Partition mesh one time
		const int numElemParitioned = partitionMeshOneTime();
		std::cout << "Number of elements partitioned : " << numElemParitioned
			<< ", Number of nodes : " << m_nodeCoordinates.size() 
			<< ", Number of elements : " << m_elemInfo.size() << std::endl;

		if( numElemParitioned < 1 ){
			std::cout << "Finishing of partition" << std::endl;
			break;
		}
	}

}

// Partition mesh one time
int MeshGenerationDHexa::partitionMeshOneTime(){

	typedef std::multimap<int, std::pair<int, double> > MMIID;
	std::set<int> elementsToBePartitioned;

	// Calculate the horizontal maximum edge length at each node
	std::multimap<int, std::pair<int, double> > nodeToElemsAndEdgeLenthHorizontal;
	calcMaxEdgeLengthHorizontalAtNode(nodeToElemsAndEdgeLenthHorizontal);
	for( MMIID::const_iterator itr = nodeToElemsAndEdgeLenthHorizontal.begin(); itr != nodeToElemsAndEdgeLenthHorizontal.end(); ++itr ){
		const int nodeIndex = itr->first;
		const int elemIndex = itr->second.first;
		const double maxLength = itr->second.second;
		const CommonParameters::XYZ coord = m_nodeCoordinates[nodeIndex];
		double upperLimit = m_observingSiteList.calcMaximumLengthHorizontal(coord);
		for (std::vector<Ellipsoids*>::const_iterator itr = m_ellipsoids.begin(); itr != m_ellipsoids.end(); ++itr) {
			const double dbuf = (*itr)->calcMaximumEdgeLengthHorizontal(coord);
			if (dbuf < upperLimit) {
				upperLimit = dbuf;
			}
		}
		for (std::vector<Cuboids*>::const_iterator itr = m_cuboids.begin(); itr != m_cuboids.end(); ++itr) {
			const double dbuf = (*itr)->calcMaximumEdgeLengthHorizontal(coord);
			if (dbuf < upperLimit) {
				upperLimit = dbuf;
			}
		}
		if( maxLength > upperLimit ){
			elementsToBePartitioned.insert(elemIndex);
		}
	}

	// Calculate the vertical maximum edge length at each node
	if( getPartitioningType() == FULL_PARTITIONING ){
		std::multimap<int, std::pair<int, double> > nodeToElemsAndEdgeLenthVertical;
		calcMaxEdgeLengthVerticalAtNode(nodeToElemsAndEdgeLenthVertical);
		for( MMIID::const_iterator itr = nodeToElemsAndEdgeLenthVertical.begin(); itr != nodeToElemsAndEdgeLenthVertical.end(); ++itr ){
			const int nodeIndex = itr->first;
			const int elemIndex = itr->second.first;
			const double maxLength = itr->second.second;
			const CommonParameters::XYZ coord = m_nodeCoordinates[nodeIndex];
			double upperLimit = m_observingSiteList.calcMaximumLengthVertical(coord);
			for (std::vector<Ellipsoids*>::const_iterator itr = m_ellipsoids.begin(); itr != m_ellipsoids.end(); ++itr) {
				const double dbuf = (*itr)->calcMaximumEdgeLengthVertical(coord);
				if (dbuf < upperLimit) {
					upperLimit = dbuf;
			}
		}
			for (std::vector<Cuboids*>::const_iterator itr = m_cuboids.begin(); itr != m_cuboids.end(); ++itr) {
				const double dbuf = (*itr)->calcMaximumEdgeLengthVertical(coord);
				if (dbuf < upperLimit) {
					upperLimit = dbuf;
				}
			}
			if( maxLength > upperLimit ){
				elementsToBePartitioned.insert(elemIndex);
			}
		}
	}
	
	const int iterationMax = 10000;
	int counter(0);
	bool continueLoop(true);
	while( continueLoop ){
		continueLoop = false;
		++counter;
		if( counter > iterationMax ){
			std::cerr << "Rearch max iteration number : " << iterationMax << std::endl;
			exit(1);
		}
		
		bool retrial = true;
		while(retrial){
			retrial = false;
			const std::set<int> elementsToBePartitionedPre = elementsToBePartitioned;
			for( std::set<int>::const_iterator itr = elementsToBePartitionedPre.begin(); itr != elementsToBePartitionedPre.end(); ++itr ){
				const int elemIndex = *itr;
				const int level = m_elemInfo[elemIndex].level;
				for( int iNeib = 0; iNeib < 6; ++iNeib ){
					const int levelNeib = m_elemInfo[elemIndex].levelNeib[iNeib];
					if( level <= levelNeib ){
						// Add element whose level is lower than the current owner element
						continue;
					}
					//Because the neigobor element has lower level, only one element is treated for each face
					const int elemIndexNeib = m_elemInfo[elemIndex].neib[iNeib][0];
					assert( m_elemInfo[elemIndex].neib[iNeib][1] < 0 &&
							m_elemInfo[elemIndex].neib[iNeib][2] < 0 &&
							m_elemInfo[elemIndex].neib[iNeib][3] < 0 );
					// Return to neighbors
					if( elemIndexNeib < 0 ){
						// Skip boundary plane
						continue;
					}
					if( !isActive(elemIndexNeib) ){
						// Skip inactive element
						continue;
					}
					if( elementsToBePartitioned.find(elemIndexNeib) == elementsToBePartitioned.end() ){
						// Insert only if the element has not been inserted
						elementsToBePartitioned.insert(elemIndexNeib);
						retrial = true;
						continueLoop = true;
					}
				}
			}
		}

		retrial = true;
		while(retrial){
			retrial = false;
			const std::set<int> elementsToBePartitionedPre = elementsToBePartitioned;
			for( std::set<int>::const_iterator itr = elementsToBePartitionedPre.begin(); itr != elementsToBePartitionedPre.end(); ++itr ){
				const int elemIndex = *itr;
				const int level = m_elemInfo[elemIndex].level;
				for( int iNeib = 0; iNeib < 6; ++iNeib ){
					if( getPartitioningType() == HORIZONTAL_ONLY && iNeib <= 3 ){
						continue;
					}
					const int levelNeib = m_elemInfo[elemIndex].levelNeib[iNeib];
					if( level < levelNeib ){
						// Add element whose level is lower than or equal to the current owner element
						continue;
					}
					//Because the neigobor element has lower or the same level, only one element is treated for each face
					const int elemIndexNeib = m_elemInfo[elemIndex].neib[iNeib][0];
					assert( m_elemInfo[elemIndex].neib[iNeib][1] < 0 &&
							m_elemInfo[elemIndex].neib[iNeib][2] < 0 &&
							m_elemInfo[elemIndex].neib[iNeib][3] < 0 );
					if( elemIndexNeib < 0 ){
						// Skip boundary plane
						continue;
					}
					// Neigbors of neighbors
					for( int iNeibNeib = 0; iNeibNeib < 6; ++iNeibNeib ){
						if( iNeibNeib == iNeib ){
							// Skip if the search direction is the same
							continue;
						}
						const int levelNeibNeib = m_elemInfo[elemIndexNeib].levelNeib[iNeibNeib];
						if( level <= levelNeibNeib ){
							// Add element whose level is lower than the current owner element
							continue;
						}
						//Because the neigobor element has lower level, only one element is treated for each face
						const int elemIndexNeibNeib = m_elemInfo[elemIndexNeib].neib[iNeibNeib][0];
						assert( m_elemInfo[elemIndexNeib].neib[iNeibNeib][1] < 0 &&
								m_elemInfo[elemIndexNeib].neib[iNeibNeib][2] < 0 &&
								m_elemInfo[elemIndexNeib].neib[iNeibNeib][3] < 0 );
						if( elemIndex == elemIndexNeibNeib  ){
							// Skip if element index is the same with the owner element
							continue;
						}
						if( elemIndexNeibNeib < 0 ){
							// Skip boundary plane
							continue;
						}
						if( !isActive(elemIndexNeibNeib) ){
							// Skip inactive element
							continue;
						}
						if( elementsToBePartitioned.find(elemIndexNeibNeib) == elementsToBePartitioned.end() ){
							// Insert only if the element has not been inserted
							elementsToBePartitioned.insert(elemIndexNeibNeib);
							retrial = true;
							continueLoop = true;
						}
					}
				}
			}
		}
	}

	const int numElemToBePartitioned = static_cast<int>( elementsToBePartitioned.size() );

	counter = 0;// Zero clear
	while( !elementsToBePartitioned.empty() ){
		++counter;
		if( counter > iterationMax ){
			std::cerr << "Rearch max iteration number" << std::endl;
			exit(1);
		}
		const std::set<int> elementsToBePartitionedPre = elementsToBePartitioned;
		for( std::set<int>::const_iterator itr = elementsToBePartitionedPre.begin(); itr != elementsToBePartitionedPre.end(); ++itr ){
			const int elemIndex = *itr;
			//const int level = m_elemInfo[elemIndex].level;
			//bool skip(false);
			//for( int i = 0; i < 6; ++i ){
			//	const std::vector<int>& neibElems = m_elemInfo[elemIndex].neib[i];
			//	for( std::vector<int>::const_iterator itrNeib = neibElems.begin(); itrNeib != neibElems.end(); ++itrNeib ){
			//		const int elemIndexNeib = *itrNeib;
			//		if( elemIndexNeib < 0 ){
			//			continue;
			//		}
			//		if( !isActive(elemIndexNeib) ){
			//			// Skip inactive element
			//			continue;
			//		}
			//		const int levelNeib = m_elemInfo[elemIndex].level;
			//		if( levelNeib < level ){
			//			skip = true;
			//			break;
			//		}
			//	}
			//}
			//if(skip){
			//	// Skip if the differene of the level become more than or equal to two.
			//	continue;
			//}
			const int level = m_elemInfo[elemIndex].level;
			bool skip(false);
			for( int iNeib = 0; iNeib < 6; ++iNeib ){
				const int elemIndexNeib = m_elemInfo[elemIndex].neib[iNeib][0];
				if( elemIndexNeib < 0 ){
					// Skip boundary plane
					continue;
				}
				if( !isActive(elemIndexNeib) ){
					// Skip inactive element
					continue;
				}
				if( m_elemInfo[elemIndex].level > m_elemInfo[elemIndex].levelNeib[iNeib] ){
					// Skip if the differene of the level become more than or equal to two.
					skip = true;
					break;
				}
			}
			if(skip){
				// Skip if the differene of the level become more than or equal to two.
				continue;
			}

			// Partitioning
			if( getPartitioningType() == FULL_PARTITIONING ){
				partitionAnElementFull(elemIndex);
			}else{
				partitionAnElementHorizontalDirectionOnly(elemIndex);
			}
			elementsToBePartitioned.erase(elemIndex);
		}
	}

	return numElemToBePartitioned;

}

// Partition an element
void MeshGenerationDHexa::partitionAnElementFull( const int elemIDParent ){

	// Partitioning
	const int newElemIDStart = static_cast<int>( m_elemInfo.size() );
	m_elemInfo[elemIDParent].isActive = false;
	for( int iChild = 0; iChild < 8; ++iChild ){
		m_elemInfo[elemIDParent].childs[iChild] = newElemIDStart + iChild;
	}

	const ElementInfo elemInfoParent = m_elemInfo[elemIDParent];
	CommonParameters::XYZ coordsPre[8];
	for( int i = 0; i < 8; ++i ){
		const int nodeID = elemInfoParent.nodes[i];
		coordsPre[i].X = m_nodeCoordinates[nodeID].X;
		coordsPre[i].Y = m_nodeCoordinates[nodeID].Y;
		coordsPre[i].Z = m_nodeCoordinates[nodeID].Z;
	}

	const int newNodeIDStart = static_cast<int>( m_nodeCoordinates.size() );
	CommonParameters::XYZ coordsNew[27];
	coordsNew[0]  = coordsPre[0];
	coordsNew[2]  = coordsPre[1];
	coordsNew[6]  = coordsPre[3];
	coordsNew[8]  = coordsPre[2];
	coordsNew[18] = coordsPre[4];
	coordsNew[20] = coordsPre[5];
	coordsNew[24] = coordsPre[7];
	coordsNew[26] = coordsPre[6];
	coordsNew[1]  = Util::averageTwoCoords( coordsPre[0], coordsPre[1] );
	coordsNew[3]  = Util::averageTwoCoords( coordsPre[0], coordsPre[3] );
	coordsNew[5]  = Util::averageTwoCoords( coordsPre[1], coordsPre[2] );
	coordsNew[7]  = Util::averageTwoCoords( coordsPre[3], coordsPre[2] );
	coordsNew[19] = Util::averageTwoCoords( coordsPre[4], coordsPre[5] );
	coordsNew[21] = Util::averageTwoCoords( coordsPre[4], coordsPre[7] );
	coordsNew[23] = Util::averageTwoCoords( coordsPre[5], coordsPre[6] );
	coordsNew[25] = Util::averageTwoCoords( coordsPre[7], coordsPre[6] );
	coordsNew[4]  = Util::averageTwoCoords( coordsNew[1],  coordsNew[7] );
	coordsNew[22] = Util::averageTwoCoords( coordsNew[19], coordsNew[25] );
	for( int i = 9; i < 18; ++i ){
		coordsNew[i] = Util::averageTwoCoords( coordsNew[i-9], coordsNew[i+9] );
	}
	int nodeIDNew[27];
#ifdef _MERGE
	nodeIDNew[0]  = elemInfoParent.nodes[0];
	nodeIDNew[2]  = elemInfoParent.nodes[1];
	nodeIDNew[6]  = elemInfoParent.nodes[3];
	nodeIDNew[8]  = elemInfoParent.nodes[2];
	nodeIDNew[18] = elemInfoParent.nodes[4];
	nodeIDNew[20] = elemInfoParent.nodes[5];
	nodeIDNew[24] = elemInfoParent.nodes[7];
	nodeIDNew[26] = elemInfoParent.nodes[6];
	nodeIDNew[1]  = newNodeIDStart;
	nodeIDNew[3]  = newNodeIDStart + 1;
	nodeIDNew[4]  = newNodeIDStart + 2;
	nodeIDNew[5]  = newNodeIDStart + 3;
	nodeIDNew[7]  = newNodeIDStart + 4;
	nodeIDNew[19] = newNodeIDStart + 14;
	nodeIDNew[21] = newNodeIDStart + 15;
	nodeIDNew[22] = newNodeIDStart + 16;
	nodeIDNew[23] = newNodeIDStart + 17;
	nodeIDNew[25] = newNodeIDStart + 18;
	for( int i = 0; i < 9; ++i ){
		nodeIDNew[i + 9] = newNodeIDStart + 5 + i;
	}

	m_nodeCoordinates.push_back(coordsNew[1]);
	m_nodeCoordinates.push_back(coordsNew[3]);
	m_nodeCoordinates.push_back(coordsNew[4]);
	m_nodeCoordinates.push_back(coordsNew[5]);
	m_nodeCoordinates.push_back(coordsNew[7]);
	for( int i = 0; i < 9; ++i ){
		m_nodeCoordinates.push_back(coordsNew[i+9]);
	}
	m_nodeCoordinates.push_back(coordsNew[19]);
	m_nodeCoordinates.push_back(coordsNew[21]);
	m_nodeCoordinates.push_back(coordsNew[22]);
	m_nodeCoordinates.push_back(coordsNew[23]);
	m_nodeCoordinates.push_back(coordsNew[25]);
#else
	for( int iNode = 0; iNode < 27; ++iNode ){
		nodeIDNew[iNode] = -1;
	}
	nodeIDNew[0]  = elemInfoParent.nodes[0];
	nodeIDNew[2]  = elemInfoParent.nodes[1];
	nodeIDNew[6]  = elemInfoParent.nodes[3];
	nodeIDNew[8]  = elemInfoParent.nodes[2];
	nodeIDNew[18] = elemInfoParent.nodes[4];
	nodeIDNew[20] = elemInfoParent.nodes[5];
	nodeIDNew[24] = elemInfoParent.nodes[7];
	nodeIDNew[26] = elemInfoParent.nodes[6];
	if( elemInfoParent.levelNeib[0] > elemInfoParent.level ){
		const int neibElem0 = elemInfoParent.neib[0][0];
		const int neibElem1 = elemInfoParent.neib[0][1];
		const int neibElem2 = elemInfoParent.neib[0][2];
		nodeIDNew[ 3] = m_elemInfo[neibElem0].nodes[2];
		nodeIDNew[ 9] = m_elemInfo[neibElem0].nodes[5];
		nodeIDNew[12] = m_elemInfo[neibElem0].nodes[6];
		nodeIDNew[15] = m_elemInfo[neibElem1].nodes[6];
		nodeIDNew[21] = m_elemInfo[neibElem2].nodes[6];
	}else if( elemInfoParent.neib[0][0] >= 0 ){
		const int neibElem = elemInfoParent.neib[0][0];
		const ElementInfo elemInfoNeib = m_elemInfo[neibElem];
		assert(elemInfoNeib.level == elemInfoParent.level);
		if( elemInfoNeib.levelNeib[2] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[2][1];
			nodeIDNew[ 9] = m_elemInfo[neibNeibElem].nodes[6];
		}
		if( elemInfoNeib.levelNeib[3] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[3][1];
			nodeIDNew[15] = m_elemInfo[neibNeibElem].nodes[5];
		}
		if( elemInfoNeib.levelNeib[4] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[4][1];
			nodeIDNew[ 3] = m_elemInfo[neibNeibElem].nodes[6];
		}
		if( elemInfoNeib.levelNeib[5] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[5][1];
			nodeIDNew[21] = m_elemInfo[neibNeibElem].nodes[2];
		}
	}
	if( elemInfoParent.levelNeib[1] > elemInfoParent.level ){
		const int neibElem0 = elemInfoParent.neib[1][0];
		const int neibElem1 = elemInfoParent.neib[1][1];
		const int neibElem2 = elemInfoParent.neib[1][2];
		nodeIDNew[ 5] = m_elemInfo[neibElem0].nodes[3];
		nodeIDNew[11] = m_elemInfo[neibElem0].nodes[4];
		nodeIDNew[14] = m_elemInfo[neibElem0].nodes[7];
		nodeIDNew[17] = m_elemInfo[neibElem1].nodes[7];
		nodeIDNew[23] = m_elemInfo[neibElem2].nodes[7];
	}else if( elemInfoParent.neib[1][0] >= 0 ){
		const int neibElem = elemInfoParent.neib[1][0];
		const ElementInfo elemInfoNeib = m_elemInfo[neibElem];
		assert(elemInfoNeib.level == elemInfoParent.level);
		if( elemInfoNeib.levelNeib[2] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[2][0];
			nodeIDNew[11] = m_elemInfo[neibNeibElem].nodes[7];
		}
		if( elemInfoNeib.levelNeib[3] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[3][0];
			nodeIDNew[17] = m_elemInfo[neibNeibElem].nodes[4];
		}
		if( elemInfoNeib.levelNeib[4] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[4][0];
			nodeIDNew[ 5] = m_elemInfo[neibNeibElem].nodes[7];
		}
		if( elemInfoNeib.levelNeib[5] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[5][0];
			nodeIDNew[23] = m_elemInfo[neibNeibElem].nodes[3];
		}
	}
	if( elemInfoParent.levelNeib[2] > elemInfoParent.level ){
		const int neibElem0 = elemInfoParent.neib[2][0];
		const int neibElem1 = elemInfoParent.neib[2][1];
		const int neibElem2 = elemInfoParent.neib[2][2];
		nodeIDNew[ 1] = m_elemInfo[neibElem0].nodes[2];
		nodeIDNew[ 9] = m_elemInfo[neibElem0].nodes[7];
		nodeIDNew[10] = m_elemInfo[neibElem0].nodes[6];
		nodeIDNew[11] = m_elemInfo[neibElem1].nodes[6];
		nodeIDNew[19] = m_elemInfo[neibElem2].nodes[6];
	}else if( elemInfoParent.neib[2][0] >= 0 ){
		const int neibElem = elemInfoParent.neib[2][0];
		const ElementInfo elemInfoNeib = m_elemInfo[neibElem];
		assert(elemInfoNeib.level == elemInfoParent.level);
		if( elemInfoNeib.levelNeib[4] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[4][2];
			nodeIDNew[ 1] = m_elemInfo[neibNeibElem].nodes[6];
		}
		if( elemInfoNeib.levelNeib[5] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[5][2];
			nodeIDNew[19] = m_elemInfo[neibNeibElem].nodes[2];
		}
	}
	if( elemInfoParent.levelNeib[3] > elemInfoParent.level ){
		const int neibElem0 = elemInfoParent.neib[3][0];
		const int neibElem1 = elemInfoParent.neib[3][1];
		const int neibElem2 = elemInfoParent.neib[3][2];
		nodeIDNew[ 7] = m_elemInfo[neibElem0].nodes[1];
		nodeIDNew[15] = m_elemInfo[neibElem0].nodes[4];
		nodeIDNew[16] = m_elemInfo[neibElem0].nodes[5];
		nodeIDNew[17] = m_elemInfo[neibElem1].nodes[5];
		nodeIDNew[25] = m_elemInfo[neibElem2].nodes[5];
	}else if( elemInfoParent.neib[3][0] >= 0 ){
		const int neibElem = elemInfoParent.neib[3][0];
		const ElementInfo elemInfoNeib = m_elemInfo[neibElem];
		assert(elemInfoNeib.level == elemInfoParent.level);
		if( elemInfoNeib.levelNeib[4] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[4][0];
			nodeIDNew[ 7] = m_elemInfo[neibNeibElem].nodes[5];
		}
		if( elemInfoNeib.levelNeib[5] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[5][0];
			nodeIDNew[25] = m_elemInfo[neibNeibElem].nodes[1];
		}
	}
	if( elemInfoParent.levelNeib[4] > elemInfoParent.level ){
		const int neibElem0 = elemInfoParent.neib[4][0];
		const int neibElem1 = elemInfoParent.neib[4][1];
		const int neibElem2 = elemInfoParent.neib[4][2];
		nodeIDNew[ 1] = m_elemInfo[neibElem0].nodes[5];
		nodeIDNew[ 3] = m_elemInfo[neibElem0].nodes[7];
		nodeIDNew[ 4] = m_elemInfo[neibElem0].nodes[6];
		nodeIDNew[ 5] = m_elemInfo[neibElem1].nodes[6];
		nodeIDNew[ 7] = m_elemInfo[neibElem2].nodes[6];
	}
	if( elemInfoParent.levelNeib[5] > elemInfoParent.level ){
		const int neibElem0 = elemInfoParent.neib[5][0];
		const int neibElem1 = elemInfoParent.neib[5][1];
		const int neibElem2 = elemInfoParent.neib[5][2];
		nodeIDNew[19] = m_elemInfo[neibElem0].nodes[1];
		nodeIDNew[21] = m_elemInfo[neibElem0].nodes[3];
		nodeIDNew[22] = m_elemInfo[neibElem0].nodes[2];
		nodeIDNew[23] = m_elemInfo[neibElem1].nodes[2];
		nodeIDNew[25] = m_elemInfo[neibElem2].nodes[2];
	}
	int icount(0);
	for( int iNode = 0; iNode < 27; ++iNode ){
		if( nodeIDNew[iNode] < 0 ){
			m_nodeCoordinates.push_back(coordsNew[iNode]);
			nodeIDNew[iNode] = newNodeIDStart + icount;
			++icount;
		}
	}
#endif

	int counter(0);
	for( int iz = 0; iz < 2; ++iz ){
		for( int iy = 0; iy < 2; ++iy ){
			for( int ix = 0; ix < 2; ++ix ){
				const int elemIDNew = newElemIDStart + counter;
				ElementInfo info;
				const int offset = ix + iy * 3 + iz * 9;
				info.nodes[0] = nodeIDNew[0  + offset];
				info.nodes[1] = nodeIDNew[1  + offset];
				info.nodes[2] = nodeIDNew[4  + offset];
				info.nodes[3] = nodeIDNew[3  + offset];
				info.nodes[4] = nodeIDNew[9  + offset];
				info.nodes[5] = nodeIDNew[10 + offset];
				info.nodes[6] = nodeIDNew[13 + offset];
				info.nodes[7] = nodeIDNew[12 + offset];
				info.type = elemInfoParent.type;
				info.isActive = true;
				info.parent = elemIDParent;
				for( int i = 0; i < 8; ++i ){
					info.childs[i] = -1;
				}
				info.level = elemInfoParent.level + 1;
				info.parameterCell = elemInfoParent.parameterCell;
				for( int i = 0; i < 6; ++i ){
					for( int j = 0; j < 4; ++j ){
						info.neib[i][j] = -1;
					}
					info.levelNeib[i] = elemInfoParent.levelNeib[i];
				}
				// Y-Z Plane
				if( ix == 0 ){
					int num1(0);
					const int num = ( elemInfoParent.levelNeib[0] > elemInfoParent.level ) ? iy+2*iz : 0;
					const int elemNeib = elemInfoParent.neib[0][num];
					info.neib[0][0] = elemNeib;
					if( elemNeib >= 0 ){
						const int num = ( info.level > elemInfoParent.levelNeib[0] ) ? iy+2*iz : 0;
						m_elemInfo[elemNeib].neib[1][num] = elemIDNew;
						m_elemInfo[elemNeib].levelNeib[1] = info.level;
					}
					info.neib[1][0] = elemIDNew + 1;
					info.levelNeib[1] = info.level;
				}else{
					const int num = ( elemInfoParent.levelNeib[1] > elemInfoParent.level ) ? iy+2*iz : 0;
					const int elemNeib = elemInfoParent.neib[1][num];
					info.neib[1][0] = elemNeib;
					if( elemNeib >= 0 ){
						const int num = ( info.level > elemInfoParent.levelNeib[1] ) ? iy+2*iz : 0;
						m_elemInfo[elemNeib].neib[0][num] = elemIDNew;
						m_elemInfo[elemNeib].levelNeib[0] = info.level;
					}
					info.neib[0][0] = elemIDNew - 1;
					info.levelNeib[0] = info.level;
				}
				// Z-X Plane
				if( iy == 0 ){
					const int num = ( elemInfoParent.levelNeib[2] > elemInfoParent.level ) ? ix+2*iz : 0;
					const int elemNeib = elemInfoParent.neib[2][num];
					info.neib[2][0] = elemNeib;
					if( elemNeib >= 0 ){
						const int num = ( info.level > elemInfoParent.levelNeib[2] ) ? ix+2*iz : 0;
						m_elemInfo[elemNeib].neib[3][num] = elemIDNew;
						m_elemInfo[elemNeib].levelNeib[3] = info.level;
					}
					info.neib[3][0] = elemIDNew + 2;
					info.levelNeib[3] = info.level;
				}else{
					const int num = ( elemInfoParent.levelNeib[3] > elemInfoParent.level ) ? ix+2*iz : 0;
					const int elemNeib = elemInfoParent.neib[3][num];
					info.neib[3][0] = elemNeib;
					if( elemNeib >= 0 ){
						const int num = ( info.level > elemInfoParent.levelNeib[3] ) ? ix+2*iz : 0;
						m_elemInfo[elemNeib].neib[2][num] = elemIDNew;
						m_elemInfo[elemNeib].levelNeib[2] = info.level;
					}
					info.neib[2][0] = elemIDNew - 2;
					info.levelNeib[2] = info.level;
				}
				// X-Y Plane
				if( iz == 0 ){
					const int num = ( elemInfoParent.levelNeib[4] > elemInfoParent.level ) ? ix+2*iy : 0;
					const int elemNeib = elemInfoParent.neib[4][num];
					info.neib[4][0] = elemNeib;
					if( elemNeib >= 0 ){
						const int num = ( info.level > elemInfoParent.levelNeib[4] ) ? ix+2*iy : 0;
						m_elemInfo[elemNeib].neib[5][num] = elemIDNew;
						m_elemInfo[elemNeib].levelNeib[5] = info.level;
					}
					info.neib[5][0] = elemIDNew + 4;
					info.levelNeib[5] = info.level;
				}else{
					const int num = ( elemInfoParent.levelNeib[5] > elemInfoParent.level ) ? ix+2*iy : 0;
					const int elemNeib = elemInfoParent.neib[5][num];
					info.neib[5][0] = elemNeib;
					if( elemNeib >= 0 ){
						const int num = ( info.level > elemInfoParent.levelNeib[5] ) ? ix+2*iy : 0;
						m_elemInfo[elemNeib].neib[4][num] = elemIDNew;
						m_elemInfo[elemNeib].levelNeib[4] = info.level;
					}
					info.neib[4][0] = elemIDNew - 4;
					info.levelNeib[4] = info.level;
				}
#ifdef _DEBUG_WRITE
				checkWhetherFourEdgesAreVertical(info);
#endif
				m_elemInfo.push_back(info);
				++counter;
			}
		}
	}

	changeElementsOnBoundaryPlanesFull( elemIDParent, newElemIDStart );

}

// Partition an element only in horizontal direction
void MeshGenerationDHexa::partitionAnElementHorizontalDirectionOnly( const int elemIDParent ){

	// Partitioning
	const int newElemIDStart = static_cast<int>( m_elemInfo.size() );
	m_elemInfo[elemIDParent].isActive = false;
	for( int iChild = 0; iChild < 4; ++iChild ){
		m_elemInfo[elemIDParent].childs[iChild] = newElemIDStart + iChild;
	}

	const ElementInfo elemInfoParent = m_elemInfo[elemIDParent];
	CommonParameters::XYZ coordsPre[8];
	for( int i = 0; i < 8; ++i ){
		const int nodeID = elemInfoParent.nodes[i];
		coordsPre[i].X = m_nodeCoordinates[nodeID].X;
		coordsPre[i].Y = m_nodeCoordinates[nodeID].Y;
		coordsPre[i].Z = m_nodeCoordinates[nodeID].Z;
	}

	const int newNodeIDStart = static_cast<int>( m_nodeCoordinates.size() );
	CommonParameters::XYZ coordsNew[18];
	coordsNew[0]  = coordsPre[0];
	coordsNew[2]  = coordsPre[1];
	coordsNew[6]  = coordsPre[3];
	coordsNew[8]  = coordsPre[2];
	coordsNew[9]  = coordsPre[4];
	coordsNew[11] = coordsPre[5];
	coordsNew[15] = coordsPre[7];
	coordsNew[17] = coordsPre[6];
	coordsNew[1]  = Util::averageTwoCoords( coordsPre[0], coordsPre[1] );
	coordsNew[3]  = Util::averageTwoCoords( coordsPre[0], coordsPre[3] );
	coordsNew[5]  = Util::averageTwoCoords( coordsPre[1], coordsPre[2] );
	coordsNew[7]  = Util::averageTwoCoords( coordsPre[3], coordsPre[2] );
	coordsNew[10] = Util::averageTwoCoords( coordsPre[4], coordsPre[5] );
	coordsNew[12] = Util::averageTwoCoords( coordsPre[4], coordsPre[7] );
	coordsNew[14] = Util::averageTwoCoords( coordsPre[5], coordsPre[6] );
	coordsNew[16] = Util::averageTwoCoords( coordsPre[7], coordsPre[6] );
	coordsNew[4]  = Util::averageTwoCoords( coordsNew[1],  coordsNew[7] );
	coordsNew[13] = Util::averageTwoCoords( coordsNew[10], coordsNew[16] );
	int nodeIDNew[18];
#ifdef _MERGE
	nodeIDNew[0]  = elemInfoParent.nodes[0];
	nodeIDNew[2]  = elemInfoParent.nodes[1];
	nodeIDNew[6]  = elemInfoParent.nodes[3];
	nodeIDNew[8]  = elemInfoParent.nodes[2];
	nodeIDNew[9]  = elemInfoParent.nodes[4];
	nodeIDNew[11] = elemInfoParent.nodes[5];
	nodeIDNew[15] = elemInfoParent.nodes[7];
	nodeIDNew[17] = elemInfoParent.nodes[6];
	nodeIDNew[1]  = newNodeIDStart;
	nodeIDNew[3]  = newNodeIDStart + 1;
	nodeIDNew[4]  = newNodeIDStart + 2;
	nodeIDNew[5]  = newNodeIDStart + 3;
	nodeIDNew[7]  = newNodeIDStart + 4;
	nodeIDNew[10] = newNodeIDStart + 5;
	nodeIDNew[12] = newNodeIDStart + 6;
	nodeIDNew[13] = newNodeIDStart + 7;
	nodeIDNew[14] = newNodeIDStart + 8;
	nodeIDNew[16] = newNodeIDStart + 9;

	m_nodeCoordinates.push_back(coordsNew[1]);
	m_nodeCoordinates.push_back(coordsNew[3]);
	m_nodeCoordinates.push_back(coordsNew[4]);
	m_nodeCoordinates.push_back(coordsNew[5]);
	m_nodeCoordinates.push_back(coordsNew[7]);
	m_nodeCoordinates.push_back(coordsNew[10]);
	m_nodeCoordinates.push_back(coordsNew[12]);
	m_nodeCoordinates.push_back(coordsNew[13]);
	m_nodeCoordinates.push_back(coordsNew[14]);
	m_nodeCoordinates.push_back(coordsNew[16]);
#else
	for( int iNode = 0; iNode < 18; ++iNode ){
		nodeIDNew[iNode] = -1;
	}
	nodeIDNew[0]  = elemInfoParent.nodes[0];
	nodeIDNew[2]  = elemInfoParent.nodes[1];
	nodeIDNew[6]  = elemInfoParent.nodes[3];
	nodeIDNew[8]  = elemInfoParent.nodes[2];
	nodeIDNew[9]  = elemInfoParent.nodes[4];
	nodeIDNew[11] = elemInfoParent.nodes[5];
	nodeIDNew[15] = elemInfoParent.nodes[7];
	nodeIDNew[17] = elemInfoParent.nodes[6];
	if( elemInfoParent.levelNeib[0] > elemInfoParent.level ){
		const int neibElem0 = elemInfoParent.neib[0][0];
		nodeIDNew[ 3] = m_elemInfo[neibElem0].nodes[2];
		nodeIDNew[12] = m_elemInfo[neibElem0].nodes[6];
	}else if( elemInfoParent.neib[0][0] >= 0 ){
		const int neibElem = elemInfoParent.neib[0][0];
		const ElementInfo elemInfoNeib = m_elemInfo[neibElem];
		assert(elemInfoNeib.level == elemInfoParent.level);
		if( elemInfoNeib.levelNeib[4] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[4][3];
			nodeIDNew[ 3] = m_elemInfo[neibNeibElem].nodes[5];
		}
		if( elemInfoNeib.levelNeib[5] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[5][3];
			nodeIDNew[12] = m_elemInfo[neibNeibElem].nodes[1];
		}
	}
	if( elemInfoParent.levelNeib[1] > elemInfoParent.level ){
		const int neibElem0 = elemInfoParent.neib[1][0];
		nodeIDNew[ 5] = m_elemInfo[neibElem0].nodes[3];
		nodeIDNew[14] = m_elemInfo[neibElem0].nodes[7];
	}else if( elemInfoParent.neib[1][0] >= 0 ){
		const int neibElem = elemInfoParent.neib[1][0];
		const ElementInfo elemInfoNeib = m_elemInfo[neibElem];
		assert(elemInfoNeib.level == elemInfoParent.level);
		if( elemInfoNeib.levelNeib[4] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[4][0];
			nodeIDNew[ 5] = m_elemInfo[neibNeibElem].nodes[7];
		}
		if( elemInfoNeib.levelNeib[5] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[5][0];
			nodeIDNew[14] = m_elemInfo[neibNeibElem].nodes[3];
		}
	}
	if( elemInfoParent.levelNeib[2] > elemInfoParent.level ){
		const int neibElem0 = elemInfoParent.neib[2][0];
		nodeIDNew[ 1] = m_elemInfo[neibElem0].nodes[2];
		nodeIDNew[10] = m_elemInfo[neibElem0].nodes[6];
	}else if( elemInfoParent.neib[2][0] >= 0 ){
		const int neibElem = elemInfoParent.neib[2][0];
		const ElementInfo elemInfoNeib = m_elemInfo[neibElem];
		assert(elemInfoNeib.level == elemInfoParent.level);
		if( elemInfoNeib.levelNeib[4] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[4][3];
			nodeIDNew[ 1] = m_elemInfo[neibNeibElem].nodes[7];
		}
		if( elemInfoNeib.levelNeib[5] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[5][3];
			nodeIDNew[10] = m_elemInfo[neibNeibElem].nodes[3];
		}
	}
	if( elemInfoParent.levelNeib[3] > elemInfoParent.level ){
		const int neibElem0 = elemInfoParent.neib[3][0];
		nodeIDNew[ 7] = m_elemInfo[neibElem0].nodes[1];
		nodeIDNew[16] = m_elemInfo[neibElem0].nodes[5];
	}else if( elemInfoParent.neib[3][0] >= 0 ){
		const int neibElem = elemInfoParent.neib[3][0];
		const ElementInfo elemInfoNeib = m_elemInfo[neibElem];
		assert(elemInfoNeib.level == elemInfoParent.level);
		if( elemInfoNeib.levelNeib[4] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[4][0];
			nodeIDNew[ 7] = m_elemInfo[neibNeibElem].nodes[5];
		}
		if( elemInfoNeib.levelNeib[5] > elemInfoParent.level ){
			const int neibNeibElem = elemInfoNeib.neib[5][0];
			nodeIDNew[16] = m_elemInfo[neibNeibElem].nodes[1];
		}
	}
	if( elemInfoParent.levelNeib[4] > elemInfoParent.level ){
		const int neibElem0 = elemInfoParent.neib[4][0];
		const int neibElem1 = elemInfoParent.neib[4][1];
		const int neibElem2 = elemInfoParent.neib[4][2];
		nodeIDNew[ 1] = m_elemInfo[neibElem0].nodes[5];
		nodeIDNew[ 3] = m_elemInfo[neibElem0].nodes[7];
		nodeIDNew[ 4] = m_elemInfo[neibElem0].nodes[6];
		nodeIDNew[ 5] = m_elemInfo[neibElem1].nodes[6];
		nodeIDNew[ 7] = m_elemInfo[neibElem2].nodes[6];
	}
	if( elemInfoParent.levelNeib[5] > elemInfoParent.level ){
		const int neibElem0 = elemInfoParent.neib[5][0];
		const int neibElem1 = elemInfoParent.neib[5][1];
		const int neibElem2 = elemInfoParent.neib[5][2];
		nodeIDNew[10] = m_elemInfo[neibElem0].nodes[1];
		nodeIDNew[12] = m_elemInfo[neibElem0].nodes[3];
		nodeIDNew[13] = m_elemInfo[neibElem0].nodes[2];
		nodeIDNew[14] = m_elemInfo[neibElem1].nodes[2];
		nodeIDNew[16] = m_elemInfo[neibElem2].nodes[2];
	}
	int icount(0);
	for( int iNode = 0; iNode < 18; ++iNode ){
		if( nodeIDNew[iNode] < 0 ){
			m_nodeCoordinates.push_back(coordsNew[iNode]);
			nodeIDNew[iNode] = newNodeIDStart + icount;
			++icount;
		}
	}
#endif

	int counter(0);
	for( int iy = 0; iy < 2; ++iy ){
		for( int ix = 0; ix < 2; ++ix ){
			const int elemIDNew = newElemIDStart + counter;
			ElementInfo info;
			const int offset = ix + iy * 3;
			info.nodes[0] = nodeIDNew[0  + offset];
			info.nodes[1] = nodeIDNew[1  + offset];
			info.nodes[2] = nodeIDNew[4  + offset];
			info.nodes[3] = nodeIDNew[3  + offset];
			info.nodes[4] = nodeIDNew[9  + offset];
			info.nodes[5] = nodeIDNew[10 + offset];
			info.nodes[6] = nodeIDNew[13 + offset];
			info.nodes[7] = nodeIDNew[12 + offset];
			info.type = elemInfoParent.type;
			info.isActive = true;
			info.parent = elemIDParent;
			for( int i = 0; i < 8; ++i ){
				info.childs[i] = -1;
			}
			info.level = elemInfoParent.level + 1;
			info.parameterCell = elemInfoParent.parameterCell;
			for( int i = 0; i < 6; ++i ){
				for( int j = 0; j < 4; ++j ){
					info.neib[i][j] = -1;
				}
				info.levelNeib[i] = elemInfoParent.levelNeib[i];
			}
			// Y-Z Plane
			if( ix == 0 ){
				int num1(0);
				const int num = ( elemInfoParent.levelNeib[0] > elemInfoParent.level ) ? iy : 0;
				const int elemNeib = elemInfoParent.neib[0][num];
				info.neib[0][0] = elemNeib;
				if( elemNeib >= 0 ){
					const int num = ( info.level > elemInfoParent.levelNeib[0] ) ? iy : 0;
					m_elemInfo[elemNeib].neib[1][num] = elemIDNew;
					m_elemInfo[elemNeib].levelNeib[1] = info.level;
				}
				info.neib[1][0] = elemIDNew + 1;
				info.levelNeib[1] = info.level;
			}else{
				const int num = ( elemInfoParent.levelNeib[1] > elemInfoParent.level ) ? iy : 0;
				const int elemNeib = elemInfoParent.neib[1][num];
				info.neib[1][0] = elemNeib;
				if( elemNeib >= 0 ){
					const int num = ( info.level > elemInfoParent.levelNeib[1] ) ? iy : 0;
					m_elemInfo[elemNeib].neib[0][num] = elemIDNew;
					m_elemInfo[elemNeib].levelNeib[0] = info.level;
				}
				info.neib[0][0] = elemIDNew - 1;
				info.levelNeib[0] = info.level;
			}
			// Z-X Plane
			if( iy == 0 ){
				const int num = ( elemInfoParent.levelNeib[2] > elemInfoParent.level ) ? ix : 0;
				const int elemNeib = elemInfoParent.neib[2][num];
				info.neib[2][0] = elemNeib;
				if( elemNeib >= 0 ){
					const int num = ( info.level > elemInfoParent.levelNeib[2] ) ? ix : 0;
					m_elemInfo[elemNeib].neib[3][num] = elemIDNew;
					m_elemInfo[elemNeib].levelNeib[3] = info.level;
				}
				info.neib[3][0] = elemIDNew + 2;
				info.levelNeib[3] = info.level;
			}else{
				const int num = ( elemInfoParent.levelNeib[3] > elemInfoParent.level ) ? ix : 0;
				const int elemNeib = elemInfoParent.neib[3][num];
				info.neib[3][0] = elemNeib;
				if( elemNeib >= 0 ){
					const int num = ( info.level > elemInfoParent.levelNeib[3] ) ? ix : 0;
					m_elemInfo[elemNeib].neib[2][num] = elemIDNew;
					m_elemInfo[elemNeib].levelNeib[2] = info.level;
				}
				info.neib[2][0] = elemIDNew - 2;
				info.levelNeib[2] = info.level;
			}
			// X-Y Plane
			{
				const int num = ( elemInfoParent.levelNeib[4] > elemInfoParent.level ) ? ix+2*iy : 0;
				const int elemNeib = elemInfoParent.neib[4][num];
				info.neib[4][0] = elemNeib;
				if( elemNeib >= 0 ){
					const int num = ( info.level > elemInfoParent.levelNeib[4] ) ? ix+2*iy : 0;
					m_elemInfo[elemNeib].neib[5][num] = elemIDNew;
					m_elemInfo[elemNeib].levelNeib[5] = info.level;
				}
			}
			{
				const int num = ( elemInfoParent.levelNeib[5] > elemInfoParent.level ) ? ix+2*iy : 0;
				const int elemNeib = elemInfoParent.neib[5][num];
				info.neib[5][0] = elemNeib;
				if( elemNeib >= 0 ){
					const int num = ( info.level > elemInfoParent.levelNeib[5] ) ? ix+2*iy : 0;
					m_elemInfo[elemNeib].neib[4][num] = elemIDNew;
					m_elemInfo[elemNeib].levelNeib[4] = info.level;
				}
			}
#ifdef _DEBUG_WRITE
			checkWhetherFourEdgesAreVertical(info);
#endif
			m_elemInfo.push_back(info);
			++counter;
		}
	}

	changeElementsOnBoundaryPlanesHorizontalDirectionOnly( elemIDParent, newElemIDStart );

}

// Change elements on boundary plane
void MeshGenerationDHexa::changeElementsOnBoundaryPlanesFull( const int elemParent, const int newElemIDStart ){


	{// YZMinus
		std::set<int>& elems = m_elementOnBoundaryPlanes[YZMinus];
		std::set<int>::const_iterator itr = elems.find(elemParent);
		if( itr != elems.end() ){
			elems.erase(itr);
			elems.insert(newElemIDStart);
			elems.insert(newElemIDStart+2); 
			elems.insert(newElemIDStart+4);
			elems.insert(newElemIDStart+6);
		}
	}

	{// YZPlus
		std::set<int>& elems = m_elementOnBoundaryPlanes[YZPlus];
		std::set<int>::const_iterator itr = elems.find(elemParent);
		if( itr != elems.end() ){
			elems.erase(itr);
			elems.insert(newElemIDStart+1);
			elems.insert(newElemIDStart+3); 
			elems.insert(newElemIDStart+5);
			elems.insert(newElemIDStart+7);
		}
	}

	{// ZXMinus
		std::set<int>& elems = m_elementOnBoundaryPlanes[ZXMinus];
		std::set<int>::const_iterator itr = elems.find(elemParent);
		if( itr != elems.end() ){
			elems.erase(itr);
			elems.insert(newElemIDStart);
			elems.insert(newElemIDStart+1); 
			elems.insert(newElemIDStart+4);
			elems.insert(newElemIDStart+5);
		}
	}

	{// ZXPlus
		std::set<int>& elems = m_elementOnBoundaryPlanes[ZXPlus];
		std::set<int>::const_iterator itr = elems.find(elemParent);
		if( itr != elems.end() ){
			elems.erase(itr);
			elems.insert(newElemIDStart+2);
			elems.insert(newElemIDStart+3); 
			elems.insert(newElemIDStart+6);
			elems.insert(newElemIDStart+7);
		}
	}

	{// XYMinus
		std::set<int>& elems = m_elementOnBoundaryPlanes[XYMinus];
		std::set<int>::const_iterator itr = elems.find(elemParent);
		if( itr != elems.end() ){
			elems.erase(itr);
			elems.insert(newElemIDStart);
			elems.insert(newElemIDStart+1); 
			elems.insert(newElemIDStart+2);
			elems.insert(newElemIDStart+3);
		}
	}

	{// XYPlus
		std::set<int>& elems = m_elementOnBoundaryPlanes[XYPlus];
		std::set<int>::const_iterator itr = elems.find(elemParent);
		if( itr != elems.end() ){
			elems.erase(itr);
			elems.insert(newElemIDStart+4);
			elems.insert(newElemIDStart+5); 
			elems.insert(newElemIDStart+6);
			elems.insert(newElemIDStart+7);
		}
	}

	{// Earth surface
		std::set<int>& elems = m_elementOfEarthSurface;
		std::set<int>::const_iterator itr = elems.find(elemParent);
		if( itr != elems.end() ){
			elems.erase(itr);
			elems.insert(newElemIDStart);
			elems.insert(newElemIDStart+1); 
			elems.insert(newElemIDStart+2);
			elems.insert(newElemIDStart+3);
			m_elementOfEarthSurfaceWithInactiveElements.insert(newElemIDStart);
			m_elementOfEarthSurfaceWithInactiveElements.insert(newElemIDStart+1); 
			m_elementOfEarthSurfaceWithInactiveElements.insert(newElemIDStart+2);
			m_elementOfEarthSurfaceWithInactiveElements.insert(newElemIDStart+3);
		}
	}

	// Sea surface
	if( m_includeSea ){
		std::set<int>& elems = m_elementOfSeaSurface;
		std::set<int>::const_iterator itr = elems.find(elemParent);
		if( itr != elems.end() ){
			elems.erase(itr);
			elems.insert(newElemIDStart);
			elems.insert(newElemIDStart+1); 
			elems.insert(newElemIDStart+2);
			elems.insert(newElemIDStart+3);
			m_elementOfSeaSurfaceWithInactiveElements.insert(newElemIDStart);
			m_elementOfSeaSurfaceWithInactiveElements.insert(newElemIDStart+1); 
			m_elementOfSeaSurfaceWithInactiveElements.insert(newElemIDStart+2);
			m_elementOfSeaSurfaceWithInactiveElements.insert(newElemIDStart+3);
		}
	}

}

// Change elements on boundary plane only in horizontal direction
void MeshGenerationDHexa::changeElementsOnBoundaryPlanesHorizontalDirectionOnly( const int elemParent, const int newElemIDStart ){

	{// YZMinus
		std::set<int>& elems = m_elementOnBoundaryPlanes[YZMinus];
		std::set<int>::const_iterator itr = elems.find(elemParent);
		if( itr != elems.end() ){
			elems.erase(itr);
			elems.insert(newElemIDStart);
			elems.insert(newElemIDStart+2); 
		}
	}

	{// YZPlus
		std::set<int>& elems = m_elementOnBoundaryPlanes[YZPlus];
		std::set<int>::const_iterator itr = elems.find(elemParent);
		if( itr != elems.end() ){
			elems.erase(itr);
			elems.insert(newElemIDStart+1);
			elems.insert(newElemIDStart+3); 
		}
	}

	{// ZXMinus
		std::set<int>& elems = m_elementOnBoundaryPlanes[ZXMinus];
		std::set<int>::const_iterator itr = elems.find(elemParent);
		if( itr != elems.end() ){
			elems.erase(itr);
			elems.insert(newElemIDStart);
			elems.insert(newElemIDStart+1); 
		}
	}

	{// ZXPlus
		std::set<int>& elems = m_elementOnBoundaryPlanes[ZXPlus];
		std::set<int>::const_iterator itr = elems.find(elemParent);
		if( itr != elems.end() ){
			elems.erase(itr);
			elems.insert(newElemIDStart+2);
			elems.insert(newElemIDStart+3); 
		}
	}

	{// XYMinus
		std::set<int>& elems = m_elementOnBoundaryPlanes[XYMinus];
		std::set<int>::const_iterator itr = elems.find(elemParent);
		if( itr != elems.end() ){
			elems.erase(itr);
			elems.insert(newElemIDStart);
			elems.insert(newElemIDStart+1); 
			elems.insert(newElemIDStart+2);
			elems.insert(newElemIDStart+3);
		}
	}

	{// XYPlus
		std::set<int>& elems = m_elementOnBoundaryPlanes[XYPlus];
		std::set<int>::const_iterator itr = elems.find(elemParent);
		if( itr != elems.end() ){
			elems.erase(itr);
			elems.insert(newElemIDStart);
			elems.insert(newElemIDStart+1); 
			elems.insert(newElemIDStart+2);
			elems.insert(newElemIDStart+3);
		}
	}

	{// Earth surface
		std::set<int>& elems = m_elementOfEarthSurface;
		std::set<int>::const_iterator itr = elems.find(elemParent);
		if( itr != elems.end() ){
			elems.erase(itr);
			elems.insert(newElemIDStart);
			elems.insert(newElemIDStart+1); 
			elems.insert(newElemIDStart+2);
			elems.insert(newElemIDStart+3);
			m_elementOfEarthSurfaceWithInactiveElements.insert(newElemIDStart);
			m_elementOfEarthSurfaceWithInactiveElements.insert(newElemIDStart+1); 
			m_elementOfEarthSurfaceWithInactiveElements.insert(newElemIDStart+2);
			m_elementOfEarthSurfaceWithInactiveElements.insert(newElemIDStart+3);
		}
	}

	// Sea surface
	if( m_includeSea ){
		std::set<int>& elems = m_elementOfSeaSurface;
		std::set<int>::const_iterator itr = elems.find(elemParent);
		if( itr != elems.end() ){
			elems.erase(itr);
			elems.insert(newElemIDStart);
			elems.insert(newElemIDStart+1); 
			elems.insert(newElemIDStart+2);
			elems.insert(newElemIDStart+3);
			m_elementOfSeaSurfaceWithInactiveElements.insert(newElemIDStart);
			m_elementOfSeaSurfaceWithInactiveElements.insert(newElemIDStart+1); 
			m_elementOfSeaSurfaceWithInactiveElements.insert(newElemIDStart+2);
			m_elementOfSeaSurfaceWithInactiveElements.insert(newElemIDStart+3);
		}
	}

}

// Calculate z-coordinate after including topography/bathymetry
double MeshGenerationDHexa::calcZCoordAfterIncludingTopography( const double zMin, const double zMax, const double zCoordSurfCur,
	const double zCoordSurfTarget, const CommonParameters::XYZ& coordCur, const bool isSea ) const{

	const double zCur = coordCur.Z;
	const double zDiffAtSurf = zCoordSurfTarget - zCoordSurfCur;

	if( isSea ){
		if( zCur >= zCoordSurfCur ){
			const double factor = (zMax - zCur) / std::max(zMax - zCoordSurfCur, 1e-6);
			return zCur + zDiffAtSurf * factor;
		}else if( zCur >= 0.0 ){
			const double zMinMod = std::max(zMin, 0.0);
			const double factor = (zCur - zMinMod) / std::max(zCoordSurfCur - zMinMod, 1e-6);
			return zCur + zDiffAtSurf * factor;
		}
	}else{
		if( zCur >= zCoordSurfCur ){
			const double factor = (zMax - zCur) / std::max(zMax - zCoordSurfCur, 1e-6);
			return zCur + zDiffAtSurf * factor;
		}else{
			const double factor = (zCur - zMin) / std::max(zCoordSurfCur - zMin, 1e-6);
			return zCur + zDiffAtSurf * factor;
		}
	}

	return zCur;

}

// Check whether all children have flat surface
bool MeshGenerationDHexa::doesAllChildrenHaveFlatSurface( const int elemIndex ) const{

	const ElementInfo& elemInfo = m_elemInfo[elemIndex];
	const double eps = 1.0e-6;
	for( int i = 0; i < 4; ++i ){
		const int nodeAtTop = elemInfo.nodes[i];
		const double depthTop = m_nodeCoordinates[nodeAtTop].Z;
		if( fabs(depthTop) > eps ){
			// The following procedure is performed only when the surface is parallel to the horizontal plane of z = 0km
			return false;
		}
	}

	for( int i = 0; i < 4; ++i ){
		const int elemIndexChild = elemInfo.childs[i];
		if( elemIndexChild < 0 ){
			continue;
		}
		if( !doesAllChildrenHaveFlatSurface(elemIndexChild) ){
			return false;
		}
	}

	return true;

}

// Insert to vertical nodes array
void MeshGenerationDHexa::insertToVerticalNodesArray( const int nodeID, std::vector< std::set<int> >& verticalNodesArray,
	std::vector< CommonParameters::XY >& horizontalCoordsArray ) const{

	const double eps = 1.0e-6;
	const CommonParameters::XY coord = { m_nodeCoordinates[nodeID].X, m_nodeCoordinates[nodeID].Y };
	bool found(false);
	int iArray(0);
	for( std::vector< std::set<int> >::iterator itrVEA = verticalNodesArray.begin();
		itrVEA != verticalNodesArray.end(); ++itrVEA, ++iArray ){
		const CommonParameters::XY coordArray = horizontalCoordsArray[iArray];
		if( fabs( coord.X - coordArray.X ) < eps && fabs( coord.Y - coordArray.Y ) < eps ){
			itrVEA->insert(nodeID);
			found = true;
			break;
		}
	}
	if( !found ){
		std::set<int> tmp;
		tmp.insert(nodeID);
		verticalNodesArray.push_back(tmp);
		horizontalCoordsArray.push_back(coord);
	}

}

// Check whether a node is connected to the element with lower level
bool MeshGenerationDHexa::isConnectedToLowerLevelElement( const std::set<int>& nodeIDs, const int level, const std::multimap<int, int>& nodeToElem ) const{

	for( std::set<int>::const_iterator itr = nodeIDs.begin(); itr != nodeIDs.end(); ++itr ){
		if( isConnectedToLowerLevelElement(*itr, level, nodeToElem) ){
			return true;
		}
	}

	return false;

}

// Check whether a node is connected to the element with lower level except top and bottom
bool MeshGenerationDHexa::isConnectedToLowerLevelElementExceptTopAndBottom( const std::set<int>& nodeIDs, const int level, const double zMin, const double zMax, 
			const std::multimap<int, int>& nodeToElem ) const{

	std::vector<double> zCoords;
	for( std::set<int>::const_iterator itr = nodeIDs.begin(); itr != nodeIDs.end(); ++itr ){
		if( isConnectedToLowerLevelElement(*itr, level, nodeToElem) ){
			const double zCoord = m_nodeCoordinates[*itr].Z;
			zCoords.push_back(zCoord);
		}
	}

	const double eps = 1.0e-6;
	if( static_cast<int>( zCoords.size() ) > 2 ){
		return true;
	}
	if( static_cast<int>( zCoords.size() ) == 2 ){
		std::sort(zCoords.begin(), zCoords.end());
		if( fabs( zCoords[0] - zMin ) < eps && fabs( zCoords[1] - zMax ) < eps ){
			return false;
		}else{
			return true;
		}
	}
	return false;

}

// Check whether a node is not connected to eight elements with the same level except top and bottom
bool MeshGenerationDHexa::isNotConnectedToEightSameLevelElementsExceptTopAndBottom( const std::set<int>& nodeIDs, const int level, const double zMin, const double zMax, 
	const std::multimap<int, int>& nodeToElem ) const{

	const double eps = 1.0e-6;
	const double xMin = m_CoordinatesXInit[0];
	const double xMax = m_CoordinatesXInit[m_numXInit];
	const double yMin = m_CoordinatesYInit[0];
	const double yMax = m_CoordinatesYInit[m_numYInit];
	std::vector<double> zCoords;
	for( std::set<int>::const_iterator itr = nodeIDs.begin(); itr != nodeIDs.end(); ++itr ){
		const int nodeID = *itr;
		const int numElemSameLevel = getNumOflElementsConnected(nodeID, level, nodeToElem);
		if (fabs(m_nodeCoordinates[nodeID].X - xMin) < eps || fabs(m_nodeCoordinates[nodeID].X - xMax) < eps ||
			fabs(m_nodeCoordinates[nodeID].Y - yMin) < eps || fabs(m_nodeCoordinates[nodeID].Y - yMax) < eps){
			// Node locates at a side boundary 
			if (numElemSameLevel != 4) {
				const double zCoord = m_nodeCoordinates[nodeID].Z;
				zCoords.push_back(zCoord);
			}
		}else{
			if (numElemSameLevel != 8) {
				const double zCoord = m_nodeCoordinates[nodeID].Z;
				zCoords.push_back(zCoord);
			}
		}
	}

	if( static_cast<int>( zCoords.size() ) > 2 ){
		return true;
	}
	if( static_cast<int>( zCoords.size() ) == 2 ){
		std::sort(zCoords.begin(), zCoords.end());
		if( fabs( zCoords[0] - zMin ) < eps && fabs( zCoords[1] - zMax ) < eps ){
			return false;
		}else{
			return true;
		}
	}
	return false;
}

// Check whether a node is connected to the element with lower level
bool MeshGenerationDHexa::isConnectedToLowerLevelElement( const int nodeID, const int level, const std::multimap<int, int>& nodeToElem ) const{

	typedef std::multimap<int, int> MMII;
	const std::pair<MMII::const_iterator, MMII::const_iterator> ret = nodeToElem.equal_range(nodeID);
	for( MMII::const_iterator itrMMII = ret.first; itrMMII != ret.second; ++itrMMII ){
		const int elemID = itrMMII->second;
		if( !isActive(elemID) ){
			continue;
		}
		if( m_elemInfo[elemID].level < level ){
			return true;
		}			
	}

	return false;

}

// Get number of element connected to the node
int MeshGenerationDHexa::getNumOflElementsConnected( const int nodeID, const int level, const std::multimap<int, int>& nodeToElem ) const{

	int counter(0);
	typedef std::multimap<int, int> MMII;
	const std::pair<MMII::const_iterator, MMII::const_iterator> ret = nodeToElem.equal_range(nodeID);
	for( MMII::const_iterator itrMMII = ret.first; itrMMII != ret.second; ++itrMMII ){
		const int elemID = itrMMII->second;
		if( m_elemInfo[elemID].level == level ){
			++counter;
		}			
	}

	return counter;

}

// Get resistivity from element
double MeshGenerationDHexa::getResistivityFromElement( const ElementInfo& info ) const{

	const int paramCellIDOrg = info.parameterCell;
	std::map<int, double>::const_iterator itr = m_parameterCellToResistivity.find(paramCellIDOrg);
	if( itr == m_parameterCellToResistivity.end() ){
		std::cerr << "Element index is wrong !!" << std::endl;
	}
	return itr->second;

}

// Auxiliary function for including topography
void MeshGenerationDHexa::includeTopographyAux( const int iLevel, const int maxLevel, const std::vector<CommonParameters::XYZ>& nodeCoordOrg, const std::multimap<int,int>& nodeToElem ){

	std::vector< std::set<int> > verticalNodesArray;
	std::vector< CommonParameters::XY > horizontalCoordsArray;
	// Earth surface
	for( std::set<int>::const_iterator itr = m_elementOfEarthSurfaceWithInactiveElements.begin();
		itr != m_elementOfEarthSurfaceWithInactiveElements.end(); ++itr ){
		const int elemID = *itr;
		const ElementInfo& info = m_elemInfo[elemID];
		if( info.level != iLevel ){
			continue;
		}
		//if( !m_includeSea && calcAverageZCoord(elemID) > 0.0 ){
		//	// If the average sea depth is positive, no topography is incorpolated
		//	// Subsequently, sea resistivity is assigned to the elements
		//	continue;
		//}
		for( int i = 0; i < 8; ++i ){
			const int nodeID = info.nodes[i];
			insertToVerticalNodesArray( nodeID, verticalNodesArray, horizontalCoordsArray );
		}
		// Search upward
		int elemIDNext = info.neib[4][0];
		while( elemIDNext >= 0 && m_elemInfo[elemIDNext].level >= iLevel ){
			if( m_elemInfo[elemIDNext].level > iLevel ){
				elemIDNext = m_elemInfo[elemIDNext].parent;
			}
			assert( m_elemInfo[elemIDNext].level == iLevel );
			for( int i = 0; i < 8; ++i ){
				const int nodeID = m_elemInfo[elemIDNext].nodes[i];
				insertToVerticalNodesArray( nodeID, verticalNodesArray, horizontalCoordsArray );
			}
			elemIDNext = m_elemInfo[elemIDNext].neib[4][0];
		}
		// Search downward
		elemIDNext = info.neib[5][0];
		while( elemIDNext >= 0 && m_elemInfo[elemIDNext].level >= iLevel ){
			if( m_elemInfo[elemIDNext].level > iLevel ){
				elemIDNext = m_elemInfo[elemIDNext].parent;
			}
			assert( m_elemInfo[elemIDNext].level == iLevel );
			for( int i = 0; i < 8; ++i ){
				const int nodeID = m_elemInfo[elemIDNext].nodes[i];
				insertToVerticalNodesArray( nodeID, verticalNodesArray, horizontalCoordsArray );
			}
			elemIDNext = m_elemInfo[elemIDNext].neib[5][0];
		}
	}
	for( std::vector< std::set<int> >::const_iterator itrArray = verticalNodesArray.begin(); itrArray != verticalNodesArray.end(); ++itrArray ){
		double zMin = 1.0e20;
		double zMax = -1.0e20;
		for( std::set<int>::const_iterator itr = itrArray->begin(); itr != itrArray->end(); ++itr ){
			const int nodeID = *itr;
			const CommonParameters::XYZ coord = m_nodeCoordinates[nodeID];
			if( coord.Z < zMin ){
				zMin = coord.Z;
			}
			if( coord.Z > zMax ){
				zMax = coord.Z;
			}
		}
		if( iLevel > 0 && isNotConnectedToEightSameLevelElementsExceptTopAndBottom( *itrArray, iLevel, zMin, zMax, nodeToElem ) ){
			continue;
		}
		bool found(false);
		const double eps = 1.0e-6;
		for( std::set<int>::const_iterator itr = itrArray->begin(); itr != itrArray->end(); ++itr ){
			const int nodeID = *itr;
			const CommonParameters::XYZ coordOrg = nodeCoordOrg[nodeID];
			const CommonParameters::XYZ coordCur = m_nodeCoordinates[nodeID];
			const double zCoordSurfOrg = m_includeSea ? m_seaDepth : 0.0;
			if( fabs( coordOrg.Z - zCoordSurfOrg ) < eps ){
				found = true;
				if( iLevel > 0 && ( fabs( coordCur.Z - zMin ) < eps || fabs( coordCur.Z - zMax ) < eps ) ){
					// Skip if the surface locates at the end of the array
					break;
				}
				const CommonParameters::XY coordXY = {coordOrg.X, coordOrg.Y};
				if( m_includeSea ){
					const double zCoordSurfTarget = std::max( m_topographyData->interpolateZCoord(coordXY), m_minimumSeaDepth );
					if (zCoordSurfTarget < zMin || zCoordSurfTarget > zMax) {
						break;
					}
					for( std::set<int>::const_iterator itrAux = itrArray->begin(); itrAux != itrArray->end(); ++itrAux ){
						const int nodeIDAux = *itrAux;
						m_nodeCoordinates[nodeIDAux].Z = calcZCoordAfterIncludingTopography(zMin, zMax, coordCur.Z, zCoordSurfTarget, m_nodeCoordinates[nodeIDAux], true);
					}
				}else{
					double zCoordSurfTarget = m_topographyData->interpolateZCoord(coordXY);
					if( zCoordSurfTarget > 0.0 ){
						// The upper limit of the depth is zero
						zCoordSurfTarget = 0.0;
					}
					if (zCoordSurfTarget < zMin || zCoordSurfTarget > zMax) {
						break;
					}
					for( std::set<int>::const_iterator itrAux = itrArray->begin(); itrAux != itrArray->end(); ++itrAux ){
						const int nodeIDAux = *itrAux;
						m_nodeCoordinates[nodeIDAux].Z = calcZCoordAfterIncludingTopography(zMin, zMax, coordCur.Z, zCoordSurfTarget, m_nodeCoordinates[nodeIDAux], false);
					}
				}
				break;
			}
		}
		if(!found){
			std::cerr << "Surface node is not found !!" << std::endl;
			exit(1);
		}
	}

	for( int i = iLevel; i < maxLevel; ++i ){
		std::set<int> dummy;
		linearInterpolationOfCoordinatesForChildren(i, false, dummy);
	}

}

// Auxiliary function for changing the sea to land
void MeshGenerationDHexa::includeTopographyAux2( const int iLevel, const int maxLevel, const std::vector<CommonParameters::XYZ>& nodeCoordOrg, const std::multimap<int,int>& nodeToElem ){

	std::vector< std::set<int> > verticalNodesArray;
	std::vector< CommonParameters::XY > horizontalCoordsArray;
	std::set<int> elementsLand;
	// Sea surface
	for( std::set<int>::const_iterator itr = m_elementOfSeaSurfaceWithInactiveElements.begin();
		itr != m_elementOfSeaSurfaceWithInactiveElements.end(); ++itr ){
		const int elemID = *itr;
		const ElementInfo& info = m_elemInfo[elemID];
		if( info.level != iLevel ){
			continue;
		}
		if( !doesAllChildrenBelongToLandOrAir(elemID) ){
			// Go to the following procedure only if all children belong to land or air
			continue;
		}
		elementsLand.insert(elemID);
		for( int i = 0; i < 8; ++i ){
			const int nodeID = info.nodes[i];
			insertToVerticalNodesArray( nodeID, verticalNodesArray, horizontalCoordsArray );
		}
		// Search upward
		int elemIDNext = info.neib[4][0];
		while( elemIDNext >= 0 && m_elemInfo[elemIDNext].level >= iLevel ){
			if( m_elemInfo[elemIDNext].level > iLevel ){
				elemIDNext = m_elemInfo[elemIDNext].parent;
			}
			assert( m_elemInfo[elemIDNext].level == iLevel );
			elementsLand.insert(elemIDNext);
			for( int i = 0; i < 8; ++i ){
				const int nodeID = m_elemInfo[elemIDNext].nodes[i];
				insertToVerticalNodesArray( nodeID, verticalNodesArray, horizontalCoordsArray );
			}
			elemIDNext = m_elemInfo[elemIDNext].neib[4][0];
		}
		// Search downward
		elemIDNext = info.neib[5][0];
		while( elemIDNext >= 0 && m_elemInfo[elemIDNext].level >= iLevel ){
			if( m_elemInfo[elemIDNext].level > iLevel ){
				elemIDNext = m_elemInfo[elemIDNext].parent;
			}
			assert( m_elemInfo[elemIDNext].level == iLevel );
			elementsLand.insert(elemIDNext);
			for( int i = 0; i < 8; ++i ){
				const int nodeID = m_elemInfo[elemIDNext].nodes[i];
				insertToVerticalNodesArray( nodeID, verticalNodesArray, horizontalCoordsArray );
			}
			elemIDNext = m_elemInfo[elemIDNext].neib[5][0];
		}
	}
	for( std::vector< std::set<int> >::const_iterator itrArray = verticalNodesArray.begin(); itrArray != verticalNodesArray.end(); ++itrArray ){
		bool goAhead(true);
		for( std::set<int>::const_iterator itr = itrArray->begin(); itr != itrArray->end(); ++itr ){
			const int nodeID = *itr;
			typedef std::multimap<int, int> MMII;
			const std::pair<MMII::const_iterator, MMII::const_iterator> ret = nodeToElem.equal_range(nodeID);
			for( MMII::const_iterator itrMMII = ret.first; itrMMII != ret.second; ++itrMMII ){
				const int elemID = itrMMII->second;
				if( m_elemInfo[elemID].level == iLevel && !doesAllChildrenBelongToLandOrAir(elemID) ){
#ifdef _DEBUG_WRITE
					const CommonParameters::XYZ coordOrg = nodeCoordOrg[nodeID];			
#endif
					// Exclude the nodes which connect to the sea element
					goAhead = false;
					break;
				}			
			}
			if(!goAhead){
				break;
			}
		}
		if(!goAhead){
			continue;
		}

		double zMin = 1.0e20;
		double zMax = -1.0e20;
		for( std::set<int>::const_iterator itr = itrArray->begin(); itr != itrArray->end(); ++itr ){
			const int nodeID = *itr;
			const CommonParameters::XYZ coord = m_nodeCoordinates[nodeID];
			if( coord.Z < zMin ){
				zMin = coord.Z;
			}
			if( coord.Z > zMax ){
				zMax = coord.Z;
			}
		}
		if( iLevel > 0 && isNotConnectedToEightSameLevelElementsExceptTopAndBottom( *itrArray, iLevel, zMin, zMax, nodeToElem ) ){
			continue;
		}
		bool found(false);
		const double eps = 1.0e-6;
		for( std::set<int>::const_iterator itr = itrArray->begin(); itr != itrArray->end(); ++itr ){
			const int nodeID = *itr;
			const CommonParameters::XYZ coordOrg = nodeCoordOrg[nodeID];
			const CommonParameters::XYZ coordCur = m_nodeCoordinates[nodeID];
			const double zCoordSurfOrg = 0.0;
			if( fabs( coordOrg.Z - zCoordSurfOrg ) < eps ){
				found = true;
				if( iLevel > 0 && ( fabs( coordCur.Z - zMin ) < eps || fabs( coordCur.Z - zMax ) < eps ) ){
					// Skip if the surface locates at the end of the array
					break;
				}
				const CommonParameters::XY coordXY = {coordOrg.X, coordOrg.Y};
				const double zCoordSurfTarget = m_topographyData->interpolateZCoord(coordXY);
				if ( zCoordSurfTarget < zMin || zCoordSurfTarget > zMax ) {
					break;
				}
				for( std::set<int>::const_iterator itrAux = itrArray->begin(); itrAux != itrArray->end(); ++itrAux ){
					const int nodeIDAux = *itrAux;
					m_nodeCoordinates[nodeIDAux].Z = calcZCoordAfterIncludingTopography(zMin, zMax, coordCur.Z, zCoordSurfTarget, m_nodeCoordinates[nodeIDAux], false);
				}
				break;
			}
		}
		if(!found){
			std::cerr << "Surface node is not found !!" << std::endl;
			exit(1);
		}
	}

	for( int i = iLevel; i < maxLevel; ++i ){
		linearInterpolationOfCoordinatesForChildren(i, true, elementsLand);
	}

}

//// Auxiliary function for including topography for land area changed from the sea
//// @todo Under construction
//void MeshGenerationDHexa::includeTopographyForLandAreaChangedFromSeaAux( const int iLevel, const std::vector<CommonParameters::XYZ>& nodeCoordOrg,
//	std::set<int>& surfElemsLand, std::set<int>& surfElemSea, const std::multimap<int,int>& nodeToElem ){
//
//	std::vector< std::set<int> > verticalNodesArray;
//	std::vector< CommonParameters::XY > horizontalCoordsArray;
//	std::set<int> nodesLandArea;
//	std::set<int> nodesSeaArea;
//	for( std::set<int>::const_iterator itr = m_elementOfSeaSurface.begin(); itr != m_elementOfSeaSurface.end(); ++itr ){
//		// Because m_elementOfSeaSurfaceWithInactiveElements was not reconstructed, it includes elements of land area changed from the sea
//		const int elemID = *itr;
//		const ElementInfo& info = m_elemInfo[elemID];
//		if( surfElemsLand.find(elemID) != surfElemsLand.end() ){
//			for( int i = 0; i < 8; ++i ){
//				const int nodeID = info.nodes[i];
//				nodesLandArea.insert(nodeID);
//			}
//		}else{
//			for( int i = 0; i < 8; ++i ){
//				const int nodeID = info.nodes[i];
//				nodesSeaArea.insert(nodeID);
//			}
//		}
//	}
//
//	for( std::set<int>::const_iterator itr = m_elementOfSeaSurfaceWithInactiveElements.begin();
//		itr != m_elementOfSeaSurfaceWithInactiveElements.end(); ++itr ){
//		// Because m_elementOfSeaSurfaceWithInactiveElements was not reconstructed, it includes elements of land area changed from the sea
//		const int elemID = *itr;
//		const ElementInfo& info = m_elemInfo[elemID];
//		if( info.level != iLevel ){
//			continue;
//		}
//		if( doChildrenHaveSea(elemID) ){
//			continue;
//		}
//		for( int i = 0; i < 8; ++i ){
//			const int nodeID = info.nodes[i];
//			insertToVerticalNodesArray( nodeID, verticalNodesArray, horizontalCoordsArray );
//		}
//		// Search upward
//		int elemIDNext = info.neib[4][0];
//		while( elemIDNext >= 0 && m_elemInfo[elemIDNext].level >= iLevel ){
//			if( m_elemInfo[elemIDNext].level > iLevel ){
//				elemIDNext = m_elemInfo[elemIDNext].parent;
//			}
//			assert( m_elemInfo[elemIDNext].level == iLevel );
//			for( int i = 0; i < 8; ++i ){
//				const int nodeID = m_elemInfo[elemIDNext].nodes[i];
//				insertToVerticalNodesArray( nodeID, verticalNodesArray, horizontalCoordsArray );
//			}
//			elemIDNext = m_elemInfo[elemIDNext].neib[4][0];
//		}
//		// Search downward
//		elemIDNext = info.neib[5][0];
//		while( elemIDNext >= 0 && m_elemInfo[elemIDNext].level >= iLevel ){
//			if( m_elemInfo[elemIDNext].level > iLevel ){
//				elemIDNext = m_elemInfo[elemIDNext].parent;
//			}
//			assert( m_elemInfo[elemIDNext].level == iLevel );
//			for( int i = 0; i < 8; ++i ){
//				const int nodeID = m_elemInfo[elemIDNext].nodes[i];
//				insertToVerticalNodesArray( nodeID, verticalNodesArray, horizontalCoordsArray );
//			}
//			elemIDNext = m_elemInfo[elemIDNext].neib[5][0];
//		}
//	}
//	for( std::vector< std::set<int> >::const_iterator itrArray = verticalNodesArray.begin(); itrArray != verticalNodesArray.end(); ++itrArray ){
//		double zMin = 1.0e20;
//		double zMax = -1.0e20;
//		for( std::set<int>::const_iterator itr = itrArray->begin(); itr != itrArray->end(); ++itr ){
//			const int nodeID = *itr;
//			const CommonParameters::XYZ coord = m_nodeCoordinates[nodeID];
//			if( coord.Z < zMin ){
//				zMin = coord.Z;
//			}
//			if( coord.Z > zMax ){
//				zMax = coord.Z;
//			}
//		}
//		if( iLevel > 0 && isNotConnectedToEightSameLevelElementsExceptTopAndBottom( *itrArray, iLevel, zMin, zMax, nodeToElem ) ){
//			continue;
//		}
//		const double eps = 1.0e-6;
//		for( std::set<int>::const_iterator itr = itrArray->begin(); itr != itrArray->end(); ++itr ){
//			const int nodeID = *itr;
//			const CommonParameters::XYZ coordOrg = nodeCoordOrg[nodeID];
//			const CommonParameters::XYZ coordCur = m_nodeCoordinates[nodeID];
//			const double zCoordSurfOrg = 0.0;
//			if( fabs( coordOrg.Z - zCoordSurfOrg ) < eps ){
//				const bool isLand = nodesLandArea.find(nodeID) != nodesLandArea.end();
//				const bool isSea = nodesSeaArea.find(nodeID) != nodesSeaArea.end();
//				if( isLand && !isSea ){
//					// Skip land-sea boundary
//					const CommonParameters::XY coordXY = {coordOrg.X, coordOrg.Y};
//					const double zCoordSurfTarget = m_topographyData->interpolateZCoord(coordXY);
//					for( std::set<int>::const_iterator itrAux = itrArray->begin(); itrAux != itrArray->end(); ++itrAux ){
//						const int nodeIDAux = *itrAux;
//						m_nodeCoordinates[nodeIDAux].Z = calcZCoordAfterIncludingTopography(zMin, zMax, coordCur.Z, zCoordSurfTarget, m_nodeCoordinates[nodeIDAux], false);
//					}
//				}
//				break;
//			}
//		}
//	}
//
//	linearInterpolationOfCoordinatesForChildren(iLevel);
//
//}

// Select elements to be changed to land from the sea
void MeshGenerationDHexa::selectElementsToBeChangedToLandFromSea( std::set<int>& elemsSeaToLand ){

	const double eps = 1.0e-6;
	
	for( std::set<int>::const_iterator itr = m_elementOfEarthSurface.begin(); itr != m_elementOfEarthSurface.end(); ++itr ){
		const int elemIndex = *itr;
		if( elemIndex < 0 ){
			continue;
		}
		const ElementInfo& info = m_elemInfo[elemIndex];
		bool searchUpward(true);
		for( int i = 0; i < 4; ++i ){
			const int nodeAtTop = info.nodes[i];
			const double depthTop = m_nodeCoordinates[nodeAtTop].Z;
			if( fabs(depthTop - m_minimumSeaDepth) > eps ){
				// The following procedure is performed only when all top depth is constrained to the minimum depth
				searchUpward = false;
			}
		}
		if( !searchUpward ){
			continue;
		}
		// Search upward
		int level = info.level;
		int elemIDNext = info.neib[4][0];
		while( m_elemInfo[elemIDNext].level > level ){
			elemIDNext = m_elemInfo[elemIDNext].parent;
		}
		while( elemIDNext >= 0 ){
			if( m_elemInfo[elemIDNext].type == AIR ){
				break;
			}
			const double seaDepth = calcAverageZCoord(elemIDNext);
			bool isLand(true);
			for( int i = 0; i < 4; ++i ){
				const int nodeAtTop = m_elemInfo[elemIDNext].nodes[i];
				const double depthTop = m_nodeCoordinates[nodeAtTop].Z;
				if( depthTop < seaDepth ){
					isLand = false;
				}
			}
			if( !isLand ){
				break;
			}
			addChildrenToLandCandidates(elemIDNext, elemsSeaToLand);
			elemIDNext = m_elemInfo[elemIDNext].neib[4][0];
			while( m_elemInfo[elemIDNext].level > level ){
				elemIDNext = m_elemInfo[elemIDNext].parent;
			}
		}
	}

	bool found(true);
	while(found){
		found = false;
		std::set<int> elemsSeaToLandOrg = elemsSeaToLand;
		for( std::set<int>::iterator itr = elemsSeaToLandOrg.begin(); itr != elemsSeaToLandOrg.end(); ++itr ){
			const int elemIndex = *itr;
			const int level = m_elemInfo[elemIndex].level;
			for( int i = 0 ; i < 4; ++i){
				// Search lower elements
				const int elemNeib = m_elemInfo[elemIndex].neib[5][i];
				if( elemNeib < 0 ){
					continue;
				}
				const int levelNeib = m_elemInfo[elemNeib].level;
				if( m_elementOfEarthSurface.find(elemNeib) != m_elementOfEarthSurface.end() ){
					// Skip if lower element belongs to surface elements
					continue;
				}
				if( elemsSeaToLand.find(elemNeib) == elemsSeaToLand.end() ){
					// Remove if there is a sea element below the element 
					elemsSeaToLand.erase(elemIndex);
					found = true;
				}
			}
		}
	}

}

// Select elements to be changed to seat from land
void MeshGenerationDHexa::selectElementsToBeChangedToSeaFromLand( std::map<int, double>& elemsLandToSea, std::map<int, double>& elemEarthSurfToSeaDepth  ) const{

	if( m_includeSea ){
		// This function is effective only when the smooth sea floor is NOT incorporated
		return;
	}

	for( std::set<int>::const_iterator itr = m_elementOfEarthSurfaceWithInactiveElements.begin(); itr != m_elementOfEarthSurfaceWithInactiveElements.end(); ++itr ){
		int elemIndex = *itr;
		if( elemIndex < 0 ){
			continue;
		}
		if( !doesAllChildrenHaveFlatSurface(elemIndex) ){
			continue;
		}
		const double depthAvg = calcAverageZCoord(elemIndex);
		if( depthAvg < CommonParameters::EPS ){
			continue;
		}
		elemEarthSurfToSeaDepth.insert( std::make_pair(elemIndex, depthAvg) );
		const ElementInfo& info = m_elemInfo[elemIndex];
		selectElementsToBeChangedToSeaFromLandAux(info.level, depthAvg, elemIndex, elemsLandToSea);
	}

}

// Auxiliary function for selecting elements to be changed to seat from land
void MeshGenerationDHexa::selectElementsToBeChangedToSeaFromLandAux( const int level, const double depthAvg, int elemIndex, std::map<int, double>& elemsLandToSea ) const{

	while( elemIndex >= 0 ){
		assert( checkWhetherElementHasFlatTopAndBottomSurfaces(m_elemInfo[elemIndex]) );
		assert( m_elemInfo[elemIndex].level == level );
		if( m_elemInfo[elemIndex].type != LAND ){
			break;
		}
		const int nodeAtTop = m_elemInfo[elemIndex].nodes[0];
		const int nodeAtBot = m_elemInfo[elemIndex].nodes[7];
		const double depthTop = m_nodeCoordinates[nodeAtTop].Z;
		const double depthBot = m_nodeCoordinates[nodeAtBot].Z;
		if( depthTop > depthAvg ){
			// The element is below the sea floor
			break;
		}
		if( depthAvg > depthBot ){
			// The whole part of the element is the sea
			elemsLandToSea.insert( std::make_pair( elemIndex, m_seaResistivity ) );
		}else{
			const double elemHeight = fabs( depthBot - depthTop );
			const double d1 = ( depthAvg - depthTop ) / elemHeight;
			const double d2 = ( depthBot - depthAvg ) / elemHeight;
			if( d1 < 0.01 ){
				break;
			}
			const double conductivitySea = 1.0 / m_seaResistivity;
			const double conductivityOrg = 1.0 / getResistivityFromElement(m_elemInfo[elemIndex]);
			const double conductivityNew = conductivitySea * d1 + conductivityOrg * d2;
			elemsLandToSea.insert( std::make_pair( elemIndex, 1.0 / conductivityNew ) );
			break;
		}
		// Search downward
		int elemIDNext = m_elemInfo[elemIndex].neib[5][0];
		if( m_elemInfo[elemIDNext].level < level ){
			break;
		}else if( m_elemInfo[elemIDNext].level > level ){
			const int levelNext = m_elemInfo[elemIDNext].level;
			for( int i = 0; i < 4; ++i ){
				const int elemIDNext = m_elemInfo[elemIndex].neib[5][i];
				selectElementsToBeChangedToSeaFromLandAux(levelNext, depthAvg, elemIDNext, elemsLandToSea);
			}
			elemIDNext = m_elemInfo[elemIDNext].parent;
		}
		elemIndex = elemIDNext;
	}

	return;

}

// Calculated average Z coordinate in an element
double MeshGenerationDHexa::calcAverageZCoord( const int elemIndex ) const{

	const int node0 = m_elemInfo[elemIndex].nodes[0];
	const int node2 = m_elemInfo[elemIndex].nodes[2];
	const double xMin = m_nodeCoordinates[node0].X;
	const double xMax = m_nodeCoordinates[node2].X;
	const double yMin = m_nodeCoordinates[node0].Y;
	const double yMax = m_nodeCoordinates[node2].Y;
	return m_topographyData->calcAverageZCoord(xMin, xMax, yMin, yMax);

}

// Change sea to land
bool MeshGenerationDHexa::changeSeaToLand( const int elemIndex, const std::set<int>& elemsSeaToLand, const int paramCellIDNew ){
	
	ElementInfo& info = m_elemInfo[elemIndex];
	if( elemsSeaToLand.find(elemIndex) != elemsSeaToLand.end() ){		
		const int paramCellIDOrg = info.parameterCell;
		const double resistivity = m_parameterCellToResistivity[paramCellIDOrg];
		const int fixFlag = m_parameterCellToFixFlag[paramCellIDOrg];
		info.parameterCell = paramCellIDNew;
		std::map<int, double> dummy;
		giveSameParamCellIDToChildren(elemIndex, paramCellIDNew, dummy);
		info.type = LAND;
		giveSameTypeToChildren(elemIndex, LAND);
		return true;
	}else{
		bool found(false);
		for( int i = 0; i < 8; ++i ){
			const int elemIndexChild = info.childs[i];
			if( elemIndexChild < 0 ){
				continue;
			}
			if( changeSeaToLand( elemIndexChild, elemsSeaToLand, paramCellIDNew ) ){
				found = true;
			}
		}
		return found;
	}

}

// Output Earth surface depth by vtk
void MeshGenerationDHexa::outputEarthSurfaceDepthByVtk( const std::map<int, double>& elemEarthSurfToSeaDepth ) const{

	std::map<int,double> surfaceNodesToDepth;
	for( std::set<int>::const_iterator itr = m_elementOfEarthSurface.begin(); itr != m_elementOfEarthSurface.end(); ++itr ){
		const int elemIndex = *itr;
		const ElementInfo& info = m_elemInfo[elemIndex];
		for( int i = 0; i < 4; ++i ){
			const int nodeID = info.nodes[i];
			double depth(0.0);
			depth = m_nodeCoordinates[nodeID].Z;
			surfaceNodesToDepth.insert( std::make_pair(nodeID, depth) );
		}
	}

	std::map<int,int> surfaceNodesToIndex;
	int index(0);
	for( std::map<int,double>::const_iterator itr = surfaceNodesToDepth.begin(); itr != surfaceNodesToDepth.end(); ++itr, ++index ){
		surfaceNodesToIndex.insert( std::make_pair( itr->first, index ) );
	}

	const std::string fileName = "depth.vtk";
	std::ofstream ofsVTK( fileName.c_str() );
	if( !ofsVTK ) {
		std::cerr << "Cannot open file " << fileName.c_str() << std::endl;
		exit(1);
	}

	ofsVTK << "# vtk DataFile Version 2.0" << std::endl;
	ofsVTK << "SurfaceTriangles" << std::endl;
	ofsVTK << "ASCII" << std::endl;
	ofsVTK << "DATASET UNSTRUCTURED_GRID" << std::endl;
	
	ofsVTK.precision(6);
	ofsVTK << std::fixed;

	const int numNodes = static_cast<int>( surfaceNodesToIndex.size() );
	ofsVTK << "POINTS " << numNodes<< " float" << std::endl;
	for( std::map<int,int>::const_iterator itr = surfaceNodesToIndex.begin(); itr != surfaceNodesToIndex.end(); ++itr ){
		const CommonParameters::XYZ coord = m_nodeCoordinates[itr->first];
		ofsVTK << std::setw(15) << std::scientific << coord.X
			   << std::setw(15) << std::scientific << coord.Y
			   << std::setw(15) << std::scientific << coord.Z << std::endl;
	}

	const int numElems = static_cast<int>( m_elementOfEarthSurface.size() );
	ofsVTK << "CELLS " << numElems << " " << numElems * 5 << std::endl;
	for( std::set<int>::const_iterator itr = m_elementOfEarthSurface.begin(); itr != m_elementOfEarthSurface.end(); ++itr ){
		const int elemIndex = *itr;
		const ElementInfo& info = m_elemInfo[elemIndex];
		ofsVTK << std::setw(10) << 4;
		for( int i = 0; i < 4; ++i ){
			const int nodeID = info.nodes[i];
			const int nodeIndex = surfaceNodesToIndex[nodeID];
			ofsVTK << std::setw(10) << nodeIndex;
		}
		ofsVTK << std::endl;
	}

	ofsVTK << "CELL_TYPES " << numElems << std::endl;
	for( int i = 0; i < numElems; ++i ){
		ofsVTK << std::setw(10) << 9 << std::endl;// QUAD
	}

	ofsVTK << "CELL_DATA " << numElems << std::endl;
	ofsVTK << "SCALARS Resistivity[Ohm-m] float" <<  std::endl;
	ofsVTK << "LOOKUP_TABLE default" <<  std::endl;
	for( std::set<int>::const_iterator itr = m_elementOfEarthSurface.begin(); itr != m_elementOfEarthSurface.end(); ++itr ){
		const int elemIndex = *itr;
		const int paramCell = m_elemInfo[elemIndex].parameterCell;
		std::map<int, double>::const_iterator itrParamCell = m_parameterCellToResistivity.find(paramCell);
		if( itrParamCell == m_parameterCellToResistivity.end() ){
			std::cerr << "Parameter cell " << paramCell << " is not found in m_parameterCellToResistivity." << std::endl;
			exit(1);
		}
		ofsVTK << itrParamCell->second << std::endl;
	}

	if( !m_includeSea ){
		ofsVTK << "SCALARS SeaDepth float" <<  std::endl;
		ofsVTK << "LOOKUP_TABLE default" <<  std::endl;
		for( std::set<int>::const_iterator itr = m_elementOfEarthSurface.begin(); itr != m_elementOfEarthSurface.end(); ++itr ){
			const int elemIndex = *itr;
			std::map<int, double>::const_iterator itrElemEarthSurfToSeaDepth = elemEarthSurfToSeaDepth.find(elemIndex);
			if( itrElemEarthSurfToSeaDepth == elemEarthSurfToSeaDepth.end() ){
				ofsVTK << 0.0 << std::endl;
			}else{
				ofsVTK << itrElemEarthSurfToSeaDepth->second << std::endl;
			}
		}
	}

	ofsVTK << "SCALARS Level int" <<  std::endl;
	ofsVTK << "LOOKUP_TABLE default" <<  std::endl;
	for( std::set<int>::const_iterator itr = m_elementOfEarthSurface.begin(); itr != m_elementOfEarthSurface.end(); ++itr ){
		const int elemIndex = *itr;
		const int level = m_elemInfo[elemIndex].level;
		ofsVTK << level << std::endl;
	}

	ofsVTK << "POINT_DATA " << numNodes << std::endl;
	ofsVTK << "SCALARS Depth float" << std::endl;
	ofsVTK << "LOOKUP_TABLE default" << std::endl;
	for( std::map<int,double>::const_iterator itr = surfaceNodesToDepth.begin(); itr != surfaceNodesToDepth.end(); ++itr ){
		ofsVTK << std::setw(15) << std::scientific << static_cast<float>(itr->second) << std::endl;
	}

	ofsVTK << "SCALARS NodeID int" << std::endl;
	ofsVTK << "LOOKUP_TABLE default" << std::endl;
	for( std::map<int,int>::const_iterator itr = surfaceNodesToIndex.begin(); itr != surfaceNodesToIndex.end(); ++itr ){
		ofsVTK << std::setw(15) << itr->first << std::endl;
	}

	ofsVTK << "SCALARS NodeIndex int" << std::endl;
	ofsVTK << "LOOKUP_TABLE default" << std::endl;
	for( std::map<int,int>::const_iterator itr = surfaceNodesToIndex.begin(); itr != surfaceNodesToIndex.end(); ++itr ){
		ofsVTK << std::setw(15) << itr->second << std::endl;
	}

	ofsVTK.close();

}

// Calculate center coordinate
CommonParameters::XY MeshGenerationDHexa::calcCenterCoordinate( const int elemIndex ) const{

	const ElementInfo& elem = m_elemInfo[elemIndex];
	const int nodesIDs[4] = {
		elem.nodes[0],
		elem.nodes[1],
		elem.nodes[2],
		elem.nodes[3] };
	const CommonParameters::XYZ coords[4] = {
		m_nodeCoordinates[nodesIDs[0]],
		m_nodeCoordinates[nodesIDs[1]],
		m_nodeCoordinates[nodesIDs[2]],
		m_nodeCoordinates[nodesIDs[3]] };
	const CommonParameters::XY coordXY = {
		0.25 * ( coords[0].X + coords[1].X + coords[2].X + coords[3].X ),
		0.25 * ( coords[0].Y + coords[1].Y + coords[2].Y + coords[3].Y ) 
	};

	return coordXY;
}

// Calculate center coordinate of element face
CommonParameters::XYZ MeshGenerationDHexa::calculateCenterCoordOfElemFace( const std::vector<ElementInfo>::const_iterator& itrElem, const int iFace ) const{

	int nodeIndex[4] = {-1, -1, -1, -1};
	switch (iFace){
		case 0:
			nodeIndex[0] = itrElem->nodes[0];
			nodeIndex[1] = itrElem->nodes[3];
			nodeIndex[2] = itrElem->nodes[7];
			nodeIndex[3] = itrElem->nodes[4];
			break;
		case 1:
			nodeIndex[0] = itrElem->nodes[1];
			nodeIndex[1] = itrElem->nodes[2];
			nodeIndex[2] = itrElem->nodes[6];
			nodeIndex[3] = itrElem->nodes[5];
			break;
		case 2:
			nodeIndex[0] = itrElem->nodes[0];
			nodeIndex[1] = itrElem->nodes[1];
			nodeIndex[2] = itrElem->nodes[5];
			nodeIndex[3] = itrElem->nodes[4];
			break;
		case 3:
			nodeIndex[0] = itrElem->nodes[3];
			nodeIndex[1] = itrElem->nodes[2];
			nodeIndex[2] = itrElem->nodes[6];
			nodeIndex[3] = itrElem->nodes[7];
			break;
		case 4:
			nodeIndex[0] = itrElem->nodes[0];
			nodeIndex[1] = itrElem->nodes[1];
			nodeIndex[2] = itrElem->nodes[2];
			nodeIndex[3] = itrElem->nodes[3];
			break;
		case 5:
			nodeIndex[0] = itrElem->nodes[4];
			nodeIndex[1] = itrElem->nodes[5];
			nodeIndex[2] = itrElem->nodes[6];
			nodeIndex[3] = itrElem->nodes[7];
			break;
		default:
			std::cerr << "Face number is wrong : " << iFace << std::endl;
			break;
	}

	CommonParameters::XYZ centerCoord = {0.0, 0.0, 0.0};
	
	for( int i = 0; i < 4; ++i ){
		const CommonParameters::XYZ cornerCoord = m_nodeCoordinates[nodeIndex[i]];
		centerCoord.X += cornerCoord.X;
		centerCoord.Y += cornerCoord.Y;
		centerCoord.Z += cornerCoord.Z;
	}
	
	centerCoord.X *= 0.25;
	centerCoord.Y *= 0.25;
	centerCoord.Z *= 0.25;

	return centerCoord;

}

// Give same parameter cell ID to children
void MeshGenerationDHexa::giveSameParamCellIDToChildren( const int elemIndex, const int paramCellID, const std::map<int, double>& elemsLandToSea  ){

	ElementInfo& elemInfo = m_elemInfo[elemIndex];
	if( elemsLandToSea.find(elemIndex) == elemsLandToSea.end() ){
		// Skip when the element is converted to the sea
		elemInfo.parameterCell = paramCellID;
	}
	for( int i = 0; i < 8; ++i ){
		if( elemInfo.childs[i] < 0 ){
			continue;
		}
		giveSameParamCellIDToChildren(elemInfo.childs[i], paramCellID, elemsLandToSea);
	}

}

// Give same type to children
void MeshGenerationDHexa::giveSameTypeToChildren( const int elemIndex, const int type ){

	ElementInfo& elemInfo = m_elemInfo[elemIndex];
	elemInfo.type = type;
	for( int i = 0; i < 8; ++i ){
		if( elemInfo.childs[i] < 0 ){
			continue;
		}
		giveSameTypeToChildren(elemInfo.childs[i], type);
	}

}

// Check whether children have sea area
bool MeshGenerationDHexa::doChildrenHaveSea( const int elemIndex ) const{

	const ElementInfo& elemInfo = m_elemInfo[elemIndex];
	if( elemInfo.type == SEA ){
		return true;
	}
	bool haveSea(false);
	for( int i = 0; i < 8; ++i ){
		if( elemInfo.childs[i] < 0 ){
			continue;
		}
		if(doChildrenHaveSea(elemInfo.childs[i])){
			haveSea = true;
		}
	}
	return haveSea;

}

// Linear interpolation of coordinates for children
void MeshGenerationDHexa::linearInterpolationOfCoordinatesForChildren( const int level, const bool selectOnlyLand, const std::set<int>& elementsLand ){

	int elemIndex(0);
	for( std::vector<ElementInfo>::const_iterator itr = m_elemInfo.begin(); itr != m_elemInfo.end(); ++itr, ++elemIndex ){
		if( itr->level != level ){
			continue;
		}
		if( selectOnlyLand && elementsLand.find(elemIndex) == elementsLand.end() ){
			// The element does not belong to land
			continue;
		}
		double coordZsParent[8];
		for( int iNode = 0; iNode < 8; ++iNode ){
			const int nodeID = itr->nodes[iNode];
			coordZsParent[iNode] = m_nodeCoordinates[nodeID].Z;
		}
		if( getPartitioningType() == FULL_PARTITIONING ){
			// Interpolate coordinates of childs
			double coordZsChild[27];
			coordZsChild[0]  = coordZsParent[0];
			coordZsChild[2]  = coordZsParent[1];
			coordZsChild[6]  = coordZsParent[3];
			coordZsChild[8]  = coordZsParent[2];
			coordZsChild[18] = coordZsParent[4];
			coordZsChild[20] = coordZsParent[5];
			coordZsChild[24] = coordZsParent[7];
			coordZsChild[26] = coordZsParent[6];
			coordZsChild[1]  = 0.5 * ( coordZsParent[0] + coordZsParent[1] );
			coordZsChild[3]  = 0.5 * ( coordZsParent[0] + coordZsParent[3] );
			coordZsChild[5]  = 0.5 * ( coordZsParent[1] + coordZsParent[2] );
			coordZsChild[7]  = 0.5 * ( coordZsParent[3] + coordZsParent[2] );
			coordZsChild[19] = 0.5 * ( coordZsParent[4] + coordZsParent[5] );
			coordZsChild[21] = 0.5 * ( coordZsParent[4] + coordZsParent[7] );
			coordZsChild[23] = 0.5 * ( coordZsParent[5] + coordZsParent[6] );
			coordZsChild[25] = 0.5 * ( coordZsParent[7] + coordZsParent[6] );
			coordZsChild[4]  = 0.5 * ( coordZsChild[1]  + coordZsChild[7] );
			coordZsChild[22] = 0.5 * ( coordZsChild[19] + coordZsChild[25] );
			for( int i = 9; i < 18; ++i ){
				coordZsChild[i] = 0.5 * ( coordZsChild[i-9] + coordZsChild[i+9] );
			}
			int iChild(0);
			for( int iz = 0; iz < 2; ++iz ){
				for( int iy = 0; iy < 2; ++iy ){
					for( int ix = 0; ix < 2; ++ix ){
						const int elemIDChild = itr->childs[iChild];
						if( elemIDChild < 0 ){
							continue;
						}
						int nodeIDChild[8];
						for( int iNode = 0; iNode < 8; ++iNode ){
							nodeIDChild[iNode] = m_elemInfo[elemIDChild].nodes[iNode];
						}
						const int offset = ix + iy * 3 + iz * 9;
						m_nodeCoordinates[ nodeIDChild[0] ].Z = coordZsChild[0  + offset];
						m_nodeCoordinates[ nodeIDChild[1] ].Z = coordZsChild[1  + offset];
						m_nodeCoordinates[ nodeIDChild[2] ].Z = coordZsChild[4  + offset];
						m_nodeCoordinates[ nodeIDChild[3] ].Z = coordZsChild[3  + offset];
						m_nodeCoordinates[ nodeIDChild[4] ].Z = coordZsChild[9  + offset];
						m_nodeCoordinates[ nodeIDChild[5] ].Z = coordZsChild[10 + offset];
						m_nodeCoordinates[ nodeIDChild[6] ].Z = coordZsChild[13 + offset];
						m_nodeCoordinates[ nodeIDChild[7] ].Z = coordZsChild[12 + offset];
						++iChild;
					}
				}
			}
		}
		else{
			// Interpolate coordinates of childs
			double coordZsChild[18];
			coordZsChild[0]  = coordZsParent[0];
			coordZsChild[2]  = coordZsParent[1];
			coordZsChild[6]  = coordZsParent[3];
			coordZsChild[8]  = coordZsParent[2];
			coordZsChild[9]  = coordZsParent[4];
			coordZsChild[11] = coordZsParent[5];
			coordZsChild[15] = coordZsParent[7];
			coordZsChild[17] = coordZsParent[6];
			coordZsChild[1]  = 0.5 * ( coordZsParent[0] + coordZsParent[1] );
			coordZsChild[3]  = 0.5 * ( coordZsParent[0] + coordZsParent[3] );
			coordZsChild[5]  = 0.5 * ( coordZsParent[1] + coordZsParent[2] );
			coordZsChild[7]  = 0.5 * ( coordZsParent[3] + coordZsParent[2] );
			coordZsChild[10] = 0.5 * ( coordZsParent[4] + coordZsParent[5] );
			coordZsChild[12] = 0.5 * ( coordZsParent[4] + coordZsParent[7] );
			coordZsChild[14] = 0.5 * ( coordZsParent[5] + coordZsParent[6] );
			coordZsChild[16] = 0.5 * ( coordZsParent[7] + coordZsParent[6] );
			coordZsChild[4]  = 0.5 * ( coordZsChild[1]  + coordZsChild[7] );
			coordZsChild[13] = 0.5 * ( coordZsChild[10] + coordZsChild[16] );
			int iChild(0);
			for( int iy = 0; iy < 2; ++iy ){
				for( int ix = 0; ix < 2; ++ix ){
					const int elemIDChild = itr->childs[iChild];
					if( elemIDChild < 0 ){
						continue;
					}
					int nodeIDChild[8];
					for( int iNode = 0; iNode < 8; ++iNode ){
						nodeIDChild[iNode] = m_elemInfo[elemIDChild].nodes[iNode];
					}
					const int offset = ix + iy * 3;
					m_nodeCoordinates[ nodeIDChild[0] ].Z = coordZsChild[0  + offset];
					m_nodeCoordinates[ nodeIDChild[1] ].Z = coordZsChild[1  + offset];
					m_nodeCoordinates[ nodeIDChild[2] ].Z = coordZsChild[4  + offset];
					m_nodeCoordinates[ nodeIDChild[3] ].Z = coordZsChild[3  + offset];
					m_nodeCoordinates[ nodeIDChild[4] ].Z = coordZsChild[9  + offset];
					m_nodeCoordinates[ nodeIDChild[5] ].Z = coordZsChild[10 + offset];
					m_nodeCoordinates[ nodeIDChild[6] ].Z = coordZsChild[13 + offset];
					m_nodeCoordinates[ nodeIDChild[7] ].Z = coordZsChild[12 + offset];
					++iChild;
				}
			}
		}
	}

}

// Check whether all children belong to land or air
bool MeshGenerationDHexa::doesAllChildrenBelongToLandOrAir( const int elemIndex ) const{

	if( isActive(elemIndex) ){
		if( m_elemInfo[elemIndex].type == LAND || m_elemInfo[elemIndex].type == AIR ){
			return true;
		}else{
			return false;
		}
	}else{
		const ElementInfo& info = m_elemInfo[elemIndex];
		for( int i = 0; i < 8; ++i ){
			const int elemIndexChild = info.childs[i];
			if( elemIndexChild < 0 ){
				continue;
			}
			if( !doesAllChildrenBelongToLandOrAir(elemIndexChild) ){
				return false;
			}
		}
		return true;
	}

}

// Check whether four edges are vertical
void MeshGenerationDHexa::checkWhetherFourEdgesAreVertical( const ElementInfo& info ) const{

	const double eps = 1.0e-6;
	for( int i = 0; i < 4; ++i ){
		const int nodeUpper = info.nodes[i];
		const int nodeLower = info.nodes[i+4];
		const CommonParameters::XYZ coordUpper = m_nodeCoordinates[nodeUpper];
		const CommonParameters::XYZ coordLower = m_nodeCoordinates[nodeLower];
		if( fabs(coordUpper.X - coordLower.X) > eps || fabs(coordUpper.Y - coordLower.Y) > eps ){
			std::cerr << "Error : Edge connectiong nodes " << nodeUpper << " and " << nodeLower << " is not vertical." << std::endl;
			exit(1);
		}
	}

}

// Check whether element has flat top and bottom surfaces
bool MeshGenerationDHexa::checkWhetherElementHasFlatTopAndBottomSurfaces( const ElementInfo& info ) const{

	for( int iTopBot = 0; iTopBot < 2; ++iTopBot ){
		const int nodeRef = info.nodes[iTopBot*4];
		const double zRef = m_nodeCoordinates[nodeRef].Z;
		const double eps = 1.0e-6;
		for( int i = 1; i < 4; ++i ){
			const int node = info.nodes[i+iTopBot*4];
			const double z = m_nodeCoordinates[node].Z;
			if( fabs(z - zRef) > eps ){
				return false;
			}
		}
	}

	return true;

}

// Check element volumes
void MeshGenerationDHexa::checkElementVolumes(){

	int elemIndex(0);
	for (std::vector<ElementInfo>::const_iterator itr = m_elemInfo.begin(); itr != m_elemInfo.end(); ++itr, ++elemIndex) {
		const double volume = calculateVolume(elemIndex);
		if (volume < CommonParameters::EPS) {
			std::cout << "Warning : Too small element volume! " << volume << std::endl;
			std::cout << "          Element index: " << elemIndex << std::endl;
			std::cout << "          Element type : " << itr->type << (itr->isActive ? "" : " (not active)") << std::endl;
			std::cout << "          Level: " << itr->level << std::endl;
			for (int iNode = 0; iNode < 8; ++iNode) {
				std::cout << "          Node#" << iNode << ": (x,y,z)=(";
				const int nodeIndex = itr->nodes[iNode];
				std::cout << m_nodeCoordinates[nodeIndex].X << ",";
				std::cout << m_nodeCoordinates[nodeIndex].Y << ",";
				std::cout << m_nodeCoordinates[nodeIndex].Z << ")" << std::endl;
			}
		}
	}

}

// Calculate volume
double MeshGenerationDHexa::calculateVolume(const int elemIndex) const {

	double volume(0.0);
	for (int ip = 0; ip < 8; ++ip) {
		const double xi = m_integralPointXi[ip];
		const double eta = m_integralPointEta[ip];
		const double zeta = m_integralPointZeta[ip];
		volume += calcDeterminantOfJacobianMatrix(elemIndex, xi, eta, zeta) * m_weights[ip];
	}
	return volume;

}

// Calculate determinant of jacobian matrix of the elements
double MeshGenerationDHexa::calcDeterminantOfJacobianMatrix(const int elemIndex, const double xi, const double eta, const double zeta) const {

	double xCoord[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double yCoord[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double zCoord[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	for (int i = 0; i < 8; ++i) {
		const int nodeID = m_elemInfo[elemIndex].nodes[i];
		xCoord[i] = m_nodeCoordinates[nodeID].X;
		yCoord[i] = m_nodeCoordinates[nodeID].Y;
		zCoord[i] = m_nodeCoordinates[nodeID].Z;
	}

	// Zero clear
	DoubleMatrix3x3 JacobMat;
	JacobMat.comp11 = 0.0;
	JacobMat.comp12 = 0.0;
	JacobMat.comp13 = 0.0;
	JacobMat.comp21 = 0.0;
	JacobMat.comp22 = 0.0;
	JacobMat.comp23 = 0.0;
	JacobMat.comp31 = 0.0;
	JacobMat.comp32 = 0.0;
	JacobMat.comp33 = 0.0;
	for (int i = 0; i < 8; ++i) {
		const double xiNode = m_xiAtNode[i];
		const double etaNode = m_etaAtNode[i];
		const double zetaNode = m_zetaAtNode[i];
		const double tmp1 = 0.125 * xiNode * (1.0 + etaNode * eta) * (1.0 + zetaNode * zeta);
		const double tmp2 = 0.125 * etaNode * (1.0 + zetaNode * zeta) * (1.0 + xiNode * xi);
		const double tmp3 = 0.125 * zetaNode * (1.0 + xiNode * xi) * (1.0 + etaNode * eta);
		JacobMat.comp11 += tmp1 * xCoord[i];
		JacobMat.comp12 += tmp1 * yCoord[i];
		JacobMat.comp13 += tmp1 * zCoord[i];
		JacobMat.comp21 += tmp2 * xCoord[i];
		JacobMat.comp22 += tmp2 * yCoord[i];
		JacobMat.comp23 += tmp2 * zCoord[i];
		JacobMat.comp31 += tmp3 * xCoord[i];
		JacobMat.comp32 += tmp3 * yCoord[i];
		JacobMat.comp33 += tmp3 * zCoord[i];
	}

	const double determinant = JacobMat.comp11 * JacobMat.comp22 * JacobMat.comp33
		+ JacobMat.comp12 * JacobMat.comp23 * JacobMat.comp31
		+ JacobMat.comp13 * JacobMat.comp21 * JacobMat.comp32
		- JacobMat.comp13 * JacobMat.comp22 * JacobMat.comp31
		- JacobMat.comp12 * JacobMat.comp21 * JacobMat.comp33
		- JacobMat.comp11 * JacobMat.comp23 * JacobMat.comp32;

	return determinant;

}

#ifdef _LAYERS
// Read data of layers
void MeshGenerationDHexa::readLayers(std::ifstream& infile){

	if( m_hasLayersRead == true ){
		std::cerr << "Data of layers of the earth has already been read." << std::endl;	
		exit(1);
	}
	else{
		std::cout << "Read data of layers of the earth." << std::endl;	
	}

	infile >> m_numLayers;

	std::cout << m_numLayers << " ";

	m_hasLayerRead = new bool[ m_numLayers ];
	for( int i = 0; i < m_numLayers; ++i ){
		m_hasLayerRead[i] = false;
	}

	m_elemGroupingZ = new int[ m_numLayers + 1 ];
	m_elemGroupingZ[0] = 0;
	for( int i = 1; i < m_numLayers + 1; ++i ){
		infile >> m_elemGroupingZ[i];
		std::cout << m_elemGroupingZ[i] << " " ;
	}
	std::cout << std::endl;

	m_numElemGroupX = new int[m_numLayers];
	m_numElemGroupY = new int[m_numLayers];
	m_elemGroupingX = new int*[m_numLayers];
	m_elemGroupingY = new int*[m_numLayers];

	// The air layer
	m_numElemGroupX[0] = 1;
	m_elemGroupingX[0] = new int[2];
	m_elemGroupingX[0][0] = 0;
	m_elemGroupingX[0][1] = m_numXInit;

	m_numElemGroupY[0] = 1;
	m_elemGroupingY[0] = new int[2];
	m_elemGroupingY[0][0] = 0;
	m_elemGroupingY[0][1] = m_numYInit;
	m_hasLayerRead[0] = true;

	//Data check
	if( m_elemGroupingZ[m_numLayers] != m_numZInit ){
		std::cerr << "The Last number of the data of layers must be equal to total division number of Z direction." << std::endl;
		exit(1);
	}

	m_hasLayersRead = true;

}

// Read data of each layer of the earth
void MeshGenerationDHexa::readLayer(std::ifstream& infile, const int layerIDStart, const int layerIDEnd){

	if( layerIDStart < 2 || layerIDEnd > m_numLayers ){
		std::cerr << "IDs of layers must be and more than 1 and less than or equal to "
			<< m_numLayers << ". : StartID = " << layerIDStart << ", EndID = " << layerIDEnd << std::endl;
		exit(1);
	}

	const int layerIDZeroStart = layerIDStart - 1;
	const int layerIDZeroEnd = layerIDEnd - 1;

	std::cout << "Read data of the layers from " << layerIDStart << " to " << layerIDEnd << "." << std::endl;
	for( int i = layerIDZeroStart; i < layerIDZeroEnd; ++i ){
		if( m_hasLayerRead[i] == true ){
			std::cerr << "Data of layer " << i+1 << " has already been read." << std::endl;
			exit(1);
		}
	}

	// Read data of X direction
	int nGroupX(0);
	infile >> nGroupX;
	int* elemGroupingXWork = NULL; 
	if( nGroupX < 0 ){
		nGroupX = m_numXInit;
		elemGroupingXWork = new int[nGroupX+1];
		elemGroupingXWork[0] = 0;
		for( int i = 1; i < nGroupX + 1; ++i ){
			elemGroupingXWork[i] = i;
		}
	}else{
		elemGroupingXWork = new int[nGroupX+1];
		elemGroupingXWork[0] = 0;
		for( int i = 1; i < nGroupX + 1; ++i ){
			infile >> elemGroupingXWork[i];
		}
	}

	// Read data of Y direction
	int nGroupY(0);
	infile >> nGroupY;
	int* elemGroupingYWork = NULL; 
	if( nGroupY < 0 ){
		nGroupY = m_numYInit;
		elemGroupingYWork = new int[nGroupY+1];
		elemGroupingYWork[0] = 0;
		for( int i = 1; i < nGroupY + 1; ++i ){
			elemGroupingYWork[i] = i;
		}
	}else{
		elemGroupingYWork = new int[nGroupY+1];
		elemGroupingYWork[0] = 0;
		for( int i = 1; i < nGroupY + 1; ++i ){
			infile >> elemGroupingYWork[i];
		}
	}

	for( int iLayer = layerIDZeroStart; iLayer <= layerIDZeroEnd; ++iLayer ){
		std::cout << "Layer# " << iLayer+1 << std::endl;

		// X direction
		std::cout << "X direction" << std::endl;
		std::cout << nGroupX << " ";

		m_numElemGroupX[iLayer] = nGroupX;
		m_elemGroupingX[iLayer] = new int[nGroupX+1];

		m_elemGroupingX[iLayer][0] = 0;
		for( int i = 1; i < m_numElemGroupX[iLayer] + 1; ++i ){
			m_elemGroupingX[iLayer][i] = elemGroupingXWork[i]; 
		}

		for( int i = 1; i < m_numElemGroupX[iLayer] + 1; ++i ){
			std::cout << m_elemGroupingX[iLayer][i] << " ";
		}

		std::cout << std::endl;

		// Y direction
		std::cout << "Y direction" << std::endl;
		std::cout << nGroupY << " ";

		m_numElemGroupY[iLayer] = nGroupY;
		m_elemGroupingY[iLayer] = new int[nGroupY+1];

		m_elemGroupingY[iLayer][0] = 0;
		for( int i = 1; i < m_numElemGroupY[iLayer] + 1; ++i ){
			m_elemGroupingY[iLayer][i] = elemGroupingYWork[i]; 
		}

		for( int i = 1; i < m_numElemGroupY[iLayer] + 1; ++i ){
			std::cout << m_elemGroupingY[iLayer][i] << " ";
		}

		std::cout << std::endl;

		m_hasLayerRead[iLayer] = true;

	}

	delete [] elemGroupingXWork;
	delete [] elemGroupingYWork;

}
#endif
