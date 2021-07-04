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
#ifndef DBLDEF_MESH_GENERATION_DHEXA
#define DBLDEF_MESH_GENERATION_DHEXA

#include "Ellipsoids.h"
#include "Cuboids.h"
#include "CommonParameters.h"
#include "ObservingSiteList.h"
#include "TopographyData.h"
#include <fstream>
#include <vector>
#include <map>
#include <set>

enum BoundaryPlanes{
	YZMinus = 0,
	YZPlus,
	ZXMinus,
	ZXPlus,
	XYMinus,
	XYPlus,
};

enum ElementType{
	LAND = 0,
	SEA,
	AIR,
};

enum PartitioningType{
	NO_PARTITIONING = 0,
	HORIZONTAL_ONLY,
	FULL_PARTITIONING
};

struct DataOfAnomaly{
	double XStart;
	double XEnd;
	double YStart;
	double YEnd;
	double ZStart;
	double ZEnd;
	double resistivityValue;
	bool fixResistivityValue;
};

struct ElementAndFace{
	int element;
	int face;
};

struct ElementInfo{
	int type;
	int parameterCell;
	int level;
	bool isActive;
	int nodes[8];
	int levelNeib[6];
	int neib[6][4];
	int parent;
	int childs[8];
};

class MeshGenerationDHexa{
	public:
		// Constructer
		explicit MeshGenerationDHexa();

		// Destructer
		~MeshGenerationDHexa();
		
		// Read input file
		void readInputFile();

		// Caluculate initial mesh data
		void calcInitialMeshData();

		// Caluculate resistivity distribution for initial mesh
		void calcResisivityDistributionForInitialMesh();

		// Reconstruct resistivity distribution to assign different parameter cell to each subsurface element
		void reconstructResisivityDistribution( const std::set<int>& elemsSeaToLand, const std::map<int, double>& elemsLandToSea );

		// Reconstruct the array of elements just below the Earth surface
		void reconstructElementOfEarthSurface();

		// Partition mesh
		void partitionMesh();

#ifdef _LAYERS
		// Caluculate accumulated numbers of resisivity blocks 
		void calcNumResisivityBlockAccumulated();

		// Caluculate ID of resisivity block
		int calcResisivityBlockID( const int ix, const int iy, const int iz );
#endif
		// Output mesh data
		void outputMeshData() const;

		// Output VTK file
		void outputVTK( const std::string& fileName ) const;

		// Output resistivity data
		void outputResistivityData() const;

		// Remove inactive elements and merge nodes
		void removeInactiveElementsAndMergeNodes();

		// Include topography
		void includeTopography( std::set<int>& elemsSeaToLand );

		//// Include topography for land area changed from the sea
		//void includeTopographyForLandAreaChangedFromSea( std::set<int>& surfElemsLand, std::set<int>& surfElemSea );

		// Merge nodes
		void mergeNodes();

		// Remove inactive elements
		void removeInactiveElements();

		// Type of partitioning
		int getPartitioningType() const;

		// Check whether each parameter cell contains at least one active element
		void checkWhetherEachParameterCellContainsAtLeastOneActiveElement() const;

	private:
		// Ellipsoids used for specifing edge lentgh
		Ellipsoids m_ellipsoids;

		// Cuboids used for specifing edge lentgh
		Cuboids m_cuboids;

		// Observation site list for specifing edge lentgh
		ObservingSiteList m_observingSiteList;

		static double distanceConversion;

		// Initial division number (X direction)
		int m_numXInit;

		// Initial division number (Y direction)
		int m_numYInit;

		// Initial division number (Z direction)
		int m_numZInit;

		// X coordinates of initial edges
		double* m_CoordinatesXInit;

		// Y coordinates of initial edges
		double* m_CoordinatesYInit;

		// Z coordinates of initial edges
		double* m_CoordinatesZInit;

#ifdef _LAYERS
		// Total number of layers
		int m_numLayers;

		// Element grouping of each layer (Z direction, Accumulated values))
		int* m_elemGroupingZ;

		// Number of element group of each layer (X direction)
		int* m_numElemGroupX;

		// Number of element group of each layer (Y direction)
		int* m_numElemGroupY;

		// Element grouping of each layer (X direction, Accumulated values))
		int** m_elemGroupingX;

		// Element grouping of each layer (Y direction, Accumulated values))
		int** m_elemGroupingY;

		// Caluculate accumulated numbers of resisivity blocks 
		int* m_numResisivityBlockAccumulated;
#endif

		// Information of elements
		std::vector<ElementInfo> m_elemInfo;

		// Array converting parameter cell index to resistivity values
		std::map<int, double> m_parameterCellToResistivity;

		// Array converting parameter cell index to flag specifing whether the resistivity is fixed or not
		std::map<int, int> m_parameterCellToFixFlag;

		// Location of nodes
		std::vector<CommonParameters::XYZ> m_nodeCoordinates;
	
		// Array of faces of elements corresponding to the boundary planes
		//   m_facesOfElementsBoundaryPlanes[0] : Y-Z Plane ( Minus Side )
		//   m_facesOfElementsBoundaryPlanes[1] : Y-Z Plane ( Plus Side  )
		//   m_facesOfElementsBoundaryPlanes[2] : Z-X Plane ( Minus Side )
		//   m_facesOfElementsBoundaryPlanes[3] : Z-X Plane ( Plus Side  )
		//   m_facesOfElementsBoundaryPlanes[4] : X-Y Plane ( Minus Side )
		//   m_facesOfElementsBoundaryPlanes[5] : X-Y Plane ( Plus Side  )
		//std::vector<ElementAndFace> m_elementAndFaceBoundaryPlanes[6];
		std::set<int> m_elementOnBoundaryPlanes[6];

		// Array of elements just below the Earth surface
		std::set<int> m_elementOfEarthSurface;

		// Array of elements just below the Earth surface (including inactive elements)
		std::set<int> m_elementOfEarthSurfaceWithInactiveElements;

		// Array of elements just below the sea surface
		std::set<int> m_elementOfSeaSurface;

		// Array of elements just below the sea surface (including inactive elements)
		std::set<int> m_elementOfSeaSurfaceWithInactiveElements;

		// Initial resistivity value
		double m_initialResistivity;

		// Resistivity value of the air
		double m_airResistivity;

		// Resistivity value of the sea
		double m_seaResistivity;

		// Sea depth
		double m_seaDepth;

		// Number of resistivity anomalies
		int m_numResisivityAnomalies;

		// Anomaly data
		struct DataOfAnomaly *m_anomalyData;

		// Index of the edges corresponding to the sea surface
		int m_edgeIndexOfSeaSurface;

		// Index of the edges corresponding to the Earth's surface
		int m_edgeIndexOfEarthSurface;

		//// Flag specifing wheter initial parameter cell distribution is saved
		//bool m_saveinitialParameterCellDistribution;

		// Upper limit of the level parameter cell partitioning
		int m_levelLimitParameterCellPartitioning;

		// Whether or not "DIVISION NUMBERS" has already read
		bool m_hasDivisionNumberRead;

		// Whether or not "X COORDINATES" has already read
		bool m_hasXCoordinatesRead;

		// Whether or not "Y COORDINATES" has already read
		bool m_hasYCoordinatesRead;

		// Whether or not "Z COORDINATES" has already read
		bool m_hasZCoordinatesRead;

		// Whether or not "INITIAL RESISTIVITY" has already read
		bool m_hasInitialResistivityRead;

		// Whether or not "AIR RESISTIVITY" has already read
		bool m_hasAirResistivityRead;
	
		// Whether or not "SEA RESISTIVITY" has already read
		bool m_hasSeaResistivityRead;

		// Whether or not "SEA DEPTH" has already read
		bool m_hasSeaDepthRead;

		// Whether or not "ANOMALIES" has already read
		bool m_hasAnomaliesRead;

#ifdef _LAYERS
		// Whether or not "LAYERS" has already read
		bool m_hasLayersRead;

		// Whether or not "LAYER" of each layer has already read
		bool* m_hasLayerRead;
#endif

		// Whether or not the sea is included
		bool m_includeSea;
		
		// Get flag specifing whether inputted element is active 
		bool m_incorporateTopo;

		// Flag specifing whether the threshold of the distance between the level boundary and observation stations is used
		bool m_isUsedThresholdsLimitDistanceBetweenLevelBoundaryAndObsSites;

		// Type of partitioning
		int m_partitioningType;

		// Minimum  sea depth
		double m_minimumSeaDepth;

		// Instance of the class TopographyData
		TopographyData* m_topographyData;

		// Get whether all children are converted to the sea
		bool allChildrenAreConvertedToSea( const int elemIndex, const std::map<int, double>& elemsLandToSea ) const;

		// Get flag specifing whether inputted element is active 
		bool isActive( const int elemID ) const;

		// Get whether the level of neighbor element is higher
		bool doesNeigoborHaveHigherLevel( const int elemID, const int iNeib ) const;

		// Search edge indexes of the sea surface and the Earth's surface
		void searchEdgeIndexesOfSeaAndEarthSurface();

		// Calculate the maximum horizontal edge length at each node;
		void calcMaxEdgeLengthHorizontalAtNode( std::multimap<int, std::pair<int, double> >& nodeToElemsAndEdgeLenth ) const;

		// Calculate the maximum vertical edge length at each node;
		void calcMaxEdgeLengthVerticalAtNode( std::multimap<int, std::pair<int, double> >& nodeToElemsAndEdgeLenth ) const;

		//// Calculate length between inputted two coordinates
		//double calcLength( const CommonParameters::XYZ& coord1, const CommonParameters::XYZ& coord2 ) const;

		// Partition mesh one time
		int partitionMeshOneTime();

		// Partition an element
		void partitionAnElementFull( const int elemIDParent );

		// Partition an element only in horizontal direction
		void partitionAnElementHorizontalDirectionOnly( const int elemIDParent );

		// Change elements on boundary plane
		void changeElementsOnBoundaryPlanesFull( const int elemParent, const int newElemIDStart );

		// Change elements on boundary plane only in horizontal direction
		void changeElementsOnBoundaryPlanesHorizontalDirectionOnly( const int elemParent, const int newElemIDStart );

		// Calculate z-coordinate after including topography/bathymetry
		double calcZCoordAfterIncludingTopography( const double zMin, const double zMax, const double zCoordSurfCur, const double zCoordSurfTarget,
			const CommonParameters::XYZ& coordCur, const bool isSea ) const;

		// Insert to vertical nodes array
		void insertToVerticalNodesArray( const int elemID, std::vector< std::set<int> >& verticalElementsArray,
			std::vector< CommonParameters::XY >& horizontalCoordsArray ) const;

		// Check whether a node is connected to the element with lower level
		bool isConnectedToLowerLevelElement( const std::set<int>& nodeIDs, const int level, const std::multimap<int, int>& nodeToElem ) const;

		// Check whether a node is connected to the element with lower level except top and bottom
		bool isConnectedToLowerLevelElementExceptTopAndBottom( const std::set<int>& nodeIDs, const int level, const double zMin, const double zMax, 
			const std::multimap<int, int>& nodeToElem ) const;

		// Check whether a node is not connected to eight elements with the same level except top and bottom
		bool isNotConnectedToEightSameLevelElementsExceptTopAndBottom( const std::set<int>& nodeIDs, const int level, const double zMin, const double zMax, 
			const std::multimap<int, int>& nodeToElem ) const;

		// Check whether a node is connected to the element with lower level
		bool isConnectedToLowerLevelElement( const int nodeID, const int level, const std::multimap<int, int>& nodeToElem ) const;

		// Get number of element connected to the node
		int getNumOflElementsConnected( const int nodeID, const int level, const std::multimap<int, int>& nodeToElem ) const;

		// Get resistivity from element
		double getResistivityFromElement( const ElementInfo& info ) const;

		// Auxiliary function for including topography
		void includeTopographyAux( const int iLevel, const int maxLevel, const std::vector<CommonParameters::XYZ>& nodeCoordOrg, const std::multimap<int,int>& nodeToElems );

		// Auxiliary function for including topography for land are changed from the sea
		void includeTopographyAux2( const int iLevel, const int maxLevel, const std::vector<CommonParameters::XYZ>& nodeCoordOrg, const std::multimap<int,int>& nodeToElems );

		//// Auxiliary function for including topography for land area changed from the sea
		//void includeTopographyForLandAreaChangedFromSeaAux( const int iLevel, const std::vector<CommonParameters::XYZ>& nodeCoordOrg,
		//	std::set<int>& surfElemsLand, std::set<int>& surfElemSea, const std::multimap<int,int>& nodeToElems );

		// Select elements to be changed to land from the sea
		void selectElementsToBeChangedToLandFromSea( std::set<int>& elemsSeaToLand );

		// Select elements to be changed to seat from land
		void selectElementsToBeChangedToSeaFromLand( std::map<int, double>& elemsLandToSea, std::map<int, double>& elemEarthSurfToSeaDepth ) const;

		// Auxiliary function for selecting elements to be changed to seat from land
		void selectElementsToBeChangedToSeaFromLandAux( const int elemIndex, const double resistivity, std::map<int, double>& elemsLandToSea ) const;

		// Calculated average Z coordinate in an element
		double calcAverageZCoord( const int elemIndex ) const;

		// Change sea to land
		bool changeSeaToLand( const int elemIndex, const std::set<int>& elemsSeaToLand, const int paramCellIDNew );

		// Output Earth surface depth by vtk
		void outputEarthSurfaceDepthByVtk( const std::map<int, double>& elemEarthSurfToSeaDepth ) const;

		// Calculate center coordinate
		CommonParameters::XY calcCenterCoordinate( const int elemIndex ) const;

		// Calculate center coordinate of element face
		CommonParameters::XYZ calculateCenterCoordOfElemFace( const std::vector<ElementInfo>::const_iterator& itrElem, const int iFace ) const;

		// Give same parameter cell ID to children
		void giveSameParamCellIDToChildren( const int elemIndex, const int paramCellID, const std::map<int, double>& elemsLandToSea  );

		// Give same type to children
		void giveSameTypeToChildren( const int elemIndex, const int type );

		// Check whether children have sea area
		bool doChildrenHaveSea( const int elemIndex ) const;

		// Linear interpolation of coordinates for children
		void linearInterpolationOfCoordinatesForChildren( const int level, const bool selectOnlyLand, const std::set<int>& elementsLand );

		// Check whether all children belong to land or air
		bool doesAllChildrenBelongToLandOrAir( const int elemIndex ) const;

		// Check whether four edges are vertical
		void checkWhetherFourEdgesAreVertical( const ElementInfo& info ) const;

		// Check whether element has flat top and bottom surfaces
		bool checkWhetherElementHasFlatTopAndBottomSurfaces( const ElementInfo& info ) const;

#ifdef _LAYERS
		// Read data of layers
		void readLayers(std::ifstream& infile);

		// Read data of each layer
		void readLayer(std::ifstream& infile, const int layerIDStart, const int layerIDEnd);
#endif


};
#endif

