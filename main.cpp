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
#include <iostream>
#include "CommonParameters.h"
#include "MeshGenerationDHexa.h"

int main(){

	std::cout << "Start Mesh Generation. Version " << CommonParameters::version << std::endl;

	MeshGenerationDHexa* meshGen = new MeshGenerationDHexa;
	meshGen->readInputFile(); // Read data from input file

	meshGen->calcInitialMeshData(); // Generate mesh data
	meshGen->calcResisivityDistributionForInitialMesh();
//#ifdef _DEBUG_WRITE
//	meshGen->outputVTK("MeshData_Init.vtk"); // Output VTK file
//#endif

	meshGen->partitionMesh();
#ifdef _MERGE
	meshGen->mergeNodes();
#endif
	std::set<int> elemsSeaToLand;
	meshGen->includeTopography(elemsSeaToLand);
	meshGen->removeInactiveElements();

	meshGen->outputMeshData(); // Output mesh data
	meshGen->outputResistivityData();
	meshGen->outputVTK("MeshData.vtk"); // Output VTK file
	delete meshGen;

	std::cout << "End Mesh Generation." << std::endl;

}
