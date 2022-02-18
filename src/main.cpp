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
#include "CommonParameters.h"
#include "MeshGenerationDHexa.h"

int main(){

	std::cout << "Start Mesh Generation. Version " << CommonParameters::version << std::endl;

	MeshGenerationDHexa* meshGen = new MeshGenerationDHexa;
	meshGen->readInputFile(); // Read data from input file

	meshGen->calcInitialMeshData(); // Generate mesh data
	meshGen->calcResisivityDistributionForInitialMesh();
#ifdef _DEBUG_WRITE
	meshGen->outputVTK("MeshData_Init.vtk"); // Output VTK file
#endif

	meshGen->partitionMesh();
#ifdef _MERGE
	meshGen->mergeNodes();
#endif
	std::set<int> elemsSeaToLand;
	meshGen->includeTopography(elemsSeaToLand);
	meshGen->removeInactiveElements();

	meshGen->checkWhetherEachParameterCellContainsAtLeastOneActiveElement();

	meshGen->outputMeshData(); // Output mesh data
	meshGen->outputResistivityData();
	meshGen->outputVTK("MeshData.vtk"); // Output VTK file
	delete meshGen;

	std::cout << "End Mesh Generation." << std::endl;

}
