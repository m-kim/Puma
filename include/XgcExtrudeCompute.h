#ifndef XGCEXTRUDECOMPUTE_H
#define XGCEXTRUDECOMPUTE_H
#include "XgcExtrude.h"

#include <adios2.h>
#include <mpi.h>
#include <vtkm/cont/DataSet.h>

class XgcExtrudeCompute : public XgcExtrude
{
public:
    std::unique_ptr<adios2::IO> fileIO;
    std::unique_ptr<adios2::IO> meshIO;

    std::unique_ptr<adios2::ADIOS> adios;
    std::unique_ptr<adios2::ADIOS> mesh;
    std::unique_ptr<adios2::Engine> fileReader, meshReader;
    bool running = true;

    int numNodes, numTris, numPhi = 4, numTimeSteps;

    void initializeReaders(std::string meshName, std::string diagName);

    void readMesh();

    void readValues();

    void openADIOS(std::string filename);


    vtkm::cont::DataSet ds;
    vtkm::cont::CellSetExtrude cells;
    vtkm::cont::ArrayHandleExtrudeCoords<double> coords;
};
#endif

