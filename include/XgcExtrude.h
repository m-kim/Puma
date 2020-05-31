#ifndef XGCEXTRUDE_H
#define XGCEXTRUDE_H
#include <string>

#include <vtkm/cont/DataSet.h>
#include <adios2.h>

class XgcExtrude
{
public:
    adios2::IO fileIO, meshIO, diagIO;
    adios2::ADIOS adios, mesh, diag;
    adios2::Engine fileReader, meshReader, diagReader;
    bool running = true;

    int numNodes, numTris, numPhi = 4, numTimeSteps;
    vtkm::cont::DataSet ds;
    
    virtual void initializeReaders(std::string meshName, std::string diagName) = 0;

    virtual void readMesh() = 0;

    virtual void openADIOS(std::string filename);
    virtual vtkm::cont::ArrayHandle<double>
        GetiTurbulence(vtkm::cont::ArrayHandle<double> &temperature);

    virtual void readValues() = 0;

};
#endif