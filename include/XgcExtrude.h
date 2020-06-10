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
    
    std::string filename, meshname, diagname;
    virtual void initializeReaders(std::string meshPath, std::string meshName, 
                                    std::string diagPath, std::string diagName,
                                    MPI_Comm comm = MPI_COMM_WORLD) = 0;

    virtual void readMesh() = 0;

    virtual void openADIOS(std::string fileName, std::string filePath, MPI_Comm comm = MPI_COMM_WORLD);
    virtual vtkm::cont::ArrayHandle<double>
        GetiTurbulence(vtkm::cont::ArrayHandle<double> &temperature);

    virtual void readValues() = 0;

    virtual adios2::StepStatus beginStep();
    virtual void endStep();
    virtual void close();

};
#endif