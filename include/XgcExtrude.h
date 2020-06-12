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

    virtual adios2::Engine kopen(adios2::IO io, std::string fp, std::string fn, adios2::Mode mode, MPI_Comm comm){
        //return io.Open(fp+fn+".bp", adios2::Mode::Read);

        return kittie::open(fn, fp+fn + ".bp", adios2::Mode::Read, comm);
    }

    virtual adios2::IO kdeclare_io(adios2::ADIOS &ad, std::string name)
    {
        // auto io = ad.DeclareIO("BP");
        // io.SetEngine("BP");
        // return io;

        return kittie::declare_io(name);
    }

    virtual void kinit(MPI_Comm comm)
    {
        kittie::initialize(comm, adios2::DebugON);
        // adios = adios2::ADIOS(comm);
        // mesh = adios2::ADIOS(comm);
        // diag = adios2::ADIOS(comm);

    }

};
#endif