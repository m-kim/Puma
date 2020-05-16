#ifndef XGCEXTRUDEMESH_H
#define XGCEXTRUDEMESH_H

#include <vtkm/cont/DeviceAdapter.h>
#include <vtkm/cont/testing/MakeTestDataSet.h>
#include <vtkm/cont/testing/Testing.h>
#include <vtkm/rendering/Actor.h>
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/MapperWireframer.h>
#include <vtkm/rendering/MapperRayTracer.h>
#include <vtkm/rendering/Scene.h>
#include <vtkm/rendering/View3D.h>
#include <vtkm/rendering/testing/RenderTest.h>
#include <adios2.h>
#include <mpi.h>

#include <memory.h>

class XgcExtrudeMesh
{
public:

    const int phiMultiplier = 8;
    int numNodes, numTris, numPhi, numTimeSteps;

    bool running = true;
    std::unique_ptr<adios2::IO> fileIO;
    std::unique_ptr<adios2::IO> meshIO;

    std::unique_ptr<adios2::Engine> fileReader, meshReader;

    vtkm::cont::ArrayHandle<vtkm::Id> wedgeConn;
    vtkm::cont::ArrayHandle<vtkm::Vec<double,3>> points;

    void initializeReaders(std::string meshName);

    vtkm::cont::DataSet readMesh(
                std::unique_ptr<adios2::IO> &meshIO,
                std::unique_ptr<adios2::Engine> &meshReader);

    vtkm::cont::DataSet readValues();
};

#endif