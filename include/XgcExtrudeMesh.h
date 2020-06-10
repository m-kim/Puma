#ifndef XGCEXTRUDEMESH_H
#define XGCEXTRUDEMESH_H

#include "XgcExtrude.h"

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

class XgcExtrudeMesh final : public XgcExtrude
{
public:
    vtkm::cont::ArrayHandle<vtkm::Id> wedgeConn;
    vtkm::cont::ArrayHandle<vtkm::Vec<double,3>> points;

    void initializeReaders(std::string mp, std::string mn, 
                            std::string dp, std::string dn,
                            MPI_Comm comm = MPI_COMM_WORLD);

    void readMesh();

    void readValues();

    void close();

    vtkm::cont::ArrayHandle<double> coords;
    std::string meshName, diagName;
};

#endif
