#ifndef XGCEXTRUDECOMPUTE_H
#define XGCEXTRUDECOMPUTE_H
#include "XgcExtrude.h"

#include <adios2.h>
#include <mpi.h>
#include <vtkm/cont/DataSet.h>

class XgcExtrudeCompute : public XgcExtrude
{
public:
    void initializeReaders(std::string meshName, std::string diagName);

    void readMesh();

    void readValues();

    vtkm::cont::CellSetExtrude cells;
    vtkm::cont::ArrayHandleExtrudeCoords<double> coords;
};
#endif

