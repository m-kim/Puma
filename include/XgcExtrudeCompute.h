#ifndef XGCEXTRUDECOMPUTE_H
#define XGCEXTRUDECOMPUTE_H
#include "XgcExtrude.h"

#include <adios2.h>
#include <mpi.h>
#include <vtkm/cont/DataSet.h>

class XgcExtrudeCompute final : public XgcExtrude
{
public:
    void initializeReaders(std::string mp, std::string mn,
                            std::string dp, std::string dn,
                            MPI_Comm comm = MPI_COMM_WORLD);

    void readMesh();

    void readValues();

    void close();
    vtkm::cont::CellSetExtrude cells;
    vtkm::cont::ArrayHandleExtrudeCoords<double> coords;
};
#endif

