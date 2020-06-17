#ifndef XGCEXTRUDECOMPUTE_H
#define XGCEXTRUDECOMPUTE_H
#include "XgcExtrude.h"

#include <adios2.h>
#include <mpi.h>
#include <vtkm/cont/DataSet.h>

namespace internal
{
template<typename T>
struct Clamp
{
    T Min, Max;
    VTKM_EXEC_CONT
    Clamp(){}
    VTKM_EXEC_CONT
    Clamp(vtkm::Float32 min, vtkm::Float32 max):Min(min), Max(max)
    {}
    VTKM_EXEC_CONT 
    T operator ()( T x) const
    {
        if (x > Max)
            return Max;
        if (x < Min)
            return Min;
        return x;
    }
};
}

class XgcExtrudeCompute final : public XgcExtrude
{
public:
    void initializeReaders(std::string mp, std::string mn,
                            std::string dp, std::string dn,
                            MPI_Comm comm = MPI_COMM_WORLD);

    void readMesh();

    void readValues();

    void close();

    void clip(vtkm::cont::ArrayHandle<double>);
    vtkm::cont::ArrayHandleTransform<vtkm::cont::ArrayHandle<double>, internal::Clamp<double>> 
     clamp(vtkm::cont::ArrayHandle<double> output);

    vtkm::cont::CellSetExtrude cells;
    vtkm::cont::ArrayHandleExtrudeCoords<double> coords;
};
#endif

