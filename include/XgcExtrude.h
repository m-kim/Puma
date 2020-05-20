#ifndef XGCEXTRUDE_H
#define XGCEXTRUDE_H
#include <string>

#include <vtkm/cont/DataSet.h>

class XgcExtrude
{
public:

    virtual void initializeReaders(std::string meshName, std::string diagName) = 0;

    virtual vtkm::cont::DataSet readMesh() = 0;

    virtual void openADIOS(std::string filename) = 0;
    virtual void readValues(vtkm::cont::DataSet &ds) = 0;

};
#endif