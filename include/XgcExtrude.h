#ifndef XGCEXTRUDE_H
#define XGCEXTRUDE_H
#include <string>

#include <vtkm/cont/DataSet.h>

class XgcExtrude
{
public:

    virtual void initializeReaders(std::string meshName, std::string diagName) = 0;

    virtual void readMesh() = 0;

    virtual void openADIOS(std::string filename) = 0;
    virtual void readValues() = 0;

};
#endif