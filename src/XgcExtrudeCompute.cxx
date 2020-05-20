#include "XgcExtrudeCompute.h"
#include <vtkm/cont/ArrayHandleExtrudeCoords.h>
#include <vtkm/cont/CellSetExtrude.h>

void XgcExtrudeCompute::initializeReaders(std::string meshName, std::string diagName)
{
    fileReader->BeginStep(adios2::StepMode::Read, 0.0f);

    meshReader = std::make_unique<adios2::Engine>(meshIO->Open(meshName, adios2::Mode::Read));
    
    adios2::Variable<int> nVar = meshIO->InquireVariable<int>("n_n");
    adios2::Variable<int> triVar = meshIO->InquireVariable<int>("n_t");
    adios2::Variable<int> phiVar = fileIO->InquireVariable<int>("nphi");

    if (nVar){
        meshReader->Get(nVar, &numNodes, adios2::Mode::Sync);

    }
    if (triVar)
        meshReader->Get(triVar, &numTris, adios2::Mode::Sync);

    if (phiVar){
        fileReader->Get(phiVar, &numPhi, adios2::Mode::Sync);
        std::cout << "phi: " << numPhi << std::endl;
    }
}

vtkm::cont::DataSet XgcExtrudeCompute::readMesh()

{
   std::cout << "numNodes: " << numNodes << ", numTris " << numTris << ", numPhi " << numPhi << std::endl;
   adios2::Variable<double> coordVar = meshIO->InquireVariable<double>("/coordinates/values");

  std::vector<double> buff;

  //const int newPhi = phiMultiplier * numPhi;
  int newPhi = numPhi/2 + 1; //+1 for the picket fence problem

  meshReader->Get(coordVar, buff, adios2::Mode::Sync);

  //TODO: go back to how it was before
  //auto coords = vtkm::cont::make_ArrayHandleExtrudeCoords(buff, newPhi, false,vtkm::Pi()/(newPhi-1), vtkm::CopyFlag::On);
  auto coords = vtkm::cont::make_ArrayHandleExtrudeCoords(buff, newPhi, false, vtkm::CopyFlag::On);
  std::vector<int> ibuffc, ibuffn;

  //vtkDataArray *conn = NULL, *nextNode = NULL;
  // meshFile->ReadScalarData("/cell_set[0]/node_connect_list", timestate, &conn);
  // if (!meshFile->ReadScalarData("/nextnode", timestate, &nextNode))
  //     meshFile->ReadScalarData("nextnode", timestate, &nextNode);
  auto nodeConnectorVar = meshIO->InquireVariable<int>("/cell_set[0]/node_connect_list");
  auto nextNodeVar = meshIO->InquireVariable<int>("nextnode");
  if (!nodeConnectorVar || !nextNodeVar)
      return vtkm::cont::DataSet();

  meshReader->Get(nodeConnectorVar, ibuffc,adios2::Mode::Sync);

  meshReader->Get(nextNodeVar, ibuffn, adios2::Mode::Sync);
  if (ibuffn.size() < 1)
      return vtkm::cont::DataSet();

  auto connectivity = vtkm::cont::make_ArrayHandle(ibuffc, vtkm::CopyFlag::On);
  auto nextNode = vtkm::cont::make_ArrayHandle(ibuffn, vtkm::CopyFlag::On);
  auto cells = vtkm::cont::make_CellSetExtrude(connectivity, coords, nextNode, false);

  vtkm::cont::DataSet ds;
  ds.AddCoordinateSystem(vtkm::cont::CoordinateSystem("coords", coords));
  ds.SetCellSet(cells);

  return ds;

}

void XgcExtrudeCompute::openADIOS(std::string filename)
{
  adios = std::make_unique<adios2::ADIOS>(MPI_COMM_WORLD);
  mesh = std::make_unique<adios2::ADIOS>(MPI_COMM_WORLD);

  int numTimeSteps;

  fileIO = std::make_unique<adios2::IO>(adios->DeclareIO("SST"));
  fileIO->SetEngine("SST");
  meshIO = std::make_unique<adios2::IO>(mesh->DeclareIO("BP"));
  meshIO->SetEngine("BP");

   fileReader = std::make_unique<adios2::Engine>(fileIO->Open(filename, adios2::Mode::Read));
  std::cout << "Open " << filename << std::endl;
  std::cout << __FILE__ << " " << __LINE__ << std::endl;

  const auto variables = fileIO->AvailableVariables();
  std::cout << variables.size() << std::endl;

  for (const auto variablePair : variables) {
    std::cout << "Name: " << variablePair.first;

    for (const auto &parameter : variablePair.second) {
      std::cout << "\t" << parameter.first << ": " << parameter.second << "\n";
    }
  }

}


void XgcExtrudeCompute::readValues(vtkm::cont::DataSet &ds)
{

    auto var = fileIO->InquireVariable<double>("dpot");
    std::vector<double> buff;
    fileReader->Get(var, buff,adios2::Mode::Sync);
    auto dpot = vtkm::cont::make_ArrayHandle(buff, vtkm::CopyFlag::On );

    ds.AddField(vtkm::cont::make_FieldPoint("pointvar",  dpot));
    
}
