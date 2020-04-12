#include <adios2.h>
#include <mpi.h>

#include <memory.h>

#include <vtkm/cont/CoordinateSystem.h>
#include <vtkm/cont/CoordinateSystem.hxx>
#include <vtkm/cont/ArrayHandleExtrudeCoords.h>
#include <vtkm/cont/ArrayHandleExtrudeField.h>
#include <vtkm/cont/CellSetExtrude.h>
#include <vtkm/cont/testing/MakeTestDataSet.h>
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/View3D.h>
#include <vtkm/rendering/Scene.h>
#include <vtkm/rendering/MapperRayTracer.h>
#include <vtkm/rendering/Camera.h>
#include <sstream>

std::unique_ptr<adios2::IO> fileIO;
std::unique_ptr<adios2::IO> meshIO;
std::unique_ptr<adios2::ADIOS> adios;
std::unique_ptr<adios2::ADIOS> mesh;
std::unique_ptr<adios2::Engine> fileReader, meshReader;
bool running = true;

int numNodes, numTris, numPhi = 4, numTimeSteps;

const auto
parse(int argc, char **argv){
  int x = 128;
  int y = 128;
  std::string meshpathname, filepathname, filename, meshname;
  filepathname = meshpathname = std::string("/home/adios/adiosvm/Tutorial/xgc/totalf_itg_tiny/");
  meshname = std::string("xgc.mesh.bp");
  filename = std::string("xgc.3d.bp");

  for (int i=1; i<argc; i++){
    if (!strcmp(argv[i], "-x")){
      if (i+1 < argc)
      {
        x = atoi(argv[i+1]);
        i += 1;
      }

    }
    else if (!strcmp(argv[i], "-y"))
    {
      if (i+1 < argc){
        y = atoi(argv[i+1]);
        i += 1;
      }
    }
    else if (!strcmp(argv[i], "-meshname"))
    {
      if (i+1 < argc){
          meshname = std::string(argv[i+1]);
          i++;
      }
    }
    else if (!strcmp(argv[i], "-meshpathname"))
    {
      if (i+1 < argc){
          meshpathname = std::string(argv[i+1]);
          i++;
      }
    }
    else if (!strcmp(argv[i], "-filename"))
    {
      if (i+1 < argc){
          filename = std::string(argv[i+1]);
          i++;
      }
    }
    else if (!strcmp(argv[i], "-filepathname"))
    {
      if (i+1 < argc){
          filepathname = std::string(argv[i+1]);
          i++;
      }
    }
  }
  return std::make_tuple(x,y, meshpathname, meshname, filepathname, filename);
}

void initializeReaders(std::string meshName)
{
    fileReader->BeginStep(adios2::StepMode::NextAvailable, 0.0f);

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
vtkm::cont::DataSet readMesh()
{
   std::cout << "numNodes: " << numNodes << ", numTris " << numTris << ", numPhi " << numPhi << std::endl;
   adios2::Variable<double> coordVar = meshIO->InquireVariable<double>("/coordinates/values");

  std::vector<double> buff;

  //const int newPhi = phiMultiplier * numPhi;
  int newPhi = numPhi/2 + 1; //+1 for the picket fence problem

  meshReader->Get(coordVar, buff, adios2::Mode::Sync);

  auto coords = vtkm::cont::make_ArrayHandleExtrudeCoords(buff, newPhi, false,vtkm::Pi()/(newPhi-1), vtkm::CopyFlag::On);
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
void openADIOS(std::string filename)
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
void
readValues(vtkm::cont::DataSet &ds)
{

    auto var = fileIO->InquireVariable<double>("dpot");
    std::vector<double> buff;
    fileReader->Get(var, buff,adios2::Mode::Sync);
    auto dpot = vtkm::cont::make_ArrayHandle(buff, vtkm::CopyFlag::On );

    ds.AddField(vtkm::cont::make_FieldPoint("pointvar",  dpot));
    
}

inline void SetCamera(vtkm::rendering::Camera& camera,
                       const vtkm::Bounds& coordBounds,
                       const vtkm::cont::Field&)
{

    vtkm::Bounds b = coordBounds;
  b.Z.Min = 0;
  b.Z.Max = 4;
  camera = vtkm::rendering::Camera();
  camera.ResetToBounds(b);
  camera.Azimuth(static_cast<vtkm::Float32>(45.0));
  camera.Elevation(static_cast<vtkm::Float32>(45.0));
}
void display(int x, int y)
{
    vtkm::cont::DataSet ds;
    if (fileReader->BeginStep() == adios2::StepStatus::OK){
        ds = readMesh();
        fileReader->EndStep();
    }
    try{
        int cnt = 0;
        while(fileReader->BeginStep() ==adios2::StepStatus::OK && running){
            readValues(ds);
            vtkm::rendering::Camera camera;
            vtkm::cont::ColorTable colorTable("inferno");

            std::string fieldNm("pointvar");
            vtkm::rendering::MapperRayTracer mapper;
            vtkm::rendering::Color background(1.0f, 1.0f, 1.0f, 1.0f);
            vtkm::rendering::Color foreground(0.0f, 0.0f, 0.0f, 1.0f);
            vtkm::rendering::CanvasRayTracer canvas(x,y);
            vtkm::rendering::Scene scene;
            scene.AddActor(vtkm::rendering::Actor(ds.GetCellSet(),
                                                  ds.GetCoordinateSystem(),
                                                  ds.GetField(fieldNm),
                                                  colorTable));
            SetCamera(camera, ds.GetCoordinateSystem().GetBounds(),
                      ds.GetField(fieldNm));

            vtkm::rendering::View3D view(scene, mapper, canvas, camera, background, foreground );
            view.Initialize();
          view.Paint();
          std::stringstream sstr;
          sstr << "output-";
          sstr << cnt << ".pnm";

              view.SaveAs(sstr.str());

            //renderer->Display(ds, canvas, fieldNm);
            fileReader->EndStep();
            cnt++;
        }
        fileReader->Close();
        MPI_Finalize();
    }
    catch(int e){
        fileReader->EndStep();
        fileReader->Close();
        MPI_Finalize();
    }

}


int main(int argc, char **argv)
{
  auto tups = parse(argc, argv);
  auto meshopen = std::get<2>(tups) + std::get<3>(tups);
  auto fileopen = std::get<4>(tups) + std::get<5>(tups);
//    renderer = std::make_unique<VTKmXeusRender>();
    MPI_Init(NULL,NULL);
    openADIOS(fileopen);
    initializeReaders(meshopen);
    display(std::get<0>(tups), std::get<1>(tups));
    fileReader->Close();
    MPI_Finalize();

}

