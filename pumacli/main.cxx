
#include <memory.h>


#include <sstream>
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

#include "XgcExtrudeMesh.h"
#include "XgcExtrudeCompute.h"
#include "TurbulenceWorklets.h"

std::unique_ptr<XgcExtrude> exrt;


const auto
parse(int argc, char **argv){
  int x = 128;
  int y = 128;
  std::string meshpathname, filepathname, filename, meshname;
  std::string diagpathname, diagname, extrudetype;
  diagpathname = filepathname = meshpathname = std::string("/home/adios/adiosvm/Tutorial/xgc/totalf_itg_tiny/");
  meshname = std::string("xgc.mesh.bp");
  filename = std::string("xgc.3d.bp");
  diagname = std::string("xgc.oneddiag.bp");
  extrudetype = std::string("compute");

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
    else if (!strcmp(argv[i], "-diagname"))
    {
      if (i+1 < argc){
          diagname = std::string(argv[i+1]);
          i++;
      }
    }
    else if (!strcmp(argv[i], "-diagpathname"))
    {
      if (i+1 < argc){
          diagpathname = std::string(argv[i+1]);
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
    else if (!strcmp(argv[i], "-extrudetype"))
    {
      if (i+1 < argc){
          extrudetype = std::string(argv[i+1]);
          i++;
      }
    }
  }

  if (extrudetype == "compute"){
    return std::make_tuple(x,y, meshpathname, meshname, filepathname, filename, diagpathname, diagname, 0);
  }
  else
    return std::make_tuple(x,y, meshpathname, meshname, filepathname, filename, diagpathname, diagname, 1);

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
  camera.Elevation(static_cast<vtkm::Float32>(-90.0));
}
void display(int x, int y)
{
    // if (exrt->fileReader->BeginStep() == adios2::StepStatus::OK){
         exrt->readMesh();
    //     exrt->fileReader->EndStep();
    // }
    try{
        int cnt = 0;
        while(exrt->fileReader->BeginStep() ==adios2::StepStatus::OK){
            exrt->readValues();
            vtkm::rendering::Camera camera;
            vtkm::cont::ColorTable colorTable("inferno");

            std::string fieldNm("pointvar");
            vtkm::rendering::MapperRayTracer mapper;
            vtkm::rendering::Color background(1.0f, 1.0f, 1.0f, 1.0f);
            vtkm::rendering::Color foreground(0.0f, 0.0f, 0.0f, 1.0f);
            vtkm::rendering::CanvasRayTracer canvas(x,y);
            vtkm::rendering::Scene scene;
            scene.AddActor(vtkm::rendering::Actor(exrt->ds.GetCellSet(),
                                                  exrt->ds.GetCoordinateSystem(),
                                                  exrt->ds.GetField(fieldNm),
                                                  colorTable));
            SetCamera(camera, exrt->ds.GetCoordinateSystem().GetBounds(),
                      exrt->ds.GetField(fieldNm));

            vtkm::rendering::View3D view(scene, mapper, canvas, camera, background, foreground );
            view.Initialize();
          view.Paint();
          std::stringstream sstr;
          sstr << "output-";
          sstr << cnt << ".pnm";

              view.SaveAs(sstr.str());

            //renderer->Display(ds, canvas, fieldNm);
            exrt->fileReader->EndStep();
            cnt++;
        }
        exrt->fileReader->Close();
        MPI_Finalize();
    }
    catch(int e){
        exrt->fileReader->EndStep();
        exrt->fileReader->Close();
        MPI_Finalize();
    }

}



int main(int argc, char **argv)
{
  auto tups = parse(argc, argv);
  auto meshopen = std::get<2>(tups) + std::get<3>(tups);
  auto fileopen = std::get<4>(tups) + std::get<5>(tups);
  auto diagopen = std::get<6>(tups) + std::get<7>(tups);
//    renderer = std::make_unique<VTKmXeusRender>();
   MPI_Init(NULL,NULL);
  if (std::get<8>(tups))
    exrt = std::make_unique<XgcExtrudeMesh>();
  else
  exrt = std::make_unique<XgcExtrudeCompute>();

  exrt->openADIOS(fileopen);
  //TODO: need diag
  exrt->initializeReaders(meshopen, diagopen);
  display(std::get<0>(tups), std::get<1>(tups));

}

