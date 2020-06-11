
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
#include <vtkm/cont/ColorTable.h>

#include <kittie.h>

#include "XgcExtrudeMesh.h"
#include "XgcExtrudeCompute.h"
#include "TurbulenceWorklets.h"

std::unique_ptr<XgcExtrude> exrt;

vtkm::cont::ColorTable makeColorTableHotDesaturated()
{
  /* Hot desaturated */
static const float ct_hot_desaturated[] = {
    0.0f,   0.28f, 0.28f, 0.86f,
    0.143f, 0.f,   0.f,   0.36f,
    0.285f, 0.f,   1.f,   1.f,
    0.429f, 0.f,   0.5f,  0.f,
    0.571f, 1.f,   1.f,   0.f,
    0.714f, 1.f,   0.38f, 0.f,
    0.857f, 0.42f, 0.f,   0.f,
    1.0f,   0.88f, 0.3f,  0.3f,
    };
  std::vector<double> rgb(8*3);
  std::vector<double> alpha(8);

  for (int i=0; i<alpha.size(); i++){
    rgb[i*3] = ct_hot_desaturated[i*4];
    rgb[i*3+1] = ct_hot_desaturated[i*4+1];
    rgb[i*3+2] = ct_hot_desaturated[i*4+2];
    alpha[i] = ct_hot_desaturated[i*4+3];
  }
  
  

  return vtkm::cont::ColorTable("ct_hot_desaturated",
                            vtkm::cont::ColorSpace::RGB,
                            vtkm::Vec<double,3>(0,0,0),
                            rgb, alpha);


}
const auto
parse(int argc, char **argv){
  int x = 128;
  int y = 128;
  std::string meshpathname, filepathname, filename, meshname;
  std::string diagpathname, diagname, extrudetype;
  diagpathname = filepathname = meshpathname = std::string("/home/adios/adiosvm/Tutorial/xgc/totalf_itg_tiny/");
  meshname = std::string("xgc.mesh");
  filename = std::string("xgc");
  diagname = std::string("xgc.oneddiag");
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
        while(true){
          auto status = exrt->beginStep();

          if (status == adios2::StepStatus::NotReady)
            continue;
          else if (status != adios2::StepStatus::OK)
            break;

          exrt->readValues();
          vtkm::rendering::Camera camera;
          vtkm::cont::ColorTable colorTable = makeColorTableHotDesaturated();

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
            exrt->endStep();
            cnt++;
        }
        exrt->close();
        kittie::finalize();
    }
    catch(int e){
        exrt->endStep();
        exrt->close();
        kittie::finalize();
    }
}



int main(int argc, char **argv)
{
  auto tups = parse(argc, argv);

//    renderer = std::make_unique<VTKmXeusRender>();
   MPI_Init(NULL,NULL);
  if (std::get<8>(tups))
    exrt = std::make_unique<XgcExtrudeMesh>();
  else
  exrt = std::make_unique<XgcExtrudeCompute>();

  exrt->openADIOS(std::get<4>(tups), std::get<5>(tups));
  //TODO: need diag

  exrt->initializeReaders(std::get<2>(tups), std::get<3>(tups), 
                          std::get<6>(tups), std::get<7>(tups));

  display(std::get<0>(tups), std::get<1>(tups));

}

