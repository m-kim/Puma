
#include <memory.h>


#include <sstream>

#include "XgcExtrudeMesh.h"
#include "XgcExtrudeCompute.h"

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

XgcExtrudeCompute exrt;

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
void display()
{
    vtkm::cont::DataSet ds;
    exrt.initializeReaders("/home/adios/adiosvm/Tutorial/xgc/totalf_itg_tiny/xgc.mesh.bp");
    if (exrt.fileReader->BeginStep() == adios2::StepStatus::OK){
        ds = exrt.readMesh();
        exrt.fileReader->EndStep();
    }
    try{
        int cnt = 0;
        while(exrt.fileReader->BeginStep() ==adios2::StepStatus::OK){
            exrt.readValues(ds);
            vtkm::rendering::Camera camera;
            vtkm::cont::ColorTable colorTable("inferno");

            std::string fieldNm("pointvar");
            vtkm::rendering::MapperRayTracer mapper;
            vtkm::rendering::Color background(1.0f, 1.0f, 1.0f, 1.0f);
            vtkm::rendering::Color foreground(0.0f, 0.0f, 0.0f, 1.0f);
            vtkm::rendering::CanvasRayTracer canvas(128,128);
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
            exrt.fileReader->EndStep();
            cnt++;
        }
        exrt.fileReader->Close();
        MPI_Finalize();
    }
    catch(int e){
        exrt.fileReader->EndStep();
        exrt.fileReader->Close();
        MPI_Finalize();
    }

}


int main()
{
//    renderer = std::make_unique<VTKmXeusRender>();
    MPI_Init(NULL,NULL);
    exrt.openADIOS();
    display();
    exrt.fileReader->Close();
    MPI_Finalize();

}

