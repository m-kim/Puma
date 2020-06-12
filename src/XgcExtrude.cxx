#include "XgcExtrude.h"

#include <kittie.h>
#include "TurbulenceWorklets.h"

void XgcExtrude::openADIOS(std::string fp, std::string fn, MPI_Comm comm)
{
  filename = fn;

  kinit(comm);

  fileIO = kdeclare_io(adios, filename);




   fileReader = kopen(fileIO, fp, fn, adios2::Mode::Read, comm);
  std::cout << "Open " << filename << std::endl;
  std::cout << __FILE__ << " " << __LINE__ << std::endl;

  const auto variables = fileIO.AvailableVariables();
  std::cout << variables.size() << std::endl;

  for (const auto variablePair : variables) {
    std::cout << "Name: " << variablePair.first;

    for (const auto &parameter : variablePair.second) {
      std::cout << "\t" << parameter.first << ": " << parameter.second << "\n";
    }
  }

}

vtkm::cont::ArrayHandle<double>
XgcExtrude::GetiTurbulence(vtkm::cont::ArrayHandle<double> &temperature)
{

    //mark: assume that we don't have numphi, numNodes
    numPhi = -1; numNodes = -1; numTris = -1;
    auto var = fileIO.InquireVariable<int>("nphi");
    if (var){
        fileReader.Get(var, &numPhi, adios2::Mode::Sync);
    }
    var = meshIO.InquireVariable<int>("n_n");
    if (var){
        meshReader.Get(var, &numNodes,adios2::Mode::Sync);

    }
    var = meshIO.InquireVariable<int>("n_t");
    if (var){
        meshReader.Get(var, &numTris,adios2::Mode::Sync);
    }

    std::vector<double> buff;

    vtkm::cont::ArrayHandle<double> pot0, potm0, dpot, psi;
  // if (kittie_found)
  // {
  //   auto vard = fileIO->InquireVariable<double>("dpot");
  //   fileReader->Get(vard, buff, adios2::Mode::Sync);
  //   dpot = vtkm::cont::make_ArrayHandle(buff, vtkm::CopyFlag::On);
    
  //  auto varP0 = fileIO->InquireVariable<double>("pot0");
  //  auto varPm0 = fileIO->InquireVariable<double>("potm0");
  //  auto vpsi = diagIO->InquireVariable<double>("psi");
  //  auto vmks = diagIO->InquireVariable<double>("psi_mks");
  //  auto vdens = diagIO->InquireVariable<double>("i_gc_density_1d");
  //  auto vtemp1 = diagIO->InquireVariable<double>("i_perp_temperature_df_1d");
  //  auto vtemp2 = diagIO->InquireVariable<double>("i_parallel_mean_en_df_1d");
  // }

    auto vard = fileIO.InquireVariable<double>("dpot");
    fileReader.Get(vard, buff, adios2::Mode::Sync);
    dpot = vtkm::cont::make_ArrayHandle(buff, vtkm::CopyFlag::On);
    
   auto varP0 = fileIO.InquireVariable<double>("pot0");
   auto varPm0 = fileIO.InquireVariable<double>("potm0");
   auto vpsi = meshIO.InquireVariable<double>("psi");
   auto vmks = diagIO.InquireVariable<double>("psi_mks");
   auto vdens = diagIO.InquireVariable<double>("i_gc_density_1d");
   auto vtemp1 = diagIO.InquireVariable<double>("i_perp_temperature_df_1d");
   auto vtemp2 = diagIO.InquireVariable<double>("i_parallel_mean_en_df_1d");


    if (!varP0 || !varPm0 || !vpsi || !vmks || !vdens || !vtemp1 || !vtemp2){
        return dpot;
    }
    
  // if (kittie_found)
  // {
  //  fileReader->Get(varP0, buff, adios2::Mode::Sync);
  //  pot0 = vtkm::cont::make_ArrayHandle(buff, vtkm::CopyFlag::On);

  //  fileReader->Get(varPm0, buff, adios2::Mode::Sync);
  //  potm0 = vtkm::cont::make_ArrayHandle(buff, vtkm::CopyFlag::On);

  //  diagReader->Get(vpsi, buff, adios2::Mode::Sync);
  //  psi = vtkm::cont::make_ArrayHandle(buff, vtkm::CopyFlag::On);

  //  vtkm::cont::ArrayHandle<double> psid, dens, temp1, temp2;
  //  diagReader->Get(vmks, buff, adios2::Mode::Sync);
  //  psid = vtkm::cont::make_ArrayHandle(buff, vtkm::CopyFlag::On);

  //  diagReader->Get(vdens, buff, adios2::Mode::Sync);
  //  dens = vtkm::cont::make_ArrayHandle(buff, vtkm::CopyFlag::On);

  //  diagReader->Get(vtemp1, buff, adios2::Mode::Sync);
  //  temp1 = vtkm::cont::make_ArrayHandle(buff, vtkm::CopyFlag::On);

  //  diagReader->Get(vtemp2, buff, adios2::Mode::Sync);
  //  temp2 = vtkm::cont::make_ArrayHandle(buff, vtkm::CopyFlag::On);
  // }

   fileReader.Get(varP0, buff, adios2::Mode::Sync);
   pot0 = vtkm::cont::make_ArrayHandle(buff), vtkm::CopyFlag::On;

   fileReader.Get(varPm0, buff, adios2::Mode::Sync);
   potm0 = vtkm::cont::make_ArrayHandle(buff, vtkm::CopyFlag::On);

   meshReader.Get(vpsi, buff, adios2::Mode::Sync);
   psi = vtkm::cont::make_ArrayHandle(buff, vtkm::CopyFlag::On);

   vtkm::cont::ArrayHandle<double> psid, dens, temp1, temp2;
   diagReader.Get(vmks, buff, adios2::Mode::Sync);
   psid = vtkm::cont::make_ArrayHandle(buff, vtkm::CopyFlag::On);

   diagReader.Get(vdens, buff, adios2::Mode::Sync);
   dens = vtkm::cont::make_ArrayHandle(buff, vtkm::CopyFlag::On);

   diagReader.Get(vtemp1, buff, adios2::Mode::Sync);
   temp1 = vtkm::cont::make_ArrayHandle(buff, vtkm::CopyFlag::On);

   diagReader.Get(vtemp2, buff, adios2::Mode::Sync);
   temp2 = vtkm::cont::make_ArrayHandle(buff, vtkm::CopyFlag::On);



   vtkm::cont::ArrayHandle<double> temp;
   temp.Allocate(temp1.GetNumberOfValues());
   // temp.PrepareForInPlace(DeviceAdapter());
   // temp1.PrepareForInPlace(DeviceAdapter());
   // temp2.PrepareForInPlace(DeviceAdapter());
   typedef typename vtkm::worklet::DispatcherMapField<Evaluate<AvgTemp<double, double>>>
     AvgTempWorkletDispatchType;
   AvgTemp<double, double> evalAvg(2.0,3.0);
   Evaluate<AvgTemp<double,double>> AvgWorklet(evalAvg);
   AvgTempWorkletDispatchType avgDispatch(AvgWorklet);
   avgDispatch.Invoke(temp1,temp2, temp);



   DoInterpolate<double> doInterp;
   auto te = doInterp.Run(psid,temp,psi);
   auto de = doInterp.Run(psid, dens, psi);

   DoMean<double> doMean;


   DoiTurbulence<double> doTurb;
   auto arr = doTurb.Run(dpot,
       pot0,
       potm0,
       te,
       de,
       numPhi,
       numNodes);

   // vtkm::cont::DataSet ds;
   // vtkm::cont::DataSetBuilderExplicit builder;
   // ds = builder.Create(arr);
   // return ds;

   int nt = arr.GetNumberOfValues();
   int nPlane = te.GetNumberOfValues();
   temperature.Allocate(nt);
   for (int i = 0; i < nt; i++)
       temperature.GetPortalControl().Set(i, te.GetPortalControl().Get(i%nPlane));
   return arr;
}

adios2::StepStatus XgcExtrude::beginStep()
{
  auto status = kittie::Couplers[filename]->begin_step();
  fileReader = kittie::Couplers[filename]->engine;

  return status;
}

void XgcExtrude::endStep()
{
  kittie::Couplers[filename]->end_step();
}
void XgcExtrude::close()
{
  kittie::Couplers[filename]->close();
}
