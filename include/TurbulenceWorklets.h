#ifndef TURBULENCEWORKLETS_H
#define TURBULENCEWORKLETS_H

#include <vtkm/Types.h>
#include <vtkm/Math.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/Field.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/cont/DeviceAdapterAlgorithm.h>
#include <vtkm/cont/ArrayHandleTransform.h>
#include <vtkm/cont/ArrayHandleCast.h>

template <typename T>
  struct BiasFunctor
  {
    VTKM_EXEC_CONT
    BiasFunctor(T bias = T(0))
      : Bias(bias)
    {
    }

    VTKM_EXEC_CONT
    T operator()(T x) const { return x - this->Bias; }

    T Bias;
  };

template <typename T>
struct NanFunctor
{
  VTKM_EXEC_CONT
  NanFunctor(T val = T(0))
      : Val(val)
  {
  }

  VTKM_EXEC_CONT
  T operator()(T x) const { if (x==x) return x; else return Val; }

  T Val;
};

template <typename T>
struct NanFlagFunctor
{
  VTKM_EXEC_CONT
  NanFlagFunctor()
  {
  }

  VTKM_EXEC_CONT
  bool operator()(T x) const { return x==x;}
};


template<class FieldType, class DeviceAdapter>
class DoRecenter
{
public:
    typedef vtkm::cont::ArrayHandleTransform<vtkm::cont::ArrayHandle<FieldType>,NanFunctor<FieldType>>
            AHTAverageType;
    typedef vtkm::cont::ArrayHandle<FieldType> OutType;
  DoRecenter()
  {

  }

  OutType
  Run(vtkm::cont::ArrayHandle<FieldType> &vals)
  {
    typedef typename vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter> DeviceAlgorithms;


    const FieldType initVal = 0.0;
    auto nanzeroed_out =
        vtkm::cont::make_ArrayHandleTransform<vtkm::cont::ArrayHandle<FieldType>>(
          vals, NanFunctor<FieldType>(0.0));

    //get the sum of the array, ignoring NaN
    FieldType sum = DeviceAlgorithms::Reduce(nanzeroed_out, initVal, vtkm::Sum());
    //Count the non-NaN.
    auto nanflagged_out =
        vtkm::cont::make_ArrayHandleTransform<vtkm::cont::ArrayHandle<FieldType>>(
          vals, NanFlagFunctor<FieldType>());

    auto nanflagged_sum = vtkm::cont::make_ArrayHandleCast<FieldType>(nanflagged_out);

    vtkm::Id cnt = DeviceAlgorithms::Reduce(nanflagged_sum, initVal, vtkm::Sum());


    //Set the NaN values to something large to be ignored
    AHTAverageType nan_out =
        vtkm::cont::make_ArrayHandleTransform<vtkm::cont::ArrayHandle<FieldType>>(
          vals, NanFunctor<FieldType>(-1e30));
    //subtract the average out

    auto handle_out =
            vtkm::cont::make_ArrayHandleTransform<AHTAverageType>(
              nan_out, BiasFunctor<FieldType>(sum/FieldType(cnt)));

    //CastAndCall2 does not support ArrayHandleTransform?
    OutType true_out;
    true_out.Allocate(vals.GetNumberOfValues());
    DeviceAlgorithms::Copy(handle_out, true_out);
    return true_out;
  }
};

class iTurbulence : public vtkm::worklet::WorkletMapField
{
public:

  iTurbulence(vtkm::Id _numPhi,
             vtkm::Id _numNodes
             ) :
    numPhi(_numPhi),
    numNodes(_numNodes)
  {

  }

  typedef void ControlSignature(
                                FieldIn,
                                WholeArrayIn,
                                FieldIn,
                                FieldIn,
                                FieldIn,
                                FieldIn,
                                WholeArrayOut);


  typedef void ExecutionSignature(
    _1,
    _2,
    _3,
    _4,
    _5,
    _6,
    _7);


  template<typename FieldType, typename FieldArrayType, typename FieldArrayOutType>
  VTKM_EXEC
  void operator()(vtkm::Id &j,
                  FieldArrayType &dpot,
                  FieldType &pot0,
                  FieldType &potm0,
                  FieldType &te,
                  FieldType &de,
                  FieldArrayOutType &out
                ) const
  {
    //for (int j = 0; j < numNodes; j++)
      for (int i=0; i<numPhi; i++)
    {
        int idx = i*numNodes+j;
        double v1 = (dpot.Get(idx) - (potm0-pot0));
        v1 = v1 / te;

        out.Set(idx, v1);
    }
  }
  vtkm::Id numNodes, numPhi;
};
template<typename FieldType, typename DeviceAdapter>
class DoiTurbulence
{
public:
  DoiTurbulence()
  {

  }
  vtkm::cont::ArrayHandle<FieldType>
  Run(
    vtkm::cont::ArrayHandle<FieldType> &dpot,
    vtkm::cont::ArrayHandle<FieldType> &pot0,
    vtkm::cont::ArrayHandle<FieldType> &potm0,
    vtkm::cont::ArrayHandle<FieldType> &te,
    vtkm::cont::ArrayHandle<FieldType> &de,
    vtkm::Id numPhi,
    vtkm::Id numNodes)
  {
    typedef typename vtkm::worklet::DispatcherMapField<iTurbulence>
      TurbulenceWorkletDispatchType;

    vtkm::cont::ArrayHandleCounting<vtkm::Id> idx;
    idx = vtkm::cont::make_ArrayHandleCounting<vtkm::Id>(0, 1, numNodes);

    iTurbulence turbulenceWorklet(numPhi, numNodes);
    TurbulenceWorkletDispatchType dispatch(turbulenceWorklet);

    vtkm::cont::ArrayHandle<FieldType> arrOut;
    arrOut.Allocate(dpot.GetNumberOfValues());

    dpot.PrepareForInPlace(DeviceAdapter());
    pot0.PrepareForInPlace(DeviceAdapter());
    potm0.PrepareForInPlace(DeviceAdapter());
    te.PrepareForInPlace(DeviceAdapter());
    de.PrepareForInPlace(DeviceAdapter());
    arrOut.PrepareForInPlace(DeviceAdapter());

    dispatch.Invoke(idx,
      dpot,
      pot0,
      potm0,
      te, de,
    arrOut);

    return arrOut;

  }
};
class Turbulence : public vtkm::worklet::WorkletMapField
{
public:
  
  Turbulence(vtkm::Id _numPhi,
             vtkm::Id _numNodes
             ) :
    numPhi(_numPhi),
    numNodes(_numNodes)
  {

  }

  typedef void ControlSignature(
                                FieldIn,
                                WholeArrayIn,
                                FieldIn,
                                FieldIn,
                                FieldIn,
                                FieldIn,
                                WholeArrayIn,
                                FieldIn,
                                WholeArrayOut);


  typedef void ExecutionSignature(
    _1, 
    _2, 
    _3,
    _4,
    _5,
    _6,
    _7,
    _8, 
    _9);



  template<typename FieldType, typename FieldArrayType, typename FieldArrayOutType>
  VTKM_EXEC
  void operator()(vtkm::Id &j,
                  FieldArrayType &dpot,
                  FieldType &pot0,
                  FieldType &potm0,
                  FieldType &te,
                  FieldType &de,
                  FieldArrayType &eden,
                  FieldType &meanEden,
                  FieldArrayOutType &out
                ) const
  {
    //for (int j = 0; j < numNodes; j++)
      for (int i=0; i<numPhi; i++)
    {
        int idx = i*numNodes+j;
        double v1 = (dpot.Get(idx) - (potm0-pot0));
        v1 = v1 / te;

        double v2 = eden.Get(idx) - meanEden;
        v2 = v2 / de;

        out.Set(idx, v1 + v2);
    }
  }
  vtkm::Id numNodes, numPhi;
};
template<typename FieldType, typename DeviceAdapter>
class DoTurbulence
{
public:
  DoTurbulence()
  {

  }
  vtkm::cont::ArrayHandle<FieldType> 
  Run(
    vtkm::cont::ArrayHandle<FieldType> &dpot,
    vtkm::cont::ArrayHandle<FieldType> &pot0,
    vtkm::cont::ArrayHandle<FieldType> &potm0,
    vtkm::cont::ArrayHandle<FieldType> &te,
    vtkm::cont::ArrayHandle<FieldType> &de,
    vtkm::cont::ArrayHandle<FieldType> &eden,
    vtkm::cont::ArrayHandle<FieldType> &meanEden,
    vtkm::Id numPhi,
    vtkm::Id numNodes)
  {
    typedef typename vtkm::worklet::DispatcherMapField<Turbulence>
      TurbulenceWorkletDispatchType;

    vtkm::cont::ArrayHandleCounting<vtkm::Id> idx;
    idx = vtkm::cont::make_ArrayHandleCounting<vtkm::Id>(0, 1, numNodes);
    
    Turbulence turbulenceWorklet(numPhi, numNodes);
    TurbulenceWorkletDispatchType dispatch(turbulenceWorklet);

    vtkm::cont::ArrayHandle<FieldType> arrOut;
    arrOut.Allocate(eden.GetNumberOfValues());
    
    dpot.PrepareForInPlace(DeviceAdapter());
    pot0.PrepareForInPlace(DeviceAdapter());
    potm0.PrepareForInPlace(DeviceAdapter());
    te.PrepareForInPlace(DeviceAdapter());
    de.PrepareForInPlace(DeviceAdapter());
    eden.PrepareForInPlace(DeviceAdapter());
    meanEden.PrepareForInPlace(DeviceAdapter());
    arrOut.PrepareForInPlace(DeviceAdapter());

    dispatch.Invoke(idx, 
      dpot, 
      pot0,
      potm0,
      te, de,
      eden,
      meanEden,
    arrOut);

    return arrOut;

  }
};


class Mean : public vtkm::worklet::WorkletMapField
{
public:
  
  Mean(vtkm::Id _numPhi,
        vtkm::Id _numNodes) :
    numPhi(_numPhi),
    numNodes(_numNodes)
  {

  }

  typedef void ControlSignature(
                                FieldIn,
                                WholeArrayIn,
                                FieldInOut);


  typedef void ExecutionSignature(_1, _2, _3);



  template<typename FieldType, typename FieldArrayType>
  VTKM_EXEC
  void operator()(vtkm::Id &idx,
                  FieldArrayType &eden,
                  FieldType &meanEden
                ) const
  {
//    for (int i = 0; i < numPhi; i++)
//        for (int j = 0; j < numNodes; j++)
//            meanEden[j] += eden->GetTuple1(i*numNodes + j);
//    for (int i = 0; i < numNodes; i++)
//        meanEden[i] /= (double)numPhi;

    meanEden = 0.0;
    for (int j=0; j<numPhi; j++){
      meanEden += eden.Get(j * numNodes + idx);
    }
    meanEden /= (FieldType)numPhi;
  }


  vtkm::Id numPhi, numNodes;
};

template<typename FieldType, typename DeviceAdapter>
class DoMean
{
public:
  DoMean()
  {

  }
    vtkm::cont::ArrayHandle<FieldType> 
    Run(vtkm::cont::ArrayHandle<FieldType> &eden,
    vtkm::Id numPhi,
    vtkm::Id numNodes)
  {
    typedef typename vtkm::worklet::DispatcherMapField<Mean>
      MeanWorkletDispatchType;

    vtkm::cont::ArrayHandleCounting<vtkm::Id> idx;
    idx = vtkm::cont::make_ArrayHandleCounting<vtkm::Id>(0, 1, numNodes);
    
    Mean meanWorklet(numPhi, numNodes);
    MeanWorkletDispatchType dispatch(meanWorklet);

    vtkm::cont::ArrayHandle<FieldType> meanEden;
    meanEden.Allocate(numNodes);
    
    eden.PrepareForInPlace(DeviceAdapter());
    meanEden.PrepareForInPlace(DeviceAdapter());
    dispatch.Invoke(idx, eden, meanEden);

    return meanEden;

  }
};

class Interpolate : public vtkm::worklet::WorkletMapField
{
public:
  
  Interpolate(vtkm::Id _n): n(_n)
  {

  }

  typedef void ControlSignature(
                                FieldIn,
                                FieldIn,
                                WholeArrayIn,
                                WholeArrayIn,
                                FieldInOut);


  typedef void ExecutionSignature(_1, _2, _3, _4, _5);



  template<typename FieldType, typename FieldProcessedType, typename FieldArrayType>
  VTKM_EXEC
  void operator()(
                  FieldType &xi,
                  FieldProcessedType &tag,
                  FieldArrayType &x,
                  FieldArrayType &y,
                  FieldType &yi
                ) const
  {
    if (!tag){
      for (int j = 0; j < n-1; j++)
      {
        if (xi >= x[j] && xi <= x[j+1])
        {
          double dy = y[j+1]-y[j];
          double t = (xi-x[j])/(x[j+1]-x[j]);
          yi = y[j] + t*dy;
          break;

        }
      }
    }
}

  vtkm::Id n;
};
template<typename FT, typename FieldProcessedType>
class AvgTemp
{
public:
  typedef  FT FieldType;
  VTKM_CONT
  AvgTemp(){}

  VTKM_CONT
  AvgTemp(FieldType _n, FieldType _d) :
    numerator(_n),
    denominator(_d)
  {

  }

  VTKM_EXEC
  void Evaluate(FieldType &temp1,
    FieldType &temp2,
    FieldProcessedType &temp
    ) const
  {
    temp = numerator*(temp1 + temp2 ) / denominator;
  }

  FieldType numerator, denominator;
};

template<typename FT, typename FieldProcessedType>
class SetMin
{
public:
  typedef  FT FieldType;
  VTKM_CONT
  SetMin(){}

  VTKM_CONT
  SetMin(FieldType _x0, FieldType _y0) :
    x0(_x0),
    y0(_y0)
  {

  }
  
  VTKM_EXEC
  void Evaluate(FieldType &xi,
    FieldType &yi,
    FieldProcessedType &tag
    ) const
  {
    yi = 0.0;
    tag = xi <= x0;
    if (tag)
      yi = y0;
  }

  FieldType x0, y0;
};
template<typename FT,typename FieldProcessedType>
class SetMax
{
public:
  typedef  FT FieldType;
  VTKM_CONT
  SetMax(){}

  VTKM_CONT
  SetMax(FieldType _xmax, FieldType _ymax) :
    xmax(_xmax),
    ymax(_ymax)
  {

  }
  
  VTKM_EXEC
  void Evaluate(FieldType &xi,
    FieldType &yi,
    FieldProcessedType &tag
    ) const
  {
    if (!tag && xi >= xmax){
      tag = true;
      yi = ymax;
    }
  }

  FieldType xmax, ymax;
};

template< typename EvaluatorType>
class Evaluate : public vtkm::worklet::WorkletMapField
{
public:
  
  Evaluate(EvaluatorType _eval):
          eval(_eval)
  {

  }

  typedef void ControlSignature(
                                FieldInOut,
                                FieldInOut,
                                FieldInOut);


  typedef void ExecutionSignature(_1, _2, _3);



  template<typename FieldType, typename FieldProcessedType>
  VTKM_EXEC
  void operator()(
                  FieldType &xi,
                  FieldType &yi,
                  FieldProcessedType &tag
                  ) const
  {
    // yi = 0.0;
    // tag = xi <= x0;
    // if (tag)
    //   yi = y0;
    eval.Evaluate(xi,yi,tag);
  }

  EvaluatorType eval;
};

template<typename FieldType, typename DeviceAdapter>
class DoInterpolate
{
public:
  DoInterpolate()
  {

  }

  vtkm::cont::ArrayHandle<FieldType>
    Run(
      vtkm::cont::ArrayHandle<FieldType> &x,
      vtkm::cont::ArrayHandle<FieldType> &y,
      vtkm::cont::ArrayHandle<FieldType> &xi) {
    typedef typename vtkm::worklet::DispatcherMapField<Evaluate<SetMin<FieldType, bool>>>
      MinWorkletDispatchType;
    typedef typename vtkm::worklet::DispatcherMapField<Evaluate<SetMax<FieldType, bool>>>
      MaxWorkletDispatchType;
    typedef typename vtkm::worklet::DispatcherMapField<Interpolate>
      InterpWorkletDispatchType;

    vtkm::cont::ArrayHandle<FieldType> yi;
    FieldType x0, y0, xmax, ymax;

    yi.Allocate(xi.GetNumberOfValues());
    x0 = x.GetPortalConstControl().Get(0);
    y0 = y.GetPortalConstControl().Get(0);
    xmax = x.GetPortalConstControl().Get(x.GetNumberOfValues() - 1);
    ymax = y.GetPortalConstControl().Get(x.GetNumberOfValues() - 1);

    SetMin<FieldType, bool> evalMin(x0,y0);
    SetMax<FieldType, bool> evalMax(xmax, ymax);
    Evaluate<SetMin<FieldType, bool>> MinWorklet(evalMin);
    Evaluate<SetMax<FieldType, bool>> MaxWorklet(evalMax);
    MinWorkletDispatchType minDispatch(MinWorklet);
    MaxWorkletDispatchType maxDispatch(MaxWorklet);
    xi.PrepareForInPlace(DeviceAdapter());
    yi.PrepareForInPlace(DeviceAdapter());
    x.PrepareForInPlace(DeviceAdapter());
    y.PrepareForInPlace(DeviceAdapter());

    vtkm::cont::ArrayHandle<bool> tag;
    tag.Allocate(xi.GetNumberOfValues());

    minDispatch.Invoke(xi,yi,tag);
    maxDispatch.Invoke(xi,yi, tag);
  
    Interpolate interpWorklet(x.GetNumberOfValues());
    InterpWorkletDispatchType interpDispatch(interpWorklet);
    interpDispatch.Invoke(xi, tag, x, y, yi);
    return yi;

  }
};
#endif
