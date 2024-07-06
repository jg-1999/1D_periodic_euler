#include "FluxFunctions.h"

template<class T>
inline FluxFunction<T>::FluxFunction()
{
  
}

template<class T>
inline LinearFlux<T>::LinearFlux()
{
  c = 1.0;
}

template<class T>
inline void LinearFlux<T>::computeFlux(DataStruct<T> &U, DataStruct<T> &F)
{
  T *dataU = U.getData();
  T *dataF = F.getData();

  for(int n = 0; n < U.getSize(); n++)
  {
    dataF[n] = c * dataU[n];
  }
}

template<class T>
inline T LinearFlux<T>::computeFlux(const T &Ui)
{
  return c * Ui;
}

template class FluxFunction<float>;
template class FluxFunction<double>;

template class LinearFlux<float>;
template class LinearFlux<double>;
