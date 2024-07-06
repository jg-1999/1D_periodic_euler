#include "rk4.h"
#include <cassert>

template<class T>
RungeKutta4<T>::RungeKutta4(DataStruct<T> &_Un) : Un(_Un)
{
  nSteps = 4;
  currentStep = 0;

  coeffs = new T[2];
  coeffs[0] = 0.5;
  coeffs[1] = 0.5;

  currentFi.setSize(Un.getSize());
  accumulatedFi.setSize(Un.getSize());

  Ui.setSize(Un.getSize());
}

template<class T>
RungeKutta4<T>::~RungeKutta4()
{
  delete[] coeffs;
}

template<class T>
int RungeKutta4<T>::getNumSteps()
{
  return nSteps;
}

template<class T>
void RungeKutta4<T>::initRK()
{
  currentStep = 0;
  accumulatedFi.zero();  
}

template<class T>
void RungeKutta4<T>::stepUi(T dt)
{
  assert(currentStep < nSteps);

  T *dataUi = Ui.getData();
  const T *dataU = Un.getData();

  if (currentStep == 0)
  {
    for (int n = 0; n < Ui.getSize(); n++)
    {
      dataUi[n] = dataU[n];
    }
  }
  else
  {
    T *dataFi = currentFi.getData();

    for (int n = 0; n < Ui.getSize(); n++)
    {
      dataUi[n] = dataU[n] + coeffs[0] * dt * dataFi[n];
    }
  }
}

template<class T>
void RungeKutta4<T>::finalizeRK(const T dt)
{
  T *dataUn = Un.getData();
  T *dataUi = Ui.getData();
  const T b1 = 2.0;
  const T b2 = 1.0;

  for (int n = 0; n < Ui.getSize(); n++)
  {
    dataUi[n] = b1 * currentFi.getData()[n] + b2 * accumulatedFi.getData()[n];
  }

  const T oneDiv6 = 1. / 6.;
  for (int n = 0; n < Ui.getSize(); n++)
  {
    dataUn[n] += dt * oneDiv6 * dataUi[n];
  }
}

template<class T>
void RungeKutta4<T>::setFi(DataStruct<T> &_F)
{
  T *dataFi = currentFi.getData();
  const T *dataF = _F.getData();

  for (int n = 0; n < Ui.getSize(); n++)
  {
    dataFi[n] = dataF[n];
    accumulatedFi.getData()[n] += coeffs[1] * dataF[n];
  }

  currentStep++;
}

template<class T>
DataStruct<T> * RungeKutta4<T>::currentU()
{
  return &Ui;
}

template class RungeKutta4<float>;
template class RungeKutta4<double>;
