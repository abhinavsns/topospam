#ifndef ENERGY_MINIMIZER_GSL_H
#define ENERGY_MINIMIZER_GSL_H

#include <vector>

// faster gsl_vector_set, gsl_vector_get
#include <gsl/gsl_multimin.h>

// #include "VertexSetAbstract.h"
// #include "EnergyMinimizerAbstract.h"

using namespace std;

// class Frame;
// class Vertex;

class EnergyMinimizerGsl {
public:
  EnergyMinimizerGsl();
  virtual ~EnergyMinimizerGsl();

protected:
  int resetMinimizer();
  Status minimizerStep();

private:
  // GSL minimizer data:
  gsl_multimin_function_fdf gslFunctionDesc;
  gsl_multimin_fdfminimizer *gslState;
  gsl_vector *gslVariables;
  double GslInitialStepSize, currentLineTolerance, cutoffForAbort, GslForceTolerancePerDof, GslForceToleranceTotal;
  int num_degrees_freedom;

  int cleanUp();
  void setCutoffForAbort();
};

#endif
