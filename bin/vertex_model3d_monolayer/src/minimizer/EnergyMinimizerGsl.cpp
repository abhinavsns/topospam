#include <math.h>
#include <vector>

#include "EnergyMinimizerGsl.h"

using namespace std;
//---------------------- MINIMIZATION ----------------------

EnergyMinimizerGsl::EnergyMinimizerGsl() {
  gslState = NULL;
  gslVariables = NULL;
  currentLineTolerance = 0.05;
  GslInitialStepSize = 0.01;
  GslForceToleranceTotal = 0.0001;
  GslForceTolerancePerDof = 1.0e-5;
  setCutoffForAbort();
}

EnergyMinimizerGsl::~EnergyMinimizerGsl() {
  cleanUp();
}

void EnergyMinimizerGsl::setCutoffForAbort() {
  if(GslForceToleranceTotal>0) {
    cutoffForAbort = GslForceToleranceTotal;
  } else {
    cutoffForAbort = sqrt(num_degrees_freedom)*GslForceTolerancePerDof;
  }
}

int EnergyMinimizerGsl::resetMinimizer() {
  cleanUp();

  // ---- reset gsl minimizer ----

  // allocate the minimizer
  const gsl_multimin_fdfminimizer_type *T;
  //INITIALIZE THE MINIMIZER: (1) SET METHOD, (2) ALLOCATE MEMORY (3) COMPUTE GRADIENT AND SET INITIAL DIRECTION
  //T = gsl_multimin_fdfminimizer_conjugate_fr;
  T = gsl_multimin_fdfminimizer_conjugate_pr;
  //T = gsl_multimin_fdfminimizer_vector_bfgs2;
  //T = gsl_multimin_fdfminimizer_steepest_descent;
  gslState = gsl_multimin_fdfminimizer_alloc(T, num_degrees_freedom);

  //ALLOCATE AND INITIALIZE THE STARTING VALUES
  gslVariables = gsl_vector_alloc(num_degrees_freedom);
  // store the current network coordinates in gsl vector
  GslVector gslVarWrapper(gslVariables);
  vertexSet->fillCoordinateVector(gslVarWrapper);

  // initialize the minimizer
  gslFunctionDesc.n = num_degrees_freedom;
  gslFunctionDesc.f = &gslEnergy;
  gslFunctionDesc.df = &gslNegForce;
  gslFunctionDesc.fdf = &gslEnergyAndNegForce;
  gslFunctionDesc.params = vertexSet;

  gsl_multimin_fdfminimizer_set(gslState, &gslFunctionDesc, gslVariables, GslInitialStepSize, currentLineTolerance);

  // abort condition:
  setCutoffForAbort();

  return 0;
}

int EnergyMinimizerGsl::cleanUp() {
  if(gslState) {
    gsl_multimin_fdfminimizer_free(gslState);
    gslState = NULL;
  }
  if(gslVariables) {
    gsl_vector_free(gslVariables);
    gslVariables = NULL;
  }
  initialized = false;
  return 0;
}


inline EnergyMinimizerGsl::Status EnergyMinimizerGsl::minimizerStep() {
  //test result
  int status = gsl_multimin_test_gradient(gslState->gradient, cutoffForAbort);
  if(status == GSL_SUCCESS){
    // copy best estimate to vertex set:
    vertexSet->loadCoordinateVector(GslVectorConst(gsl_multimin_fdfminimizer_x(gslState)));
    // done
    return Finished;
  } else if (status != GSL_CONTINUE) {
    throw MyException("EnergyMinimizer::step: gsl_multimin_test_gradient delivered unexpected status: %d!", status);
    return Unknown;
  }

  //ITERATE MINIMIZER
  status = gsl_multimin_fdfminimizer_iterate(gslState);

  if(status){
    //AN ERROR OCCURRED DURING MINIMIZATION
    if(currentLineTolerance>0.0){
      //TRY MINIMIZING WITH EXACT LINE MINIMIZATION
      currentLineTolerance = 0.0;
      // have to reinitialize:
      initialized = false;
//       return Error;  // for debugging
      return ExactLine;
    } else {
      if((status==GSL_ENOPROG) && Parameters::GslStartAgainIfNoProgress) {
        // have to reinitialize:
        initialized = false;
        return StartAgain;
      } else {
//           vertexSet->debugOnError();
        Log(Log::MinimizerImportant, "GSL error while minimizing, status: %d!\n", status);
        return Error;
      }
    }
  }

  return Success;
}
