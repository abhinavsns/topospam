#include <math.h>
#include <stdio.h>

#include "../MyException.h"
#include "../Log.h"
#include "../Parameters.h"

#include "EnergyMinimizerAbstract.h"

using namespace std;

void EnergyMinimizerAbstract::reset() {
  vertexSet->prepareForMinimization();
  resetMinimizer();
  maxIterations = Parameters::MaxIterationsPerDof*vertexSet->numEquations();
  initialized = true;
}

EnergyMinimizerAbstract::Status EnergyMinimizerAbstract::step() {
  // not initialized
  if(!initialized) {
    reset();
  }
  
  // NOTE: I can assume that distances etc are already up to date, here, since:
  //  1. upon initialization, distances are updated
  //  2. after each change of vertex/frame coordinates via loadCoordinateVector, distances are updated
  //  3. after each change of vertex/frame coordinates via frame routines, distances are updated TODO

  // check if vertex set enforces finish
  if(vertexSet->forceFinished()) {
    return Finished;
  }

  // check for topological updates
  if(vertexSet->updateTopology()){
    // topology changed -> have to reinitialize
    initialized = false;
    if(vertexSet->numEquations()==0) {
      return NothingLeft;
    }
    return TopologicalChanges;
  }

  return minimizerStep();
}


EnergyMinimizerAbstract::Status EnergyMinimizerAbstract::relax(bool stopOnTopologicalChanges) {
  //PERFORM MINIMIZATION
  long iter = 0;
  int retry = 0;
  int status;
  if(Parameters::DebugPrintForcesBeforeFirstStep) {
    vertexSet->debugOnError();
  }
  do {
    status = step();
    if(Parameters::DebugPrintForcesAfterEachStep) {
      vertexSet->debugOnError();
    }
    switch(status) {
      case Success:
        ++iter;
        break;
      case TopologicalChanges:
        if(stopOnTopologicalChanges) {
          Log(Log::Minimizer, "Stopping minimization due to topological changes after %d iterations...\n", iter);
          return TopologicalChanges;
        }
        iter = 0;
        break;
      case ExactLine:
        Log(Log::MinimizerImportant, "Starting exact line minimization after %d iterations...\n", iter);
        iter = 0;
        break;
      case StartAgain:
        if(retry<Parameters::NumberOfMinimizerRetries) {
          ++retry;
          Log(Log::MinimizerImportant, "Minimization retry #%d after %d iterations...\n", retry, iter);
        } else {
          if(Parameters::DebugOnMinimizerError) {
            vertexSet->debugOnError();
          }
          if(Parameters::ExceptionIfMinimizationError) {
            throw MyException("EnergyMinimizerAbstract::relax: No progress during iteration %d after %d retries.", iter, retry);
          } else {
            Log(Log::MinimizerImportant, "EnergyMinimizerAbstract::relax: No progress during iteration %d after %d retries.\n", iter, retry);
          }
        }
        iter = 0;
        break;
      case Error:
        if(Parameters::DebugOnMinimizerError) {
          vertexSet->debugOnError();
        }
        if(Parameters::ExceptionIfMinimizationError) {
          throw MyException("EnergyMinimizerAbstract::relax: Error during iteration %d.", iter);
        } else {
          Log(Log::MinimizerImportant, "EnergyMinimizerAbstract::relax: Error during iteration %d.\n", iter);
        }
        return Error;
      case Unknown:
        return Unknown;
      case NothingLeft:
        return NothingLeft;
    }

    if(iter >= maxIterations) {
//       vertexSet->debugOnError();
      Log(Log::MinimizerImportant, "Maximum number of iterations %ld reached.\n", maxIterations);
      return MaxIterations;
    }

  } while(status != Finished);

  vertexSet->debugPrintMaxForce();
  Log(Log::Minimizer, "Minimum for %d degrees of freedom found after %d iterations, mechanical energy: %e\n", (int)vertexSet->numEquations(), iter, vertexSet->mechanicalEnergy());
  return Finished;
}
