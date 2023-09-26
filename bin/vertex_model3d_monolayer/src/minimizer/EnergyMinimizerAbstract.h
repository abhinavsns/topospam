#ifndef ENERGY_MINIMIZER_ABSTRACT_H
#define ENERGY_MINIMIZER_ABSTRACT_H

#include "VertexSetAbstract.h"

using namespace std;

class EnergyMinimizerAbstract {
public:
  enum Status {
    Success, Finished, TopologicalChanges, ExactLine, StartAgain, Error, Unknown, MaxIterations, NothingLeft
  };

  EnergyMinimizerAbstract(VertexSetAbstract *vset) : vertexSet(vset), maxIterations(0), initialized(false) { };
  virtual ~EnergyMinimizerAbstract() {};

  void reset();
  Status step();
  Status relax(bool stopOnTopologicalChanges=false);

protected:
  VertexSetAbstract *vertexSet;
  long maxIterations;
  bool initialized;

  virtual int resetMinimizer() = 0;
  virtual Status minimizerStep() = 0;
};

#endif
