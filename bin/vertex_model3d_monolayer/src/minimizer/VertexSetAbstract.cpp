#include <vector>

#include "../MyException.h"
#include "../Log.h"

#include "../Parameters.h"    // needed for CACHE_DISTANCES

#include "../Vertex.h"
#include "../Cell.h"
#include "../Bond.h"
#include "../Frame.h"

#include "VertexSetAbstract.h"

#include "../Frame-inline.h"
#include "../Bond-inline.h"

using namespace std;

VertexSetAbstract::VertexSetAbstract(Frame *frame) : _frame(frame) {
  if(!_frame) {
    throw MyException("VertexSetAbstract::VertexSetAbstract: _frame==NULL!");
  }
  prepareForMinimization();
}

VertexSetAbstract::~VertexSetAbstract() {
  cleanUp();
}

void VertexSetAbstract::prepareForMinimization() {
  resetFlags();
  verticesForMinimization.clear();
  for(unsigned int i = 0; i < vertices.size(); ++i) {
    if(vertices[i]->get_to_be_minimized()){
      verticesForMinimization.push_back(vertices[i]);
    }
  }
}

void VertexSetAbstract::loadCoordinateVector(const AbstractVectorConst &v) {
  // copy from position vector
  for(unsigned int i = 0; i < verticesForMinimization.size(); ++i) {
    verticesForMinimization[i]->x = v.get(2*i);
    verticesForMinimization[i]->y = v.get(2*i+1);
  }

  int nextIndex = 2*verticesForMinimization.size();
  if(hasFrame) {
    if(hasFrameX) {
      _frame->setXWidth(v.get(nextIndex++));
    }
    if(hasFrameY) {
      _frame->setYWidth(v.get(nextIndex++));
    }
    if(hasFrameShear) {
      _frame->setSkewVariable(v.get(nextIndex++));
    }
  }
}

void VertexSetAbstract::fillCoordinateVector(AbstractVector &v) {
  // copy to position vector
  for(unsigned int i = 0; i < verticesForMinimization.size(); ++i) {
    v.set(2*i, verticesForMinimization[i]->x);
    v.set(2*i+1, verticesForMinimization[i]->y);
  }

  int nextIndex = 2*verticesForMinimization.size();
  if(hasFrame) {
    if(hasFrameX) {
      v.set(nextIndex++, _frame->xWidth());
    }
    if(hasFrameY) {
      v.set(nextIndex++, _frame->yWidth());
    }
    if(hasFrameShear) {
      v.set(nextIndex++, _frame->skewVariable());
    }
  }
}

void VertexSetAbstract::fillForceVector(AbstractVector &v) {
  // copy to force vector
  for(unsigned int i = 0; i < verticesForMinimization.size(); ++i) {
    v.set(2*i, verticesForMinimization[i]->negForceX);
    v.set(2*i+1, verticesForMinimization[i]->negForceY);
  }

  int nextIndex = 2*verticesForMinimization.size();
  if(hasFrame) {
    if(hasFrameX) {
      v.set(nextIndex++, _frame->effectiveNegForceX());
    }
    if(hasFrameY) {
      v.set(nextIndex++, _frame->effectiveNegForceY());
    }
    if(hasFrameShear) {
      v.set(nextIndex++, _frame->effectiveNegForceShear());
    }
  }
}

int VertexSetAbstract::recomputeDistances() {
#if CACHE_DISTANCES
  // distances:
//   #pragma omp parallel for
  for(unsigned int i = 0; i < bonds.size(); ++i) {
    //    if(!bonds[i]->isFrozen()) {
    bonds[i]->recomputeLength();
    //    }
  }
#endif
  return 0;
}

int VertexSetAbstract::recomputeEnergies() {
  recomputeDistances();

  // areas, perimeters, energies:
  for(unsigned int i = 0; i < cells.size(); ++i) {
    cells[i]->recomputeAreaElongationPerimeterEnergy();
  }
  return 0;
}

int VertexSetAbstract::recomputeEnergiesAndForces() {
  recomputeEnergies();

  // recompute forces:
  for(unsigned int i = 0; i < vertices.size(); ++i) {
    vertices[i]->resetForce();
  }
  for(unsigned int i = 0; i < cells.size(); ++i) {
    cells[i]->computeNegForcesOnVertices();
  }
#if SCALING_FRAME
  for(unsigned int i = 0; i < vertices.size(); ++i) {
    vertices[i]->negForceX *= _frame->xWidth();
    vertices[i]->negForceY *= _frame->yWidth();
  }
#endif

  if(hasFrame) {
    _frame->resetForce();
    for(unsigned int i = 0; i < cells.size(); ++i) {
      cells[i]->computeNegForcesOnFrame();
    }
  }

  return 0;
}

int VertexSetAbstract::recomputeForces() {
  return recomputeEnergiesAndForces();
}

int VertexSetAbstract::addVertex(Vertex *v) {
  vertices.push_back(v);

  bool found;
  for(unsigned int i = 0; i < v->bonds.size(); ++i) {
    Bond *b = v->bonds[i];
    // check if we already have added the bond:
    found = false;
    for(unsigned int j = 0; j < bonds.size(); ++j) {
      if(bonds[j] == b) {
        found = true;
        break;
      }
    }
    if(!found) {
      // not found -> add bond together with its conjugated
      bonds.push_back(b);
      bonds.push_back(b->conjBond);
    }

    Cell *c = b->cell;
    // check if we already have added the cell:
    found = false;
    for(unsigned int j = 0; j < cells.size(); ++j) {
      if(cells[j] == c) {
        found = true;
        break;
      }
    }
    if(!found) {
      // not found -> add cell
      cells.push_back(c);
    }
  }
  return 0;
}

int VertexSetAbstract::cleanUp() {
  bonds.clear();
  cells.clear();
  vertices.clear();
  verticesForMinimization.clear();
  return 0;
}

int VertexSetAbstract::debugPrintMaxForce() {
  if(Parameters::DebugThresholdMaxForceForOutput < 0) {
      return 0;
  }
  recomputeEnergiesAndForces();

  double totalForceSq = 0.0, curForce, maxForce = -1.0;
  int maxForceIndex = -2;

  for(unsigned int i = 0; i < vertices.size(); ++i) {
      curForce = fabs(vertices[i]->negForceX);
      totalForceSq += curForce*curForce;
      if(curForce > maxForce) {
          maxForce = curForce;
          maxForceIndex = vertices[i]->id();
      }
      curForce = fabs(vertices[i]->negForceY);
      totalForceSq += curForce*curForce;
      if(curForce > maxForce) {
          maxForce = curForce;
          maxForceIndex = vertices[i]->id();
      }
  }

  if(hasFrameX) {
      curForce = fabs(_frame->negForceX);
      totalForceSq += curForce*curForce;
      if(curForce > maxForce) {
          maxForce = curForce;
          maxForceIndex = -1;
      }
  }
  if(hasFrameY) {
      curForce = fabs(_frame->negForceY);
      totalForceSq += curForce*curForce;
      if(curForce > maxForce) {
          maxForce = curForce;
          maxForceIndex = -1;
      }
  }

  if(totalForceSq > Parameters::DebugThresholdMaxForceForOutput) {
      Log(Log::MinimizerDebug, "\n");
      Log(Log::MinimizerDebug, "Max force (on vertex %2d): %g\n", maxForceIndex, maxForce);
      Log(Log::MinimizerDebug, "             Max forceSq: %g\n", maxForce * maxForce);
      Log(Log::MinimizerDebug, "           Total forceSq: %g\n", totalForceSq);
  }
  return 0;
}

int VertexSetAbstract::debugOnError() {
    Log(Log::MinimizerDebug, "--------\n\n");

    recomputeForces();

    double e0, e1, numericalX, numericalY;
    double deltaDeriv = Parameters::DebugNumericalDerivativeDx;
    for(unsigned int i = 0; i < verticesForMinimization.size(); ++i) {
        Vertex *v = verticesForMinimization[i];

        v->x -= 0.5 * deltaDeriv;
        recomputeEnergies();
        e0 = mechanicalEnergy();
        v->x += deltaDeriv;
        recomputeEnergies();
        e1 = mechanicalEnergy();
        v->x -= 0.5 * deltaDeriv;
        numericalX = (e1 - e0) / deltaDeriv;

        v->y -= 0.5 * deltaDeriv;
        recomputeEnergies();
        e0 = mechanicalEnergy();
        v->y += deltaDeriv;
        recomputeEnergies();
        e1 = mechanicalEnergy();
        v->y -= 0.5 * deltaDeriv;
        numericalY = (e1 - e0) / deltaDeriv;

        recomputeEnergies();

        double thresDiff = Parameters::DebugThresholdForceDiffExactNumericalForOutput;
        double thresAbs = Parameters::DebugThresholdForceAbsoluteForOutput;
        if(    (fabs(v->negForceX - numericalX) > thresDiff)
            || (fabs(v->negForceY - numericalY) > thresDiff)
            || (fabs(v->negForceX) > thresAbs)
            || (fabs(v->negForceY) > thresAbs)) {
            Log(Log::MinimizerDebug, "Vertex %d: (%e, %e)\n", v->id(), v->x, v->y);
            Log(Log::MinimizerDebug, "NegForceX on vertex %d: %g (numerical: %g)\n", v->id(), v->negForceX, numericalX);
            Log(Log::MinimizerDebug, "NegForceY on vertex %d: %g (numerical: %g)\n", v->id(), v->negForceY, numericalY);
            for(unsigned int j=0; j<v->bonds.size(); ++j) {
                Log(Log::MinimizerDebug, "length %d-%d: %e, Bond Force: %e, %e\n", v->id(), v->bonds[j]->leftVertex->id(), v->bonds[j]->length(), v->bonds[j]->negForceX(), v->bonds[j]->negForceY());
                Log(Log::MinimizerDebug, "cell %d has area %e.\n", v->bonds[j]->cell->id(), v->bonds[j]->cell->area);
            }
        }
        for(unsigned int j = 0; j < v->bonds.size(); ++j) {
            if(v->bonds[j]->length() < Parameters::T1Cutoff * 100.0) {
                Log(Log::MinimizerDebug, "bond %d(%d): %d-%d, length: %g, deltaX: %g, deltaY: %g, force: (%g,%g), per: (%d,%d)\n", v->bonds[j]->id(), v->bonds[j]->conjBond->id(), v->id(), v->bonds[j]->leftVertex->id(), v->bonds[j]->length(), v->bonds[j]->deltaX(), v->bonds[j]->deltaY(), v->bonds[j]->negForceX(), v->bonds[j]->negForceY(), v->bonds[j]->xPeriodicity(), v->bonds[j]->yPeriodicity());
//                Log(Log::MinimizerDebug, "Vertex %d: (%.20e, %.20e)\n", v->id(), v->x, v->y);
//                Log(Log::MinimizerDebug, "NegForceX on vertex %d: %g (numerical: %g)\n", v->id(), v->negForceX, numericalX);
//                Log(Log::MinimizerDebug, "NegForceY on vertex %d: %g (numerical: %g)\n", v->id(), v->negForceY, numericalY);
                break;
            }
        }
    }

    if(hasFrameX) {
        _frame->setXWidth(_frame->xWidth() - 0.5 * deltaDeriv);
        recomputeEnergies();
        e0 = mechanicalEnergy();
        _frame->setXWidth(_frame->xWidth() + deltaDeriv);
        recomputeEnergies();
        e1 = mechanicalEnergy();
        _frame->setXWidth(_frame->xWidth() - 0.5 * deltaDeriv);
        numericalX = (e1 - e0) / deltaDeriv;
        Log(Log::MinimizerDebug, "NegForceX on frame: %g (numerical: %g)\n", _frame->negForceX, numericalX);
        recomputeEnergies();
    }

    if(hasFrameY) {
        _frame->setYWidth(_frame->yWidth() - 0.5 * deltaDeriv);
        recomputeEnergies();
        e0 = mechanicalEnergy();
        _frame->setYWidth(_frame->yWidth() + deltaDeriv);
        recomputeEnergies();
        e1 = mechanicalEnergy();
        _frame->setYWidth(_frame->yWidth() - 0.5 * deltaDeriv);
        numericalY = (e1 - e0) / deltaDeriv;
        Log(Log::MinimizerDebug, "NegForceY on frame: %g (numerical: %g)\n", _frame->negForceY, numericalY);
        recomputeEnergies();
    }
    
    if(hasFrameShear) {
        _frame->setSkewVariable(_frame->skewVariable() - 0.5 * deltaDeriv);
        recomputeEnergies();
        e0 = mechanicalEnergy();
        _frame->setSkewVariable(_frame->skewVariable() + deltaDeriv);
        recomputeEnergies();
        e1 = mechanicalEnergy();
        _frame->setSkewVariable(_frame->skewVariable() - 0.5 * deltaDeriv);
        double numericalShear = (e1 - e0) / deltaDeriv;
        Log(Log::MinimizerDebug, "negForceShear on frame: %g (numerical: %g)\n", _frame->negForceShear, numericalShear);
        recomputeEnergies();
    }
    
    Log(Log::MinimizerDebug, "Mechanical Energy: %.20e\n", mechanicalEnergy());

    // ---------

    recomputeForces();
    double totalForceSq = 0.0, curForce, maxForce = -1.0;
    int maxForceIndex = -2, maxForceI = -1;

    for(unsigned int i = 0; i < verticesForMinimization.size(); ++i) {
        curForce = fabs(verticesForMinimization[i]->negForceX);
        totalForceSq += curForce*curForce;
        if(curForce > maxForce) {
            maxForce = curForce;
            maxForceIndex = verticesForMinimization[i]->id();
        }
        curForce = fabs(verticesForMinimization[i]->negForceY);
        totalForceSq += curForce*curForce;
        if(curForce > maxForce) {
            maxForce = curForce;
            maxForceIndex = verticesForMinimization[i]->id();
            maxForceI = i;
        }
    }

    Log(Log::MinimizerDebug, "\n");
    if(hasFrameX) {
        curForce = fabs(_frame->negForceX);
        totalForceSq += curForce*curForce;
        if(curForce > maxForce) {
            maxForce = curForce;
            maxForceIndex = -1;
        }
    }
    if(hasFrameY) {
        curForce = fabs(_frame->negForceY);
        totalForceSq += curForce*curForce;
        if(curForce > maxForce) {
            maxForce = curForce;
            maxForceIndex = -1;
        }
    }
    if(hasFrameShear) {
        curForce = fabs(_frame->negForceShear);
        totalForceSq += curForce*curForce;
        if(curForce > maxForce) {
            maxForce = curForce;
            maxForceIndex = -1;
        }
    }

    Log(Log::MinimizerDebug, "Max force (on vertex %2d): %g\n", maxForceIndex, maxForce);
    Log(Log::MinimizerDebug, "             Total force: %g\n", sqrt(totalForceSq));
    Log(Log::MinimizerDebug, "\n");
    if(maxForceI>=0) {
      Vertex *v = verticesForMinimization[maxForceI];
      for(unsigned int j = 0; j < v->bonds.size(); ++j) {
        Bond *b = v->bonds[j];
        Log(Log::MinimizerDebug, "bond %d(%d): %d-%d, length: %g, deltaX: %g, deltaY: %g, force: (%g,%g), per: (%d,%d)\n", b->id(), b->conjBond->id(), v->id(), b->leftVertex->id(), b->length(), b->deltaX(), b->deltaY(), b->negForceX(), b->negForceY(), b->xPeriodicity(), b->yPeriodicity());
        Cell *c = b->cell;
        Log(Log::MinimizerDebug, "cell %d: area: %g\n", c->id(), c->area);
      }
    }
    Log(Log::MinimizerDebug, "\n");

    if(0) {
        double relPosChange = deltaDeriv / sqrt(totalForceSq);

        for(unsigned int i = 0; i < vertices.size(); ++i) {
            Vertex *v = vertices[i];
            v->x -= relPosChange * v->negForceX;
            v->y -= relPosChange * v->negForceY;
        }
        if(_frame) {
            _frame->setXWidth(_frame->xWidth() - _frame->negForceX * relPosChange);
            _frame->setYWidth(_frame->yWidth() - _frame->negForceY * relPosChange);
        }

        // ---------

        recomputeForces();
        for(unsigned int i = 0; i < verticesForMinimization.size(); ++i) {
            Vertex *v = verticesForMinimization[i];

            v->x -= 0.5 * deltaDeriv;
            recomputeEnergies();
            e0 = mechanicalEnergy();
            v->x += deltaDeriv;
            recomputeEnergies();
            e1 = mechanicalEnergy();
            v->x -= 0.5 * deltaDeriv;
            numericalX = (e1 - e0) / deltaDeriv;

            v->y -= 0.5 * deltaDeriv;
            recomputeEnergies();
            e0 = mechanicalEnergy();
            v->y += deltaDeriv;
            recomputeEnergies();
            e1 = mechanicalEnergy();
            v->y -= 0.5 * deltaDeriv;
            numericalY = (e1 - e0) / deltaDeriv;
            recomputeEnergies();

            double thresFactor = 1000;
            if((fabs(v->negForceX - numericalX) > thresFactor * deltaDeriv) || (fabs(v->negForceY - numericalY) > thresFactor * deltaDeriv)) {
                // 			Log(Log::MinimizerDebug, "Vertex %d: (%e, %e)\n", v->id(), v->x, v->y);
                // 			Log(Log::MinimizerDebug, "NegForceX on vertex %d: %g (numerical: %g)\n", v->id(), v->negForceX, numericalX);
                // 			Log(Log::MinimizerDebug, "NegForceY on vertex %d: %g (numerical: %g)\n", v->id(), v->negForceY, numericalY);
                /*			for(unsigned int j=0; j<v->bonds.size(); ++j) {
                                                Log(Log::MinimizerDebug, "length %d-%d: %e, lineTensionForce: %e\n", v->id(), v->bonds[j]->leftVertex->id(), v->bonds[j]->length(), v->bonds[j]->negForce());
                                                Log(Log::MinimizerDebug, "cell %d has area %e.\n", v->bonds[j]->cell->id(), v->bonds[j]->cell->area);
                                        }*/
            }
            for(unsigned int j = 0; j < v->bonds.size(); ++j) {
                if(v->bonds[j]->length() < Parameters::T1Cutoff * 2.0) {
                    Log(Log::MinimizerDebug, "bond %d-%d, length : %g, force: %g,%g\n", v->id(), v->bonds[j]->leftVertex->id(), v->bonds[j]->length(), v->bonds[j]->negForceX(), v->bonds[j]->negForceY());
                    Log(Log::MinimizerDebug, "Vertex %d: (%e, %e)\n", v->id(), v->x, v->y);
                    Log(Log::MinimizerDebug, "NegForceX on vertex %d: %g (numerical: %g)\n", v->id(), v->negForceX, numericalX);
                    Log(Log::MinimizerDebug, "NegForceY on vertex %d: %g (numerical: %g)\n", v->id(), v->negForceY, numericalY);
                }
            }
        }

        if(_frame) {
            _frame->setXWidth(_frame->xWidth() - 0.5 * deltaDeriv);
            recomputeEnergies();
            e0 = mechanicalEnergy();
            _frame->setXWidth(_frame->xWidth() + deltaDeriv);
            recomputeEnergies();
            e1 = mechanicalEnergy();
            _frame->setXWidth(_frame->xWidth() - 0.5 * deltaDeriv);
            numericalX = (e1 - e0) / deltaDeriv;

            _frame->setYWidth(_frame->yWidth() - 0.5 * deltaDeriv);
            recomputeEnergies();
            e0 = mechanicalEnergy();
            _frame->setYWidth(_frame->yWidth() + deltaDeriv);
            recomputeEnergies();
            e1 = mechanicalEnergy();
            _frame->setYWidth(_frame->yWidth() - 0.5 * deltaDeriv);
            numericalY = (e1 - e0) / deltaDeriv;
            recomputeEnergies();

            Log(Log::MinimizerDebug, "NegForceX on frame: %g (numerical: %g)\n", _frame->negForceX, numericalX);
            Log(Log::MinimizerDebug, "NegForceY on frame: %g (numerical: %g)\n", _frame->negForceY, numericalY);
        }
        Log(Log::MinimizerDebug, "Mechanical Energy: %.20e\n", mechanicalEnergy());
    }

    return 0;
}

