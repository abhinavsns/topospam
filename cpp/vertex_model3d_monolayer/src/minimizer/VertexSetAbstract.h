#ifndef VERTEX_SET_ABSTRACT_H
#define VERTEX_SET_ABSTRACT_H

#include "../Tissue.h"
// #include "../Cell.h"
// #include "../Frame-inline.h"

// class Vertex;
// class Bond;

class AbstractVectorConst {
public:
  virtual double get(int i) const = 0;
};

class AbstractVector : public AbstractVectorConst {
public:
  virtual void set(int i, double value) = 0;
};

class VertexSetAbstract {
public:
  VertexSetAbstract(Frame *frame);
  virtual ~VertexSetAbstract();

  int numEquations() { return 2*verticesForMinimization.size() + (hasFrame? ((hasFrameX?1:0) + (hasFrameY?1:0) + (hasFrameShear?1:0)) : 0); }
  void loadCoordinateVector(const AbstractVectorConst &v);
  void fillCoordinateVector(AbstractVector &v);
  void fillForceVector(AbstractVector &v);

  int recomputeDistances();
  int recomputeEnergies();
  int recomputeForces();
  int recomputeEnergiesAndForces();
  double mechanicalEnergy() const {
    double totalEnergy = 0.0;
    for(unsigned int i=0; i<cells.size(); ++i) {
      totalEnergy += cells[i]->mechanicalEnergy;
    }
    if(_frame) {
      return totalEnergy + _frame->energy();
    } else {
      return totalEnergy;
    }
  }

  virtual void prepareForMinimization();
  virtual void resetFlags() {
    hasFrame = true;
    hasFrameX = !_frame->xIsochoric();
    hasFrameY = !_frame->yIsochoric();
    hasFrameShear = !_frame->simpleShearIsochoric();
  }
  virtual bool forceFinished() { return false; }
  virtual bool updateTopology() { return false; }

  // DEBUG
  int debugOnError();
  int debugPrintMaxForce();

protected:
  Frame *_frame;
  SearchVector<Vertex *> vertices;
  SearchVector<Vertex *> verticesForMinimization;
  SearchVector<Cell *> cells;
  SearchVector<Bond *> bonds;

  bool hasFrame, hasFrameX, hasFrameY, hasFrameShear;

  virtual int addVertex(Vertex *v);
  int cleanUp();
};

#endif
