#include "Edge.hpp"

//Vertex constructor
Edge::Edge(Vertex *tail_, Vertex *head_, Edge *conj_, int id_) : tail(tail_),head(head_), conj(conj_), id(id_){
  cent_l = Vector3d(0,0,0);

}

// void Vertex::setNetwork(Network *n) {
//   _net = n;
// }

//Vertex destructor
Edge::~Edge(){
}
