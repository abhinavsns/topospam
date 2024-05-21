#include "Cell.hpp"


//Vertex constructor
Cell::Cell(std::list<Vertex*> Vlist_, std::list<Edge*> Elist_, int id_) : Vlist(Vlist_), Elist(Elist_), id(id_){
  mag_F_a = 0.;// default: cells are not active
  cent_a = Vector3d(0,0,0);
  cent_b = Vector3d(0,0,0);
  centroid = Vector3d(0,0,0);
  polarity = Vector3d(0,0,0);
}

// void Vertex::setNetwork(Network *n) {
//   _net = n;
// }

//Vertex destructor
Cell::~Cell(){
}
