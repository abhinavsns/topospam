#include "Vertex.hpp"

//Vertex constructor
Vertex::Vertex( Vector3d pos_a_, Vector3d pos_b_, int id_) : pos_a(pos_a_),pos_b(pos_b_),id(id_){

}

// void Vertex::setNetwork(Network *n) {
//   _net = n;
// }

//Vertex destructor
Vertex::~Vertex(){
}
