#ifndef EDGE_HPP_
#define EDGE_HPP_

#include "Triangle.hpp"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

using Eigen::Matrix3d;
using Eigen::Vector3d;

class Vertex;
class Cell;

class Edge
{
public:
  // ---- CREATION & DESTRUCTION ----
  Edge();
  // Vertex(bool is_3D);      //constructor
  virtual ~Edge();     //destructor

  // --- GEOMETRY ---
  Vertex *tail;                // coordinates of the vertices (apical side)
  Vertex *head;                // coordinates of the vertices (apical side)

  Edge *conj;
  Edge *next;
  Cell *cell;

  int id;

  Vector3d l_a_vec, l_b_vec;
  double l_a, l_b;

  double area_lateral; // lateral area; used for the 3D model
  Vector3d cent_l; // lateral centroid of the bond

  Triangle *tri_a, *tri_b;

  /// constructor that takes the coordinates of the apical and basal point
  Edge(Vertex *tail, Vertex *head, Edge *conj, int id);

};

#endif//EDGE_HPP_
