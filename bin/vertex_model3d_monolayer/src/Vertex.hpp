#ifndef VERTEX_HPP_
#define VERTEX_HPP_

#include <vector>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

using Eigen::Matrix3d;
using Eigen::Vector3d;

class Vertex
{
public:
  // ---- CREATION & DESTRUCTION ----
  Vertex();
  // Vertex(bool is_3D);      //constructor
  virtual ~Vertex();     //destructor

  // --- GEOMETRY ---
  Vector3d pos_a;                // coordinates of the vertices (apical side)
  Vector3d pos_b;                // coordinates of the vertices (apical side)
  double lateral_distance;

  Vector3d pos_a_prev;     // coordinates of the vertices at the previous time step (apical side)... used for dynamic simulations

  Vector3d dW_a;         // gradient of energy with respect to apical poistion dW/d pos_a
  Vector3d dW_b;         // gradient of energy with respect to basal poistion dW/d pos_b
  Vector3d F_active_a;

  int id;

  /// constructor that takes the coordinates of the apical and basal point
  Vertex(Vector3d pos_a, Vector3d pos_b, int id);

};

#endif//VERTEX_HPP_
