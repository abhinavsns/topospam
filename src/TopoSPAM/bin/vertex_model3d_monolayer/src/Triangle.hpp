#ifndef TRIANGLE_HPP_
#define TRIANGLE_HPP_
#include "Common-inline.h"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

using Eigen::Matrix3d;
using Eigen::Vector3d;

class Triangle;

class TriVertex
{
public:
  // ---- CREATION & DESTRUCTION ----
  TriVertex(){connected_triangles.reserve(16);};
  TriVertex(Vector3d pos_, int id_) : pos(pos_), id(id_){connected_triangles.reserve(16);};
  ~TriVertex(){};    //destructor

  // --- GEOMETRY ---
  Vector3d pos;                // coordinates of the vertices (apical side)

  int id;

  std::vector<Triangle*> connected_triangles;

  Vector3d vec_field = Vector3d(0,0,0); // vector field to be decomposed
  Vector3d pos_sph = Vector3d(0,0,0); // vertex position given in spherical coordinates (r, theta, phi)
  Vector3d r_hat = Vector3d(0,0,0); // computing once, and keeping the spherical unit vectors to speeds up
  Vector3d theta_hat = Vector3d(0,0,0);
  Vector3d phi_hat = Vector3d(0,0,0);
  double cos_theta = 0; // keeping a copy of cos(pos_sph[1]) due to it's common use

};

class TriEdge
{
public:
  // ---- CREATION & DESTRUCTION ----
  TriEdge(){};
  TriEdge(TriVertex *tail_, TriVertex *head_, TriEdge* conj_, int id_) : tail(tail_), head(head_), conj(conj_), id(id_){};
  ~TriEdge(){};    //destructor

  // --- GEOMETRY ---
  TriVertex *tail;
  TriVertex *head;

  TriEdge *conj;

  int id;

  Triangle *tri = nullptr;

  Vector3d l_vec = Vector3d(0,0,0); // head - tail
  double l=-1; //edge length
  // double a_ratio=0;// = tri->area/( tri->area+conj->tri->area )
  // double alpha_n1n2=0; // acos( tri->unit_normal_vec .dot. conj->tri->unit_normal_vec )
  // Vector3d n_a = Vector3d(0,0,0);
  // Vector3d n_i = Vector3d(0,0,0);
  Matrix3d C_int = Matrix3d::Zero();
};

class Triangle
{
public:
  // ---- CREATION & DESTRUCTION ----
  Triangle();
  Triangle(TriVertex *v1_, TriVertex *v2_, TriVertex *v3_, TriEdge *e1_, TriEdge *e2_, TriEdge *e3_, int id_);
  // Vertex(bool is_3D);      //constructor
  virtual ~Triangle();     //destructor

  // --- GEOMETRY ---
  // vertices ordered in clock-wise order
  TriVertex *v1;
  TriVertex *v2;
  TriVertex *v3;

  TriEdge *e1; // v1->v2
  TriEdge *e2; // v2->v3
  TriEdge *e3; // v3->v1

  int id;
  //-----
  Vector3d com; // triangle center of mass

  double area;

  Vector3d normal_vec, unit_normal_vec;

  Vector3d polarity;

  double angle_1;//angle between triangle bonds at v1
  double angle_2;//angle between triangle bonds at v2
  double angle_3;//angle between triangle bonds at v3

  double excess;

  Matrix3d Q_xyplane; // elongation tensor rotated into x-y plane

  Matrix3d C_int;

  double update_excess();
  void update_corner_angles();

};

#endif//TRIANGLE_HPP_
