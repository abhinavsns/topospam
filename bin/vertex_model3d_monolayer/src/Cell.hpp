#ifndef CELL_HPP_
#define CELL_HPP_

#include <iostream>
#include <list>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

using Eigen::Matrix3d;
using Eigen::Vector3d;

class Vertex;
class Edge;

class Cell
{
public:
  // ---- CREATION & DESTRUCTION ----
  Cell();
  // Vertex(bool is_3D);      //constructor
  virtual ~Cell();     //destructor

  // --- GEOMETRY ---
  std::list<Vertex*> Vlist;
  std::list<Edge*> Elist;
  int id;

  /// constructor that takes the coordinates of the apical and basal point
  Cell(std::list<Vertex*> Vlist, std::list<Edge*> Elist, int id);

  Vector3d cent_a, cent_b, centroid;
  double area_a, area_b; // note: lateral area is assigned to each cell edge
  double perim_a;//apical perimeter
  double perim_b;//basal perimeter
  double vol; // cell volume for the 3D model

  Vector3d polarity, avg_polarity_j;
  Vector3d cent_a_prev;

  double mag_F_a; // magnitude of the traction force directed by polarity

  Vector3d unit_normal_vec;

  Matrix3d Q_3d;

 // Vector3d update_cell_apical_centroid();
  // double get_apical_area();

};

#endif//CELL_HPP_
