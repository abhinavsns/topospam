#ifndef TRIANGULATION_HPP_
#define TRIANGULATION_HPP_

#include "Common-inline.h"
#include "PostProcess-inline.h"
#include "Bond-inline.h"
#include "Triangle.hpp"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

using Eigen::Matrix3d;
using Eigen::Vector3d;

struct Patch {
  std::vector<Triangle*> triangles;
  double area; // sum area of contained triangles
  Vector3d unit_normal_vec;// area weighted average of normal vectors of its triangles
  Matrix3d Q_3D; // coarse grained Q tensor in patch plane
  Matrix3d C_int; // coarse grained curvature tensor in patch plane
  Vector3d com; // center of mass of the patch
  unsigned id;
  Vector3d Q_eigvec; //will hold a unit vector
  double Q_norm;
  double rotated_theta; // speherical coordinate after organoid is aligned
  double rotated_phi;
  double tilt_angle; // in range [0, 180); also know as beta in our organoid paper
  Patch(unsigned id_)
  {
    id = id_;
    triangles.reserve(20);
    area=0;
    unit_normal_vec = Vector3d(0,0,0);
    com = Vector3d(0,0,0);
    Q_3D = Matrix3d::Zero();
    C_int = Matrix3d::Zero();
    Q_eigvec = Vector3d(0,0,0);
    Q_norm = 0.;
    rotated_theta = 0.;
    rotated_phi = 0.;
    tilt_angle = 0.;
  }
};


class Triangulation
{
public:
  // ---- CREATION & DESTRUCTION ----
  Triangulation();
  Triangulation(bool);
  // Vertex(bool is_3D);      //constructor
  virtual ~Triangulation();     //destructor

  // --- GEOMETRY ---
  std::vector<TriVertex> Vertices;
  std::vector<TriEdge> Edges;
  std::vector<Triangle> Triangles;
  std::vector<Patch> Patches;
  bool VSH_analysis;
  // ceep a copy of vertex coordinarted in spherical coords
  //std::vector<Vector3d> Vertices_sph;

  TriVertex* append_vertex(Vector3d pos_);
  TriEdge* append_edge(TriVertex *tail_, TriVertex *head_);
  void update_Vertices_sph();
  void update_Vertices_vec_field(); // v->pos - v->rHat
  void update_Vertices_vec_field(std::vector<Vector3d>);
  //void update_Vertices_vec_field(std::string vec_filed_name); // loading vec_field from a file
  void update_Edges_geometry();
  void update_Triangles_area_and_normal_vec();
  void update_Triangles_geomtery();
  void update_all_geometries();
  void update_all_geometries(std::vector<Vector3d>);
  void update_all_triangles_integrated_curvature_tensor();
  void update_all_Q_3D_xyplane();
  std::tuple<double, double, double> get_Elm_coefficients(int l, int m);//return real E^r, E^1, E^2 for mode l,m
  // returning both real and complex modes
  std::tuple<double, double, double, std::complex<double>, std::complex<double>, std::complex<double> > get_Elm_Clm_coefficients(int l, int m);
  VSH_MODES_REAL get_VectorSphericalHarmonicsModes(int Lmax_);
  VSH_MODES_REAL_COMPLEX get_VectorSphericalHarmonicsModes_real_and_complex(int Lmax_);
  std::vector<Vector3d> reconstruct_vec_field_from_modes(std::vector<double> Er, std::vector<double> E1, std::vector<double> E2, int Lmax_, bool add_rhat_);
  std::vector<Vector3d> get_purely_rotational_part_from_E2(std::vector<double> E2);
  Matrix3d compute_Q3d_Patch(unsigned patch_id_);
  void update_all_Patches_Q_3D();
  void update_all_elongation_tilt_angles_if_already_aligned();
  //Vector3d compute_Qeigvec_Patch(unsigned patch_id_);
  Matrix3d compute_Cint_Patch(unsigned patch_id_);
  void update_all_Patches_Cint();
};

#endif//TRIANGULATION_HPP_
