#ifndef BOND_INLINE_H
#define BOND_INLINE_H
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "tools.hpp"

using Eigen::Matrix2d;
using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::Vector2d;
using Eigen::EigenSolver;

#ifndef NORM_EPS
#define NORM_EPS 1e-15
#endif

inline std::pair<Vector3d,Vector3d> triangle_side_vectors(Vector3d v1, Vector3d v2, Vector3d v3)
{
  Vector3d e12 = v2-v1;
  Vector3d e13 = v3-v1;
  return std::make_pair(e12, e13);
}

inline std::tuple<Vector3d,Vector3d,double> triangle_normal_and_area(Vector3d e12, Vector3d e13)
{
  Vector3d normal_vec = e12.cross(e13);
  double normal_vec_norm = normal_vec.norm();
  return std::make_tuple(normal_vec, normal_vec/(normal_vec_norm+NORM_EPS), normal_vec_norm/2.);
}

inline std::pair<double,double> normal_ThetaX_ThetaY(Vector3d unit_normal)
{
  //Express the orientation of the normal vector in two angles, for a rotation around the x and y axis.
  //These angles is how much to rotate along the respective axis to rotate the z unit vector to the triangle normal.
  //definition: Define Theta_x, Theta_y such that: unit_N = RotX_ThetaX . RotY_ThetaY . unit_z
  // Rot_x and Rot_y are the rotation matrices along the x- and y- axis and z is the unit normal along the z-axis.
  double normal_ThetaX = - std::atan2(unit_normal[1],unit_normal[2]);
  double normal_ThetaY = std::atan2(unit_normal[0], std::sqrt(1.-unit_normal[0]*unit_normal[0]));

  return std::make_pair(normal_ThetaX,normal_ThetaY);
}

inline Matrix3d R_x(double rot_ang_X)
{
  double cosX = std::cos(rot_ang_X);
  double sinX = std::sin(rot_ang_X);
  Matrix3d rotation_x;
  rotation_x(0,0)=1.; rotation_x(0,1)=0.; rotation_x(0,2)=0.;
  rotation_x(1,0)=0.; rotation_x(1,1)=cosX; rotation_x(1,2)=-sinX;
  rotation_x(2,0)=0.; rotation_x(2,1)=sinX; rotation_x(2,2)=cosX;
  return rotation_x;
}

inline Matrix3d R_y(double rot_ang_Y)
{
  double cosY = std::cos(rot_ang_Y);
  double sinY = std::sin(rot_ang_Y);
  Matrix3d rotation_y;
  rotation_y(0,0)=cosY; rotation_y(0,1)=0.; rotation_y(0,2)=sinY;
  rotation_y(1,0)=0.; rotation_y(1,1)=1.; rotation_y(1,2)=0.;
  rotation_y(2,0)=-sinY; rotation_y(2,1)=0.; rotation_y(2,2)=cosY;
  return rotation_y;
}

inline Matrix3d R_z(double rot_ang_Z)
{
  double cs = std::cos(rot_ang_Z);
  double sn = std::sin(rot_ang_Z);
  Matrix3d rotation_z;
  rotation_z(0,0)=cs; rotation_z(0,1)=-sn; rotation_z(0,2)=0.;
  rotation_z(1,0)=sn; rotation_z(1,1)=cs;  rotation_z(1,2)=0.;
  rotation_z(2,0)=0.; rotation_z(2,1)=0.;  rotation_z(2,2)=1.;
  return rotation_z;
}


inline Matrix3d get_inverse_Ry_Rx(double rot_ang_X, double rot_ang_Y)
{
  // Rotate the triangle vectors into the XY plain.
  //Minus signs for angles because you rotate normal vector towards z-unit vector.
  Matrix3d inv_rotation_x = R_x(-rot_ang_X);
  Matrix3d inv_rotation_y = R_y(-rot_ang_Y);
  //Define (Rot_x . Rot_y)^-1 rotation matrix.
  return inv_rotation_y * inv_rotation_x;
}

inline std::pair<Vector2d,Vector2d> rotate_triangles_into_xy_plane(double rot_ang_X, double rot_ang_Y, Vector3d e12, Vector3d e13)
{
  // Rotate the triangle vectors into the XY plain.
  //Minus signs for angles because you rotate normal vector towards z-unit vector.
  Matrix3d inv_Ry_Rx = get_inverse_Ry_Rx(rot_ang_X, rot_ang_Y);
  Vector3d e12_xy_3d = inv_Ry_Rx * e12;
  Vector3d e13_xy_3d = inv_Ry_Rx * e13;
  return std::make_pair(Vector2d(e12_xy_3d[0],e12_xy_3d[1]), Vector2d(e13_xy_3d[0],e13_xy_3d[1]));
}

//Split a given tensor T into trace, anti-symmetric part and a nematic part: T = A + B_nematic
//Return angles and norms of both tensors.
inline std::tuple<double,double,double,double> splitTensor_angleNorm(Matrix2d T)
{
  // T = A + B
  // A = [[a,-b],
  //      [b, a]]
  // B = [[c, d],
  //      [d,-c]]
  double a = 0.5 * (T(0, 0) + T(1, 1));
  double b = 0.5 * (T(1, 0) - T(0, 1));
  double c = 0.5 * (T(0, 0) - T(1, 1));
  double d = 0.5 * (T(0, 1) + T(1, 0));

  double theta = std::atan2(b, a);
  double AabsSq = a * a + b * b;

  double twophi = theta + std::atan2(d, c);
  double BabsSq = c * c + d * d;

  return std::make_tuple(theta, AabsSq, twophi, BabsSq);
}

// decompose the triangle state tensor in x-y plane into anti-symmetric part and nematic parts
//Return angles and norms of both tensors.
inline std::tuple<double,double,double,double> calculate_triangle_state_tensor_symmetric_antisymmetric(Vector2d e12_xy, Vector2d e13_xy)
{
  double A0=1.; double l = std::sqrt(4. * A0 /std::sqrt(3.));
  Matrix2d equilat_tri;
  equilat_tri(0,0) = l; equilat_tri(0,1)=l/2.;
  equilat_tri(1,0) =0.; equilat_tri(1,1)=std::sqrt(3.)*l/2.;
  Matrix2d inv_equilat_tri = equilat_tri.inverse();
  Matrix2d tri_state_R;
  tri_state_R(0,0) = e12_xy[0]; tri_state_R(0,1) = e13_xy[0];
  tri_state_R(1,0) = e12_xy[1]; tri_state_R(1,1) = e13_xy[1];
  Matrix2d tri_state_S = tri_state_R * inv_equilat_tri;
  // double theta, AabsSq, twophi, BabsSq;
  return splitTensor_angleNorm(tri_state_S);
}

// Given the norm of the symmetric and anti-symmetric part of the
// decomposed state tensor, calculate the elongation norm.
// For zero-area triangles, set elongation to high value 100.
inline double Qnorm(double AabsSq, double BabsSq)
{
  if ( (AabsSq - BabsSq)<NORM_EPS )
  {
    return 0.;
  }
  return std::asinh( std::sqrt(BabsSq)/std::sqrt(AabsSq - BabsSq) );
}


inline Matrix3d calculate_triangle_elongation_tensor_xy_plane(double rot_ang_X, double rot_ang_Y, Vector3d e12, Vector3d e13)
{
  // Generate the 2d elongation tensors in the xy-plane, written in 3d matrix form.
  Vector2d e12_xy, e13_xy;
  std::tie(e12_xy, e13_xy) = rotate_triangles_into_xy_plane(rot_ang_X, rot_ang_Y, e12, e13);

  double rotation_angle_z, AabsSq, twophi, BabsSq;
  std::tie(rotation_angle_z, AabsSq, twophi, BabsSq) =  calculate_triangle_state_tensor_symmetric_antisymmetric(e12_xy, e13_xy);

  double triangle_elongation_norm = Qnorm(AabsSq, BabsSq);
  // calculate the components of the elongation tensor qxx and qxy.
  double cos_two_phi = std::cos(twophi);
  double sin_two_phi = std::sin(twophi);
  double qxx = triangle_elongation_norm * cos_two_phi;
  double qxy = triangle_elongation_norm * sin_two_phi;

  // Generate the 2d elongation tensors in the xy-plane, written in 3d matrix form.
  Matrix3d triangle_q_2d;
  triangle_q_2d(0,0) = qxx; triangle_q_2d(0,1) = qxy;  triangle_q_2d(0,2) = 0.;
  triangle_q_2d(1,0) = qxy; triangle_q_2d(1,1) = -qxx; triangle_q_2d(1,2) = 0.;
  triangle_q_2d(2,0) = 0.;  triangle_q_2d(2,1) = 0.;   triangle_q_2d(2,2) = 0.;

  return triangle_q_2d;
}



//Given an array of 3d matrices describing a quantity in the xy-plane,
//function returns matrix transformed to the plane of the triangle.
inline Matrix3d transform_matrices_from_xy_plane_to_triangle_plane(Matrix3d triangle_q_2d, double rot_ang_X, double rot_ang_Y)
{
  /*
  Rotate triangle_q_2d such that the normal vector coincides with the triangle normal.

  Parameters
  ----------
  */
  // Rotate the triangle vectors from the XY plane back to it's plane .
  Matrix3d rotation_x = R_x(rot_ang_X);
  Matrix3d rotation_y = R_y(rot_ang_Y);

  Matrix3d RotX_RotY = rotation_x * rotation_y;
  Matrix3d inverse_RotX_RotY = RotX_RotY.inverse();

  Matrix3d triangle_q_3d = RotX_RotY * (triangle_q_2d * inverse_RotX_RotY);
  //Matrix3d triangle_q_3d = inverse_RotX_RotY * (triangle_q_2d * RotX_RotY);

  return triangle_q_3d;
}


inline Matrix3d calculate_triangle_elongation_tensor(double rot_ang_X, double rot_ang_Y, Vector3d e12, Vector3d e13)
{
  // Generate the 2d elongation tensors in the xy-plane, written in 3d matrix form.
  Matrix3d triangle_q_2d = calculate_triangle_elongation_tensor_xy_plane(rot_ang_X, rot_ang_Y, e12, e13);

  Matrix3d triangle_q_3d = transform_matrices_from_xy_plane_to_triangle_plane(triangle_q_2d, rot_ang_X, rot_ang_Y);

  return triangle_q_3d;
}

inline std::tuple<Matrix3d, Vector3d, double> calculate_triangle_elongation_tensor_xy_plane(Vector3d v1, Vector3d v2, Vector3d v3)
{
  Vector3d e12, e13;
  std::tie(e12, e13) = triangle_side_vectors(v1, v2, v3);
  Vector3d normal, unit_normal;
  double area;
  std::tie(normal, unit_normal,area) = triangle_normal_and_area(e12, e13);
  double rot_ang_X, rot_ang_Y;
  std::tie(rot_ang_X, rot_ang_Y) = normal_ThetaX_ThetaY(unit_normal);

  Matrix3d triangle_q_2d = calculate_triangle_elongation_tensor_xy_plane(rot_ang_X, rot_ang_Y, e12, e13);

  return std::make_tuple(triangle_q_2d, unit_normal, area);
}

inline std::tuple<Matrix3d, Vector3d, double> calculate_triangle_elongation_tensor(Vector3d v1, Vector3d v2, Vector3d v3)
{
  Vector3d e12, e13;
  std::tie(e12, e13) = triangle_side_vectors(v1, v2, v3);
  Vector3d normal, unit_normal;
  double area;
  std::tie(normal, unit_normal,area) = triangle_normal_and_area(e12, e13);
  double rot_ang_X, rot_ang_Y;
  std::tie(rot_ang_X, rot_ang_Y) = normal_ThetaX_ThetaY(unit_normal);

  Matrix3d triangle_q_3d = calculate_triangle_elongation_tensor(rot_ang_X, rot_ang_Y, e12, e13);

  return std::make_tuple(triangle_q_3d, unit_normal, area);
}

/*
- Get tensor eigenvalues and eigenvectors with specified return_type.
- return_type=1: return eigenvalue/vec with highest value
- return_type=0      : return eigenvalue/vec with smallest norm
*/
inline std::pair<double,Vector3d> calculate_3x3_tensor_eigenvalues_and_eigenvector(Matrix3d tens, int return_type)
{
  EigenSolver<Matrix3d> es(tens);
  Vector3d eig_vals = es.eigenvalues().real();
  Matrix3d eig_vecs = es.eigenvectors().real();
  double val=0; Vector3d vec(0,0,0);
  if(return_type==1)
  {
    Eigen::Index max_index;
    eig_vals.maxCoeff(&max_index);
    val = eig_vals[max_index];
    vec = eig_vecs.col(max_index);
  }
  else if(return_type==0)
  {
    Eigen::Index min_abs_index;
    eig_vals.array().abs().minCoeff(&min_abs_index);
    val = eig_vals[min_abs_index];
    vec = eig_vecs.col(min_abs_index);
  }
  else
  {
    std::cout << "must specify return type: 1=eigen vec for max eigenvalue;     0=eigen vec for min eigenvalue;" << std::endl;
  }
  return std::make_pair(val,vec);
}

/*
Rotate rank-2 tensors into the xy plane, such that the eigenvector with eigenvalue closest to zero is paralell with the system z-axis.

Parameters:
-----------
nematic_tensors : Matrix3d
    rank2 tensor  on a triangle
surface_normal : Vector3d
    Normal vector on a triangle (Required to determine orientation of surface).

Returns:
--------
nematic_tensors_inplane : Matrix3d
    For every triangle a 3x3 matrix, representing the rank 2 tensor rotated into the xy-plane.
    */
inline Matrix3d rotate_tensors_in_plane(Matrix3d tens, Vector3d surface_normal)
{
  Vector3d nematic_tensors_normal;
  std::tie(std::ignore, nematic_tensors_normal) = calculate_3x3_tensor_eigenvalues_and_eigenvector(tens, 0);//eigenvector with eigenvalue closest to zero
  double nematic_tensors_normal_orientation = sign_func(surface_normal.dot(nematic_tensors_normal));
  nematic_tensors_normal *= nematic_tensors_normal_orientation;

  double nematic_tensors_to_z_angle = std::acos(nematic_tensors_normal[2]);
  Vector3d nematic_tensors_to_z_rotdir = nematic_tensors_normal.cross(Vector3d(0.,0.,1.));
  nematic_tensors_to_z_rotdir /= (nematic_tensors_to_z_rotdir.norm()+NORM_EPS);
  return rotation(tens, nematic_tensors_to_z_angle, nematic_tensors_to_z_rotdir);
}

#endif //ifndef
