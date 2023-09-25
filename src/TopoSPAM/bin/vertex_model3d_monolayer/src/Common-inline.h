#ifndef COMMON_INLINE_H
#define COMMON_INLINE_H

#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

using Eigen::Matrix3d;
using Eigen::Vector3d;

#ifndef NORM_EPS
#define NORM_EPS 1e-14
#endif


// double inline SIGN(double a,double b){return (b>0 ? std::abs(a):-1.0*std::abs(a));}
// inline void shft3(double &a, double &b, double &c, const double d)
// {
// a=b;
// b=c;
// c=d;
// }
// inline void SWAP(double &x, double &y) { double temp=x; x=y; y=temp; }

template<typename T>
inline double sign_func(T x)
{
  if(x>0) return +1.;
  else if(x<0) return -1.;
  else return 0.;
}


inline Vector3d cartesian_to_spherical(Vector3d vec)
{
  double xy = vec[0]*vec[0] + vec[1]*vec[1];
  double rr = vec.norm();
  double tta = std::atan2(std::sqrt(xy), vec[2]);
  double phi = std::atan2(vec[1], vec[0]);
  // std::atan2(std::hypot(vec[0],vec[1]), vec[1]);
  return Vector3d(rr, tta, phi);
}

inline std::tuple<Vector3d,Vector3d,Vector3d> get_spherical_coord_unit_vectors(double theta, double phi)
{
  double st = std::sin(theta); double ct = std::cos(theta);
  double sp = std::sin(phi); double cp = std::cos(phi);
  Vector3d r_hat(st*cp, st*sp, ct);
  Vector3d theta_hat(cp*ct, sp*ct, -st);
  Vector3d phi_hat(-sp, cp, 0.);
  return std::make_tuple(r_hat, theta_hat, phi_hat);
}

inline std::tuple<double,double,double> get_tilt_angle(Vector3d q_vec, Vector3d cell_cent)
{
  Vector3d q_vec_unit = q_vec/q_vec.norm();
  Vector3d cent_sph_coord = cartesian_to_spherical(cell_cent);
  auto [r_hat, thet_hat, phi_hat] = get_spherical_coord_unit_vectors(cent_sph_coord[1], cent_sph_coord[2]);
  (void) r_hat;
  if(std::acos(q_vec_unit.dot(thet_hat)) < M_PI/2.) q_vec_unit = -q_vec_unit;
  return std::make_tuple(std::acos(q_vec_unit.dot(phi_hat)), cent_sph_coord[1], cent_sph_coord[2]);
}

inline std::pair<Vector3d, Matrix3d> get_ang_momentum_inertia_tensor(std::vector<Vector3d> pos_t1, std::vector<Vector3d> pos_t2, double dt)
{
  int n_points = pos_t1.size();
  std::vector<Vector3d> velocity; velocity.reserve(n_points);
  Vector3d RR(0,0,0), VV(0,0,0);
  for(int i=0; i<n_points; i++)
  {
    RR += pos_t1[i];
    Vector3d vi = (pos_t2[i]-pos_t1[i])/dt;
    velocity.emplace_back(vi);
    VV += vi;
  }
  RR = RR/n_points; VV = VV/n_points;
  //------------
  Vector3d Gamma(0,0,0); Matrix3d MM  = Matrix3d::Zero();
  Matrix3d II = Matrix3d::Identity();
  for(int i=0; i<n_points; i++)
  {
    Vector3d dr_i = pos_t1[i] - RR;
    Vector3d dv_i = velocity[i] - VV;
    Gamma += dr_i.cross(dv_i);
    Matrix3d M_i = II* dr_i.dot(dr_i) - dr_i * dr_i.transpose();
    MM += M_i;
  }
  Gamma = Gamma/n_points; MM = MM/n_points;
  //std::cout << MM/dt << "\n\n" << std::endl;
  return std::make_pair(Gamma, MM);
}

inline Vector3d project_on_tangent_plane(Vector3d nn, Vector3d p)
{
  /*  projecting a vector p onto a tangent plane
    defined by normal vector nn */
  nn = nn/nn.norm();
  Matrix3d II = Matrix3d::Identity();
  Matrix3d nxn = nn * nn.transpose();
  Matrix3d proj_op = II - nxn;
  Vector3d tt = proj_op * p;
  return tt;
}


inline Vector3d rotation(Vector3d vec, double angle, Vector3d axis)
{
  /*    Return the vector 'vec' rotated by an angle 'angle' around the axis 'axis' in 3d
    and return also the associated rotation matrix 'rot'.
    The computation goes as follow:
    let 'u' be the unit vector along 'axis', i.e. u = axis/norm(axis)
    and A = I × u be the skew-symmetric matrix associated to 'u', i.e.
    the cross product of the identity matrix with 'u'
    then rot = exp(tta A) is the 3d rotation matrix.*/
  if(angle<NORM_EPS || axis.norm()<NORM_EPS) return vec;

  Vector3d unit_axis = axis/axis.norm();
  Vector3d n1(1.,0,0), n2(0,1.,0), n3(0,0,1.);
  Vector3d IU1 = n1.cross(unit_axis*angle) ;
  Vector3d IU2 = n2.cross(unit_axis*angle) ;
  Vector3d IU3 = n3.cross(unit_axis*angle) ;
  Matrix3d ttaA;
  ttaA(0,0) = IU1(0); ttaA(0,1) = IU1(1); ttaA(0,2) = IU1(2);
  ttaA(1,0) = IU2(0); ttaA(1,1) = IU2(1); ttaA(1,2) = IU2(2);
  ttaA(2,0) = IU3(0); ttaA(2,1) = IU3(1); ttaA(2,2) = IU3(2);
  Matrix3d rot = ttaA.exp();
  Vector3d rotated_vec = rot*vec;
  return rotated_vec;
  // Vector3d rotated_vec;
  // rotated_vec(0) = rot(0,0)*vec(0) + rot(0,1)*vec(1) + rot(0,2)*vec(2);
  // rotated_vec(1) = rot(1,0)*vec(0) + rot(1,1)*vec(1) + rot(1,2)*vec(2);
  // rotated_vec(2) = rot(2,0)*vec(0) + rot(2,1)*vec(1) + rot(2,2)*vec(2);
  // return rotated_vec;
}

inline Matrix3d rotation(Matrix3d tens, double angle, Vector3d axis)
{
  /*    Return the vector 'vec' rotated by an angle 'angle' around the axis 'axis' in 3d
    and return also the associated rotation matrix 'rot'.
    The computation goes as follow:
    let 'u' be the unit vector along 'axis', i.e. u = axis/norm(axis)
    and A = I × u be the skew-symmetric matrix associated to 'u', i.e.
    the cross product of the identity matrix with 'u'
    then rot = exp(tta A) is the 3d rotation matrix.*/
  if(angle<NORM_EPS || axis.norm()<NORM_EPS) return tens;

  Vector3d unit_axis = axis/axis.norm();
  Vector3d n1(1.,0,0), n2(0,1.,0), n3(0,0,1.);
  Vector3d IU1 = n1.cross(unit_axis*angle) ;
  Vector3d IU2 = n2.cross(unit_axis*angle) ;
  Vector3d IU3 = n3.cross(unit_axis*angle) ;
  Matrix3d ttaA;
  ttaA(0,0) = IU1(0); ttaA(0,1) = IU1(1); ttaA(0,2) = IU1(2);
  ttaA(1,0) = IU2(0); ttaA(1,1) = IU2(1); ttaA(1,2) = IU2(2);
  ttaA(2,0) = IU3(0); ttaA(2,1) = IU3(1); ttaA(2,2) = IU3(2);
  Matrix3d rot = ttaA.exp();
  return rot*tens*rot.transpose();
}

inline double signed_terahedron_volume(Vector3d a, Vector3d b, Vector3d d, Vector3d cent)
{
  // treat cent as centroid of the closed surface (i.e. a cell)
  // b is before a, in CCW order, and d is the face center
  Vector3d cs = (a-b).cross(d-b);
  return cs.dot(b-cent)/6.;
}

inline Vector3d gradient_face_area(Vector3d a, Vector3d b, Vector3d c, Vector3d Mf, double Nv)
{
  Vector3d tt = (1. - Nv) * b + Mf;
  Vector3d uu = (1. - Nv) * c + Mf;
  Vector3d rr = Mf.cross(b) + a.cross(tt);
  Vector3d ss = Mf.cross(c) + a.cross(uu);
  Vector3d grad_A = tt.cross(rr)/(2.*Nv*rr.norm()+NORM_EPS) + uu.cross(ss)/(2.*Nv*ss.norm()+NORM_EPS); // gradient of area
  return grad_A;
}

inline Vector3d gradient_volume(Vector3d b, Vector3d c, Vector3d M_cell, Vector3d Mf, double Nvc, double Nvf)
{
  Vector3d kk = (2.*Nvc-1.)*Mf - (Nvf-1.)*M_cell;
  return (c-b).cross(kk)/(12.*Nvc*Nvf);
}

inline std::pair<std::vector<Vector3d>, std::vector<std::vector<int>>> dodecahedron_verts_faces()
{
  // a regular dodecahedron centered at origin; it has 20 vertices and 12 pentagonal faces
  std::vector<Vector3d> verts{
  Vector3d(-1.3763819204711736, 0., 0.2628655560595668),
  Vector3d(1.3763819204711736, 0., -0.2628655560595668),
  Vector3d(-0.42532540417601994, -1.3090169943749475, 0.2628655560595668),
  Vector3d(-0.42532540417601994, 1.3090169943749475, 0.2628655560595668),
  Vector3d(1.1135163644116066, -0.8090169943749475, 0.2628655560595668),
  Vector3d(1.1135163644116066, 0.8090169943749475, 0.2628655560595668),
  Vector3d(-0.2628655560595668, -0.8090169943749475, 1.1135163644116066),
  Vector3d(-0.2628655560595668, 0.8090169943749475, 1.1135163644116066),
  Vector3d(-0.6881909602355868, -0.5, -1.1135163644116068),
  Vector3d(-0.6881909602355868, 0.5, -1.1135163644116068),
  Vector3d(0.6881909602355868, -0.5, 1.1135163644116066),
  Vector3d(0.6881909602355868, 0.5, 1.1135163644116066),
  Vector3d(0.85065080835204, 0., -1.1135163644116068),
  Vector3d(-1.1135163644116068, -0.8090169943749475, -0.2628655560595668),
  Vector3d(-1.1135163644116068, 0.8090169943749475, -0.2628655560595668),
  Vector3d(-0.8506508083520399, 0., 1.1135163644116066),
  Vector3d(0.2628655560595668, -0.8090169943749475, -1.1135163644116068),
  Vector3d(0.2628655560595668, 0.8090169943749475, -1.1135163644116068),
  Vector3d(0.42532540417601994, -1.3090169943749475, -0.2628655560595668),
  Vector3d(0.42532540417601994, 1.3090169943749475, -0.2628655560595668)};

  std::vector<std::vector<int>> face_inds{
    {14, 9, 8, 13, 0},
    {1, 5, 11, 10, 4},
    {4, 10, 6, 2, 18},
    {10, 11, 7, 15, 6},
    {11, 5, 19, 3, 7},
    {5, 1, 12, 17, 19},
    {1, 4, 18, 16, 12},
    {3, 19, 17, 9, 14},
    {17, 12, 16, 8, 9},
    {16, 18, 2, 13, 8},
    {2, 6, 15, 0, 13},
    {15, 7, 3, 14, 0}};

  return std::make_pair(verts, face_inds);
}

#endif //ifndef
