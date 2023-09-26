#include "Triangle.hpp"

Triangle::Triangle()
{
  polarity = Vector3d(0,0,0);
  v1->connected_triangles.emplace_back(this);
  v2->connected_triangles.emplace_back(this);
  v3->connected_triangles.emplace_back(this);
  Q_xyplane = Matrix3d::Zero();
  C_int = Matrix3d::Zero();
}

//Triangle(TriVertex *v1_, TriVertex *v2_, TriVertex *v3_, int id_);
Triangle::Triangle(TriVertex *v1_, TriVertex *v2_, TriVertex *v3_, TriEdge *e1_, TriEdge *e2_, TriEdge *e3_, int id_) :
        v1(v1_), v2(v2_), v3(v3_), e1(e1_), e2(e2_), e3(e3_), id(id_)
      {
        polarity = Vector3d(0,0,0);
        v1->connected_triangles.emplace_back(this);
        v2->connected_triangles.emplace_back(this);
        v3->connected_triangles.emplace_back(this);
        Q_xyplane = Matrix3d::Zero();
        C_int = Matrix3d::Zero();
      }
// void Vertex::setNetwork(Network *n) {
//   _net = n;
// }

//Vertex destructor
Triangle::~Triangle(){
}

double Triangle::update_excess()
{
  //double m_pi = 3.141592653589793238462643383279502884197169399;
  Vector3d ra = v1->pos/(v1->pos).norm();//Point(sin(v1->sphCoord.y)*cos(v1->sphCoord.z), sin(v1->sphCoord.y)*sin(v1->sphCoord.z), v1->cosTheta);
  Vector3d rb = v2->pos/(v2->pos).norm();//Point(sin(v2->sphCoord.y)*cos(v2->sphCoord.z), sin(v2->sphCoord.y)*sin(v2->sphCoord.z), v2->cosTheta);
  Vector3d rc = v3->pos/(v3->pos).norm();//Point(sin(v3->sphCoord.y)*cos(v3->sphCoord.z), sin(v3->sphCoord.y)*sin(v3->sphCoord.z), v3->cosTheta);

  double cosa = rb.dot(rc); double cosb = ra.dot(rc); double cosc = ra.dot(rb);
  double sina = std::sqrt(1.0-cosa*cosa); double sinb = std::sqrt(1.0-cosb*cosb); double sinc = std::sqrt(1.0-cosc*cosc);
  double A = std::acos( (cosa-cosb*cosc)/(sinb*sinc) );
  double B = std::acos( (cosb-cosa*cosc)/(sina*sinc) );
  double C = std::acos( (cosc-cosa*cosb)/(sina*sinb) );
  excess = (A+B+C-M_PI);
  return excess;
}

void Triangle::update_corner_angles()
{
  Vector3d unit_l1 = e1->l_vec/(e1->l+NORM_EPS);
  Vector3d unit_l2 = e2->l_vec/(e2->l+NORM_EPS);
  Vector3d unit_l3 = e3->l_vec/(e3->l+NORM_EPS);
  angle_1 = std::acos(-unit_l1.dot(unit_l3));
  angle_2 = std::acos(-unit_l1.dot(unit_l2));
  angle_3 = std::acos(-unit_l2.dot(unit_l3));
}
// double Triangle::update_excess()
// {
//   //double m_pi = 3.141592653589793238462643383279502884197169399;
//   Point ra = Point(sin(v1->sphCoord.y)*cos(v1->sphCoord.z), sin(v1->sphCoord.y)*sin(v1->sphCoord.z), v1->cosTheta);
//   Point rb = Point(sin(v2->sphCoord.y)*cos(v2->sphCoord.z), sin(v2->sphCoord.y)*sin(v2->sphCoord.z), v2->cosTheta);
//   Point rc = Point(sin(v3->sphCoord.y)*cos(v3->sphCoord.z), sin(v3->sphCoord.y)*sin(v3->sphCoord.z), v3->cosTheta);
//
//   double cosa = rb.dot(rc); double cosb = ra.dot(rc); double cosc = ra.dot(rb);
//   double sina = sqrt(1.0-cosa*cosa); double sinb = sqrt(1.0-cosb*cosb); double sinc = sqrt(1.0-cosc*cosc);
//   double A = acos( (cosa-cosb*cosc)/(sinb*sinc) );
//   double B = acos( (cosb-cosa*cosc)/(sina*sinc) );
//   double C = acos( (cosc-cosa*cosb)/(sina*sinb) );
//   excess = (A+B+C-PI);
//   return excess;
// }
