// #include <type_traits>
// #include <typeinfo>
//#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include "../accessories/include/json.hpp"
#include "Tissue.hpp"
#include <functional>
#include <gsl/gsl_multimin.h>

//#define assertm(exp, msg) assert(((void)msg, exp))

// template <typename T> std::string type_name();

// template< typename F >
//   class gsl_function_pp : public gsl_function {
//   public:
//   gsl_function_pp(const F& func) : _func(func) {
//     function = &gsl_function_pp::invoke;
//     params=this;
//   }
//   private:
//   // const F& _func;
//   const F& _func;
//   static double invoke(const gsl_vector * x, void *params) {
//     return static_cast<gsl_function_pp*>(params)->_func(x);
//   }
// };

typedef double ( * GSLMultiMinFuncPointer ) ( const gsl_vector *, void *);
typedef void   ( * GSLMultiMinDfPointer )   ( const gsl_vector *, void *, gsl_vector *);
typedef void   ( * GSLMultiMinFdfPointer ) ( const gsl_vector *, void *, double *, gsl_vector *);

Tissue::Tissue()
{
  ListVertex.reserve(2000);
  ListEdge.reserve(6000);
  ListCell.reserve(1000);
}

Tissue::Tissue(int N0) // this argument helps with a reasonable estimation of required memory
{
  ListVertex.reserve(4*N0);
  ListEdge.reserve(12*N0);
  ListCell.reserve(2*N0);
}

//Tissue destructor
Tissue::~Tissue(){
  ListVertex.clear();
  ListEdge.clear();
  ListCell.clear();
}

void Tissue::set_CellMechProps3D(double Kc_, double V0c_, double Ta_, double Tl_, double Tb_)
{
  CellMechProp.Kc = Kc_;
  CellMechProp.V0c = V0c_;
  CellMechProp.Ta = Ta_;
  CellMechProp.Tl = Tl_;
  CellMechProp.Tb = Tb_;
}

void Tissue::set_SpheroidMechProp(double P0_, double Kl_, double V0l_)
{
  SphMechProp.P0 = P0_;
  SphMechProp.Kl = Kl_;
  SphMechProp.V0l = V0l_;
}

void Tissue::set_CellMechProps2D(double Kc_, double A0c_, double BondT_, double perim_elasticity_)
{
  CellMechProp.Kc = Kc_;
  CellMechProp.A0c = A0c_;
  CellMechProp.BondT = BondT_;
  CellMechProp.perim_elasticity = perim_elasticity_;
}

// void Tissue::set_minimization_params()
// {
//   minPar.MinimizationMode = 2;
//   minPar.ftol = 1e-12;
//   minPar.GTOL = 1e-8;
//   minPar.rateT1check = 0;
//   minPar.tol = 1e-9;
//   minPar.ITMAX = 10000;
//   minPar.EPS = 1e-18;
// }

Vertex* Tissue::similar_or_new_vertex_2D(Vector3d pos)
{
  // iterate through all the vertices in the system and check if the vertex has the same coordinates
  std::vector<Vertex>::iterator it = ListVertex.begin();
  for(;it!=ListVertex.end();it++)
  {
    Vertex* v = &(*it);
    if((v->pos_a-pos).norm() < 1e-10) return v;
  }
  //if it makes it this far, the vertes does not exist, and needs to be made
  ListVertex.emplace_back(Vertex(pos, pos, ListVertex.size()));

  return &(ListVertex.back());
}

Vertex* Tissue::similar_or_new_vertex_3D(Vector3d _pos_a, Vector3d _pos_b)
{
  // iterate through all the vertices in the system and check if the vertex has the same coordinates
  std::vector<Vertex>::iterator it = ListVertex.begin();
  for(;it!=ListVertex.end();it++)
  {
    Vertex* v = &(*it);
    if((v->pos_a-_pos_a).norm() < 1e-10 && (v->pos_b-_pos_b).norm() < 1e-10) return v;
  }
  //if it makes it this far, the vertes does not exist, and needs to be made
  ListVertex.emplace_back(Vertex(_pos_a, _pos_b, ListVertex.size()));

  return &(ListVertex.back());
}

Edge* Tissue::similar_or_new_edge(Vertex* v_tail, Vertex* v_head)
{
  Edge *e_conj = nullptr;
  // iterate through all the edges in the system and check if another Edge has the same vertices
  std::vector<Edge>::iterator it = ListEdge.begin();
  for(;it!=ListEdge.end();it++)
  {
    Edge* e = &(*it);
    if(e->tail == v_tail && e->head == v_head) return e;
    else if(e->head == v_tail && e->tail == v_head) e_conj = e;
  }
  //if it makes it this far, the edge does not exist, and needs to be made
  ListEdge.emplace_back(Edge(v_tail, v_head, e_conj, ListEdge.size()));
  return &(ListEdge.back());
}

Vector3d Tissue::get_cell_centroid(Cell *c)
{
  c->centroid = Vector3d(0,0,0);
  for(auto v : c->Vlist){c->centroid += (v->pos_a + v->pos_b);}

  c->centroid /= (2. * c->Vlist.size());

  return c->centroid;
}

double Tissue::get_cell_volume(Cell *c)
{
  c->vol = 0.;
  for(auto e : c->Elist) // 6 tetrahedrons per each edge
  {
    // sum of 1 apical, a basal and 4 latera tehtrahedron volumes
    c->vol += signed_terahedron_volume(e->head->pos_a, e->tail->pos_a, c->cent_a, c->centroid);//apical
    c->vol += signed_terahedron_volume(e->tail->pos_b, e->head->pos_b, c->cent_b, c->centroid);//basal
    c->vol += signed_terahedron_volume(e->head->pos_b, e->tail->pos_b, e->cent_l, c->centroid); //lateral
    c->vol += signed_terahedron_volume(e->head->pos_a, e->head->pos_b, e->cent_l, c->centroid); //lateral
    c->vol += signed_terahedron_volume(e->tail->pos_a, e->head->pos_a, e->cent_l, c->centroid); //lateral
    c->vol += signed_terahedron_volume(e->tail->pos_b, e->tail->pos_a, e->cent_l, c->centroid); //lateral
    //std::cout << c->vol << std::endl;
  }
  return c->vol;
}

double Tissue::get_lumen_volume()
{
  lumen_volume =0.;

  #if PRALLEL
  auto func = [this](auto&& it){this->lumen_volume += signed_terahedron_volume((&it)->tail->pos_a, (&it)->head->pos_a, (&it)->cell->cent_a, Vector3d(0.,0.,0.) ); };
  // __gnu_parallel::for_each(
  //   ListEdge.begin(),
  //   ListEdge.end(),
  //   func);
  std::for_each(
    std::execution::par_unseq,
    ListEdge.begin(),
    ListEdge.end(),
    func);
  #else
  Vector3d lumen_cent = Vector3d(0.,0.,0.); // this point can be any arbitrary point, since signed volumes are added for a closed surface
  for(auto e : ListEdge)
  {
    lumen_volume += signed_terahedron_volume(e.tail->pos_a, e.head->pos_a, e.cell->cent_a, lumen_cent);
  }
  #endif
  return lumen_volume;
}

Vector3d Tissue::get_cell_apical_centroid(Cell *c)
{
  c->cent_a = Vector3d(0,0,0);
  for(std::list<Vertex*>::const_iterator it = c->Vlist.begin(); it != c->Vlist.end(); it++)
  {
  auto &v = *it;
  c->cent_a += v->pos_a;
  }
  c->cent_a /= c->Vlist.size();
  return c->cent_a;
}

Vector3d Tissue::get_cell_basal_centroid(Cell *c)
{
  c->cent_b = Vector3d(0,0,0);
  for(std::list<Vertex*>::const_iterator it = c->Vlist.begin(); it != c->Vlist.end(); it++)
  {
  auto &v = *it;
  c->cent_b += v->pos_b;
  }
  c->cent_b /= c->Vlist.size();
  return c->cent_b;
}

double Tissue::get_cell_apical_area(Cell *c)
{
  //Vector3d area_vec(0,0,0);
  //double a_scalar=0;
  c->area_a = 0.;
  for(std::list<Edge*>::const_iterator it = c->Elist.begin(); it != c->Elist.end(); it++)
  {
  auto &e = *it;
  Vector3d u1 = e->head->pos_a - c->cent_a;
  Vector3d u2 = e->tail->pos_a - c->cent_a;
  Vector3d avec = u1.cross(u2);
  c->area_a += avec.norm();
  //std::cout << (u1.cross(u2)).dot(c->cent_a) << std::endl;
  //a_scalar += (u1.cross(u2)).norm()/2.;
  }
  c->area_a /= 2.;
  //std::cout << "****   " << c->area_a/a_scalar << std::endl;
  return c->area_a;
}

double Tissue::get_cell_basal_area(Cell *c)
{
  //Vector3d area_vec(0,0,0);
  //double a_scalar=0;
  c->area_b = 0.;
  for(std::list<Edge*>::const_iterator it = c->Elist.begin(); it != c->Elist.end(); it++)
  {
  auto &e = *it;
  Vector3d u1 = e->head->pos_b - c->cent_b;
  Vector3d u2 = e->tail->pos_b - c->cent_b;
  Vector3d avec = u1.cross(u2);
  c->area_b += avec.norm();
  //std::cout << (u1.cross(u2)).dot(c->cent_a) << std::endl;
  //a_scalar += (u1.cross(u2)).norm()/2.;
  }
  c->area_b /= 2.;
  //std::cout << "****   " << c->area_a/a_scalar << std::endl;
  return c->area_b;
}

double Tissue::get_cell_apical_perimeter(Cell *c)
{
  c->perim_a = 0.;
  for(std::list<Edge*>::const_iterator it = c->Elist.begin(); it != c->Elist.end(); it++)
  {
    auto &e = *it;
    c->perim_a += e->l_a;
  }
  return c->perim_a;
}

double Tissue::get_cell_basal_perimeter(Cell *c)
{
  c->perim_b = 0.;
  for(std::list<Edge*>::const_iterator it = c->Elist.begin(); it != c->Elist.end(); it++)
  {
    auto &e = *it;
    c->perim_b += e->l_b;
  }
  return c->perim_b;
}

std::pair<std::vector<Vector3d>, std::vector<Vector3d>> Tissue::get_vertex_pos_changes_of_cell(Cell *c)
{
  std::vector<Vector3d> pos_t1, pos_t2; pos_t1.reserve(c->Vlist.size()); pos_t2.reserve(c->Vlist.size());

  Vector3d t1_Bc(0,0,0), t2_Bc(0,0,0);

  for(std::list<Vertex*>::const_iterator it = c->Vlist.begin(); it != c->Vlist.end(); it++)
  {
    auto &v = *it;
    t1_Bc += v->pos_a_prev;
    t2_Bc += v->pos_a;
    pos_t2.emplace_back(v->pos_a);
    pos_t1.emplace_back(v->pos_a_prev);
  }

  c->cent_a_prev = t1_Bc/pos_t1.size();
  c->cent_a = t2_Bc/pos_t2.size();
  return std::make_pair(pos_t1, pos_t2);
}

Matrix3d Tissue::get_apical_elongation_tensor(Cell *c)
{
  Matrix3d coarsegrained_Q2d;
  coarsegrained_Q2d(0,0)=0.; coarsegrained_Q2d(0,1)=0.; coarsegrained_Q2d(0,2)=0.;
  coarsegrained_Q2d(1,0)=0.; coarsegrained_Q2d(1,1)=0.; coarsegrained_Q2d(1,2)=0.;
  coarsegrained_Q2d(2,0)=0.; coarsegrained_Q2d(2,1)=0.; coarsegrained_Q2d(2,2)=0.;

  double sum_area=0.;
  Vector3d avg_normals(0,0,0);
  for(std::list<Edge*>::const_iterator it = c->Elist.begin(); it != c->Elist.end(); it++)
  {
    auto &e = *it;
    auto [tri_Q2d, tri_n, tri_area] = calculate_triangle_elongation_tensor_xy_plane(c->cent_a, e->head->pos_a, e->tail->pos_a);
    coarsegrained_Q2d += tri_Q2d * tri_area;
    avg_normals += tri_n * tri_area;
    sum_area += tri_area;
  }
  avg_normals /= sum_area;
  coarsegrained_Q2d /= sum_area;

  double rot_ang_X, rot_ang_Y;
  std::tie(rot_ang_X, rot_ang_Y) = normal_ThetaX_ThetaY(avg_normals);

  c->Q_3d = transform_matrices_from_xy_plane_to_triangle_plane(coarsegrained_Q2d, rot_ang_X, rot_ang_Y);

  return c->Q_3d;
}

double Tissue::get_vertex_lateral_distance(Vertex *v)
{
  v->lateral_distance = (v->pos_b - v->pos_a).norm();
  return v->lateral_distance;
}

Vector3d Tissue::get_edge_lateral_centroid(Edge *e)
{
  e->cent_l = (e->tail->pos_a + e->tail->pos_b + e->head->pos_a + e->head->pos_b)/4.;
  return e->cent_l;
}

double Tissue::get_lateral_area(Edge *e)
{
  //e->area_lateral = 0.;
  double a1 = ((e->tail->pos_a - e->cent_l).cross(e->head->pos_a - e->cent_l)).norm();
  double a2 = ((e->tail->pos_b - e->cent_l).cross(e->head->pos_b - e->cent_l)).norm();
  double a3 = ((e->tail->pos_a - e->cent_l).cross(e->tail->pos_b - e->cent_l)).norm();
  double a4 = ((e->head->pos_a - e->cent_l).cross(e->head->pos_b - e->cent_l)).norm();
  e->area_lateral = (a1+a2+a3+a4)/2.;
  return e->area_lateral;
}

double Tissue::get_apical_bond_length(Edge *e)
{
  e->l_a_vec = e->head->pos_a - e->tail->pos_a;
  e->l_a = e->l_a_vec.norm();
  return e->l_a;
}

double Tissue::get_basal_bond_length(Edge *e)
{
  e->l_b_vec = e->head->pos_b - e->tail->pos_b;
  e->l_b = e->l_b_vec.norm();
  return e->l_b;
}

void Tissue::update_all_apical_bond_lengths()
{
  #if PRALLEL
  auto func = [this](auto&& it){this->get_apical_bond_length(&it);};
  // __gnu_parallel::for_each(
  //   ListEdge.begin(),
  //   ListEdge.end(),
  //   func);
  std::for_each(
    std::execution::par_unseq,
    ListEdge.begin(),
    ListEdge.end(),
    func);
  #else
  for(std::vector<Edge>::iterator it = ListEdge.begin(); it != ListEdge.end(); it++)
  {
    get_apical_bond_length(&(*it));
  }
  #endif
}

void Tissue::update_all_basal_bond_lengths()
{
  #if PRALLEL
  auto func = [this](auto&& it){this->get_basal_bond_length(&it);};
  std::for_each(
    std::execution::par_unseq,
    ListEdge.begin(),
    ListEdge.end(),
    func);
  #else
  for(std::vector<Edge>::iterator it = ListEdge.begin(); it != ListEdge.end(); it++)
  {
    get_basal_bond_length(&(*it));
  }
  #endif
}


void Tissue::rescale_spheroid_apical_size()
{
  double R_sp = std::sqrt ( ListCell.size()* CellMechProp.A0c/(4.*M_PI) );
  Vector3d t_cent(0,0,0);
  #if PRALLEL
  // auto func = [t_cent](auto&& it){t_cent[0] += ((&it)->pos_a)[0]; t_cent[1] += ((&it)->pos_a)[1]; t_cent[2] += ((&it)->pos_a)[2];};
  // std::for_each(
  //   std::execution::par_unseq,
  //   ListVertex.begin(),
  //   ListVertex.end(),
  //   func);
  for(std::vector<Vertex>::iterator it = ListVertex.begin(); it != ListVertex.end(); it++) {t_cent += it->pos_a;}
  t_cent /= ListVertex.size();

  auto func2 = [t_cent,R_sp](auto&& it){(&it)->pos_a -= t_cent; (&it)->pos_a = (&it)->pos_a *(R_sp/((&it)->pos_a.norm()) );};
  std::for_each(
    std::execution::par_unseq,
    ListVertex.begin(),
    ListVertex.end(),
    func2);
  #else
  for(std::vector<Vertex>::iterator it = ListVertex.begin(); it != ListVertex.end(); it++) {t_cent += it->pos_a;}
  t_cent /= ListVertex.size();

  for(std::vector<Vertex>::iterator it = ListVertex.begin(); it != ListVertex.end(); it++)
  {
  it->pos_a -= t_cent;
  it->pos_a = it->pos_a *(R_sp/(it->pos_a.norm()) );
  }
  #endif
}

void Tissue::build_2D_spherical_tissue_from_file(const char* fileName)
{
  std::ifstream myfile;
  myfile.open(fileName);

  int Ncells;
  double v1x, v1y, v1z;
  double vx_this, vy_this, vz_this;
  Vertex *Ver1=nullptr, *thisVer=nullptr, *VerOld=nullptr;

  if (myfile.is_open())
  {
    myfile >> Ncells; // read the 1st for the number of cells

    for(int cid=0; cid<Ncells; cid++)
    {
      std::list<Vertex*> vert_list;
      std::list<Edge*> edge_list;

      int nv; // number of cell vertices
      myfile >> nv;

      myfile >> v1x >> v1y >> v1z;//getting the first apical vertex coordinate for this cell
      Vector3d pos_a_1(v1x, v1y, v1z);
      Ver1 = similar_or_new_vertex_2D(pos_a_1);
      vert_list.push_back(Ver1);// adding pointer to list for c->Vlist

      VerOld = Ver1;
      for(int itv=0;itv<nv-1;itv++) //iterating over the remaining nv-1 cell vertices
      {
        myfile >> vx_this >> vy_this >> vz_this;
        Vector3d pos_a_this(vx_this, vy_this, vz_this);
        thisVer = similar_or_new_vertex_2D(pos_a_this);
        vert_list.push_back(thisVer);// adding pointer to list for c->Vlist

        Edge *e = similar_or_new_edge(VerOld, thisVer);
        edge_list.push_back(e);
        VerOld=thisVer;
      }

      Edge *e = similar_or_new_edge(thisVer, Ver1);
      edge_list.push_back(e); // finished with making cell edges
      Cell &c = ListCell.emplace_back(Cell(vert_list, edge_list, cid));// now the cell can be made

      // now we set the next edge of each cell edge
      std::list<Edge*>::iterator itr = edge_list.begin();
      e->next = *itr;
      //(*itr)->next = e;
      while(itr != --edge_list.end())
      {
        Edge *ee = *itr;
        itr++;
        Edge *en = *itr;
        ee->next = en;
        //en->next = ee;
      }

      // get_cell_apical_centroid(&c);
      // double ca = get_cell_apical_area(&c);
      // std::cout << "area:  " << ca << std::endl;
      // //std::cout << cent[0] << " " << cent[1] << " " << cent[2] << std::endl;

      for(std::list<Edge*>::iterator ite=edge_list.begin();ite!=edge_list.end();ite++) (*ite)->cell = &c;
    }
    //  at this point only half of conj bonds are set. let set the other half
    for(std::vector<Edge>::iterator ite=ListEdge.begin();ite!=ListEdge.end();ite++)
    {
      if(ite->conj != nullptr) ite->conj->conj = &(*ite);
    }
  }

  myfile.close();
}


void Tissue::build_2D_dodecahedron_tissue()
{
  std::vector<Vector3d> dodecaheron_verts;
  std::vector<std::vector<int>> dodecaheron_faces;
  std::tie(dodecaheron_verts, dodecaheron_faces) = dodecahedron_verts_faces();

  //std::assertm(dodecaheron_faces.size()==12, "Reading from Dodecahedron did not give 2 faces!");
  Vertex *Ver1=nullptr, *thisVer=nullptr, *VerOld=nullptr;
  for(int cid=0; cid<12; cid++)
  {
    std::list<Vertex*> vert_list;
    std::list<Edge*> edge_list;

    int vid = dodecaheron_faces[cid][0];
    Vector3d pos_a_1 = dodecaheron_verts[vid];
    Ver1 = similar_or_new_vertex_2D(pos_a_1);
    vert_list.push_back(Ver1);// adding pointer to list for c->Vlist

    VerOld = Ver1;
    int nv = 5;
    for(int itv=0;itv<nv-1;itv++) //iterating over the remaining nv-1 cell vertices
    {
      vid = dodecaheron_faces[cid][itv+1];
      Vector3d pos_a_this = dodecaheron_verts[vid];
      thisVer = similar_or_new_vertex_2D(pos_a_this);
      vert_list.push_back(thisVer);// adding pointer to list for c->Vlist

      Edge *e = similar_or_new_edge(VerOld, thisVer);
      edge_list.push_back(e);
      VerOld=thisVer;
    }

    Edge *e = similar_or_new_edge(thisVer, Ver1);
    edge_list.push_back(e); // finished with making cell edges
    Cell &c = ListCell.emplace_back(Cell(vert_list, edge_list, cid));// now the cell can be made

    // now we set the next edge of each cell edge
    std::list<Edge*>::iterator itr = edge_list.begin();
    e->next = *itr;
    while(itr != --edge_list.end())
    {
      Edge *ee = *itr;
      itr++;
      Edge *en = *itr;
      ee->next = en;
    }

    for(std::list<Edge*>::iterator ite=edge_list.begin();ite!=edge_list.end();ite++) (*ite)->cell = &c;
  }

  //  at this point only half of conj bonds are set. let set the other half
  for(std::vector<Edge>::iterator ite=ListEdge.begin();ite!=ListEdge.end();ite++)
  {
    if(ite->conj != nullptr) ite->conj->conj = &(*ite);
  }
}


void Tissue::build_3D_spherical_tissue_from_file(const char* fileName, double R_sph, double Rb_Ra_ratio)
{
  std::ifstream myfile;
  myfile.open(fileName);

  int Ncells;
  double v1x, v1y, v1z;
  double vx_this, vy_this, vz_this;
  Vertex *Ver1=nullptr, *thisVer=nullptr, *VerOld=nullptr;

  double r0 = 1.;

  if (myfile.is_open())
  {
    myfile >> Ncells; // read the 1st for the number of cells

    for(int cid=0; cid<Ncells; cid++)
    {
      std::list<Vertex*> vert_list;
      std::list<Edge*> edge_list;

      int nv; // number of cell vertices
      myfile >> nv;

      myfile >> v1x >> v1y >> v1z;//getting the first apical vertex coordinate for this cell
      r0 = std::sqrt( v1x*v1x + v1y*v1y + v1z*v1z );
      Vector3d pos_a_1(R_sph*v1x/r0, R_sph*v1y/r0, R_sph*v1z/r0);
      Vector3d pos_b_1 = pos_a_1 * Rb_Ra_ratio;
      Ver1 = similar_or_new_vertex_3D(pos_a_1, pos_b_1);
      vert_list.push_back(Ver1);// adding pointer to list for c->Vlist

      VerOld = Ver1;
      for(int itv=0;itv<nv-1;itv++) //iterating over the remaining nv-1 cell vertices
      {
        myfile >> vx_this >> vy_this >> vz_this;
        r0 = std::sqrt( vx_this*vx_this + vy_this*vy_this + vz_this*vz_this );
        Vector3d pos_a_this(R_sph*vx_this/r0, R_sph*vy_this/r0, R_sph*vz_this/r0);
        Vector3d pos_b_this = pos_a_this * Rb_Ra_ratio;
        thisVer = similar_or_new_vertex_3D(pos_a_this, pos_b_this);
        vert_list.push_back(thisVer);// adding pointer to list for c->Vlist

        Edge *e = similar_or_new_edge(VerOld, thisVer);
        edge_list.push_back(e);
        VerOld=thisVer;
      }

      Edge *e = similar_or_new_edge(thisVer, Ver1);
      edge_list.push_back(e); // finished with making cell edges
      Cell &c = ListCell.emplace_back(Cell(vert_list, edge_list, cid));// now the cell can be made

      // now we set the next edge of each cell edge
      std::list<Edge*>::iterator itr = edge_list.begin();
      e->next = *itr;
      //(*itr)->next = e;
      while(itr != --edge_list.end())
      {
        Edge *ee = *itr;
        itr++;
        Edge *en = *itr;
        ee->next = en;
        //en->next = ee;
      }

      // get_cell_apical_centroid(&c);
      // double ca = get_cell_apical_area(&c);
      // std::cout << "area:  " << ca << std::endl;
      // //std::cout << cent[0] << " " << cent[1] << " " << cent[2] << std::endl;

      for(std::list<Edge*>::iterator ite=edge_list.begin();ite!=edge_list.end();ite++) (*ite)->cell = &c;
    }
    //  at this point only half of conj bonds are set. let set the other half
    for(std::vector<Edge>::iterator ite=ListEdge.begin();ite!=ListEdge.end();ite++)
    {
      if(ite->conj != nullptr) ite->conj->conj = &(*ite);
    }
  }

  myfile.close();
}

void Tissue::build_3D_spherical_tissue_from_json_file(std::string json_file)
{
  using json = nlohmann::json;

  std::ifstream f(json_file);
  json json_data = json::parse(f);

  int Nv=json_data["Num_vertices"];
  int Ne=json_data["Num_edges"];
  int Nc=json_data["Num_cells"];

  auto vertex_coord_a = json_data["vertex_pos_a"];
  auto vertex_prev_coord_a = json_data["vertex_prev_pos_a"];
  auto vertex_coord_b = json_data["vertex_pos_b"];
  auto cell_polarity = json_data["cells_polarity"];
  auto cell_Factive_mag = json_data["cell_Factive_mag"];
  auto cell_vertex_ids = json_data["cell_vertex_ids"];

  ListVertex.reserve(Nv);
  ListCell.reserve(Nc);
  ListEdge.reserve(Ne);

  for(int i=0; i<Nv; i++)
  {
    Vector3d posa = Vector3d(vertex_coord_a[i][0], vertex_coord_a[i][1], vertex_coord_a[i][2]);
    Vector3d posb = Vector3d(vertex_coord_b[i][0], vertex_coord_b[i][1], vertex_coord_b[i][2]);
    Vector3d posa_prev = Vector3d(vertex_prev_coord_a[i][0], vertex_prev_coord_a[i][1], vertex_prev_coord_a[i][2]);
    Vertex *vv = similar_or_new_vertex_3D(posa, posb);
    vv->pos_a_prev = posa_prev;
  }

  for(int cid=0; cid<Nc; cid++)
  {
    auto cell_vids = cell_vertex_ids[cid];
    int nv = cell_vids.size();

    std::list<Vertex*> vert_list;
    std::list<Edge*> edge_list;

    Vertex *Ver1 = &ListVertex[cell_vids[0]];
    vert_list.push_back(Ver1);// adding pointer to list for c->Vlist

    Vertex *VerOld = Ver1;
    Vertex *thisVer = nullptr;
    for(int itv=0;itv<nv-1;itv++) //iterating over the remaining nv-1 cell vertices
    {
      thisVer = &ListVertex[cell_vids[itv+1]];
      vert_list.push_back(thisVer);// adding pointer to list for c->Vlist

      Edge *e = similar_or_new_edge(VerOld, thisVer);
      edge_list.push_back(e);
      VerOld=thisVer;
    }

    Edge *e = similar_or_new_edge(thisVer, Ver1);
    edge_list.push_back(e); // finished with making cell edges
    Cell &c = ListCell.emplace_back(Cell(vert_list, edge_list, cid));// now the cell can be made
    c.polarity = Vector3d(cell_polarity[cid][0], cell_polarity[cid][1], cell_polarity[cid][2]);
    c.mag_F_a = cell_Factive_mag[cid];

    // now we set the next edge of each cell edge
    std::list<Edge*>::iterator itr = edge_list.begin();
    e->next = *itr;
    //(*itr)->next = e;
    while(itr != --edge_list.end())
    {
      Edge *ee = *itr;
      itr++;
      Edge *en = *itr;
      ee->next = en;
      //en->next = ee;
    }

    for(std::list<Edge*>::iterator ite=edge_list.begin();ite!=edge_list.end();ite++) (*ite)->cell = &c;
  }

  //  at this point only half of conj bonds are set. let set the other half
  for(std::vector<Edge>::iterator ite=ListEdge.begin();ite!=ListEdge.end();ite++)
  {
    if(ite->conj != nullptr) ite->conj->conj = &(*ite);
  }

}

void Tissue::build_3D_spherical_tissue_from_output_files(const char* vertices_file, const char* cells_file, const char* cell_data_file)
{
  std::vector<int> cells_nv; cells_nv.reserve(1000);

  std::ifstream cdata_file;
  cdata_file.open(cell_data_file);
  int nv_cell; // number of cell vertices
  double px, py, pz, F_mag; // cell polarity and traction force magnitude
  int cid=0;
  if(cdata_file.is_open())
  while(cdata_file.good())
  {
    cdata_file >> nv_cell >> px >> py >> pz >> F_mag;
    if(cdata_file.eof()) break;
    //std::cout << "line: " << nv_cell << " " << px << " " << py << " " << pz << " " << F_mag << std::endl;
    cells_nv.emplace_back(nv_cell);
    std::list<Vertex*> vert_list;
    std::list<Edge*> edge_list;
    Cell *c = &ListCell.emplace_back(Cell(vert_list, edge_list, cid));// vert_list, edge_list will be filled later
    c->polarity = Vector3d(px, py, pz); c->mag_F_a = F_mag; //c->cent_a_prev = Vector3d(cent_prev_x, cent_prev_y, cent_prev_z);
    cid++;
    //std::cout << cid << std::endl;
  }
  cdata_file.close();
  for(size_t ss=0; ss<cells_nv.size(); ss++) {std::cout << cells_nv[ss] << " ";}

  std::cout << std::endl;
  //----------------------
  double v1x, v1y, v1z, v1bx, v1by, v1bz;
  //double vx_this, vy_this, vz_this;

  std::ifstream vert_file;
  vert_file.open(vertices_file);
  double pos_a_prev_x, pos_a_prev_y, pos_a_prev_z;
  if(vert_file.is_open())
  while (vert_file.good())
  {
    vert_file >> v1x >> v1y >> v1z >> v1bx >> v1by >> v1bz >> pos_a_prev_x >> pos_a_prev_y >> pos_a_prev_z;
    //if(vert_file.eof()) break;
    Vector3d pos_a(v1x, v1y, v1z); Vector3d pos_b(v1bx, v1by, v1bz);
    Vertex *vv = similar_or_new_vertex_3D(pos_a, pos_b);
    vv->pos_a_prev = Vector3d(pos_a_prev_x, pos_a_prev_y, pos_a_prev_z);// the code will update cell->cent_a_prev accordingly
  }
  vert_file.close();
  //-------------------------------
  int Ncells = cells_nv.size();
  std::cout << "Ncells:  " << Ncells << std::endl;
  std::ifstream c_file; // assign vertices and edges to each cell
  c_file.open(cells_file);
  if(c_file.is_open())
  while (c_file.good())
  {
    for(cid=0; cid<Ncells; cid++)
    {
      Cell *c = &ListCell[cid];
      c->Vlist.clear(); c->Elist.clear(); //SOMETHIG GOES WRONG IN READING FILE
      int v_idx; // vertex index
      for(int v_count=0; v_count<cells_nv[cid]; v_count++)
      {
        c_file >> v_idx;
        //if(c_file.eof()) break;
        if(cid==0) std::cout << "writing v_idx: " <<v_idx << std::endl;
        c->Vlist.emplace_back(&ListVertex[v_idx]);// filling each cells Vlist
      }
      if(cid==0) for(auto v : c->Vlist) std:: cout << "double chec: " <<  v->id << std::endl;
      //if(c_file.eof()) break;
      Vertex *thisVer=nullptr, *VerOld=nullptr;
      std::list<Vertex*>::iterator itv = c->Vlist.begin();
      Vertex *Ver1 = *itv;
      VerOld = Ver1;
      for(int v_count=0; v_count<cells_nv[cid]-1; v_count++) //iterating over the remaining nv-1 cell vertices
      {
        thisVer = *(itv++);
        Edge *e = similar_or_new_edge(VerOld, thisVer);
        c->Elist.push_back(e);
        VerOld=thisVer;
      }

      Edge *e = similar_or_new_edge(thisVer, Ver1);
      c->Elist.push_back(e); // finished with making cell edges

      // now we set the next edge of each cell edge
      std::list<Edge*>::iterator itr = c->Elist.begin();
      e->next = *itr; // e is still the last edge of the cell
      //(*itr)->next = e;
      while(itr != --c->Elist.end())
      {
        Edge *ee = *itr;
        itr++;
        Edge *en = *itr;
        ee->next = en;
        //en->next = ee;
      }

      for(std::list<Edge*>::iterator ite=c->Elist.begin();ite!=c->Elist.end();ite++) (*ite)->cell = c;

    }
    //  at this point only half of conj bonds are set. let set the other half
    for(std::vector<Edge>::iterator ite=ListEdge.begin();ite!=ListEdge.end();ite++)
    {
      if(ite->conj != nullptr) ite->conj->conj = &(*ite);
    }
  }
  c_file.close();
}

void Tissue::update_apical_geometric_quantities()
{
  update_all_apical_bond_lengths();

  #if PRALLEL
  auto func = [this](auto&& it){this->get_cell_apical_centroid(&it); this->get_cell_apical_area(&it); this->get_cell_apical_perimeter(&it);};
  // __gnu_parallel::for_each(
  //   ListCell.begin(),
  //   ListCell.end(),
  //   func);
  std::for_each(
    std::execution::par_unseq,
    ListCell.begin(),
    ListCell.end(),
    func);
  #else
  for(std::vector<Cell>::iterator it=ListCell.begin();it!=ListCell.end();it++)
  {
    get_cell_apical_centroid(&(*it));
    get_cell_apical_area(&(*it));
    get_cell_apical_perimeter(&(*it));
  }
  #endif
}

void Tissue::update_basal_geometric_quantities()
{
  update_all_basal_bond_lengths();

  #if PRALLEL
  auto func = [this](auto&& it){this->get_cell_basal_centroid(&it); this->get_cell_basal_area(&it); this->get_cell_basal_perimeter(&it);};
  std::for_each(
    std::execution::par_unseq,
    ListCell.begin(),
    ListCell.end(),
    func);
  #else
  for(std::vector<Cell>::iterator it=ListCell.begin();it!=ListCell.end();it++)
  {
    get_cell_basal_centroid(&(*it));
    get_cell_basal_area(&(*it));
    get_cell_basal_perimeter(&(*it));
  }
  #endif
}

void Tissue::update_lateral_geometric_quantities()
{
  #if PRALLEL
  auto func = [this](auto&& it){this->get_vertex_lateral_distance(&it);};
  std::for_each(
    std::execution::par_unseq,
    ListVertex.begin(),
    ListVertex.end(),
    func);
  #else
  for(std::vector<Vertex>::iterator it=ListVertex.begin(); it!=ListVertex.end(); it++) get_vertex_lateral_distance(&(*it));
  #endif

  #if PRALLEL
  auto func2 = [this](auto&& it){this->get_edge_lateral_centroid(&it); this->get_lateral_area(&it);};
  std::for_each(
    std::execution::par_unseq,
    ListEdge.begin(),
    ListEdge.end(),
    func2);
  #else
  for(std::vector<Edge>::iterator it=ListEdge.begin(); it!=ListEdge.end(); it++)
  {
    get_edge_lateral_centroid(&(*it));
    get_lateral_area(&(*it));
  }
  #endif
}

void Tissue::update_3D_geometric_quantities()
{
  update_apical_geometric_quantities();
  update_basal_geometric_quantities();
  update_lateral_geometric_quantities();
  for(std::vector<Cell>::iterator it=ListCell.begin();it!=ListCell.end();it++) {get_cell_centroid(&(*it)); get_cell_volume(&(*it));}

  if(SphMechProp.Kl>0 || SphMechProp.P0 !=0) get_lumen_volume();
}

std::vector<int> Tissue::find_short_bonds(double bond_T1_cutoff, bool random_shufflin)
{
  std::vector<int> short_indices;
  for(size_t i=0; i<ListEdge.size(); i++) if(ListEdge[i].id < ListEdge[i].conj->id && ListEdge[i].l_a<bond_T1_cutoff) short_indices.emplace_back(i);

  if(random_shufflin)
  {
    auto rd = std::random_device {};
    auto rng = std::default_random_engine { rd() };
    std::shuffle(std::begin(short_indices), std::end(short_indices), rng);
  }

  return short_indices;
}

bool Tissue::flip_edge(Edge *e, double opening_length)
{
  Vertex *vt = e->tail;
  Vertex *vh = e->head;

  Cell *c1 = e->cell; // right hand side of and owning the helf-edge
  Cell *c2 = e->conj->cell; // left hand side of the edge, andowning it's conj
  Cell *c3 = e->conj->next->conj->cell; // connected to tail
  Cell *c4 = e->next->conj->cell; // connected to head

  //std::cout << "c3 old nv: " << c3->Vlist.size() << "," << c3->Elist.size() << "    c4 old nv: " << c4->Vlist.size() << "," << c4->Elist.size() << std::endl;

  if(c1->Vlist.size() < 4 || c2->Vlist.size() < 4) return 0; //triangle can not lose another side

  Vector3d c1_cent_a = c1->cent_a, c2_cent_a = c2->cent_a;
  Vector3d new_e_hat_a = c2_cent_a - c1_cent_a; new_e_hat_a/=new_e_hat_a.norm();
  Vector3d e_mid_a = (e->tail->pos_a + e->head->pos_a)/2.;

  #if THREE_DIMENSIONS
  // basal setCoordinates
  Vector3d c1_cent_b = c1->cent_b, c2_cent_b = c2->cent_b;
  Vector3d new_e_hat_b = c2_cent_b - c1_cent_b; new_e_hat_b/=new_e_hat_b.norm();
  Vector3d e_mid_b = (e->tail->pos_b + e->head->pos_b)/2.;
  // Vector3d c3_cent_a = c3->cent_a, c4_cent_a = c4->cent_a;
  #endif

  // std::cout << c1->id << " " << c2->id << " " << c3->id << " " << c4->id << std::endl;
  //language of indices:  n=next,   s=conj
  Edge *e_s = e->conj;
  Edge *e_sn = e->conj->next;
  Edge *e_sns = e->conj->next->conj;
  Edge *e_snsn = e->conj->next->conj->next;
  Edge *e_snsns = e->conj->next->conj->next->conj;
  //...
  Edge *e_n = e->next;
  Edge *e_ns = e->next->conj;
  Edge *e_nsn = e->next->conj->next;
  Edge *e_nsns = e->next->conj->next->conj;
  //... adding and removing vertices from involved cells ...
  std::list<Vertex*>::iterator it=c2->Vlist.begin(); // c2 will lose the e->tail
  while(it!=c2->Vlist.end()) { if((*it)==vt){c2->Vlist.erase(it++); break;} else {++it;} }
  it=c1->Vlist.begin(); // c1 will lose the e->head
  while(it!=c1->Vlist.end()) { if((*it)==vh){c1->Vlist.erase(it++); break;} else {++it;} }
  it=c4->Vlist.begin(); // c4 will gain the e->tail
  while(it!=c4->Vlist.end()) { if((*it)==vh){c4->Vlist.insert(it,vt); break;} else {++it;} }
  it=c3->Vlist.begin(); // c3 will gain the e->head
  while(it!=c3->Vlist.end()) { if((*it)==vt){c3->Vlist.insert(it,vh); break;} else {++it;} }

  // ... adding and removing the flipping bond from the involved cells ...
  std::list<Edge*>::iterator ite=c1->Elist.begin(); // c1 will lose the e
  while(ite!=c1->Elist.end()) { if((*ite)==e){c1->Elist.erase(ite++); break;} else {++ite;} }
  ite=c2->Elist.begin(); // c2 will lose the e->conj
  while(ite!=c2->Elist.end()) { if((*ite)==e_s){c2->Elist.erase(ite++); break;} else {++ite;} }
  ite=c4->Elist.begin(); // c4 will gain the e
  while(ite!=c4->Elist.end()) { if((*ite)==e_nsn){c4->Elist.insert(ite,e); break;} else {++ite;} }
  ite=c3->Elist.begin(); // c3 will gain the e->conj
  while(ite!=c3->Elist.end()) { if((*ite)==e_snsn){c3->Elist.insert(ite,e_s); break;} else {++ite;} }

  //... now, fixing the "next" pointer of each affected edge ..
  e_sns->next = e_s; e_s->next = e_snsn; // in c3
  e_snsns->next = e_n;  // in c1
  e_nsns->next = e_sn; // in c2
  e_ns->next = e; e->next = e_nsn; // in c4
  e_n->tail = vt; e_ns->head = vt;
  e_sn->tail = vh; e_sns->head = vh;
  e->cell = c4; e_s->cell = c3; // e->cell changes to c4, and e->conj->cell to c3
  // === now topology change is done ====
  e->head->pos_a = e_mid_a + new_e_hat_a*opening_length/2.;
  e->tail->pos_a = e_mid_a - new_e_hat_a*opening_length/2.;

  #if THREE_DIMENSIONS
  e->head->pos_b = e_mid_b + new_e_hat_b*opening_length/2.;
  e->tail->pos_b = e_mid_b - new_e_hat_b*opening_length/2.;
  update_3D_geometric_quantities();
  #else
  update_apical_geometric_quantities(); //TODO: do for only 10 involved bonds and 4 involved cells. Do the same for basal
  #endif

  //std::cout << "c3 new nv: " << c3->Vlist.size() << "," << c3->Elist.size() << "    c4 new nv: " << c4->Vlist.size() << "," << c4->Elist.size() << std::endl;

  return 1;
}

bool Tissue::flip_short_bonds(double bond_T1_cutoff, double opening_length, bool random_shuffling)
{
  std::vector<int> short_bonds =  find_short_bonds(bond_T1_cutoff, random_shuffling);

  if (short_bonds.size()<1) return false;
  else
  {
    for(auto ei : short_bonds) flip_edge(&ListEdge[ei], opening_length);
  }
  return true;
}


bool Tissue::divide_cell_random_axis(Cell *M)
{
  // M is the mother cell, and let's name the new-born X
  int nv = M->Vlist.size();
  //std::cout << "M inialz: " << M->Vlist.size() << " " << M->Vlist.size() << std::endl;
  // Seed with a real random value, if available
  std::random_device r;
  // Choose a random mean between 0 and nv
  std::default_random_engine egn(r());
  std::uniform_int_distribution<int> uniform_dist(0, nv-1);
  int e1_i = uniform_dist(egn);
  int forw_seq = nv/2; // to find the edge that cuts cell almost in half
  int e2_i = (e1_i + forw_seq)%nv;
  //std::cout << "nv: " << nv << "    Randomly-chosen b: " << e1_i << "  " << e2_i << '\n';

  Edge *e1=nullptr, *e2=nullptr; //pointers to cutting edges
  int e_ct=0;
  for(std::list<Edge*>::iterator it = M->Elist.begin(); it != M->Elist.end(); it++)
  {
  if(e_ct==e1_i) e1 = *it;
  else if(e_ct==e2_i) e2 = *it;

  e_ct++;
  }

  Cell *c1 = e1->conj->cell;
  Cell *c2 = e2->conj->cell;
  //std::cout <<"initials C1: " << c1->Vlist.size() << " " << c1->Elist.size() << "   intilaz C2: " << c2->Vlist.size() << " " << c2->Elist.size() << std::endl;
  //if(e1==nullptr || e2==nullptr) std::cout << " Warn!" << std::endl;
  Vector3d e1_mid_a = (e1->tail->pos_a + e1->head->pos_a)/2.;
  //Vector3d e1_mid_b = (e1->tail->pos_b + e1->head->pos_b)/2.;
  Vector3d e2_mid_a = (e2->tail->pos_a + e2->head->pos_a)/2.;
  //Vector3d e2_mid_b = (e2->tail->pos_b + e2->head->pos_b)/2.;

  Vertex *v_e1m = similar_or_new_vertex_2D(e1_mid_a); // vertex in the middle of e1
  Vertex *v_e2m = similar_or_new_vertex_2D(e2_mid_a); // vertex in the middle of e2

  Edge *e1_sn = e1->conj->next;
  //Edge *e1_sns = e1->conj->next->conj;
  //Edge *e1_snsn = e1->conj->next->conj->next;
  Edge *e1_snsns = e1->conj->next->conj->next->conj;

  Edge *e2_n = e2->next;
  //Edge *e2_ns = e2->next->conj;
  //Edge *e2_nsn = e2->next->conj->next;
  Edge *e2_nsns = e2->next->conj->next->conj;

  Edge *eM = similar_or_new_edge(v_e2m, v_e1m); eM->cell = M;
  Edge *eM_s = similar_or_new_edge(v_e1m, v_e2m);
  eM->conj = eM_s; eM_s->conj = eM;

  Edge *e1X = similar_or_new_edge(e1->tail, v_e1m);
  Edge *e1X_s = similar_or_new_edge(v_e1m, e1->tail); e1X_s->cell = c1;
  e1X->conj = e1X_s; e1X_s->conj = e1X;

  Edge *e2X = similar_or_new_edge(v_e2m, e2->head);
  Edge *e2X_s = similar_or_new_edge(e2->head, v_e2m); e2X_s->cell = c2;
  e2X->conj = e2X_s; e2X_s->conj = e2X;

  std::list<Edge*> M_Elist={eM, e1}, X_Elist={e1X, eM_s, e2X};//edge pointers of mother and the new cell
  std::list<Vertex*> M_Vlist={v_e2m, v_e1m}, X_Vlist={v_e1m, v_e2m, e2X->head}; //vertex pointers of mother and the new cell

  Edge *e_it=e1; //filling in the rest
  for(int i=0; i<forw_seq; i++)
  {
    M_Vlist.push_back(e_it->head);
    e_it = e_it->next;
    M_Elist.push_back(e_it);
  }

  Edge *e_it2 = e2;
  for(int i=0; i<(nv-forw_seq-1); i++)
  {
    e_it2 = e_it2->next;
    X_Elist.push_back(e_it2);
    X_Vlist.push_back(e_it2->head);
  }

  Cell *X = &(ListCell.emplace_back(Cell(X_Vlist, X_Elist, ListCell.size())));// now the cell X can be born
  M->Elist = M_Elist; M->Vlist = M_Vlist; // update mother cell edges and vertices

  for(std::list<Edge*>::iterator it=X->Elist.begin(); it!=X->Elist.end(); it++) (*it)->cell = X;
  for(std::list<Edge*>::iterator it=M->Elist.begin(); it!=M->Elist.end(); it++) (*it)->cell = M;

//  std::cout << M->Vlist.size() << " " << M->Elist.size() << " " << X->Vlist.size() << " " << X->Elist.size() << std::endl;

  //... adding new vertices to cells c1 and c2 ...
  std::list<Vertex*>::iterator it=c1->Vlist.begin(); // c1 will gain the v_e1m
  while(it!=c1->Vlist.end()) { if((*it)==e1X->tail){c1->Vlist.insert(it,v_e1m); break;} else {++it;} }
  it=c2->Vlist.begin(); // c2 will gain the v_e2m
  while(it!=c2->Vlist.end()) { if((*it)==e2->tail){c2->Vlist.insert(it,v_e2m); break;} else {++it;} }

  // ... adding new bonds to c1 and c2 ...
  std::list<Edge*>::iterator ite=c1->Elist.begin(); // c1 will lgain e1X_s
  while(ite!=c1->Elist.end()) { if((*ite)==e1_sn){c1->Elist.insert(ite,e1X_s); break;} else {++ite;} }
  ite=c2->Elist.begin(); // c2 will gain e2X_s
  while(ite!=c2->Elist.end()) { if((*ite)==e2->conj){c2->Elist.insert(ite,e2X_s); break;} else {++ite;} }

  //some necessary re-wiring
  e2->head = v_e2m; e2->conj->tail = v_e2m;
  e1->tail = v_e1m; e1->conj->head = v_e1m;
  e2_nsns->next = e2X_s; e2X_s->next = e2->conj;
  e1->conj->next = e1X_s; e1X_s->next = e1_sn;
  e1_snsns->next = e1X; e1X->next=eM_s; eM_s->next = e2X; e2X->next = e2_n;
  e2->next=eM; eM->next=e1;

  update_apical_geometric_quantities(); //TODO: do for only involved bonds and involved cells. Do the same for basal
  //std::cout <<"finals C1: " << c1->Vlist.size() << " " << c1->Elist.size() << "   finals C2: " << c2->Vlist.size() << " " << c2->Elist.size() << std::endl;

  return 1;
}


bool Tissue::divide_N_random_cell_with_random_axis(int N_d)
{
  auto rd = std::random_device {};
  auto rng = std::default_random_engine { rd() };

  if(N_d<1) return 0;

  std::vector<int> cell_indices(ListCell.size());
  for(size_t i=0; i<ListCell.size(); i++) cell_indices[i]=i;
  std::shuffle(std::begin(cell_indices), std::end(cell_indices), rng);

  int div_count = 0;
  size_t list_idx=0;
  while (div_count<N_d)
  {
    int cid = cell_indices[list_idx];
    bool flag = divide_cell_random_axis(&ListCell[cid]);
    if(flag) div_count++;

    list_idx++;

    if(list_idx==cell_indices.size())
    {
      cell_indices.resize(ListCell.size());
      list_idx = 0;
      for(size_t i=0; i<ListCell.size(); i++) cell_indices[i]=i;
      std::shuffle(std::begin(cell_indices), std::end(cell_indices), rng);
    }
  }

  return 1;
}

double Tissue::W_2D()
{
  W_frame = 0.;

  #if PRALLEL
  auto func = [this](auto&& it){double dA = ((&it)->area_a - this->CellMechProp.A0c); this->W_frame += dA*dA* this->CellMechProp.Kc/2.; };
  // __gnu_parallel::for_each(
  //   ListCell.begin(),
  //   ListCell.end(),
  //   func);
  std::for_each(
    std::execution::par_unseq,
    ListCell.begin(),
    ListCell.end(),
    func);
  #else
  for(std::vector<Cell>::iterator it=ListCell.begin();it!=ListCell.end();it++)
  {
    double dA = (it->area_a - CellMechProp.A0c);
    W_frame += dA*dA*CellMechProp.Kc/2.;
  }
  #endif
  // now the edge energy
  #if PRALLEL
  auto func2 = [this](auto&& it){this->W_frame += this->CellMechProp.BondT * (&it)->l_a/2.;};
  // __gnu_parallel::for_each(
  //   ListEdge.begin(),
  //   ListEdge.end(),
  //   func2);
  std::for_each(
    std::execution::par_unseq,
    ListEdge.begin(),
    ListEdge.end(),
    func2);
  #else
  for(std::vector<Edge>::iterator it=ListEdge.begin();it!=ListEdge.end();it++) W_frame += CellMechProp.BondT * it->l_a/2.;// half factor due to half edges
  #endif

  return W_frame;
}


double Tissue::W_3D()
{
  W_frame = 0.;
  for(std::vector<Cell>::iterator it=ListCell.begin();it!=ListCell.end();it++)
  {
    double dV = (it->vol - CellMechProp.V0c);
    W_frame += dV*dV*CellMechProp.Kc/2.;
    W_frame += CellMechProp.Ta * it->area_a + CellMechProp.Tb * it->area_b;
  }

  for(std::vector<Edge>::iterator it=ListEdge.begin();it!=ListEdge.end();it++) W_frame += CellMechProp.Tl * it->area_lateral/2.;// half factor due to half edges

  if(SphMechProp.P0 != 0 )
  {
    W_frame -= SphMechProp.P0 * lumen_volume; //  -P0 * V_l
  }
  else if(SphMechProp.Kl >0 )
  {
    W_frame += SphMechProp.Kl * (lumen_volume - SphMechProp.V0l) * (lumen_volume - SphMechProp.V0l);
  }

  return W_frame;
}

void Tissue::update_dW_a_2D()
{
  #if PRALLEL
  auto func = [](auto&& it){(&it)->dW_a = Vector3d(0.,0.,0.);};
  // __gnu_parallel::for_each(
  //   ListVertex.begin(),
  //   ListVertex.end(),
  //   func);
  std::for_each(
    std::execution::par_unseq,
    ListVertex.begin(),
    ListVertex.end(),
    func);
  #else
  for(std::vector<Vertex>::iterator itv=ListVertex.begin();itv!=ListVertex.end();itv++) itv->dW_a = Vector3d(0.,0.,0.);//first reset gradients
  #endif

  //next, go around cells and update gradients on all vertices
  //double M_EPS = 1e-24;
  for(std::vector<Cell>::iterator it=ListCell.begin();it!=ListCell.end();it++)
  {
    double KcxdA = CellMechProp.Kc * (it->area_a - CellMechProp.A0c);
    double Nv = (double)it->Vlist.size();
    for(std::list<Edge*>::const_iterator ite = it->Elist.begin(); ite != it->Elist.end(); ite++)
    {
      auto &e = *ite;
      Vector3d aa = e->head->pos_a;
      Vector3d bb = e->tail->pos_a; //prev vertex
      Vector3d cc = e->next->head->pos_a; //next vertex
      Vector3d mm = Nv * it->cent_a - aa; // X_cent - aa
      // Vector3d tt = (1. - Nv) * bb + mm;
      // Vector3d uu = (1. - Nv) * cc + mm;
      // Vector3d rr = mm.cross(bb) + aa.cross(tt);
      // Vector3d ss = mm.cross(cc) + aa.cross(uu);
      Vector3d grad_A = gradient_face_area(aa, bb, cc, mm, Nv);
      //Vector3d grad_A = tt.cross(rr)/(2.*Nv*rr.norm()+eps) + uu.cross(ss)/(2.*Nv*ss.norm()+eps); // gradient of area
      Vector3d l1 = e->l_a_vec/(e->l_a+NORM_EPS);
      Vector3d l2 = e->next->l_a_vec/(e->next->l_a+NORM_EPS);
      Vector3d grad_l = l1 - l2;
      e->head->dW_a += KcxdA * grad_A + CellMechProp.BondT * grad_l + CellMechProp.perim_elasticity * e->cell->perim_a * grad_l;

    }

  }
}

void Tissue::update_dW_3D_spherical()// will take lumen arguments later
{
  for(std::vector<Vertex>::iterator itv=ListVertex.begin();itv!=ListVertex.end();itv++)
  {
    itv->dW_a = Vector3d(0.,0.,0.); itv->dW_b = Vector3d(0.,0.,0.);
  }//first reset gradients
  for(std::vector<Cell>::iterator it=ListCell.begin();it!=ListCell.end();it++)
  {
    double KcxdV = CellMechProp.Kc * (it->vol - CellMechProp.V0c);
    double Nv = (double)it->Vlist.size();
    for(std::list<Edge*>::const_iterator ite = it->Elist.begin(); ite != it->Elist.end(); ite++)
    {
      auto &e = *ite;
      Vector3d aa = e->head->pos_a;
      Vector3d bb = e->tail->pos_a; //prev vertex
      Vector3d cc = e->next->head->pos_a; //next vertex
      Vector3d aa_prime = e->head->pos_b;
      Vector3d bb_prime = e->tail->pos_b; //prev vertex
      Vector3d cc_prime = e->next->head->pos_b; //next vertex
      Vector3d MM_c = 2.*Nv * it->centroid - aa; // centroid - aa
      Vector3d MM_c_prime = 2.*Nv * it->centroid - aa_prime; // centroid - aa_pr
      Vector3d mm_apical = Nv * it->cent_a - aa;
      Vector3d mm_basal = Nv * it->cent_b - aa_prime;
      Vector3d mm_p = 4. * e->cent_l - aa;
      Vector3d mm_p_prime = 4. * e->cent_l - aa_prime;
      Vector3d mm_q = 4. * e->next->cent_l - aa;
      Vector3d mm_q_prime = 4. * e->next->cent_l - aa_prime;
      Vector3d grad_apical_area = gradient_face_area(aa, bb, cc, mm_apical, Nv);
      Vector3d grad_basal_area = gradient_face_area(aa_prime, bb_prime, cc_prime, mm_basal, Nv);
      // gradient with respect to apical point a
      Vector3d grad_lateral_area_a =  gradient_face_area(aa, bb, aa_prime, mm_p, 4.)+gradient_face_area(aa, aa_prime, cc, mm_q, 4.);
      // gradient with respect to basal point a_prime
      Vector3d grad_lateral_area_b =  gradient_face_area(aa_prime, bb_prime, aa, mm_p_prime, 4.)+gradient_face_area(aa_prime, aa, cc_prime, mm_q_prime, 4.);
      // gradients of cell volume
      Vector3d grad_vol_a1 = gradient_volume(bb, cc, MM_c, mm_apical, Nv, Nv);
      Vector3d grad_vol_a2 = gradient_volume(aa_prime, bb, MM_c, mm_p, Nv, 4.);
      Vector3d grad_vol_a3 = gradient_volume(cc, aa_prime, MM_c, mm_q, Nv, 4.);
      Vector3d grad_vol_a = grad_vol_a1 + grad_vol_a2 + grad_vol_a3;
      Vector3d grad_vol_b = gradient_volume(cc_prime, bb_prime, MM_c_prime, mm_basal, Nv, Nv)+gradient_volume(bb_prime, aa, MM_c_prime, mm_p_prime, Nv, 4.)
       + gradient_volume(aa, cc_prime, MM_c_prime, mm_q_prime, Nv, 4.);
       //adding volume elasticity terms
      e->head->dW_a += KcxdV * grad_vol_a;
      e->head->dW_b += KcxdV * grad_vol_b;
      //adding surface tension terms
      e->head->dW_a += CellMechProp.Ta * grad_apical_area + 1/2. * CellMechProp.Tl * grad_lateral_area_a;
      e->head->dW_b += CellMechProp.Tb * grad_basal_area + 1/2. * CellMechProp.Tl * grad_lateral_area_b;
      if(SphMechProp.P0 != 0 )
      {
        e->head->dW_a += SphMechProp.P0 * grad_vol_a1; // gradient of -P0 * V_l; note that grad_vol_a1 comes with negative signe, thus the term is added with + sign
      }
      else if(SphMechProp.Kl >0 )
      {
        e->head->dW_a -= SphMechProp.Kl * (lumen_volume - SphMechProp.V0l) * grad_vol_a1; // grad of 1/2 Kl (Vl - V0l)^2
      }
      // Vector3d l1 = e->l_a_vec/(e->l_a+eps);
      // Vector3d l2 = e->next->l_a_vec/(e->next->l_a+eps);
      // Vector3d grad_l = l1 - l2;
      // e->head->dW_a += KcxdA * grad_A
      // + 0.5 * CellMechProp.BondT * grad_l
      // + CellMechProp.perim_elasticity * e->cell->perim_a * grad_l;

    }

  }
}

void Tissue::evolve_vertex_positions(double xi, double dt, bool spherical_constrain)
{
  #if PRALLEL
  auto func = [spherical_constrain,xi,dt](auto&& it){
    (&it)->pos_a_prev = (&it)->pos_a;
    Vector3d F_net = (&it)->F_active_a - (&it)->dW_a;
    if(spherical_constrain) F_net = project_on_tangent_plane((&it)->pos_a, F_net);
    (&it)->pos_a += F_net * (dt/xi);
  };
  std::for_each(
    std::execution::par_unseq,
    ListVertex.begin(),
    ListVertex.end(),
    func);
  #else
  for(std::vector<Vertex>::iterator itv=ListVertex.begin();itv!=ListVertex.end();itv++)
  {
    itv->pos_a_prev = itv->pos_a;
    Vector3d F_net = itv->F_active_a - itv->dW_a;
    if(spherical_constrain) F_net = project_on_tangent_plane(itv->pos_a, F_net);
    itv->pos_a += F_net * (dt/xi);
  }
  #endif
}


void Tissue::set_apical_points_on_sphere_at_radius(double R_sp)
{
  //double R_sp = std::sqrt ( ListCell.size()* CellMechProp.A0c/(4.*PI) );
  Vector3d T_cent = Vector3d(0,0,0);
  for(std::vector<Vertex>::iterator it = ListVertex.begin(); it != ListVertex.end(); it++)
  {
  T_cent += it->pos_a;
  }
  T_cent /= ListVertex.size();

  for(std::vector<Vertex>::iterator it = ListVertex.begin(); it != ListVertex.end(); it++)
  {
    Vector3d pt = it->pos_a - T_cent;
    it->pos_a = pt *(R_sp/(pt.norm()) );
  }
}

// void Tissue::update()
// {
//   update_apical_geometric_quantities();
//   W();
//   update_dW_a_2D();
// }

// void Tissue::set_pos()
// {
//   if(vpos.size()!=3*ListVertex.size()) throw;
//
//   for(int i=0; i<ListVertex.size(); i++)
//   {
//     ListVertex[i].pos_a[0] = vpos[3*i];
//     vpos[3*i]=ListVertex[i].pos_a[0];
//     ListVertex[i].pos_a[1] = vpos[3*i+1];
//     vpos[3*i+1]=ListVertex[i].pos_a[1];
//     ListVertex[i].pos_a[2] = vpos[3*i+2];
//     pos[3*i+2]=ListVertex[i].pos_a[2];
//   }
//
// };
//
// void Tissue::update_pos()
// {
//   if(vpos.size()!=3*ListVertex.size()) throw;
//
//   for(int i=0; i<ListVertex.size(); i++)
//   {
//     vpos[3*i]=ListVertex[i].pos_a[0];
//     vpos[3*i+1]=ListVertex[i].pos_a[1];
//     vpos[3*i+2]=ListVertex[i].pos_a[2];
//   }
//
// };
//
// void Tissue::update_pos_and_der()
// {
//   vpos = std::vector<double>(3*ListVertex.size());
//   vder = std::vector<double>(3*ListVertex.size());
//
//   update_dW_a_2D();
//
//   for(int i=0; i<ListVertex.size(); i++)
//   {
//     vpos[3*i]=ListVertex[i].pos_a[0];
//     vpos[3*i+1]=ListVertex[i].pos_a[1];
//     vpos[3*i+2]=ListVertex[i].pos_a[2];
//
//     vder[3*i]=ListVertex[i].dW_a[0];
//     vder[3*i+1]=ListVertex[i].dW_a[1];
//     vder[3*i+2]=ListVertex[i].dW_a[2];
//   }
//
// };
//
// // function to bracket the minimum (minimize Energy(T->pos + u*dir) for u)
// void Tissue::bracket(std::vector<double> &dir, double &ax, double &bx, double &cx, double &fb, double &current_u)
// {
//     const double GOLD = 1.618034; // default ratio by which successive intervals are magnified
//     const double GLIMIT = 100.0; //maximum magnification allowed for a parabolic-fit step
//     const double TINY = 1.0e-20;
//
//     double fa,fc,fu,u;
//
//     // evaluate energy at T->pos+ax*dir and at T->pos+bx*dir
//     for(unsigned long i=0;i<vpos.size();i++) vpos[i] += (ax-current_u)*dir[i];
//
//     set_pos();
//     W();
//     fa= W_frame;
//     current_u=ax;
//
//     for(unsigned long i=0;i<vpos.size();i++) vpos[i] += (bx-current_u)*dir[i];
//
//     set_pos();
//
//     fb=W();;
//     current_u=bx;
//
//     if(fb>fa) { // switch roles of a and b so that we can go downhill in the direction from a to b
//         SWAP(ax,bx);
//         SWAP(fa,fb);
//     }
//
//     cx = bx + GOLD*(bx-ax); // first guess for c - now bx is golden section between ax and cx
//     // evaluate energy at T->pos+cx*dir
//     for(unsigned long i=0;i<vpos.size();i++) vpos[i] += (cx-current_u)*dir[i];
//
//     set_pos();
//     fc=W();;
//     current_u=cx;
//
//     while(fb>fc) { // keep returning here until we bracket
//         double r=(bx-ax)*(fb-fc);
//         double q=(bx-cx)*(fb-fa);
//         u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0*SIGN(std::max(std::abs(q-r),TINY),q-r));
//         double ulim=bx+GLIMIT*(cx-bx);
//
//         // now test various possibilities
//         if((bx-u)*(u-cx)>0.0) { // parabolic u is between b and c: try it
//             // fu=f(u):
//             for(unsigned long i=0;i<vpos.size();i++) vpos[i] += (u-current_u)*dir[i];
//
//             set_pos();
//             fu=W();;
//             current_u=u;
//             if(fu<fc){ // got a minumum between b and c
//                 ax=bx;
//                 bx=u;
//                 fa=fb;
//                 fb=fu;
//                 return;
//             } else if(fu>fb){ // got a minimum between a and u
//                 cx=u;
//                 fc=fu;
//                 return;
//             }
//
//             u=cx+GOLD*(cx-bx); // parabolic fit was no use; use default magnification
//             for(unsigned long i=0;i<vpos.size();i++) vpos[i] += (u-current_u)*dir[i];
//
//             set_pos();
//             fu=W();
//             current_u=u;
//         } else if((cx-u)*(u-ulim)>0.0) { // parabolic fit is between c and its allowed limit
//             // fu=f(u):
//             for(unsigned long i=0;i<vpos.size();i++) vpos[i] += (u-current_u)*dir[i];
//
//             set_pos();
//             fu=W();
//             current_u=u;
//             if(fu<fc) {
//                 shft3(bx,cx,u,u+GOLD*(u-cx));
//                 fb=fc;
//                 fc=fu;
//                 for(unsigned long i=0;i<vpos.size();i++) vpos[i] += (u-current_u)*dir[i];
//
//                 set_pos();
//                 fu=W();
//                 current_u=u;
//             }
//
//         } else if((u-ulim)*(ulim-cx) >= 0.0) { // limit parabolic u to maximum allowed value
//             u=ulim;
//             for(unsigned long i=0;i<vpos.size();i++) vpos[i] += (u-current_u)*dir[i];
//
//             set_pos();
//             fu=W();
//             current_u=u;
//
//         } else { // reject parabolic u, use default magnification
//             u=cx+GOLD*(cx-bx);
//             for(unsigned long i=0;i<vpos.size();i++) vpos[i] += (u-current_u)*dir[i];
//
//             set_pos();
//             fu=W();
//             current_u=u;
//         }
//
//         shft3(ax,bx,cx,u); // eliminate oldest point and continue
//         shft3(fa,fb,fc,fu);
//
//     }
//
// }
//
// // minimize the energy along the line p+lambda*dir, the result is written into T->pos
// double Tissue::line_minimization_Brent(std::vector<double> &dir)
// {
//     const int ITMAX = minPar.ITMAX; // max number of iterations before the algorithm stops
//     const double CGOLD = 0.3819660; // golden ratio
//     const double ZEPS = std::numeric_limits<double>::epsilon()*10e-3; // machine limit
//     double a,b,d(0),etemp,fu,fv,fw,fx;
//     double p,q,r,tol1,tol2,u,v,w,x,xm;
//     double e=0; // this will be the distance on the step moved before last
//     double tol = 1e-14; //TO CHECK
//     update_pos(); //TODO; //update_pos(); //TODO
//     std::vector<double> pos_init = vpos; //TODO
//
//     std::vector<double> dir_vec = dir; // take copy of the direction vector to be able to modify it
//
//     // call Brents bracket method to enclose the minumum along the gradient
//     double ax(-0.0000001),bx(0.0000001), cx;
//     bracket(dir,ax,bx,cx,fx,u);
//
//     a = (ax<cx ? ax:cx); // a and b must be in ascending order
//     b = (ax>cx ? ax:cx);
//     x=w=v=bx; // initializations...
//
//     double u_previous=u;
//
//     fw = fv = fx; // energy at middle point of interval (bx)
//
//     for(int iter = 0; iter<ITMAX;iter++){ // main program loop
//
//         xm=0.5*(a+b);
//         tol2=2.0*(tol1=tol*std::abs(x)+ZEPS);
//         if(std::abs(x-xm)<=(tol2-0.5*(b-a))){ // test for done here
//             return 1;
//         }
//         if(std::abs(e)>tol1){ // construct a trial parabolic fit
//             r=(x-w)*(fx-fv);
//             q=(x-v)*(fx-fw);
//             p=(x-v)*q - (x-w)*r;
//             q=2.0*(q-r);
//             if(q>0) p=-p;
//             q=std::abs(q);
//             etemp=e;
//             e=d;
//             if(std::abs(p)>=std::abs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x)) d=CGOLD*(e=(x>=xm ? a-x:b-x));
//             // the above conditions determine the acceptability of the parabolic fit. Here we take the golden section step into the larger of the two segments
//             else{
//                 d=p/q;
//                 u=x+d;
//                 if(u-a<tol2 || b-u<tol2) d = SIGN(tol1,xm-x);
//             }
//         } else {
//             d=CGOLD*(e=(x>=xm ? a-x:b-x));
//         }
//
//         u= (std::abs(d)>=tol1 ? x+d:x+SIGN(tol1,d));
//
//         // now evaluate W for the only time in the iteration
//         for(unsigned long i=0;i<vpos.size();i++) { // move the tissue along dir
//             vpos[i] += (u-u_previous)*dir[i];
//         }
//
//         u_previous = u;
//
//         set_pos();
//         fu = W();
//
//         //std::cout << "a= " << a << ", u= " << u  << ", b=" << b << ", fu= " << fu << std::endl;
//
//         if(fu<=fx){
//             if(u>=x) a=x; else b=x;
//             shft3(v,w,x,u);
//             shft3(fv,fw,fx,fu);
//         } else {
//             if(u<x) a=u; else b=u;
//             if(fu<=fw || w==x) {
//                 v=w;
//                 w=u;
//                 fv=fw;
//                 fw=fu;
//             } else if(fu<=fv || v==x || v==w){
//                 v=u;
//                 fv=fu;
//             }
//         }
//     }
//     std::cout << "Too many iterations in brent" << std::endl;
//     throw("Too many iterations in brent");
// }
//
// // multidimensional minimization by the Fletcher-Reese-Polak-Ribiere method; algorithm following "Numerical Recipes"
// int Tissue::my_cg_minimizer()
// {
//   set_minimization_params();
//
//     unsigned step = 0;
//
//     const int ITMAX = minPar.ITMAX;
//     const double EPS = minPar.EPS;
//     const double GTOL = minPar.GTOL; // maximal gradient size as stop criterion
//     const double ftol = minPar.ftol;//1.0e-8; // relative max change in energy as stop criterion
//
//     double gg, dgg;
//     int iter;
//     double fret;
//
//
//     int rateT1check= minPar.rateT1check; // number of steps after each which the system is checked for possible T1s
//
//     update();
//     update_pos_and_der();
//     int n = vpos.size();
//     std::vector<double> g(n), h(n);
//
//     double fp = W_frame; // energy of starting point
//     std::vector<double> xi = vder; // derivative at starting point of minimization
//     std::vector<double> pp = vpos; // position of the starting point
//
//
//     for(int j=0; j<n; j++){ // initial values for gradient calculation
//         g[j]=-xi[j];
//         xi[j]=h[j]=g[j];
//     }
//
//
//
//     for(int its=0;its<ITMAX;its++){// loop over iterations
//
//         iter=its; // current iteration
//         line_minimization_Brent(xi); // do a line minimization along vector xi and the Tissue will be in this optimal position afterwards
//         update(); // update the derivatives
//         update_pos_and_der(); // save the derivatives in der
//
//         fret = W_frame;
//
//         if(2.0*std::abs(fret-fp) < ftol){ // one possible return - check absolute change
//         //if(2.0*std::abs(fret-f p) < ftol*(std::abs(fret)+std::abs(fp)+EPS)){ // one possible return - check for relative change
//             std::cout << "Minimization successful - difference between successive energies sufficiently low." << std::endl;
//             return 1;
//         }
//
//         fp=fret; // save previous energy fp
//
//         xi = vder;
//         double test = 0.0; // test for convergence on zero gradient
//         double den = std::max(fp, 1.0);
//         for(int j=0; j<n; j++){
//             double temp = std::abs(xi[j])*std::max(std::abs(vpos[j]),1.0)/den; // Standard way
//           //  double temp = std::abs(xi[j])*(n/6.)/den; // here: seems more useful to set something like: change/(energy per vertex) < GTOL   as a criterion since the tissue behaves translational invariant
//             if(temp>test) test=temp;
//         }
//
//         if(test<GTOL){ // the other possible return
//             std::cout << "Minimization successful - gradient short enough." << std::endl;
//             return 1;
//         }
//
//         dgg=gg=0;
//         for(int j=0; j<n; j++){ // calculate the conjugate gradient according to Fletcher-Reese or Polak-Ribiere
//             gg += g[j]*g[j];
//             if(minPar.MinimizationMode == 1) {
//                 dgg += xi[j]*xi[j]; // statement for Fletcher-Reese
//             } else if(minPar.MinimizationMode == 2) {
//                 dgg += (xi[j]+g[j])*xi[j]; // statement for Polak-Ribiere
//             }
//         }
//
//         if(gg==0){ // unlikely: if gradient is exactly zero, then we are already done
//         std::cout << "Minimization successful - gradient zero." << std::endl;
//             return 1;
//         }
//
//         double gam = dgg/gg;
//
//         for(int j=0; j<n; j++) {
//             g[j]=-xi[j];
//             xi[j]=h[j]=g[j]+gam*h[j];
//         }
//
//         step++;
//         // // check for T1 transitions
//         // if(rateT1check>0 && its>0 && its%rateT1check==0)
//         // {
//         //     bool a = T->popAllEdges();
//         //     bool b = T->unpopAllVertices();
//         //     if(a | b) return 2;
//         // }
//     }
//
//     //update();
//
//     std::cout << "Minimization did not succeed in predefined maximal number of iterations!" << std::endl;
//     return 0; // too many iterations in frpr - method
//
// }
//
//

double Tissue::gsl_f3D (const gsl_vector *v)
{
  for(size_t i=0; i<ListVertex.size(); i++)
  {
    ListVertex[i].pos_a[0] = gsl_vector_get(v, 6*i);
    ListVertex[i].pos_a[1] = gsl_vector_get(v, 6*i+1);
    ListVertex[i].pos_a[2] = gsl_vector_get(v, 6*i+2);
    ListVertex[i].pos_b[0] = gsl_vector_get(v, 6*i+3);
    ListVertex[i].pos_b[1] = gsl_vector_get(v, 6*i+4);
    ListVertex[i].pos_b[2] = gsl_vector_get(v, 6*i+5);
  }

  update_3D_geometric_quantities();

  W_3D();

  return W_frame;
}


double Tissue::gsl_f (const gsl_vector *v)
{
  // for(int i=0; i<ListVertex.size(); i++)
  // {
  //   gsl_vector_get(v, ListVertex[i].pos_a[0]);
  //   gsl_vector_get(v, ListVertex[i].pos_a[1]);
  //   gsl_vector_get(v, ListVertex[i].pos_a[2]);
  // }

  // std::vector<Vector3d> old_pos(ListVertex.size());
  // for(int i=0; i<ListVertex.size(); i++) old_pos[i]=ListVertex[i].pos_a;

  for(size_t i=0; i<ListVertex.size(); i++)
  {
    ListVertex[i].pos_a[0] = gsl_vector_get(v, 3*i);
    ListVertex[i].pos_a[1] = gsl_vector_get(v, 3*i+1);
    ListVertex[i].pos_a[2] = gsl_vector_get(v, 3*i+2);
  }

  update_apical_geometric_quantities();

  W_2D();

  // for(int i=0; i<ListVertex.size(); i++) ListVertex[i].pos_a=old_pos[i];
  // update_apical_geometric_quantities();

  return W_frame;
}
//

static double static_gsl_f3D(const gsl_vector* v, void* params) // proxy function
{
    Tissue* object = static_cast<Tissue*>(params);
    return object->gsl_f3D(v);
}

static double static_gsl_f(const gsl_vector* v, void* params) // proxy function
{
    Tissue* object = static_cast<Tissue*>(params);
    return object->gsl_f(v);
}
//
// //

void Tissue::gsl_df3D (const gsl_vector *v, gsl_vector *df)
{

  for(size_t i=0; i<ListVertex.size(); i++)
  {
    ListVertex[i].pos_a[0] = gsl_vector_get(v, 6*i);
    ListVertex[i].pos_a[1] = gsl_vector_get(v, 6*i+1);
    ListVertex[i].pos_a[2] = gsl_vector_get(v, 6*i+2);
    ListVertex[i].pos_b[0] = gsl_vector_get(v, 6*i+3);
    ListVertex[i].pos_b[1] = gsl_vector_get(v, 6*i+4);
    ListVertex[i].pos_b[2] = gsl_vector_get(v, 6*i+5);
  }

  update_3D_geometric_quantities();

  update_dW_3D_spherical();

  //std::cout << "df clled:  " << gsl_vector_get(df, 0) << " " << gsl_vector_get(df, 1) << " " << gsl_vector_get(df, 2) << std::endl;
  for(size_t i=0; i<ListVertex.size(); i++)
  {
    //if(i<4) std::cout << ListVertex[i].dW_a[0] << " " << ListVertex[i].dW_a[1] << " " << ListVertex[i].dW_a[2] << std::endl;
    gsl_vector_set(df, 6*i, ListVertex[i].dW_a[0]);
    gsl_vector_set(df, 6*i+1, ListVertex[i].dW_a[1]);
    gsl_vector_set(df, 6*i+2, ListVertex[i].dW_a[2]);
    gsl_vector_set(df, 6*i+3, ListVertex[i].dW_b[0]);
    gsl_vector_set(df, 6*i+4, ListVertex[i].dW_b[1]);
    gsl_vector_set(df, 6*i+5, ListVertex[i].dW_b[2]);
  }

}


// /* The gradient of f, df = (df/dx, df/dy). */
void Tissue::gsl_df (const gsl_vector *v, gsl_vector *df)
{

  // std::vector<Vector3d> old_pos(ListVertex.size());
  // for(int i=0; i<ListVertex.size(); i++) old_pos[i]=ListVertex[i].pos_a;

  for(size_t i=0; i<ListVertex.size(); i++)
  {
    ListVertex[i].pos_a[0] = gsl_vector_get(v, 3*i);
    ListVertex[i].pos_a[1] = gsl_vector_get(v, 3*i+1);
    ListVertex[i].pos_a[2] = gsl_vector_get(v, 3*i+2);
  }

  update_apical_geometric_quantities();

  update_dW_a_2D();

  //std::cout << "df clled:  " << gsl_vector_get(df, 0) << " " << gsl_vector_get(df, 1) << " " << gsl_vector_get(df, 2) << std::endl;
  for(size_t i=0; i<ListVertex.size(); i++)
  {
    //if(i<4) std::cout << ListVertex[i].dW_a[0] << " " << ListVertex[i].dW_a[1] << " " << ListVertex[i].dW_a[2] << std::endl;
    gsl_vector_set(df, 3*i, ListVertex[i].dW_a[0]);
    gsl_vector_set(df, 3*i+1, ListVertex[i].dW_a[1]);
    gsl_vector_set(df, 3*i+2, ListVertex[i].dW_a[2]);
  }

  // for(int i=0; i<ListVertex.size(); i++) ListVertex[i].pos_a=old_pos[i];
  // update_apical_geometric_quantities();
 // std::cout << "df clled:  " << gsl_vector_get(df, 0) << " " << gsl_vector_get(df, 1) << " " << gsl_vector_get(df, 2) << std::endl;
}
//
static void static_gsl_df3D(const gsl_vector* v, void* params, gsl_vector *df) // proxy function
{
    Tissue* object = static_cast<Tissue*>(params);
    return object->gsl_df3D(v, df);
}

static void static_gsl_df(const gsl_vector* v, void* params, gsl_vector *df) // proxy function
{
    Tissue* object = static_cast<Tissue*>(params);
    return object->gsl_df(v, df);
}
// // //
// /* Compute both f and df together. */

void Tissue::gsl_fdf3D (const gsl_vector *v, double *f, gsl_vector *df)
{
  *f = static_gsl_f3D(v, this);
  static_gsl_df3D(v, this, df);
}


void Tissue::gsl_fdf (const gsl_vector *v, double *f, gsl_vector *df)
{
  *f = static_gsl_f(v, this);
  static_gsl_df(v, this, df);
  // *f = gsl_f(x);
  // gsl_df(x, df);
}
//
static void static_gsl_fdf3D(const gsl_vector* v, void* params, double *f, gsl_vector *df) // proxy function
{
    Tissue* object = static_cast<Tissue*>(params);
    object->gsl_fdf3D(v, f, df);
}

static void static_gsl_fdf(const gsl_vector* v, void* params, double *f, gsl_vector *df) // proxy function
{
    Tissue* object = static_cast<Tissue*>(params);
    object->gsl_fdf(v, f, df);
}
// //

void Tissue::minimize_gsl3D()
{
  update_dW_3D_spherical();
  size_t iter = 0;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  T = gsl_multimin_fdfminimizer_conjugate_fr;
  gsl_vector *x;
  gsl_multimin_fdfminimizer *s;
  gsl_multimin_function_fdf w_func;

 GSLMultiMinFuncPointer fptr = static_gsl_f3D;
 GSLMultiMinDfPointer dfptr = static_gsl_df3D;
 GSLMultiMinFdfPointer fdfptr = static_gsl_fdf3D;

  w_func.n = 6 * ListVertex.size(); // x, y, and z Coordinates
  w_func.f = fptr; //static_gsl_f;
  w_func.df = dfptr;//&DW_PTR; //std::invoke(&Tissue::gsl_df);
  w_func.fdf = fdfptr;//&WDW_PTR; //std::invoke(&Tissue::gsl_fdf);
  w_func.params = this;

  /* Starting point, x = (5,7) */
  x = gsl_vector_alloc (6*ListVertex.size());
  for(size_t i=0; i<ListVertex.size(); i++)
  {
    gsl_vector_set (x, 6*i, ListVertex[i].pos_a[0]);
    gsl_vector_set (x, 6*i+1, ListVertex[i].pos_a[1]);
    gsl_vector_set (x, 6*i+2, ListVertex[i].pos_a[2]);
    gsl_vector_set (x, 6*i+3, ListVertex[i].pos_b[0]);
    gsl_vector_set (x, 6*i+4, ListVertex[i].pos_b[1]);
    gsl_vector_set (x, 6*i+5, ListVertex[i].pos_b[2]);
  }

  s = gsl_multimin_fdfminimizer_alloc (T, 6*ListVertex.size());

  gsl_multimin_fdfminimizer_set (s, &w_func, x, 0.01, 1e-16);

  do
    {
      update_3D_geometric_quantities();
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);

      //std::cout << "here " << s->f << " " << s->gradient << std::endl;

      if (status)
        break;

      status = gsl_multimin_test_gradient (s->gradient, 1e-8);

      std::cout << std::setprecision(16) <<  iter << ":  " << s->f << std::endl;
      if (status == GSL_SUCCESS)
      {
        printf ("Minimum found at \n");
        std::cout << iter << " " << s->f << std::endl;
      }

    }
  while (status == GSL_CONTINUE && iter < 2000);

  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);

}


void Tissue::minimize_gsl()
{
  update_dW_a_2D();
  size_t iter = 0;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  T = gsl_multimin_fdfminimizer_conjugate_fr;
  gsl_vector *x;
  gsl_multimin_fdfminimizer *s;
  gsl_multimin_function_fdf w_func;

  // const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  // gsl_multimin_fminimizer *s = NULL;
  // gsl_vector *ss, *x;
  // gsl_multimin_function minex_func;

  // gsl_function F;
  // F.function = gsl_f;
  // F.params = NULL:

  // Tissue* ptr0 = this;
  // std::function<double(gsl_vector *, void *)> ptr = [=](gsl_vector *v, void *){return ptr0->gsl_f(v);};
  // //std::function<double(const gsl_vector *)> ff = ptr;
  // gsl_function F = {
  //   [](gsl_vector *v, void *vf)->double {auto &f = *static_cast<std::function<double(gsl_vector *, void *)>*>(vf); return f(*v);},
  //    &ptr};

  //double (*F)(const gsl_vector *, void *) = &ptr;

  //gsl_function_pp<decltype(ptr)> Fp(ptr);
  // gsl_function_pp<double> Fp(ptr);
  // auto &F = *static_cast<gsl_function*>(&Fp);

 GSLMultiMinFuncPointer fptr = static_gsl_f;
 GSLMultiMinDfPointer dfptr = static_gsl_df;
 GSLMultiMinFdfPointer fdfptr = static_gsl_fdf;

//const gsl_vector* v, void* params, double *f, gsl_vector *df
// std::function<void(const gsl_vector* , void* , gsl_vector *)> df_func = static_gsl_df;
// std::function<void(const gsl_vector* , void* , double *, gsl_vector *)> fdf_func = static_gsl_fdf;

  w_func.n = 3 * ListVertex.size(); // x, y, and z Coordinates
  w_func.f = fptr; //static_gsl_f;
  w_func.df = dfptr;//&DW_PTR; //std::invoke(&Tissue::gsl_df);
  w_func.fdf = fdfptr;//&WDW_PTR; //std::invoke(&Tissue::gsl_fdf);
  w_func.params = this;

  /* Starting point, x = (5,7) */
  x = gsl_vector_alloc (3*ListVertex.size());
  for(size_t i=0; i<ListVertex.size(); i++)
  {
    gsl_vector_set (x, 3*i, ListVertex[i].pos_a[0]);
    //gsl_vector_set (df, 3*i, ListVertex[i].dW_a[0]);
    gsl_vector_set (x, 3*i+1, ListVertex[i].pos_a[1]);
    //gsl_vector_set (df, 3*i+1, ListVertex[i].dW_a[1]);
    gsl_vector_set (x, 3*i+2, ListVertex[i].pos_a[2]);
    //gsl_vector_set (df, 3*i+2, ListVertex[i].dW_a[2]);
  }

  s = gsl_multimin_fdfminimizer_alloc (T, 3*ListVertex.size());

  gsl_multimin_fdfminimizer_set (s, &w_func, x, 0.1, 1e-16);

  do
    {
      update_apical_geometric_quantities();
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);

      //std::cout << "here " << s->f << " " << s->gradient << std::endl;

      if (status)
        break;

      status = gsl_multimin_test_gradient (s->gradient, 1e-8);

      std::cout << std::setprecision(16) <<  iter << ":  " << s->f << std::endl;
      if (status == GSL_SUCCESS)
      {
        printf ("Minimum found at \n");
        std::cout << iter << " " << s->f << std::endl;
      }

    }
  while (status == GSL_CONTINUE && iter < 2000);

  // for(int i=0; i<ListVertex.size(); i++)
  // {
  //   ListVertex[i].pos_a[0] = gsl_vector_get(x, 3*i);
  //   ListVertex[i].pos_a[1] = gsl_vector_get(x, 3*i+1);
  //   ListVertex[i].pos_a[2] = gsl_vector_get(x, 3*i+2);
  // }
  // update_apical_geometric_quantities();

  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);

}
//
//
// void Tissue::minimize_gsl2()
// {
//
//   const gsl_multimin_fminimizer_type *T =
//     gsl_multimin_fminimizer_nmsimplex2;
//   gsl_multimin_fminimizer *s = NULL;
//   gsl_vector *ss, *x;
//   gsl_multimin_function minex_func;
//
//   size_t iter = 0;
//   int status;
//   double size;
//
//   /* Starting point */
//   x = gsl_vector_alloc (3*ListVertex.size());
//   for(int i=0; i<ListVertex.size(); i++)
//   {
//     gsl_vector_set (x, 3*i, ListVertex[i].pos_a[0]);
//     gsl_vector_set (x, 3*i+1, ListVertex[i].pos_a[1]);
//     gsl_vector_set (x, 3*i+2, ListVertex[i].pos_a[2]);
//   }
//
//   /* Set initial step sizes to .1 */
//   ss = gsl_vector_alloc (3*ListVertex.size());
//   gsl_vector_set_all (ss, 0.1);
//
//
//   // Tissue* ptr2 = this;
//   // auto ptr = [=](const gsl_vector * x)->double{return ptr2->gsl_f(x);};
//
//
//   // std::cout << type_name<decltype(ptr)>() << std::endl;
//
//   //gsl_function_pp<decltype(ptr)> Fp(ptr);
//   // gsl_function_pp<const gsl_vector> Fp(ptr);
//
//   // gsl_function *F = static_cast<gsl_function*>(&Fp);
//
// /* Initialize method and iterate */
//   minex_func.n = 3*ListVertex.size();
//   minex_func.f = static_gsl_f;//F->function;//static_gsl_f;
//   minex_func.params = &(*this);
//
//   s = gsl_multimin_fminimizer_alloc (T, 3*ListVertex.size());
//   gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
//
//   do
//     {
//       update_apical_geometric_quantities();
//       iter++;
//       status = gsl_multimin_fminimizer_iterate(s);
//
//       if (status)
//         break;
//
//       size = gsl_multimin_fminimizer_size (s);
//       status = gsl_multimin_test_size (size, 1e-14);
//
//       // if (status == GSL_SUCCESS)
//       //   {
//       //     printf ("converged to minimum at\n");
//       //   }
//       //
//       // printf ("%5d f() = %7.3f size = %.3f\n",
//       //         iter, s->fval, size);
//     }
//   while (status == GSL_CONTINUE && iter < 10000);
//
//   printf ("%5d f() = %7.3f size = %.8f\n",
//           iter, s->fval, size);
//
//   for(int i=0; i<ListVertex.size(); i++)
//   {
//     ListVertex[i].pos_a[0] = gsl_vector_get(s->x, 3*i);
//     ListVertex[i].pos_a[1] = gsl_vector_get(s->x, 3*i+1);
//     ListVertex[i].pos_a[2] = gsl_vector_get(s->x, 3*i+2);
//   }
//
//   gsl_vector_free(x);
//   gsl_vector_free(ss);
//   gsl_multimin_fminimizer_free (s);
//
//
// }

void Tissue::randomize_cell_planar_polarity()
{
	std::cout << "random polarity direction is set!" << std::endl;
	std::random_device rd;  // Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis_tta(0., M_PI);
	std::uniform_real_distribution<> dis_phi(0., 2.*M_PI);
	std::uniform_real_distribution<> dis_a(-1., 1.);
	std::uniform_real_distribution<> dis_b(-1., 1.);

	for(std::vector<Cell>::iterator it=ListCell.begin();it!=ListCell.end();it++)
	{
		double a = dis_a(gen); double b = dis_b(gen);
		double theta = dis_tta(gen); double phi = dis_phi(gen);

		//...
		Vector3d theta_hat( cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta) );
		Vector3d phi_hat( -sin(phi), cos(phi), 0.0 );
		//..

		it->polarity = a*theta_hat + b*phi_hat ;
    //it->polarity = project_on_tangent_plane(it->cent_a, it->polarity);
    Vector3d p_t = project_on_tangent_plane(it->cent_a, it->polarity);
    it->polarity = p_t/(p_t.norm()+NORM_EPS);
	}

	for(std::vector<Vertex>::iterator itv=ListVertex.begin();itv!=ListVertex.end();itv++)
	{
		itv->pos_a_prev = itv->pos_a;
	}

}

void Tissue::set_axissymmetric_polarity(double a0, double b0)
{
  std::cout << "axissyymetric polarity direction is set!  a0="<< a0 << "   b0=" << b0 << std::endl;

  for(std::vector<Cell>::iterator it=ListCell.begin();it!=ListCell.end();it++)
  {
    Vector3d cent_sph_coord = cartesian_to_spherical(it->cent_a);
    double theta = cent_sph_coord[1]; double phi = cent_sph_coord[2];
    //...
    Vector3d theta_hat( cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta) );
    Vector3d phi_hat( -sin(phi), cos(phi), 0.0 );
    //..

    double a_theta = a0 + b0 * std::log( (1.-std::cos(theta))/(1.+std::cos(theta)+NORM_EPS) );

    it->polarity = std::sin(a_theta)*theta_hat + std::cos(a_theta)*phi_hat ;
    //it->polarity = project_on_tangent_plane(it->cent_a, it->polarity);
    Vector3d p_t = project_on_tangent_plane(it->cent_a, it->polarity);
    it->polarity = p_t/(p_t.norm()+NORM_EPS);
  }

  for(std::vector<Vertex>::iterator itv=ListVertex.begin();itv!=ListVertex.end();itv++)
  {
    itv->pos_a_prev = itv->pos_a;
  }

}

void Tissue::update_average_nearest_neighbor_polarity(bool _elastic_polarity)
{
  #if PRALLEL
  auto func = [_elastic_polarity](auto&& it){
    (&it)->avg_polarity_j=Vector3d(0,0,0);
    if(_elastic_polarity)
    {
      Vector3d p_t = project_on_tangent_plane((&it)->cent_a, (&it)->polarity);
      (&it)->polarity = p_t/(p_t.norm()+NORM_EPS);
    }
  };
  std::for_each(
    std::execution::par_unseq,
    ListCell.begin(),
    ListCell.end(),
    func);

  auto func2 = [](auto&& it){
    if ( (&it)->cell!=NULL && (&it)->conj->cell!=NULL )
    {
      (&it)->cell->avg_polarity_j += (&it)->conj->cell->polarity/(&it)->cell->Elist.size();
    }
    };
  std::for_each(
    std::execution::par_unseq,
    ListEdge.begin(),
    ListEdge.end(),
      func2);

  auto func3 = [](auto&& it){ (&it)->avg_polarity_j = project_on_tangent_plane((&it)->cent_a, (&it)->avg_polarity_j);};
  std::for_each(
    std::execution::par_unseq,
    ListCell.begin(),
    ListCell.end(),
    func3);
  #else
	for(std::vector<Cell>::iterator itc=ListCell.begin();itc!=ListCell.end();itc++)
	{
		itc->avg_polarity_j=Vector3d(0,0,0);
		//---------------- also make sure cell polarities are tangential---------------
    if(_elastic_polarity)
    {
      Vector3d p_t = project_on_tangent_plane(itc->cent_a, itc->polarity);
      itc->polarity = p_t/(p_t.norm()+NORM_EPS);
    }
	}

	for(std::vector<Edge>::iterator ite = ListEdge.begin();ite!=ListEdge.end();ite++)
	{
		if(ite->cell!=NULL && ite->conj->cell!=NULL)
		{
			ite->cell->avg_polarity_j += ite->conj->cell->polarity/ite->cell->Elist.size();
		}
	}

  for(std::vector<Cell>::iterator itc=ListCell.begin();itc!=ListCell.end();itc++)
  {
    itc->avg_polarity_j = project_on_tangent_plane(itc->cent_a, itc->avg_polarity_j);
    //itc->avg_polarity_j = p_t/(p_t.norm()+NORM_EPS);
  }
  #endif
}

void Tissue::evolve_polarity(double _P_nu, double _P_gamma, double diffD, double dt, double _P_A, bool _elastic_polarity)
{
  std::random_device rd;
  std::default_random_engine generator;
  generator.seed( rd() ); //Now this is seeded differently each time.
  std::normal_distribution<double> distribution(0.0,1.0);
  //...............................................................
  double Ap = _P_A * dt;
  //...............................................................
	update_average_nearest_neighbor_polarity(_elastic_polarity);

	for(std::vector<Cell>::iterator itc=ListCell.begin();itc!=ListCell.end();itc++)
	{
    //... noise vector on the polarity ....
    double eta_noise = distribution(generator);

    Vector3d t_perp = itc->cent_a.cross(itc->polarity);
    t_perp /= (t_perp.norm()+NORM_EPS);
    Vector3d noise_polarity = t_perp * (eta_noise * std::sqrt(2. * diffD * dt) );
    //...  cell center veloicty for flow alignment ....
		std::vector<Vector3d> pos_t1, pos_t2;
		std::tie(pos_t1, pos_t2) = get_vertex_pos_changes_of_cell(&(*itc));
    //... rotating polarity with cell's solid-body rotation axis Omega
		Vector3d Gamma, Omega;
		Matrix3d MM;
		std::tie(Gamma, MM) = get_ang_momentum_inertia_tensor(pos_t1, pos_t2, dt);
		Omega = MM.inverse() * Gamma;
		//if(Omega.norm()>0.1) std::cout << "Omega: " << Omega.norm()/PI << std::endl;
		itc->polarity = rotation(itc->polarity, -Omega.norm()*dt, Omega); // notice the negative sign of angle
    //if(itc->polarity.dot(itc->cent_a)>1e-8) {std::cout << itc->polarity.dot(itc->cent_a) << std::endl; }
    //Vector3d p1 = itc->polarity;
		itc->polarity +=  itc->avg_polarity_j * (_P_gamma * dt) + noise_polarity + (itc->cent_a - itc->cent_a_prev) * _P_nu;
    itc->polarity = project_on_tangent_plane(itc->cent_a, itc->polarity);
    double p_norm = itc->polarity.norm();

    if(_elastic_polarity){itc->polarity += Ap * (1 - p_norm*p_norm) * itc->polarity;}
    else{itc->polarity /= (p_norm+NORM_EPS);}
	}
}

void Tissue::update_active_forces()
{
  #if PRALLEL
  std::for_each( std::execution::par_unseq, ListVertex.begin(), ListVertex.end(), [](auto &&it) {(&it)->F_active_a = Vector3d(0,0,0);});

  auto func = [](auto&& it){double mc = (double)(&it)->cell->Elist.size(); (&it)->head->F_active_a += (&it)->cell->polarity * ((&it)->cell->mag_F_a/mc);};
  std::for_each(
    std::execution::par_unseq,
    ListEdge.begin(),
    ListEdge.end(),
    func);

  std::for_each( std::execution::par_unseq, ListVertex.begin(), ListVertex.end(),
  [](auto &&it) {(&it)->F_active_a = project_on_tangent_plane((&it)->pos_a, (&it)->F_active_a);});
  #else
  for(std::vector<Vertex>::iterator itv=ListVertex.begin();itv!=ListVertex.end();itv++) {itv->F_active_a = Vector3d(0,0,0);}

  for(std::vector<Edge>::iterator it=ListEdge.begin();it!=ListEdge.end();it++)
  {
    double mc = (double)it->cell->Elist.size();
    it->head->F_active_a += it->cell->polarity * (it->cell->mag_F_a/mc);
  }

  for(std::vector<Vertex>::iterator itv=ListVertex.begin();itv!=ListVertex.end();itv++){itv->F_active_a = project_on_tangent_plane(itv->pos_a, itv->F_active_a);}
  #endif


}

void Tissue::evolve_organoid_one_step(double _P_nu, double _P_gamma, double _diffD, double _time_step, double _xi, double _P_A, bool _elastic_polarity, bool _spherical_constrain=true)
{
    update_apical_geometric_quantities();
    evolve_polarity(_P_nu, _P_gamma, _diffD, _time_step, _P_A, _elastic_polarity);
    update_dW_a_2D();
    update_active_forces();
    evolve_vertex_positions(_xi, _time_step, _spherical_constrain);
    rescale_spheroid_apical_size();
}

double Tissue::get_tissue_polarity_elastic_energy()
{
  update_average_nearest_neighbor_polarity(true);

  double polrity_elastic_energy = 0.;

  #if PRALLEL
  for(std::vector<Cell>::iterator itc=ListCell.begin();itc!=ListCell.end();itc++)
  {
    polrity_elastic_energy += (itc->avg_polarity_j).dot(itc->polarity);
  }
  // auto func = [polrity_elastic_energy](auto&& it){ polrity_elastic_energy += (&it)->avg_polarity_j.dot((&it)->polarity);};
  // std::for_each(
  //   std::execution::par_unseq,
  //   ListCell.begin(),
  //   ListCell.end(),
  //   func);
  #else
	for(std::vector<Cell>::iterator itc=ListCell.begin();itc!=ListCell.end();itc++)
	{
		polrity_elastic_energy += (itc->avg_polarity_j).dot(itc->polarity);
	}
  #endif

  return polrity_elastic_energy/2.;

}

void Tissue::write_to_json_file(std::string output_folder, int frame)
{
  using json = nlohmann::json;
  std::string json_file = output_folder + "frame_" + std::to_string(frame) + ".json";
  std::ofstream json_out(json_file);
  json json_data;
  //   preparing data ...
  std::vector<Vector3d> vertex_coord_a(ListVertex.size());
  std::vector<Vector3d> vertex_prev_coord_a(ListVertex.size());
  std::vector<Vector3d> vertex_coord_b(ListVertex.size());
  std::vector<Vector3d> cell_polarity(ListCell.size());
  std::vector<double> cell_Factive_mag(ListCell.size());
  std::vector<std::vector<int>> cell_vertex_ids(ListCell.size());
  std::vector<std::vector<double>> cells_Q3d; cells_Q3d.reserve(ListCell.size());
  for(size_t i=0; i<ListVertex.size(); i++)
  {
    vertex_coord_a[i]=ListVertex[i].pos_a; vertex_coord_b[i]=ListVertex[i].pos_b; vertex_prev_coord_a[i]=ListVertex[i].pos_a_prev;
  }
  for(size_t i=0; i<ListCell.size(); i++)
  {
    Matrix3d q3d = get_apical_elongation_tensor( &ListCell[i] );
    std::vector<double> cell_q3d_vec={q3d(0,0), q3d(0,1), q3d(0,2), q3d(1,0), q3d(1,1), q3d(1,2), q3d(2,0), q3d(2,1), q3d(2,2)};
    cells_Q3d.emplace_back( cell_q3d_vec );
    cell_polarity[i]=ListCell[i].polarity; cell_Factive_mag[i]=ListCell[i].mag_F_a;
    std::vector<int> cell_v_ids; cell_v_ids.reserve(ListCell[i].Vlist.size());
    for(auto v : ListCell[i].Vlist) cell_v_ids.emplace_back(v->id);
    cell_vertex_ids[i] = cell_v_ids;
  }
  // .. writing them json
  json_data["Num_cells"] = ListCell.size();
  json_data["Num_edges"] = ListEdge.size();
  json_data["Num_vertices"] = ListVertex.size();
  json_data["vertex_pos_a"] = vertex_coord_a;
  json_data["vertex_pos_b"] = vertex_coord_b;
  json_data["vertex_prev_pos_a"] = vertex_prev_coord_a;
  json_data["cells_polarity"] = cell_polarity;
  json_data["cell_Factive_mag"] = cell_Factive_mag;
  json_data["cell_vertex_ids"] = cell_vertex_ids;
  json_data["cell_Q3d_tensor"] = cells_Q3d;
  json_out  << json_data << std::endl;

  json_out.close();
}


void Tissue::write_to_txt_files(std::string output_folder, int frame)
{
	// if(!boost::algorithm::ends_with(output_folder,"/")) output_folder = output_folder + "/";
	// std::cout << "frame:  " << frame  << std::endl;
	std::string vertices_file = output_folder + "vertices_" + std::to_string(frame) + ".txt";
	std::string cells_file = output_folder + "cells_" + std::to_string(frame) + ".txt";
  std::string cell_data_file = output_folder + "cell_data_" + std::to_string(frame) + ".txt";

	std::ofstream vertex_out; vertex_out.open(vertices_file);
	std::ofstream cells_out; cells_out.open(cells_file);
  std::ofstream cell_data_out; cell_data_out.open(cell_data_file);

	//-----------
	if(vertex_out.is_open())
	{
		for(std::vector<Vertex>::iterator itv=ListVertex.begin();itv!=ListVertex.end();itv++)
		{
			vertex_out << itv->pos_a[0] << " " << itv->pos_a[1] << " " << itv->pos_a[2] <<
			" " << itv->pos_b[0] << " " << itv->pos_b[1] << " " << itv->pos_b[2] <<
      " " << itv->pos_a_prev[0] << " " << itv->pos_a_prev[1] << " " << itv->pos_a_prev[2] << std::endl;
		}
	}
	//-----------
	if(cells_out.is_open() && cell_data_out.is_open())
	{
		for(std::vector<Cell>::iterator itc=ListCell.begin();itc!=ListCell.end();itc++)
		{
      cell_data_out << itc->Vlist.size() << " " << itc->polarity[0] << " " << itc->polarity[1] << " " << itc->polarity[2] << " " << itc->mag_F_a << std::endl;
      // now writing cell coordinates:
			for(std::list<Vertex*>::iterator itv=itc->Vlist.begin(); itv!=itc->Vlist.end(); itv++)
			{
				cells_out << (*itv)->id << " ";
			}
			cells_out << std::endl;
		}
	}
	//-----------
	vertex_out.close(); cells_out.close(); cell_data_out.close();
}
