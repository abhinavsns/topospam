#include "Triangulation.hpp"

//Triangulation constructor
Triangulation::Triangulation(){}

Triangulation::Triangulation(bool VSH_analysis_) : VSH_analysis(VSH_analysis_){
}

// void Vertex::setNetwork(Network *n) {
//   _net = n;
// }

//Vertex destructor
Triangulation::~Triangulation(){
}

TriVertex* Triangulation::append_vertex(Vector3d pos_)
{
  // iterate through all the vertices in the system and check if the vertex has the same coordinates
  std::vector<TriVertex>::iterator it = Vertices.begin();
  for(;it!=Vertices.end();it++)
  {
    TriVertex* v = &(*it);
    if((v->pos-pos_).norm() < 1e-10) return v;
  }
  //if it makes it this far, the vertes does not exist, and needs to be made
  Vertices.emplace_back(TriVertex(pos_, Vertices.size()));

  return &(Vertices.back());
}


TriEdge* Triangulation::append_edge(TriVertex *tail_, TriVertex *head_)
{
  TriEdge *e_conj = nullptr;
  // iterate through all the edges in the system and check if another Edge has the same vertices
  std::vector<TriEdge>::iterator it = Edges.begin();
  for(;it!=Edges.end();it++)
  {
    TriEdge* e = &(*it);
    if(e->tail == tail_ && e->head == head_) return e;
    else if(e->head == tail_ && e->tail == head_) e_conj = e;
  }
  //if it makes it this far, the edge does not exist, and needs to be made
  //int eid = (int)(Edges.size());
  Edges.emplace_back(TriEdge(tail_, head_, e_conj, Edges.size()) );
  return &(Edges.back());
}


void Triangulation::update_Vertices_sph()
{
  for(size_t i=0; i<Vertices.size(); i++)
  {
    Vertices[i].pos_sph = cartesian_to_spherical(Vertices[i].pos);
    Vertices[i].cos_theta = std::cos((Vertices[i].pos_sph)[1]);
    std::tie(Vertices[i].r_hat, Vertices[i].theta_hat, Vertices[i].phi_hat) = get_spherical_coord_unit_vectors((Vertices[i].pos_sph)[1], (Vertices[i].pos_sph)[2]);
  }
}

void Triangulation::update_Vertices_vec_field()
{
  for(size_t i=0; i<Vertices.size(); i++)
  {
    Vertices[i].vec_field = Vertices[i].pos - Vertices[i].r_hat;
  }
}

void Triangulation::update_Vertices_vec_field(std::vector<Vector3d> given_field)
{
  for(size_t i=0; i<Vertices.size(); i++)
  {
    Vertices[i].vec_field = given_field[i];
  }
}

void Triangulation::update_Edges_geometry()
{
  for(size_t i=0; i<Edges.size(); i++)
  {
    Edges[i].l_vec = Edges[i].head->pos - Edges[i].tail->pos;
    Edges[i].l = (Edges[i].l_vec).norm();
  }
}

void Triangulation::update_Triangles_area_and_normal_vec()
{
  for(size_t i=0; i<Triangles.size(); i++)
  {
    Triangles[i].normal_vec = 0.5 * (Triangles[i].e1->l_vec).cross(Triangles[i].e3->l_vec);
    Triangles[i].area = (Triangles[i].normal_vec).norm();
    Triangles[i].unit_normal_vec = Triangles[i].normal_vec/(Triangles[i].area+NORM_EPS);
  }
}

void Triangulation::update_Triangles_geomtery()
{
  update_Triangles_area_and_normal_vec();
  for(size_t i=0; i<Triangles.size(); i++)
  {
    Triangles[i].com = (Triangles[i].v1->pos + Triangles[i].v2->pos + Triangles[i].v3->pos)/3.;
    Triangles[i].update_excess();
    Triangles[i].update_corner_angles();
  }
}

void Triangulation::update_all_geometries()
{
  update_Vertices_sph();
  update_Vertices_vec_field();
  update_Edges_geometry();
  update_Triangles_geomtery();
}

void Triangulation::update_all_geometries(std::vector<Vector3d> given_field)
{
  update_Vertices_sph();
  update_Vertices_vec_field(given_field);
  update_Edges_geometry();
  update_Triangles_geomtery();
}

std::tuple<double,double,double> Triangulation::get_Elm_coefficients(int l, int m)//return E^r, E^1, E^2 for mode l,m
{
  double Er=0.; double E1=0.; double E2=0.;
  Vector3d PSIlm, PHIlm;
  for(size_t i=0; i<Triangles.size(); i++)
  {
    double triExcess = Triangles[i].excess;
    //vertex v1
    Vector3d EE = Triangles[i].v1->vec_field;
    double Ylm = realYlm(l,m, Triangles[i].v1->cos_theta, (Triangles[i].v1->pos_sph)[2]);
    std::tie(PSIlm, PHIlm) = PSIlmPHIlm(l, m, Triangles[i].v1->cos_theta, (Triangles[i].v1->pos_sph)[2], Triangles[i].v1->theta_hat, Triangles[i].v1->phi_hat);
    Er +=  EE.dot(Triangles[i].v1->r_hat) * Ylm * triExcess;
    E1 += EE.dot(PSIlm) * triExcess;
    E2 += EE.dot(PHIlm) * triExcess;
    //vertex v2
    EE = Triangles[i].v2->vec_field;
    Ylm = realYlm(l,m, Triangles[i].v2->cos_theta, (Triangles[i].v2->pos_sph)[2]);
    std::tie(PSIlm, PHIlm) = PSIlmPHIlm(l, m, Triangles[i].v2->cos_theta, (Triangles[i].v2->pos_sph)[2], Triangles[i].v2->theta_hat, Triangles[i].v2->phi_hat);
    Er +=  EE.dot(Triangles[i].v2->r_hat) * Ylm * triExcess;
    E1 += EE.dot(PSIlm) * triExcess;
    E2 += EE.dot(PHIlm) * triExcess;
    //vertex v3
    EE = Triangles[i].v3->vec_field;
    Ylm = realYlm(l,m, Triangles[i].v3->cos_theta, (Triangles[i].v3->pos_sph)[2]);
    std::tie(PSIlm, PHIlm) = PSIlmPHIlm(l, m, Triangles[i].v3->cos_theta, (Triangles[i].v3->pos_sph)[2], Triangles[i].v3->theta_hat, Triangles[i].v3->phi_hat);
    Er +=  EE.dot(Triangles[i].v3->r_hat) * Ylm * triExcess;
    E1 += EE.dot(PSIlm) * triExcess;
    E2 += EE.dot(PHIlm) * triExcess;
  }
  Er = Er/3.;//divide by three to divide excess to all vertices equally
  if(l==0) {E1=0; E2=0;}
  else{ E1 = E1/(3.*l*(l+1.)); E2 = E2/(3.*l*(l+1.));}

  return std::make_tuple(Er,E1,E2);
}

std::tuple<double, double, double, std::complex<double>, std::complex<double>, std::complex<double> > Triangulation::get_Elm_Clm_coefficients(int l, int m)
{
  double Er=0.; double E1=0.; double E2=0.;
  std::complex<double> Cr(0,0), C1(0,0), C2(0,0);
  std::complex<double> pf;//to fill in the real to complex convertor pre-factor
  double temp_r, temp_1, temp_2;

  Vector3d PSIlm, PHIlm;
  for(size_t i=0; i<Triangles.size(); i++)
  {
    double triExcess = Triangles[i].excess;
    //vertex v1
    Vector3d EE = Triangles[i].v1->vec_field;
    double Ylm = realYlm(l,m, Triangles[i].v1->cos_theta, (Triangles[i].v1->pos_sph)[2]);
    std::tie(PSIlm, PHIlm) = PSIlmPHIlm(l, m, Triangles[i].v1->cos_theta, (Triangles[i].v1->pos_sph)[2], Triangles[i].v1->theta_hat, Triangles[i].v1->phi_hat);
    pf = real_to_complex_prefactor(m, (Triangles[i].v1->pos_sph)[2]);
    temp_r = EE.dot(Triangles[i].v1->r_hat) * Ylm * triExcess;
    Er +=  temp_r; Cr += temp_r * pf;
    temp_1 = EE.dot(PSIlm) * triExcess;;
    E1 += temp_1; C1 += temp_1 * pf;
    temp_2 = EE.dot(PHIlm) * triExcess;
    E2 += temp_2; C2 += temp_2*pf;
    //vertex v2
    EE = Triangles[i].v2->vec_field;
    Ylm = realYlm(l,m, Triangles[i].v2->cos_theta, (Triangles[i].v2->pos_sph)[2]);
    std::tie(PSIlm, PHIlm) = PSIlmPHIlm(l, m, Triangles[i].v2->cos_theta, (Triangles[i].v2->pos_sph)[2], Triangles[i].v2->theta_hat, Triangles[i].v2->phi_hat);
    pf = real_to_complex_prefactor(m, (Triangles[i].v2->pos_sph)[2]);
    temp_r = EE.dot(Triangles[i].v2->r_hat) * Ylm * triExcess;
    Er += temp_r; Cr += temp_r*pf;
    temp_1 = EE.dot(PSIlm) * triExcess;
    E1 += temp_1; C1 += temp_1*pf;
    temp_2 = EE.dot(PHIlm) * triExcess;
    E2 += temp_2; C2 += temp_2*pf;
    //vertex v3
    EE = Triangles[i].v3->vec_field;
    Ylm = realYlm(l,m, Triangles[i].v3->cos_theta, (Triangles[i].v3->pos_sph)[2]);
    std::tie(PSIlm, PHIlm) = PSIlmPHIlm(l, m, Triangles[i].v3->cos_theta, (Triangles[i].v3->pos_sph)[2], Triangles[i].v3->theta_hat, Triangles[i].v3->phi_hat);
    pf = real_to_complex_prefactor(m, (Triangles[i].v3->pos_sph)[2]);
    temp_r = EE.dot(Triangles[i].v3->r_hat) * Ylm * triExcess;
    Er += temp_r; Cr += temp_r*pf;
    temp_1 = EE.dot(PSIlm) * triExcess;
    E1 += temp_1; C1 += temp_1*pf;
    temp_2 = EE.dot(PHIlm) * triExcess;
    E2 += temp_2; C2 += temp_2*pf;
  }
  Er = Er/3.; Cr =Cr/3.;//divide by three to divide excess to all vertices equally
  if(l==0) {E1=0; E2=0; C1={0,0}; C2={0,0};}
  else{ E1 = E1/(3.*l*(l+1.)); E2 = E2/(3.*l*(l+1.)); C1 = C1/(3.*l*(l+1.)); C2 = C2/(3.*l*(l+1.));}


  return make_tuple(Er,E1,E2,Cr,C1,C2);
}

VSH_MODES_REAL Triangulation::get_VectorSphericalHarmonicsModes(int Lmax_)
{
  size_t nModes=(Lmax_+1)*(Lmax_+1);

  std::vector<double> Er_lm, E1_lm, E2_lm;// MSQDiff_l;
  Er_lm.reserve(nModes); E1_lm.reserve(nModes); E2_lm.reserve(nModes);

  double Er,E1,E2;
  for(int l=0; l<=Lmax_; l++)
    for(int m=-l; m<=l; m++)
    {
      std::tie(Er, E1, E2) = get_Elm_coefficients(l, m);
      Er_lm.emplace_back(Er);
      E1_lm.emplace_back(E1);
      E2_lm.emplace_back(E2);
    }
  return std::make_tuple(Er_lm, E1_lm, E2_lm);
}

VSH_MODES_REAL_COMPLEX Triangulation::get_VectorSphericalHarmonicsModes_real_and_complex(int Lmax_)
{
  size_t nModes=(Lmax_+1)*(Lmax_+1);

  std::vector<double> Er_lm, E1_lm, E2_lm;// MSQDiff_l;
  COMP_VEC Cr_lm, C1_lm, C2_lm;
  Er_lm.reserve(nModes); E1_lm.reserve(nModes); E2_lm.reserve(nModes);
  Cr_lm.reserve(nModes); C1_lm.reserve(nModes); C2_lm.reserve(nModes);

  double Er,E1,E2;
  std::complex<double> Cr, C1, C2;
  for(int l=0; l<=Lmax_; l++)
    for(int m=-l; m<=l; m++)
    {
      std::tie(Er, E1, E2, Cr, C1, C2) = get_Elm_Clm_coefficients(l, m);
      Er_lm.emplace_back(Er);
      E1_lm.emplace_back(E1);
      E2_lm.emplace_back(E2);
      Cr_lm.emplace_back(Cr);
      C1_lm.emplace_back(C1);
      C2_lm.emplace_back(C2);
    }
  return std::make_tuple(Er_lm, E1_lm, E2_lm, Cr_lm, C1_lm, C2_lm);
}

std::vector<Vector3d> Triangulation::reconstruct_vec_field_from_modes(std::vector<double> Er, std::vector<double> E1, std::vector<double> E2, int Lmax_, bool add_rhat_)
{
  std::vector<Vector3d> reconst_field; reconst_field.resize(Vertices.size());
  if(add_rhat_)
  {
    for(size_t i=0; i<Vertices.size(); i++) reconst_field[i] = Vertices[i].r_hat;
  }
  else
  {
    for(size_t i=0; i<Vertices.size(); i++) reconst_field[i] = Vector3d(0,0,0);
  }

  //---
  // Vector3d PSIlm, PHIlm;
  // double Ylm;
  for(size_t i=0; i<Vertices.size(); i++)
  {
    unsigned mc=0;//to move along the given vectors
    for(int l=0; l<=Lmax_; l++)
      for(int m=-l; m<=l; m++)
      {
        reconst_field[i] += build_vector_from_one_mode(l, m, Er[mc], E1[mc], E2[mc], Vertices[i].cos_theta, (Vertices[i].pos_sph)[2], Vertices[i].r_hat, Vertices[i].theta_hat, Vertices[i].phi_hat);
        // Ylm = realYlm(l,m, Vertices[i].cos_theta, (Vertices[i].pos_sph)[2]);
        // std::tie(PSIlm, PHIlm) = PSIlmPHIlm(l, m, Vertices[i].cos_theta, (Vertices[i].pos_sph)[2], Vertices[i].theta_hat, Vertices[i].phi_hat);
        // reconst_field[i] += ( Vertices[i].r_hat * (Ylm*Er[mc]) + PSIlm*E1[mc] + PHIlm*E2[mc] );
        mc++;
      }
  }

  return reconst_field;
}

void Triangulation::update_all_triangles_integrated_curvature_tensor()
{
  for(size_t i=0; i<Edges.size(); i++)
  if(Edges[i].id<Edges[i].conj->id)
  {
    double ar = Edges[i].tri->area / ( Edges[i].tri->area + Edges[i].conj->tri->area );
    //...
    double cc =  (Edges[i].tri->unit_normal_vec).dot(Edges[i].conj->tri->unit_normal_vec);
    if(cc>1) cc=1.;
    double alpha = std::acos(cc);
    Vector3d convex = Edges[i].tri->com + Edges[i].tri->unit_normal_vec - ( Edges[i].conj->tri->com + Edges[i].conj->tri->unit_normal_vec );
    Vector3d concav = Edges[i].tri->com - Edges[i].tri->unit_normal_vec - ( Edges[i].conj->tri->com - Edges[i].conj->tri->unit_normal_vec );
    if(convex.norm()<concav.norm()) alpha *= -1;
    //...
    Vector3d n_a = Edges[i].tri->unit_normal_vec + Edges[i].conj->tri->unit_normal_vec;
    n_a = n_a/(n_a.norm()+NORM_EPS);
    //...
    Vector3d ee = Edges[i].l_vec/(Edges[i].l+NORM_EPS);
    Vector3d n_i = n_a.cross(ee);
    //... now C_int for each bond...
    for(size_t k=0; k<3; k++)
    for(size_t l=0; l<3; l++)
    {
      double t1 = ( 2. * ar * alpha + std::sin(alpha) + std::sin(alpha - 2. * ar * alpha) ) * n_a[k]*n_a[l];
      double t2 = ( 2. * ar * alpha - std::sin(alpha) - std::sin(alpha - 2. * ar * alpha) ) * n_i[k]*n_i[l];
      double t3 = ( 2. * std::cos(ar * alpha) * std::cos(alpha - ar * alpha) ) * ( n_a[k]*n_i[l] + n_a[l]*n_i[k] );
      Edges[i].C_int(k,l) += Edges[i].l/4. * (t1 + t2 + t3);
      Edges[i].conj->C_int(k,l) += Edges[i].l/4. * (t1 + t2 - t3);
    }
  }
  //..... now we compute the triangles C_int ....
  for(size_t i=0; i<Triangles.size(); i++)
  {
    Triangles[i].C_int += Triangles[i].e1->C_int;
    Triangles[i].C_int += Triangles[i].e2->C_int;
    Triangles[i].C_int += Triangles[i].e3->C_int;
  }
}


void Triangulation::update_all_Q_3D_xyplane()
{
  for(size_t i=0; i<Triangles.size(); i++)
  {
    double rot_ang_X, rot_ang_Y;
    std::tie(rot_ang_X, rot_ang_Y) = normal_ThetaX_ThetaY(Triangles[i].unit_normal_vec);
    //std::cout << rot_ang_X << " " << rot_ang_Y << std::endl;
    Triangles[i].Q_xyplane = calculate_triangle_elongation_tensor_xy_plane(rot_ang_X, rot_ang_Y, -(Triangles[i].e3->l_vec), Triangles[i].e1->l_vec);
    //Triangles[i].Q_xyplane = calculate_triangle_elongation_tensor_xy_plane(rot_ang_X, rot_ang_Y, Triangles[i].e1->l_vec, -(Triangles[i].e3->l_vec));
    // std::cout << "e1:\n" <<Triangles[i].e1->l_vec << "\ne3:\n" << -(Triangles[i].e3->l_vec) << std::endl;
    //std::cout << "Q_xyplane:\n" << Triangles[i].Q_xyplane << std::endl;
    //std::cout << Triangles[i].e1->l_vec << std::endl;
    //std::cout << Triangles[i].Q_xyplane << std::endl;
  }
}

Matrix3d Triangulation::compute_Q3d_Patch(unsigned patch_id_)
{
  double rot_ang_X, rot_ang_Y;
  Patches[patch_id_].Q_3D = Matrix3d::Zero();
  Patches[patch_id_].area = 0; Patches[patch_id_].unit_normal_vec = Vector3d(0,0,0);
  Patches[patch_id_].com = Vector3d(0,0,0);
  //Matrix3d Q_xy = Matrix3d::Zero();
  for(auto tt : Patches[patch_id_].triangles)
  {
    Patches[patch_id_].com += tt->com;
    Patches[patch_id_].area += tt->area;
    Patches[patch_id_].unit_normal_vec += tt->normal_vec; // for now, it has area in it
    std::tie(rot_ang_X, rot_ang_Y) = normal_ThetaX_ThetaY( tt->unit_normal_vec );
    Patches[patch_id_].Q_3D += tt->area * transform_matrices_from_xy_plane_to_triangle_plane(tt->Q_xyplane, rot_ang_X, rot_ang_Y);
    //Q_xy += tt->Q_xyplane * tt->area;
  }
  Patches[patch_id_].com /= Patches[patch_id_].triangles.size();
  //Q_xy /= Patches[patch_id_].area;
  Patches[patch_id_].unit_normal_vec /= ( (Patches[patch_id_].unit_normal_vec).norm()+NORM_EPS );
  //std::cout << "computed Q3d one::\n" << Patches[patch_id_].Q_3D << std::endl;
  Patches[patch_id_].Q_3D /= Patches[patch_id_].area;

  auto [qval, q_vec] = calculate_3x3_tensor_eigenvalues_and_eigenvector(Patches[patch_id_].Q_3D, 1);

  Patches[patch_id_].Q_eigvec = q_vec;
  Patches[patch_id_].Q_norm = qval;

  return Patches[patch_id_].Q_3D;
}

// Vector3d Triangulation::compute_Qeigvec_Patch(unsigned patch_id_)
// {
//   //Patches[patch_id_].Q_3D = Matrix3d::Zero();
//   Patches[patch_id_].area = 0; Patches[patch_id_].unit_normal_vec = Vector3d(0,0,0);
//   Patches[patch_id_].com = Vector3d(0,0,0);
//   Matrix3d Q_xy = Matrix3d::Zero();
//   for(auto tt : Patches[patch_id_].triangles)
//   {
//     Patches[patch_id_].com += tt->com;
//     Patches[patch_id_].area += tt->area;
//     Patches[patch_id_].unit_normal_vec += tt->normal_vec; // for now, it has area in it
//     Q_xy += tt->Q_xyplane * tt->area;
//   }
//   Patches[patch_id_].com /= Patches[patch_id_].triangles.size();
//   Q_xy /= Patches[patch_id_].area;
//   Patches[patch_id_].unit_normal_vec /= ( (Patches[patch_id_].unit_normal_vec).norm()+NORM_EPS );
//   double rot_ang_X, rot_ang_Y;
//   std::tie(rot_ang_X, rot_ang_Y) = normal_ThetaX_ThetaY( Patches[patch_id_].unit_normal_vec );
//
//   double eigval = std::sqrt( Q_xy(0,0)*Q_xy(0,0) + Q_xy(0,1)*Q_xy(0,1) );
//   double v0 = (eigval+Q_xy(0,0))/Q_xy(0,1);
//   Vector3d vxy = Vector3d(v0, 1, 0)/std::sqrt(2.);
//   Matrix3d Rx = R_x(rot_ang_X);
//   Matrix3d Ry = R_y(rot_ang_Y);
//   Matrix3d RxRy = Rx * Ry;
//   Patches[patch_id_].Q_eigvec = RxRy*vxy;
//   return Patches[patch_id_].Q_eigvec;
// }

void Triangulation::update_all_Patches_Q_3D()
{
  update_all_Q_3D_xyplane();
  for(size_t pid=0; pid<Patches.size(); pid++){compute_Q3d_Patch(pid);}
}

void Triangulation::update_all_elongation_tilt_angles_if_already_aligned()
{
  for(size_t pid=0; pid<Patches.size(); pid++)
  {
    auto [beta, tta, phi] = get_tilt_angle(Patches[pid].Q_eigvec, Patches[pid].com);
    Patches[pid].tilt_angle = beta;
    Patches[pid].rotated_theta = tta;
    Patches[pid].rotated_phi = phi;
  }
}

Matrix3d Triangulation::compute_Cint_Patch(unsigned patch_id_)
{
  Matrix3d C_int = Matrix3d::Zero();
  for(auto tt : Patches[patch_id_].triangles) C_int += tt->C_int;

  Patches[patch_id_].C_int = C_int;
  return Patches[patch_id_].C_int;
}

void Triangulation::update_all_Patches_Cint()
{
  update_all_triangles_integrated_curvature_tensor();

  for(size_t pid=0; pid<Patches.size(); pid++) compute_Cint_Patch(pid);
}
