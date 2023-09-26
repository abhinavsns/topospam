#include "Writer.hpp"

//Vertex constructor
// Writer::Writer(){}

Writer::Writer(unsigned frame_, std::string output_path_) : frame(frame_), output_path(output_path_){
  //create output_path id not there ...
  if (!output_path.empty() && output_path.back() != '/')
    output_path += '/';

  //apical_tri_surface = std::unique_ptr<Triangulation>();
  //apical_tri_surface = Triangulation();
  //
  // namespace fs = std::filesystem;
  // const fs::path out_path{output_path};
  // fs::create_directory(out_path);
}

void Writer::rotate_tissue_to_align_z_with_omega(Vector3d omega_)
{
  auto [rot_ang_X, rot_ang_Y] = normal_ThetaX_ThetaY(omega_/omega_.norm());
  auto inv_RyRx = get_inverse_Ry_Rx(rot_ang_X, rot_ang_Y);

  for(size_t i=0; i<ListVertex.size(); i++)
  {
    ListVertex[i].pos_a = inv_RyRx * ListVertex[i].pos_a;
    ListVertex[i].pos_b = inv_RyRx * ListVertex[i].pos_b;
    ListVertex[i].pos_a_prev = inv_RyRx * ListVertex[i].pos_a_prev;
  }
  //...
  for(size_t i=0; i<ListCell.size(); i++)
  {
    ListCell[i].polarity = inv_RyRx * ListCell[i].polarity;
    ListCell[i].avg_polarity_j = inv_RyRx * ListCell[i].avg_polarity_j;
    ListCell[i].cent_a = inv_RyRx * ListCell[i].cent_a;
  }
}

void Writer::rotate_tissue_around_z_axis(double angle_)
{
  auto Rz = R_z(angle_);
  for(size_t i=0; i<ListVertex.size(); i++)
  {
    ListVertex[i].pos_a = Rz * ListVertex[i].pos_a;
    ListVertex[i].pos_b = Rz * ListVertex[i].pos_b;
    ListVertex[i].pos_a_prev = Rz * ListVertex[i].pos_a_prev;
  }
  //...
  for(size_t i=0; i<ListCell.size(); i++)
  {
    ListCell[i].polarity = Rz * ListCell[i].polarity;
    ListCell[i].avg_polarity_j = Rz * ListCell[i].avg_polarity_j;
    ListCell[i].cent_a = Rz * ListCell[i].cent_a;
  }
}

void Writer::build_apical_tissue_triangulation(bool update_apical_geometry)
{
  if (update_apical_geometry) update_apical_geometric_quantities();
  //go through edges and build triangles with basic info..
  apical_tri_surface.Triangles.reserve(ListEdge.size());
  apical_tri_surface.Vertices.reserve(ListVertex.size()+ListCell.size());
  apical_tri_surface.Edges.reserve(3*ListEdge.size());
  apical_tri_surface.Patches.reserve(ListCell.size());

  // constructing enough Patches -- used for coarse graining
  // I first contruct enough with a unique id, and then fill them up with their triangulation list
  for(size_t i=0; i<ListCell.size(); i++) {ListCell[i].id =i; apical_tri_surface.Patches.emplace_back(Patch(i));}

  int tid=0;
  for(auto te : ListEdge)
  {
    TriVertex *v1 = apical_tri_surface.append_vertex(te.cell->cent_a);
    TriVertex *v2 = apical_tri_surface.append_vertex(te.tail->pos_a);
    TriVertex *v3 = apical_tri_surface.append_vertex(te.head->pos_a);
    //....
    TriEdge *e1 = apical_tri_surface.append_edge(v1, v2);
    TriEdge *e2 = apical_tri_surface.append_edge(v2, v3);
    TriEdge *e3 = apical_tri_surface.append_edge(v3, v1);
    //...
    apical_tri_surface.Triangles.emplace_back(v1, v2, v3, e1, e2, e3, tid);
    Triangle *tri_ = &(apical_tri_surface.Triangles.back());
    e1->tri = tri_; e2->tri = tri_; e3->tri = tri_;
    tri_->polarity = te.cell->polarity;
    apical_tri_surface.Patches[te.cell->id].triangles.emplace_back(tri_);
    tid++;
  }

  //  at this point only half of conj bonds are set. let set the other half
  for(std::vector<TriEdge>::iterator ite=apical_tri_surface.Edges.begin();ite!=apical_tri_surface.Edges.end();ite++)
  {
    if(ite->conj != nullptr) ite->conj->conj = &(*ite);
  }
}

void Writer::write_apical_surface_vector_spherical_harmonics_modes_to_json(int Lmax_)
{
  apical_tri_surface.update_all_geometries();
  std::vector<double> Er, E1, E2;
  COMP_VEC Cr, C1, C2;
  //std::tie(Er, E1, E2) = apical_tri_surface.get_VectorSphericalHarmonicsModes(Lmax_);
  std::tie(Er, E1, E2, Cr, C1, C2) = apical_tri_surface.get_VectorSphericalHarmonicsModes_real_and_complex(Lmax_);
  size_t nModes = Er.size();
  std::vector<double> Cr_real(nModes), Cr_imag(nModes), C1_real(nModes), C1_imag(nModes), C2_real(nModes), C2_imag(nModes);
  for(size_t i=0; i<nModes; i++)
  {
    Cr_real[i]=Cr[i].real(); Cr_imag[i] = Cr[i].imag();
    C1_real[i]=C1[i].real(); C1_imag[i] = C1[i].imag();
    C2_real[i]=C2[i].real(); C2_imag[i] = C2[i].imag();
  }
  //----
  std::string json_file = output_path + "VSH_modes_apicalsurface_fr"+std::to_string(frame)+".json";
  std::ofstream json_out(json_file);
  json json_data;
  json_data["Er"] = Er; json_data["E1"] = E1; json_data["E2"] = E2;
  json_data["Cr_realpart"] = Cr_real; json_data["Cr_imagpart"] = Cr_imag;
  json_data["C1_realpart"] = C1_real; json_data["C1_imagpart"] = C1_imag;
  json_data["C2_realpart"] = C2_real; json_data["C2_imagpart"] = C2_imag;
  json_data["Lmax"] = Lmax_;
  json_out  << json_data << std::endl;
  json_out.close();
  //=----------------
}

void Writer::write_apical_polarity_vector_spherical_harmonics_modes_to_json(int Lmax_)
{
  std::vector<Vector3d> pol; pol.resize(apical_tri_surface.Vertices.size());

  for(size_t i=0; i<apical_tri_surface.Vertices.size(); i++)
  {
    Vector3d avg_polarity = Vector3d(0,0,0);
    for(auto tt : apical_tri_surface.Vertices[i].connected_triangles) avg_polarity += tt->polarity;

    avg_polarity /= apical_tri_surface.Vertices[i].connected_triangles.size();

    apical_tri_surface.Vertices[i].vec_field = project_on_tangent_plane(apical_tri_surface.Vertices[i].r_hat, avg_polarity);
    pol[i] = apical_tri_surface.Vertices[i].vec_field;
  }

  apical_tri_surface.update_all_geometries(pol);
  std::vector<double> Er, E1, E2;
  COMP_VEC Cr, C1, C2;
  //std::tie(Er, E1, E2) = apical_tri_surface.get_VectorSphericalHarmonicsModes(Lmax_);
  std::tie(Er, E1, E2, Cr, C1, C2) = apical_tri_surface.get_VectorSphericalHarmonicsModes_real_and_complex(Lmax_);
  size_t nModes = Er.size();
  std::vector<double> Cr_real(nModes), Cr_imag(nModes), C1_real(nModes), C1_imag(nModes), C2_real(nModes), C2_imag(nModes);
  for(size_t i=0; i<nModes; i++)
  {
    Cr_real[i]=Cr[i].real(); Cr_imag[i] = Cr[i].imag();
    C1_real[i]=C1[i].real(); C1_imag[i] = C1[i].imag();
    C2_real[i]=C2[i].real(); C2_imag[i] = C2[i].imag();
  }
  //----
  std::string json_file = output_path + "VSH_modes_apicalpolarity_fr"+std::to_string(frame)+".json";
  std::ofstream json_out(json_file);
  json json_data;
  json_data["Er"] = Er; json_data["E1"] = E1; json_data["E2"] = E2;
  json_data["Cr_realpart"] = Cr_real; json_data["Cr_imagpart"] = Cr_imag;
  json_data["C1_realpart"] = C1_real; json_data["C1_imagpart"] = C1_imag;
  json_data["C2_realpart"] = C2_real; json_data["C2_imagpart"] = C2_imag;
  json_data["Lmax"] = Lmax_;
  json_out  << json_data << std::endl;
  json_out.close();
  //=----------------

}

void Writer::write_apical_or_basal_polygonal_surface_to_vtk(std::string surf_type)
{
  std::string write_name = output_path + surf_type + "poly_" + std::to_string(frame) + ".vtk";

  std::ofstream write_out; write_out.open(write_name);
  if (write_out.is_open())
  {
    write_out << "# vtk DataFile Version 1.0\n2D Unstructured Grid of Polygons\nASCII\n\nDATASET POLYDATA\n";
    write_out << "POINTS " + std::to_string(ListVertex.size()) + " float\n" ;
    if(surf_type=="apical") for(auto v : ListVertex){write_out << v.pos_a[0] << " " << v.pos_a[1] << " " << v.pos_a[2] << std::endl;}
    else if(surf_type=="basal") for(auto v : ListVertex){write_out << v.pos_b[0] << " " << v.pos_b[1] << " " << v.pos_b[2] << std::endl;}
    else {std::cout << " ERROR:  couldn't write the polygonal mesh to vtk. the surface type not recognized!" << std::endl;}
    int n_tot_elements =0;
    for(auto c : ListCell) n_tot_elements += c.Vlist.size()+1;
    write_out << "\nPOLYGONS "+ std::to_string(ListCell.size()) + " " + std::to_string(n_tot_elements);
    for(auto c : ListCell)
    {
      write_out << "\n" << c.Vlist.size();//+std::to_string(c.Vlist.size());
      for(auto v : c.Vlist) write_out << " " << v->id;// + std::to_string(v->id) ;
    }
    write_out << "\n\nCELL_DATA "+ std::to_string(ListCell.size())+"\n";
    write_out << "SCALARS neighbors float 1\n";
    write_out << "LOOKUP_TABLE default\n";
    for(auto c : ListCell) write_out << std::to_string(c.Vlist.size()) + " ";
    write_out << std::endl;
    //........
    write_out << "SCALARS beta float 1\n";
    write_out << "LOOKUP_TABLE default\n";
    for(size_t i=0; i<apical_tri_surface.Patches.size(); i++) write_out << apical_tri_surface.Patches[i].tilt_angle * (180/M_PI) << " ";
    write_out << std::endl;
    //........
    write_out << "SCALARS theta float 1\n";
    write_out << "LOOKUP_TABLE default\n";
    for(size_t i=0; i<apical_tri_surface.Patches.size(); i++) write_out << apical_tri_surface.Patches[i].rotated_theta * (180/M_PI) << " ";
    write_out << std::endl;
  }
  write_out.close();
}

void Writer::write_apical_or_basal_triangulated_surface_to_vtk(std::string surf_type)
{
  std::string write_name = output_path + surf_type + "tri_" + std::to_string(frame) + ".vtk";
  std::ofstream write_out; write_out.open(write_name);
  if (write_out.is_open())
  {
    write_out << "# vtk DataFile Version 1.0\n2D Unstructured Grid of Linear Triangles\nASCII\n\nDATASET POLYDATA\n";
    write_out << "POINTS " + std::to_string(ListVertex.size()+ListCell.size()) + " float\n" ;
    if(surf_type=="apical")
    {
      update_apical_geometric_quantities();
      for(auto v : ListVertex) {write_out << v.pos_a[0] << " " << v.pos_a[1] << " " << v.pos_a[2] << std::endl;}
      for(auto c : ListCell) {write_out << c.cent_a[0] << " " << c.cent_a[1] << " " << c.cent_a[2]<<std::endl;}
    }
    else if(surf_type=="basal")
    {
      update_basal_geometric_quantities();
      for(auto v : ListVertex) {write_out << v.pos_b[0] << " " << v.pos_b[1] << " " << v.pos_b[2] << std::endl;}
      for(auto c : ListCell) {write_out << c.cent_b[0] << " " << c.cent_b[1] << " " << c.cent_b[2]<<std::endl;}
    } else {std::cout << " ERROR:  couldn't write the polygonal mesh to vtk. the surface type not recognized!" << std::endl;}

    size_t n_triangle=0;
    for(auto c : ListCell) n_triangle += c.Vlist.size();
    std::vector<std::tuple<int,int,int>> triangles; triangles.reserve(n_triangle);
    for(auto c : ListCell){
      for(auto e : c.Elist) triangles.emplace_back(std::make_tuple(ListVertex.size()+c.id, e->tail->id, e->head->id));
    };
    if(triangles.size() != n_triangle) std::cout << "number of triangles did not match!!!" << n_triangle <<" "<< triangles.size()<< std::endl;

    write_out << "\nPOLYGONS "+ std::to_string(n_triangle) + " " + std::to_string(4*n_triangle) + "\n";
    for(auto tri : triangles)
    {
      int a, b, c;
      std::tie(a,b,c) = tri;
      write_out << "3 " << a << " " << b << " " << c << std::endl;
    }
    write_out << "\n\nCELL_DATA "+ std::to_string(n_triangle)+"\n";
    write_out << "SCALARS neighbors float 1\n";
    write_out << "LOOKUP_TABLE default\n";
    for(size_t i_tri=0; i_tri<n_triangle; i_tri++) write_out << 3 << " ";
    write_out << std::endl;
  }
  write_out.close();
}

void Writer::write_full_3D_triangulation_to_vtk()
{
  std::string write_name = output_path + "3Dtri_" + std::to_string(frame) + ".vtk";
  std::ofstream write_out; write_out.open(write_name);
  update_apical_geometric_quantities();
  update_basal_geometric_quantities();
  update_lateral_geometric_quantities();
  int Nv = ListVertex.size();
  int Nc = ListCell.size();
  int Ne = ListEdge.size();
  //std::cout << Nv << " " << Nc << " " << Ne << " " << ListEdge.size() << std::endl;
  size_t n_triangle= 6 * Ne;
  //for(auto c : ListCell) n_triangle += 2 * c.Vlist.size();
  std::vector<std::tuple<int,int,int>> triangles; triangles.reserve(n_triangle);
  std::vector<int> p_class; p_class.reserve(n_triangle);//lateral=4, apical=basal=c.Vlist.size()
  for(int i=0; i<Nv; i++) ListVertex[i].id=i;
  for(int i=0; i<Nc; i++) ListCell[i].id=i;
  //int e_counter=0;
  for(size_t i=0; i<ListEdge.size(); i++) {
    ListEdge[i].id=i;
    //if(ListEdge[i].id<ListEdge[i].conj->id) {ListEdge[i].id=e_counter; ListEdge[i].conj->id=Ne+e_counter; e_counter++;}
  }
  if (write_out.is_open())
  {
    write_out << "# vtk DataFile Version 1.0\n2D Unstructured Grid of Linear Triangles\nASCII\n\nDATASET POLYDATA\n";
    //write_out << "POINTS " + std::to_string(2*Nv+2*Nc) + " float\n" ;
    write_out << "POINTS " + std::to_string(2*Nv+2*Nc+Ne) + " float\n" ;
    for(auto v : ListVertex) {
      write_out << v.pos_a[0] << " " << v.pos_a[1] << " " << v.pos_a[2] << std::endl;
      write_out << v.pos_b[0] << " " << v.pos_b[1] << " " << v.pos_b[2] << std::endl;
    }
    for(auto c : ListCell) {
      write_out << c.cent_a[0] << " " << c.cent_a[1] << " " << c.cent_a[2]<<std::endl;
      write_out << c.cent_b[0] << " " << c.cent_b[1] << " " << c.cent_b[2]<<std::endl;
    }
    for(auto e : ListEdge) /*{if(e.id < e.conj->id)*/{write_out << e.cent_l[0] << " " << e.cent_l[1] << " " << e.cent_l[2] << std::endl;}// avoid repeating conj lateral triangles

    for(auto c : ListCell){ // fill in the triangle coordinates
      int apical_center = 2*Nv+2*c.id;
      int basal_center = 2*Nv+2*c.id+1;
      for(auto e : c.Elist)
      {
        //apical_triangle
        triangles.emplace_back(std::make_tuple(apical_center, 2*e->tail->id, 2*e->head->id));
        p_class.emplace_back(c.Vlist.size());
        //basal triangle
        triangles.emplace_back(std::make_tuple(basal_center, 2*e->tail->id+1, 2*e->head->id+1));
        p_class.emplace_back(c.Vlist.size());
        //4 lateral triangles
        //if(e->id < e->conj->id) // avoid repeating lateral triangles
        //{
          int e_lat_cent = 2*Nv + 2*Nc + e->id;
  //  std::cout << "*******  e_id:    " << e->id << std::endl;
          triangles.emplace_back(std::make_tuple(e_lat_cent, 2*e->tail->id, 2*e->head->id)); p_class.emplace_back(4);
          triangles.emplace_back(std::make_tuple(e_lat_cent, 2*e->tail->id+1, 2*e->head->id+1)); p_class.emplace_back(4);
          triangles.emplace_back(std::make_tuple(e_lat_cent, 2*e->tail->id, 2*e->tail->id+1)); p_class.emplace_back(4);
          triangles.emplace_back(std::make_tuple(e_lat_cent, 2*e->head->id, 2*e->head->id+1)); p_class.emplace_back(4);
        //}
      }
    };
    if(triangles.size() != n_triangle) std::cout << "number of triangles did not match!!!" << n_triangle <<" "<< triangles.size()<< std::endl;

    write_out << "\nPOLYGONS "+ std::to_string(n_triangle) + " " + std::to_string(4*n_triangle) + "\n";
    for(auto tri : triangles)
    {
      int a, b, c;
      std::tie(a,b,c) = tri;
      write_out << "3 " << a << " " << b << " " << c << std::endl;
    }
    write_out << "\nCELL_DATA "+ std::to_string(n_triangle)+"\n";
    write_out << "SCALARS neighbors float 1\n";
    write_out << "LOOKUP_TABLE default\n";
    for(size_t i_tri=0; i_tri<n_triangle; i_tri++) write_out << p_class[i_tri] << " ";
    write_out << std::endl;
  }
  write_out.close();
}

void Writer::write_polarity_field_to_vtk()
{
  update_apical_geometric_quantities();
  std::string write_name = output_path + "polarity_" + std::to_string(frame) + ".vtk";
  std::ofstream write_out; write_out.open(write_name);
  if (write_out.is_open())
  {
    write_out << "# vtk DataFile Version 1.0\npolarity field\nASCII\n\nDATASET POLYDATA\n";
    write_out << "POINTS " + std::to_string(ListCell.size()) + " float\n" ;
    for(auto c : ListCell) {write_out << c.cent_a[0] << " " << c.cent_a[1] << " " << c.cent_a[2]<<std::endl;}

    write_out << "\nPOINT_DATA "+ std::to_string(ListCell.size()) + "\n";
    write_out << "SCALARS F_mag float\nLOOKUP_TABLE default\n";
    for(auto c : ListCell) {write_out << c.mag_F_a << " ";}

    write_out << "\nVECTORS polarity float\n";
    for(auto c : ListCell) {write_out << c.polarity[0] << " " << c.polarity[1] << " " << c.polarity[2]<<std::endl;}
  }
  write_out.close();
}

void Writer::write_nematicity_field_to_vtk()
{
  size_t n_patches = apical_tri_surface.Patches.size();
  std::vector<Vector3d> patches_com; patches_com.reserve(n_patches);
  std::vector<Vector3d> eigenvecs; eigenvecs.reserve(n_patches);
  // std::vector<Vector3d> eigenvecs_min; eigenvecs_min.reserve(n_patches);
  std::vector<double> eigenvals; eigenvals.reserve(n_patches);
  for(size_t pid=0; pid<n_patches; pid++)
  {
    Vector3d pcom = apical_tri_surface.Patches[pid].com * 1.01;
    //Matrix3d pQ3d = apical_tri_surface.Patches[pid].Q_3D;
    //auto [qval, q_vec] = calculate_3x3_tensor_eigenvalues_and_eigenvector(pQ3d, 1);
    patches_com.emplace_back(pcom);
    eigenvecs.emplace_back(apical_tri_surface.Patches[pid].Q_eigvec);
    eigenvals.emplace_back(apical_tri_surface.Patches[pid].Q_norm);
  }

  std::string write_name = output_path + "elongation_" + std::to_string(frame) + ".vtk";
  std::ofstream write_out; write_out.open(write_name);
  if (write_out.is_open())
  {
    write_out << "# vtk DataFile Version 1.0\nelongation field\nASCII\n\nDATASET POLYDATA\n";
    write_out << "POINTS " << n_patches << " float\n" ;
    for(size_t pid=0; pid<n_patches; pid++) {write_out << patches_com[pid][0] << " " << patches_com[pid][1] << " " << patches_com[pid][2]<<std::endl;}

    write_out << "\nPOINT_DATA " << n_patches << "\n";
    write_out << "SCALARS q float\nLOOKUP_TABLE default\n";
    for(size_t pid=0; pid<n_patches; pid++) {write_out << eigenvals[pid] << " ";}

    write_out << "\nVECTORS elongation float\n";
    for(size_t pid=0; pid<n_patches; pid++) {write_out << eigenvecs[pid][0] << " " << eigenvecs[pid][1] << " " << eigenvecs[pid][2]<<std::endl;}
  }
  write_out.close();
}

void Writer::write_data_to_vtk()
{
  if(write_apical_polygonal_surface) write_apical_or_basal_polygonal_surface_to_vtk("apical");
  else if(write_basal_polygonal_surface) write_apical_or_basal_polygonal_surface_to_vtk("basal");

  if(write_apical_triangulated_surface) write_apical_or_basal_triangulated_surface_to_vtk("apical");
  else if(write_basal_triangulated_surface) write_apical_or_basal_triangulated_surface_to_vtk("basal");

  if(write_full_3D_triangulated_surface) write_full_3D_triangulation_to_vtk();

  if(write_polarity_field) write_polarity_field_to_vtk();

  if(write_nematic_field) write_nematicity_field_to_vtk();
}

void Writer::write_Patches_data()
{
  apical_tri_surface.update_all_Patches_Q_3D();
  apical_tri_surface.update_all_elongation_tilt_angles_if_already_aligned();
  apical_tri_surface.update_all_Patches_Cint();

  size_t n_patches = apical_tri_surface.Patches.size();
  std::vector<int> patch_ids; patch_ids.reserve(n_patches);
  std::vector<double> patches_area; patches_area.reserve(n_patches);
  std::vector<double> patches_beta; patches_beta.reserve(n_patches);
  std::vector<double> patches_rotated_theta; patches_rotated_theta.reserve(n_patches);
  std::vector<double> patches_rotated_phi; patches_rotated_phi.reserve(n_patches);
  std::vector<std::vector<double>> patches_com; patches_com.reserve(n_patches);
  std::vector<std::vector<double>> patches_Q3d; patches_Q3d.reserve(n_patches);
  std::vector<std::vector<double>> patches_Cint; patches_Cint.reserve(n_patches);
  for(size_t pid=0; pid<n_patches; pid++)
  {
    //std::cout << "num of tri: " <<  apical_tri_surface.Patches[pid].triangles.size() << std::endl;
    patches_area.emplace_back(apical_tri_surface.Patches[pid].area);
    patches_beta.emplace_back(apical_tri_surface.Patches[pid].tilt_angle);
    patches_rotated_theta.emplace_back(apical_tri_surface.Patches[pid].rotated_theta);
    patches_rotated_phi.emplace_back(apical_tri_surface.Patches[pid].rotated_phi);
    Vector3d pcom = apical_tri_surface.Patches[pid].com;
    Matrix3d pQ3d = apical_tri_surface.Patches[pid].Q_3D;
    Matrix3d Cint = apical_tri_surface.Patches[pid].C_int;
    patch_ids.emplace_back( pid );
    std::vector<double> com_vec = {pcom[0], pcom[1], pcom[2]};
    std::vector<double> Q3d_vec = {pQ3d(0,0), pQ3d(0,1), pQ3d(0,2), pQ3d(1,0), pQ3d(1,1), pQ3d(1,2), pQ3d(2,0), pQ3d(2,1), pQ3d(2,2)};
    std::vector<double> Cint_vec = {Cint(0,0), Cint(0,1), Cint(0,2), Cint(1,0), Cint(1,1), Cint(1,2), Cint(2,0), Cint(2,1), Cint(2,2)};
    patches_com.emplace_back(com_vec);
    patches_Q3d.emplace_back(Q3d_vec);
    patches_Cint.emplace_back(Cint_vec);
  }

  std::string json_file = output_path + "Patches_fr"+std::to_string(frame)+".json";
  std::ofstream json_out(json_file);
  json json_data;
  json_data["patch_ids"] = patch_ids; json_data["patches_com"] = patches_com;
  json_data["patches_Q3d"] = patches_Q3d; json_data["patches_Cint"] = patches_Cint;
  json_data["patches_area"] = patches_area;
  json_data["patches_beta"] = patches_beta;
  json_data["patches_rotated_theta"] = patches_rotated_theta;
  json_data["patches_rotated_phi"] = patches_rotated_phi;
  json_out  << json_data << std::endl;
  json_out.close();
}

//Write destructor
Writer::~Writer(){
}
