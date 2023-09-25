#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>
#include "Tissue.hpp"

#include "../accessories/include/json.hpp"

namespace fs = std::filesystem;

#define PI 3.14159265359

int main(int argc, char *argv[])
{
  (void) argc;

  int Nc = std::stoi(argv[1]);
  std::string realisation = std::string(argv[2]);

  Tissue *tissue = new Tissue(Nc); // heap memory allocation

  double bond_T1_cutoff=0.04, opening_length=0.045;

  std::string input_file = "accessories/misc/input_files/random/"+std::to_string(Nc)+"/random" + std::string(argv[1])+"v"+std::string(argv[2])+".txt";
  //std::string input_file = "accessories/misc/input_files/thomson/thomson" + std::string(argv[1])+".txt";

  double Kc=1., V0c = 1., Ta=0.015, Tb=0.015, Tl=0.05;
  double lumen_P0=0, lumen_Kl=0.005, lumen_V0l=300.;

  double Heq = pow(2*V0c/sqrt(3), 1.0/3.0) * pow((Ta+Tb)/Tl, 2.0/3.0);
  double t1=9.0*Heq*Nc*V0c/PI - 3.0*pow(Heq, 4.0);
  double R_sph = sqrt(t1)/(6.0*Heq) - Heq/2.0;
  std::cout << "R_sph: " << R_sph << "  t1: " << t1 << std::endl;
  if(t1<=0 || R_sph<0.1) {R_sph = 4.0; Heq=0.25;}

  double Rb_Ra_ratio = 1. + Heq/R_sph;

  tissue->build_3D_spherical_tissue_from_file(input_file.c_str(), R_sph, Rb_Ra_ratio);

  tissue->set_CellMechProps3D(Kc, V0c, Ta, Tl, Tb);

  tissue->set_SpheroidMechProp(lumen_P0, lumen_Kl, lumen_V0l);

  std::string output_folder = "/data/biophys/aamiri/y2023/test/" ;
  //
  // tissue->rescale_spheroid_apical_size();

  // ==== write paramters -----
  const fs::path out_path{output_folder};
  fs::create_directory(out_path);
    using json = nlohmann::json;
    std::string json_file = output_folder + "parameters.json";
    std::ofstream json_out(json_file);
    json json_data;
    json_data["Nc"] = Nc; json_data["input_file"] = input_file;
    json_data["lumen_P0"] = lumen_P0; json_data["lumen_Kl"] = lumen_Kl;
    json_data["Kc"]=Kc; json_data["V0c"]=V0c; json_data["Ta"]=Ta; json_data["Tb"]=Tb; json_data["Tl"]=Tl;
    // json_data["output_folder"] = output_folder; json_data["xi"]=xi; json_data["dt"]=dt;
    // json_data["bond_T1_cutoff"] = bond_T1_cutoff; json_data["opening_length"] = opening_length;
    // json_data["Kc"]=Kc; json_data["A0c"]=A0c; json_data["BondT"]=BondT; json_data["perim_elasticity"]=perim_elasticity;
    json_out  << json_data << std::endl;
    json_out.close();
    //=----------------
  //

  Vector3d tissue_cent = Vector3d(0,0,0);
  for(size_t i=0; i<tissue->ListVertex.size(); i++) {tissue_cent += (tissue->ListVertex[i].pos_a + tissue->ListVertex[i].pos_b); }
  tissue_cent /= (2.*tissue->ListVertex.size());

  for(size_t i=0; i<tissue->ListVertex.size(); i++) {tissue->ListVertex[i].pos_a -= tissue_cent; tissue->ListVertex[i].pos_b -= tissue_cent; }
  //std::cout << "tissue_cent:   " << tissue_cent[0] << " " << tissue_cent[1] << " " << tissue_cent[2] << std::endl;

  tissue->update_3D_geometric_quantities();
  tissue->update_dW_3D_spherical();

  std::cout << "get_lumen_volume():  " << tissue->get_lumen_volume() << std::endl;


  // Vector3d ccent = tissue->get_cell_centroid(&(tissue->ListCell[0]));
  // std::cout << ccent[0] << " " << ccent[1] << " " << ccent[2] << " " << tissue->get_cell_volume(&tissue->ListCell[0]) << std::endl;
  //
  // ccent = tissue->get_cell_centroid(&(tissue->ListCell[1]));
  // std::cout << ccent[0] << " " << ccent[1] << " " << ccent[2] << " " << tissue->get_cell_volume(&tissue->ListCell[1]) << std::endl;

  // std::cout << tissue->ListCell[0].vol << " " << tissue->ListCell[1].vol << " " << tissue->ListCell[2].vol << std::endl;
  // std::cout << tissue->ListVertex[0].dW_a[0] << " " << tissue->ListVertex[0].dW_a[1] << " " << tissue->ListVertex[0].dW_a[2] << std::endl;
  //   std::cout << tissue->ListVertex[0].dW_b[0] << " " << tissue->ListVertex[0].dW_b[1] << " " << tissue->ListVertex[0].dW_b[2] << std::endl;
  // for(int i=0; i<tissue->ListCell.size(); i++) {std::cout << tissue->ListCell[i].vol << " "; }
  // std::cout << std::endl;

  //Vector3d net_grad(0,0,0);
  //for(size_t i=0; i<tissue->ListVertex.size(); i++) {net_grad += tissue->ListVertex[i].dW_a + tissue->ListVertex[i].dW_b;}
  //std::cout << "net_grad: " << net_grad[0] << " " << net_grad[1] << " " << net_grad[2] << std::endl;


  tissue->write_to_json_file(output_folder, 0);

  // for(std::vector<Cell>::iterator itc=tissue->ListCell.begin();itc!=tissue->ListCell.end();itc++)
  // {
  //  if(itc->id < 1) tissue->divide_cell_random_axis(&(*itc));
  // //   std::cout << itc->area_a << std::endl;
  // }
  double tt=0, dt=0.01;
  int fr=1;
  while(tissue->SphMechProp.V0l>200)
  {
    tt=0;
    while(tt<5.)
    {
      tt += dt;
      tissue->update_3D_geometric_quantities();
      tissue->update_dW_3D_spherical();
      for(size_t i=0; i<tissue->ListVertex.size(); i++)
      {
        tissue->ListVertex[i].pos_a -= tissue->ListVertex[i].dW_a * dt;
        tissue->ListVertex[i].pos_b -= tissue->ListVertex[i].dW_b * dt;
      }
    }
    //tissue->minimize_gsl3D();
    tissue->write_to_json_file(output_folder, fr);
    tissue->flip_short_bonds(bond_T1_cutoff, opening_length, true);

    tissue->SphMechProp.V0l -= 0.1;
    fr++;
  }

    std::cout << "get_lumen_volume():  " << tissue->get_lumen_volume() << std::endl;

    tissue->SphMechProp.Kl = 0;
    tissue->minimize_gsl3D();
    tissue->write_to_json_file(output_folder, 1001);


  //
  // net_grad = Vector3d(0,0,0);
  // for(int i=0; i<tissue->ListVertex.size(); i++) {net_grad += tissue->ListVertex[i].dW_a + tissue->ListVertex[i].dW_b;}
  // std::cout << "net_grad: " << net_grad[0] << " " << net_grad[1] << " " << net_grad[2] << std::endl;
  // //
  // // for(int i=0; i<tissue->ListCell.size(); i++) {std::cout << tissue->ListCell[i].vol << " "; }
  // // std::cout << std::endl;
  // //
  // // for(int i=0; i<tissue->ListEdge.size(); i++) std::cout << "area_lateral: " << tissue->ListEdge[i].area_lateral <<
  // // "  cent_l: " << tissue->ListEdge[i].cent_l[0] << " " << tissue->ListEdge[i].cent_l[1] << " " << tissue->ListEdge[i].cent_l[2] << std::endl;
  // //
  // //
  // std::cout << "  2.   get_lumen_volume():  " << tissue->get_lumen_volume() << std::endl;

  //tissue->minimize_gsl3D();
  //tissue->write_to_txt_files(output_name, 1);
  //
  // std::cout << "W0: " << tissue->W() << std::endl;

  // int ct=0;
  // for(std::vector<Edge>::iterator it=tissue->ListEdge.begin();it!=tissue->ListEdge.end();it++)
  // {
  //   if (it->id==1) {tissue->flip_edge(&(*it), 0.1);}
  //   // if ((it->head->pos_a - it->next->tail->pos_a).norm() >1e-10) {ct++; std::cout << "Warn!   " << (it->head->pos_a - it->next->tail->pos_a).norm() << std::endl;}
  //   // std::cout << "bond length: " << it->l_a << "  ids: " << it->id << ", " << it->next->id << std::endl;
  // }
  // // std::cout << "ct   " << ct << std::endl;
  // //
  // // std::cout << tissue->ListVertex.size() << " " << tissue->ListEdge.size() << " " << tissue->ListCell.size() << std::endl;
  //
  // tissue->write_to_txt_files(output_name, 1);
  //
  // for(std::vector<Cell>::iterator itc=tissue->ListCell.begin();itc!=tissue->ListCell.end();itc++)
  // {
  //  if(itc->id == 28) tissue->divide_cell_random_axis(&(*itc));
  // //   std::cout << itc->area_a << std::endl;
  // }

  // for(std::vector<Edge>::iterator it=tissue->ListEdge.begin();it!=tissue->ListEdge.end();it++)
  // {
  //   if (it->id == 1) {tissue->flip_edge(&(*it), 0.1); break;}
  // }

  // tissue->write_to_txt_files(output_name, 1);
  // double xi=1., dt=1e-4;
  // for(int i=0; i<10000; i++)
  // {
  //   tissue->update_apical_geometric_quantities();
  //   tissue->update_dW_a_2D();
  //   tissue->evolve_vertex_positions(xi, dt);
  //
  //   // double R_sp = std::sqrt ( tissue->ListCell.size()* tissue->CellMechProp.A0c/(4.*PI) );
  //   // tissue->set_apical_points_on_sphere_at_radius(R_sp);
  //
  //   for(std::vector<Edge>::iterator it=tissue->ListEdge.begin();it!=tissue->ListEdge.end();it++)
  //   {
  //     if (it->l_a< 0.04) {tissue->flip_edge(&(*it), 0.05);}
  //   }
  //
  // }
  //
  // tissue->update_apical_geometric_quantities();
  // std::cout << "W1: " << tissue->W() << std::endl;
  //
  // tissue->write_to_txt_files(output_name, 2);
  //
  // //tissue->minimize_gsl();
  // tissue->my_cg_minimizer();
  // // //tissue->minimize_gsl2();
  // tissue->update_apical_geometric_quantities();
  // tissue->write_to_txt_files(output_name, 3);
  // std::cout << "W2: " << tissue->W() << std::endl;

  delete tissue;
  return 0;
}
