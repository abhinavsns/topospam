#include <iostream>
#include <filesystem>
#include <vector>
// #include <algorithm>
// #include <execution>
#include "Tissue.hpp"

namespace fs = std::filesystem;

int main(int argc, char *argv[])
{
  int Nc = std::stoi(argv[1]);
  std::string realisation = std::string(argv[2]);
  std::string output_name = std::string(argv[3]);

  //std::string input_file = "/data/biophys/aamiri/stableFiles/inputInitConfig/random/"+std::string(argv[1])+"/random" + std::string(argv[1])+"v"+std::string(argv[2])+".txt";
  std::string input_file = "/data/biophys/aamiri/stableFiles/inputInitConfig/thomson/thomson" + std::string(argv[1])+".txt";
  ///data/biophys/aamiri/stableFiles/inputInitConfig/thomson/

  Tissue *tissue = new Tissue(Nc); // heap memory allocation

  tissue->build_2D_spherical_tissue_from_file(input_file.c_str());

  double Kc=1., A0c = 1., BondT=0.1, perim_elasticity = 0.;

  tissue->set_CellMechProps2D(Kc, A0c, BondT, perim_elasticity);

  tissue->rescale_spheroid_apical_size();

  tissue->update_apical_geometric_quantities();

  // auto func = [tissue](auto it){std::cout << it.id << " " << tissue->get_apical_bond_length(&it) << std::endl;};
  // std::for_each(
  //     std::execution::par_unseq,
  //     tissue->ListEdge.begin(),
  //     tissue->ListEdge.end(),
  //     func);

  const fs::path out_path{output_name};
  fs::create_directory(out_path);

  tissue->write_to_txt_files(output_name, 0);

  std::cout << "W0: " << tissue->W_2D() << std::endl;

  tissue->update_dW_a_2D();

  Vector3d net_froce = Vector3d(0.,0.,0.);
  for(std::vector<Vertex>::iterator it=tissue->ListVertex.begin();it!=tissue->ListVertex.end();it++) net_froce += it->dW_a;
  // net_froce /= tissue->ListVertex.size();

  std::cout << "net_froce:   " << net_froce[0] << " " << net_froce[1] << " " << net_froce[2] <<std::endl;

  tissue->write_to_json_file(output_name, 0);

  tissue->flip_short_bonds(0.04, 0.045, true);

  tissue->minimize_gsl();

  tissue->write_to_json_file(output_name, 1);

  double W1 = tissue->W_2D();

  net_froce = Vector3d(0.,0.,0.);
  for(std::vector<Vertex>::iterator it=tissue->ListVertex.begin();it!=tissue->ListVertex.end();it++) net_froce += it->dW_a;
  // net_froce /= tissue->ListVertex.size();

  std::cout << "net_froce2:   " << net_froce[0] << " " << net_froce[1] << " " << net_froce[2] <<std::endl;


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
  // for(int i=0; i<200000; i++)
  // {
  //   tissue->update_apical_geometric_quantities();
  //   tissue->update_dW_a_2D();
  //   tissue->evolve_vertex_positions(xi, dt, 0);
  //
  //   //double R_sp = std::sqrt ( tissue->ListCell.size()* tissue->CellMechProp.A0c/(4.*PI) );
  //   //tissue->set_apical_points_on_sphere_at_radius(R_sp);
  //
  //   for(std::vector<Edge>::iterator it=tissue->ListEdge.begin();it!=tissue->ListEdge.end();it++)
  //   {
  //     if (it->l_a< 0.04) {tissue->flip_edge(&(*it), 0.05);}
  //   }
  //
  // }
  //
  // tissue->update_apical_geometric_quantities();
  // double W2 = tissue->W();
  //
  // tissue->write_to_json_file(output_name, 2);
  //
  // tissue->minimize_gsl();
  //
  // tissue->write_to_json_file(output_name, 3);
  //
  // double W3 = tissue->W();
  //
  // std::cout << "W1: " << W1 << "  W2: " << W2 << "  W3:  " << W3 << std::endl;
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
