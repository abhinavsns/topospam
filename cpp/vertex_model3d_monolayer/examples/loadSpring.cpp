#include <iostream>
#include <vector>
#include "Tissue.hpp"

#define PI 3.14159265359

int main(int argc, char *argv[])
{

  Tissue *tissue = new Tissue(); // heap memory allocation

  int Nc = std::stoi(argv[1]);
  std::string realisation = std::string(argv[2]);
  double v_spring = std::stod(argv[3]);
  double k_spring = 0.01;
  double gamma = 0.;//std::stod(argv[4]);
  double lambda=0.;
  double diffD=0.;
  double xi=1., dt=0.0005;

  std::string input_file = "/data/biophys/aamiri/stableFiles/inputInitConfig/random/"+std::string(argv[1])+"/random" + std::string(argv[1])+"v"+std::string(argv[2])+".txt";

  tissue->build_2D_spherical_tissue_from_file(input_file.c_str());

  double Kc=1., A0c = 1., BondT=0.1, perim_elasticity = 0.;

  tissue->set_CellMechProps2D(Kc, A0c, BondT, perim_elasticity);

  tissue->rescale_spheroid_apical_size();

  tissue->update_apical_geometric_quantities();

  std::string output_name = "output/" ;

  bool spherical_constrain = true;

  //std::vector<double> cells_F_mag(tissue->ListCell.size(), 0.); // first zero traction
  for(int i=0; i<tissue->ListCell.size(); i++) tissue->ListCell[i].mag_F_a = 0.;
  tissue->update_active_forces(); // first pass mag_F_a=0 for passive evolution

  tissue->randomize_cell_planar_polarity();

  for(int i=0; i<100000; i++)//10x of this safer
  {
    tissue->update_apical_geometric_quantities();
    tissue->update_dW_a_2D();
    tissue->evolve_vertex_positions(xi, dt, spherical_constrain);
    tissue->rescale_spheroid_apical_size();

    for(std::vector<Edge>::iterator it=tissue->ListEdge.begin();it!=tissue->ListEdge.end();it++)
    {
      if (it->l_a< 0.04) {tissue->flip_edge(&(*it), 0.05);}
    }

  }
  tissue->write_to_txt_files(output_name, 0);

  // === now the active model; The spring extension ensemble ...
  double vt = 0; int vt_id = 0;
  while(vt<20)//10x of this safer
  {
    vt += v_spring * dt;
    tissue->update_apical_geometric_quantities();
    tissue->evolve_polarity(lambda, gamma, diffD, dt);
    for(int i=0; i<tissue->ListCell.size(); i++)
    {
      Vector3d cell_delta_cent = tissue->ListCell[i].cent_a - tissue->ListCell[i].cent_a_prev;
      double pu = tissue->ListCell[i].polarity.dot(cell_delta_cent);
      tissue->ListCell[i].mag_F_a = tissue->ListCell[i].mag_F_a - k_spring * (pu - v_spring * dt);
    }
    tissue->update_dW_a_2D();
    tissue->update_active_forces();
    tissue->evolve_vertex_positions(xi, dt, spherical_constrain);
    tissue->rescale_spheroid_apical_size();

    for(std::vector<Edge>::iterator it=tissue->ListEdge.begin();it!=tissue->ListEdge.end();it++)
    {
      if (it->l_a< 0.04) {tissue->flip_edge(&(*it), 0.05);}
    }

    if((vt_id+1)%10000==0)tissue->write_to_txt_files(output_name, int(vt_id/10000)+1);

    vt_id++;
  }



  delete tissue;
  return 0;
}
