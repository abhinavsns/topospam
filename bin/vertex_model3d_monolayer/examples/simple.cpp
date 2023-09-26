#include <iostream>
#include <vector>
#include "Tissue.hpp"

int main(int argc, char *argv[])
{

  Tissue *tissue = new Tissue(); // heap memory allocation

  int Nc = std::stoi(argv[1]);
  std::string realisation = std::string(argv[2]);

  std::string input_file = "/data/biophys/aamiri/stableFiles/inputInitConfig/random/"+std::string(argv[1])+"/random" + std::string(argv[1])+"v"+std::string(argv[2])+".txt";

  tissue->build_2D_spherical_tissue_from_file(input_file.c_str());

  double Kc=1., A0c = 1., BondT=0.1;

  tissue->set_CellMechProps2D(Kc, A0c, BondT);

  tissue->rescale_spheroid_apical_size();

  tissue->update_apical_geometric_quantities();

  std::string output_name = "output/" ;

  tissue->write_to_txt_files(output_name, 0);

  int ct=0;
  for(std::vector<Edge>::iterator it=tissue->ListEdge.begin();it!=tissue->ListEdge.end();it++)
  {
    if (it->id==1) {tissue->flip_edge(&(*it), 0.1);}
    // if ((it->head->pos_a - it->next->tail->pos_a).norm() >1e-10) {ct++; std::cout << "Warn!   " << (it->head->pos_a - it->next->tail->pos_a).norm() << std::endl;}
    // std::cout << "bond length: " << it->l_a << "  ids: " << it->id << ", " << it->next->id << std::endl;
  }
  // std::cout << "ct   " << ct << std::endl;
  //
  // std::cout << tissue->ListVertex.size() << " " << tissue->ListEdge.size() << " " << tissue->ListCell.size() << std::endl;

  tissue->write_to_txt_files(output_name, 1);

  for(std::vector<Cell>::iterator itc=tissue->ListCell.begin();itc!=tissue->ListCell.end();itc++)
  {
   if(itc->id == 28) tissue->divide_cell_random_axis(&(*itc));
  //   std::cout << itc->area_a << std::endl;
  }

  tissue->write_to_txt_files(output_name, 2);

  delete tissue;
  return 0;
}
