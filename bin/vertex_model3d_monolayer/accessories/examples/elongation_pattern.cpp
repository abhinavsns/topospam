#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include "../include/json.hpp"
using json = nlohmann::json;

//#include "PostProcess-inline.h"
#include "Common-inline.h"
#include "Bond-inline.h"


#include "Tissue.hpp"

namespace fs = std::filesystem;

int main(int argc, char *argv[])
{
  (void) argc;

  std::string data_folder = std::string(argv[1]);
  int fr_begin = std::stoi(argv[2]);
  int fr_end = std::stoi(argv[3]);
  //std::string out_folder = std::string(argv[4]);

  // ==== make post_process folder -----
  std::string analysis_folder = data_folder+"analysis";
  const fs::path out_path{analysis_folder};
  fs::create_directory(out_path);

  for(int fr=fr_begin; fr<=fr_end; fr++)
  {
    std::string json_param_file = data_folder + "parameters.json";
    std::ifstream parfile(json_param_file);
    json param_data = json::parse(parfile);
    double dt = param_data["dt"];


    std::string json_file = data_folder + "frame_" + std::to_string(fr) + ".json";
    std::ifstream f(json_file);
    json json_data = json::parse(f);

    int Nv=json_data["Num_vertices"];
    // int Ne=json_data["Num_edges"];
    int Nc=json_data["Num_cells"];

    auto vertex_coord_a = json_data["vertex_pos_a"];
    auto vertex_prev_coord_a = json_data["vertex_prev_pos_a"];
    auto cells_Q_tensor = json_data["cell_Q3d_tensor"];

    auto cells_polarity_js = json_data["cells_polarity"];

    //=========================
    std::vector<Vector3d> cell_centers; cell_centers.resize(Nc);
    Tissue *tissue = new Tissue(Nc); // heap memory allocation

    tissue->build_3D_spherical_tissue_from_json_file(json_file);
    tissue->update_apical_geometric_quantities();
    for(int i=0; i<Nc; i++) cell_centers[i]=tissue->ListCell[i].cent_a;
    delete tissue;
    //=========================

    std::vector<Vector3d> cells_polarity; cells_polarity.resize(Nc);
    for(int i=0; i<Nc; i++) {cells_polarity[i] = Vector3d(cells_polarity_js[i][0], cells_polarity_js[i][1], cells_polarity_js[i][2]);}

    std::vector<Vector3d> pos_t1(Nv), pos_t2(Nv);
    for(int i=0; i<Nv; i++)
    {
      pos_t1[i] = Vector3d(vertex_prev_coord_a[i][0], vertex_prev_coord_a[i][1], vertex_prev_coord_a[i][2]);
      pos_t2[i] = Vector3d(vertex_coord_a[i][0], vertex_coord_a[i][1], vertex_coord_a[i][2]);
    }
    // //std::pair<Vector3d, Matrix3d>
    auto [Gamma, MM] = get_ang_momentum_inertia_tensor(pos_t1, pos_t2, dt);

    auto Omega = MM.inverse() * Gamma;

    auto [rot_ang_X, rot_ang_Y] = normal_ThetaX_ThetaY(Omega/Omega.norm());

    auto inv_RyRx = get_inverse_Ry_Rx(rot_ang_X, rot_ang_Y);
    //auto RyRx = inv_RyRx.inverse();

    //auto Omega_rotated = inv_RyRx*Omega; //sanity check; it must be along z axis

    //std::cout << "Omega rotated: \n"<< Omega_rotated << "\nRatio:\n"<< Omega.norm()/Omega_rotated.norm() << std::endl;

    std::vector<Vector3d> cell_centers_rot; cell_centers_rot.resize(Nc);
    for(int i=0; i<Nc; i++) {cell_centers_rot[i]= inv_RyRx*cell_centers[i];}

    // std::cout << "Gamma:\n" << Gamma << "\n MM:\n " << MM << "\nOmega:\n" << Omega << std::endl;
    std::vector<double> cells_beta; cells_beta.resize(Nc);
    std::vector<double> cells_qval; cells_qval.resize(Nc);
    std::vector<double> rot_cells_theta; rot_cells_theta.resize(Nc);
    std::vector<double> rot_cells_phi; rot_cells_phi.resize(Nc);
    std::vector<Vector3d> rot_cells_q_vec; rot_cells_q_vec.resize(Nc);
    std::vector<Vector3d> rot_cells_p; rot_cells_p.resize(Nc);
    for(int i=0; i<Nc; i++)
    {
      auto cq3d = cells_Q_tensor[i];
      Matrix3d Q3d_eigen;
      Q3d_eigen(0,0)= cq3d[0]; Q3d_eigen(0,1)= cq3d[1]; Q3d_eigen(0,2)= cq3d[2];
      Q3d_eigen(1,0)= cq3d[3]; Q3d_eigen(1,1)= cq3d[4]; Q3d_eigen(1,2)= cq3d[5];
      Q3d_eigen(2,0)= cq3d[6]; Q3d_eigen(2,1)= cq3d[7]; Q3d_eigen(2,2)= cq3d[8];
      //Matrix3d rot_Q = transform_matrices_from_xy_plane_to_triangle_plane(Q3d_eigen, rot_ang_X, rot_ang_Y);

      // Matrix3d rot_Q2 = RyRx*Q3d_eigen*RyRx.transpose();
      // std::cout << rot_Q2-rot_Q << std::endl;

      auto [qval, q_vec] = calculate_3x3_tensor_eigenvalues_and_eigenvector(Q3d_eigen, 1);
      cells_qval[i] = qval;
      auto rot_q_vec = inv_RyRx*q_vec;
      auto rot_pol = inv_RyRx*cells_polarity[i];

      rot_cells_q_vec[i] = rot_q_vec;
      rot_cells_p[i] = rot_pol;
      // std::cout << (cell_centers[i]/cell_centers[i].norm()).dot(q_vec/q_vec.norm()) << " "
      // << (cell_centers_rot[i]/cell_centers_rot[i].norm()).dot(rot_q_vec/rot_q_vec.norm()) << std::endl;
      auto [beta, tta, phi] = get_tilt_angle(rot_q_vec, cell_centers_rot[i]);
      cells_beta[i] = beta*180/M_PI;
      rot_cells_theta[i] = tta;
      rot_cells_phi[i] = phi;
    }

    std::string json_analysis_file = analysis_folder + "/cells_beta_fr"+std::to_string(fr)+".json";
    std::ofstream json_out_analysis(json_analysis_file);
    json json_analysis_data;
    json_analysis_data["frame"] = fr;
    json_analysis_data["cells_beta"] = cells_beta;
    json_analysis_data["cells_qval"] = cells_qval;
    json_analysis_data["rot_cells_q_vec"] = rot_cells_q_vec;
    json_analysis_data["rot_cells_pol"] = rot_cells_p;
    json_analysis_data["Omega"] = Omega;
    json_analysis_data["rot_cells_theta"] = rot_cells_theta;
    json_analysis_data["rot_cells_phi"] = rot_cells_phi;
    json_out_analysis  << json_analysis_data << std::endl;
    json_out_analysis.close();
    //=----------------
    // std::cout << Omega_rotated.norm() << " " << b_45 << " " << b_135 << std::endl;

  }

  return 0;
}
