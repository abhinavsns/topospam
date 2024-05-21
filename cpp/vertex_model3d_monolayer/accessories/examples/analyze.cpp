#include "header.hpp"
#include "analyze_options.hpp"
#include "tools.hpp"
#include "error_msg.hpp"

#include <chrono>

#include "Writer.hpp"

#define T1TRANISTION 1


unsigned verbose = 2;
/** name of the output dir */
std::string outname;
// /** data dir */
std::string data_dir;
/** write any output? */
bool no_write = false;
//extern bool ,, , , align_before_writing;
bool analyze_elongation = false;
bool compute_curvature_tensor = false;
bool write_to_vtk = true;
/** skip runtime warnings? */
bool analyze_rotation = false;

bool analyze_apical_surface_VSH = false;

bool analyze_cell_polarity_VSH = false;

bool align_before_writing = false;

bool analysis_in_corotating_frame = false;

bool write_apical_polygonal_surface = false;
bool write_basal_polygonal_surface = false;
bool write_apical_triangulated_surface = false;
bool write_basal_triangulated_surface = false;
bool write_full_3D_triangulated_surface = false;
bool write_polarity_field = false;
bool write_nematic_field = false;
//...... dtaa frame range
unsigned first_frame, last_frame;
// ......Vector pshericla harmmonics modes related ...
unsigned Lmax;

//---------- Defining the algorithms for post-processing -------------
void Write_VTK_Files(std::unique_ptr<Writer> &VM_writer)
{
  VM_writer->write_apical_polygonal_surface = write_apical_polygonal_surface;
  VM_writer->write_basal_polygonal_surface = write_basal_polygonal_surface;
  VM_writer->write_apical_triangulated_surface = write_apical_triangulated_surface;
  VM_writer->write_basal_triangulated_surface = write_basal_triangulated_surface;
  VM_writer->write_full_3D_triangulated_surface = write_full_3D_triangulated_surface;
  VM_writer->write_polarity_field = write_polarity_field;
  VM_writer->write_nematic_field = write_nematic_field;

  VM_writer->write_data_to_vtk();
}

std::pair<Vector3d,double> GetAngularVelocity(unsigned fr_number, bool use_prev_frame_nr=true)
{
  // ... let's obtain the time step used in the siulation ....
  std::string json_param_file = data_dir + "parameters.json";
  std::ifstream parfile(json_param_file);
  json param_data = json::parse(parfile);
  double dt = param_data["dt"];
  // .........
  std::string json_file = data_dir + "frame_" + std::to_string(fr_number) + ".json";
  std::ifstream f(json_file);
  json json_data = json::parse(f);
  int Nv=json_data["Num_vertices"];
  // int Nc=json_data["Num_cells"];
  auto vertex_coord_a = json_data["vertex_pos_a"];
  std::vector<std::vector<double>> vertex_prev_coord_a;
  if(use_prev_frame_nr)
  {
    // int ninfo = json_data["output_freq"];
    // dt *= (double)ninfo;
    dt = 1;
    std::string json_file0 = data_dir + "frame_" + std::to_string(fr_number-1) + ".json";
    std::ifstream f0(json_file0);
    json json_data = json::parse(f0);
    vertex_prev_coord_a = json_data["vertex_pos_a"];
  }
  else
  {
    vertex_prev_coord_a = json_data["vertex_prev_pos_a"];
  }


  std::vector<Vector3d> pos_t1(Nv), pos_t2(Nv);
  for(int i=0; i<Nv; i++)
  {
    pos_t1[i] = Vector3d(vertex_prev_coord_a[i][0], vertex_prev_coord_a[i][1], vertex_prev_coord_a[i][2]);
    pos_t2[i] = Vector3d(vertex_coord_a[i][0], vertex_coord_a[i][1], vertex_coord_a[i][2]);
  }
  auto [Gamma, MM] = get_ang_momentum_inertia_tensor(pos_t1, pos_t2, dt);
  Vector3d Omega = MM.inverse() * Gamma;

  return std::make_pair(Omega, dt);
}

void Algorithm()
{
  // //.......... now run the actual simulation steps with active forces included
  // for(size_t i=0; i<organoid->ListCell.size(); i++) organoid->ListCell[i].mag_F_a = F_mag;
  //
  // for(unsigned fr=0; fr<nframes; fr++)
  // {
  //   if(fr==noiseoff) P_noise=0;
  //
  //   for(unsigned i=0; i<ninfo; i++)
  //   {
  //     organoid->evolve_organoid_one_step(P_nu, P_gamma, P_noise, time_step, xi, P_A, elastic_polarity, spherical_constrain);
  //     #if T1TRANISTION
  //       organoid->flip_short_bonds(bond_T1_cutoff, bond_T1_opening, true);
  //     #endif
  //   }
  //   if(!no_write) organoid->write_to_json_file(outname, fr+1);
  // }

}
//=====================================================================
// ======================== MAIN FILE =================================
//=====================================================================
int main(int argc, char *argv[])
{
  try
  {
    if (argc<2) throw error_msg("You need to identify the arguments. Type -h for help.");
    ParseProgramOptions(argc, argv);
    ProcessProgramOptions();
    PrintParameters();
    if(!no_write) write_parameters_to_jsonfile();

    const auto start = std::chrono::steady_clock::now();
    for(unsigned fr=first_frame; fr<=last_frame; fr++)
    {
      std::cout << fr << std::endl;
      // ...................
      std::unique_ptr<Writer> writer_obj(new Writer(fr, outname)); // a unique ptr

      std::string jf = data_dir + "frame_"+std::to_string(fr)+".json";

      writer_obj->build_3D_spherical_tissue_from_json_file(jf);

      if(align_before_writing)
      {
        Vector3d omega_fr; double dt_;
        std::tie(omega_fr, dt_) = GetAngularVelocity(fr, 1);
        writer_obj->rotate_tissue_to_align_z_with_omega(omega_fr);
      }

      // if(align_before_writing)
      // {
      //   Vector3d omega_fr; double dt_;
      //   std::tie(omega_fr, dt_) = GetAngularVelocity(fr, 1);
      //   writer_obj->rotate_tissue_around_z_axis( - dt_ * omega_fr.norm());
      // }

      writer_obj->build_apical_tissue_triangulation(true);
      // if(analyze_elongation) // these will be done in write_Patches_data()
      // {
      //   writer_obj->apical_tri_surface.update_all_Patches_Q_3D();
      //   writer_obj->apical_tri_surface.update_all_elongation_tilt_angles_if_already_aligned();
      // }
      if(analyze_apical_surface_VSH) writer_obj->write_apical_surface_vector_spherical_harmonics_modes_to_json(Lmax);
      if(analyze_cell_polarity_VSH) writer_obj->write_apical_polarity_vector_spherical_harmonics_modes_to_json(Lmax);
      if(!no_write) writer_obj->write_Patches_data();

      if(!no_write) Write_VTK_Files(writer_obj);

      //double p_energy = writer_obj->get_tissue_polarity_elastic_energy();
      //std::cout << p_energy << std::endl;

      //std::cout << "Nc:  " << writer_obj->ListCell.size()  << "  nv:  " << writer_obj->ListVertex.size()<< std::endl;
    }

    const auto duration = std::chrono::steady_clock::now() - start;

    if(verbose)
    {
      std::cout << "Analysis run time :                    "
           << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count()
              /1000. << " s" << std::endl;
    }

    //Tissue *tissue = new Tissue(Nc); // heap memory allocation


    // std::unique_ptr<Tissue> tissue(new Tissue(Nc)); // a unique ptr to avoid the memory leak duto an accidental crash
    // tissue->verbose = verbose;
    //
    // const auto prep_start = std::chrono::steady_clock::now();
    // PrepareInitialTissue(tissue);
    // const auto prep_duration = std::chrono::steady_clock::now() - prep_start;
    //
    // const auto start = std::chrono::steady_clock::now();
    // Algorithm(tissue);
    // const auto duration = std::chrono::steady_clock::now() - start;
    //
    // if(verbose) std::cout << "done" <<std::endl;
    //
    // if(verbose)
    // {
    //   std::cout << std::endl << "Statistics" << std::endl << std::string(width, '=') << std::endl;
    //   std::cout << "Preparation run time :                    "
    //        << std::chrono::duration_cast<std::chrono::milliseconds>(prep_duration).count()
    //           /1000. << " s" << std::endl;
    //
    //   std::cout << "Algorithm run time :                    "
    //        << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count()
    //           /1000. << " s" << std::endl;
    // }

  }
  // custom small messages
  catch(const error_msg& e) {
    std::cerr << argv[0] << ": error: " << e.what() << std::endl;
    return 1;
  }
  // bad alloc (possibly from memory())
  catch(const std::bad_alloc& ba) {
    std::cerr << argv[0] << ": error initializing memory: " << ba.what() << std::endl;
    return 1;
  }
  // all the rest (mainly from boost)
  catch(const std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    return 1;
  }


  return 0;
}
