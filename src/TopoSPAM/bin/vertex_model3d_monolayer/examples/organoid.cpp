#include "header.hpp"
#include "options.hpp"
#include "tools.hpp"
#include "error_msg.hpp"

#include <chrono>

#include "Tissue.hpp"

#define T1TRANISTION 1


unsigned verbose = 2;
/** name of the output dir */
std::string outname;
// /** Output dir (or tmp dir before moving files to the archive) */
// std::string output_dir;
/** write any output? */
bool no_write = false;
bool P0_axissymmetric = false;
bool read_initsurface = true;
/** skip runtime warnings? */
bool no_warning = false;

bool spherical_constrain = true;

bool elastic_polarity = true;

unsigned nthreads;
/** Total number of time steps */
unsigned nframes;
/** Time interval between data outputs */
unsigned ninfo;
/** Time at which to start the output */
unsigned noiseoff;
/** number of cells */
unsigned Nc, run_id;

double Kc, A0c, bond_Tension, P_noise, P_A, noise, perim_Elasticity, F_mag, P_gamma, time_step, xi, P_nu;
double bond_T1_cutoff, bond_T1_opening;
double P0_a0_by_PI, P0_b0_by_PI;
std::string inputfileadress = "/data/biophys/aamiri/y2023/VM/v1/2/vertex_model3d_monolayer/accessories/misc/input_files/";

//---------- Defining the algorithms for organoid dynamics -------------
void PrepareInitialTissue(std::unique_ptr<Tissue> &organoid)
{
  if(read_initsurface)
  {
    std::string input_file = inputfileadress + "random/"+std::to_string(Nc)+"/random" + std::to_string(Nc)+"v"+std::to_string(run_id)+".txt";
    organoid->build_2D_spherical_tissue_from_file(input_file.c_str());
  }
  else
  {
    organoid->build_2D_dodecahedron_tissue();
    organoid->divide_N_random_cell_with_random_axis(Nc-12);
  }

  organoid->set_CellMechProps2D(Kc, A0c, bond_Tension, perim_Elasticity);

  organoid->rescale_spheroid_apical_size();

  organoid->update_apical_geometric_quantities();

  if(P0_axissymmetric) organoid->set_axissymmetric_polarity(P0_a0_by_PI*M_PI, P0_b0_by_PI*M_PI);
  else organoid->randomize_cell_planar_polarity();
  //............first prepare the tissue without active forces -----
  for(size_t i=0; i<organoid->ListCell.size(); i++) organoid->ListCell[i].mag_F_a = 0.;

  double time_scale = xi * std::sqrt(A0c) / bond_Tension;
  int steps_for_one_time_scale = (int)(time_scale/time_step);
  for(int i=0; i<2*steps_for_one_time_scale; i++)
  {
    organoid->evolve_organoid_one_step(P_nu, P_gamma, P_noise, 2.0*time_step, xi, P_A, elastic_polarity, spherical_constrain);
    #if T1TRANISTION
      organoid->flip_short_bonds(bond_T1_cutoff, bond_T1_opening, true);
    #endif
  }

  if(!no_write) organoid->write_to_json_file(outname, 0);

}

void Algorithm(std::unique_ptr<Tissue> &organoid)
{
  //.......... now run the actual simulation steps with active forces included
  for(size_t i=0; i<organoid->ListCell.size(); i++) organoid->ListCell[i].mag_F_a = F_mag;

  for(unsigned fr=0; fr<nframes; fr++)
  {
    if(fr==noiseoff) P_noise=0;

    for(unsigned i=0; i<ninfo; i++)
    {
      organoid->evolve_organoid_one_step(P_nu, P_gamma, P_noise, time_step, xi, P_A, elastic_polarity, spherical_constrain);
      #if T1TRANISTION
        organoid->flip_short_bonds(bond_T1_cutoff, bond_T1_opening, true);
      #endif
    }
    if(!no_write) organoid->write_to_json_file(outname, fr+1);
  }

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

    //Tissue *tissue = new Tissue(Nc); // heap memory allocation


    std::unique_ptr<Tissue> tissue(new Tissue(Nc)); // a unique ptr to avoid the memory leak duto an accidental crash
    tissue->verbose = verbose;

    const auto prep_start = std::chrono::steady_clock::now();
    PrepareInitialTissue(tissue);
    const auto prep_duration = std::chrono::steady_clock::now() - prep_start;

    const auto start = std::chrono::steady_clock::now();
    Algorithm(tissue);
    const auto duration = std::chrono::steady_clock::now() - start;

    if(verbose) std::cout << "done" <<std::endl;

    if(verbose)
    {
      std::cout << std::endl << "Statistics" << std::endl << std::string(width, '=') << std::endl;
      std::cout << "Preparation run time :                    "
           << std::chrono::duration_cast<std::chrono::milliseconds>(prep_duration).count()
              /1000. << " s" << std::endl;

      std::cout << "Algorithm run time :                    "
           << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count()
              /1000. << " s" << std::endl;
    }

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
