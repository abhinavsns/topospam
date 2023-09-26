#include "header.hpp"
#include "options.hpp"
#include "tools.hpp"
#include "error_msg.hpp"
//#include "random.hpp"

using namespace std;
//namespace opt = boost::program_options;

//================================
extern double Kc, A0c, bond_Tension, P_noise, P_A, noise, perim_Elasticity, F_mag, P_gamma, time_step, xi, P_nu;
extern double bond_T1_cutoff, bond_T1_opening, P0_a0_by_PI, P0_b0_by_PI;
//using option_list = std::vector<opt::options_description>;
// =============================================================================
// Options

// declare variables externally
// parameters
extern unsigned verbose, nthreads, nframes, ninfo, noiseoff, Nc, run_id;

extern bool no_write, P0_axissymmetric, read_initsurface, spherical_constrain, elastic_polarity;

extern string outname;

/** the variables map used to collect options */
opt::variables_map vm;
/** name of the runcarde file */
string runcardname = "";

bool force_delete;

unsigned seed;

std::vector<opt::options_description> GetOptions()
{
  // model specific options
  opt::options_description model_options("Model options");
  model_options.add_options()
    ("Kc", opt::value<double>(&Kc),
     "Area stiffness")
    ("A0c", opt::value<double>(&A0c),
     "cell prefereed area")
    ("bond_Tension", opt::value<double>(&bond_Tension),
     "Cell bond tension")
    ("perim_Elasticity", opt::value<double>(&perim_Elasticity),
     "Cell perimeter elasticity")
    ("F_mag", opt::value<double>(&F_mag),
      "Cell traction force magnitude")
    ("P_gamma", opt::value<double>(&P_gamma),
      "Rate of cell polarity alignment with nieighbors")
    ("time_step", opt::value<double>(&time_step),
      "size of each time step")
    ("xi", opt::value<double>(&xi),
      "friction coefficient of vertex with extern environemnt.")
    ("P_nu", opt::value<double>(&P_nu),
      "rate of alignment of polarity with cell velocity.")
    ("P_noise", opt::value<double>(&P_noise),
      "Diffusion coefficient of polarity rotational noise")
    ("P_A", opt::value<double>(&P_A),
      "stiffness of polarity norm.");

  opt::options_description config_options("Initial configuration options");
  config_options.add_options()
    ("noise", opt::value<double>(&noise),
     "size of initial variations")
    ("bond_T1_cutoff", opt::value<double>(&bond_T1_cutoff),
      "T1 cutoff bond lenght")
    ("bond_T1_opening", opt::value<double>(&bond_T1_opening),
      "T1 bond opening length")
    ("Nc", opt::value<unsigned>(&Nc),
      "Number of cells at t=0")
    ("run_id", opt::value<unsigned>(&run_id),
      "run id for each randomize surface at t=0")
    ("P0_a0_by_PI", opt::value<double>(&P0_a0_by_PI),
      "a_0 for axissymetric polarity pattern at t=0")
    ("P0_b0_by_PI", opt::value<double>(&P0_b0_by_PI),
      "b_0 for axissymetric polarity pattern at t=0")
    ("P0_axissymmetric", opt::bool_switch(&P0_axissymmetric),
      "if true, axissetric polarity initialization.")
    ("read_initsurface", opt::bool_switch(&read_initsurface)->default_value(true),
      "true=reads from given inputsurface address; false: generats random surface with Nc cells.")
    ("spherical_constrain", opt::bool_switch(&spherical_constrain)->default_value(true),
      "if true, vertices contrained on spherical geomtery. [default=true]")
    ("elastic_polarity", opt::bool_switch(&elastic_polarity)->default_value(true),
      "if fasle polarity norm is contrianed to 1. [default=true]");

  return { model_options, config_options };
}
//===================================================


void ParseProgramOptions(int ac, char **av)
{
  // options allowed only in the command line
  opt::options_description generic("Generic options");
  generic.add_options()
    ("help,h",
     "produce help message")
    ("verbose,v", opt::value<unsigned>(&verbose)->implicit_value(2),
     "verbosity level (0=none, 1=little, 2=normal, 3=debug)")
    ("runcard,i", opt::value<string>(&runcardname),
     "runcard file")
    ("force-delete,f", opt::bool_switch(&force_delete)->default_value(true),
     "force deletion of existing output file. [default=true]")
    ("threads,t",
     opt::value<unsigned>(&nthreads)->default_value(0)->implicit_value(1),
     "number of threads (0=no multithreading, 1=OpenMP default, "
     ">1=your favorite number)")
    ("no-write", opt::bool_switch(&no_write)->default_value(false),
     "disable file output (for testing purposes). [default=false]");

  // options allowed both in the command line and config file
  opt::options_description config("Program options");
  config.add_options()
    ("output,o", opt::value<string>(&outname),
     "output name")
    ("nframes", opt::value<unsigned>(&nframes),
     "iterate untill it writes this many frames in total")
    ("ninfo", opt::value<unsigned>(&ninfo),
     "save frame every so many steps")
    ("seed", opt::value<unsigned>(&seed),
     "set seed for random number generation (random if unset)")
    ("noiseoff", opt::value<unsigned>(&noiseoff)->default_value(std::numeric_limits<unsigned>::max()),
     "starting from this frame noise turns off");

  // command line options
  opt::options_description cmdline_options;
  cmdline_options.add(generic).add(config);
  opt::options_description config_file_options;
  config_file_options.add(config);

  // first unnamed argument is the runcard file
  opt::positional_options_description p;
  p.add("runcard", 1);

  // reintialize vm in case we run this function twice
  vm = opt::variables_map();

  // parse first the cmd line to get model name (no throw)
  opt::store(
    opt::command_line_parser(ac, av)
    .options(cmdline_options)
    .positional(p)
    .allow_unregistered()
    .run(), vm);
  opt::notify(vm);

  // print help msg and exit
  if(vm.count("help") and runcardname.empty())
  {
    cout << cmdline_options << endl;
    exit(0);
  }

  // parse runcard file (values are not erased, such that cmd line args
  // are 'stronger') (no throw)
  if(not runcardname.empty())
  {
    std::fstream file(runcardname.c_str(), std::fstream::in);
    if(!file.good()) throw error_msg("can not open runcard file ", runcardname);
    opt::store(opt::parse_config_file(file, config_file_options, true), vm);
    opt::notify(vm);
  }

  // print help msg and exit
  if(vm.count("help"))
  {
    cout << cmdline_options << endl;
    exit(0);
  }


auto model_options = GetOptions();

for(const auto& o : model_options)
{
  cmdline_options.add(o);
  config_file_options.add(o);
}

// parse the cmd line a second time with model options (throw)
opt::store(
  opt::command_line_parser(ac, av)
  .options(cmdline_options)
  .positional(p)
  .run(), vm);
opt::notify(vm);

// print help msg and exit
if(vm.count("help"))
{
  cout << cmdline_options << endl;
  exit(0);
}

// parse runcard file (values are not erased, such that cmd line args
// are 'stronger')
if(runcardname.empty())
  throw error_msg("please provide an runcard file / type -h for help.");
else
{
  std::fstream file(runcardname.c_str(), std::fstream::in);
  if(!file.good()) throw error_msg("can not open runcard file ", runcardname);
  opt::store(opt::parse_config_file(file, config_file_options), vm);
  opt::notify(vm);
}

}

void ProcessProgramOptions()
{
  if(vm.count("output")==0)
  {
    outname = "output/";
  }
  if (!outname.empty() && outname.back() != '/')
    outname += '/';

  namespace fs = std::filesystem;
  const fs::path out_path{outname};
  fs::create_directory(out_path);

  if(force_delete)
  {
    for (const auto& entry : std::filesystem::directory_iterator(out_path))
      std::filesystem::remove_all(entry.path());
  }
}
//
// /** Print variables from variables_map
//   *
//   * from: https://gist.github.com/gesquive/8673796
//   */
void print_vm(const opt::variables_map& vm, unsigned padding)
{
  for (opt::variables_map::const_iterator it = vm.begin(); it != vm.end(); ++it)
  {
    // pass if defaulted or empty
    if (vm[it->first].defaulted() || it->second.defaulted()) continue;
    if (((boost::any)it->second.value()).empty()) continue;

    std::cout << std::left << std::setw(floor(padding/2)) << it->first;

    /*if (((boost::any)it->second.value()).empty()) {
      std::cout << "(empty)";
    }
    if (vm[it->first].defaulted() || it->second.defaulted()) {
      std::cout << "(default)";
    }*/

    std::cout << std::right << std::setw(ceil(padding/2));

    bool is_char;
    try {
      boost::any_cast<const char*>(it->second.value());
      is_char = true;
    } catch (const boost::bad_any_cast &) {
      is_char = false;
    }
    bool is_str;
    try {
      boost::any_cast<std::string>(it->second.value());
      is_str = true;
    } catch (const boost::bad_any_cast &) {
      is_str = false;
    }

    if (((boost::any)it->second.value()).type() == typeid(int)) {
      std::cout << vm[it->first].as<int>() << std::endl;
    } else if (((boost::any)it->second.value()).type() == typeid(unsigned)) {
      std::cout << vm[it->first].as<unsigned>() << std::endl;
    } else if (((boost::any)it->second.value()).type() == typeid(size_t)) {
      std::cout << vm[it->first].as<size_t>() << std::endl;
    } else if (((boost::any)it->second.value()).type() == typeid(bool)) {
      std::cout << (vm[it->first].as<bool>() ? "true" : "false") << std::endl;
    } else if (((boost::any)it->second.value()).type() == typeid(double)) {
      std::cout << vm[it->first].as<double>() << std::endl;
    } else if (((boost::any)it->second.value()).type()
               == typeid(vector<double>)) {
      std::cout << vec2str(vm[it->first].as<vector<double>>()) << std::endl;
    } else if (((boost::any)it->second.value()).type()
               == typeid(vector<unsigned>)) {
      std::cout << vec2str(vm[it->first].as<vector<unsigned>>()) << std::endl;
    } else if (is_char) {
      std::cout << vm[it->first].as<const char *>() << std::endl;
    } else if (is_str) {
      std::string temp = vm[it->first].as<std::string>();
      if (temp.size()) {
        std::cout << temp << std::endl;
      } else {
        std::cout << "true" << std::endl;
      }
    } else { // Assumes that the only remainder is vector<string>
      try {
        auto vect = vm[it->first].as<std::vector<std::string> >();
        unsigned int i = 0;
        for (auto oit=vect.begin();
            oit != vect.end(); oit++, ++i) {
          std::cout << "\r> " << it->first
                    << "[" << i << "]=" << (*oit) << std::endl;
        }
      } catch (const boost::bad_any_cast &) {
        std::cout << "UnknownType("
                  << ((boost::any)it->second.value()).type().name() << ")" << std::endl;
      }
    }
  }
}
//
void PrintParameters()
{
  // print the simulation parameters
  if(verbose)
  {
    cout << "Run parameters" << endl;
    cout << string(width, '=') << endl;
    print_vm(vm, width);
    cout << string(width, '.') << endl;
  }
}


void write_parameters_to_jsonfile()
{
    using json = nlohmann::json;
    std::string json_file = outname + "parameters.json";
    std::ofstream json_out(json_file);
    json json_data;
    json_data["Nc"] = Nc; json_data["runcard"] = runcardname;
    json_data["F_mag"] = F_mag; json_data["P_gamma"] = P_gamma;
    json_data["P_nu"]=P_nu; json_data["P_noise"] = P_noise;
    json_data["output_folder"] = outname; json_data["xi"]=xi; json_data["dt"]=time_step;
    json_data["output_freq"] = ninfo; json_data["total_frames"] = nframes; json_data["noise_off_frame"]=noiseoff;
    json_data["bond_T1_cutoff"] = bond_T1_cutoff; json_data["bond_T1_opening"] = bond_T1_opening;
    json_data["Kc"]=Kc; json_data["A0c"]=A0c; json_data["bond_Tension"]=bond_Tension; json_data["perim_Elasticity"]=perim_Elasticity;
    json_data["P0_axissymmetric"] = P0_axissymmetric; json_data["P0_a0_by_PI"]=P0_a0_by_PI; json_data["P0_b0_by_PI"]=P0_b0_by_PI;
    json_data["run_id"]=run_id; json_data["read_initsurface"]=read_initsurface;
    json_data["spherical_constrain"] = spherical_constrain; json_data["elastic_polarity"]=elastic_polarity; json_data["P_A"]=P_A;
    json_out  << json_data << std::endl;
    json_out.close();
    //=----------------
}
