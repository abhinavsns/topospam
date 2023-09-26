#include "header.hpp"
#include "analyze_options.hpp"
#include "tools.hpp"
#include "error_msg.hpp"
//#include "random.hpp"

using namespace std;
//namespace opt = boost::program_options;

//================================
//extern double A0c, bond_Tension, P_noise, P_A, perim_Elasticity, P_gamma, time_step, xi, P_nu;
//using option_list = std::vector<opt::options_description>;
// =============================================================================
// Options

// declare variables externally
// parameters
extern unsigned verbose, Lmax, first_frame, last_frame;

extern bool no_write, analyze_elongation, compute_curvature_tensor, write_to_vtk, analyze_rotation, align_before_writing, analysis_in_corotating_frame;
extern bool analyze_apical_surface_VSH, analyze_cell_polarity_VSH;

extern bool write_apical_polygonal_surface, write_basal_polygonal_surface;
extern bool write_apical_triangulated_surface, write_basal_triangulated_surface;
extern bool write_full_3D_triangulated_surface;
extern bool write_polarity_field, write_nematic_field;

extern string outname, data_dir;

/** the variables map used to collect options */
opt::variables_map vm;
/** name of the runcarde file */
string runcardname = "";

bool force_delete;

std::vector<opt::options_description> GetOptions()
{
  // model specific options
  opt::options_description model_options("Model options");
  model_options.add_options()
    ("Lmax", opt::value<unsigned>(&Lmax),
     "Maxiumum l mode decomposition for Vector Spherical Harmonics.")
    ("first_frame", opt::value<unsigned>(&first_frame)->default_value(0),
      "Fist frame number to be analyzed. [default=0]")
    ("last_frame", opt::value<unsigned>(&last_frame)->default_value(0),
      "Last frame number to be analyzed. [default=0]");

  opt::options_description config_options("Initial configuration options");
  config_options.add_options()
    ("analyze_elongation", opt::bool_switch(&analyze_elongation)->default_value(false),
      "if ture, cell elongation patter is analyzed and output-ed. [default=false]")
    ("compute_curvature_tensor", opt::bool_switch(&compute_curvature_tensor)->default_value(false),
      "if ture, coarse-grained curvature tensor is analyzed and output-ed. [default=false]")
    ("write_to_vtk", opt::bool_switch(&write_to_vtk)->default_value(false),
      "if true, vtk fiels will be produced for visualisation. [default=false]")
    ("analyze_rotation", opt::bool_switch(&analyze_rotation)->default_value(false),
      "if true, rotation angular velocity is computed and output-ed. [default=false]")
    ("analyze_apical_surface_VSH", opt::bool_switch(&analyze_apical_surface_VSH)->default_value(false),
      "if true, apical surface is decomposed to VSH and output-ed. [default=false]")
    ("analyze_cell_polarity_VSH", opt::bool_switch(&analyze_cell_polarity_VSH)->default_value(false),
      "if true, cell polarity is decomposed to VSH and output-ed. [default=false]")
    ("align_before_writing", opt::bool_switch(&align_before_writing)->default_value(false),
      "if true, aligns z axis with the rotation axis before output. [default=false]")
    ("analysis_in_corotating_frame", opt::bool_switch(&analysis_in_corotating_frame)->default_value(false),
      "if true, analysis will be done in corotating frame. [default=false]")
    ("write_apical_polygonal_surface", opt::bool_switch(&write_apical_polygonal_surface)->default_value(false),
      "write_apical_polygonal_surface. [default=false]")
    ("write_basal_polygonal_surface", opt::bool_switch(&write_basal_polygonal_surface)->default_value(false),
      "write_basal_polygonal_surface. [default=false]")
    ("write_apical_triangulated_surface", opt::bool_switch(&write_apical_triangulated_surface)->default_value(false),
      "write_apical_triangulated_surface. [default=false]")
    ("write_basal_triangulated_surface", opt::bool_switch(&write_basal_triangulated_surface)->default_value(false),
      "write_basal_triangulated_surface. [default=false]")
    ("write_full_3D_triangulated_surface", opt::bool_switch(&write_full_3D_triangulated_surface)->default_value(false),
      "write_full_3D_triangulated_surface. [default=false]")
    ("write_polarity_field", opt::bool_switch(&write_polarity_field)->default_value(false),
      "write_polarity_field. [default=false]")
    ("write_nematic_field", opt::bool_switch(&write_nematic_field)->default_value(false),
      "write_nematic_field. [default=false]");

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
    ("no-write", opt::bool_switch(&no_write)->default_value(false),
     "disable file output (for testing purposes). [default=false]");

  // options allowed both in the command line and config file
  opt::options_description config("Program options");
  config.add_options()
    ("output,o", opt::value<string>(&outname),
     "output name")
    ("data_dir,d", opt::value<string>(&data_dir),
      "directory path to the data to be analyzed.");

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
  if(data_dir.empty())
    throw error_msg("please provide path to data_dir / type -h for help.");
  else if(data_dir.back() != '/'){
    data_dir += '/';
  }

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

  if(no_write)
  {
    write_apical_polygonal_surface = false;
    write_basal_polygonal_surface = false;
    write_apical_triangulated_surface = false;
    write_basal_triangulated_surface = false;
    write_full_3D_triangulated_surface = false;
    write_polarity_field = false;
    write_nematic_field = false;
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
    json_data["runcard"] = runcardname; json_data["align_before_writing"]=align_before_writing;
    json_data["Lmax"] = Lmax;
    json_data["first_frame"] = first_frame; json_data["last_frame"]=last_frame;
    json_data["output_folder"] = outname;
    json_data["analyze_elongation"]=analyze_elongation;
    json_data["analysis_in_corotating_frame"] = analysis_in_corotating_frame;
    json_data["analyze_apical_surface_VSH"] = analyze_apical_surface_VSH;
    json_data["analyze_cell_polarity_VSH"] = analyze_cell_polarity_VSH;
    json_data["compute_curvature_tensor"] = compute_curvature_tensor;
    json_data["write_to_vtk"] = write_to_vtk; json_data["analyze_rotation"]=analyze_rotation;
    json_data["write_apical_polygonal_surface"] = write_apical_polygonal_surface;
    json_data["write_basal_polygonal_surface"] = write_basal_polygonal_surface;
    json_data["write_apical_triangulated_surface"] = write_apical_triangulated_surface;
    json_data["write_basal_triangulated_surface"] = write_basal_triangulated_surface;
    json_data["write_full_3D_triangulated_surface"] = write_full_3D_triangulated_surface;
    json_data["write_polarity_field"] = write_polarity_field;
    json_data["write_nematic_field"] = write_nematic_field;
    json_out  << json_data << std::endl;
    json_out.close();
    //=----------------
}
