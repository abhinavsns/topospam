#ifndef TISSUE_HPP_
#define TISSUE_HPP_

// #include <boost/math/constants/constants.hpp>

#include "Vertex.hpp"
#include "Edge.hpp"
#include "Cell.hpp"
#include "Common-inline.h"
#include "Bond-inline.h"

#include <gsl/gsl_multimin.h>

// #define M_PI boost::math::constants::pi<double>()

#define THREE_DIMENSIONS 0

//#define PRALLEL 0

#if PRALLEL
	#include <parallel/algorithm>
	#include <algorithm>
	#include <execution>
#endif

#ifndef NORM_EPS
#define NORM_EPS 1e-15
#endif

struct CellMechanicalProperties
{
	double Ta; // apical surface tension
  double Tb; // apical surface tension
  double Tl; // apical surface tension
	double V0c; // Cell preferred volume
	double Kc; // Cell bulk modulus
	double A0c;//preferred area. used for 2D curved tissue
  double BondT; // Bond tension
	double perim_elasticity; // perimeter elasticity, mostly used for 2D model
};

struct SpheroidMechanicalProperties
{
	double P0; // pressue difference accross spheroid
  double Kl; // lumen volume elasticity
	double V0l; // lumen preferred volume
};

struct MinimizationParameters
{
    int MinimizationMode; // 1 - Fletcher Reese, 2 - Polak Ribiere
    double ftol; // convergence criterion - relative error in successive energy steps
    double GTOL; // convergence criterion - max gradient contribution
    int rateT1check; // number of steps after which the system is checked for T1s
    double tol; // tolerance in line minimization
    int ITMAX; // maximum number of iterations before abort
    double EPS; //
};

class Tissue
{
public:
  // ---- CREATION & DESTRUCTION ----
  Tissue();

	Tissue(int N0);
  // Vertex(bool is_3D);      //constructor
  virtual ~Tissue();     //destructor

  std::vector<Vertex> ListVertex;
  std::vector<Edge> ListEdge;
  std::vector<Cell> ListCell;
  /// constructor that takes the coordinates of the apical and basal point
  // Vertex(Vector3d pos_a, Vector3d pos_b);

  CellMechanicalProperties CellMechProp;
	SpheroidMechanicalProperties SphMechProp;
	MinimizationParameters minPar;

	double W_frame;
	double lumen_volume;

	std::vector<double> vpos;
	std::vector<double> vder;

	unsigned verbose;

  void set_CellMechProps3D(double Kc_, double V0c_, double Ta_, double Tl_, double Tb_);
	void set_SpheroidMechProp(double P0_, double Kl_, double V0l_);
  void set_CellMechProps2D(double Kc_, double A0c_, double BondT_, double perim_elasticity_);

  void build_2D_spherical_tissue_from_file(const char* fileName);
	void build_2D_dodecahedron_tissue();

  Vertex* similar_or_new_vertex_2D(Vector3d pos);
  Edge* similar_or_new_edge(Vertex* v_tail, Vertex* v_head);

	double get_cell_volume(Cell *c);
	Vector3d get_cell_centroid(Cell *c);
  Vector3d get_cell_apical_centroid(Cell *c);
	Vector3d get_cell_basal_centroid(Cell *c);
  double get_cell_apical_area(Cell *c);
	double get_cell_basal_area(Cell *c);
	double get_cell_apical_perimeter(Cell *c);
	double get_cell_basal_perimeter(Cell *c);

	Matrix3d get_apical_elongation_tensor(Cell *c);

	std::pair<std::vector<Vector3d>, std::vector<Vector3d>> get_vertex_pos_changes_of_cell(Cell *c);

  void rescale_spheroid_apical_size();

	double get_lumen_volume();

	double get_apical_bond_length(Edge *e);
	double get_basal_bond_length(Edge *e);

	double get_vertex_lateral_distance(Vertex *v);
	Vector3d get_edge_lateral_centroid(Edge *e);
	double get_lateral_area(Edge *e);

	void update_all_apical_bond_lengths();
	void update_all_basal_bond_lengths();

	void update_apical_geometric_quantities();
	void update_basal_geometric_quantities();
	void update_lateral_geometric_quantities();
	void update_3D_geometric_quantities();

	std::vector<int> find_short_bonds(double bond_T1_cutoff, bool random_shuffling);

	bool flip_edge(Edge *e, double opening_length);

	bool flip_short_bonds(double bond_T1_cutoff, double opening_length, bool random_shuffling=true);

	bool divide_cell_random_axis(Cell *M);

	bool divide_N_random_cell_with_random_axis(int N_d);

	double W_2D();
	void update_dW_a_2D();
	double W_3D();
	void update_dW_3D_spherical();// will take lumen arguments later

	void evolve_vertex_positions(double xi, double dt, bool spherical_constrain);

	void set_apical_points_on_sphere_at_radius(double R_sp);

	// void set_minimization_params();
	// void update();
	// void set_pos();
	// void update_pos();
	// void update_pos_and_der();
	// void bracket(std::vector<double> &dir, double &ax, double &bx, double &cx, double &fb, double &current_u);
	// double line_minimization_Brent(std::vector<double> &dir);
	// int my_cg_minimizer();
	//
	double gsl_f (const gsl_vector *v);
	double gsl_f3D (const gsl_vector *v);
	void gsl_df (const gsl_vector *v, gsl_vector *df);
	void gsl_df3D (const gsl_vector *v, gsl_vector *df);
	void gsl_fdf (const gsl_vector *x, double *f, gsl_vector *df);
	void gsl_fdf3D (const gsl_vector *x, double *f, gsl_vector *df);
	//
	void minimize_gsl3D();
	void minimize_gsl();
	// void minimize_gsl2();

	void randomize_cell_planar_polarity();
	void set_axissymmetric_polarity(double a0, double b0);
	void update_average_nearest_neighbor_polarity(bool _elastic_polarity);
	void evolve_polarity(double lambda, double gamma, double diffD, double dt, double _P_A, bool _elastic_polarity);
	void update_active_forces();
	void evolve_organoid_one_step(double nu, double P_gamma, double diffD, double dt, double xi, double _P_A, bool _elastic_polarity, bool spherical_constrain);
	double get_tissue_polarity_elastic_energy();

	// ====== 3D stuff -------------------
	Vertex* similar_or_new_vertex_3D(Vector3d _pos_a, Vector3d _pos_b);
	void build_3D_spherical_tissue_from_file(const char* fileName, double R_sph, double Rb_Ra_ratio);
	void build_3D_spherical_tissue_from_json_file(std::string json_file);
	void build_3D_spherical_tissue_from_output_files(const char* vertices_file, const char* cells_file, const char* cell_data_file);


  void write_to_txt_files(std::string output_folder, int frame);

	void write_to_json_file(std::string output_folder, int frame);
};

#endif//TISSUE_HPP_
