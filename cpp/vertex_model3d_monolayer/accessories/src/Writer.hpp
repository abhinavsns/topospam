  #ifndef WRITER_HPP_
  #define WRITER_HPP_

  #include <iostream>
  #include <fstream>
  #include <filesystem>
  // #include <vector>
  #include "json.hpp"
  using json = nlohmann::json;

  namespace fs = std::filesystem;

  #include "Tissue.hpp"
  #include "Triangulation.hpp"

  /** Base class for writing post prcessed data.
    *
    * Should take care of:
    */
  class Writer : public Tissue
  {
  public:
    //
    // Writer() = default;

    Writer(unsigned, std::string);

    virtual ~Writer();
    /** explain
      *
      * This function is called
      * */
    unsigned frame = 0; // frame number
    std::string output_path = "postprocessed/"; // the output directory path

    bool write_apical_polygonal_surface = false;
    bool write_basal_polygonal_surface = false;
    bool write_apical_triangulated_surface = false;
    bool write_basal_triangulated_surface = false;
    bool write_full_3D_triangulated_surface = false;
    bool write_polarity_field = false;
    bool write_nematic_field = false;

    //std::unique_ptr<Triangulation> apical_tri_surface;//(new Triangulation());
    Triangulation apical_tri_surface = Triangulation();

    // now, a bunch of memebr functions:...
    void rotate_tissue_to_align_z_with_omega(Vector3d omega_);

    void rotate_tissue_around_z_axis(double angle_);

    void build_apical_tissue_triangulation(bool update_apical_geometry);

    void write_apical_surface_vector_spherical_harmonics_modes_to_json(int Lmax_);

    void write_apical_polarity_vector_spherical_harmonics_modes_to_json(int Lmax_);

    void write_apical_or_basal_polygonal_surface_to_vtk(std::string surf_type);

    void write_apical_or_basal_triangulated_surface_to_vtk(std::string surf_type);

    void write_full_3D_triangulation_to_vtk();

    void write_polarity_field_to_vtk();

    void write_nematicity_field_to_vtk();

    void write_data_to_vtk();

    void write_Patches_data();

  };

  #endif//WRITE_HPP_
