#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include "../include/json.hpp"
using json = nlohmann::json;

namespace fs = std::filesystem;

#include "Tissue.hpp"

class VtkWriter : public Tissue {
  public:
    std::string write_folder;
    //std::string write_name;
    int m;
    VtkWriter(int _m, std::string _write_folder):m(_m), write_folder(_write_folder){

    }
    // VtkWriter(std::string _write_name): write_name(_write_name){}
    void write_apical_or_basal_polygonal_surface(std::string surf_type)
    {
      std::string write_name = write_folder + surf_type + "_" + std::to_string(m) + ".vtk";
      std::ofstream write_out; write_out.open(write_name);
      if (write_out.is_open())
      {
        write_out << "# vtk DataFile Version 1.0\n2D Unstructured Grid of Polygons\nASCII\n\nDATASET POLYDATA\n";
        write_out << "POINTS " + std::to_string(ListVertex.size()) + " float\n" ;
        if(surf_type=="apical") for(auto v : ListVertex){write_out << v.pos_a[0] << " " << v.pos_a[1] << " " << v.pos_a[2] << std::endl;}
        else if(surf_type=="basal") for(auto v : ListVertex){write_out << v.pos_b[0] << " " << v.pos_b[1] << " " << v.pos_b[2] << std::endl;}
        else {std::cout << " ERROR:  couldn't write the polygonal mesh to vtk. the surface type not recognized!" << std::endl;}
        int n_tot_elements =0;
        for(auto c : ListCell) n_tot_elements += c.Vlist.size()+1;
        write_out << "\nPOLYGONS "+ std::to_string(ListCell.size()) + " " + std::to_string(n_tot_elements);
        for(auto c : ListCell)
        {
          write_out << "\n"+std::to_string(c.Vlist.size());
          for(auto v : c.Vlist) write_out << " " + std::to_string(v->id) ;
        }
        write_out << "\n\nCELL_DATA "+ std::to_string(ListCell.size())+"\n";
        write_out << "SCALARS neighbors float 1\n";
        write_out << "LOOKUP_TABLE default\n";
        for(auto c : ListCell) write_out << std::to_string(c.Vlist.size()) + " ";
        write_out << std::endl;
      }
      write_out.close();
    }

    void write_apical_or_basal_triangulated_surface(std::string surf_type)
    {
      std::string write_name = write_folder + surf_type + "tri_" + std::to_string(m) + ".vtk";
      std::ofstream write_out; write_out.open(write_name);
      if (write_out.is_open())
      {
        write_out << "# vtk DataFile Version 1.0\n2D Unstructured Grid of Linear Triangles\nASCII\n\nDATASET POLYDATA\n";
        write_out << "POINTS " + std::to_string(ListVertex.size()+ListCell.size()) + " float\n" ;
        if(surf_type=="apical")
        {
          update_apical_geometric_quantities();
          for(auto v : ListVertex) {write_out << v.pos_a[0] << " " << v.pos_a[1] << " " << v.pos_a[2] << std::endl;}
          for(auto c : ListCell) {write_out << c.cent_a[0] << " " << c.cent_a[1] << " " << c.cent_a[2]<<std::endl;}
        }
        else if(surf_type=="basal")
        {
          update_basal_geometric_quantities();
          for(auto v : ListVertex) {write_out << v.pos_b[0] << " " << v.pos_b[1] << " " << v.pos_b[2] << std::endl;}
          for(auto c : ListCell) {write_out << c.cent_b[0] << " " << c.cent_b[1] << " " << c.cent_b[2]<<std::endl;}
        } else {std::cout << " ERROR:  couldn't write the polygonal mesh to vtk. the surface type not recognized!" << std::endl;}

        int n_triangle=0;
        for(auto c : ListCell) n_triangle += c.Vlist.size();
        std::vector<std::tuple<int,int,int>> triangles; triangles.reserve(n_triangle);
        for(auto c : ListCell){
          for(auto e : c.Elist) triangles.emplace_back(std::make_tuple(ListVertex.size()+c.id, e->tail->id, e->head->id));
        };
        if(triangles.size() != n_triangle) std::cout << "number of triangles did not match!!!" << n_triangle <<" "<< triangles.size()<< std::endl;

        write_out << "\nPOLYGONS "+ std::to_string(n_triangle) + " " + std::to_string(4*n_triangle) + "\n";
        for(auto tri : triangles)
        {
          int a, b, c;
          std::tie(a,b,c) = tri;
          write_out << "3 " << a << " " << b << " " << c << std::endl;
        }
        write_out << "\n\nCELL_DATA "+ std::to_string(n_triangle)+"\n";
        write_out << "SCALARS neighbors float 1\n";
        write_out << "LOOKUP_TABLE default\n";
        for(int i_tri=0; i_tri<n_triangle; i_tri++) write_out << 3 << " ";
        write_out << std::endl;
      }
      write_out.close();
    }

    void write_full_3D_triangulation()
    {
      std::string write_name = write_folder + "3Dtri_" + std::to_string(m) + ".vtk";
      std::ofstream write_out; write_out.open(write_name);
      update_apical_geometric_quantities();
      update_basal_geometric_quantities();
      update_lateral_geometric_quantities();
      int Nv = ListVertex.size();
      int Nc = ListCell.size();
      int Ne = ListEdge.size();
      std::cout << Nv << " " << Nc << " " << Ne << " " << ListEdge.size() << std::endl;
      int n_triangle= 6 * Ne;
      //for(auto c : ListCell) n_triangle += 2 * c.Vlist.size();
      std::vector<std::tuple<int,int,int>> triangles; triangles.reserve(n_triangle);
      std::vector<int> p_class; p_class.reserve(n_triangle);//lateral=4, apical=basal=c.Vlist.size()
      for(int i=0; i<Nv; i++) ListVertex[i].id=i;
      for(int i=0; i<Nc; i++) ListCell[i].id=i;
      int e_counter=0;
      for(int i=0; i<ListEdge.size(); i++) {
	ListEdge[i].id=i;
        //if(ListEdge[i].id<ListEdge[i].conj->id) {ListEdge[i].id=e_counter; ListEdge[i].conj->id=Ne+e_counter; e_counter++;}
      }
      if (write_out.is_open())
      {
        write_out << "# vtk DataFile Version 1.0\n2D Unstructured Grid of Linear Triangles\nASCII\n\nDATASET POLYDATA\n";
        //write_out << "POINTS " + std::to_string(2*Nv+2*Nc) + " float\n" ;
        write_out << "POINTS " + std::to_string(2*Nv+2*Nc+Ne) + " float\n" ;
        for(auto v : ListVertex) {
          write_out << v.pos_a[0] << " " << v.pos_a[1] << " " << v.pos_a[2] << std::endl;
          write_out << v.pos_b[0] << " " << v.pos_b[1] << " " << v.pos_b[2] << std::endl;
        }
        for(auto c : ListCell) {
          write_out << c.cent_a[0] << " " << c.cent_a[1] << " " << c.cent_a[2]<<std::endl;
          write_out << c.cent_b[0] << " " << c.cent_b[1] << " " << c.cent_b[2]<<std::endl;
        }
        for(auto e : ListEdge) /*{if(e.id < e.conj->id)*/{write_out << e.cent_l[0] << " " << e.cent_l[1] << " " << e.cent_l[2] << std::endl;}// avoid repeating conj lateral triangles

        for(auto c : ListCell){ // fill in the triangle coordinates
          int apical_center = 2*Nv+2*c.id;
          int basal_center = 2*Nv+2*c.id+1;
          for(auto e : c.Elist)
          {
            //apical_triangle
            triangles.emplace_back(std::make_tuple(apical_center, 2*e->tail->id, 2*e->head->id));
            p_class.emplace_back(c.Vlist.size());
            //basal triangle
            triangles.emplace_back(std::make_tuple(basal_center, 2*e->tail->id+1, 2*e->head->id+1));
            p_class.emplace_back(c.Vlist.size());
            //4 lateral triangles
            //if(e->id < e->conj->id) // avoid repeating lateral triangles
            //{
              int e_lat_cent = 2*Nv + 2*Nc + e->id;
	    //  std::cout << "*******  e_id:    " << e->id << std::endl;
              triangles.emplace_back(std::make_tuple(e_lat_cent, 2*e->tail->id, 2*e->head->id)); p_class.emplace_back(4);
              triangles.emplace_back(std::make_tuple(e_lat_cent, 2*e->tail->id+1, 2*e->head->id+1)); p_class.emplace_back(4);
              triangles.emplace_back(std::make_tuple(e_lat_cent, 2*e->tail->id, 2*e->tail->id+1)); p_class.emplace_back(4);
              triangles.emplace_back(std::make_tuple(e_lat_cent, 2*e->head->id, 2*e->head->id+1)); p_class.emplace_back(4);
            //}
          }
        };
        if(triangles.size() != n_triangle) std::cout << "number of triangles did not match!!!" << n_triangle <<" "<< triangles.size()<< std::endl;

        write_out << "\nPOLYGONS "+ std::to_string(n_triangle) + " " + std::to_string(4*n_triangle) + "\n";
        for(auto tri : triangles)
        {
          int a, b, c;
          std::tie(a,b,c) = tri;
          write_out << "3 " << a << " " << b << " " << c << std::endl;
        }
        write_out << "\nCELL_DATA "+ std::to_string(n_triangle)+"\n";
        write_out << "SCALARS neighbors float 1\n";
        write_out << "LOOKUP_TABLE default\n";
        for(int i_tri=0; i_tri<n_triangle; i_tri++) write_out << p_class[i_tri] << " ";
        write_out << std::endl;
      }
      write_out.close();
    }

    void write_polarity_field()
    {
      update_apical_geometric_quantities();
      std::string write_name = write_folder + "polarity_" + std::to_string(m) + ".vtk";
      std::ofstream write_out; write_out.open(write_name);
      if (write_out.is_open())
      {
        write_out << "# vtk DataFile Version 1.0\npolarity field\nASCII\n\nDATASET POLYDATA\n";
        write_out << "POINTS " + std::to_string(ListCell.size()) + " float\n" ;
        for(auto c : ListCell) {write_out << c.cent_a[0] << " " << c.cent_a[1] << " " << c.cent_a[2]<<std::endl;}

        write_out << "\nPOINT_DATA "+ std::to_string(ListCell.size()) + "\n";
        write_out << "SCALARS F_mag float\nLOOKUP_TABLE default\n";
        for(auto c : ListCell) {write_out << c.mag_F_a << " ";}

        write_out << "\nVECTORS polarity float\n";
        for(auto c : ListCell) {write_out << c.polarity[0] << " " << c.polarity[1] << " " << c.polarity[2]<<std::endl;}
      }
      write_out.close();
    }

};


int main(int argc, char *argv[])
{
  std::string data_folder = std::string(argv[1]);
  int fr_begin = std::stoi(argv[2]);
  int fr_end = std::stoi(argv[3]);
  std::string out_folder = std::string(argv[4]);

  const fs::path out_path{out_folder};// build the directory if does not exist
  fs::create_directory(out_path);

  //Tissue *tissue0 = new Tissue(); // heap memory allocation

  for(int fr=fr_begin; fr<=fr_end; fr++)
  {
    //Tissue *tissue = new Tissue(); // "placement" of the new tissue into previously allocated heap memory
    std::unique_ptr<Tissue> tissue(new Tissue()); // a unique ptr

    std::string jf = data_folder + "frame_"+std::to_string(fr)+".json";

    tissue->build_3D_spherical_tissue_from_json_file(jf);

    std::cout << "Nc:  " << tissue->ListCell.size()  << "  nv:  " << tissue->ListVertex.size()<< std::endl;
    //
    VtkWriter vtkWriter(fr, out_folder);
    //
    vtkWriter.ListCell = tissue->ListCell;
    vtkWriter.ListEdge = tissue->ListEdge;
    vtkWriter.ListVertex = tissue->ListVertex;
    //

    vtkWriter.write_apical_or_basal_polygonal_surface("apical");
    vtkWriter.write_apical_or_basal_triangulated_surface("apical");
    vtkWriter.write_polarity_field();

    //vtkWriter.write_full_3D_triangulation();
    //delete tissue;
  }


  return 0;
}
