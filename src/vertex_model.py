import glob,os,subprocess,shutil
import pyvista as pv

class vertex_model:
    class Params:
        def __init__(self):
            self.time_step    = 0.0002
            # ...number of time steps after which is frame will be written to output
            self.ninfo    = 100000
            # ...the simulation will continue for nframes frames.
            # ...therefore, total steps = nframes * ninfo
            self.nframes    = 3000
            # ...if set, after this many frames the noise on polarity turns off
            # ...if commented out, it gets default value: std::numeric_limits<unsigned>::max()
            self.noiseoff    = 2800
            # ...by default vertices are constrained on a spherical surface.
            # ...to relax this constrain and allow for radial deformations, un-comment:
            #self.spherical_constrain       = 0

            # ==============================  model parameters     ===================================
            # ...number of cells
            self.Nc       = 200
            # ...cell area stiffness
            self.Kc       = 1.0
            # ...cell preferred area
            self.A0c       = 1.0
            # ...cell bond tension
            self.bond_Tension   = 0.1
            # ...cell perimeter elasticity
            self.perim_Elasticity   = 0.0
            # ...cut-off and opening lengths of bonds for T1
            self.bond_T1_cutoff   = 0.04
            self.bond_T1_opening   = 0.045

            # ...vertex friction coefficient with external environment
            self.xi       = 1.0

            # ================== initial patterns of cell polarity (default=randomized) =============
            self.P0_axissymmetric       = 1
            # ...if P0_axissymmetric is set to 1, the code uses the parameters:
            self.P0_a0_by_PI       = 0
            self.P0_b0_by_PI       = 0

            # ================== changes initial cellular network, for ensemble averages =============
            self.run_id       = 1

            # ...if you do not want to initialize with a predefined surface, set the following to zero
            # ... if set to zero, a random tissue will be generated (default=1)
            self.read_initsurface       = 1

            #  ============================  parameters for active terms =============================
            # ...cell traction force magnitude
            self.F_mag   = 0.02
            # ...rate of cell polarity alignment with neighbors
            self.P_gamma   = 0.005
            # ...rate of cell polarity alignment with cell velocity
            self.P_nu   = 0.0
            # ...strength of rotational noise in the polarity
            self.P_noise   = 0.001
            # ...polarity norm constrains
            self.elastic_polarity       = 1
            # ...if elastic_polarity is set to 1, instead of a hard constrain |p|=1, the molecular field
            # ...in polarity p dynamics will include the term:  P_A * (1. |p|^2 ) * p
            self.P_A       = 0.001

            #  ======================  setting seed and initial noise used for debugging =============
            self.noise   = 0.01
            self.seed   = 1
    class AnalyzeParams:
        def __init__(self):
            #no_write       = 1
            # ...should I write the data to vtk files?
            write_to_vtk       = 1
            # ...should I analyze the cell elongation patterns?
            analyze_elongation       = 1
            # ...should I decompose apical surface to vector spherical harmonics modes
            analyze_apical_surface_VSH       = 1
            # ...should I decompose cell polarity field to vector spherical harmonics modes
            analyze_cell_polarity_VSH       = 1
            # ...should I analyze the coarse-grained curvature tensor on defined patches
            compute_curvature_tensor       = 1
            # ...should I analyze tissue rotation, angular velocity and residual from solid body?
            analyze_rotation       = 1
            # ...should I align such that rotation axis points to z-direction?
            align_before_writing       = 1
            # ...should I analyze data in co-rotating fram?
            analysis_in_corotating_frame       = 1
            # ... what kind of data should be written to vtk?
            write_apical_polygonal_surface       = 1
            write_basal_polygonal_surface       = 0
            write_apical_triangulated_surface       = 1
            write_basal_triangulated_surface       = 0
            write_full_3D_triangulated_surface       = 0
            write_polarity_field       = 1
            write_nematic_field       = 1

            # ==============================  model parameters     ===================================
            # ...the maximum l mode for vector spherical harmonics mode decomposition?
            Lmax       = 4
            # ... first frame number inside data_dir to be analyzed
            first_frame       = 1
            # ... last frame number inside data_dir to be analyzed
            last_frame       = 100
        
    def __init__(self,path):
        self.params = self.Params()
        self.analyzeParams = self.AnalyzeParams()
        self.multiblock = None
        self.multiblock2 = None
        self.repo_path=path


    def RunSimulation(self, nThreads=1):
        SimParams=self.params
        with open(os.path.join(self.repo_path, 'cpp', 'VertexModel.dat'), 'w') as f:
            for attr in dir(SimParams):
                if not attr.startswith('__'):
                    f.write(f"{attr}={getattr(SimParams, attr)}\n")
        output_dir = os.path.join(self.repo_path, 'cpp', 'VertexModelOutput')

        # Check if the directory exists
        if os.path.exists(output_dir):
            # Clear the directory
            for file in os.listdir(output_dir):
                file_path = os.path.join(output_dir, file)
                try:
                    if os.path.isfile(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                except Exception as e:
                    print(f"Failed to delete {file_path}. Reason: {e}")
        else:
            # Create the directory
            os.makedirs(output_dir)

        os.popen('[ -d'+self.repo_path+'/cpp/VertexModelOutput/ ] && rm -rf '+self.repo_path +
                '/cpp/VertexModelOutput/* || mkdir -p '+self.repo_path+'/cpp/VertexModelOutput/')
        cmd = [
            "/bin/bash", "-c",
            "../vertex_model3d_monolayer/organoid ../VertexModel.dat"
        ]
        try:
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, cwd=output_dir)
            stdout, stderr = p.communicate()
            # Check for errors
            if p.returncode != 0:
                print(f"Command failed with error code {p.returncode}")
                if stderr:
                    print(f"Error message: {stderr.decode('utf-8')}")
            else:
                # Split the output into lines and print the last two lines
                output_lines = stdout.decode('utf-8').splitlines()
                print(output_lines[-2])
                print(output_lines[-1])

        except Exception as e:
            print(f"Failed to run command: {e}")
        return 
    
    def load_sim_data(self):
        output_dir = os.path.join(self.repo_path, 'cpp', 'VertexModelOutput')
        file_pattern = os.path.join(output_dir,"./outputAnalyze/apicaltri*.vtk")
        vtk_files = glob.glob(file_pattern)
        file_pattern2 = os.path.join(output_dir,"./outputAnalyze/polarity*.vtk")
        vtk_files2 = glob.glob(file_pattern2)

        # Create a MultiBlock dataset
        self.multiblock = pv.MultiBlock()
        self.multiblock2 = pv.MultiBlock()

        # Load each file into the MultiBlock dataset
        for i, vtk_file in enumerate(vtk_files):
            mesh = pv.read(vtk_file)
            self.multiblock[str(i)] = mesh
        for i, vtk_file in enumerate(vtk_files2):
            mesh = pv.read(vtk_file)
            self.multiblock2[str(i)] = mesh        
        return 
    
    def AnalyzeData(self, nThreads=1):
        AnaParams=self.analyzeParams
        with open(os.path.join(self.repo_path, 'cpp', 'VertexModelAnalyze.dat'), 'w') as f:
            for attr in dir(AnaParams):
                if not attr.startswith('__'):
                    f.write(f"{attr}={getattr(AnaParams, attr)}\n")
        output_dir = os.path.join(self.repo_path, 'cpp', 'VertexModelOutput')
    
        # Check if the directory exists
        if not os.path.exists(output_dir+'/output/frame_0.json'):
            print("Please run the simulation first")
            return

        cmd = [
            "/bin/bash", "-c",
            "../vertex_model3d_monolayer/accessories/analyze ../VertexModelAnalyze.dat -d ./output -o ./outputAnalyze"
        ]
        try:
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, cwd=output_dir)
            stdout, stderr = p.communicate()
            # Check for errors
            if p.returncode != 0:
                print(f"Command failed with error code {p.returncode}")
                if stderr:
                    print(f"Error message: {stderr.decode('utf-8')}")
            else:
                # Split the output into lines and print the last two lines
                output_lines = stdout.decode('utf-8').splitlines()
                print(output_lines[-2])
                print(output_lines[-1])

        except Exception as e:
            print(f"Failed to run command: {e}")
        
        self.load_sim_data()
        return

    
    def VizualizeIteration(self, iteration="0",edges=False):
        if self.multiblock is None:
            print("Please analyze the data first")
            return

        # Ensure the iteration index is within bounds
        if iteration >= len(self.multiblock) or iteration < 0:
            print("Invalid iteration index")
            return
        
        iteration = str(iteration)
        # Get the array names of the specified block
        #array_names = self.multiblock[iteration].array_names

        # Print array names for debugging or analysis
        #print(f"All arrays: {array_names}")

        # Visualize the specified iteration using PyVista plotter
        plotter = pv.Plotter()
        # Add glyphs on top of the surface
        glyph = self.multiblock2[iteration].glyph(orient='polarity', factor=30,tolerance=0.001)
        plotter.add_mesh(glyph, lighting=False, render_points_as_spheres=True, point_size=10, color='firebrick')
        # Add the mesh as surface with edges
        plotter.add_mesh(self.multiblock[iteration],show_edges=edges, color='lightgrey')
        plotter.add_text("Time:"+iteration, color="black",position='upper_edge',font_size=14)
        plotter.show()

    def VizualizeAnimate(self,gif_name,edges=False):
        if self.multiblock is None:
            print("Please analyze the data first")
            return
        
        plotter = pv.Plotter()
        # Open a gif
        plotter.open_gif(gif_name)

        for iteration in range(len(self.multiblock)):
            iter=str(iteration)
            glyph = self.multiblock2[iter].glyph(orient='polarity', factor=30,tolerance=0.001)
            plotter.add_mesh(glyph, lighting=False, render_points_as_spheres=True, point_size=10, color='firebrick')
            plotter.add_mesh(self.multiblock[iter],show_edges=edges, color='lightgrey')
            plotter.add_text("Time:"+iter, color="black",position='upper_edge',font_size=14)
            plotter.write_frame()
            plotter.clear()
        plotter.close()  
        return 
