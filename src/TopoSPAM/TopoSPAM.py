import numpy as np
import os
import subprocess
import time
import vtk
import matplotlib.pyplot as plt
from vtk.util.numpy_support import vtk_to_numpy
from vtk.numpy_interface import dataset_adapter as dsa
from IPython.display import HTML, display
from ipywidgets import widgets, interact_manual, FloatProgress
from matplotlib.animation import FuncAnimation
from matplotlib.offsetbox import AnchoredText
import linecache
from TopoSPAM.mesh_methods import *
repo_path = None


def set_repo_path(path):
    """
    Set the repository path and check if there's a cpp folder with a Makefile inside.

    Args:
    - path (str): The path to the repository.

    Returns:
    - bool: True if there's a cpp folder with a makefile inside, False otherwise.
    """
    global repo_path
    # Ensure the provided path exists
    if not os.path.exists(path):
        print(f"Error: The path '{path}' does not exist.")
        return

    # Check for the bin folder
    bin_path = os.path.join(path, 'cpp')
    if not os.path.exists(bin_path):
        print(f"Error: The bin folder does not exist in '{path}'.")
        return

    # Check for the Makefile inside the bin folder
    makefile_path = os.path.join(bin_path, 'makefile')
    if not os.path.exists(makefile_path):
        print(f"Error: makefile not found in the bin folder of '{path}'.")
        return

    print(f"Success: The path '{path}' contains the TopoSPAM repository.")
    repo_path = path
    return


class ActiveFluid2D:
    class Params:
        def __init__(self):
            self.Gd_Sz = 41
            self.Ks = 1.0
            self.Kb = 1.0
            self.dmu = 0.0
            self.nu = 0.0
            self.zeta = 1.0
            self.lamda = 1.0
            self.eta = 1.0
            self.gama = 1.0
            self.tf = 1
            self.dt = 1e-2
            self.wrat = 1
            self.Vtol = 1e-2
            self.absttol = 1e-3
            self.relttol = 1e-3
        def get(self):
            return np.array([self.Gd_Sz, self.Ks, self.Kb, self.dmu, self.nu, self.zeta, self.lamda, self.eta, self.gama, self.tf, self.dt, self.wrat, self.Vtol, self.absttol, self.relttol])

    def __init__(self):
        self.params = self.Params()
    
    def RunSimulation(self,nCores=1):
        global repo_path
        SimParams=self.params.get()
        np.savetxt(repo_path+'/cpp/active2d.csv', SimParams, delimiter=' ')
        output_dir = os.path.join(repo_path, 'cpp', 'Active2dOutput')

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

        os.popen('[ -d'+repo_path+'/cpp/Active2dOutput/ ] && rm -rf '+repo_path +
                '/cpp/Active2dOutput/* || mkdir -p '+repo_path+'/cpp/Active2dOutput/')
        cmd = [
            "/bin/bash", "-c",
            "mpirun -np " +
            str(nCores)+" ../Active2d ../active2d.csv"
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
        return int(output_lines[-1].split()[-1])


    def RunSimulationWithProgress(self, nCores=1):
        global repo_path
        SimParams=self.params
        np.savetxt(repo_path+'/cpp/active2d.csv', SimParams, delimiter=' ')
        output_dir = os.path.join(repo_path, 'cpp', 'Active2dOutput')
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
        
        cmd = [
            "/bin/bash", "-c",
            "mpirun -np " +
            str(nCores)+" ../Active2d ../active2d.csv"
        ]    
        p = subprocess.Popen(cmd,
                            stdout=subprocess.PIPE, cwd=repo_path+'/cpp/Active2dOutput', shell=True)
        f = FloatProgress(value=0.0, min=0.0, max=self.params.tf)  # instantiate the bar
        print("Simulation Progress (Based on current time value)")
        display(f)  # display the bar
        ctr = 0
        outp = "0"
        for line in iter(p.stdout.readline, b''):
            outp = line.rstrip().decode('utf-8')
            if (outp.split()[0] == "Total" or outp.split()[0] == "The"):
                print(outp)
            else:
                f.value = float(line.rstrip().decode('utf-8'))
            ctr += 1
        print("Simulation finished (reached t="+str(self.params.tf)+")")
        return int(outp.split()[-1])


    def getSimData(self,iteration):
        global repo_path
        PropNames = ['00-Polarization', '01-Velocity', '02-Vorticity_0_0', '02-Vorticity_0_1', '02-Vorticity_1_0', '02-Vorticity_1_1', '03-ExternalForce', '04-Pressure', '05-StrainRate_0_0', '05-StrainRate_0_1', '05-StrainRate_1_0', '05-StrainRate_1_1', '06-Stress_0_0',
                    '06-Stress_0_1', '06-Stress_1_0', '06-Stress_1_1', '07-MolecularField', '08-DPOL', '09-DV', '10-VRHS', '11-f1', '12-f2', '13-f3', '14-f4', '15-f5', '16-f6', '17-V_T', '18-DIV', '19-DELMU', '20-HPB', '21-FrankEnergy', '22-R', 'SubsetNumber', 'domain']
        if not os.path.exists(repo_path+'/cpp/Active2dOutput/Polar_'+str(iteration)+'.pvtp'):
            print("Iteration does not exists")
            return None
        reader = vtk.vtkXMLPPolyDataReader()
        reader.SetFileName(
            repo_path+'/cpp/Active2dOutput/Polar_'+str(iteration)+'.pvtp')
        reader.Update()
        nodes = vtk_to_numpy(dsa.WrapDataObject(reader.GetOutput()).Points)
        t = float(linecache.getline(repo_path+'/cpp/Active2dOutput/Polar_' +
                str(iteration)+'.pvtp', 5, module_globals=None))
        Pol = vtk_to_numpy(
            reader.GetOutput().GetPointData().GetArray(PropNames[0]))
        Vel = vtk_to_numpy(
            reader.GetOutput().GetPointData().GetArray(PropNames[1]))
        FE = vtk_to_numpy(
            reader.GetOutput().GetPointData().GetArray('21-FrankEnergy'))
        x, y = nodes[:, 0], nodes[:, 1]
        return x, y, Pol, Vel, FE, t


    def VizualizeIteration(self,iteration=0):
        x, y, Pol, Vel, FE, t = self.getSimData(iteration)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        q1 = ax1.quiver(x, y, Pol[:, 0], Pol[:, 1], FE, cmap=plt.cm.jet)
        fig.colorbar(q1, cmap=plt.cm.jet, label='Frank Energy Density', ax=ax1)
        q2 = ax2.quiver(x, y, Vel[:, 0], Vel[:, 1], np.linalg.norm(
            Vel, axis=1), cmap=plt.cm.viridis)
        ax1.set_title("Polarity Vectors")
        ax2.set_title("Velocity Vectors")
        fig.colorbar(q2, cmap=plt.cm.viridis, label='Velocity Magnitude', ax=ax2)
        fig.suptitle('Time:'+str(round(t, 2)))
        # ax1.annotate('Time:'+str(round(t,2)), xy=(-0.5, 11), xycoords='data', annotation_clip=False,size=15)
        plt.show()
        return plt


    def VizualizeAnimate(self,finalStep, jump=1):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        x, y, Pol, Vel, FE, t = self.getSimData(finalStep)
        q1 = ax1.quiver(x, y, Pol[:, 0], Pol[:, 1], FE,
                        cmap=plt.cm.jet, width=0.005, headlength=4)
        q2 = ax2.quiver(x, y, Vel[:, 0], Vel[:, 1], np.linalg.norm(
            Vel, axis=1), cmap=plt.cm.viridis, width=0.005, headlength=4)
        ax1.set_title("Polarity Vectors")
        ax2.set_title("Velocity Vectors")
        c1 = fig.colorbar(q1, cmap=plt.cm.jet,
                        label='Frank Energy Density', ax=ax1)
        c2 = fig.colorbar(q2, cmap=plt.cm.viridis,
                        label='Velocity Magnitude', ax=ax2)
        tittle = fig.suptitle("Time:"+str(round(t, 2)))
        # anno=ax1.annotate('Time:'+str(t), xy=(-0.5, 10), xycoords='data', annotation_clip=False,size=15)

        def animate(i):
            x, y, Pol, Vel, FE, t = self.getSimData(i)
            q1.set_UVC(Pol[:, 0], Pol[:, 1], FE)
            q2.set_UVC(Vel[:, 0], Vel[:, 1], np.linalg.norm(Vel, axis=1))
            tittle.set_text("Time:"+str(round(t, 2)))
        anim = FuncAnimation(fig, animate, np.arange(
            1, finalStep, jump), interval=150)
        display(HTML(anim.to_html5_video()))
        plt.close()


class SpringLattice:
    class Params:
        def __init__(self):
            self.thickness = 1.0
            self.lambda_tensor_diagonal = [1.5, 1/1.5, 1]
            self.nematic_coordinates = "polar"
            self.mesh_geometry = "circle"
            self.balls = None
            self.springs = None
            self.init_positions = None

    def __init__(self):
        self.params = self.Params()

    def load_mesh(self, thick_mesh_file=None):  # find a better way to refer to the mesh file

        global repo_path
        if self.params.mesh_geometry == "square":
            if thick_mesh_file is None:
                thick_mesh_file = os.path.join(
                    repo_path, "ExampleNotebooks/meshes/square.pkl")
            # replace the line below to either create a mesh or read a mesh
            [balls_df, springs_df] = pickle.load(open(thick_mesh_file, 'rb'))
            # setting thickness
            balls_df["z"][len(balls_df)//2:] = self.params.thickness
            # offsetting the mesh so that it starts from y = 0
            balls_df["y"] = balls_df["y"] - balls_df["y"].min()
            # update springs
            springs_df = update_springs(springs_df, balls_df[['x', 'y', 'z']])

        elif self.params.mesh_geometry == "circle":
            if thick_mesh_file is None:
                thick_mesh_file = os.path.join(
                    repo_path, "ExampleNotebooks/meshes/circle.pkl")
            [balls_df, springs_df] = pickle.load(open(thick_mesh_file, 'rb'))
            # setting thickness
            balls_df["z"][len(balls_df)//2:] = self.params.thickness
            # if any of the balls have x = 0 then offset by small amount
            ind = balls_df[balls_df["x"] == 0].index
            balls_df.loc[ind, "x"] = 0.00001
            # update springs
            springs_df = update_springs(springs_df, balls_df[['x', 'y', 'z']])
            # compute the polar coordinates of the balls
            balls_df["r"] = np.sqrt(balls_df["x"]**2 + balls_df["y"]**2)
            balls_df["theta"] = np.arctan2(balls_df["y"], balls_df["x"])

        # springs_df = update_springs(springs_df, balls_df[['x', 'y', 'z']]) #, compute_lo = True
        self.params.balls = balls_df
        self.params.springs = springs_df
        self.params.init_positions = balls_df[["x", "y", "z"]]

        return


    def load_strain_pattern(self):
        balls_df = self.params.balls
        springs_df = self.params.springs
        lambda_tensor_diagonal = self.params.lambda_tensor_diagonal
        nematic_coordinates = self.params.nematic_coordinates
        # compute the rest length of each spring
        springs_df['l0'] = springs_df.apply(get_rest_length, axis=1, lambda_tensor_diagonal=lambda_tensor_diagonal,
                                            nematic_coordinates=nematic_coordinates, balls_df=balls_df)
        springs_df['l0_target'] = springs_df['l0']
        springs_df['l1_initial'] = springs_df['l1']
        self.params.springs = springs_df

        return


    def get_lambda_tensor(self,pos_vector, lambda_tensor_diagonal, nematic_coordinates="polar", ):
        """
        This function computes the director field for a given point in space
        and the magnitude of the lambda coefficient at that point
        and it returns the lambda tensor at that point in cartesian coordinates
        This can be made more general by allowing the user to input the director field
        """
        [lambda1, lambda2, lambda3] = lambda_tensor_diagonal
        if nematic_coordinates == "polar":
            theta = np.arctan2(pos_vector[1], pos_vector[0])
            r = np.sqrt(pos_vector[0]**2 + pos_vector[1]**2)
            director1 = np.array([pos_vector[0]/r, pos_vector[1]/r, 0])
            director2 = np.array([-pos_vector[1]/r, pos_vector[0]/r, 0])
            director3 = np.array([0, 0, 1])
            lambda_tensor = lambda1(r, theta)*np.tensordot(director1, director1, axes=0) \
                + lambda2(r, theta)*np.tensordot(director2, director2, axes=0) \
                + lambda3(r, theta)*np.tensordot(director3, director3, axes=0)

        elif nematic_coordinates == "cartesian":
            director1 = [1, 0, 0]
            director2 = [0, 1, 0]
            director3 = [0, 0, 1]
            [x, y, z] = pos_vector
            lambda_tensor = lambda1(x, y)*np.tensordot(director1, director1, axes=0) \
                + lambda2(x, y)*np.tensordot(director2, director2, axes=0) \
                + lambda3(x, y)*np.tensordot(director3, director3, axes=0)

        return (lambda_tensor)


    def get_rest_length(self,row, lambda_tensor_diagonal=None, nematic_coordinates="polar", balls_df=None):
        # this function will take a row of the springs_df and return the rest length

        # get lambda tensor at one endpoint in cartesian coordinates
        lambda_alpha = self.get_lambda_tensor(
            balls_df.iloc[row["ball1"]][["x", "y", "z"]].values,
            lambda_tensor_diagonal,
            nematic_coordinates=nematic_coordinates
        )

        # get lambda tensor at other endpoint
        lambda_beta = self.get_lambda_tensor(
            balls_df.iloc[row["ball2"]][["x", "y", "z"]].values,
            lambda_tensor_diagonal,
            nematic_coordinates=nematic_coordinates
        )
        # average the lambda tensors
        lambda_avg = (lambda_alpha + lambda_beta)/2
        # Multiply the average lambda tensor with the original spring vector
        # get the rest length by taking the norm of the above vector
        spring_vector = np.array(
            [row['x1'] - row['x2'], row['y1'] - row['y2'], row['z1'] - row['z2']])
        new_rest_length = np.linalg.norm(
            np.matmul(lambda_avg, spring_vector))

        return (new_rest_length)


    def add_noise(self, noise=0.1):
        # add noise to the initial positions
        self.params.balls[["x", "y", "z"]] = self.params.balls[["x", "y", "z"]] + \
            np.abs(noise*np.random.randn(*self.params.balls[["x", "y", "z"]].shape))
        self.params.springs = update_springs(
            self.params.springs, self.params.balls[['x', 'y', 'z']])
        return


    def visualize(self, ax=None, fig=None, mode="discrete", x='x', y='y',
                title="Strain Pattern", color_min=0.7, color_max=1.3, state="final",
                tick_fontsize=10, label_fontsize=16, cbar_name=r'$\frac{l_{rest}}{l_{init}}$', title_fontsize=20,
                xlim_min=-1.2, xlim_max=1.2, ylim_min=-1.2, ylim_max=1.2, plot_only_top=True,
                ):
        Params=self.params

        if mode == "continuous":
            pass

        elif mode == "discrete":

            balls_df = Params.balls
            springs_df = Params.springs

            if state == "initial":
                balls_df[['x', 'y', 'z']] = Params.init_positions.values
            elif state == "final":
                balls_df[['x', 'y', 'z']] = Params.final_positions.values
            springs_df = update_springs(springs_df, balls_df[['x', 'y', 'z']])

            if ax is None:
                fig, ax = plt.subplots()

            fig, ax = plot_shell_on_given_ax(balls_df, springs_df, x=x, y=y, ax=ax, fig=fig,
                                            xlim_min=xlim_min, xlim_max=xlim_max, ylim_min=ylim_min, ylim_max=ylim_max,
                                            title=title, color_min=color_min, color_max=color_max, cmap="RdBu_r", return_ax=True,
                                            show_cbar=True, tick_fontsize=tick_fontsize, cbar_name=cbar_name, label_fontsize=label_fontsize,
                                            plot_only_top=plot_only_top
                                            )
            # title
            ax.set_title(title, fontsize=title_fontsize)
            return fig, ax

        # return ax


    def RunSimulation(self,dt=0.01, tol=1e-6, csv_t_save=500):
        global repo_path
        dirname = repo_path+'/cpp/SpringLatticeOutput/'
        bin_dir = repo_path+'/cpp/'

        os.makedirs(dirname, exist_ok=True)
        os.makedirs(dirname+'runfiles/', exist_ok=True)
        os.makedirs(dirname+'sim_output/', exist_ok=True)

        [balls_df, springs_df] = initialize_cpp_simulation(
            self.params.balls, self.params.springs, dt=dt, csv_t_save=csv_t_save, tol=tol, path=dirname)

        # filelist = ['Makefile', 'SpringLattice.cpp']
        # for file in filelist:
        #    shutil.copy(bin_dir + file, dirname)

        # running the simulation
        print('$$$$$$$ Running openfpm $$$$$$$')
        os.system("cd " + bin_dir +
                " && source /usr/local/openfpm/source/openfpm_vars && make SpringLattice")
        os.system("cd " + dirname +
                " && source /usr/local/openfpm/source/openfpm_vars && ../SpringLattice")
        print('$$$$ Exit OpenFPM $$$$')

        # access the output files
        final_balls_df = pd.read_csv(dirname + 'files/final_0_0.csv')
        self.params.final_positions = pd.DataFrame(
            final_balls_df[['x[0]', 'x[1]', 'x[2]']].values, columns=['x', 'y', 'z'])
        self.params.balls[['x', 'y', 'z']] = self.params.final_positions
        self.params.springs = update_springs(
            self.params.springs, self.params.balls[['x', 'y', 'z']])

        return 


class VertexModel:
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
        
    def __init__(self):
        self.params = self.Params()
        self.analyzeParams = self.AnalyzeParams()


    def RunSimulation(self, nThreads=1):
        global repo_path
        SimParams=self.params
        with open(os.path.join(repo_path, 'cpp', 'VertexModel.dat'), 'w') as f:
            for attr in dir(SimParams):
                if not attr.startswith('__'):
                    f.write(f"{attr}={getattr(SimParams, attr)}\n")
        output_dir = os.path.join(repo_path, 'cpp', 'VertexModelOutput')

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

        os.popen('[ -d'+repo_path+'/cpp/VertexModelOutput/ ] && rm -rf '+repo_path +
                '/cpp/VertexModelOutput/* || mkdir -p '+repo_path+'/cpp/VertexModelOutput/')
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
    
    def AnalyzeData(self, nThreads=1):
        global repo_path
        AnaParams=self.analyzeParams
        with open(os.path.join(repo_path, 'cpp', 'VertexModelAnalyze.dat'), 'w') as f:
            for attr in dir(AnaParams):
                if not attr.startswith('__'):
                    f.write(f"{attr}={getattr(AnaParams, attr)}\n")
        output_dir = os.path.join(repo_path, 'cpp', 'VertexModelOutput')
    
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
        return 

