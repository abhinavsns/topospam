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
from TopoSPAM.methods import *
repo_path = None


def set_repo_path(path):
    """
    Set the repository path and check if there's a bin folder with a Makefile inside.

    Args:
    - path (str): The path to the repository.

    Returns:
    - bool: True if there's a bin folder with a makefile inside, False otherwise.
    """
    global repo_path
    # Ensure the provided path exists
    if not os.path.exists(path):
        print(f"Error: The path '{path}' does not exist.")
        return

    # Check for the bin folder
    bin_path = os.path.join(path, 'bin')
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


class ActiveFluidParameters:
    Gd_Sz = 41
    Ks = 1.0
    Kb = 1.0
    dmu = 0.0
    nu = 0.0
    zeta = 1.0
    lamda = 1.0
    eta = 1.0
    gama = 1.0
    tf = 1
    dt = 1e-2
    wrat = 1
    Vtol = 1e-2
    absttol = 1e-3
    relttol = 1e-3


def RunActiveFluidSimulation(Params, nCores=1):
    global repo_path
    data = np.array([Params.Gd_Sz,
                     Params.Ks,
                     Params.Kb,
                     Params.dmu,
                     Params.nu,
                     Params.zeta,
                     Params.lamda,
                     Params.eta,
                     Params.gama,
                     Params.tf,
                     Params.dt,
                     Params.wrat,
                     Params.Vtol,
                     Params.absttol,
                     Params.relttol
                     ])
    np.savetxt(repo_path+'/bin/active2d.csv', data, delimiter=' ')
    output_dir = os.path.join(repo_path, 'bin', 'Active2dOutput')

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
    os.popen('[ -d'+repo_path+'/bin/Active2dOutput/ ] && rm -rf '+repo_path +
             '/bin/Active2dOutput/* || mkdir -p '+repo_path+'/bin/Active2dOutput/')
    cmd = [
        "/bin/bash", "-c",
        "source /usr/local/openfpm/source/openfpm_vars && mpirun -np " +
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


def RunActiveFluidSimulationWithProgress(Params, nCores=1):
    global repo_path
    data = np.array([Params.Gd_Sz,
                     Params.Ks,
                     Params.Kb,
                     Params.dmu,
                     Params.nu,
                     Params.zeta,
                     Params.lamda,
                     Params.eta,
                     Params.gama,
                     Params.tf,
                     Params.dt,
                     Params.wrat,
                     Params.Vtol,
                     Params.absttol,
                     Params.relttol
                     ])
    np.savetxt(repo_path+'/bin/active2d.csv', data, delimiter=' ')
    output_dir = os.path.join(repo_path, 'bin', 'Active2dOutput')
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
    p = subprocess.Popen('source  /usr/local/openfpm/source/openfpm_vars && ../Active2d ../active2d.csv',
                         stdout=subprocess.PIPE, cwd=repo_path+'/bin/Active2dOutput', shell=True)
    f = FloatProgress(value=0.0, min=0.0, max=Params.tf)  # instantiate the bar
    print("Simulation Progress (Based on current time value)")
    display(f)  # display the bar
    ctr = 0
    for line in iter(p.stdout.readline, b''):
        if (ctr > 5):
            outp = line.rstrip().decode('utf-8')
            if (outp.split()[0] == "Total" or outp.split()[0] == "The"):
                print(outp)
            else:
                f.value = float(line.rstrip().decode('utf-8'))
        ctr += 1
    print("Simulation finished (reached t="+str(Params.tf)+")")
    return int(outp.split()[-1])


def getSimData(iteration):
    global repo_path
    PropNames = ['00-Polarization', '01-Velocity', '02-Vorticity_0_0', '02-Vorticity_0_1', '02-Vorticity_1_0', '02-Vorticity_1_1', '03-ExternalForce', '04-Pressure', '05-StrainRate_0_0', '05-StrainRate_0_1', '05-StrainRate_1_0', '05-StrainRate_1_1', '06-Stress_0_0',
                 '06-Stress_0_1', '06-Stress_1_0', '06-Stress_1_1', '07-MolecularField', '08-DPOL', '09-DV', '10-VRHS', '11-f1', '12-f2', '13-f3', '14-f4', '15-f5', '16-f6', '17-V_T', '18-DIV', '19-DELMU', '20-HPB', '21-FrankEnergy', '22-R', 'SubsetNumber', 'domain']
    if not os.path.exists(repo_path+'/bin/Active2dOutput/Polar_'+str(iteration)+'.pvtp'):
        print("Iteration does not exists")
        return None
    reader = vtk.vtkXMLPPolyDataReader()
    reader.SetFileName(
        repo_path+'/bin/Active2dOutput/Polar_'+str(iteration)+'.pvtp')
    reader.Update()
    nodes = vtk_to_numpy(dsa.WrapDataObject(reader.GetOutput()).Points)
    t = float(linecache.getline(repo_path+'/bin/Active2dOutput/Polar_' +
              str(iteration)+'.pvtp', 5, module_globals=None))
    Pol = vtk_to_numpy(
        reader.GetOutput().GetPointData().GetArray(PropNames[0]))
    Vel = vtk_to_numpy(
        reader.GetOutput().GetPointData().GetArray(PropNames[1]))
    FE = vtk_to_numpy(
        reader.GetOutput().GetPointData().GetArray('21-FrankEnergy'))
    x, y = nodes[:, 0], nodes[:, 1]
    return x, y, Pol, Vel, FE, t


def VizualizeIteration(iteration=0):
    x, y, Pol, Vel, FE, t = getSimData(iteration)
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


def VizualizeAnimate(finalStep, jump=1):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    x, y, Pol, Vel, FE, t = getSimData(finalStep)
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
        x, y, Pol, Vel, FE, t = getSimData(i)
        q1.set_UVC(Pol[:, 0], Pol[:, 1], FE)
        q2.set_UVC(Vel[:, 0], Vel[:, 1], np.linalg.norm(Vel, axis=1))
        tittle.set_text("Time:"+str(round(t, 2)))
    anim = FuncAnimation(fig, animate, np.arange(
        1, finalStep, jump), interval=150)
    display(HTML(anim.to_html5_video()))
    plt.close()


class SpringLatticeParameters:

    thickness = 1.0
    Geometry = "flat-square"
    length = 10.0
    width = 10.0
    nematic_direction = [1, 0, 0]  # some vector function of x,y
    lambda_n = 1.5  # some scalar function of x,y
    lambda_n_perp = 0.5  # some scalar function of x,y
    lambda_h = 1  # some scalar function of x,y

    def load_mesh(self, thick_mesh_file=None):  # find a better way to refer to the mesh file

        if self.mesh_geometry == "square":
            if thick_mesh_file is None:
                thick_mesh_file = os.path.join(repo_path, "src/meshes/square.pkl")
            # replace the line below to either create a mesh or read a mesh
            [balls_df, springs_df] = pickle.load(open(thick_mesh_file, 'rb'))
            # setting thickness
            balls_df["z"][len(balls_df)//2:] = self.thickness
            # offsetting the mesh so that it starts from y = 0
            balls_df["y"] = balls_df["y"] - balls_df["y"].min()
            # update springs
            springs_df = update_springs(springs_df, balls_df[['x', 'y', 'z']])

        elif self.mesh_geometry == "circle":
            if thick_mesh_file is None:
                thick_mesh_file = os.path.join(repo_path, "meshes/circle.pkl")
            [balls_df, springs_df] = pickle.load(open(thick_mesh_file, 'rb'))
            # setting thickness
            balls_df["z"][len(balls_df)//2:] = self.thickness
            # if any of the balls have x = 0 then offset by small amount
            ind = balls_df[balls_df["x"] == 0].index
            balls_df.loc[ind, "x"] = 0.00001
            # update springs
            springs_df = update_springs(springs_df, balls_df[['x', 'y', 'z']])
            # compute the polar coordinates of the balls
            balls_df["r"] = np.sqrt(balls_df["x"]**2 + balls_df["y"]**2)
            balls_df["theta"] = np.arctan2(balls_df["y"], balls_df["x"])

        self.balls = balls_df
        self.springs = springs_df
        self.init_positions = balls_df[["x", "y", "z"]]

    def load_strain_pattern(self, theta_o=0.2, lambda1=1.25, lambda2=1.25**(-0.5), lambda3=1, step_size=0.1):

        def get_lambda_tensor(pos_vector, theta_o=np.pi/4, step_size=0.01, lambda1=1.25, lambda2=1.25**(-0.5), lambda3=1,):
            """
            This function computes the director field for a given point in space
            and the magnitude of the lambda coefficient at that point
            """
            # compute theta
            if self.nematic_coordinates == "zig-zag":
                theta = theta_o if (
                    pos_vector[1]//step_size) % 2 == 0 else np.pi-theta_o
                director1 = np.array([np.cos(theta), np.sin(theta), 0])
                director2 = np.array([-np.sin(theta), np.cos(theta), 0])
                director3 = np.array([0, 0, 1])
                lambda_tensor = lambda1*np.tensordot(director1, director1, axes=0) + lambda2*np.tensordot(
                    director2, director2, axes=0) + lambda3*np.tensordot(director3, director3, axes=0)

            elif self.nematic_coordinates == "polar":
                theta = np.arctan2(pos_vector[1], pos_vector[0])
                r = np.sqrt(pos_vector[0]**2 + pos_vector[1]**2)
                director1 = np.array([pos_vector[0]/r, pos_vector[1]/r, 0])
                director2 = np.array([-pos_vector[1]/r, pos_vector[0]/r, 0])
                director3 = np.array([0, 0, 1])
                lambda_tensor = lambda1(r)*np.tensordot(director1, director1, axes=0) \
                    + lambda2(r)*np.tensordot(director2, director2, axes=0) \
                    + lambda3(r)*np.tensordot(director3, director3, axes=0)

            elif self.nematic_coordinates == "cartesian":
                director1 = [1, 0, 0]
                director2 = [0, 1, 0]
                director3 = [0, 0, 1]
                [x, y, z] = pos_vector
                lambda_tensor = lambda1(x, y)*np.tensordot(director1, director1, axes=0) \
                    + lambda2(x, y)*np.tensordot(director2, director2, axes=0) \
                    + lambda3(x, y)*np.tensordot(director3, director3, axes=0)

            return (lambda_tensor)

        # function to compute the rest length of a single spring
        def get_rest_length(row):
            # this function will take a row of the springs_df and return the rest length

            # get lambda tensor at one endpoint
            lambda_alpha = get_lambda_tensor(
                balls_df.iloc[row["ball1"]][["x", "y", "z"]].values,
                theta_o=theta_o, step_size=step_size, lambda1=lambda1, lambda2=lambda2, lambda3=lambda3,
            )

            # get lambda tensor at other endpoint
            lambda_beta = get_lambda_tensor(
                balls_df.iloc[row["ball2"]][["x", "y", "z"]].values,
                theta_o=theta_o, step_size=step_size, lambda1=lambda1, lambda2=lambda2, lambda3=lambda3,
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

        balls_df = self.balls
        springs_df = self.springs
        # compute the rest length of each spring
        springs_df['l0'] = springs_df.apply(get_rest_length, axis=1)
        springs_df['l0_target'] = springs_df['l0']
        springs_df['l1_initial'] = springs_df['l1']
        self.springs = springs_df

    def add_noise(self, noise=0.1):
        # add noise to the initial positions
        self.balls[["x", "y", "z"]] = self.balls[["x", "y", "z"]] + \
            np.abs(noise*np.random.randn(*self.balls[["x", "y", "z"]].shape))
        self.springs = update_springs(
            self.springs, self.balls[['x', 'y', 'z']])

    def visualize(self, ax=None, fig=None, mode="discrete", x='x', y='y',
                  title="Strain Pattern", color_min=0.7, color_max=1.3, state="final",
                  tick_fontsize=10, label_fontsize=16, cbar_name=r'$\frac{l_{rest}}{l_{init}}$', title_fontsize=20,
                  xlim_min=-1.2, xlim_max=1.2, ylim_min=-1.2, ylim_max=1.2,
                  ):

        if mode == "continuous":
            pass

        elif mode == "discrete":

            balls_df = self.balls
            springs_df = self.springs

            if state == "initial":
                balls_df[['x', 'y', 'z']] = self.init_positions.values
            elif state == "final":
                balls_df[['x', 'y', 'z']] = self.final_positions.values
            springs_df = update_springs(springs_df, balls_df[['x', 'y', 'z']])

            if ax is None:
                plot_shell(balls_df, springs_df, x=x, y=y,  # filename=dirname + 'sim_output/top_view_init.pdf',
                           cbar_name=cbar_name, title=title, color_min=color_min, color_max=color_min, cmap="RdBu_r", norm="linear",)
            else:
                fig, ax = plot_shell_on_given_ax(balls_df, springs_df, x=x, y=y, ax=ax, fig=fig,
                                                 xlim_min=xlim_min, xlim_max=xlim_max, ylim_min=ylim_min, ylim_max=ylim_max,
                                                 title=title, color_min=color_min, color_max=color_max, cmap="RdBu_r", return_ax=True,
                                                 show_cbar=True, tick_fontsize=tick_fontsize, cbar_name=cbar_name, label_fontsize=label_fontsize,
                                                 )
                # title
                ax.set_title(title, fontsize=title_fontsize)
                return fig, ax

        # return ax

    def RunSpringLatticeSimulation(self, dt=0.01, tol=1e-6, csv_t_save=500):
        global repo_path
        dirname = repo_path+'/bin/SpringLatticeOutput/'
        bin_dir = repo_path+'/bin/'

        os.makedirs(dirname, exist_ok=True)
        os.makedirs(dirname+'runfiles/', exist_ok=True)
        os.makedirs(dirname+'sim_output/', exist_ok=True)

        [balls_df, springs_df] = initialize_cpp_simulation(
            self.balls, self.springs, dt=dt, csv_t_save=csv_t_save, tol=tol, path=dirname)

        filelist = ['Makefile', 'main.cpp']
        for file in filelist:
            shutil.copy(bin_dir + file, dirname)

        # running the simulation
        print('$$$$$$$ Running openfpm $$$$$$$')
        os.system("cd " + dirname +
                  " && /usr/local/openfpm/source/openfpm_vars && make && SpringLattice")
        print('$$$$ Exit OpenFPM $$$$')

        # access the output files
        final_balls_df = pd.read_csv(dirname + 'files/final_0_0.csv')
        self.final_positions = pd.DataFrame(
            final_balls_df[['x[0]', 'x[1]', 'x[2]']].values, columns=['x', 'y', 'z'])
        self.balls[['x', 'y', 'z']] = self.final_positions
        self.springs = update_springs(
            self.springs, self.balls[['x', 'y', 'z']])

        return (self)
