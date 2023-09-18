import numpy as np
import os, subprocess, time
import vtk
import matplotlib.pyplot as plt
from vtk.util.numpy_support import vtk_to_numpy
from vtk.numpy_interface import dataset_adapter as dsa
from IPython.display import HTML,display
from ipywidgets import widgets,interact_manual,FloatProgress
from matplotlib.animation import FuncAnimation
from matplotlib.offsetbox import AnchoredText
import linecache

from TopoSPAM.methods import *

class SpringLatticeParameters:

    thickness = 1.0
    Geometry = "flat-square"
    length = 10.0
    width = 10.0
    nematic_direction = [1,0,0]#some vector function of x,y
    lambda_n = 1.5 #some scalar function of x,y
    lambda_n_perp = 0.5 #some scalar function of x,y
    lambda_h = 1 #some scalar function of x,y

    def load_mesh(self, thick_mesh_file = None): #find a better way to refer to the mesh file

        if self.mesh_geometry == "square":
            if thick_mesh_file is None: thick_mesh_file = "../src/TopoSPAM/meshes/square.pkl"
            #replace the line below to either create a mesh or read a mesh
            [balls_df, springs_df] = pickle.load(open(thick_mesh_file, 'rb'))
            # setting thickness
            balls_df["z"][len(balls_df)//2:] = self.thickness
            # offsetting the mesh so that it starts from y = 0
            balls_df["y"] = balls_df["y"] - balls_df["y"].min()
            #update springs
            springs_df = update_springs(springs_df, balls_df[['x', 'y', 'z']])

        elif self.mesh_geometry == "circle":
            if thick_mesh_file is None: thick_mesh_file = "../src/TopoSPAM/meshes/circle.pkl"
            [balls_df, springs_df] = pickle.load(open(thick_mesh_file, 'rb'))
            # setting thickness
            balls_df["z"][len(balls_df)//2:] = self.thickness
            # if any of the balls have x = 0 then offset by small amount
            ind = balls_df[balls_df["x"] == 0].index
            balls_df.loc[ind, "x"] = 0.00001
            #update springs
            springs_df = update_springs(springs_df, balls_df[['x', 'y', 'z']])
            #compute the polar coordinates of the balls
            balls_df["r"] = np.sqrt(balls_df["x"]**2 + balls_df["y"]**2)
            balls_df["theta"] = np.arctan2(balls_df["y"], balls_df["x"])

        self.balls = balls_df
        self.springs = springs_df
        self.init_positions = balls_df[["x", "y", "z"]]

    def load_strain_pattern(self, theta_o = 0.2, lambda1 = 1.25, lambda2 = 1.25**(-0.5), lambda3=1, step_size = 0.1):

        def get_lambda_tensor(pos_vector, theta_o=np.pi/4, step_size=0.01, lambda1=1.25, lambda2=1.25**(-0.5), lambda3=1,):
            """
            This function computes the director field for a given point in space
            and the magnitude of the lambda coefficient at that point
            """
            # compute theta
            if self.nematic_coordinates == "zig-zag":
                theta = theta_o if (pos_vector[1]//step_size) % 2 == 0 else np.pi-theta_o
                director1 = np.array([np.cos(theta), np.sin(theta), 0])
                director2 = np.array([-np.sin(theta), np.cos(theta), 0])
                director3 = np.array([0, 0, 1])
                lambda_tensor = lambda1*np.tensordot(director1, director1, axes=0) + lambda2*np.tensordot(
                    director2, director2, axes=0) + lambda3*np.tensordot(director3, director3, axes=0)

            elif self.nematic_coordinates == "polar":
                theta = np.arctan2(pos_vector[1], pos_vector[0]) 
                r = np.sqrt(pos_vector[0]**2 + pos_vector[1]**2)
                director1 = np.array([pos_vector[0]/r, pos_vector[1]/r,0])
                director2 = np.array([-pos_vector[1]/r, pos_vector[0]/r,0])
                director3 = np.array([0, 0, 1])
                lambda_tensor = lambda1(r)*np.tensordot(director1, director1, axes=0) \
                + lambda2(r)*np.tensordot(director2, director2, axes=0) \
                + lambda3(r)*np.tensordot(director3, director3, axes=0)
            
            elif self.nematic_coordinates == "cartesian":
                director1 = [1,0,0]
                director2 = [0,1,0]
                director3 = [0,0,1]
                [x,y,z] = pos_vector
                lambda_tensor = lambda1(x,y)*np.tensordot(director1, director1, axes=0) \
                + lambda2(x,y)*np.tensordot(director2, director2, axes=0) \
                + lambda3(x,y)*np.tensordot(director3, director3, axes=0)

            return (lambda_tensor)

        #function to compute the rest length of a single spring
        def get_rest_length(row):
            # this function will take a row of the springs_df and return the rest length

            # get lambda tensor at one endpoint
            lambda_alpha = get_lambda_tensor(
                balls_df.iloc[row["ball1"]][["x", "y", "z"]].values,
                theta_o=theta_o, step_size=step_size, lambda1=lambda1, lambda2=lambda2, lambda3 = lambda3, 
            )
            
            # get lambda tensor at other endpoint
            lambda_beta = get_lambda_tensor(
                balls_df.iloc[row["ball2"]][["x", "y", "z"]].values,
                theta_o=theta_o, step_size=step_size, lambda1=lambda1, lambda2=lambda2, lambda3 = lambda3, 
            )
            # average the lambda tensors
            lambda_avg = (lambda_alpha + lambda_beta)/2
            # Multiply the average lambda tensor with the original spring vector
            # get the rest length by taking the norm of the above vector
            spring_vector = np.array(
                [row['x1'] - row['x2'], row['y1'] - row['y2'], row['z1'] - row['z2']])
            new_rest_length = np.linalg.norm(np.matmul(lambda_avg, spring_vector))

            return (new_rest_length)
        
        balls_df = self.balls
        springs_df = self.springs
        # compute the rest length of each spring
        springs_df['l0'] = springs_df.apply(get_rest_length, axis=1)
        springs_df['l0_target'] = springs_df['l0']
        springs_df['l1_initial'] = springs_df['l1']
        self.springs = springs_df

    def add_noise(self, noise = 0.1):
        # add noise to the initial positions
        self.balls[["x", "y", "z"]] = self.balls[["x", "y", "z"]] + np.abs(noise*np.random.randn(*self.balls[["x", "y", "z"]].shape))
        self.springs = update_springs(self.springs, self.balls[['x', 'y', 'z']])

    def visualize(self, ax = None, fig = None, mode = "discrete", x = 'x', y = 'y', 
                  title = "Strain Pattern", color_min = 0.7, color_max = 1.3, state = "final", 
                  tick_fontsize = 10, label_fontsize = 16, cbar_name=r'$\frac{l_{rest}}{l_{init}}$', title_fontsize = 20,
                  xlim_min = -1.2, xlim_max = 1.2, ylim_min = -1.2, ylim_max = 1.2,
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
                plot_shell(balls_df, springs_df, x=x, y=y, #filename=dirname + 'sim_output/top_view_init.pdf',
                        cbar_name=cbar_name, title=title, color_min=color_min, color_max=color_min, cmap="RdBu_r", norm="linear",)
            else:
                fig,ax = plot_shell_on_given_ax( balls_df, springs_df, x=x, y=y, ax = ax, fig = fig, 
                                                xlim_min = xlim_min, xlim_max = xlim_max, ylim_min = ylim_min, ylim_max = ylim_max,
                                       title = title, color_min=color_min, color_max=color_max, cmap="RdBu_r", return_ax=True, 
                                       show_cbar=True, tick_fontsize=tick_fontsize, cbar_name=cbar_name, label_fontsize=label_fontsize,
                )
                #title
                ax.set_title(title, fontsize=title_fontsize)
                return fig,ax

        #return ax

    def RunSpringLatticeSimulation(self, dt = 0.01, tol = 1e-6, csv_t_save = 500,
                                   dirname = "../only_local/test_run/", bin_dir = "../bin/") :

        os.makedirs(dirname, exist_ok=True)
        os.makedirs(dirname+'runfiles/', exist_ok=True)
        os.makedirs(dirname+'sim_output/', exist_ok=True)

        [balls_df, springs_df] = initialize_cpp_simulation(
            self.balls, self.springs, dt=dt, csv_t_save=csv_t_save, tol=tol, path=dirname)
        
        filelist = ['Makefile', 'main.cpp']
        for file in filelist: shutil.copy(bin_dir + file, dirname)

        # running the simulation
        print('$$$$$$$ Running openfpm $$$$$$$')
        os.system("cd " + dirname + " && source ~/openfpm_vars && make && grid")
        print('$$$$ Exit OpenFPM $$$$')

        #access the output files
        final_balls_df = pd.read_csv(dirname + 'files/final_0_0.csv')
        self.final_positions = pd.DataFrame(final_balls_df[['x[0]', 'x[1]', 'x[2]']].values, columns=['x', 'y', 'z'])
        self.balls[['x', 'y', 'z']] = self.final_positions
        self.springs = update_springs(self.springs, self.balls[['x', 'y', 'z']])

        return(self) 



    
    @classmethod
    def helloworld(cls):
        print("hello world")


    
