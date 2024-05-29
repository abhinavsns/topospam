import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mesh_methods import *

class spring_lattice:
    class Params:
        def __init__(self):
            self.thickness = 1.0
            self.lambda_tensor_diagonal = [1.5, 1/1.5, 1]
            self.nematic_coordinates = "polar"
            self.mesh_geometry = "./mesh_utils_spring_lattice/circle.pkl"
            self.balls = None
            self.springs = None
            self.init_positions = None

    def __init__(self,path):
        self.params = self.Params()
        self.repo_path=path

    def load_mesh(self, thick_mesh_file=None):  # find a better way to refer to the mesh file
        if self.params.mesh_geometry == "square":
            if thick_mesh_file is None:
                thick_mesh_file = os.path.join(
                    self.repo_path, "examples/meshe_utils_spring_lattice/square.pkl")
            # replace the line below to either create a mesh or read a mesh
            [balls_df, springs_df] = pickle.load(open(thick_mesh_file, 'rb'))
            # setting thickness
            balls_df["z"][len(balls_df)//2:] = self.params.thickness
            # offsetting the mesh so that it starts from y = 0
            balls_df["y"] = balls_df["y"] - balls_df["y"].min()
            # update springs
            springs_df = update_springs(springs_df, balls_df[['x', 'y', 'z']])

        else:
            thick_mesh_file = self.params.mesh_geometry
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


    def initialize_strain_pattern(self):
        balls_df = self.params.balls
        springs_df = self.params.springs
        lambda_tensor_diagonal = self.params.lambda_tensor_diagonal
        nematic_coordinates = self.params.nematic_coordinates
        # compute the rest length of each spring
        springs_df['l0'] = springs_df.apply(self.get_rest_length, axis=1, lambda_tensor_diagonal=lambda_tensor_diagonal,
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


    def add_noise(self, noise=0.1, seed=0):
        # add noise to the initial positions
        np.random.seed(seed)
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
        if mode == "continuous":
            pass

        elif mode == "discrete":

            balls_df = self.params.balls
            springs_df = self.params.springs

            if state == "initial":
                balls_df[['x', 'y', 'z']] = self.params.init_positions.values
            elif state == "final":
                balls_df[['x', 'y', 'z']] = self.params.final_positions.values
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
            fig.show()
            return fig, ax

        # return ax


    def RunSimulation(self,dt=0.01, tol=1e-6, csv_t_save=500):
        dirname = self.repo_path+'/cpp/SpringLatticeOutput/'
        bin_dir = self.repo_path+'/cpp/'

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
                "make SpringLattice")
        os.system("cd " + dirname +
                " && ../SpringLattice")
        print('$$$$ Exit OpenFPM $$$$')

        # access the output files
        final_balls_df = pd.read_csv(dirname + 'files/final_0_0.csv')
        self.params.final_positions = pd.DataFrame(
            final_balls_df[['x[0]', 'x[1]', 'x[2]']].values, columns=['x', 'y', 'z'])
        self.params.balls[['x', 'y', 'z']] = self.params.final_positions
        self.params.springs = update_springs(
            self.params.springs, self.params.balls[['x', 'y', 'z']])
        return 