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

    def load_mesh(self, thick_mesh_file = "../src/TopoSPAM/meshes/square_boundary.pkl"): #find a better way to refer to the mesh file
        if self.mesh_geometry == "flat-square":

            #replace the line below to either create a mesh or read a mesh
            [balls_df, springs_df] = pickle.load(open(thick_mesh_file, 'rb'))
            # setting thickness
            balls_df["z"][len(balls_df)//2:] = self.thickness
            # offsetting the mesh so that it starts from y = 0
            balls_df["y"] = balls_df["y"] - balls_df["y"].min()
            #update springs
            #springs_df = update_springs(springs_df, balls_df[['x', 'y', 'z']])
            self.balls = balls_df
            self.springs = springs_df
        else:
            pass

    def load_strain_pattern(self, theta_o = 0.2, lambda1 = 1.25, lambda2 = 1.25**(-0.5), step_size = 0.1):

        def get_lambda_tensor(pos_vector, theta_o=np.pi/4, step_size=0.01, lambda1=1.25, lambda2=1, lambda3=1):
            """
            This function computes the director field for a given point in space
            and the magnitude of the lambda coefficient at that point
            """
            # compute theta
            theta = theta_o if (pos_vector[1]//step_size) % 2 == 0 else np.pi-theta_o
            # compute director1
            director1 = np.array([np.cos(theta), np.sin(theta), 0])
            # compute director2
            director2 = np.array([-np.sin(theta), np.cos(theta), 0])
            # compute director3
            director3 = np.array([0, 0, 1])
            #combine to form lambda tensor
            lambda_tensor = lambda1*np.tensordot(director1, director1, axes=0) + lambda2*np.tensordot(
                director2, director2, axes=0) + lambda3*np.tensordot(director3, director3, axes=0)

            return (lambda_tensor)

        #function to compute the rest length of a single spring
        def get_rest_length(row):
            # this function will take a row of the springs_df and return the rest length

            # get lambda tensor at one endpoint
            lambda_alpha = get_lambda_tensor(
                balls_df.iloc[row["ball1"]][["x", "y", "z"]].values,
                theta_o=theta_o, step_size=step_size, lambda1=lambda1, lambda2=lambda2,
            )
            
            # get lambda tensor at other endpoint
            lambda_beta = get_lambda_tensor(
                balls_df.iloc[row["ball2"]][["x", "y", "z"]].values,
                theta_o=theta_o, step_size=step_size, lambda1=lambda1, lambda2=lambda2,
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

    def visualize(Params, mode = "continuous", fig = None, ax = None):

        if ax == None:
            fig,ax = plt.subplots(1,1,figsize = (10,10))
        #draw the shape of the lattice
        if Params.Geometry == "flat-square":
            #draw a square of length Params.length and width Params.width
            ax.plot([0,Params.length,Params.length,0,0],[0,0,Params.width,Params.width,0],color = "black")  

        if mode == "continuous":
            pass
            
        elif mode == "discrete":
            pass

        return ax
    
    @classmethod
    def helloworld(cls):
        print("hello world")


    
