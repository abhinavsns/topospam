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

class SpringLatticeParameters:

    thickness = 1.0
    Geometry = "flat-square"
    length = 10.0
    width = 10.0
    nematic_direction = [1,0,0]#some vector function of x,y
    lambda_n = np.poly1d([1.5,1])#some scalar function of x,y
    lambda_n_perp = np.poly1d([0.5, 1])#some scalar function of x,y
    lambda_h = np.poly1d([1])#some scalar function of x,y

    def visualize(Params, mode = "continuous"):

        fig,ax = plt.subplots(1,1,figsize = (10,10))
        #draw the shape of the lattice
        if Params.Geometry == "flat-square":
            #draw a square of length Params.length and width Params.width
            ax.plot([0,Params.length,Params.length,0,0],[0,0,Params.width,Params.width,0],color = "black")  

        if mode == "continuous":
            pass
            
        elif mode == "discrete":
            pass
        