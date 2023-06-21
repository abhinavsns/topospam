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

class ActiveFluidParameters:
    Gd_Sz=41
    Ks=1.0
    Kb=1.0
    dmu=0.0
    nu=0.0
    zeta=1.0
    lamda=1.0
    eta=1.0
    gama=1.0
    tf=1
    dt=1e-2
    wrat=1
    Vtol=1e-2
    absttol=1e-3
    relttol=1e-3
def RunActiveFluidSimulation(Params,nCores=1):
    data=np.array([Params.Gd_Sz,
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
    np.savetxt('bin/active2d.csv', data, delimiter=' ')
    os.popen('rm -r bin/Active2dOutput/* && mkdir -p bin/Active2dOutput')
    stream=os.popen('cd ./bin/Active2dOutput && source ~/openfpm_vars && mpirun -np '+str(nCores)+' ../Active2D ../active2d.csv')
    output = stream.readlines()
    print(output[-2])
    print(output[-1])
    return int(output[-1].split()[-1])
def RunActiveFluidSimulationWithProgress(Params,nCores=1):
    data=np.array([Params.Gd_Sz,
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
    np.savetxt('bin/active2d.csv', data, delimiter=' ')
    os.popen('rm -rf bin/Active2dOutput/*')
    time.sleep(1)
    os.popen('mkdir -p bin/Active2dOutput')
    p=subprocess.Popen('source ~/openfpm_vars && ../Active2d ../active2d.csv', stdout=subprocess.PIPE,cwd=os.getcwd()+'/bin/Active2dOutput',shell=True) 
    f = FloatProgress(value=0.0,min=0.0, max=Params.tf) # instantiate the bar
    print("Simulation Progress (Based on current time value)")
    display(f) # display the bar
    ctr=0
    for line in iter(p.stdout.readline, b''):
        if(ctr>5):
            outp=line.rstrip().decode('utf-8')
            if(outp.split()[0]=="Total" or outp.split()[0]=="The"):
                print(outp)
            else:
                f.value = float(line.rstrip().decode('utf-8'))
        ctr+=1
    print("Simulation finished (reached t="+str(Params.tf)+")")            
    return int(outp.split()[-1])
    
def getSimData(iteration):
    PropNames = ['00-Polarization', '01-Velocity', '02-Vorticity_0_0', '02-Vorticity_0_1', '02-Vorticity_1_0', '02-Vorticity_1_1', '03-ExternalForce', '04-Pressure', '05-StrainRate_0_0', '05-StrainRate_0_1', '05-StrainRate_1_0', '05-StrainRate_1_1', '06-Stress_0_0', '06-Stress_0_1', '06-Stress_1_0', '06-Stress_1_1', '07-MolecularField', '08-DPOL', '09-DV', '10-VRHS', '11-f1', '12-f2', '13-f3', '14-f4', '15-f5', '16-f6', '17-V_T', '18-DIV', '19-DELMU', '20-HPB', '21-FrankEnergy', '22-R', 'SubsetNumber', 'domain']
    if not os.path.exists('bin/Active2dOutput/Polar_'+str(iteration)+'.pvtp'):
        print("Iteration does not exists")
        return None 
    reader = vtk.vtkXMLPPolyDataReader()
    reader.SetFileName('bin/Active2dOutput/Polar_'+str(iteration)+'.pvtp')
    reader.Update()
    nodes = vtk_to_numpy(dsa.WrapDataObject(reader.GetOutput()).Points)
    t=float(linecache.getline('bin/Active2dOutput/Polar_'+str(iteration)+'.pvtp', 5, module_globals=None))
    Pol = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(PropNames[0]))
    Vel = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(PropNames[1]))
    FE = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray('21-FrankEnergy'))
    x,y= nodes[:,0] , nodes[:,1]
    return x,y,Pol,Vel,FE,t

def VizualizeIteration(iteration=0):
    x,y,Pol,Vel,FE,t=getSimData(iteration)
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(16, 6))
    q1=ax1.quiver(x,y,Pol[:,0],Pol[:,1],FE,cmap=plt.cm.jet)
    fig.colorbar(q1, cmap=plt.cm.jet,label='Frank Energy Density',ax=ax1)
    q2=ax2.quiver(x,y,Vel[:,0],Vel[:,1],np.linalg.norm(Vel,axis=1),cmap=plt.cm.viridis)
    ax1.set_title("Polarity Vectors")
    ax2.set_title("Velocity Vectors")
    fig.colorbar(q2, cmap=plt.cm.viridis,label='Velocity Magnitude',ax=ax2)
    fig.suptitle('Time:'+str(round(t,2)))
    #ax1.annotate('Time:'+str(round(t,2)), xy=(-0.5, 11), xycoords='data', annotation_clip=False,size=15)
    plt.show()
    return plt
def VizualizeAnimate(finalStep,jump=1):
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(10, 4))
    x,y,Pol,Vel,FE,t=getSimData(finalStep)
    q1=ax1.quiver(x,y,Pol[:,0],Pol[:,1],FE,cmap=plt.cm.jet,width=0.005,headlength=4)
    q2=ax2.quiver(x,y,Vel[:,0],Vel[:,1],np.linalg.norm(Vel,axis=1),cmap=plt.cm.viridis,width=0.005,headlength=4)
    ax1.set_title("Polarity Vectors")
    ax2.set_title("Velocity Vectors")
    c1=fig.colorbar(q1, cmap=plt.cm.jet,label='Frank Energy Density',ax=ax1)
    c2=fig.colorbar(q2, cmap=plt.cm.viridis,label='Velocity Magnitude',ax=ax2)
    tittle=fig.suptitle("Time:"+str(round(t,2)))
    #anno=ax1.annotate('Time:'+str(t), xy=(-0.5, 10), xycoords='data', annotation_clip=False,size=15)
    def animate(i):
        x,y,Pol,Vel,FE,t=getSimData(i)
        q1.set_UVC(Pol[:,0],Pol[:,1],FE)
        q2.set_UVC(Vel[:,0],Vel[:,1],np.linalg.norm(Vel,axis=1))
        tittle.set_text("Time:"+str(round(t,2)))
    anim = FuncAnimation(fig, animate,np.arange(1,finalStep,jump),interval=150)
    display(HTML(anim.to_html5_video()))
    plt.close()



class SpringLatticeParameters:
    seed = 0
    nb_iterations = 20
    thickness = 0.1
    R_initial = 1
    theta_DV = 0.1931
    DV_present = True
    #'outDV_gradient':True,
    #'volume_conservation':False,
    k_type = 'k_c'
    theta_max =  0.8662, #32*np.pi/180
    theta_ref =  0.8662
    overwrite_old_simulation = True
    #'tol':1e-6,
    dt = 0.01
    angle_of_rotation = 0.1
    mesh_refine_factor = 30
    tol_by_dt = 1e-8
    isotropic_contribution = "all"
    anisotropic_contribution = "all"
    height_contribution = "no"
    genotype = "ecadGFPnbG4",#"fit_lambdas_df_ecadGFPnbG4.pkl", #fit_lambdas_df_ecadGFPnbG4myoVI.pk
    disc_name = 'isotropic_homogeneous'
    lambda_anisotropic_coeffs = : [0, 1.1]
def RunSpringLatticeSimulation(Params,nCores=1):
    data=np.array([Params.Gd_Sz,
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
    np.savetxt('bin/active2d.csv', data, delimiter=' ')
    os.popen('rm -r bin/Active2dOutput/* && mkdir -p bin/Active2dOutput')
    stream=os.popen('cd ./bin/Active2dOutput && source ~/openfpm_vars && mpirun -np '+str(nCores)+' ../Active2D ../active2d.csv')
    output = stream.readlines()
    print(output[-2])
    print(output[-1])
    return int(output[-1].split()[-1])
def RunSpringLatticeSimulationWithProgress(Params,nCores=1):
    data=np.array([Params.Gd_Sz,
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
    np.savetxt('bin/active2d.csv', data, delimiter=' ')
    os.popen('rm -rf bin/Active2dOutput/*')
    time.sleep(1)
    os.popen('mkdir -p bin/Active2dOutput')
    p=subprocess.Popen('source ~/openfpm_vars && ../Active2d ../active2d.csv', stdout=subprocess.PIPE,cwd=os.getcwd()+'/bin/Active2dOutput',shell=True) 
    f = FloatProgress(value=0.0,min=0.0, max=Params.tf) # instantiate the bar
    print("Simulation Progress (Based on current time value)")
    display(f) # display the bar
    ctr=0
    for line in iter(p.stdout.readline, b''):
        if(ctr>5):
            outp=line.rstrip().decode('utf-8')
            if(outp.split()[0]=="Total" or outp.split()[0]=="The"):
                print(outp)
            else:
                f.value = float(line.rstrip().decode('utf-8'))
        ctr+=1
    print("Simulation finished (reached t="+str(Params.tf)+")")            
    return int(outp.split()[-1])
    
def getSpringLatticeSimData(iteration):
    PropNames = ['00-Polarization', '01-Velocity', '02-Vorticity_0_0', '02-Vorticity_0_1', '02-Vorticity_1_0', '02-Vorticity_1_1', '03-ExternalForce', '04-Pressure', '05-StrainRate_0_0', '05-StrainRate_0_1', '05-StrainRate_1_0', '05-StrainRate_1_1', '06-Stress_0_0', '06-Stress_0_1', '06-Stress_1_0', '06-Stress_1_1', '07-MolecularField', '08-DPOL', '09-DV', '10-VRHS', '11-f1', '12-f2', '13-f3', '14-f4', '15-f5', '16-f6', '17-V_T', '18-DIV', '19-DELMU', '20-HPB', '21-FrankEnergy', '22-R', 'SubsetNumber', 'domain']
    if not os.path.exists('bin/Active2dOutput/Polar_'+str(iteration)+'.pvtp'):
        print("Iteration does not exists")
        return None 
    reader = vtk.vtkXMLPPolyDataReader()
    reader.SetFileName('bin/Active2dOutput/Polar_'+str(iteration)+'.pvtp')
    reader.Update()
    nodes = vtk_to_numpy(dsa.WrapDataObject(reader.GetOutput()).Points)
    t=float(linecache.getline('bin/Active2dOutput/Polar_'+str(iteration)+'.pvtp', 5, module_globals=None))
    Pol = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(PropNames[0]))
    Vel = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(PropNames[1]))
    FE = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray('21-FrankEnergy'))
    x,y= nodes[:,0] , nodes[:,1]
    return x,y,Pol,Vel,FE,t

def SpringLatticeVizualizeIteration(iteration=0):
    x,y,Pol,Vel,FE,t=getSimData(iteration)
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(16, 6))
    q1=ax1.quiver(x,y,Pol[:,0],Pol[:,1],FE,cmap=plt.cm.jet)
    fig.colorbar(q1, cmap=plt.cm.jet,label='Frank Energy Density',ax=ax1)
    q2=ax2.quiver(x,y,Vel[:,0],Vel[:,1],np.linalg.norm(Vel,axis=1),cmap=plt.cm.viridis)
    ax1.set_title("Polarity Vectors")
    ax2.set_title("Velocity Vectors")
    fig.colorbar(q2, cmap=plt.cm.viridis,label='Velocity Magnitude',ax=ax2)
    fig.suptitle('Time:'+str(round(t,2)))
    #ax1.annotate('Time:'+str(round(t,2)), xy=(-0.5, 11), xycoords='data', annotation_clip=False,size=15)
    plt.show()
    return plt
def SpringLatticeVizualizeAnimate(finalStep,jump=1):
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(10, 4))
    x,y,Pol,Vel,FE,t=getSimData(finalStep)
    q1=ax1.quiver(x,y,Pol[:,0],Pol[:,1],FE,cmap=plt.cm.jet,width=0.005,headlength=4)
    q2=ax2.quiver(x,y,Vel[:,0],Vel[:,1],np.linalg.norm(Vel,axis=1),cmap=plt.cm.viridis,width=0.005,headlength=4)
    ax1.set_title("Polarity Vectors")
    ax2.set_title("Velocity Vectors")
    c1=fig.colorbar(q1, cmap=plt.cm.jet,label='Frank Energy Density',ax=ax1)
    c2=fig.colorbar(q2, cmap=plt.cm.viridis,label='Velocity Magnitude',ax=ax2)
    tittle=fig.suptitle("Time:"+str(round(t,2)))
    #anno=ax1.annotate('Time:'+str(t), xy=(-0.5, 10), xycoords='data', annotation_clip=False,size=15)
    def animate(i):
        x,y,Pol,Vel,FE,t=getSimData(i)
        q1.set_UVC(Pol[:,0],Pol[:,1],FE)
        q2.set_UVC(Vel[:,0],Vel[:,1],np.linalg.norm(Vel,axis=1))
        tittle.set_text("Time:"+str(round(t,2)))
    anim = FuncAnimation(fig, animate,np.arange(1,finalStep,jump),interval=150)
    display(HTML(anim.to_html5_video()))
    plt.close()