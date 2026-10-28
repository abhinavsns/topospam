import numpy as np
import os,shutil
import subprocess,vtk,linecache
import matplotlib.pyplot as plt
from vtk.util.numpy_support import vtk_to_numpy
from vtk.numpy_interface import dataset_adapter as dsa
from IPython.display import HTML, display
from ipywidgets import widgets, interact_manual, FloatProgress
from matplotlib.animation import FuncAnimation
from topospam_utils import animation_to_html, find_mpirun, require_binary, run_backend


def _mpirun_command(nCores):
    cmd = [find_mpirun(), "-np", str(nCores)]
    if os.geteuid() == 0:  # e.g. a container started as root
        cmd.append("--allow-run-as-root")
    return cmd + ["../Active2d", "../active2d.csv"]


class active_fluid_2d:
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

    def __init__(self,path):
        self.params = self.Params()
        self.repo_path=path

    def RunSimulation(self,nCores=1):
        SimParams=self.params.get()
        np.savetxt(self.repo_path+'/cpp/active2d.csv', SimParams, delimiter=' ')
        output_dir = os.path.join(self.repo_path, 'cpp', 'Active2dOutput')

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

        require_binary(os.path.join(self.repo_path, 'cpp', 'Active2d'))
        output_lines = run_backend(_mpirun_command(nCores), cwd=output_dir)
        # Print the timing and step-count lines written at the end of the run
        summary = [line for line in output_lines
                   if line.startswith(("The simulation took", "Total number of Time steps"))]
        if not summary or not summary[-1].startswith("Total number of Time steps"):
            raise RuntimeError("Active2d did not report the number of time steps. Last output:\n"
                               + "\n".join(output_lines[-30:]))
        for line in summary:
            print(line)
        return int(summary[-1].split()[-1])


    def RunSimulationWithProgress(self, nCores=1):
        SimParams=self.params.get()
        np.savetxt(self.repo_path+'/cpp/active2d.csv', SimParams, delimiter=' ')
        output_dir = os.path.join(self.repo_path, 'cpp', 'Active2dOutput')
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
        
        require_binary(os.path.join(self.repo_path, 'cpp', 'Active2d'))
        p = subprocess.Popen(_mpirun_command(nCores),
                            stdout=subprocess.PIPE, cwd=output_dir)
        f = FloatProgress(value=0.0, min=0.0, max=self.params.tf)  # instantiate the bar
        print("Simulation Progress (Based on current time value)")
        display(f)  # display the bar
        ctr = 0
        outp = "0"
        for line in iter(p.stdout.readline, b''):
            outp = line.rstrip().decode('utf-8')
            try:
                f.value = float(outp)
            except ValueError:
                print(outp)
            ctr += 1
        if p.wait() != 0:
            raise RuntimeError(f"Active2d failed with exit code {p.returncode}")
        print("Simulation finished (reached t="+str(self.params.tf)+")")
        return int(outp.split()[-1])


    def getSimData(self,iteration):
        PropNames = ['00-Polarization', '01-Velocity', '02-Vorticity_0_0', '02-Vorticity_0_1', '02-Vorticity_1_0', '02-Vorticity_1_1', '03-ExternalForce', '04-Pressure', '05-StrainRate_0_0', '05-StrainRate_0_1', '05-StrainRate_1_0', '05-StrainRate_1_1', '06-Stress_0_0',
                    '06-Stress_0_1', '06-Stress_1_0', '06-Stress_1_1', '07-MolecularField', '08-DPOL', '09-DV', '10-VRHS', '11-f1', '12-f2', '13-f3', '14-f4', '15-f5', '16-f6', '17-V_T', '18-DIV', '19-DELMU', '20-HPB', '21-FrankEnergy', '22-R', 'SubsetNumber', 'domain']
        if not os.path.exists(self.repo_path+'/cpp/Active2dOutput/Polar_'+str(iteration)+'.pvtp'):
            print("Iteration does not exists")
            return None
        reader = vtk.vtkXMLPPolyDataReader()
        reader.SetFileName(
            self.repo_path+'/cpp/Active2dOutput/Polar_'+str(iteration)+'.pvtp')
        reader.Update()
        nodes = vtk_to_numpy(dsa.WrapDataObject(reader.GetOutput()).Points)
        linecache.checkcache()  # the file is rewritten when the simulation is run again
        t = float(linecache.getline(self.repo_path+'/cpp/Active2dOutput/Polar_' +
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
        display(HTML(animation_to_html(anim)))
        plt.close()