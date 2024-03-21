def RunSpringLatticeSimulation(Params, dt=0.01, tol=1e-6, csv_t_save=500, nCores = 1):
    global repo_path
    dirname = repo_path+'/bin/SpringLatticeOutput/'
    bin_dir = repo_path+'/bin/'

    os.makedirs(dirname, exist_ok=True)
    os.makedirs(dirname+'runfiles/', exist_ok=True)
    os.makedirs(dirname+'sim_output/', exist_ok=True)

    [balls_df, springs_df] = initialize_cpp_simulation(
        Params.balls, Params.springs, dt=dt, csv_t_save=csv_t_save, tol=tol, path=dirname)

    filelist = ['Makefile', 'SpringLattice.cpp']
    for file in filelist:
        shutil.copy(bin_dir + file, dirname)

    # running the simulation
    print('$$$$$$$ Running openfpm $$$$$$$')

    cmd = [
        "/bin/bash", "-c",
        "source /usr/local/openfpm/source/openfpm_vars && mpirun -np " +
        str(nCores)+" ../SpringLattice"
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

    print('$$$$ Exit OpenFPM $$$$')

    # access the output files
    final_balls_df = pd.read_csv(dirname + 'files/final_0_0.csv')
    Params.final_positions = pd.DataFrame(
        final_balls_df[['x[0]', 'x[1]', 'x[2]']].values, columns=['x', 'y', 'z'])
    Params.balls[['x', 'y', 'z']] = Params.final_positions
    Params.springs = update_springs(
        Params.springs, Params.balls[['x', 'y', 'z']])

    return (Params)
