#define SE_CLASS1 // always use it
// #define PRINT_STACKTRACE
// #define STOP_ON_ERROR

#include "Vector/vector_dist.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "hash_map/hopscotch_map.h"
#include "hash_map/hopscotch_set.h"
#include <sys/stat.h>
#include <string>

// main program

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

constexpr int mass = 0;
constexpr int velocity = 1;
constexpr int force = 2;
constexpr int pid = 3;
constexpr int nbs = 4;
constexpr int onbs = 5;
constexpr int length = 6;
constexpr int n_length = 7;
constexpr int k_s = 8;
constexpr int n_s = 9;
constexpr int active = 10;
constexpr int spring_active = 11;
constexpr int n_nbs = 12;
constexpr int ext_force = 13;
constexpr int n_tri = 14;
constexpr int tri_V1 = 15;
constexpr int tri_V2 = 16;

constexpr int nbs_max = 22;

int n_particles = 0;
double n_springs = 0;

struct ModelCustom
{
    template <typename Decomposition, typename vector>
    inline void addComputation(Decomposition &dec,
                               vector &particles,
                               size_t cell_id,
                               size_t p)
    {
        dec.addComputationCost(cell_id, 300);
    }

    template <typename Decomposition>
    inline void applyModel(Decomposition &dec, size_t v)
    {
    }

    double distributionTol()
    {
        return 1.01;
    }
};

struct node_sim
{
    typedef boost::fusion::vector<float[3]> type;

    type data;

    struct attributes
    {
        static const std::string name[];
    };

    typedef float s_type;

    node_sim(){};

    static const unsigned int max_prop = 1;
};

const std::string node_sim::attributes::name[] = {"x"};

template <typename particle_type>
void write_connection_vtk(std::string file, particle_type &part, int n)
{
    Graph_CSR<node_sim, aggregate<int>> graph;

    auto it = part.getDomainAndGhostIterator();

    while (it.isNext())
    {
        auto p = it.get();

        graph.addVertex();
        graph.vertex(p.getKey()).template get<0>()[x] = part.getPos(p)[x];
        graph.vertex(p.getKey()).template get<0>()[y] = part.getPos(p)[y];
        graph.vertex(p.getKey()).template get<0>()[z] = part.getPos(p)[z];

        ++it;
    }

    auto it2 = part.getDomainIterator();

    while (it2.isNext())
    {
        auto p = it2.get();

        for (int j = 0; j < part.template getProp<n_nbs>(p); j++)
        {
            auto nnp = part.template getProp<onbs>(p)[j];

            graph.addEdge(p.getKey(), nnp);
        }

        ++it2;
    }

    VTKWriter<Graph_CSR<node_sim, aggregate<int>>, VTK_GRAPH> gv2(graph);
    gv2.write(file + "_" + std::to_string(create_vcluster().rank()) + "_" + std::to_string(n) + ".vtk");
}

template <typename particles_type>
void reconnect(particles_type &part)
{
    part.template ghost_get<pid>();

    auto it2 = part.getDomainAndGhostIterator();

    tsl::hopscotch_map<int, int> map;

    while (it2.isNext())
    {
        auto p = it2.get();

        map[part.template getProp<pid>(p)] = p.getKey();

        ++it2;
    }

    auto it = part.getDomainIterator();

    while (it.isNext())
    {
        auto p = it.get();

        for (int j = 0; j < part.template getProp<n_nbs>(p); j++)
        {
            int pid = part.template getProp<nbs>(p)[j];

            auto fnd = map.find(pid);
            if (fnd == map.end())
            {
                std::cout << "RECONNECTION FAILED " << pid << std::endl;
                part.write("crash");
                exit(1);
            }
            else
            {
                part.template getProp<onbs>(p)[j] = fnd->second;
            }
        }

        ++it;
    }
}

int main(int argc, char *argv[])
{

    struct stat st;
    std::string dir = "files/";

    if (stat("files", &st) != 0 || !S_ISDIR(st.st_mode))
    {
        mkdir("files", 0777);
    }

    // read simulation parameters from .csv
    std::ifstream param_file("runfiles/sim_params.csv");
    std::string param_line;

    std::getline(param_file, param_line);
    std::stringstream ss(param_line);

    double dt, tf, dim;
    double save_csv_d, save_vtk_d;
    double tol;
    // int nbs_m;

    ss >> dt;
    // ss >> tf;
    ss >> save_csv_d;
    // ss >> save_vtk_d;
    ss >> dim;
    ss >> tol;
    // ss >> folder;

    int save_csv = save_csv_d / dt;
    int save_vtk = save_vtk_d / dt;

    // string dir = folder;
    // ss >> nbs_m;

    // const int nbs_max  = nbs_m;

    std::cout << "dt " << dt << "\n";
    std::cout << "tf " << tf << "\n";
    std::cout << "save_csv " << save_csv << "\n";
    std::cout << "save_vtk " << save_vtk << "\n";
    std::cout << "dim " << dim << "\n";
    std::cout << "tol " << tol << "\n";

    size_t psz = 20;
    const size_t sz[3] = {psz, psz, 2};
    Box<3, double> box({-dim, -dim, -dim}, {dim, dim, dim});
    size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
    double spacing = 1.0;
    Ghost<3, double> ghost(spacing * 3);

    grid_sm<3, void> g(sz);

    size_t cell_sz[] = {32, 32, 1};
    grid_sm<3, void> g_cell(cell_sz);

    openfpm_init(&argc, &argv);
    //                                 mass   velocity force    pid    nbs         onbs          length          n_length        k_s             n_s          active  spr_act n_nbs  ext_frc  n_tri    tri_V1       tri_V2
    vector_dist<3, double, aggregate<double, double[3], double[3], int, int[nbs_max], int[nbs_max], double[nbs_max], double[nbs_max], double[nbs_max], double[nbs_max], int, int[nbs_max], int, double[3], int, int[nbs_max], int[nbs_max]>> particles(0, box, bc, ghost, g_cell);

    auto &v_cl = create_vcluster<>();
    particles.setPropNames({"mass", "velocity", "force", "pid", "nbs", "onbs", "length", "n_length", "k_s", "n_s", "active", "spring_active", "n_nbs", "ext_force", "n_tri", "tri_V1", "tri_V2"});

    if (v_cl.rank() == 0)
    {
        // Create an input filestream
        std::ifstream ball_file("runfiles/balls.csv");
        std::ifstream neighbours_file("runfiles/neighbours.csv");
        std::ifstream l0_file("runfiles/neigh_l0.csv");
        std::ifstream k_file("runfiles/neigh_k.csv");
        std::ifstream ns_file("runfiles/neigh_viscoelastic_coeff.csv");
        std::ifstream triV1_file("runfiles/triangles_Vertex1.csv");
        std::ifstream triV2_file("runfiles/triangles_Vertex2.csv");

        std::string line;
        std::string line_neigh;
        std::string line_l0;
        std::string line_k;
        std::string line_ns;
        std::string line_triV1;
        std::string line_triV2;

        while (std::getline(ball_file, line))
        {
            // std::cout << line << "\n";
            std::getline(neighbours_file, line_neigh);
            std::getline(l0_file, line_l0);
            std::getline(k_file, line_k);
            std::getline(ns_file, line_ns);
            std::getline(triV1_file, line_triV1);
            std::getline(triV2_file, line_triV2);

            std::stringstream ss(line);
            std::stringstream ss_neigh(line_neigh);
            std::stringstream ss_l0(line_l0);
            std::stringstream ss_k(line_k);
            std::stringstream ss_ns(line_ns);
            std::stringstream ss_triV1(line_triV1);
            std::stringstream ss_triV2(line_triV2);

            particles.add();

            double g_id;
            ss >> g_id;

            particles.getLastProp<pid>() = g_id;

            for (int i = 0; i < 3; i++)
            {
                double coord;
                ss >> coord;
                // if(i>1){
                //     std::cout << "Z : " << coord << "\n";
                // }
                particles.getLastPos()[i] = coord;
            }

            double n_neigh;
            ss >> n_neigh;
            particles.getLastProp<n_nbs>() = n_neigh;
            if (nbs_max < n_neigh)
            {
                std::cout << "INCREASE NBS_MAX to " << n_neigh << "\n";
            }

            double act;
            ss >> act;
            particles.getLastProp<active>() = (int)act;

            for (int i = 0; i < 3; i++)
            {
                double ext_force_input;
                ss >> ext_force_input;
                // ext_force_input = ext_force_x; // for now x = y =z not input from csv
                particles.getLastProp<ext_force>()[i] = ext_force_input;
            }

            double n_triangles; // number of triangles
            ss >> n_triangles;
            particles.getLastProp<n_tri>() = n_triangles;

            // std::cout << particles.template getLastProp<pid>() << " " << act << " "  << particles.getLastProp<active>() << "\n";

            n_springs = n_springs + n_neigh;

            for (int i = 0; i < n_neigh; i++)
            {

                double neigh, l0, k, ns;
                int act;
                ss_neigh >> neigh;
                ss_l0 >> l0;
                ss_k >> k;
                ss_ns >> ns;

                // std::cout << particles.getLastProp<pid>() << " " << neigh << " " << k << " " << ns << "\n";

                particles.getLastProp<nbs>()[i] = neigh;
                particles.getLastProp<n_length>()[i] = l0;
                particles.getLastProp<k_s>()[i] = k;
                particles.getLastProp<n_s>()[i] = ns;
            }

            for (int i = 0; i < n_triangles; i++)
            {

                int triangle_Vertex1, triangle_Vertex2;
                ss_triV1 >> triangle_Vertex1;
                ss_triV2 >> triangle_Vertex2;
                particles.getLastProp<tri_V1>()[i] = triangle_Vertex1;
                particles.getLastProp<tri_V2>()[i] = triangle_Vertex2;
            }

            particles.getLastProp<velocity>()[x] = 0;
            particles.getLastProp<velocity>()[y] = 0;
            particles.getLastProp<velocity>()[z] = 0;

            ++n_particles;
        }
    }

    particles.map();
    particles.ghost_get<pid>();

    reconnect(particles);
    //

    // write_connection_vtk(dir + "Test",particles,0);

    ModelCustom md;
    // particles.addComputationCosts(md);
    // particles.getDecomposition().decompose();
    particles.map();

    // particles.write("DEBUG");
    particles.getDecomposition().write(dir + "Decomposition");
    reconnect(particles);

    timer timer;
    timer.start();

    auto itp = particles.getDomainIterator();

    while (itp.isNext())
    {

        auto p = itp.get();
        Point<3, double> xp = particles.getPos(p);

        double fx = 0 + particles.template getProp<ext_force>(p)[x];
        double fy = 0 + particles.template getProp<ext_force>(p)[y];
        double fz = 0 + particles.template getProp<ext_force>(p)[z];
        double x1 = 0;
        double y1 = 0;
        double z1 = 0;

        if (particles.template getProp<active>(p) != 1 && particles.template getProp<active>(p) != 0)
        {
            std::cout << particles.template getProp<pid>(p) << " " << particles.template getProp<active>(p) << "\n";
        }

        for (int j = 0; j < particles.getProp<n_nbs>(p); j++)
        {

            auto nnp = particles.getProp<onbs>(p)[j];
            Point<3, double> xq = particles.getPos(nnp);

            if (xp != xq)
            {

                x1 = xq[x] - xp[x];
                y1 = xq[y] - xp[y];
                z1 = xq[z] - xp[z];

                double xm = sqrt(x1 * x1 + y1 * y1 + z1 * z1);

                particles.template getProp<length>(p)[j] = xm;

                fx = fx + particles.template getProp<k_s>(p)[j] * (particles.template getProp<length>(p)[j] - particles.template getProp<n_length>(p)[j]) * x1 / xm; //-particles.template getProp<n_s>(p)[j]*particles.template getProp<velocity>(p)[x];
                fy = fy + particles.template getProp<k_s>(p)[j] * (particles.template getProp<length>(p)[j] - particles.template getProp<n_length>(p)[j]) * y1 / xm; //-particles.template getProp<n_s>(p)[j]*particles.template getProp<velocity>(p)[y];
                fz = fz + particles.template getProp<k_s>(p)[j] * (particles.template getProp<length>(p)[j] - particles.template getProp<n_length>(p)[j]) * z1 / xm; //-particles.template getProp<n_s>(p)[j]*particles.template getProp<velocity>(p)[z];
                // count = count + 1;

                // Adding viscoelastic forces
                // fx=fx+particles.template getProp<n_s>(p)[j]*(particles.template getProp<velocity>(nnp)[x]-particles.template getProp<velocity>(p)[x])* x1/xm;
                // fy=fy+particles.template getProp<n_s>(p)[j]*(particles.template getProp<velocity>(nnp)[y]-particles.template getProp<velocity>(p)[y])* y1/xm;
                // fz=fz+particles.template getProp<n_s>(p)[j]*(particles.template getProp<velocity>(nnp)[z]-particles.template getProp<velocity>(p)[z])* z1/xm;
            }
        }

        particles.template getProp<force>(p)[x] = fx;
        particles.template getProp<force>(p)[y] = fy;
        particles.template getProp<force>(p)[z] = fz;

        // std::cout << fx << " " << fy << " " << fz << "\n";

        ++itp;
    }

    // openfpm_finalize();
    // return 0;

    // particles.write_frame("Forces",0);

    double t = 0;
    // particles.write_frame("Spring",0,CSV_WRITER);

    int ctr = 1;
    // Writing first frame
    // particles.write_frame(dir + "Spring",ctr,CSV_WRITER);

    double avg_movement = 99999; // to sum the displacement of each particle in order to check whether the simulation is moving

    // while(t<=tf){
    while (avg_movement > tol)
    {

        // double f_av = 0;
        // double d_av = 0;

        avg_movement = 0;

        auto it3 = particles.getDomainIterator();
        while (it3.isNext())
        {
            auto p = it3.get();

            // f_av += sqrt(particles.template getProp<force>(p)[x]*particles.template getProp<force>(p)[x] + particles.template getProp<force>(p)[y]*particles.template getProp<force>(p)[y] + particles.template getProp<force>(p)[z]*particles.template getProp<force>(p)[z]);

            // particles.template getProp<velocity>(p)[x]=particles.template getProp<velocity>(p)[x]+dt*particles.template getProp<force>(p)[x];
            // particles.template getProp<velocity>(p)[y]=particles.template getProp<velocity>(p)[y]+dt*particles.template getProp<force>(p)[y];
            // particles.template getProp<velocity>(p)[z]=particles.template getProp<velocity>(p)[z]+dt*particles.template getProp<force>(p)[z];

            // add hard surface
            // Changed here
            // if(particles.getPos(p)[z] == 0.0){
            //     if(particles.template getProp<velocity>(p)[z] < 0){
            //         particles.template getProp<velocity>(p)[z] = 0.0;
            //     }
            // }

            // double dx = dt*particles.template getProp<velocity>(p)[x];
            // double dy = dt*particles.template getProp<velocity>(p)[y];
            // double dz = dt*particles.template getProp<velocity>(p)[z];

            // d_av += sqrt(dx*dx + dy*dy + dz*dz);

            // particles.getPos(p)[x]=particles.getPos(p)[x]+dt*particles.template getProp<velocity>(p)[x];
            // particles.getPos(p)[y]=particles.getPos(p)[y]+dt*particles.template getProp<velocity>(p)[y];
            // particles.getPos(p)[z]=particles.getPos(p)[z]+dt*particles.template getProp<velocity>(p)[z];

            // Overdamped dynamics

            double dx = dt * particles.template getProp<force>(p)[x];
            double dy = dt * particles.template getProp<force>(p)[y];
            double dz = dt * particles.template getProp<force>(p)[z];

            particles.getPos(p)[x] = particles.getPos(p)[x] + dx;
            particles.getPos(p)[y] = particles.getPos(p)[y] + dy;
            particles.getPos(p)[z] = particles.getPos(p)[z] + dz;

            avg_movement += sqrt(dx * dx + dy * dy + dz * dz);

            // std::cout << t << " " << particles.template getProp<pid>(p) << " x_coord: ";
            // std::cout << particles.getPos(p)[x] << "\n";

            // add hard surface
            // Changed here
            // if(particles.getPos(p)[z] < 0){
            //     particles.getPos(p)[z] = 0.0;
            // }

            ++it3;
        }

        avg_movement = avg_movement / n_particles;
        // particles.write_frame("Debug",ctr);

        // only gives position, if we use any other properties
        // need to do ghost_get<prop1,prop2>
        if (ctr % 1000 == 0)
        {
            ModelCustom md;
            particles.addComputationCosts(md);
            particles.getDecomposition().decompose();
            particles.map();
            particles.ghost_get<>();

            reconnect(particles);
        }

        particles.ghost_get<>(SKIP_LABELLING);

        // double strain = 0;

        auto it2 = particles.getDomainIterator();
        while (it2.isNext())
        {
            auto p = it2.get();
            Point<3, double> xp = particles.getPos(p);
            // size_t count  = 0;
            double fx = 0 + particles.template getProp<ext_force>(p)[x];
            double fy = 0 + particles.template getProp<ext_force>(p)[y];
            double fz = 0 + particles.template getProp<ext_force>(p)[z];
            double x1 = 0;
            double y1 = 0;
            double z1 = 0;
            double xm = 0;
            double scalar_viscoelastic = 0;

            for (int j = 0; j < particles.getProp<n_nbs>(p); j++)
            {

                auto nnp = particles.getProp<onbs>(p)[j];
                Point<3, double> xq = particles.getPos(nnp);
                if (xp != xq)
                {

                    x1 = xq[x] - xp[x];
                    y1 = xq[y] - xp[y];
                    z1 = xq[z] - xp[z];

                    xm = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
                    particles.template getProp<length>(p)[j] = xm;

                    // std::cout << particles.template getProp<pid>(p) << " " << particles.template getProp<pid>(nnp) << ": ";
                    // std::cout << particles.template getProp<length>(p)[j] << " " << particles.template getProp<n_length>(p)[j] << "\n";

                    // Adding elastic forces
                    fx = fx + particles.template getProp<k_s>(p)[j] * (particles.template getProp<length>(p)[j] - particles.template getProp<n_length>(p)[j]) * x1 / xm; //-particles.template getProp<n_s>(p)[j]*particles.template getProp<velocity>(p)[x];
                    fy = fy + particles.template getProp<k_s>(p)[j] * (particles.template getProp<length>(p)[j] - particles.template getProp<n_length>(p)[j]) * y1 / xm; //-particles.template getProp<n_s>(p)[j]*particles.template getProp<velocity>(p)[y];
                    fz = fz + particles.template getProp<k_s>(p)[j] * (particles.template getProp<length>(p)[j] - particles.template getProp<n_length>(p)[j]) * z1 / xm; //-particles.template getProp<n_s>(p)[j]*particles.template getProp<velocity>(p)[z];
                    // count = count + 1;

                    // Adding viscoelastic forces

                    // scalar viscoelastic force
                    // scalar_viscoelastic = (particles.template getProp<velocity>(nnp)[x]-particles.template getProp<velocity>(p)[x])* x1/xm + (particles.template getProp<velocity>(nnp)[y]-particles.template getProp<velocity>(p)[y])* y1/xm + (particles.template getProp<velocity>(nnp)[z]-particles.template getProp<velocity>(p)[z])* z1/xm;

                    // fx=fx+particles.template getProp<n_s>(p)[j]*scalar_viscoelastic*x1/xm;
                    // fy=fy+particles.template getProp<n_s>(p)[j]*scalar_viscoelastic*y1/xm;
                    // fz=fz+particles.template getProp<n_s>(p)[j]*scalar_viscoelastic*z1/xm;

                    // fx=fx+particles.template getProp<n_s>(p)[j]*(particles.template getProp<velocity>(nnp)[x]-particles.template getProp<velocity>(p)[x])* x1/xm;
                    // fy=fy+particles.template getProp<n_s>(p)[j]*(particles.template getProp<velocity>(nnp)[y]-particles.template getProp<velocity>(p)[y])* y1/xm;
                    // fz=fz+particles.template getProp<n_s>(p)[j]*(particles.template getProp<velocity>(nnp)[z]-particles.template getProp<velocity>(p)[z])* z1/xm;
                }
            }

            // are we not updating the force on the particles?
            particles.template getProp<force>(p)[x] = fx;
            particles.template getProp<force>(p)[y] = fy;
            particles.template getProp<force>(p)[z] = fz;

            // std::cout << "t= " << t << "; " << particles.template getProp<pid>(p) << "; f= ";
            // std::cout << fx << " " << fy << " " << fz << "\n";

            ++it2;
        }

        // std::cout << "average f: " << f_av/392 << std::endl;
        // std::cout << "average displacement: " << d_av/392 << std::endl;

        // int save_step = 500/dt;

        // Uncomment this to outputt csv
        if ((ctr - 1) % save_csv == 0)
        {
            std::cout << "TIME step: " << t << std::endl;
            std::cout << "Avg movement: " << avg_movement << std::endl;
            particles.write_frame(dir + "Spring", ctr, CSV_WRITER);
        }

        // if ((ctr - 1)%save_vtk == 0) {
        //     particles.ghost_get<force>(SKIP_LABELLING);
        // particles.write_frame(dir + "Spring",ctr);
        // write_connection_vtk(dir + "connections",particles,ctr);
        // std::cout << "TIME in second: " << t.getwct() << std::endl;
        //}

        ctr += 1;
        t += dt;
    }

    particles.write_frame(dir + "final", 0, CSV_WRITER);

    std::ofstream outf("sim_output/sim_params.txt");

    outf << "Simulation completed with the following parameters:"
         << "\n";
    outf << "number of balls: " << n_particles << "\n";
    outf << "number of springs: " << n_springs << "\n";
    outf << "time: " << tf << "\n";
    outf << "dt: " << dt << "\n";
    outf << "mass: " << 1 << "\n";
    outf << "box dimensions: " << dim << "\n";
    outf << "running time: " << timer.getwct() << "\n";
    outf << "number of iterations: " << ctr << "\n";

    openfpm_finalize();
    return 0;
}
