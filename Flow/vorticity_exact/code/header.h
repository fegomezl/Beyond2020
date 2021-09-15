#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include "mfem.hpp"

struct Config{
    //Constructor
    Config(bool master, int nproc);

    //Passing parameters
    bool master;
    int nproc;

    int order;
    int refinements;
    bool last;

    double T_f;
    double invDeltaT;
    double viscosity;
    double epsilon_eta;
};

using namespace std;
using namespace mfem;

class Artic_sea{
    public:
        Artic_sea(Config config);
        void run(const char *mesh_file);
        ~Artic_sea();
    private:
        void make_grid(const char *mesh_file);
        void assemble_system();
        void solve_system();
        void output_results();

        //Global parameters
        Config config;

        //Output parameters
        int dim;
        int serial_refinements;
        HYPRE_BigInt size;
        HYPRE_BigInt size_v;
        double h_min;

        //Mesh objects
        ParMesh *pmesh;

        FiniteElementCollection *fec;
        FiniteElementCollection *fec_v;

        ParFiniteElementSpace *fespace;
        ParFiniteElementSpace *fespace_v;

        //System objects
        ParGridFunction *w;
        ParGridFunction *psi;
        ParGridFunction *v;

        //Solver objects
        BlockVector X;
        BlockVector B;

        HypreParMatrix *M;
        HypreParMatrix *D;
        HypreParMatrix *C;
        HypreParMatrix *Ct;

        //Convergence analysis
        double error_w;
        double error_psi;
};

extern double height;
extern double int_rad;
extern double out_rad;

extern double epsilon_r;

extern double exact_w(const Vector &x);
extern double exact_psi(const Vector &x);
