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
        HYPRE_BigInt size_w;
        HYPRE_BigInt size_psi;
        HYPRE_BigInt size_v;
        double h_min;
        double v_error;
        double p_error;

        //Mesh objects
        ParMesh *pmesh;

        FiniteElementCollection *fec_w;
        FiniteElementCollection *fec_psi;
        FiniteElementCollection *fec_v;

        ParFiniteElementSpace *fespace_w;
        ParFiniteElementSpace *fespace_psi;
        ParFiniteElementSpace *fespace_v;

        Array<int> block_offsets;
        Array<int> block_true_offsets;

        //System objects
        ParGridFunction *w;
        ParGridFunction *psi;
        ParGridFunction *v;
        BlockVector x;
        BlockVector b;
        ParLinearForm *g;
        ParLinearForm *f;
        ParBilinearForm *m;
        ParBilinearForm *d;
        ParMixedBilinearForm *c;

        //Solver objects
        BlockVector X;
        BlockVector B;
        HypreParMatrix *M;
        HypreParMatrix *D;
        HypreParMatrix *C;
};

extern double height;
extern double int_rad;
extern double out_rad;
