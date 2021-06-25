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
        HYPRE_BigInt size_rt;
        HYPRE_BigInt size_l2;
        double h_min;
        double v_error;
        double p_error;

        //Mesh objects
        ParMesh *pmesh;
        FiniteElementCollection *fec_rt;
        FiniteElementCollection *fec_l2;
        ParFiniteElementSpace *fespace_rt;
        ParFiniteElementSpace *fespace_l2;
        Array<int> block_offsets;
        Array<int> block_true_offsets;

        //System objects
        BlockVector x;
        BlockVector b;
        ParLinearForm *f;
        ParLinearForm *g;
        ParBilinearForm *m;
        ParMixedBilinearForm *c;
        FunctionCoefficient p_exact_coeff;
        VectorFunctionCoefficient v_exact_coeff;

        //Solver objects
        BlockVector X;
        BlockVector B;
        HypreParMatrix *M;
        HypreParMatrix *C;
        TransposeOperator *Ct;
        BlockOperator *A;
        ParGridFunction *v;
        ParGridFunction *p;
};

extern double p_exact(const Vector &x);
extern void v_exact(const Vector &x, Vector &f);

extern double height;
extern double int_rad;
extern double out_rad;

extern double k;
