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
        HYPRE_Int size_w;
        HYPRE_Int size_psi;
        HYPRE_Int size_v;
        double h_min;
        double l2_error;

        //Mesh objects
        ParMesh *pmesh;

        Array<int> ess_tdof_list_w;
        Array<int> ess_tdof_list_psi;

        Array<int> ess_bdr_w;
        Array<int> ess_bdr_psi;

        FiniteElementCollection *fec_w;
        FiniteElementCollection *fec_psi;
        FiniteElementCollection *fec_v;

        ParFiniteElementSpace *fespace_w;
        ParFiniteElementSpace *fespace_psi;
        ParFiniteElementSpace *fespace_v;

        //System objects
        ParBilinearForm *a_w;
        ParLinearForm *b_w;

        ParBilinearForm *a_psi;
        ParLinearForm *b_psi;

        ParGridFunction *x_w;
        ParGridFunction *x_psi;
        ParGridFunction *x_v;

        FunctionCoefficient r;

        //Solver objects
        HypreParMatrix A_w;
        HypreParMatrix A_psi;

        HypreParVector *B_w;
        HypreParVector *B_psi;

        HypreParVector *X_w;
        HypreParVector *X_psi;
};

extern double rf(const Vector &x);

extern double exact(const Vector &x);

extern double height;
extern double int_rad;
extern double out_rad;
