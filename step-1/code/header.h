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
        HYPRE_Int size;
        double h_min;
        double l2_error;

        //Mesh objects
        ParMesh *pmesh;
        FiniteElementCollection *fec;
        ParFiniteElementSpace *fespace;

        //System objects
        ParBilinearForm *a;
        ParLinearForm *b;
        ParGridFunction *x;
        FunctionCoefficient *u;

        //Solver objects
        HypreParMatrix A;
        Vector B;
        Vector X;
};

extern double height;
extern double int_rad;
extern double out_rad;
