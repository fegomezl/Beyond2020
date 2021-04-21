#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "mfem.hpp"

using namespace std;
using namespace mfem;

struct Config{
    //Constructor
    Config(bool master, int nproc);
    //Passing parameters
    bool master;
    int nproc;
    int order;
    int refinements;
    double dt;
    double t_init;
    double t_final;
};

class Artic_sea{
    public:
        Artic_sea(Config config);
        void run(const char *mesh_file);
        ~Artic_sea();
    private:
        void make_grid(const char *mesh_file);
        void assemble_system();
        void time_step();
        void output_results();

        //Global parameters
        Config config;

        //Iteration parameters
        int iteration;
        double t;
        bool last;

        //Output parameters
        int dim;
        int serial_refinements;
        double h_min;
        HYPRE_Int size;

        //Mesh objects
        ParMesh *pmesh;
        FiniteElementCollection *fec;
        ParFiniteElementSpace *fespace;

        //System objects
        ParGridFunction *x;
        Vector X;
        FunctionCoefficient function;

        //Extra
        ParaViewDataCollection *paraview_out;
};

extern double Function(const Vector &x, double t);

extern double T_f;     //Fusion temperature
extern double T_i;     //Initial temperature

extern double alpha_l; //Liquid thermal conduction
extern double alpha_s; //Solid thermal conduction
extern double lambda;   //s(t) = sqrt(4*lamda*(alpha_s+alpha_l)*t)
