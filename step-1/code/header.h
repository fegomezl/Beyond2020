#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include "mfem.hpp"

using namespace std;
using namespace mfem;

class Artic_sea{
    public:
        Artic_sea(bool master, int order, int refinements, bool last);
        void run(const char *mesh_file);
        ~Artic_sea();
    private:
        void make_grid(const char *mesh_file);
        void assemble_system();
        void solve_system();
        void output_results();

        //Global parameters
        bool master;
        bool last;
        int order;
        int dim;
        int refinements;
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

        //Extra
        bool delete_fec;
};

//Compute functions
double rhs(const Vector &x);        //Right hand side of the equation
double boundary(const Vector &x);   //for Neumann condition
double exact(const Vector &x);      //Exact solution for the equation

extern double height;
extern double int_rad;
extern double out_rad;
