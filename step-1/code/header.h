#pragma once

#include <iostream>
#include <fstream>
#include "mfem.hpp"

using namespace std;
using namespace mfem;

class Artic_sea{
    public:
        Artic_sea(bool master, int order, int refinements);
        void run(const char *mesh_file);
    private:
        void make_grid(const char *mesh_file);
        void assemble_system();
        void solve_system();
        void output_results();

        //Global parameters
        bool master;
        int order;
        int dim;
        int refinements;
        int serial_refinements;
        HYPRE_Int size;

        //Mesh objects
        ParMesh *pmesh;
        FiniteElementCollection *fec;
        ParFiniteElementSpace *fespace;

        //System objects
        ParBilinearForm *a;
        ParLinearForm *b;
        ParGridFunction *x;

        //Solver objects
        HypreParMatrix A;
        Vector B;
        Vector X;

        //Extra
        bool delete_fec;
};
