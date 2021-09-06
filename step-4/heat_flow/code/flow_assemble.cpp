#include "header.h"

//Boundary values
double boundary_w(const Vector &x);
double boundary_psi(const Vector &x);

Flow_Operator::Flow_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes, const Vector &Theta):
    config(config),
    fespace(fespace),
    block_true_offsets(3),
    ess_bdr_psi(attributes), ess_bdr_w(attributes),
    f(NULL), g(NULL),
    m(NULL), d(NULL), c(NULL), ct(NULL),
    M(NULL), D(NULL), C(NULL), Ct(NULL),
    hBlocks(2,2), blockCoeff(2,2), H(NULL),
    superlu(MPI_COMM_WORLD), SLU_A(NULL),
    psi(&fespace), w(&fespace), v(&fespace_v),
    theta(&fespace), theta_eta(&fespace), 
    psi_grad(&fespace_v), theta_dr(&fespace),
    r(r_f), r_inv_hat(dim, r_inv_hat_f),
    w_coeff(boundary_w), psi_coeff(boundary_psi), 
    grad(&fespace, &fespace_v),
    rot(dim, rot_f), Psi_grad(&psi_grad),
    rot_Psi_grad(rot, Psi_grad)
{
    //Create the block offsets
    block_true_offsets[0] = 0;
    block_true_offsets[1] = fespace.TrueVSize();
    block_true_offsets[2] = fespace.TrueVSize();
    block_true_offsets.PartialSum();
  
    //Initialize the corresponding vectors
    X.Update(block_true_offsets);
    B.Update(block_true_offsets);
  
    //Define essential boundary conditions
    //
    //                  1
    //            /------------\
    // (w,psi=0)  |            |
    //           2|            |3
    //            |            |
    //            \------------/
    //                  0
    //
    ess_bdr_w[0] = 1; ess_bdr_w[1] = 1;
    ess_bdr_w[2] = 1; ess_bdr_w[3] = 1;
  
    ess_bdr_psi[0] = 1; ess_bdr_psi[1] = 1;
    ess_bdr_psi[2] = 1; ess_bdr_psi[3] = 1;

    //Create gradient interpolator
    grad.AddDomainIntegrator(new GradientInterpolator);
    grad.Assemble();
    grad.Finalize();

    //Initialize solver
    hBlocks = NULL;
    blockCoeff = 1.;

    superlu.SetPrintStatistics(false);
    superlu.SetSymmetricPattern(true);
    superlu.SetColumnPermutation(superlu::PARMETIS);
    superlu.SetIterativeRefine(superlu::SLU_DOUBLE);

    //Set initial system
    SetParameters(Theta);
}

//Boundary values
double boundary_w(const Vector &x){
    return 0.;
}

double boundary_psi(const Vector &x){
    double vel = 0.;
    return -0.5*vel*pow(x(0), 2);
}
