#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    transport_oper->SetParameters(X, *r_Velocity);
    ode_solver->Step(X, t, dt);

    flow_oper->SetParameters(X);
    flow_oper->Solve(Y, *Velocity, *r_Velocity);

    //Update visualization steps
    vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

    if (last || vis_steps <= vis_iteration){
        //Update parameters
        vis_iteration = 0;
        vis_impressions++;

        //Update information
        enthalphy->Distribute(X.GetBlock(0));
        salinity->Distribute(X.GetBlock(1));
        vorticity->Distribute(Y.GetBlock(0));
        stream->Distribute(Y.GetBlock(1));
        velocity->Distribute(Velocity);
        r_velocity->Distribute(r_Velocity);

        //Graph
        paraview_out->SetCycle(vis_impressions);
        paraview_out->SetTime(t);
        paraview_out->Save();
    }

    //Print the system state
    double percentage = 100*t/config.t_final;
    string progress = to_string((int)percentage)+"%";
    if (config.master){
        cout.precision(4);
        cout << left << setw(12)
             << iteration << setw(12)
             << dt << setw(12)
             << t  << setw(12)
             << progress << "\r";
        cout.flush();
    }
}

void Transport_Operator::SetParameters(const BlockVector &X, const Vector &r_Velocity){
    //Recover actual information
    ParGridFunction enthalphy(&fespace_H1);
    ParGridFunction salinity(&fespace_H1);
    ParGridFunction r_velocity(&fespace_ND);

    enthalphy.Distribute(X.GetBlock(0));
    salinity.Distribute(X.GetBlock(1));
    r_velocity.Distribute(r_Velocity);

    //Create the auxiliar grid functions
    ParGridFunction heat_diffusivity(&fespace_H1);
    ParGridFunction salt_diffusivity(&fespace_H1);
    ParGridFunction fusion_diffusivity(&fespace_H1);

    //Associate the values of each auxiliar function
    for (int ii = 0; ii < enthalphy.Size(); ii++){
        heat_diffusivity(ii) = HeatDiffusivity(enthalphy(ii));
        salt_diffusivity(ii) = SaltDiffusivity(enthalphy(ii));
        fusion_diffusivity(ii) = FusionDiffusivity(enthalphy(ii), salinity(ii));
    }

    //Set the grid coefficients
    GridFunctionCoefficient coeff_A0(&heat_diffusivity);
    GridFunctionCoefficient coeff_A1(&salt_diffusivity);
    GridFunctionCoefficient coeff_A2(&fusion_diffusivity);

    ProductCoefficient coeff_rA0(coeff_r, coeff_A0);
    ProductCoefficient coeff_rA1(coeff_r, coeff_A1);
    ProductCoefficient coeff_rA2(coeff_r, coeff_A2);

    GradientGridFunctionCoefficient coeff_grad_salinity(&salinity);
    ScalarVectorProductCoefficient coeff_fusion(coeff_rA2, coeff_grad_salinity);
    ScalarVectorProductCoefficient coeff_neg_fusion(-1., coeff_fusion);

    VectorGridFunctionCoefficient coeff_rV(&r_velocity);

    //Create transport matrix
    if (k0) delete k0;
    if (K0) delete K0;
    k0 = new ParBilinearForm(&fespace_H1);
    k0->AddDomainIntegrator(new DiffusionIntegrator(coeff_rA0));
    k0->AddDomainIntegrator(new ConvectionIntegrator(coeff_rV));
    k0->Assemble();
    k0->Finalize();
    K0 = k0->ParallelAssemble();    

    if (k1) delete k1;
    if (K1) delete K1;
    k1 = new ParBilinearForm(&fespace_H1);
    k1->AddDomainIntegrator(new DiffusionIntegrator(coeff_rA1));
    k1->AddDomainIntegrator(new ConvectionIntegrator(coeff_rV));
    k1->Assemble();
    k1->Finalize();
    K1 = k1->ParallelAssemble();    

    //Create RHS
    if (b0) delete b0;
    if (B0) delete B0;
    b0 = new ParLinearForm(&fespace_H1);
    b0->AddDomainIntegrator(new DomainLFGradIntegrator(coeff_fusion));
    b0->AddBoundaryIntegrator(new BoundaryTangentialLFIntegrator(coeff_neg_fusion));
    b0->Assemble();
    B0 = b0->ParallelAssemble();

    if (b1) delete b1;
    if (B1) delete B1;
    b1 = new ParLinearForm(&fespace_H1);
    ConstantCoefficient zero(0.);
    b1->AddDomainIntegrator(new DomainLFIntegrator(zero));
    b1->Assemble();
    B1 = b1->ParallelAssemble();
}

void Flow_Operator::SetParameters(const BlockVector &X){
    //Recover actual information
    ParGridFunction enthalphy(&fespace_H1);
    ParGridFunction salinity(&fespace_H1);

    enthalphy.Distribute(X.GetBlock(0));
    salinity.Distribute(X.GetBlock(1));

    //Create the auxiliar grid functions
    ParGridFunction inverse_permeability(&fespace_H1);
    //ParGridFunction expansivity_enthalphy(&fespace_H1);
    //ParGridFunction expansivity_salinity(&fespace_H1);i
    ParGridFunction buoyancy(&fespace_H1);
    //ParGridFunction enthalphy_dr(&fespace_H1);
    //ParGridFunction salinity_dr(&fespace_H1);
    ParGridFunction buoyancy_dr(&fespace_H1);

    //Associate the values of each auxiliar function
    for (int ii = 0; ii < enthalphy.Size(); ii++){
        inverse_permeability(ii) = InversePermeability(enthalphy(ii));
        //expansivity_enthalphy(ii) = ExpansivityEnthalphy(enthalphy(ii), salinity(ii));
        //expansivity_salinity(ii) = ExpansivitySalinity(enthalphy(ii), salinity(ii));
        buoyancy(ii) = Buoyancy(enthalphy(ii), salinity(ii));
    }

    //enthalphy.GetDerivative(1, 0, enthalphy_dr);
    //salinity.GetDerivative(1, 0, salinity_dr);
    buoyancy.GetDerivative(1, 0, buoyancy_dr);

    

    //Set the grid coefficients
    GridFunctionCoefficient coeff_gamma(&inverse_permeability);

    //GridFunctionCoefficient coeff_expansivity_enthalphy(&expansivity_enthalphy);
    //GridFunctionCoefficient coeff_expansivity_salinity(&expansivity_salinity);
    GridFunctionCoefficient coeff_buoyancy(&buoyancy_dr);
    //GridFunctionCoefficient coeff_enthalphy_dr(&enthalphy_dr);
    //GridFunctionCoefficient coeff_salinity_dr(&salinity_dr);

    ProductCoefficient coeff_neg_gamma(-1., coeff_gamma);
    ScalarVectorProductCoefficient coeff_neg_r_inv_gamma_hat(coeff_neg_gamma, coeff_r_inv_hat);

    //ProductCoefficient coeff_buoyancy_enthalphy(coeff_expansivity_enthalphy, coeff_enthalphy_dr);
    //ProductCoefficient coeff_buoyancy_salinity(coeff_expansivity_salinity, coeff_salinity_dr);
    //ProductCoefficient coeff_r_buoyancy_enthalphy(coeff_r, coeff_buoyancy_enthalphy);
    //ProductCoefficient coeff_r_buoyancy_salinity(coeff_r, coeff_buoyancy_salinity);
    ProductCoefficient coeff_r_buoyancy(coeff_r, coeff_buoyancy);
    
    //Define non-constant bilinear forms of the system
    if(A11) delete A11;
    ParBilinearForm a11(&fespace_H1);
    a11.AddDomainIntegrator(new DiffusionIntegrator(coeff_neg_gamma));
    a11.AddDomainIntegrator(new ConvectionIntegrator(coeff_neg_r_inv_gamma_hat));
    a11.Assemble();
    a11.Finalize();
    A11 = a11.ParallelAssemble();
    A11_e = A11->EliminateRowsCols(ess_tdof_1);

    //Define the non-constant RHS
    ParLinearForm b1(&fespace_H1);
    //b1.AddDomainIntegrator(new DomainLFIntegrator(coeff_r_buoyancy_enthalphy));
    //b1.AddDomainIntegrator(new DomainLFIntegrator(coeff_r_buoyancy_salinity));
    b1.AddDomainIntegrator(new DomainLFIntegrator(coeff_r_buoyancy));
    b1.Assemble();
    b1.ParallelAssemble(B1);
    A10_e->Mult(Vorticity, B1, -1., 1.);
    EliminateBC(*A11, *A11_e, ess_tdof_1, Stream, B1);
}
