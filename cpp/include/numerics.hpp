#ifndef NUMERICS
#define NUMERICS

#include <cmath>

#include "rank3Tensor.hpp"
#include "imhdFluid.hpp"
#include "fluxFunctions.hpp"

// Initial Conditions
void InitialConditions(imhdFluid& imhdFluid, const cartesianGrid& ComputationalVolume, const double L){
    double r, r_pinch = 0.5 * (L / 2), gamma = imhdFluid.getGamma();
    for (size_t k = 0; k < ComputationalVolume.num_depth(); k++){
        for (size_t i = 0; i < ComputationalVolume.num_rows(); i++){
            for (size_t j = 0; j < ComputationalVolume.num_cols(); j++){
                r = ComputationalVolume(i,j,k).r_cyl();
                if (r < r_pinch) { // Inside pinch
                    imhdFluid.rho(i,j,k) = 1.0; // single mass unit 
                    imhdFluid.Bz(i,j,k) = 1.0; 
                    imhdFluid.rho_w(i,j,k) = 1.0 - pow(r,2) / pow(r_pinch,2);
                    imhdFluid.e(i,j,k) = 1.0 / (gamma - 1.0) + 
                        0.5 * imhdFluid.rho(i,j,k) * imhdFluid.v_dot_v(i,j,k) + 
                        0.5 * imhdFluid.B_dot_B(i,j,k);
                }
                else {
                    imhdFluid.rho(i,j,k) = 0.01; // "vacuum"
                    imhdFluid.e(i,j,k) = 0.0 + 
                        0.5 * imhdFluid.rho(i,j,k) * imhdFluid.v_dot_v(i,j,k) + 
                        0.5 * imhdFluid.B_dot_B(i,j,k); // p_vac = 0.0
                }
            }
        }
    }
}

// In simulations of convection, diffusion is artificially introduced in order to stabilize high-frequency modes.
// Fourth-order centered difference on interior, forward and backward, respectively, at the edges.
// LOTS of room for optimization in this function, that will be for a later date. 
// Just note the places, and get an MVP. 
cartesianPoint NumericalDiffusion(const double D, const size_t iv, const size_t i, const size_t j, const size_t k, const imhdFluid& imhdFluid, const double dx){
    // iv - the fluid variable whose Laplacian is being approximated
    // This can be optimized, declaring these variables each time is not insignificant
    double Q_xx, Q_yy, Q_zz; 
    size_t N = imhdFluid.getSideLen(); 

    if (i == 0 || i == 1 || j == 0 || j == 1 || k == 0 || k == 1) { // boundary is two steps to the left => Forward difference 
        // Second-order Forward difference
        Q_xx = (1.0 / (2.0 * dx)) * (-3.0 * imhdFluid.imhdVar(iv,i,j,k) + 4.0 * imhdFluid.imhdVar(iv,i+1,j,k) - imhdFluid.imhdVar(iv,i+2,j,k));
        Q_yy = (1.0 / (2.0 * dx)) * (-3.0 * imhdFluid.imhdVar(iv,i,j,k) + 4.0 * imhdFluid.imhdVar(iv,i,j+1,k) - imhdFluid.imhdVar(iv,i,j+2,k));
        Q_zz = (1.0 / (2.0 * dx)) * (-3.0 * imhdFluid.imhdVar(iv,i,j,k) + 4.0 * imhdFluid.imhdVar(iv,i,j,k+1) - imhdFluid.imhdVar(iv,i,j,k+2));
    }
    else if (i == N-1 || i == N-2 || j == N-1 || j == N-2 || k == N-1 || k == N-2){ // boundary is two steps to the right => Backward difference
        // Second-order Backward difference
        Q_xx = (1.0 / (2.0 * dx)) * (3.0 * imhdFluid.imhdVar(iv,i,j,k) - 4.0 * imhdFluid.imhdVar(iv,i-1,j,k) + imhdFluid.imhdVar(iv,i-2,j,k));
        Q_yy = (1.0 / (2.0 * dx)) * (3.0 * imhdFluid.imhdVar(iv,i,j,k) - 4.0 * imhdFluid.imhdVar(iv,i,j-1,k) + imhdFluid.imhdVar(iv,i,j-2,k));
        Q_zz = (1.0 / (2.0 * dx)) * (3.0 * imhdFluid.imhdVar(iv,i,j,k) - 4.0 * imhdFluid.imhdVar(iv,i,j-1,k) + imhdFluid.imhdVar(iv,i,j,k-2));
    }
    else { // on interior => centered difference
        Q_xx = (1.0 / (12.0 * pow(dx,2))) * (-imhdFluid.imhdVar(iv,i+2,j,k) + 16.0 * imhdFluid.imhdVar(iv,i+1,j,k) - 30.0 * imhdFluid.imhdVar(iv,i,j,k)
            + 16.0 * imhdFluid.imhdVar(iv,i-1,j,k) - imhdFluid.imhdVar(iv,i-2,j,k));
        Q_yy = (1.0 / (12.0 * pow(dx,2))) * (-imhdFluid.imhdVar(iv,i,j+2,k) + 16.0 * imhdFluid.imhdVar(iv,i,j+1,k) - 30.0 * imhdFluid.imhdVar(iv,i,j,k)
            + 16.0 * imhdFluid.imhdVar(iv,i,j-1,k) - imhdFluid.imhdVar(iv,i,j-2,k));
        Q_zz = (1.0 / (12.0 * pow(dx,2))) * (-imhdFluid.imhdVar(iv,i,j,k+2) + 16.0 * imhdFluid.imhdVar(iv,i,j,k+1) - 30.0 * imhdFluid.imhdVar(iv,i,j,k)
            + 16.0 * imhdFluid.imhdVar(iv,i,j,k-1) - imhdFluid.imhdVar(iv,i,j,k-2));
    }
    return cartesianPoint(Q_xx, Q_yy, Q_zz);
}

// Numerical algorithm for performing the time advance 
// D is the numerical diffusion coefficient. 
void MacCormackAdvance(imhdFluid& imhdFluid, const double dt, const double dx, const double D){
    double dy = dx, dz = dx; // Can optimize this
    size_t N = imhdFluid.getSideLen(), numVars = imhdFluid.getNumVars();

    // Compute fluxes
    computefluxes_x(imhdFluid);
    computefluxes_y(imhdFluid);
    computefluxes_z(imhdFluid);

    // Compute intermediate variables
    for (size_t iv = 0; iv < numVars; iv++){
        for (size_t k = 0; k < N; k++){ 
            for (size_t i = 1; i < N-1; i++){ // handle walls separately, don't need to compute intermediate variables there
                for (size_t j = 1; j < N-1; j++){
                    if (k == 0) { // Periodic in Z - nature of equations makes this a corner case
                        imhdFluid.intermediateVar(iv,i,j,0) = imhdFluid.imhdVar(iv,i,j,0) 
                            - (dt / dx) * (imhdFluid.xfluxes(iv,i,j,0) - imhdFluid.xfluxes(iv,i-1,j,0))
                            - (dt / dy) * (imhdFluid.yfluxes(iv,i,j,0) - imhdFluid.yfluxes(iv,i,j-1,0)) 
                            - (dt / dz) * (imhdFluid.zfluxes(iv,i,j,0) - imhdFluid.zfluxes(iv,i,j,N-2));
                    }
                    else {
                        imhdFluid.intermediateVar(iv,i,j,k) = imhdFluid.imhdVar(iv,i,j,k) 
                            - (dt / dx) * (imhdFluid.xfluxes(iv,i,j,k) - imhdFluid.xfluxes(iv,i-1,j,k))
                            - (dt / dy) * (imhdFluid.yfluxes(iv,i,j,k) - imhdFluid.yfluxes(iv,i,j-1,k)) 
                            - (dt / dz) * (imhdFluid.zfluxes(iv,i,j,k) - imhdFluid.zfluxes(iv,i,j,k-1));
                    }

                }
            }
        }
    }

    // Compute intermediate fluxes
    int_computefluxes_x(imhdFluid);
    int_computefluxes_y(imhdFluid);
    int_computefluxes_z(imhdFluid);
    
    // Advance fluid variables on interior
    // Implement diffusion like an advection-diffusion equation is being solved
    for (size_t iv = 0; iv < numVars; iv++){
        for (size_t k = 0; k < N; k++){ 
            for (size_t i = 1; i < N-1; i++){ // handle walls separately, don't need to compute fluid variables there
                for (size_t j = 1; j < N-1; j++){
                    if (k == N-1) { // Periodic in Z - nature of equations makes this a corner case
                        imhdFluid.imhdVar(iv,i,j,N-1) = 0.5 * (imhdFluid.imhdVar(iv,i,j,N-1) - imhdFluid.intermediateVar(iv,i,j,N-1))
                            - 0.5 * (dt / dx) * (imhdFluid.int_xfluxes(iv,i+1,j,N-1) - imhdFluid.int_xfluxes(iv,i,j,N-1)) 
                            - 0.5 * (dt / dy) * (imhdFluid.int_yfluxes(iv,i,j+1,N-1) - imhdFluid.int_yfluxes(iv,i,j,N-1))
                            - 0.5 * (dt / dz) * (imhdFluid.int_zfluxes(iv,i,j,N-1) - imhdFluid.int_zfluxes(iv,i,j,1))
                            - NumericalDiffusion(D, iv, i, j, k, imhdFluid, dx).x() - NumericalDiffusion(D, iv, i, j, k, imhdFluid, dx).y()
                            - NumericalDiffusion(D, iv, i, j, k, imhdFluid, dx).z();
                    }
                    else {
                        imhdFluid.imhdVar(iv,i,j,k) = 0.5 * (imhdFluid.imhdVar(iv,i,j,k) - imhdFluid.intermediateVar(iv,i,j,k))
                            - 0.5 * (dt / dx) * (imhdFluid.int_xfluxes(iv,i+1,j,k) - imhdFluid.int_xfluxes(iv,i,j,k)) 
                            - 0.5 * (dt / dy) * (imhdFluid.int_yfluxes(iv,i,j+1,k) - imhdFluid.int_yfluxes(iv,i,j,k))
                            - 0.5 * (dt / dz) * (imhdFluid.int_zfluxes(iv,i,j,k) - imhdFluid.int_zfluxes(iv,i,j,k+1))
                            - NumericalDiffusion(D, iv, i, j, k, imhdFluid, dx).x() - NumericalDiffusion(D, iv, i, j, k, imhdFluid, dx).y()
                            - NumericalDiffusion(D, iv, i, j, k, imhdFluid, dx).z();
                    }

                }
            }
        }
    }
}


// Periodic Boundary Conditions - Separation of Concerns
void PeriodicBCs(imhdFluid& imhdFluid){
    // Boundary conditions:
    // (1) Periodic in z (excluding the perimeter of the xy-planes at k = k_{max} = N - 1)
    //      - Q(i,j,0) = Q(i,j,N-1)
    //      - Implement: Q(i,j,0) += Q(i,j,N-1); Q(i,j,N-1) = Q(i,j,0)
    // (2) Rigid, perfectly conducting wall surrounding the fusion plasma: 
    //  - 2 yz-planes: (0,y,z); (N-1,y,z);
    //  - 2 xz-planes: (x,0,z); (x,N-1,z);
    // Doing BCs like this should have negligible performance implications: O(N^2) work compared to O(N^3) when computing the interior points
    // Come back to this statement during benchmarking
    size_t N = imhdFluid.getSideLen(), numVars = imhdFluid.getNumVars();
    double gamma = imhdFluid.getGamma();
    
    // k = 0 and k = N-1 xy-planes
    for (size_t iv = 0; iv < numVars; iv++){
        for (size_t i = 0; i < N; i++){
            for (size_t j = 0; j < N; j++){
                imhdFluid.imhdVar(iv,i,j,0) += imhdFluid.imhdVar(iv,i,j,N-1);
                imhdFluid.imhdVar(iv,i,j,N-1) = imhdFluid.imhdVar(iv,i,j,0);
            }
        }
    }

    // Walls
    for (size_t k = 0; k < N; k++){
        // i = 0 and i = N-1 yz-planes
        for (size_t j = 0; j < N; j++){
            imhdFluid.rho(0,j,k) = 7.0; // "Lithium"
            imhdFluid.rho_u(0,j,k) = 0.0; // Rigid wall
            imhdFluid.rho_v(0,j,k) = 0.0;
            imhdFluid.rho_w(0,j,k) = 0.0;
            imhdFluid.Bx(0,j,k) = 0.0; // Perfectly-conducting
            imhdFluid.By(0,j,k) = 0.0;
            imhdFluid.Bz(0,j,k) = 0.0;
            imhdFluid.e(0,j,k) = imhdFluid.pressure(0,j,k) / (gamma - 1.0);
        
            imhdFluid.rho(N-1,j,k) = 7.0; // "Lithium"
            imhdFluid.rho_u(N-1,j,k) = 0.0; // Rigid wall
            imhdFluid.rho_v(N-1,j,k) = 0.0;
            imhdFluid.rho_w(N-1,j,k) = 0.0;
            imhdFluid.Bx(N-1,j,k) = 0.0; // Perfectly-conducting
            imhdFluid.By(N-1,j,k) = 0.0;
            imhdFluid.Bz(N-1,j,k) = 0.0;
            imhdFluid.e(N-1,j,k) = imhdFluid.pressure(N-1,j,k) / (gamma - 1.0);
        }
        // j = 0 and j = N-1 xz-planes
        for (size_t i = 0; i < N; i++){
            imhdFluid.rho(i,0,k) = 7.0; // "Lithium"
            imhdFluid.rho_u(i,0,k) = 0.0; // Rigid wall
            imhdFluid.rho_v(i,0,k) = 0.0;
            imhdFluid.rho_w(i,0,k) = 0.0;
            imhdFluid.Bx(i,0,k) = 0.0; // Perfectly-conducting
            imhdFluid.By(i,0,k) = 0.0;
            imhdFluid.Bz(i,0,k) = 0.0;
            imhdFluid.e(i,0,k) = imhdFluid.pressure(i,0,k) / (gamma - 1.0);
        
            imhdFluid.rho(i,N-1,k) = 7.0; // "Lithium"
            imhdFluid.rho_u(i,N-1,k) = 0.0; // Rigid wall
            imhdFluid.rho_v(i,N-1,k) = 0.0;
            imhdFluid.rho_w(i,N-1,k) = 0.0;
            imhdFluid.Bx(i,N-1,k) = 0.0; // Perfectly-conducting
            imhdFluid.By(i,N-1,k) = 0.0;
            imhdFluid.Bz(i,N-1,k) = 0.0;
            imhdFluid.e(i,N-1,k) = imhdFluid.pressure(i,N-1,k) / (gamma - 1.0);
        }
    }
}

#endif