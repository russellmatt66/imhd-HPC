#ifndef NUMERICS
#define NUMERICS

#include "rank3Tensor.hpp"
#include "imhdFluid.hpp"
#include "fluxFunctions.hpp"

// Timestep
void MacCormackAdvance(imhdFluid& imhdFluid, const double dt, const double dx){
    double dy = dx, dz = dx;
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
                    if (k == 0) { // Periodic in Z
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
    for (size_t iv = 0; iv < numVars; iv++){
        for (size_t k = 0; k < N; k++){ 
            for (size_t i = 1; i < N-1; i++){ // handle walls separately, don't need to compute intermediate variables there
                for (size_t j = 1; j < N-1; j++){
                    if (k == N-1) {
                        imhdFluid.imhdVar(iv,i,j,N-1) = 0.5 * (imhdFluid.imhdVar(iv,i,j,N-1) - imhdFluid.intermediateVar(iv,i,j,N-1))
                            - 0.5 * (dt / dx) * (imhdFluid.int_xfluxes(iv,i+1,j,N-1) - imhdFluid.int_xfluxes(iv,i,j,N-1)) 
                            - 0.5 * (dt / dy) * (imhdFluid.int_yfluxes(iv,i,j+1,N-1) - imhdFluid.int_yfluxes(iv,i,j,N-1))
                            - 0.5 * (dt / dz) * (imhdFluid.int_zfluxes(iv,i,j,N-1) - imhdFluid.int_zfluxes(iv,i,j,1));
                    }
                    else {
                        imhdFluid.imhdVar(iv,i,j,k) = 0.5 * (imhdFluid.imhdVar(iv,i,j,k) - imhdFluid.intermediateVar(iv,i,j,k))
                            - 0.5 * (dt / dx) * (imhdFluid.int_xfluxes(iv,i+1,j,k) - imhdFluid.int_xfluxes(iv,i,j,k)) 
                            - 0.5 * (dt / dy) * (imhdFluid.int_yfluxes(iv,i,j+1,k) - imhdFluid.int_yfluxes(iv,i,j,k))
                            - 0.5 * (dt / dz) * (imhdFluid.int_zfluxes(iv,i,j,k) - imhdFluid.int_zfluxes(iv,i,j,k+1));
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
            // Come back to this during benchmarking

}

#endif