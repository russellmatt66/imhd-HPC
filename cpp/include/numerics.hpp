#ifndef NUMERICS
#define NUMERICS

#include "rank3Tensor.hpp"
#include "imhdFluid.hpp"
#include "fluxFunctions.hpp"

// Timestep
void MacCormackAdvance(imhdFluid& imhdFluid, const double dt, const double dx){
    double dy = dx, dz = dx;
    size_t N = imhdFluid.getSideLen(), numVars = imhdFluid.getNumVars();

    // Compute intermediate variables
    for (size_t iv = 0; iv < numVars; iv++){
        for (size_t k = 0; k < N-1; k++){
            for (size_t i = 1; i < N-1; i++){
                for (size_t j = 1; j < N-1; j++){

                }
            }
        }
    }
    // Compute intermediate fluxes

    // Advance fluid variables on interior


}

// Periodic Boundary Conditions - Separation of Concerns
void PeriodicBCs(imhdFluid& imhdFluid){
            // Boundary conditions:
            // (1) Periodic in z (excluding the perimeter of the xy-planes at k = k_{max} = N - 1)
            // (2) Rigid, perfectly conducting wall surrounding the fusion plasma: 
            //  - 2 yz-planes: (0,y,z); (N-1,y,z);
            //  - 2 xz-planes: (x,0,z); (x,N-1,z);
            // Doing BCs like this should have negligible performance implications: O(N^2) work compared to O(N^3) when computing the interior points
            // Caveat to the above: memory access pattern is not optimized for some of the walls 
            // Closing Remark: Performance degradation of ~75 MFlops observed compared to MacCormackAdvance
}

#endif