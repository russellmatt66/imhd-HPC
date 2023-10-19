#ifndef NUMERICS
#define NUMERICS

#include "rank3Tensor.hpp"
#include "imhdFluid.hpp"

void computefluxes_x(imhdFluid& imhdFluid) {
    double p, BdotB, Bdotv;
    double rho, rho_u, rho_v, rho_w, B_x, B_y, B_z, e;
    size_t N = imhdFluid.getSideLen();
    for (size_t k = 0; k < N; k++){ // This is the fastest way to scan through rank3Tensor class
        for (size_t i = 0; i < N; i++){
            for (size_t j = 0; j < N; j++){
                // Hoisting reduces the number of memory accesses each loop, and increases readability
                p = imhdFluid.pressure(i,j,k);
                BdotB = imhdFluid.B_dot_B(i,j,k);
                Bdotv = imhdFluid.B_dot_v(i,j,k);
                rho = imhdFluid.rho(i,j,k);
                rho_u = imhdFluid.rho_u(i,j,k);
                rho_v = imhdFluid.rho_v(i,j,k);
                rho_w = imhdFluid.rho_w(i,j,k);
                B_x = imhdFluid.Bx(i,j,k);
                B_y = imhdFluid.By(i,j,k);
                B_z = imhdFluid.Bz(i,j,k);
                e = imhdFluid.e(i,j,k);
                
                // x-directed flux of mass density
                xfluxes_[0](i,j,k) = rho_u; 
                
                // x-directed flux of x-momentum density
                xfluxes_[1](i,j,k) = (rho_u * rho_u) / rho - B_x * B_x + p + BdotB / 2.0;
                
                // x-directed flux of y-momentum density
                xfluxes_[2](i,j,k) = (rho_v * rho_u) / rho - B_y * B_x;
                
                // x-directed flux of z-momentum density
                xfluxes_[3](i,j,k) = (rho_w * rho_u) / rho - B_x * B_z;
                
                // x-directed flux of B_{x} 
                xfluxes_[4](i,j,k) = 0.0; 
                
                // x-directed flux of B_{y}
                xfluxes_[5](i,j,k) = (rho_v / rho) * B_x - rho_u / rho * B_y;
                
                // x-directed flux of B_{z}
                xfluxes_[6](i,j,k) = rho_w / rho * B_x - rho_u / rho * B_z;
                
                // x-directed flux of energy  
                xfluxes_[7](i,j,k) = (rho_u / rho) * (e + p + BdotB / 2.0) - Bdotv*B_x;  
            }
        }
    }
}

void computefluxes_y(imhdFluid& imhdFluid){
    
}

// Compute intermediate fluxes

#endif
