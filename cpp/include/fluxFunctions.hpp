#ifndef FLUX_FUNCTION
#define FLUX_FUNCTION

#include "rank3Tensor.hpp"
#include "imhdFluid.hpp"

// Flux functions
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
                imhdFluid.xflux_rho(i,j,k) = rho_u; 
                
                // x-directed flux of x-momentum density
                imhdFluid.xflux_rhou(i,j,k) = (rho_u * rho_u) / rho - B_x * B_x + p + BdotB / 2.0;
                
                // x-directed flux of y-momentum density
                imhdFluid.xflux_rhov(i,j,k) = (rho_v * rho_u) / rho - B_y * B_x;
                
                // x-directed flux of z-momentum density
                imhdFluid.xflux_rhow(i,j,k) = (rho_w * rho_u) / rho - B_x * B_z;
                
                // x-directed flux of B_{x} 
                imhdFluid.xflux_Bx(i,j,k) = 0.0; 
                
                // x-directed flux of B_{y}
                imhdFluid.xflux_By(i,j,k) = (rho_v / rho) * B_x - rho_u / rho * B_y;
                
                // x-directed flux of B_{z}
                imhdFluid.xflux_Bz(i,j,k) = rho_w / rho * B_x - rho_u / rho * B_z;
                
                // x-directed flux of energy  
                imhdFluid.xflux_e(i,j,k) = (rho_u / rho) * (e + p + BdotB / 2.0) - Bdotv*B_x;  
            }
        }
    }
}

void computefluxes_y(imhdFluid& imhdFluid){
    
}

void computefluxes_z(imhdFluid& imhdFluid){

}

// Intermediate flux functions
void int_computefluxes_x(imhdFluid& imhdFluid){

}

void int_computefluxes_y(imhdFluid& imhdFluid){
    
}

void int_computefluxes_z(imhdFluid& imhdFluid){
    
}
#endif