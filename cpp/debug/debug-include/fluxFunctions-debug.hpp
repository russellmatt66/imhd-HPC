#ifndef FLUX_FUNCTION_DEBUG
#define FLUX_FUNCTION_DEBUG

#include "rank3Tensor-debug.hpp"
#include "imhdFluid-debug.hpp"

// Flux functions
void computefluxes_x(imhdFluid& imhdFluid) {
    double p, BdotB, Bdotv;
    double rho, rho_u, rho_v, rho_w, B_x, B_y, B_z, e;
    size_t N = imhdFluid.getSideLen();
    
    // This nested loop ordering is the fastest way to scan through rank3Tensor class
    for (size_t k = 0; k < N; k++){ 
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
    double p, BdotB, Bdotv;
    double rho, rho_u, rho_v, rho_w, B_x, B_y, B_z, e;
    size_t N = imhdFluid.getSideLen();

    for (size_t k = 0; k < N; k++){ 
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
                
                // y-directed flux of mass density
                imhdFluid.yflux_rho(i,j,k) = rho_v;
                
                // y-directed flux of x-momentum density
                imhdFluid.yflux_rhou(i,j,k) = (rho_u * rho_v) / rho - B_x * B_y;
                
                // y-directed flux of y-momentum density
                imhdFluid.yflux_rhov(i,j,k) = (rho_v * rho_v / rho) - B_y * B_y + p + BdotB /2.0;
                
                // y-directed flux of z-momentum density
                imhdFluid.yflux_rhow(i,j,k) = (rho_w * rho_v / rho) - B_y * B_z;
                
                // y-directed flux of B_{x}
                imhdFluid.yflux_Bx(i,j,k) = (rho_u / rho) * B_y - (rho_v / rho) * B_x;
                
                // y-directed flux of B_{y}
                imhdFluid.yflux_By(i,j,k) = 0.0; // technically a waste to even calculate this - should be initialized to 0 anyways
                
                // y-directed flux of B_{z}
                imhdFluid.yflux_Bz(i,j,k) = (rho_w / rho) * B_y - (rho_v / rho) * B_z;
                
                // y-directed flux of energy 
                imhdFluid.yflux_e(i,j,k) = (e + p + BdotB / 2.0)*(rho_v / rho) - Bdotv * B_y;
            }
        }
    }
}    

void computefluxes_z(imhdFluid& imhdFluid){
    double p, BdotB, Bdotv;
    double rho, rho_u, rho_v, rho_w, B_x, B_y, B_z, e;
    size_t N = imhdFluid.getSideLen();

    for (size_t k = 0; k < N; k++){ 
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
                
                // z-directed flux of mass density 
                imhdFluid.zflux_rho(i,j,k) = rho_w;
                
                // z-directed flux of x-momentum density
                imhdFluid.zflux_rhou(i,j,k) = (rho_u * rho_w / rho) - (B_x * B_z);
                
                // z-directed flux of y-momentum density
                imhdFluid.zflux_rhov(i,j,k) = ((rho_v * rho_w) / rho) - B_y * B_z;
                
                // z-directed flux of z-momentum density
                imhdFluid.zflux_rhow(i,j,k) = (rho_w*rho_w) / rho - B_z*B_z + p + BdotB / 2.0;
                
                // z-directed flux of B_{x}
                imhdFluid.zflux_Bx(i,j,k) = (rho_u / rho)*B_z - (rho_w / rho) * B_x;
                
                // z-directed flux of B_{y}
                imhdFluid.zflux_By(i,j,k) = (rho_v / rho) * B_z - (rho_w / rho) * B_y;
                
                // z-directed flux of B_{z}
                imhdFluid.zflux_Bz(i,j,k) = 0.0; // a waste to even calculate this - should be initialized to 0 anyways
                
                // z-directed flux of energy
                imhdFluid.zflux_e(i,j,k) = (e + p + BdotB / 2.0)*(rho_w / rho) - Bdotv*B_z;
            }
        }
    }
}

// Intermediate flux functions
void int_computefluxes_x(imhdFluid& imhdFluid){
    double p, BdotB, Bdotv;
    double rho, rho_u, rho_v, rho_w, B_x, B_y, B_z, e;
    size_t N = imhdFluid.getSideLen();

    for (size_t k = 0; k < N; k++){ 
        for (size_t i = 0; i < N; i++){
            for (size_t j = 0; j < N; j++){
                // Hoisting reduces the number of memory accesses each loop, and increases readability
                p = imhdFluid.pressure(i,j,k);
                BdotB = imhdFluid.B_dot_B(i,j,k);
                Bdotv = imhdFluid.B_dot_v(i,j,k);
                rho = imhdFluid.int_rho(i,j,k);
                rho_u = imhdFluid.int_rhou(i,j,k);
                rho_v = imhdFluid.int_rhov(i,j,k);
                rho_w = imhdFluid.int_rhow(i,j,k);
                B_x = imhdFluid.int_Bx(i,j,k);
                B_y = imhdFluid.int_By(i,j,k);
                B_z = imhdFluid.int_Bz(i,j,k);
                e = imhdFluid.int_e(i,j,k);
                
                // x-directed flux of mass density
                imhdFluid.int_xflux_rho(i,j,k) = rho_u; 
                
                // x-directed flux of x-momentum density
                imhdFluid.int_xflux_rhou(i,j,k) = (rho_u * rho_u) / rho - B_x * B_x + p + BdotB / 2.0;
                
                // x-directed flux of y-momentum density
                imhdFluid.int_xflux_rhov(i,j,k) = (rho_v * rho_u) / rho - B_y * B_x;
                
                // x-directed flux of z-momentum density
                imhdFluid.int_xflux_rhow(i,j,k) = (rho_w * rho_u) / rho - B_x * B_z;
                
                // x-directed flux of B_{x} 
                imhdFluid.int_xflux_Bx(i,j,k) = 0.0; // technically a waste to even calculate this - should be initialized to 0 anyways
                
                // x-directed flux of B_{y}
                imhdFluid.int_xflux_By(i,j,k) = (rho_v / rho) * B_x - rho_u / rho * B_y;
                
                // x-directed flux of B_{z}
                imhdFluid.int_xflux_Bz(i,j,k) = rho_w / rho * B_x - rho_u / rho * B_z;
                
                // x-directed flux of energy  
                imhdFluid.int_xflux_e(i,j,k) = (rho_u / rho) * (e + p + BdotB / 2.0) - Bdotv*B_x;   
            }
        }
    }
}


void int_computefluxes_y(imhdFluid& imhdFluid){
    double p, BdotB, Bdotv;
    double rho, rho_u, rho_v, rho_w, B_x, B_y, B_z, e;
    size_t N = imhdFluid.getSideLen();

    for (size_t k = 0; k < N; k++){ 
        for (size_t i = 0; i < N; i++){
            for (size_t j = 0; j < N; j++){
                // Hoisting reduces the number of memory accesses each loop, and increases readability
                p = imhdFluid.pressure(i,j,k);
                BdotB = imhdFluid.B_dot_B(i,j,k);
                Bdotv = imhdFluid.B_dot_v(i,j,k);
                rho = imhdFluid.int_rho(i,j,k);
                rho_u = imhdFluid.int_rhou(i,j,k);
                rho_v = imhdFluid.int_rhov(i,j,k);
                rho_w = imhdFluid.int_rhow(i,j,k);
                B_x = imhdFluid.int_Bx(i,j,k);
                B_y = imhdFluid.int_By(i,j,k);
                B_z = imhdFluid.int_Bz(i,j,k);
                e = imhdFluid.int_e(i,j,k);
                
                // y-directed flux of mass density
                imhdFluid.int_yflux_rho(i,j,k) = rho_v;
                
                // y-directed flux of x-momentum density
                imhdFluid.int_yflux_rhou(i,j,k) = (rho_u * rho_v) / rho - B_x * B_y;
                
                // y-directed flux of y-momentum density
                imhdFluid.int_yflux_rhov(i,j,k) = (rho_v * rho_v / rho) - B_y * B_y + p + BdotB /2.0;
                
                // y-directed flux of z-momentum density
                imhdFluid.int_yflux_rhow(i,j,k) = (rho_w * rho_v / rho) - B_y * B_z;
                
                // y-directed flux of B_{x}
                imhdFluid.int_yflux_Bx(i,j,k) = (rho_u / rho) * B_y - (rho_v / rho) * B_x;
                
                // y-directed flux of B_{y}
                imhdFluid.int_yflux_By(i,j,k) = 0.0; // technically a waste to even calculate this - should be initialized to 0 anyways
                
                // y-directed flux of B_{z}
                imhdFluid.int_yflux_Bz(i,j,k) = (rho_w / rho) * B_y - (rho_v / rho) * B_z;
                
                // y-directed flux of energy 
                imhdFluid.int_yflux_e(i,j,k) = (e + p + BdotB / 2.0)*(rho_v / rho) - Bdotv * B_y;
            }
        }
    }
}    

void int_computefluxes_z(imhdFluid& imhdFluid){
    double p, BdotB, Bdotv;
    double rho, rho_u, rho_v, rho_w, B_x, B_y, B_z, e;
    size_t N = imhdFluid.getSideLen();

    for (size_t k = 0; k < N; k++){ 
        for (size_t i = 0; i < N; i++){
            for (size_t j = 0; j < N; j++){
                // Hoisting reduces the number of memory accesses each loop, and increases readability
                p = imhdFluid.pressure(i,j,k);
                BdotB = imhdFluid.B_dot_B(i,j,k);
                Bdotv = imhdFluid.B_dot_v(i,j,k);
                rho = imhdFluid.int_rho(i,j,k);
                rho_u = imhdFluid.int_rhou(i,j,k);
                rho_v = imhdFluid.int_rhov(i,j,k);
                rho_w = imhdFluid.int_rhow(i,j,k);
                B_x = imhdFluid.int_Bx(i,j,k);
                B_y = imhdFluid.int_By(i,j,k);
                B_z = imhdFluid.int_Bz(i,j,k);
                e = imhdFluid.int_e(i,j,k);
                
                // z-directed flux of mass density 
                imhdFluid.int_zflux_rho(i,j,k) = rho_w;
                
                // z-directed flux of x-momentum density
                imhdFluid.int_zflux_rhou(i,j,k) = (rho_u * rho_w / rho) - (B_x * B_z);
                
                // z-directed flux of y-momentum density
                imhdFluid.int_zflux_rhov(i,j,k) = ((rho_v * rho_w) / rho) - B_y * B_z;
                
                // z-directed flux of z-momentum density
                imhdFluid.int_zflux_rhow(i,j,k) = (rho_w*rho_w) / rho - B_z*B_z + p + BdotB / 2.0;
                
                // z-directed flux of B_{x}
                imhdFluid.int_zflux_Bx(i,j,k) = (rho_u / rho)*B_z - (rho_w / rho) * B_x;
                
                // z-directed flux of B_{y}
                imhdFluid.int_zflux_By(i,j,k) = (rho_v / rho) * B_z - (rho_w / rho) * B_y;
                
                // z-directed flux of B_{z}
                imhdFluid.int_zflux_Bz(i,j,k) = 0.0; // a waste to even calculate this - should be initialized to 0 anyways
                
                // z-directed flux of energy
                imhdFluid.int_zflux_e(i,j,k) = (e + p + BdotB / 2.0)*(rho_w / rho) - Bdotv*B_z;
            }
        }
    }
}
#endif