#ifndef IMHDFLUID
#define IMHDFLUID

#include <vector>

#include "rank3Tensor.hpp"

class imhdFluid {
	public:
		imhdFluid(size_t numVars, size_t N) : 
            variables_(numVars, rank3Tensor(N)), 
            Q_intermediate_(numVars, rank3Tensor(N)),
            xfluxes_(numVars, rank3Tensor(N)), 
            yfluxes_(numVars, rank3Tensor(N)), 
            zfluxes_(numVars, rank3Tensor(N)),
			int_xfluxes_(numVars, rank3Tensor(N)),
			int_yfluxes_(numVars, rank3Tensor(N)),
			int_zfluxes_(numVars, rank3Tensor(N)),
            numVars_(numVars),
            N_(N) {}
        
        // Fluid variable accessors
        double& imhdVar(size_t iv, size_t i, size_t j, size_t k) { return variables_[iv](i,j,k); }
        const double& imhdVar(size_t iv, size_t i, size_t j, size_t k) const { return variables_[iv](i,j,k); }

        double& rho(size_t i, size_t j, size_t k) { return variables_[0](i,j,k); }
        const double& rho(size_t i, size_t j, size_t k) const { return variables_[0](i,j,k); }
        
        double& rho_u(size_t i, size_t j, size_t k) { return variables_[1](i,j,k); }
        const double& rho_u(size_t i, size_t j, size_t k) const { return variables_[1](i,j,k); }
        
        double& rho_v(size_t i, size_t j, size_t k) { return variables_[2](i,j,k); }
        const double& rho_v(size_t i, size_t j, size_t k) const { return variables_[2](i,j,k); }

        double& rho_w(size_t i, size_t j, size_t k) { return variables_[3](i,j,k); }
        const double& rho_w(size_t i, size_t j, size_t k) const { return variables_[3](i,j,k); }

        double& Bx(size_t i, size_t j, size_t k) { return variables_[4](i,j,k); }
        const double& Bx(size_t i, size_t j, size_t k) const { return variables_[4](i,j,k); }
        
        double& By(size_t i, size_t j, size_t k) { return variables_[5](i,j,k); }
        const double& By(size_t i, size_t j, size_t k) const { return variables_[5](i,j,k); }

        double& Bz(size_t i, size_t j, size_t k) { return variables_[6](i,j,k); }
        const double& Bz(size_t i, size_t j, size_t k) const { return variables_[6](i,j,k); }

        double& e(size_t i, size_t j, size_t k) { return variables_[7](i,j,k); }
        const double& e(size_t i, size_t j, size_t k) const { return variables_[7](i,j,k); }

		// Fluid flux accessors
        // xfluxes
		double& xfluxes(size_t iv, size_t i, size_t j, size_t k) { return xfluxes_[iv](i,j,k); }
        const double& xfluxes(size_t iv, size_t i, size_t j, size_t k) const { return xfluxes_[iv](i,j,k); }
		
        double& xflux_rho(size_t i, size_t j, size_t k) { return xfluxes_[0](i,j,k); }
        const double& xflux_rho(size_t i, size_t j, size_t k) const { return xfluxes_[0](i,j,k); }
        
        double& xflux_rhou(size_t i, size_t j, size_t k) { return xfluxes_[1](i,j,k); }
        const double& xflux_rhou(size_t i, size_t j, size_t k) const { return xfluxes_[1](i,j,k); }

        double& xflux_rhov(size_t i, size_t j, size_t k) { return xfluxes_[2](i,j,k); }
        const double& xflux_rhov(size_t i, size_t j, size_t k) const { return xfluxes_[2](i,j,k); }

        double& xflux_rhow(size_t i, size_t j, size_t k) { return xfluxes_[3](i,j,k); }
        const double& xflux_rhow(size_t i, size_t j, size_t k) const { return xfluxes_[3](i,j,k); }

        double& xflux_Bx(size_t i, size_t j, size_t k) { return xfluxes_[4](i,j,k); }
        const double& xflux_Bx(size_t i, size_t j, size_t k) const { return xfluxes_[4](i,j,k); }

        double& xflux_By(size_t i, size_t j, size_t k) { return xfluxes_[5](i,j,k); }
        const double& xflux_By(size_t i, size_t j, size_t k) const { return xfluxes_[5](i,j,k); }

        double& xflux_Bz(size_t i, size_t j, size_t k) { return xfluxes_[6](i,j,k); }
        const double& xflux_Bz(size_t i, size_t j, size_t k) const { return xfluxes_[6](i,j,k); }

        double& xflux_e(size_t i, size_t j, size_t k) { return xfluxes_[7](i,j,k); }
        const double& xflux_e(size_t i, size_t j, size_t k) const { return xfluxes_[7](i,j,k); }

        // yfluxes
        double& yfluxes(size_t iv, size_t i, size_t j, size_t k) { return yfluxes_[iv](i,j,k); }
        const double& yfluxes(size_t iv, size_t i, size_t j, size_t k) const { return yfluxes_[iv](i,j,k); }
		
        double& yflux_rho(size_t i, size_t j, size_t k) { return yfluxes_[0](i,j,k); }
        const double& yflux_rho(size_t i, size_t j, size_t k) const { return yfluxes_[0](i,j,k); }
        
        double& yflux_rhou(size_t i, size_t j, size_t k) { return yfluxes_[1](i,j,k); }
        const double& yflux_rhou(size_t i, size_t j, size_t k) const { return yfluxes_[1](i,j,k); }

        double& yflux_rhov(size_t i, size_t j, size_t k) { return yfluxes_[2](i,j,k); }
        const double& yflux_rhov(size_t i, size_t j, size_t k) const { return yfluxes_[2](i,j,k); }

        double& yflux_rhow(size_t i, size_t j, size_t k) { return yfluxes_[3](i,j,k); }
        const double& yflux_rhow(size_t i, size_t j, size_t k) const { return yfluxes_[3](i,j,k); }

        double& yflux_Bx(size_t i, size_t j, size_t k) { return yfluxes_[4](i,j,k); }
        const double& yflux_Bx(size_t i, size_t j, size_t k) const { return yfluxes_[4](i,j,k); }

        double& yflux_By(size_t i, size_t j, size_t k) { return yfluxes_[5](i,j,k); }
        const double& yflux_By(size_t i, size_t j, size_t k) const { return yfluxes_[5](i,j,k); }

        double& yflux_Bz(size_t i, size_t j, size_t k) { return yfluxes_[6](i,j,k); }
        const double& yflux_Bz(size_t i, size_t j, size_t k) const { return yfluxes_[6](i,j,k); }

        double& yflux_e(size_t i, size_t j, size_t k) { return yfluxes_[7](i,j,k); }
        const double& yflux_e(size_t i, size_t j, size_t k) const { return yfluxes_[7](i,j,k); }

        // zfluxes
        double& zfluxes(size_t iv, size_t i, size_t j, size_t k) { return zfluxes_[iv](i,j,k); }
        const double& zfluxes(size_t iv, size_t i, size_t j, size_t k) const { return zfluxes_[iv](i,j,k); }
		
        double& zflux_rho(size_t i, size_t j, size_t k) { return zfluxes_[0](i,j,k); }
        const double& zflux_rho(size_t i, size_t j, size_t k) const { return zfluxes_[0](i,j,k); }
        
        double& zflux_rhou(size_t i, size_t j, size_t k) { return zfluxes_[1](i,j,k); }
        const double& zflux_rhou(size_t i, size_t j, size_t k) const { return zfluxes_[1](i,j,k); }

        double& zflux_rhov(size_t i, size_t j, size_t k) { return zfluxes_[2](i,j,k); }
        const double& zflux_rhov(size_t i, size_t j, size_t k) const { return zfluxes_[2](i,j,k); }

        double& zflux_rhow(size_t i, size_t j, size_t k) { return zfluxes_[3](i,j,k); }
        const double& zflux_rhow(size_t i, size_t j, size_t k) const { return zfluxes_[3](i,j,k); }

        double& zflux_Bx(size_t i, size_t j, size_t k) { return zfluxes_[4](i,j,k); }
        const double& zflux_Bx(size_t i, size_t j, size_t k) const { return zfluxes_[4](i,j,k); }

        double& zflux_By(size_t i, size_t j, size_t k) { return zfluxes_[5](i,j,k); }
        const double& zflux_By(size_t i, size_t j, size_t k) const { return zfluxes_[5](i,j,k); }

        double& zflux_Bz(size_t i, size_t j, size_t k) { return zfluxes_[6](i,j,k); }
        const double& zflux_Bz(size_t i, size_t j, size_t k) const { return zfluxes_[6](i,j,k); }

        double& zflux_e(size_t i, size_t j, size_t k) { return zfluxes_[7](i,j,k); }
        const double& zflux_e(size_t i, size_t j, size_t k) const { return zfluxes_[7](i,j,k); }
        
        // Intermediate variable accessors
		double& intermediateVar(size_t iv, size_t i, size_t j, size_t k) { return Q_intermediate_[iv](i,j,k); }
        const double& intermediateVar(size_t iv, size_t i, size_t j, size_t k) const { return Q_intermediate_[iv](i,j,k); }

        double& int_rho(size_t i, size_t j, size_t k) { return Q_intermediate_[0](i,j,k); }
        const double& int_rho(size_t i, size_t j, size_t k) const { return Q_intermediate_[0](i,j,k); }
        
        double& int_rhou(size_t i, size_t j, size_t k) { return Q_intermediate_[1](i,j,k); }
        const double& int_rhou(size_t i, size_t j, size_t k) const { return Q_intermediate_[1](i,j,k); }
        
        double& int_rhov(size_t i, size_t j, size_t k) { return Q_intermediate_[2](i,j,k); }
        const double& int_rhov(size_t i, size_t j, size_t k) const { return Q_intermediate_[2](i,j,k); }

        double& int_rhow(size_t i, size_t j, size_t k) { return Q_intermediate_[3](i,j,k); }
        const double& int_rhow(size_t i, size_t j, size_t k) const { return Q_intermediate_[3](i,j,k); }

        double& int_Bx(size_t i, size_t j, size_t k) { return Q_intermediate_[4](i,j,k); }
        const double& int_Bx(size_t i, size_t j, size_t k) const { return Q_intermediate_[4](i,j,k); }
        
        double& int_By(size_t i, size_t j, size_t k) { return Q_intermediate_[5](i,j,k); }
        const double& int_By(size_t i, size_t j, size_t k) const { return Q_intermediate_[5](i,j,k); }

        double& int_Bz(size_t i, size_t j, size_t k) { return Q_intermediate_[6](i,j,k); }
        const double& int_Bz(size_t i, size_t j, size_t k) const { return Q_intermediate_[6](i,j,k); }

        double& int_e(size_t i, size_t j, size_t k) { return Q_intermediate_[7](i,j,k); }
        const double& int_e(size_t i, size_t j, size_t k) const { return Q_intermediate_[7](i,j,k); }
		
        // Intermediate flux accessors
        // intermediate xfluxes
		double& int_xfluxes(size_t iv, size_t i, size_t j, size_t k) { return int_xfluxes_[iv](i,j,k); }
        const double& int_xfluxes(size_t iv, size_t i, size_t j, size_t k) const { return int_xfluxes_[iv](i,j,k); }
		
        double& int_xflux_rho(size_t i, size_t j, size_t k) { return int_xfluxes_[0](i,j,k); }
        const double& int_xflux_rho(size_t i, size_t j, size_t k) const { return int_xfluxes_[0](i,j,k); }
        
        double& int_xflux_rhou(size_t i, size_t j, size_t k) { return int_xfluxes_[1](i,j,k); }
        const double& int_xflux_rhou(size_t i, size_t j, size_t k) const { return int_xfluxes_[1](i,j,k); }

        double& int_xflux_rhov(size_t i, size_t j, size_t k) { return int_xfluxes_[2](i,j,k); }
        const double& int_xflux_rhov(size_t i, size_t j, size_t k) const { return int_xfluxes_[2](i,j,k); }

        double& int_xflux_rhow(size_t i, size_t j, size_t k) { return int_xfluxes_[3](i,j,k); }
        const double& int_xflux_rhow(size_t i, size_t j, size_t k) const { return int_xfluxes_[3](i,j,k); }

        double& int_xflux_Bx(size_t i, size_t j, size_t k) { return int_xfluxes_[4](i,j,k); }
        const double& int_xflux_Bx(size_t i, size_t j, size_t k) const { return int_xfluxes_[4](i,j,k); }

        double& int_xflux_By(size_t i, size_t j, size_t k) { return int_xfluxes_[5](i,j,k); }
        const double& int_xflux_By(size_t i, size_t j, size_t k) const { return int_xfluxes_[5](i,j,k); }

        double& int_xflux_Bz(size_t i, size_t j, size_t k) { return int_xfluxes_[6](i,j,k); }
        const double& int_xflux_Bz(size_t i, size_t j, size_t k) const { return int_xfluxes_[6](i,j,k); }

        double& int_xflux_e(size_t i, size_t j, size_t k) { return int_xfluxes_[7](i,j,k); }
        const double& int_xflux_e(size_t i, size_t j, size_t k) const { return int_xfluxes_[7](i,j,k); }

        // intermediate yfluxes
        double& int_yfluxes(size_t iv, size_t i, size_t j, size_t k) { return int_yfluxes_[iv](i,j,k); }
        const double& int_yfluxes(size_t iv, size_t i, size_t j, size_t k) const { return int_yfluxes_[iv](i,j,k); }
		
        double& int_yflux_rho(size_t i, size_t j, size_t k) { return int_yfluxes_[0](i,j,k); }
        const double& int_yflux_rho(size_t i, size_t j, size_t k) const { return int_yfluxes_[0](i,j,k); }
        
        double& int_yflux_rhou(size_t i, size_t j, size_t k) { return int_yfluxes_[1](i,j,k); }
        const double& int_yflux_rhou(size_t i, size_t j, size_t k) const { return int_yfluxes_[1](i,j,k); }

        double& int_yflux_rhov(size_t i, size_t j, size_t k) { return int_yfluxes_[2](i,j,k); }
        const double& int_yflux_rhov(size_t i, size_t j, size_t k) const { return int_yfluxes_[2](i,j,k); }

        double& int_yflux_rhow(size_t i, size_t j, size_t k) { return int_yfluxes_[3](i,j,k); }
        const double& int_yflux_rhow(size_t i, size_t j, size_t k) const { return int_yfluxes_[3](i,j,k); }

        double& int_yflux_Bx(size_t i, size_t j, size_t k) { return int_yfluxes_[4](i,j,k); }
        const double& int_yflux_Bx(size_t i, size_t j, size_t k) const { return int_yfluxes_[4](i,j,k); }

        double& int_yflux_By(size_t i, size_t j, size_t k) { return int_yfluxes_[5](i,j,k); }
        const double& int_yflux_By(size_t i, size_t j, size_t k) const { return int_yfluxes_[5](i,j,k); }

        double& int_yflux_Bz(size_t i, size_t j, size_t k) { return int_yfluxes_[6](i,j,k); }
        const double& int_yflux_Bz(size_t i, size_t j, size_t k) const { return int_yfluxes_[6](i,j,k); }

        double& int_yflux_e(size_t i, size_t j, size_t k) { return int_yfluxes_[7](i,j,k); }
        const double& int_yflux_e(size_t i, size_t j, size_t k) const { return int_yfluxes_[7](i,j,k); }

        // intermediate zfluxes
        double& int_zfluxes(size_t iv, size_t i, size_t j, size_t k) { return int_zfluxes_[iv](i,j,k); }
        const double& int_zfluxes(size_t iv, size_t i, size_t j, size_t k) const { return int_zfluxes_[iv](i,j,k); }
		
        double& int_zflux_rho(size_t i, size_t j, size_t k) { return int_zfluxes_[0](i,j,k); }
        const double& int_zflux_rho(size_t i, size_t j, size_t k) const { return int_zfluxes_[0](i,j,k); }
        
        double& int_zflux_rhou(size_t i, size_t j, size_t k) { return int_zfluxes_[1](i,j,k); }
        const double& int_zflux_rhou(size_t i, size_t j, size_t k) const { return int_zfluxes_[1](i,j,k); }

        double& int_zflux_rhov(size_t i, size_t j, size_t k) { return int_zfluxes_[2](i,j,k); }
        const double& int_zflux_rhov(size_t i, size_t j, size_t k) const { return int_zfluxes_[2](i,j,k); }

        double& int_zflux_rhow(size_t i, size_t j, size_t k) { return int_zfluxes_[3](i,j,k); }
        const double& int_zflux_rhow(size_t i, size_t j, size_t k) const { return int_zfluxes_[3](i,j,k); }

        double& int_zflux_Bx(size_t i, size_t j, size_t k) { return int_zfluxes_[4](i,j,k); }
        const double& int_zflux_Bx(size_t i, size_t j, size_t k) const { return int_zfluxes_[4](i,j,k); }

        double& int_zflux_By(size_t i, size_t j, size_t k) { return int_zfluxes_[5](i,j,k); }
        const double& int_zflux_By(size_t i, size_t j, size_t k) const { return int_zfluxes_[5](i,j,k); }

        double& int_zflux_Bz(size_t i, size_t j, size_t k) { return int_zfluxes_[6](i,j,k); }
        const double& int_zflux_Bz(size_t i, size_t j, size_t k) const { return int_zfluxes_[6](i,j,k); }

        double& int_zflux_e(size_t i, size_t j, size_t k) { return int_zfluxes_[7](i,j,k); }
        const double& int_zflux_e(size_t i, size_t j, size_t k) const { return int_zfluxes_[7](i,j,k); }
		
        // Misc inner products
        double v_dot_v(size_t i, size_t j, size_t k){
            return (variables_[1](i,j,k)*variables_[1](i,j,k) + variables_[2](i,j,k)*variables_[2](i,j,k) 
                + variables_[3](i,j,k)*variables_[3](i,j,k)) / (variables_[0](i,j,k) * variables_[0](i,j,k));
        }

        double B_dot_B(size_t i, size_t j, size_t k) {
            return variables_[4](i,j,k)*variables_[4](i,j,k) + variables_[5](i,j,k)*variables_[5](i,j,k) 
                    + variables_[6](i,j,k)*variables_[6](i,j,k);
        }

		double B_dot_v(size_t i, size_t j, size_t k){
            return 1.0 / (variables_[0](i,j,k)) * (variables_[1](i,j,k) * variables_[4](i,j,k) 
				+ variables_[2](i,j,k) * variables_[5](i,j,k)
                + variables_[3](i,j,k) * variables_[6](i,j,k));
        }

        double pressure(size_t i, size_t j, size_t k){
            return (gamma - 1)*(variables_[7](i,j,k) - variables_[0](i,j,k)*v_dot_v(i,j,k) / 2.0 
				- B_dot_B(i,j,k) / 2.0);
        }    

        const double getGamma() const { return gamma; }

        const size_t getNumVars() const { return numVars_; }

        const size_t getSideLen() const { return N_; }
	
	// Key for the fluid variables: 
    // 0 - rho; mass density 
    // 1 - rho_u; // x-momentum
    // 2 - rho_v; // y-momentum
    // 3 - rho_w; // z-momentum
    // 4 - B_x; // B-field in the x-direction
    // 5 - B_y; // B-field in the y-direction
    // 6 - B_z; // B-field in the z-direction
    // 7 - e; // Total energy
    private:
        size_t numVars_, N_;
        double gamma = 5.0 / 3.0; // polytropic index
        std::vector<rank3Tensor> variables_; // fluid variables
        std::vector<rank3Tensor> xfluxes_;
        std::vector<rank3Tensor> yfluxes_;
        std::vector<rank3Tensor> zfluxes_;
		std::vector<rank3Tensor> Q_intermediate_; // intermediate variables
		std::vector<rank3Tensor> int_xfluxes_;
		std::vector<rank3Tensor> int_yfluxes_;
		std::vector<rank3Tensor> int_zfluxes_;
};

#endif