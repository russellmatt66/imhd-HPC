#ifndef IMHDFLUID
#define IMHDFLUID

#include <vector>

#include "rank3Tensor.hpp"

class imhdFluid {
	public:
		imhdFluid(size_t numVars, size_t N) : 
            variables_(numVars, rank3Tensor(N)), // fluid variables
            intVariables_(numVars, rank3Tensor(N)), // intermediate variables
            xFluxes_(numVars, rank3Tensor(N)), // fluid fluxes
            yFluxes_(numVars, rank3Tensor(N)), 
            zFluxes_(numVars, rank3Tensor(N)),
			int_xFluxes_(numVars, rank3Tensor(N)), // intermediate fluxes
			int_yFluxes_(numVars, rank3Tensor(N)),
			int_zFluxes_(numVars, rank3Tensor(N)),
            numVars_(numVars),
            N_(N) // number of elements to data cube side 
            {}
        
        // Spaghetti for accessing all the member variables beyond this point
        // In scientific software development, there is a tradeoff between elegance and correctness that must sometimes be made. 
        // It is far more important that the code be easy to debug, than that it be the cleverest, most compact set of statements ever written.  
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
		double& xfluxes(size_t iv, size_t i, size_t j, size_t k) { return xFluxes_[iv](i,j,k); }
        const double& xfluxes(size_t iv, size_t i, size_t j, size_t k) const { return xFluxes_[iv](i,j,k); }
		
        double& xflux_rho(size_t i, size_t j, size_t k) { return xFluxes_[0](i,j,k); }
        const double& xflux_rho(size_t i, size_t j, size_t k) const { return xFluxes_[0](i,j,k); }
        
        double& xflux_rhou(size_t i, size_t j, size_t k) { return xFluxes_[1](i,j,k); }
        const double& xflux_rhou(size_t i, size_t j, size_t k) const { return xFluxes_[1](i,j,k); }

        double& xflux_rhov(size_t i, size_t j, size_t k) { return xFluxes_[2](i,j,k); }
        const double& xflux_rhov(size_t i, size_t j, size_t k) const { return xFluxes_[2](i,j,k); }

        double& xflux_rhow(size_t i, size_t j, size_t k) { return xFluxes_[3](i,j,k); }
        const double& xflux_rhow(size_t i, size_t j, size_t k) const { return xFluxes_[3](i,j,k); }

        double& xflux_Bx(size_t i, size_t j, size_t k) { return xFluxes_[4](i,j,k); }
        const double& xflux_Bx(size_t i, size_t j, size_t k) const { return xFluxes_[4](i,j,k); }

        double& xflux_By(size_t i, size_t j, size_t k) { return xFluxes_[5](i,j,k); }
        const double& xflux_By(size_t i, size_t j, size_t k) const { return xFluxes_[5](i,j,k); }

        double& xflux_Bz(size_t i, size_t j, size_t k) { return xFluxes_[6](i,j,k); }
        const double& xflux_Bz(size_t i, size_t j, size_t k) const { return xFluxes_[6](i,j,k); }

        double& xflux_e(size_t i, size_t j, size_t k) { return xFluxes_[7](i,j,k); }
        const double& xflux_e(size_t i, size_t j, size_t k) const { return xFluxes_[7](i,j,k); }

        // yfluxes
        double& yfluxes(size_t iv, size_t i, size_t j, size_t k) { return yFluxes_[iv](i,j,k); }
        const double& yfluxes(size_t iv, size_t i, size_t j, size_t k) const { return yFluxes_[iv](i,j,k); }
	
        double& yflux_rho(size_t i, size_t j, size_t k) { return yFluxes_[0](i,j,k); }
        const double& yflux_rho(size_t i, size_t j, size_t k) const { return yFluxes_[0](i,j,k); }
    
        double& yflux_rhou(size_t i, size_t j, size_t k) { return yFluxes_[1](i,j,k); }
        const double& yflux_rhou(size_t i, size_t j, size_t k) const { return yFluxes_[1](i,j,k); }
        
        double& yflux_rhov(size_t i, size_t j, size_t k) { return yFluxes_[2](i,j,k); }
        const double& yflux_rhov(size_t i, size_t j, size_t k) const { return yFluxes_[2](i,j,k); }

        double& yflux_rhow(size_t i, size_t j, size_t k) { return yFluxes_[3](i,j,k); }
        const double& yflux_rhow(size_t i, size_t j, size_t k) const { return yFluxes_[3](i,j,k); }
        
        double& yflux_Bx(size_t i, size_t j, size_t k) { return yFluxes_[4](i,j,k); }
        const double& yflux_Bx(size_t i, size_t j, size_t k) const { return yFluxes_[4](i,j,k); }
        
        double& yflux_By(size_t i, size_t j, size_t k) { return yFluxes_[5](i,j,k); }
        const double& yflux_By(size_t i, size_t j, size_t k) const { return yFluxes_[5](i,j,k); }

        double& yflux_Bz(size_t i, size_t j, size_t k) { return yFluxes_[6](i,j,k); }
        const double& yflux_Bz(size_t i, size_t j, size_t k) const { return yFluxes_[6](i,j,k); }

        double& yflux_e(size_t i, size_t j, size_t k) { return yFluxes_[7](i,j,k); }
        const double& yflux_e(size_t i, size_t j, size_t k) const { return yFluxes_[7](i,j,k); }

        // zfluxes
        double& zfluxes(size_t iv, size_t i, size_t j, size_t k) { return zFluxes_[iv](i,j,k); }
        const double& zfluxes(size_t iv, size_t i, size_t j, size_t k) const { return zFluxes_[iv](i,j,k); }
		
        double& zflux_rho(size_t i, size_t j, size_t k) { return zFluxes_[0](i,j,k); }
        const double& zflux_rho(size_t i, size_t j, size_t k) const { return zFluxes_[0](i,j,k); }
        
        double& zflux_rhou(size_t i, size_t j, size_t k) { return zFluxes_[1](i,j,k); }
        const double& zflux_rhou(size_t i, size_t j, size_t k) const { return zFluxes_[1](i,j,k); }

        double& zflux_rhov(size_t i, size_t j, size_t k) { return zFluxes_[2](i,j,k); }
        const double& zflux_rhov(size_t i, size_t j, size_t k) const { return zFluxes_[2](i,j,k); }

        double& zflux_rhow(size_t i, size_t j, size_t k) { return zFluxes_[3](i,j,k); }
        const double& zflux_rhow(size_t i, size_t j, size_t k) const { return zFluxes_[3](i,j,k); }

        double& zflux_Bx(size_t i, size_t j, size_t k) { return zFluxes_[4](i,j,k); }
        const double& zflux_Bx(size_t i, size_t j, size_t k) const { return zFluxes_[4](i,j,k); }

        double& zflux_By(size_t i, size_t j, size_t k) { return zFluxes_[5](i,j,k); }
        const double& zflux_By(size_t i, size_t j, size_t k) const { return zFluxes_[5](i,j,k); }

        double& zflux_Bz(size_t i, size_t j, size_t k) { return zFluxes_[6](i,j,k); }
        const double& zflux_Bz(size_t i, size_t j, size_t k) const { return zFluxes_[6](i,j,k); }

        double& zflux_e(size_t i, size_t j, size_t k) { return zFluxes_[7](i,j,k); }
        const double& zflux_e(size_t i, size_t j, size_t k) const { return zFluxes_[7](i,j,k); }
        
        // Intermediate variable accessors
		double& intermediateVar(size_t iv, size_t i, size_t j, size_t k) { return intVariables_[iv](i,j,k); }
        const double& intermediateVar(size_t iv, size_t i, size_t j, size_t k) const { return intVariables_[iv](i,j,k); }

        double& int_rho(size_t i, size_t j, size_t k) { return intVariables_[0](i,j,k); }
        const double& int_rho(size_t i, size_t j, size_t k) const { return intVariables_[0](i,j,k); }
        
        double& int_rhou(size_t i, size_t j, size_t k) { return intVariables_[1](i,j,k); }
        const double& int_rhou(size_t i, size_t j, size_t k) const { return intVariables_[1](i,j,k); }
        
        double& int_rhov(size_t i, size_t j, size_t k) { return intVariables_[2](i,j,k); }
        const double& int_rhov(size_t i, size_t j, size_t k) const { return intVariables_[2](i,j,k); }

        double& int_rhow(size_t i, size_t j, size_t k) { return intVariables_[3](i,j,k); }
        const double& int_rhow(size_t i, size_t j, size_t k) const { return intVariables_[3](i,j,k); }

        double& int_Bx(size_t i, size_t j, size_t k) { return intVariables_[4](i,j,k); }
        const double& int_Bx(size_t i, size_t j, size_t k) const { return intVariables_[4](i,j,k); }
        
        double& int_By(size_t i, size_t j, size_t k) { return intVariables_[5](i,j,k); }
        const double& int_By(size_t i, size_t j, size_t k) const { return intVariables_[5](i,j,k); }

        double& int_Bz(size_t i, size_t j, size_t k) { return intVariables_[6](i,j,k); }
        const double& int_Bz(size_t i, size_t j, size_t k) const { return intVariables_[6](i,j,k); }

        double& int_e(size_t i, size_t j, size_t k) { return intVariables_[7](i,j,k); }
        const double& int_e(size_t i, size_t j, size_t k) const { return intVariables_[7](i,j,k); }
		
        // Intermediate flux accessors
        // intermediate xfluxes
		double& int_xfluxes(size_t iv, size_t i, size_t j, size_t k) { return int_xFluxes_[iv](i,j,k); }
        const double& int_xfluxes(size_t iv, size_t i, size_t j, size_t k) const { return int_xFluxes_[iv](i,j,k); }
		
        double& int_xflux_rho(size_t i, size_t j, size_t k) { return int_xFluxes_[0](i,j,k); }
        const double& int_xflux_rho(size_t i, size_t j, size_t k) const { return int_xFluxes_[0](i,j,k); }
        
        double& int_xflux_rhou(size_t i, size_t j, size_t k) { return int_xFluxes_[1](i,j,k); }
        const double& int_xflux_rhou(size_t i, size_t j, size_t k) const { return int_xFluxes_[1](i,j,k); }

        double& int_xflux_rhov(size_t i, size_t j, size_t k) { return int_xFluxes_[2](i,j,k); }
        const double& int_xflux_rhov(size_t i, size_t j, size_t k) const { return int_xFluxes_[2](i,j,k); }

        double& int_xflux_rhow(size_t i, size_t j, size_t k) { return int_xFluxes_[3](i,j,k); }
        const double& int_xflux_rhow(size_t i, size_t j, size_t k) const { return int_xFluxes_[3](i,j,k); }

        double& int_xflux_Bx(size_t i, size_t j, size_t k) { return int_xFluxes_[4](i,j,k); }
        const double& int_xflux_Bx(size_t i, size_t j, size_t k) const { return int_xFluxes_[4](i,j,k); }

        double& int_xflux_By(size_t i, size_t j, size_t k) { return int_xFluxes_[5](i,j,k); }
        const double& int_xflux_By(size_t i, size_t j, size_t k) const { return int_xFluxes_[5](i,j,k); }

        double& int_xflux_Bz(size_t i, size_t j, size_t k) { return int_xFluxes_[6](i,j,k); }
        const double& int_xflux_Bz(size_t i, size_t j, size_t k) const { return int_xFluxes_[6](i,j,k); }

        double& int_xflux_e(size_t i, size_t j, size_t k) { return int_xFluxes_[7](i,j,k); }
        const double& int_xflux_e(size_t i, size_t j, size_t k) const { return int_xFluxes_[7](i,j,k); }

        // intermediate yfluxes
        double& int_yfluxes(size_t iv, size_t i, size_t j, size_t k) { return int_yFluxes_[iv](i,j,k); }
        const double& int_yfluxes(size_t iv, size_t i, size_t j, size_t k) const { return int_yFluxes_[iv](i,j,k); }
	
        double& int_yflux_rho(size_t i, size_t j, size_t k) { return int_yFluxes_[0](i,j,k); }
        const double& int_yflux_rho(size_t i, size_t j, size_t k) const { return int_yFluxes_[0](i,j,k); }
    
        double& int_yflux_rhou(size_t i, size_t j, size_t k) { return int_yFluxes_[1](i,j,k); }
        const double& int_yflux_rhou(size_t i, size_t j, size_t k) const { return int_yFluxes_[1](i,j,k); }

        double& int_yflux_rhov(size_t i, size_t j, size_t k) { return int_yFluxes_[2](i,j,k); }
        const double& int_yflux_rhov(size_t i, size_t j, size_t k) const { return int_yFluxes_[2](i,j,k); }

        double& int_yflux_rhow(size_t i, size_t j, size_t k) { return int_yFluxes_[3](i,j,k); }
        const double& int_yflux_rhow(size_t i, size_t j, size_t k) const { return int_yFluxes_[3](i,j,k); }

        double& int_yflux_Bx(size_t i, size_t j, size_t k) { return int_yFluxes_[4](i,j,k); }
        const double& int_yflux_Bx(size_t i, size_t j, size_t k) const { return int_yFluxes_[4](i,j,k); }

        double& int_yflux_By(size_t i, size_t j, size_t k) { return int_yFluxes_[5](i,j,k); }
        const double& int_yflux_By(size_t i, size_t j, size_t k) const { return int_yFluxes_[5](i,j,k); }

        double& int_yflux_Bz(size_t i, size_t j, size_t k) { return int_yFluxes_[6](i,j,k); }
        const double& int_yflux_Bz(size_t i, size_t j, size_t k) const { return int_yFluxes_[6](i,j,k); }

        double& int_yflux_e(size_t i, size_t j, size_t k) { return int_yFluxes_[7](i,j,k); }
        const double& int_yflux_e(size_t i, size_t j, size_t k) const { return int_yFluxes_[7](i,j,k); }

        // intermediate zfluxes
        double& int_zfluxes(size_t iv, size_t i, size_t j, size_t k) { return int_zFluxes_[iv](i,j,k); }
        const double& int_zfluxes(size_t iv, size_t i, size_t j, size_t k) const { return int_zFluxes_[iv](i,j,k); }
		
        double& int_zflux_rho(size_t i, size_t j, size_t k) { return int_zFluxes_[0](i,j,k); }
        const double& int_zflux_rho(size_t i, size_t j, size_t k) const { return int_zFluxes_[0](i,j,k); }
        
        double& int_zflux_rhou(size_t i, size_t j, size_t k) { return int_zFluxes_[1](i,j,k); }
        const double& int_zflux_rhou(size_t i, size_t j, size_t k) const { return int_zFluxes_[1](i,j,k); }

        double& int_zflux_rhov(size_t i, size_t j, size_t k) { return int_zFluxes_[2](i,j,k); }
        const double& int_zflux_rhov(size_t i, size_t j, size_t k) const { return int_zFluxes_[2](i,j,k); }

        double& int_zflux_rhow(size_t i, size_t j, size_t k) { return int_zFluxes_[3](i,j,k); }
        const double& int_zflux_rhow(size_t i, size_t j, size_t k) const { return int_zFluxes_[3](i,j,k); }

        double& int_zflux_Bx(size_t i, size_t j, size_t k) { return int_zFluxes_[4](i,j,k); }
        const double& int_zflux_Bx(size_t i, size_t j, size_t k) const { return int_zFluxes_[4](i,j,k); }

        double& int_zflux_By(size_t i, size_t j, size_t k) { return int_zFluxes_[5](i,j,k); }
        const double& int_zflux_By(size_t i, size_t j, size_t k) const { return int_zFluxes_[5](i,j,k); }

        double& int_zflux_Bz(size_t i, size_t j, size_t k) { return int_zFluxes_[6](i,j,k); }
        const double& int_zflux_Bz(size_t i, size_t j, size_t k) const { return int_zFluxes_[6](i,j,k); }

        double& int_zflux_e(size_t i, size_t j, size_t k) { return int_zFluxes_[7](i,j,k); }
        const double& int_zflux_e(size_t i, size_t j, size_t k) const { return int_zFluxes_[7](i,j,k); }
		
        // Miscellaneous inner products, pressure function, and accessors for the constants
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
    // 1 - rho_u; // x-momentum density
    // 2 - rho_v; // y-momentum density
    // 3 - rho_w; // z-momentum density
    // 4 - B_x; // B-field in the x-direction
    // 5 - B_y; // B-field in the y-direction
    // 6 - B_z; // B-field in the z-direction
    // 7 - e; // Total energy
    private:
        size_t numVars_, N_;
        double gamma = 5.0 / 3.0; // polytropic index
        std::vector<rank3Tensor> variables_; // fluid variables
        std::vector<rank3Tensor> xFluxes_;
        std::vector<rank3Tensor> yFluxes_;
        std::vector<rank3Tensor> zFluxes_;
		std::vector<rank3Tensor> intVariables_; // intermediate variables
		std::vector<rank3Tensor> int_xFluxes_;
		std::vector<rank3Tensor> int_yFluxes_;
		std::vector<rank3Tensor> int_zFluxes_;
};

#endif