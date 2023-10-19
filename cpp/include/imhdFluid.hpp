#ifndef IMHDFLUID
#define IMHDFLUID

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
		
		// Intermediate variable accessors
		
		// Intermediate flux accessors

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