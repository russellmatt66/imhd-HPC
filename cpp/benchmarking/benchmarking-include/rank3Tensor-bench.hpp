#ifndef RANK3TENS_BENCH
#define RANK3TENS_BENCH

#include<cmath>
#include<vector>


class rank3Tensor {
	public:
		rank3Tensor(size_t N): num_rows_(N), num_cols_(N), num_depth_(N), storage_(num_rows_ * num_cols_ * num_depth_, 0.0) {}

		// row-major access order
		double& operator()(size_t i, size_t j, size_t k) { return storage_[k * num_rows_* num_cols_ + i * num_cols_ + j]; }
		const double& operator()(size_t i, size_t j, size_t k) const { return storage_[k * num_rows_* num_cols_ + i * num_cols_ + j]; }	
	
		const std::vector<double>& get_storage() { return storage_; }

		const size_t num_rows() const { return num_rows_; }
		const size_t num_cols() const { return num_cols_; }
		const size_t num_depth() const { return num_depth_; }

	private:
		size_t num_rows_, num_cols_, num_depth_;
		std::vector<double> storage_;
};

class cartesianPoint{
	public:
		cartesianPoint(double x, double y, double z): x_(x), y_(y), z_(z) {}

		double& x() { return x_; }
		const double& x() const { return x_; } 
		double& y() { return y_; }
		const double& y() const { return y_; } 
		double& z() { return z_; }
		const double& z() const { return z_; }

		const double r_cyl() const { return sqrt(x_ * x_ + y_ * y_); } 
	private:
		double x_, y_, z_;
};

class cartesianGrid{ 
    public:
        cartesianGrid(size_t N) : num_rows_(N), num_cols_(N), num_depth_(N), storage_(num_rows_ * num_cols_ * num_depth_, cartesianPoint(0.0,0.0,0.0)) {}
        
        // row-major order
        cartesianPoint& operator()(size_t i, size_t j, size_t k)    { return storage_[k * num_rows_ * num_cols_ + i * num_cols_ + j]; }
        const cartesianPoint& operator()(size_t i, size_t j, size_t k) const { return storage_[k * num_rows_ * num_cols_ + i * num_cols_ + j]; }

        size_t num_rows() const { return num_rows_; }
        size_t num_cols() const { return num_cols_; }
        size_t num_depth() const { return num_depth_; }
    private:
        size_t num_rows_, num_cols_, num_depth_;
        std::vector<cartesianPoint> storage_;
};

#endif
