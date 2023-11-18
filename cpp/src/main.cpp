#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <string>
#include <experimental/filesystem> 
#include <cmath>
#include <unordered_map>
#include <variant>
#include <stdexcept>
#include <chrono>
#include <H5Cpp.h> 

#include "../include/rank3Tensor.hpp"
#include "../include/imhdFluid.hpp"
#include "../include/numerics.hpp"
#include "../include/fluxFunctions.hpp"

using std::string;
using std::ofstream;
using std::ifstream;
using std::endl;
using std::unordered_map;
using std::get;
using std::cerr;
namespace fs = std::experimental::filesystem;

using ParameterValue = std::variant<size_t, double, string>;
using timeunits = std::chrono::milliseconds;

// Function prototypes
bool createHDF5File(const imhdFluid &imhdData, const cartesianGrid &gridData, string &filename);
void clearDataDirectory(const string &directoryPath);
unordered_map<string, ParameterValue> parseInputFile(const string& filename);

// Main
int main(){
    ofstream simlog;
    simlog.open("../build/imhd.log");
    simlog << "Beginning simulation ...\n";

    // Parse input file
    simlog << "Parsing input file ...\n";
    unordered_map<string, ParameterValue> inputHash = parseInputFile("../build/imhd.inp");

    size_t N = get<size_t>(inputHash["N"]);
    size_t Nt = get<size_t>(inputHash["Nt"]);
    double dx = get<double>(inputHash["dx"]);
    double dt = get<double>(inputHash["dt"]);
    double D = get<double>(inputHash["D"]);

    simlog << "N = " << N << endl;
    simlog << "Nt = " << Nt << endl;

    // Initialize simulation
    simlog << "Initializing timestep, and grid spacing.\n";
    double CFL = dt / dx;

    simlog << "dx = " << dx << endl;
    simlog << "dt = " << dt << endl;
    simlog << "CFL = " << CFL << endl;

    simlog << "Initializing Cartesian grid.\n";
    cartesianGrid ComputationalVolume = cartesianGrid(N);
    double L = N * dx; // side length

    // Origin in the middle of the domain
    double x,y,z;
    for (size_t k = 0; k < ComputationalVolume.num_depth(); k++){
        for (size_t i = 0; i < ComputationalVolume.num_rows(); i++){
            for (size_t j = 0; j < ComputationalVolume.num_cols(); j++){
                x = L/(N-1)*i - L/2;
                y = L/(N-1)*j - L/2;
                z = L/(N-1)*k - L/2;
                ComputationalVolume(i,j,k) = cartesianPoint(x,y,z);
            }
        }
    }
    simlog << "Cartesian grid successfully initialized.\n";

    simlog << "Initializing Screw-pinch (Ideal MHD fluid).\n";
    // 8 rank-3 tensors containing the fluid variables
    // + 24 containing the fluid fluxes
    // + 8 rank-3 tensors containing the intermediate variables (required by MacCormack method)
    // + 24 containing the intermediate fluxes
    imhdFluid screwPinchSim = imhdFluid(8, N); 
    double gamma = screwPinchSim.getGamma(); 
    simlog << "Screw-pinch gamma is " << gamma << endl;

    // Initial Conditions
    InitialConditions(screwPinchSim, ComputationalVolume, L);
    simlog << "Screw-pinch successfully initialized.\n";

    // For creating HDF5 files and holding return value
    const string dataPath = "../data"; 
    string filePath; 
    bool fileFlag;

    // Clear data directory of all files before beginning simulation
    simlog << "Clearing data files from " << dataPath << endl;
    clearDataDirectory(dataPath);
    simlog << "Data files successfully cleared " << endl;

    // Write ICs
    simlog << "Writing Initial Conditions to data.\n";
    filePath = "../data/FluidVars_timestep_0.h5";
    fileFlag = createHDF5File(screwPinchSim,ComputationalVolume,filePath);
    simlog << "Initial Conditions successfully written to data.\n";

    // Fluid (MacCormack) Advance
    simlog << "Beginning time advance.\n";
    auto loop_start = std::chrono::high_resolution_clock::now();
    for (size_t it = 1; it < Nt+1; it++){
        simlog << "Beginning timestep " << it << endl;
        MacCormackAdvance(screwPinchSim,dt,dx,D);
        simlog << "Advance completed. Writing boundary conditions\n";
        PeriodicBCs(screwPinchSim);
        simlog << "Boundary conditions written.\n";
        filePath = dataPath + "/FluidVars" + "_timestep_" + std::to_string(it) + ".h5";  
        simlog << "Writing data to file " << filePath << endl;
        fileFlag = createHDF5File(screwPinchSim,ComputationalVolume,filePath);
        simlog << "Data successfully written out.\n";
    }
    auto loop_stop = std::chrono::high_resolution_clock::now();
    auto loop_time = std::chrono::duration_cast<timeunits>(loop_stop - loop_start); 
    simlog << "Simulation loop took " << loop_time.count() << " ms for " 
        << N << " side-elements, and " << Nt << " timesteps " << endl;

    // Close log and halt
    simlog << "Simulation halting successfully. Closing log and returning control.";
    simlog.close();
    return 0;
}

// HDF5 file creation with C++ API
bool createHDF5File(const imhdFluid &imhdData, const cartesianGrid &gridData, string &filename){
    // Create file with the given name and fluid variables
    try {
        // Create file
        H5::H5File datafile(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        
        // Declare dataspace dimensions and dataspace for fluid variables 
        size_t numVars = imhdData.getNumVars();
        size_t N = imhdData.getSideLen();

        hsize_t dims[3] = {N, N, N}; 
        rank3Tensor buffer(N);
        const double* dataPtr = buffer.get_storage().data();

        H5::DataSpace dspcFluidVars(3,dims);

        // Create group for fluid variables
        H5::Group grpFluidVars = datafile.createGroup("/fluidvars");

        // Write fluid variables to their respective datasets - a little spaghetti but it's worth it
        H5::DataSet dstFluidVars_rho = 
            grpFluidVars.createDataSet("rho", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);

        H5::DataSet dstFluidVars_rhou = 
            grpFluidVars.createDataSet("rho_u", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);

        H5::DataSet dstFluidVars_rhov = 
            grpFluidVars.createDataSet("rho_v", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);

        H5::DataSet dstFluidVars_rhow = 
            grpFluidVars.createDataSet("rho_w", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);

        H5::DataSet dstFluidVars_Bx = 
            grpFluidVars.createDataSet("B_x", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);

        H5::DataSet dstFluidVars_By = 
            grpFluidVars.createDataSet("B_y", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);

        H5::DataSet dstFluidVars_Bz = 
            grpFluidVars.createDataSet("B_z", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);

        H5::DataSet dstFluidVars_e = 
            grpFluidVars.createDataSet("e", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);
        
        H5::DataSet dstCartMesh_x = 
            grpFluidVars.createDataSet("cartesian_x", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);

        H5::DataSet dstCartMesh_y = 
            grpFluidVars.createDataSet("cartesian_y", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);
        
        H5::DataSet dstCartMesh_z = 
            grpFluidVars.createDataSet("cartesian_z", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);
        
        // Write the tensors one by one instead of in a continuous manner.
        for (size_t iv = 0; iv < numVars; iv++){
            for (size_t k = 0; k < N; k++){
                for (size_t i = 0; i < N; i++){
                    for (size_t j = 0; j < N; j++){
                        buffer(i,j,k) = imhdData.imhdVar(iv,i,j,k);
                    }
                }
            }
            switch (iv) {
                case 0: // rho
                    dstFluidVars_rho.write(dataPtr, H5::PredType::NATIVE_DOUBLE);
                    break;
                case 1: // rho_u
                    dstFluidVars_rhou.write(dataPtr, H5::PredType::NATIVE_DOUBLE);
                    break;
                case 2: // rho_v
                    dstFluidVars_rhov.write(dataPtr, H5::PredType::NATIVE_DOUBLE);
                    break;    
                case 3: // rho_w
                    dstFluidVars_rhow.write(dataPtr, H5::PredType::NATIVE_DOUBLE);
                    break;
                case 4: // rho_Bx
                    dstFluidVars_Bx.write(dataPtr, H5::PredType::NATIVE_DOUBLE);
                    break;
                case 5: // rho_By
                    dstFluidVars_By.write(dataPtr, H5::PredType::NATIVE_DOUBLE);
                    break;
                case 6: // rho_Bz
                    dstFluidVars_Bz.write(dataPtr, H5::PredType::NATIVE_DOUBLE);
                    break;
                case 7: // energy
                    dstFluidVars_e.write(dataPtr, H5::PredType::NATIVE_DOUBLE);
                    break;
            }
        }

        // Now write Cartesian grid - x
        for (size_t k = 0; k < N; k++){
            for (size_t i = 0; i < N; i++){
                for (size_t j = 0; j < N; j++){
                    buffer(i,j,k) = gridData(i,j,k).x();
                }
            }
        }
        dstCartMesh_x.write(dataPtr, H5::PredType::NATIVE_DOUBLE);

        // Now write Cartesian grid - y
        for (size_t k = 0; k < N; k++){
            for (size_t i = 0; i < N; i++){
                for (size_t j = 0; j < N; j++){
                    buffer(i,j,k) = gridData(i,j,k).y();
                }
            }
        }
        dstCartMesh_y.write(dataPtr, H5::PredType::NATIVE_DOUBLE);

        // Now write Cartesian grid - z
        for (size_t k = 0; k < N; k++){
            for (size_t i = 0; i < N; i++){
                for (size_t j = 0; j < N; j++){
                    buffer(i,j,k) = gridData(i,j,k).z();
                }
            }
        }
        dstCartMesh_z.write(dataPtr, H5::PredType::NATIVE_DOUBLE);

        datafile.close();

        return true;
    } catch (const H5::Exception& e) {
        cerr << "Error creating HDF5 file: " << e.getCDetailMsg() << endl;
        return false;
    }
}

// Parse input file and put values into hash
unordered_map<string, ParameterValue> parseInputFile(const string& filename){
    
    unordered_map<string, ParameterValue> parameters;
    ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        cerr << "Error opening input file: " << filename << endl;
        return parameters;
    }
    
    string line;
    while (getline(inputFile, line)) {
        size_t delimiterPos = line.find('=');
        if (delimiterPos != string::npos) {
            string paramName = line.substr(0,delimiterPos);
            string paramValueStr = line.substr(delimiterPos + 1);
            if (paramName == "N" || paramName == "Nt"){
                size_t paramValue = std::stoul(paramValueStr);
                parameters[paramName] = paramValue;
            } else if (paramName == "dx" || paramName == "dt" || paramName == "D"){
                double paramValue = std::stod(paramValueStr);
                parameters[paramName] = paramValue;
            }
        }
    }
    
    inputFile.close();
    return parameters;
}

// Intended usage is to clear ../data/ before every run
void clearDataDirectory(const string& directoryPath){
    fs::path directory(directoryPath);

    for (const auto& file : fs::directory_iterator(directory)){
        if (file.path().filename() == "README.md"){ }
        else {
        fs::remove(file);
        }
    }
}
