#define _USE_MATH_DEFINES

#include<iostream>
#include<fstream>
#include<string>
#include<experimental/filesystem>
#include<cmath>
#include<unordered_map>
#include<variant>

#include<H5Cpp.h> 

#include "../debug-include/rank3Tensor-debug.hpp"
#include "../debug-include/imhdFluid-debug.hpp"
#include "../debug-include/numerics-debug.hpp"
#include "../debug-include/fluxFunctions-debug.hpp"

// Don't want to use "using namespace std;"
using std::string;
using std::ofstream;
using std::ifstream;
using std::endl;
using std::unordered_map;
using std::get;
using std::cerr;
namespace fs = std::experimental::filesystem;

using ParameterValue = std::variant<size_t, double, string>;

// Function prototypes
bool createHDF5File(const imhdFluid &imhdData, const cartesianGrid &gridData, const string &filename, ofstream& debuglog);
void clearDataDirectory(const string &directoryPath); 
unordered_map<string, ParameterValue> parseInputFile(const string& filename);

// Main
int main(){
    // Initialize log
    ofstream simlog;
    simlog.open("../debug-build/logs/debug-ICs.log");
    simlog << "Beginning simulation ...\n";
    
    // Parse input file
    simlog << "Parsing input file ...\n";
    unordered_map<string, ParameterValue> inputHash = parseInputFile("../debug-build/debug.inp");

    size_t N = get<size_t>(inputHash["N"]);
    size_t Nt = get<size_t>(inputHash["Nt"]);
    double dx = get<double>(inputHash["dx"]);
    double dt = get<double>(inputHash["dt"]);

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

    double gamma = screwPinchSim.getGamma(); // d_pinch = L / 2
    simlog << "Screw-pinch gamma is " << gamma << endl;

    // Initial Conditions
    ofstream IClog;
    IClog.open("../debug-build/logs/ICs.log");
    InitialConditions(screwPinchSim, ComputationalVolume, L, IClog);
    simlog << "Screw-pinch successfully initialized.\n";

    // For creating HDF5 files and holding return value
    const string dataPath = "../debug-data"; 
    string filePath; 
    bool fileFlag;

    // Clear data directory of all files before beginning simulation
    simlog << "Clearing data files from " << dataPath << endl;
    clearDataDirectory(dataPath);
    simlog << "Data files successfully cleared " << endl;

    // Write Initial Conditions to HDF5 file
    simlog << "Writing Initial Conditions to data.\n";
    filePath = "../debug-data/FluidVars_timestep_0.h5";
    ofstream HDF5log;
    HDF5log.open("../debug-build/logs/hdf5.log");
    fileFlag = createHDF5File(screwPinchSim,ComputationalVolume,filePath,HDF5log);
    simlog << "HDF5 file creation returned with boolean " << fileFlag << endl; 
    simlog << "Initial Conditions successfully written to data.\n";
    return 0;
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
            } else if (paramName == "dx" || paramName == "dt"){
                double paramValue = std::stod(paramValueStr);
                parameters[paramName] = paramValue;
            }
        }
    }
    
    inputFile.close();
    return parameters;
}

// HDF5 file creation with C++ API
bool createHDF5File(const imhdFluid &imhdData, const cartesianGrid &gridData, const string &filename, ofstream& debuglog){
    // Create file with the given name and fluid variables
    debuglog << "Inside createHDF5File()" << endl;
    try {
        // Create file
        debuglog << "Creating datafile" << endl;
        H5::H5File datafile(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        debuglog << "datafile object created" << endl;

        // Declare dataspace dimensions and dataspace for fluid variables 
        size_t numVars = imhdData.getNumVars();
        size_t N = imhdData.getSideLen();

        debuglog << "N is " << N << endl;
        debuglog << "numVars is " << numVars << endl;

        hsize_t dims[3] = {N, N, N}; 
        debuglog << "dims[0] = " << dims[0] << endl;
        debuglog << "dims[1] = " << dims[1] << endl;
        debuglog << "dims[2] = " << dims[2] << endl;


        debuglog << "Creating buffer object" << endl;
        /* Problem is here! */ 
        // double buffer[N][N][N];
        rank3Tensor buffer(N);
        const double *dataPtr = buffer.get_storage().data();
        debuglog << "Buffer object created" << endl;

        debuglog << "Creating dataspace object" << endl;
        H5::DataSpace dspcFluidVars(3,dims);
        debuglog << "Dataspace object created" << endl;

        // Create group for fluid variables
        debuglog << "Creating group object" << endl;        
        H5::Group grpFluidVars = datafile.createGroup("/fluidvars");
        debuglog << "Group object created" << endl;

        // Write fluid variables to their respective datasets - a little spaghetti but it's worth it
        debuglog << "Creating dataset for rho" << endl;
        H5::DataSet dstFluidVars_rho = 
            grpFluidVars.createDataSet("rho", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);
        debuglog << "Dataset object for rho created" << endl;
        
        debuglog << "Creating dataset for rho_u" << endl;
        H5::DataSet dstFluidVars_rhou = 
            grpFluidVars.createDataSet("rho_u", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);
        debuglog << "Dataset object for rho_u created" << endl;

        debuglog << "Creating dataset for rho_v" << endl;
        H5::DataSet dstFluidVars_rhov = 
            grpFluidVars.createDataSet("rho_v", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);
        debuglog << "Dataset object for rho_v created" << endl;

        debuglog << "Creating dataset for rho_w" << endl;
        H5::DataSet dstFluidVars_rhow = 
            grpFluidVars.createDataSet("rho_w", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);
        debuglog << "Dataset object for rho_w created" << endl;

        debuglog << "Creating dataset for Bx" << endl;
        H5::DataSet dstFluidVars_Bx = 
            grpFluidVars.createDataSet("B_x", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);
        debuglog << "Dataset object for Bx created" << endl;

        debuglog << "Creating dataset for By" << endl;
        H5::DataSet dstFluidVars_By = 
            grpFluidVars.createDataSet("B_y", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);
        debuglog << "Dataset object for By created" << endl;

        debuglog << "Creating dataset for Bz" << endl;
        H5::DataSet dstFluidVars_Bz = 
            grpFluidVars.createDataSet("B_z", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);
        debuglog << "Dataset object for Bz created" << endl;

        debuglog << "Creating dataset for e" << endl;
        H5::DataSet dstFluidVars_e = 
            grpFluidVars.createDataSet("e", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);
        debuglog << "Dataset object for e created" << endl;        
        
        debuglog << "Creating dataset for Cartesian Mesh - x" << endl;
        H5::DataSet dstCartMesh_x = 
            grpFluidVars.createDataSet("cartesian_x", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);
        debuglog << "Dataset object created for Cartesian Mesh - x" << endl;

        debuglog << "Creating dataset for Cartesian Mesh - y" << endl;
        H5::DataSet dstCartMesh_y = 
            grpFluidVars.createDataSet("cartesian_y", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);
        debuglog << "Dataset object created for Cartesian Mesh - x" << endl;
        
        debuglog << "Creating dataset for Cartesian Mesh - z" << endl;        
        H5::DataSet dstCartMesh_z = 
            grpFluidVars.createDataSet("cartesian_z", H5::PredType::NATIVE_DOUBLE, dspcFluidVars);
        debuglog << "Dataset object created for Cartesian Mesh - x" << endl;
        
        // Write the tensors one by one instead of in a continuous manner.
        for (size_t iv = 0; iv < numVars; iv++){
            for (size_t k = 0; k < N; k++){
                for (size_t i = 0; i < N; i++){
                    for (size_t j = 0; j < N; j++){
                        // buffer[k][i][j] = imhdData.imhdVar(iv,i,j,k);
                        buffer(i,j,k) = imhdData.imhdVar(iv,i,j,k);
                    }
                }
            }
            switch (iv) {
                case 0: // rho
                    debuglog << "Writing buffer for into rho dataset" << endl;
                    dstFluidVars_rho.write(dataPtr, H5::PredType::NATIVE_DOUBLE);
                    debuglog << "Buffer written into rho dataset" << endl;                    
                    break;
                case 1: // rho_u
                    debuglog << "Writing buffer for into rho_u dataset" << endl;                
                    dstFluidVars_rhou.write(dataPtr, H5::PredType::NATIVE_DOUBLE);
                    debuglog << "Buffer written into rho_u dataset" << endl;                    
                    break;
                case 2: // rho_v
                    debuglog << "Writing buffer for into rho_v dataset" << endl;
                    dstFluidVars_rhov.write(dataPtr, H5::PredType::NATIVE_DOUBLE);
                    debuglog << "Buffer written into rho_v dataset" << endl;                    
                    break;    
                case 3: // rho_w
                    debuglog << "Writing buffer for into rho_w dataset" << endl;
                    dstFluidVars_rhow.write(dataPtr, H5::PredType::NATIVE_DOUBLE);
                    debuglog << "Buffer written into rho_w dataset" << endl;                    
                    break;
                case 4: // rho_Bx
                    debuglog << "Writing buffer for into Bx dataset" << endl;
                    dstFluidVars_Bx.write(dataPtr, H5::PredType::NATIVE_DOUBLE);
                    debuglog << "Buffer written into Bx dataset" << endl;                    
                    break;
                case 5: // rho_By
                    debuglog << "Writing buffer for into By dataset" << endl;
                    dstFluidVars_By.write(dataPtr, H5::PredType::NATIVE_DOUBLE);
                    debuglog << "Buffer written into By dataset" << endl;                    
                    break;
                case 6: // rho_Bz
                    debuglog << "Writing buffer for into Bz dataset" << endl;
                    dstFluidVars_Bz.write(dataPtr, H5::PredType::NATIVE_DOUBLE);
                    debuglog << "Buffer written into Bz dataset" << endl;                    
                    break;
                case 7: // energy
                    debuglog << "Writing buffer for into e dataset" << endl;
                    dstFluidVars_e.write(dataPtr, H5::PredType::NATIVE_DOUBLE);
                    debuglog << "Buffer written into e dataset" << endl;                    
                    break;
            }
        }

        // Now write Cartesian grid - x
        for (size_t k = 0; k < N; k++){
            for (size_t i = 0; i < N; i++){
                for (size_t j = 0; j < N; j++){
                    // buffer[k][i][j] = gridData(i,j,k).x();
                    buffer(i,j,k) = gridData(i,j,k).x();
                }
            }
        }
        dstCartMesh_x.write(dataPtr, H5::PredType::NATIVE_DOUBLE);

        // Now write Cartesian grid - y
        for (size_t k = 0; k < N; k++){
            for (size_t i = 0; i < N; i++){
                for (size_t j = 0; j < N; j++){
                    // buffer[k][i][j] = gridData(i,j,k).y();
                    buffer(i,j,k) = gridData(i,j,k).y();
                }
            }
        }
        dstCartMesh_y.write(dataPtr, H5::PredType::NATIVE_DOUBLE);

        // Now write Cartesian grid - z
        for (size_t k = 0; k < N; k++){
            for (size_t i = 0; i < N; i++){
                for (size_t j = 0; j < N; j++){
                    // buffer[k][i][j] = gridData(i,j,k).z();
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

// Intended usage is to clear ../data/ before every run
void clearDataDirectory(const string& directoryPath){
    fs::path directory(directoryPath);

    for (const auto& file : fs::directory_iterator(directory)){
        fs::remove(file);
    }
}