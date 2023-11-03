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
bool createHDF5File(const imhdFluid &imhdData, const cartesianGrid &gridData, string &filename);
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

    simlog << "dx = " << endl;
    simlog << "dt = " << endl;
    simlog << "CFL = " << endl;

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

// Intended usage is to clear ../data/ before every run
void clearDataDirectory(const string& directoryPath){
    fs::path directory(directoryPath);

    for (const auto& file : fs::directory_iterator(directory)){
        fs::remove(file);
    }
}