/*
This file benchmarks the entire simulation with numerical diffusion added into the MacCormack Advance.
The parts of the project that it measures are the following functions:

(1) Construction of imhdFluid class by STL
(2) Writing of initial conditions
(3) MacCormack Advance (w/numerical diffusion)
(4) Periodic BCs
(5) Full loop

Performance metrics are obtained using calculations of the number of "floating-point" operations performed by each component.
*/
#define USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <unordered_map>
#include <variant>
#include <chrono>

#include "../benchmarking-include/rank3Tensor-bench.hpp"
#include "../benchmarking-include/imhdFluid-bench.hpp"
#include "../benchmarking-include/numerics-bench.hpp"
#include "../benchmarking-include/fluxFunctions-bench.hpp"

using std::string;
using std::ofstream;
using std::ifstream;
using std::endl;
using std::unordered_map;
using std::get;
using std::cerr;

using ParameterValue = std::variant<size_t, double, string>;
using timeunits = std::chrono::milliseconds;

unordered_map<string, ParameterValue> parseInputFile(const string& filename);
void computeExecutionStatistics(std::ofstream& log, const std::vector<timeunits> execTimes_MacAdv, const std::vector<timeunits> execTimes_BCs);

int main(){
    ofstream simlog;
    simlog.open("../benchmarking-build/bench.log");
    simlog << "Beginning simulation ...\n";

    // Parse input file
    simlog << "Parsing input file ...\n";
    unordered_map<string, ParameterValue> inputHash = parseInputFile("../benchmarking-build/bench.inp");

    size_t N = get<size_t>(inputHash["N"]);
    size_t Nt = get<size_t>(inputHash["Nt"]);
    double dx = get<double>(inputHash["dx"]);
    double dt = get<double>(inputHash["dt"]);
    double D = get<double>(inputHash["D"]);

    simlog << "Total data volume = " << 64*pow(N,3) << " doubles " << endl;

    double CFL = dt / dx;

    cartesianGrid ComputationalVolume = cartesianGrid(N);
    double L = N * dx; // side length

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

    auto start_sPSinit = std::chrono::high_resolution_clock::now();
    imhdFluid screwPinchSim = imhdFluid(8, N); 
    auto stop_sPSinit = std::chrono::high_resolution_clock::now();
    auto sPSinit_duration = std::chrono::duration_cast<timeunits>(stop_sPSinit - start_sPSinit);

    simlog << "Initializing 64 rank3Tensors of " << N << " elements per side took: " << sPSinit_duration.count() << " ms" << endl;

    double gamma = screwPinchSim.getGamma(); 
    
    auto start_ICs = std::chrono::high_resolution_clock::now();
    InitialConditions(screwPinchSim, ComputationalVolume, L);
    auto stop_ICs = std::chrono::high_resolution_clock::now();
    auto ICs_duration = std::chrono::duration_cast<timeunits>(stop_ICs - start_ICs);

    simlog << "Writing Initial Conditions took: " << ICs_duration.count() << " ms" << endl; 

    std::vector<timeunits> executionTimes_MacCormack(Nt);
    std::vector<timeunits> executionTimes_BCs(Nt);

    auto full_start = std::chrono::high_resolution_clock::now();
    for (size_t it = 1; it < Nt+1; it++){
        auto start_MacAdv = std::chrono::high_resolution_clock::now();
        MacCormackAdvance(screwPinchSim,dt,dx,D);
        auto stop_MacAdv = std::chrono::high_resolution_clock::now();
        auto MacAdv_duration = std::chrono::duration_cast<timeunits>(stop_MacAdv - start_MacAdv);
        executionTimes_MacCormack[it-1] = MacAdv_duration;

        auto start_BCs = std::chrono::high_resolution_clock::now();
        PeriodicBCs(screwPinchSim);
        auto stop_BCs = std::chrono::high_resolution_clock::now();
        auto BCs_duration = std::chrono::duration_cast<timeunits>(stop_BCs - start_BCs);
        executionTimes_BCs[it-1] = BCs_duration;
    }
    auto full_stop = std::chrono::high_resolution_clock::now();

    auto full_duration = std::chrono::duration_cast<timeunits>(full_stop - full_start);
    simlog << "Time taken for " << Nt << " timesteps, with " << N << " elements per side: " << full_duration.count() << " ms " << endl; 
    return 0;
}

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

void computeExecutionStatistics(ofstream& log, const std::vector<timeunits> execTimes_MacAdv, const std::vector<timeunits> execTimes_BCs){

}