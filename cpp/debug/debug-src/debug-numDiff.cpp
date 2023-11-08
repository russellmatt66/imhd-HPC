#include<iostream>
#include<fstream>
#include<string>
#include<experimental/filesystem>
#include<cmath>
#include<unordered_map>
#include<variant>

#include "../debug-include/rank3Tensor-debug.hpp"
#include "../debug-include/imhdFluid-debug.hpp"
#include "../debug-include/numerics-debug.hpp"
#include "../debug-include/fluxFunctions-debug.hpp"

using std::string;
using std::ofstream;
using std::ifstream;
using std::endl;
using std::unordered_map;
using std::get;
using std::cerr;

using ParameterValue = std::variant<size_t, double, string>;

unordered_map<string, ParameterValue> parseInputFile(const string& filename);

int main(){
    ofstream debuglog;
    debuglog.open("../debug-build/logs/debugNumDiff.log");
    debuglog << "Beginning simulation ...\n";

    // Parse input file
    debuglog << "Parsing input file ...\n";
    unordered_map<string, ParameterValue> inputHash = parseInputFile("../debug-build/debug.inp");

    size_t N = get<size_t>(inputHash["N"]);
    size_t Nt = get<size_t>(inputHash["Nt"]);
    double dx = get<double>(inputHash["dx"]);
    double dt = get<double>(inputHash["dt"]);
    double D = get<double>(inputHash["D"]);

    imhdFluid testFluid = imhdFluid(8,N);

    size_t numVars = testFluid.getNumVars();

    cartesianPoint debugOutput = cartesianPoint(0.0,0.0,0.0);

    for (size_t iv = 0; iv < numVars; iv++){
        for (size_t k = 0; k < N; k++){
            for (size_t i = 0; i < N; i++){
                for (size_t j = 0; j < N; j++){
                    debugOutput = NumericalDiffusion(debuglog,D,iv,i,j,k,testFluid,dx);
                }
            }
        }
    }
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
            } else if (paramName == "dx" || paramName == "dt" || paramName == "D"){
                double paramValue = std::stod(paramValueStr);
                parameters[paramName] = paramValue;
            }
        }
    }
    
    inputFile.close();
    return parameters;
}