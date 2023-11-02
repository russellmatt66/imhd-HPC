#define _USE_MATH_DEFINES

#include <iostream>
#include <unordered_map>
#include <variant>
#include <fstream>
#include <typeinfo>
#include <iomanip>

using std::unordered_map;
using std::string;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::get;

using ParameterValue = std::variant<size_t, double, string>;

unordered_map<string, ParameterValue> parseInputFile(const string& filename, ofstream& outputFile);

int main() {
    ofstream simlog;
    simlog.open("../debug-build/parseInput.log");
    simlog << "Beginning simulation ...\n";

    // Parse input file
    simlog << "Parsing input file ...\n";
    unordered_map<string, ParameterValue> inputHash = parseInputFile("../debug-build/debug.inp", simlog);

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
    simlog << typeid(dx).name() << endl;
    simlog << std::setprecision(3) << dx << endl;
    simlog << "dt = " << dt << endl;
    simlog << "CFL = " << CFL << endl;
    return 0;
}

// Parse input file and put values into hash
unordered_map<string, ParameterValue> parseInputFile(const string& filename, ofstream& outputFile){
    unordered_map<string, ParameterValue> parameters;
    ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        cerr << "Error opening input file: " << filename << endl;
        return parameters;
    }

    outputFile << "Inside parseInputFile() " << endl;

    string line;
    while (getline(inputFile, line)) {
        
        outputFile << "Current line is " << line << endl; 

        size_t delimiterPos = line.find('=');
        if (delimiterPos != string::npos) {
            string paramName = line.substr(0,delimiterPos);
            string paramValueStr = line.substr(delimiterPos + 1);
            if (paramName == "N" || paramName == "Nt"){
                outputFile << "paramName is " << paramName << endl;
                outputFile << "paramValueStr is " << paramValueStr << endl;
                size_t paramValue = std::stoul(paramValueStr);
                outputFile << "paramValue is " << paramValue << endl;
                parameters[paramName] = paramValue;
                outputFile << paramName << " = " << stoul(paramValueStr) << endl; 
            } else if (paramName == "dx" || paramName == "dt"){
                outputFile << "paramName is " << paramName << endl;
                outputFile << "paramValueStr is " << paramValueStr << endl;
                double paramValue = std::stod(paramValueStr);
                outputFile << "paramValue is " << paramValue << endl;
                parameters[paramName] = paramValue;
                outputFile << paramName << " = " << paramValue << endl;
            }
        }
    }

    // for (const auto& [key, value] : parameters) {
    //     outputFile << "paramName: " << key << " paramValue: ";
        
    //     std::visit([&outputFile](const auto& paramValue) {
    //         outputFile << paramValue;
    //     }, value);

    //     outputFile << endl;
    // }

    inputFile.close();
    return parameters;
}