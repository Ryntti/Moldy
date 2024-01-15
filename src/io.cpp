#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <vector>
#include <cmath>
#include <iomanip>

#include "header.h"

void read_inputfile(const std::string &inpfname, double &t, double &dt, double &a, double &m, int &maxstep, int &n, double &k, int &dumprate, double &k_vec, double &amp) {
    /*
    Reads a text file and saves the simulation parameters therein to the variables given as arguments.
    */
    std::ifstream inputFile(inpfname);
    int i;

    if (inputFile.is_open()) {
        std::string line;
        i = 1;

        while (std::getline(inputFile, line)) {
            std::istringstream iss(line);
            std::string parname;
            std::string val;

            if (iss >> parname >> val) {
                try {
                    size_t pos;
                    int ival = std::stoi(val, &pos);
                
                if (pos == val.size()) {
                    // val converted to int, save as int
                    if (parname == "maxstep") {
                        maxstep = ival;
                    } else if (parname == "dumprate") {
                        dumprate = ival;
                    } else if (parname == "size") {
                        n = ival;
                    }
                } else {
                    // Not int, save as a double 
                    double dval = std::stod(val);

                    if (parname == "temperature") {
                        t = dval;
                    } else if (parname == "timestep") {
                        dt = dval;
                    } else if (parname == "mass") {
                        m = dval;
                    } else if (parname == "amplitude") {
                        amp = dval;
                    } else if (parname == "latticeconstant") {
                        a = dval;
                    } else if (parname == "springconstant") {
                        k = dval;
                    } else if (parname == "wavevector") {
                        k_vec = dval;
                    } else {
                        std::cerr << "Format error on line " << i << ": " << line << std::endl;
                    }
                }

                    
            } catch (const std::invalid_argument &e) {
                std::cerr << "Invalid argument: " << val << std::endl;
            }
            i += 1;   
        } else {
            std::cerr << "Format error on line " << i << ": " << line << std::endl;
        }
        i += 1;
        }
    }
    inputFile.close();

}


void write_thermo(std::string thermoname, double t, double time, int iteration, double totpoteng, double totkineng, double toteng) {
    /*
    Write thermodata, like temperature, potential energy etc. into a log file
    */
    std::ofstream outputFile;

    if (iteration == 0) {
        outputFile.open(thermoname, std::ofstream::trunc);
        if (outputFile.is_open()) {
            outputFile << "Iteration time temperature PotEng EKin totEng \n" 
            << iteration << " " << time << " " << t << " " << totpoteng << " " << totkineng << " " << toteng << std::endl;

        } else {
            std::cerr << "Error in opening thermo log file" << std::endl;
            return;
        }
    } else {
        outputFile.open(thermoname, std::ofstream::app);
        if (outputFile.is_open()) {
            outputFile << iteration << " " << time << " " << t << " " << totpoteng << " " << totkineng << " " << toteng << std::endl;
        } else {
            std::cerr << "Error in opening thermo log file" << std::endl;
            return;
        }
    }
    outputFile.close();
}

void write_data(std::string outputfname, std::vector<atom> &r, int n, double m, std::string format, double time, int iteration, double lx, double ly, int test_ind) {
    /*
    Write either dump or single frame data files in the extended xyz format. 
    */
    
    std::ofstream outputFile;
    std::string type;
    int size = n*n;
    int precision = 8;

    if (iteration == 0) {
        outputFile.open(outputfname, std::ofstream::trunc);
    } else {
        outputFile.open(outputfname, std::ofstream::app);
    }

    if (m == 63.5) {
        type = "Cu";
    } else {
        type = "NaN";
        std::cout << "Unexpected mass \n" << std::endl;
    }

    if (outputFile.is_open()) {
        
        m *= amu_conversion;
        if (format == "dump") {
            // Extended xyz format
            outputFile << size << "\n";
            outputFile << "Properties=species:S:1:id:I:1:pos:R:2:velo:R:2:force:R:2 Time=" << time << "\n";
            for (const auto &elem : r) {
                outputFile << type << " " << elem.i << " "
                << elem.x <<  " " << elem.y << " "
                << elem.vx << " " << elem.vy << " "
                << elem.ax*m << " " << elem.ay*m << " "
                << std::endl;
            }
        } else if (format == "data") {
            // Ext. xyz format again
            outputFile << size << "\n";
            outputFile << "Properties=species:S:1:id:I:1:pos:R:2 Time=" << time << "\n";
            for (const auto &elem : r) {
                outputFile << type << " " << elem.i << " " << elem.x << " " << elem.y << std::endl;
            }
        } else if (format == "x_only") {
            outputFile << time << " " << std::fixed << std::setprecision(precision) << r[test_ind].x << std::endl;
        }
    } else {
        std::cerr << "Error in opening output file" << std::endl;
        return;
    }
    outputFile.close();
}
