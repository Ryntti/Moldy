#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <vector>
#include <cmath>
#include <iomanip>

#include "header.h"

#define pi 3.14159265358979
#define amu_conversion 103.642697   // Conversion factor from Da to fs^2 * eV/Å^2
#define kb 8.617e-5                 // Boltzmann constant in units eV / K


int main(int argc, char *argv[]) {
    double dt;
    double t;
    double a;
    double k;
    double m;
    double amp;
    double k_vec;
    double totkineng = 0.0;
    double toteng = 0.0;
    double totpoteng = 0.0;
    int n;
    int dumprate;
    int thermorate;
    int maxstep;
    std::string inputfname;
    std::string dname;
    std::string thermoname;
    std::vector<atom> r;
    double lx = 0.0;
    double ly = 0.0;
    double time = 0.0;
    int step = 0;

    thermorate = 10;
    dname = "dump/dump.xyz";
    thermoname = "thermo.log";

    if (argc == 1) {
        inputfname = "inputfile";
    } else {
        inputfname = argv[1];
    }

    // Initialize parameters to some default values
    n = 0;
    dt = 0.0;
    maxstep = 0;
    t = 0.0;
    a = 0.0;
    m = 0.0;
    k = 0.0;

    // Read the proper input parameters
    read_inputfile(inputfname, t, dt, a, m, maxstep, n, k, dumprate, k_vec, amp);

    // Check that the sim parameters are right
    std::cout << "Simulation input parameters: \n"
            << "T = " << t << "\n" 
            << "timestep = " << dt << "\n" 
            << "lattice constant = " << a << "\n" 
            << "mass = " << m << "\n"
            << "steps = " << maxstep << "\n"
            << "size = " << n << "\n"
            << "spring constant = " << k << "\n" 
            << "dump rate = " << dumprate << "\n"
            << "wave number = " << k_vec << "\n"
            << "amplitude = " << amp << "\n" << std::endl;
    
// ----- Initialize ---------
    r = create_lattice(n, a, lx, ly);
    apply_longitudinal_wave(r, n, a, amp, k_vec);

    // Track the x-pos of atom 28
    int test_ind = 28;

    std::cout << "lx = " << lx << "\n"
    << "ly = " << ly << "\n" << std::endl;

    write_data("data.initial", r, n, m, "data", time, step, lx, ly, test_ind);
    set_temp(r, t, n, m);
    compute_temperature(r, t, totkineng, n, m);

    // int col = 0;
    // displace_one_column(r, n, col, a);

// ------- Main loop ---------

    for (step = 0; step < maxstep; step++) {
        time = step*dt;
        velocity_verlet(r, n, lx, ly, k, a, m, dt, time);
        compute_temperature(r, t, totkineng, n, m);
        compute_energies(r, n, totpoteng, totkineng, toteng);

        if (step%dumprate == 0) {
            write_data(dname, r, n, m, "dump", time, step, lx, ly, test_ind);
            write_data("dump/x_dump.txt", r, n, m, "x_only", time, step, lx, ly, test_ind);
        } 
        if (step%thermorate == 0) {
            write_thermo(thermoname, t, time, step, totpoteng, totkineng, toteng);
        }
    
    }

// -------------------------------
}

