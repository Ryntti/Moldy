#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <vector>
#include <cmath>
#include <iomanip>


#define pi 3.14159265358979
#define amu_conversion 103.642697   // Conversion factor from Da to fs^2 * eV/Å^2
#define kb 8.617e-5                 // Boltzmann constant in units eV / K


class atom{
    // This defines the atom class that is used in a molecular dynamics simulation. An atom has a unique id, and x and y components 
    // of position, velocity and acceleration. It also has a unique integer list of the ids of its nearest neighbors.

    // Constructors not used in the md code, they are included for clarity and illustration.
    public:
        // Attributes
        int i; 
        double x; 
        double y; 
        double vx; 
        double vy; 
        double ax; 
        double ay;
        std::vector<int> nlist;
        double u;
        double ekin;

        // Default constructor, creates the very first atom of the lattice i.e. its id, position, velocity, acceleration and neighbor list.
        atom() : i(0), x(0.0), y(0.0), vx(0.0), vy(0.0), ax(0.0), ay(0.0), nlist({2, 7, 6, 31}), u(0.0), ekin(0.0)  {}
        // Custom constructor: the attribute values taken in as constructor parameters
        atom(int ind, double a, double b, double c, double d, double e, double f, std::vector<int> neighbors, double g, double h) : i(ind), x(a), y(b), vx(c), vy(d), ax(e), ay(f), nlist(neighbors), u(g), ekin(h) {}
    private:
};

//------------------------------------------------------------------------------

std::vector<atom> create_lattice(int n, double a, double &lx, double &ly);

void apply_longitudinal_wave(std::vector<atom> &r, int n, double a, double amp, double kvec);

void displace_one_column(std::vector<atom> &r, int n, int col, double a);

void set_temp(std::vector<atom> &r, double t, int n, double m);

void compute_temperature(std::vector<atom> &r, double &t, double &totkineng, int n, double m);

void compute_energies(std::vector<atom> &r, int n, double &totpeng, double &totkineng, double &toteng);

void velocity_verlet(std::vector <atom> &r, int n, double lx, double ly, double k, double a, double m, double dt, double time);

void read_inputfile(const std::string &inpfname, double &t, double &dt, double &a, double &m, int &maxstep, int &n, double &k, int &dumprate, double &k_vec, double &amp);

void write_thermo(std::string thermoname, double t, double time, int iteration, double totpoteng, double totkineng, double toteng);

void write_data(std::string outputfname, std::vector<atom> &r, int n, double m, std::string format, double time, int iteration, double lx, double ly, int test_ind = 0);

