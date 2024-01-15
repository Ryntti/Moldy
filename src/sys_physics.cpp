#include "header.h"


std::vector<atom> create_lattice(int n, double a, double &lx, double &ly) {
    /*
    Construct a two-dimensional n x n lattice of atoms at 0 K with lattice constant a and with atom id's ranging from 0 to n*n - 1. 
    Also creates null velocities and accelerations, as well as neighbor lists for each atom in the lattice.
    The atoms are summoned into an std::vector.
    */
    int size = n*n;
    std::vector <atom> r(size);
    int i, j, k;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            k = i*n + j;

            r[k].i = k;
            r[k].x = j*a;
            r[k].y = i*a;
            r[k].vx = 0.0;
            r[k].vy = 0.0;
            r[k].ax = 0.0;
            r[k].ay = 0.0;

            r[k].nlist.clear();
            // Next, add the indices of the neighbor atoms to the nlist of the ind:th atom.
            // Periodic boundary conditions must be taken into account.
            // Right neighbor (in case of leftmost edge i.e. j+1%n = 0, wrap around to the opposite edge)
            r[k].nlist.push_back(((j + 1)%n) + i*n );
            // left neighbor (in case of rightmost edge i.e. j+1%n = 1, wrap around to the opposite edge)
            r[k].nlist.push_back(((j + n - 1)%n) + i*n);
            // lower neighbor (in case of lower edge, wrap around to the upmost edge)
            r[k].nlist.push_back(j + ((i + 1)%n)*n);
            // Upper neighbor (In case of upper edge, wrap around to the bottom)
            r[k].nlist.push_back(j + ((i + n - 1)%n)*n);
        }   
    }
    lx = r[n-1].x - r[0].x + a; 
    ly = r[size-1].y - r[n-1].y + a;
    return r;
}

void apply_longitudinal_wave(std::vector<atom> &r, int n, double a, double amp, double kvec) {
    /*
    Displace atoms in a sinusoidal wave fashion
    */ 
    // double phi = 0.25*pi;       // An optional phase angle of pi/4
    for (int i = 0; i < r.size(); i++) {

        // Modify the x-coordinate with a sine pattern
        r[i].x += amp*sin(kvec*r[i].x);
    }
}

void displace_one_column(std::vector<atom> &r, int n, int col, double a) {
    /*
    Displaces the first column of atoms by 0.3*a in the positive x-axis direction
    */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) { 
            int ind = i*n + j;
            if (ind%n == 0) {
                r[ind].x += 0.3*a;
            }
        }
    }
}


void set_temp(std::vector<atom> &r, double t, int n, double m) {
    /*
    This function is inspired by the contents of lecture 3 page 15 of Antti Kuronen's
    lecture notes for the course FYS2088 Introduction to molecular dynamics simulations.
    */
    int i;
    int size = n*n;
    double vnorm;
    double vtotx, vtoty;
    std::random_device seed;
    std::mt19937 generator(seed());
    std::normal_distribution<double> distribution(-1, 1);

    vnorm = sqrt(kb*2*t / (m*amu_conversion));    // Units Å/fs                   
                     
    vtotx = 0.0;
    vtoty = 0.0;

    for (i = 0; i < size; i++) {
        r[i].vx = vnorm * distribution(generator);
        r[i].vy = vnorm * distribution(generator);
        vtotx += r[i].vx;
        vtoty += r[i].vy;
    }

    vtotx /= size;
    vtoty /= size;

    // lastly, subtract total velocity per atom from all the atoms's velocity components
    // i.e. remove the net momentum to eliminate any bulk movement due to possible initial bias in velocity sampling
    for (i = 0; i < size; i++) {
        r[i].vx -= vtotx;
        r[i].vy -= vtoty;
    }
}

void compute_temperature(std::vector<atom> &r, double &t, double &totkineng, int n, double m) {

    int i;
    int size = n*n;
    double k_tot = 0;       // k_tot is the ensemble average of the kinetic energy per atom

// Calculate k_tot
    // first, sum the squared velocities of atoms
    for (i = 0; i < size; i++) {
        k_tot += r[i].vx*r[i].vx + r[i].vy*r[i].vy;
    }
    // then, add in the prefactor 1/2 * m
    k_tot *= 0.5*m;
    // Convert from Daltons to fs^2 eV/Å^2
    k_tot *= amu_conversion;

    // update the value of total kinetic energy based on this. Units in eV
    totkineng = k_tot;

    // Now calculate temperature according to the equipartition theorem using the mean kinetic energy
    t = 2.0*k_tot / (3.0*size*kb);
}





void compute_energies(std::vector<atom> &r, int n, double &totpeng, double &totkineng, double &toteng) {
    /*
    Compute potential energy, kinetic energy and the total energy of the system
    */
    int i;
    double pot = 0.0;
    double kin = 0.0;

    for (i = 0; i < r.size(); i++) {
        pot += r[i].u;
        kin += r[i].ekin;
    }

    totpeng = pot;
    totkineng = kin;
    toteng = pot + kin;
}