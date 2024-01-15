#include "header.h"
 
void velocity_verlet(std::vector <atom> &r, int n, double lx, double ly, double k, double a, double m, double dt, double time) {
    /*
    Velocity Verlet solver. Takes in an std::vector of atom class objects as its entries by reference and modifies it and performs
    velocity verlet integration. In addition, it takes in the amount of atoms in one side of the lattice i.e. the n in the n x n lattice .
    */
    int i, j, size, n_ind;
    double r_x, r_y, v_x, v_y, dx, dy, f_x, f_y, dr;
    double dt2 = 0.5*dt*dt;
    double poteng;

    size = n*n;
    for (i = 0; i < size; i++) {
        // update the position of atom i
        r_x = r[i].x + dt*r[i].vx + dt2 * r[i].ax;
        r_y = r[i].y + dt*r[i].vy + dt2 * r[i].ay;

        // Adjust the positions according to the boundary conditions
        if (r_x >= lx) {
            r_x -= lx;
        }
        if (r_x < 0.0) {
            r_x += lx;
        }
        if (r_y >= ly) {
            r_y -= ly;
        }
        if (r_y < 0.0) {
            r_y += ly;
        }

        r[i].x = r_x;
        r[i].y = r_y;
    }


    for (i = 0; i < size; i++) {
        // prepare the first update to the velocity of atom i using the old acceleration
        v_x = r[i].vx + 0.5*dt*r[i].ax;
        v_y = r[i].vy + 0.5*dt*r[i].ay;

        f_x = 0.0;
        f_y = 0.0;
        r[i].u = 0.0;

        // Next, use the updated positions to calculate the new acceleration of atom i from the force
        for (j = 0; j < r[i].nlist.size(); j++) {
            // loop over all the neighbor atoms of atom i
            n_ind = r[i].nlist[j];

            dx = r[i].x - r[n_ind].x;
            dy = r[i].y - r[n_ind].y;

            if (dx > lx / 2) {
                dx -= lx;
            } 
            if (dx < -lx / 2) {
                dx += lx;
            } 
            if (dy > ly/2) {
                dy -= ly;
            } 
            if (dy < -ly/2) {
                dy += ly;
            }
 
            dr = sqrt(dx*dx + dy*dy);

            // Harmonic pair potential energy
            poteng = 0.5*k* (dr-a)*(dr-a);        // Potential energy between atom i and its neighbor n_ind. This shouldn't be double counted when the outer loop iterator i reaches the current n_ind value!
            // Eliminate double counting by division by 2, since the poteng of a pair interaction is calculated separately for every atom 
            poteng *= 0.5;
            r[i].u += poteng;

            // Sum the x and y components of force knowing dr, dx and dy and the harmonic interaction potential.
            f_x += -k*(dr - a) * dx/dr;
            f_y += -k*(dr - a) * dy/dr;        
        }
        // Update accelerations
        r[i].ax = (f_x/m)*0.00964853321;
        r[i].ay = (f_y/m)*0.00964853321;

        // Make the final update to the velocity of atom i using the updated acceleration
        r[i].vx = v_x + 0.5*dt * r[i].ax;
        r[i].vy = v_y + 0.5*dt * r[i].ay;

        // Update the kinetic energy of the atom
        r[i].ekin = r[i].vx*r[i].vx + r[i].vy*r[i].vy;
        r[i].ekin *= 0.5*m;

        // Lastly, convert Da to fs^2 eV/Å^2
        r[i].ekin *= amu_conversion;
    }
}