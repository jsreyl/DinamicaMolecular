#include "particle.h"
#include <vector>
#include <random>
#include <iostream>

// some constants for this simple simul
const double G=-9.81; // m/s^2
const double RHO=21.4; // kg/m^3
const double K=1230.45; // kg/m, lastic constant
const double L=3.45; // m, domain lenght
const double RAD=0.2*L; // m, radius for example particle
const int N=1; // Number of particles
const int NSTEPS=1000; // Number of time steps
const double DT=0.01;

// function declarations
void material_properties(std::vector<Particle> & balls);
void initial_conditions(std::vector<Particle> & balls);
void compute_force(std::vector<Particle> & balls);
void start_integration(std::vector<Particle> & balls);
void integrate(std::vector<Particle> & balls);
void print_info(const std::vector<Particle> & balls, const double & time);
