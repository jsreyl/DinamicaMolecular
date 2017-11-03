#include<iostream>
#include<cstdlib>
#include<vector>
#include"mdsimul.h"

int main(int argc, char **argv){

  std::vector<Particle> balls(N); // N=Number of particles

  // t=0
  material_properties(balls);
  initial_conditions(balls);
  compute_force(balls);
  start_integration(balls);
  print_info(balls, 0);

  // evolve
  for(int istep=1; istep<NSTEPS; ++istep){
    compute_force(balls);
    integrate(balls);
    print_info(balls, istep*DT); // time = istep*DT
  }
  
  return 0;
}
