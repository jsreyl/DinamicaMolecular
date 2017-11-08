#include<iostream>
#include<cstdlib>
#include<vector>
#include"mdsimul.h"

int main(int argc, char **argv){

  std::vector<Particle> balls(N); // N=Number of particles

  // t=0
  material_properties(balls);
  //initial_conditions(balls);
  balls[0].Rx = 0.50;
  balls[0].Ry = 1.25;
  balls[1].Rx = 2.50;
  balls[1].Ry = 1.25;
  balls[0].Vx = +5.0;
  balls[1].Vx = -5.0;
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
