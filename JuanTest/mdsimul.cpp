#include "mdsimul.h"

// function implementations
void material_properties(std::vector<Particle> & balls)
{
  std::mt19937 gen; // random generator: Mersenne-Twister
  gen.seed(1); // control the seed
  std::uniform_real_distribution<> dis(0.5*RAD, 1.0*RAD); // radii between R/2 and R
  
  // set radii randomly
  for(auto & body : balls){
    body.rad=dis(gen);
    body.mass=RHO*4*M_PI*body.rad*body.rad*body.rad/3;
  }
}

void initial_conditions(std::vector<Particle> & balls)
{
  // set up uniformly on a 2d grid

  // extract max rad
  double max_rad=-1.0e300;
  for(auto & body : balls){
    if (body.rad >= max_rad) {
      max_rad=body.rad;
    }
  }

  // put on the 2d grid
  const int NX=L/(2*RAD);
  for(int id=0; id<balls.size(); ++id){
    int ix=id%NX; int iy=id%NX;
    balls[id].Rx = max_rad + ix*max_rad + 2;
    balls[id].Ry = 2*max_rad + iy*max_rad;
    balls[id].Rz = 0.0;
  }

  // initial velocity could be random, for now everything going up
  for(auto & body : balls){
    body.Vy = +1.2347;
    body.Vx = -1.5;
  }
}
void compute_force(std::vector<Particle> & balls)
{
  for(auto & body : balls){
    // reset force
    body.Fx = body.Fy = body.Fz = 0.0;

    // add gravitational force
    body.Fy += body.mass*G; // G is already negative

    // add force with walls
    double deltayb = body.rad - body.Ry; //bottom
    double deltayt = body.rad + body.Ry- L; //top
    double deltaxl = body.rad -body.Rx; // left
    double deltaxr = body.rad +body.Rx-L; //right
    if (deltayb > 0) {
      body.Fy += K*deltayb;
      //body.Fy += K*deltayb -10.1*body.Vy; //using dispersion
    }
    if (deltayt > 0) {
      body.Fy -= K*deltayt;
    }    
    if (deltaxl > 0) {
      body.Fx += K*deltaxl;
    }
    if (deltaxr > 0) {
      body.Fx -= K*deltaxr;
    }
    
    // force with other particles? other walls?
  }
}

void start_integration(std::vector<Particle> & balls)
{
  for(auto & body : balls){
    body.Vx -= 0.5*body.Fx*DT/body.mass;
    body.Vy -= 0.5*body.Fy*DT/body.mass;
    body.Vz -= 0.5*body.Fz*DT/body.mass;
  }
}

void integrate(std::vector<Particle> & balls)
{
  // use simple leap-frog
  for(auto & body : balls){
    body.Vx += body.Fx*DT/body.mass;
    body.Vy += body.Fy*DT/body.mass;
    body.Vz += body.Fz*DT/body.mass;
    body.Rx += body.Vx*DT;
    body.Ry += body.Vy*DT;
    body.Rz += body.Vz*DT;
  }
}

void print_info(const std::vector<Particle> & balls, const double & time)
{
  // here we print something useful
  // for now is just a test
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  
  std::cout << time
            << "\t" << balls[0].Rx
            << "\t" << balls[0].Ry 
            << "\t" << balls[0].Rz 
            << "\t" << balls[0].Vx 
            << "\t" << balls[0].Vy 
            << "\t" << balls[0].Vz 
            << "\n";
}
