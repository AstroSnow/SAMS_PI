/*
 *    Copyright 2025 SAMS Team
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/
#include "shared_data.h"

void simulation::controlvariables(simulationData &data) {

  data.nx=10; // Number of cells in the x-direction
  data.ny=2; // Number of cells in the y-direction
  data.nz=2; // Number of cells in the z-direction

  data.dt_multiplier = 0.8; // Default multiplier for time step
  data.dt=0.0;

  // Maximum number of iterations; if nsteps < 0, run until t_end
  data.nsteps = 2;
  data.t_end = 0.2; // One day in seconds

  // Geometry options: cartesian, cylindrical, spherical
  data.geometry = geometryType::Cartesian;

  // Domain limits
  data.x_min = 0.0;
  data.x_max = 1.0;
  data.y_min = -1.0;
  data.y_max = 1.0;
  data.z_min = -1.0;
  data.z_max = 1.0;

  // Boundary conditions
  data.xbc_min = BCType::BC_OTHER;
  data.xbc_max = BCType::BC_OTHER;
  data.ybc_min = BCType::BC_OTHER;
  data.ybc_max = BCType::BC_OTHER;
  data.zbc_min = BCType::BC_OTHER;
  data.zbc_max = BCType::BC_OTHER;

  // Grid stretching
  data.x_stretch = false;
  data.y_stretch = false;
  data.z_stretch = false;

  // Shock viscosity coefficients
  data.visc1 = 0.1;
  data.visc2 = 1.0;

  // Ratio of specific heat capacities
  data.gas_gamma = 1.4;

  // Average mass of an ion in proton masses
  data.mf = 1.2;

  // Resistive MHD options
  data.resistiveMHD = false;
  data.eta_background = 1.e-10;
  data.j_max = 1.0;
  data.eta0 = 2.e-10;

  // Remap kinetic energy correction
  data.rke = true;
  
  // Two-fluid flag
  data.two_fluid=true;

  // Output frequency and directory
  data.dt_snapshots = 0.02;
}

void simulation::initial_conditions(simulationData &data,simulationData &dataNeutral) {

  using Range = portableWrapper::Range;

  //std::cout << "Setting up initial conditions" << std::endl;
  // Set initial conditions for the simulation
  
  printf("Setting up initial conditions\n"); 
  
  portableWrapper::assign(data.vx,0.0);
  portableWrapper::assign(data.vy,0.0);
  portableWrapper::assign(data.vz,0.0);
  
  portableWrapper::assign(data.bx,0.0);
  portableWrapper::assign(data.by,0.0);
  portableWrapper::assign(data.bz,0.0);
  
  
  ////////////////////////////////////////////////////
  // Sod Shock tube
  T_dataType rho_L = 1.0;
  T_dataType P_L = 1.0;
  T_dataType vx_L = 0.0;
  
  T_dataType rho_R = 0.125;
  T_dataType P_R = 0.1;
  T_dataType vx_R = 0.0;
  ////////////////////////////////////////////////////

  portableWrapper::assign(data.rho,rho_R);
  portableWrapper::assign(data.energy_ion,P_R/rho_R/(data.gas_gamma-1.0));
  portableWrapper::assign(data.energy_electron,P_R/rho_R/(data.gas_gamma-1.0));

  //Some Neutral conditions
  if (data.two_fluid) {
    portableWrapper::assign(dataNeutral.vx,0.0);
    portableWrapper::assign(dataNeutral.vx,0.0);
    portableWrapper::assign(dataNeutral.vx,0.0);
    portableWrapper::assign(dataNeutral.rho,0.1);
    portableWrapper::assign(dataNeutral.energy_neutral,0.1);    
  }

  portableWrapper::applyKernel(
    LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) {
    
    if (data.xb(ix) < 0.5) {
      data.vx(ix, iy, iz) = vx_L;
      data.rho(ix, iy, iz) = rho_L;
      data.energy_ion(ix, iy, iz) = P_L/rho_L/(data.gas_gamma-1.0);
    } 

    //printf("%ld %f %f \n",ix,data.rho(ix,iy,iz),dataNeutral.rho(ix,iy,iz));

    },
    portableWrapper::Range(0, data.nx),
    portableWrapper::Range(0, data.ny),
    portableWrapper::Range(0, data.nz)
  );

//data.rho(1,0,0)=1.0;
//dataNeutral.rho(1,0,0)=2.0;
//printf("\n %f %f \n",data.rho(1,0,0),dataNeutral.rho(1,0,0));
//printf("\n%p",&data.rho(1,0,0));
//printf("\n%p \n",&dataNeutral.rho(1,0,0));

  std::cout << "Range of vx: " << portableWrapper::minval(data.vx) << " to " << portableWrapper::maxval(data.vx) << "\n";
  std::cout << "Range of vy: " << portableWrapper::minval(data.vy) << " to " << portableWrapper::maxval(data.vy) << "\n";
  std::cout << "Range of vz: " << portableWrapper::minval(data.vz) << " to " << portableWrapper::maxval(data.vz) << "\n";

  T_dataType bmult = 000.0;
  portableWrapper::assign(data.bx,0.00);
  portableWrapper::assign(data.by,0.00);
  portableWrapper::assign(data.bz,0.00);
  // Set the initial density field in kg/m^3
  //portableWrapper::assign(data.rho, 1.0e-6);

}
