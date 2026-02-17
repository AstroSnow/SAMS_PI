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

  data.nx=32; // Number of cells in the x-direction
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
  data.x_min = -0.5;
  data.x_max = 0.5;
  data.y_min = -1.0;
  data.y_max = 1.0;
  data.z_min = -1.0;
  data.z_max = 1.0;

  // Boundary conditions
  data.xbc_min = BCType::BC_OTHER;
  data.xbc_max = BCType::BC_OTHER;
  data.ybc_min = BCType::BC_PERIODIC;
  data.ybc_max = BCType::BC_PERIODIC;
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
  data.gas_gamma = 1.66;

  // Average mass of an ion in proton masses
  data.mf = 1.66;

  // Resistive MHD options
  data.resistiveMHD = false;
  data.eta_background = 1.e-10;
  data.j_max = 1.0;
  data.eta0 = 2.e-10;

  // Remap kinetic energy correction
  data.rke = true;
  
  // Two-fluid flag
  data.two_fluid=true;
  data.collisions=true;
  data.ion_rec=true;
  data.ion_rec_empirical=false;
  data.ion_rec_nlevel=true;
  data.alpha0=100.0;
  data.T_reference=10000.0; //Reference electron temperature (Kelvin)
  data.ne_reference=1.0e14; //Reference electron number density (cm^-3)

  // Output frequency and directory
  data.dt_snapshots = 0.02;
}

void simulation::initial_conditions(simulationData &data,simulationData &dataNeutral) {

  using Range = portableWrapper::Range;

  SAMS::cout << "Setting up initial conditions" << std::endl;

  char shock_tube_problem[8]="sod";
  //char shock_tube_problem[8]="briowu";
  
  portableWrapper::assign(data.vx,0.0);
  portableWrapper::assign(data.vy,0.0);
  portableWrapper::assign(data.vz,0.0);
  
  portableWrapper::assign(data.bx,0.0);
  portableWrapper::assign(data.by,0.0);
  portableWrapper::assign(data.bz,0.0);
 
  
  //declare some arrays. Rewritten from the shock tube declaration
  T_dataType rho_L = 1.0;
  T_dataType P_L = 1.0;
  T_dataType vx_L = 0.0;
  T_dataType vy_L = 0.0;
  T_dataType vz_L = 0.0;
  T_dataType bx_L = 0.0;
  T_dataType by_L = 0.0;
  T_dataType bz_L = 0.0; 
  T_dataType rho_R = 0.125;
  T_dataType P_R = 0.1;
  T_dataType vx_R = 0.0;
  T_dataType vy_R = 0.0;
  T_dataType vz_R = 0.0;
  T_dataType bx_R = 0.0;
  T_dataType by_R = 0.0;
  T_dataType bz_R = 0.0;
  ////////////////////////////////////////////////////
  if (std::strcmp(shock_tube_problem,"sod")==0){
      // Sod Shock tube
      printf("Sod Shock Tube \n");
      rho_L = 1.0;
      P_L = 1.0;
      vx_L = 0.0;
      
      rho_R = 0.125;
      P_R = 0.1;
      vx_R = 0.0;
      
      if (data.ion_rec_empirical){
        P_L=P_L/data.gas_gamma; //Normalise to sound speed of 1
        P_R=P_R/data.gas_gamma;
      }
  }
  if (std::strcmp(shock_tube_problem,"briowu")==0){
      // Brio & Wu Shock tube
      printf("Brio Wu Shock Tube \n");
      rho_L = 1.0;
      P_L = 1.0;
      vx_L = 0.0;
      bx_L = 0.75;
      by_L = 1.0;
      
      rho_R = 0.125;
      P_R = 0.1;
      vx_R = 0.0;
      bx_R = 0.75;
      by_R = -1.0;
  }
  
  T_dataType en_L =P_L/2.0/rho_L/(data.gas_gamma-1.0); //The half comes from electron and proton pressures
  T_dataType en_R =P_R/2.0/rho_R/(data.gas_gamma-1.0);
  
  //Assign perminent values for BCs
  data.rho_L=rho_L;
  data.vx_L=vx_L;
  data.vy_L=vy_L;
  data.vz_L=vz_L;
  data.bx_L=bx_L;
  data.by_L=by_L;
  data.bz_L=bz_L;
  data.en_L=en_L;
  data.rho_R=rho_R;
  data.vx_R=vx_R;
  data.vy_R=vy_R;
  data.vz_R=vz_R;
  data.bx_R=bx_R;
  data.by_R=by_R;
  data.bz_R=bz_R;
  data.en_R=en_R;
  
  if (data.two_fluid){
      dataNeutral.rho_L=rho_L;
      dataNeutral.vx_L=vx_L;
      dataNeutral.vy_L=vy_L;
      dataNeutral.vz_L=vz_L;
      dataNeutral.bx_L=0.0;
      dataNeutral.by_L=0.0;
      dataNeutral.bz_L=0.0;
      dataNeutral.en_L=en_L*2.0; //the factor of 2 is to pair up with total plasma pressure
      dataNeutral.rho_R=rho_R;
      dataNeutral.vx_R=vx_R;
      dataNeutral.vy_R=vy_R;
      dataNeutral.vz_R=vz_R;
      dataNeutral.bx_R=0.0;
      dataNeutral.by_R=0.0;
      dataNeutral.bz_R=0.0;
      dataNeutral.en_R=en_R*2.0;
      
      if (data.ion_rec_empirical){

        //Much of this should go elsewhere
        T_dataType T0=data.T_reference; //Reference temperature
        T_dataType n0=data.ne_reference; //Reference electron number density
        T_dataType t_ir=1.0e0; //Reference recombination timescale (relative to collisional timescale)

	    T_dataType Te_0=T0/1.1604e4; //Calculate electron temperature in eV
	    T_dataType rec_fac=2.6e-19*(n0*1.0e6)/std::sqrt(Te_0);  //reference recombination rate (n0 converted to m^-3)

	    //initial equilibrium fractions
	    T_dataType ioneq=(2.6e-19/std::sqrt(Te_0))/(2.91e-14/(0.232+13.6/Te_0)*std::pow(13.6/Te_0,0.39)*std::exp(-13.6/Te_0));
	    T_dataType f_n=ioneq/(ioneq+1.0);
	    T_dataType f_p=1.0-f_n;
	    T_dataType f_p_p=f_p/(f_n+2.0*f_p); //note the missing factor of 2 from the PIP code since electron and proton have different energies
	    T_dataType f_p_n=f_n/(f_n+2.0*f_p); 
        data.rho_L=data.rho_L*f_p;
        data.rho_R=data.rho_R*f_p; //This is not in balance!
        dataNeutral.rho_L=dataNeutral.rho_L*f_n;
        dataNeutral.rho_R=dataNeutral.rho_R*f_n; //This is not in balance!
        
        data.en_L=data.en_L*f_p_p/f_p;
        data.en_R=data.en_R*f_p_p/f_p; //This is not in balance!
        dataNeutral.en_L=dataNeutral.en_L*f_p_n/f_n;
        dataNeutral.en_R=dataNeutral.en_R*f_p_n/f_n; //This is not in balance!
      
      }
  }
  
  T_dataType w_lay=0.01;
  ////////////////////////////////////////////////////

  portableWrapper::assign(data.rho,rho_R);
  portableWrapper::assign(data.vx,vx_R);
  portableWrapper::assign(data.vy,vy_R);
  portableWrapper::assign(data.vz,vz_R);
  portableWrapper::assign(data.bx,bx_R);
  portableWrapper::assign(data.by,by_R);
  portableWrapper::assign(data.bz,bz_R);
  portableWrapper::assign(data.energy_ion,P_R/2.0/rho_R/(data.gas_gamma-1.0));
  portableWrapper::assign(data.energy_electron,P_R/2.0/rho_R/(data.gas_gamma-1.0));

  //Some Neutral conditions
  /*if (data.two_fluid) {
    portableWrapper::assign(dataNeutral.vx,vx_R);
    portableWrapper::assign(dataNeutral.vy,vy_R);
    portableWrapper::assign(dataNeutral.vz,vz_R);
    portableWrapper::assign(dataNeutral.bx,0.0);
    portableWrapper::assign(dataNeutral.by,0.0);
    portableWrapper::assign(dataNeutral.bz,0.0);
    portableWrapper::assign(dataNeutral.rho,rho_R);
    portableWrapper::assign(dataNeutral.energy_neutral,P_R/rho_R/(data.gas_gamma-1.0));    
    }
  */

  if ((std::strcmp(shock_tube_problem,"briowu")==0) || (std::strcmp(shock_tube_problem,"sod")==0)){
      portableWrapper::applyKernel(
        LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) {
        
        data.vx(ix, iy, iz)=data.vx_L+(data.vx_R-data.vx_L)*(std::tanh(data.xb(ix)/w_lay)+1.0)*0.5;
        data.vy(ix, iy, iz)=data.vy_L+(data.vy_R-data.vy_L)*(std::tanh(data.xb(ix)/w_lay)+1.0)*0.5;
        data.vz(ix, iy, iz)=data.vz_L+(data.vz_R-data.vz_L)*(std::tanh(data.xb(ix)/w_lay)+1.0)*0.5;
        data.bx(ix, iy, iz)=data.bx_L+(data.bx_R-data.bx_L)*(std::tanh(data.xb(ix)/w_lay)+1.0)*0.5;
        data.by(ix, iy, iz)=data.by_L+(data.by_R-data.by_L)*(std::tanh(data.xb(ix)/w_lay)+1.0)*0.5;
        data.bz(ix, iy, iz)=data.bz_L+(data.bz_R-data.bz_L)*(std::tanh(data.xb(ix)/w_lay)+1.0)*0.5;
        data.rho(ix, iy, iz)=data.rho_L+(data.rho_R-data.rho_L)*(std::tanh(data.xb(ix)/w_lay)+1.0)*0.5;
        data.energy_ion(ix, iy, iz)=data.en_L+(data.en_R-data.en_L)*(std::tanh(data.xb(ix)/w_lay)+1.0)*0.5;
        data.energy_electron(ix, iy, iz)=data.en_L+(data.en_R-data.en_L)*(std::tanh(data.xb(ix)/w_lay)+1.0)*0.5;
          if (data.two_fluid) {
              portableWrapper::assign(dataNeutral.bx,0.0);
              portableWrapper::assign(dataNeutral.by,0.0);
              portableWrapper::assign(dataNeutral.bz,0.0);
              dataNeutral.vx(ix, iy, iz) = dataNeutral.vx_L+(dataNeutral.vx_R-dataNeutral.vx_L)*(std::tanh(data.xb(ix)/w_lay)+1.0)*0.5;
              dataNeutral.vy(ix, iy, iz) = dataNeutral.vy_L+(dataNeutral.vy_R-dataNeutral.vy_L)*(std::tanh(data.xb(ix)/w_lay)+1.0)*0.5;
              dataNeutral.vz(ix, iy, iz) = dataNeutral.vz_L+(dataNeutral.vz_R-dataNeutral.vz_L)*(std::tanh(data.xb(ix)/w_lay)+1.0)*0.5;
              dataNeutral.rho(ix, iy, iz) = dataNeutral.rho_L+(dataNeutral.rho_R-dataNeutral.rho_L)*(std::tanh(data.xb(ix)/w_lay)+1.0)*0.5;
              dataNeutral.energy_neutral(ix, iy, iz) = dataNeutral.en_L+(dataNeutral.en_R-dataNeutral.en_L)*(std::tanh(data.xb(ix)/w_lay)+1.0)*0.5;
          }
          //printf("%ld %f %f \n",ix,data.energy_ion(ix,iy,iz),dataNeutral.energy_neutral(ix,iy,iz));
        }, 
        portableWrapper::Range(-1, data.nx+1),
        portableWrapper::Range(-1, data.ny+1),
        portableWrapper::Range(-1, data.nz+1)
      );
  }
  else{
  printf("unknown initial condition");
  }

  printf("%f %f \n",data.rho(1,1,1),data.energy_electron(1,1,1)*(data.gas_gamma-1.0)*data.rho(1,1,1));

  if (data.rke) portableWrapper::assign(data.delta_ke, 0.0);

}
