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

////////////////////////////////////////////////////////////////////////////////////////
void simulation::two_fluid_grid(simulationData &data,simulationData &dataNeutral){
     /*
     This is my work around for allocating the grid properties to the neutral data structure
     */
     dataNeutral.cv=data.cv;
     dataNeutral.dxab=data.dxab;
     dataNeutral.dyab=data.dyab;
     dataNeutral.dzab=data.dzab;
     dataNeutral.dxac=data.dxac;
     dataNeutral.dyac=data.dyac;
     dataNeutral.dzac=data.dzac;
     dataNeutral.cv1=data.cv1;
     dataNeutral.cvc=data.cvc;
     dataNeutral.xc=data.xc;
     dataNeutral.yc=data.yc;
     dataNeutral.zc=data.zc;
     dataNeutral.xb=data.xb;
     dataNeutral.yb=data.yb;
     dataNeutral.zb=data.zb;
     dataNeutral.xb_global=data.xb_global;
     dataNeutral.yb_global=data.yb_global;
     dataNeutral.zb_global=data.zb_global;
     dataNeutral.dxc=data.dxc;
     dataNeutral.dyc=data.dyc;
     dataNeutral.dzc=data.dzc;
     dataNeutral.dxb=data.dxb;
     dataNeutral.dyb=data.dyb;
     dataNeutral.dzb=data.dzb;
     dataNeutral.hy=data.hy;
     dataNeutral.hz=data.hz;
     dataNeutral.hyc=data.hyc;
     dataNeutral.hzc=data.hzc;
     dataNeutral.hz1=data.hz1;
     dataNeutral.hz2=data.hz2;
     dataNeutral.x=data.x;
     dataNeutral.y=data.y;
     dataNeutral.z=data.z;
     dataNeutral.xp=data.xp;
     dataNeutral.yp=data.yp;
     dataNeutral.zp=data.zp;
}

////////////////////////////////////////////////////////////////////////////////////////
void simulation::two_fluid_source(simulationData &data,simulationData &dataNeutral){

    //data.two_fluid_timestep=1.0;
    
    //Calculate the source terms for the two-fluid interactions
    using Range = portableWrapper::Range;
    portableWrapper::applyKernel(LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) {
        //Get Temperatures
        T_dataType temperature_ion = data.gas_gamma*data.energy_ion(ix,iy,iz)*(data.gas_gamma-1.0);
        T_dataType temperature_electron = data.gas_gamma*data.energy_electron(ix,iy,iz)*(data.gas_gamma-1.0);
        T_dataType temperature_neutral = data.gas_gamma*dataNeutral.energy_neutral(ix,iy,iz)*(data.gas_gamma-1.0);

        //This needs temeprature dependence
        T_dataType ac=data.alpha0*std::sqrt(0.5*(temperature_neutral+temperature_ion));
        
        //Apply the velocity exchange terms
        data.vx(ix,iy,iz)       +=data.dt*ac*(dataNeutral.rho(ix,iy,iz)*dataNeutral.vx(ix,iy,iz)-dataNeutral.rho(ix,iy,iz)*data.vx(ix,iy,iz));
        dataNeutral.vx(ix,iy,iz)-=data.dt*ac*(data.rho(ix,iy,iz)       *dataNeutral.vx(ix,iy,iz)-data.rho(ix,iy,iz)       *data.vx(ix,iy,iz));
        
        data.vy(ix,iy,iz)       +=data.dt*ac*(dataNeutral.rho(ix,iy,iz)*dataNeutral.vy(ix,iy,iz)-dataNeutral.rho(ix,iy,iz)*data.vy(ix,iy,iz));
        dataNeutral.vy(ix,iy,iz)-=data.dt*ac*(data.rho(ix,iy,iz)       *dataNeutral.vy(ix,iy,iz)-data.rho(ix,iy,iz)       *data.vy(ix,iy,iz));
        
        data.vz(ix,iy,iz)       +=data.dt*ac*(dataNeutral.rho(ix,iy,iz)*dataNeutral.vz(ix,iy,iz)-dataNeutral.rho(ix,iy,iz)*data.vz(ix,iy,iz));
        dataNeutral.vz(ix,iy,iz)-=data.dt*ac*(data.rho(ix,iy,iz)       *dataNeutral.vz(ix,iy,iz)-data.rho(ix,iy,iz)       *data.vz(ix,iy,iz));
        
        //Energy source terms - the 3/2 here needs fixing
        data.energy_ion(ix,iy,iz)=data.energy_ion(ix,iy,iz)+data.dt*ac*dataNeutral.rho(ix,iy,iz)*(0.5*(\
                        (dataNeutral.vx(ix,iy,iz)*dataNeutral.vx(ix,iy,iz)-data.vx(ix,iy,iz)*data.vx(ix,iy,iz))+\
                        (dataNeutral.vy(ix,iy,iz)*dataNeutral.vy(ix,iy,iz)-data.vy(ix,iy,iz)*data.vy(ix,iy,iz))+\
                        (dataNeutral.vz(ix,iy,iz)*dataNeutral.vz(ix,iy,iz)-data.vz(ix,iy,iz)*data.vz(ix,iy,iz)))\
                        + 3.0/data.gas_gamma/2.0*(temperature_neutral-temperature_ion));
        data.energy_electron(ix,iy,iz)=data.energy_electron(ix,iy,iz)+data.dt*ac*dataNeutral.rho(ix,iy,iz)*(0.5*(\
                        (dataNeutral.vx(ix,iy,iz)*dataNeutral.vx(ix,iy,iz)-data.vx(ix,iy,iz)*data.vx(ix,iy,iz))+\
                        (dataNeutral.vy(ix,iy,iz)*dataNeutral.vy(ix,iy,iz)-data.vy(ix,iy,iz)*data.vy(ix,iy,iz))+\
                        (dataNeutral.vz(ix,iy,iz)*dataNeutral.vz(ix,iy,iz)-data.vz(ix,iy,iz)*data.vz(ix,iy,iz)))\
                        + 3.0/data.gas_gamma/2.0*(temperature_neutral-temperature_ion));
        dataNeutral.energy_neutral(ix,iy,iz)=dataNeutral.energy_neutral(ix,iy,iz)-data.dt*ac*dataNeutral.rho(ix,iy,iz)*(0.5*(\
                        (dataNeutral.vx(ix,iy,iz)*dataNeutral.vx(ix,iy,iz)-data.vx(ix,iy,iz)*data.vx(ix,iy,iz))+\
                        (dataNeutral.vy(ix,iy,iz)*dataNeutral.vy(ix,iy,iz)-data.vy(ix,iy,iz)*data.vy(ix,iy,iz))+\
                        (dataNeutral.vz(ix,iy,iz)*dataNeutral.vz(ix,iy,iz)-data.vz(ix,iy,iz)*data.vz(ix,iy,iz)))\
                        + 3.0/data.gas_gamma/2.0*(temperature_neutral-temperature_ion));  
                        
        //Two-fluid time-step
        //T_dataType collisional_timestep_temp=0.3/(ac*data.rho(ix,iy,iz));
        //if (data.two_fluid_timestep < collisional_timestep_temp) printf("%f \n", collisional_timestep_temp); 
        //data.two_fluid_timestep=collisional_timestep_temp;
        //collisional_timestep_temp=0.3/(ac*dataNeutral.rho(ix,iy,iz));
        //if (data.two_fluid_timestep < collisional_timestep_temp) data.two_fluid_timestep=collisional_timestep_temp;
              
    }, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));
    
    
    //Two-fluid time-step
    //printf("%f \n",data.two_fluid_timestep);
}

////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////
