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

void set_dt_collisional(simulationData &data, simulationData &dataNeutral);
void ion_rec_rates_empirical(simulationData &data, simulationData &dataNeutral);
void set_dt_ion_rec(simulationData &data,simulationData &dataNeutral);

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
                        
    }, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));
    
    //Get the ionisation rates
    if (data.ion_rec_empirical) ion_rec_rates_empirical(data,dataNeutral);
    
    //Set the timestep for the collisions
    set_dt_collisional(data, dataNeutral);
    
    //Set the ionisation/recombination timestep
    set_dt_ion_rec(data,dataNeutral);
    
    //Two-fluid time-step
    printf("%f \n",data.two_fluid_timestep);
}

////////////////////////////////////////////////////////////////////////////////////////

void ion_rec_rates_empirical(simulationData &data, simulationData &dataNeutral){

    //Much of this should go elsewhere
    T_dataType T0=1.0e4; //Reference temperature
    T_dataType n0=1.0e14; //Reference electron number density
    T_dataType t_ir=1.0e-5; //Reference recombination timescale (relative to collisional timescale)
	//Formulation from Snow+2021 paper
	//Empirical estimates for the rates

	T_dataType Te_0=T0/1.1604e4; //Calculate electron temperature in eV
	T_dataType rec_fac=2.6e-19*(n0*1.0e6)/std::sqrt(Te_0);  //reference recombination rate (n0 converted to m^-3)

	//initial equilibrium fractions
	T_dataType ioneq=(2.6e-19/std::sqrt(Te_0))/(2.91e-14/(0.232+13.6/Te_0)*std::pow(13.6/Te_0,0.39)*std::exp(-13.6/Te_0));
	T_dataType f_n=ioneq/(ioneq+1.0);
	T_dataType f_p=1.0-f_n;
	T_dataType f_p_p=2.0*f_p/(f_n+2.0*f_p);
	
    T_dataType tfac=0.5*f_p_p/f_p; //Normalisation assumes sound speed normalisation
	

    using Range = portableWrapper::Range;
    portableWrapper::applyKernel(LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) {
        //Get Temperatures
        T_dataType temperature_electron = data.gas_gamma*data.energy_electron(ix,iy,iz)*(data.gas_gamma-1.0);
        T_dataType numberDensity_electron=data.rho(ix,iy,iz); // This isn't actually the numebr density. Neet to fix

        //Get ionisation and recomination rates
    	data.Gm_rec(ix,iy,iz)=numberDensity_electron/std::sqrt(temperature_electron)*t_ir/f_p*std::sqrt(tfac);
    	data.Gm_ion(ix,iy,iz)=2.91e-14*(n0*1.0e6)*numberDensity_electron*std::exp(-13.6/Te_0/temperature_electron*tfac)*std::pow(13.6/Te_0/temperature_electron*tfac,0.39);
    	data.Gm_ion(ix,iy,iz)=data.Gm_ion(ix,iy,iz)/(0.232+13.6/Te_0/temperature_electron*tfac)/rec_fac/f_p *t_ir;        
    }, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));


}

////////////////////////////////////////////////////////////////////////////////////////

void set_dt_collisional(simulationData &data,simulationData &dataNeutral) {

    using Range = portableWrapper::Range;

    int i0 = data.geometry == geometryType::Cartesian ? 0:1;

    //Now need to do a map and reduction
    data.two_fluid_timestep = data.dt_multiplier * 
    portableWrapper::applyReduction(LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) {
        //Get Temperatures
        T_dataType temperature_ion = data.gas_gamma*data.energy_ion(ix,iy,iz)*(data.gas_gamma-1.0);
        T_dataType temperature_electron = data.gas_gamma*data.energy_electron(ix,iy,iz)*(data.gas_gamma-1.0);
        T_dataType temperature_neutral = data.gas_gamma*dataNeutral.energy_neutral(ix,iy,iz)*(data.gas_gamma-1.0);

        //This needs temeprature dependence
        T_dataType ac=data.alpha0*std::sqrt(0.5*(temperature_neutral+temperature_ion));
        
        T_dataType collisional_timestep_plasma=0.3/(ac*data.rho(ix,iy,iz));
        T_dataType collisional_timestep_neutral=0.3/(ac*dataNeutral.rho(ix,iy,iz));
        
        T_dataType t1;
        
        if (collisional_timestep_plasma < collisional_timestep_neutral){
            t1=collisional_timestep_plasma;
        } else{
            t1=collisional_timestep_neutral;
        }
                return t1;
    }, LAMBDA(T_dataType &a, const T_dataType &b) {
        a=portableWrapper::min(a, b);
    }, data.largest_number,
    Range(i0, data.nx), Range(0, data.ny), Range(0, data.nz));

}

////////////////////////////////////////////////////////////////////////////////////////

void set_dt_ion_rec(simulationData &data,simulationData &dataNeutral) {

    using Range = portableWrapper::Range;

    int i0 = data.geometry == geometryType::Cartesian ? 0:1;

    //Now need to do a map and reduction
    T_dataType ir_timestep= data.dt_multiplier * 
    portableWrapper::applyReduction(LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) {        
        T_dataType t1=1.0/(data.rho(ix,iy,iz)*data.Gm_rec(ix,iy,iz)-dataNeutral.rho(ix,iy,iz)*data.Gm_ion(ix,iy,iz));
        return t1;
    }, LAMBDA(T_dataType &a, const T_dataType &b) {
        a=portableWrapper::min(a, b);
    }, data.largest_number,
    Range(i0, data.nx), Range(0, data.ny), Range(0, data.nz));
    
    if (ir_timestep < data.two_fluid_timestep) data.two_fluid_timestep=ir_timestep;

}
