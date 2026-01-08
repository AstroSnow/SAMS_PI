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


struct data_two_fluid_source_ir
{
    volumeArray source_mass; // mass source term
    volumeArray source_v_x; // velocity source term
    volumeArray source_v_y; // velocity source term
    volumeArray source_v_z; // velocity source term
    volumeArray source_energy; // energy source term
    volumeArray source_electron_energy; // energy source term
    volumeArray ac; //coupling coeficient
};

void get_ac(simulationData &data, simulationData &dataNeutral, data_two_fluid_source_ir &plasma_ir_source);
void set_dt_collisional(simulationData &data, simulationData &dataNeutral, data_two_fluid_source_ir &plasma_ir_source);
void get_collisional_source_terms(simulationData &data, simulationData &dataNeutral, data_two_fluid_source_ir &plasma_ir_source, data_two_fluid_source_ir &neutral_ir_source);
void ion_rec_rates_empirical(simulationData &data, simulationData &dataNeutral);
void get_ion_rec_source_terms(simulationData &data, simulationData &dataNeutral, data_two_fluid_source_ir &plasma_ir_source, data_two_fluid_source_ir &neutral_ir_source);
void set_dt_ion_rec(simulationData &data,simulationData &dataNeutral);
//void get_ac(T_dataType alpha0,T_dataType temperature_ion,T_dataType temperature_neutral);

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
     dataNeutral.isxLB=data.isxLB;
     dataNeutral.isxUB=data.isxUB;
     dataNeutral.isyLB=data.isyLB;
     dataNeutral.isyUB=data.isyUB;
     dataNeutral.iszLB=data.iszLB;
     dataNeutral.iszUB=data.iszUB;
     
}

////////////////////////////////////////////////////////////////////////////////////////
void simulation::two_fluid_source(simulationData &data,simulationData &dataNeutral,bool first_step){

    //data.two_fluid_timestep=1.0;
    
    data_two_fluid_source_ir plasma_ir_source;
    data_two_fluid_source_ir neutral_ir_source;
    portableWrapper::portableArrayManager irSourceManager;
    using Range = portableWrapper::Range;
    
    irSourceManager.allocate(plasma_ir_source.ac, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));
    irSourceManager.allocate(plasma_ir_source.source_mass, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));
    irSourceManager.allocate(plasma_ir_source.source_v_x, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));
    irSourceManager.allocate(plasma_ir_source.source_v_y, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));
    irSourceManager.allocate(plasma_ir_source.source_v_z, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));
    irSourceManager.allocate(plasma_ir_source.source_energy, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));
    irSourceManager.allocate(plasma_ir_source.source_electron_energy, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));
    irSourceManager.allocate(neutral_ir_source.source_mass, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));
    irSourceManager.allocate(neutral_ir_source.source_v_x, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));
    irSourceManager.allocate(neutral_ir_source.source_v_y, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));
    irSourceManager.allocate(neutral_ir_source.source_v_z, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));
    irSourceManager.allocate(neutral_ir_source.source_energy, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));
    
    portableWrapper::assign(plasma_ir_source.source_mass,0.0);
    portableWrapper::assign(plasma_ir_source.source_v_x,0.0);
    portableWrapper::assign(plasma_ir_source.source_v_y,0.0);
    portableWrapper::assign(plasma_ir_source.source_v_z,0.0);
    portableWrapper::assign(plasma_ir_source.source_energy,0.0);
    portableWrapper::assign(plasma_ir_source.source_electron_energy,0.0);
    portableWrapper::assign(neutral_ir_source.source_mass,0.0);
    portableWrapper::assign(neutral_ir_source.source_v_x,0.0);
    portableWrapper::assign(neutral_ir_source.source_v_y,0.0);
    portableWrapper::assign(neutral_ir_source.source_v_z,0.0);
    portableWrapper::assign(neutral_ir_source.source_energy,0.0);
    portableWrapper::assign(plasma_ir_source.ac,0.0);
    
    //Get collisional coefficient
    get_ac(data,dataNeutral,plasma_ir_source);
    
    //Get the ionisation rates
    if (data.ion_rec_empirical){        
        ion_rec_rates_empirical(data,dataNeutral);
    }
    
    //Calculate the source terms for the two-fluid interactions
    get_collisional_source_terms(data,dataNeutral,plasma_ir_source,neutral_ir_source);
    
    //Calculate the source terms for Ionisation/recombination
    if (data.ion_rec) get_ion_rec_source_terms(data,dataNeutral,plasma_ir_source,neutral_ir_source);
    
    
    // Make sure the timestep is the same in both species
    //Set dt to be the minimum of the neutral and plasma times
    if (first_step){        
        //Set the timestep for the collisions
        set_dt_collisional(data, dataNeutral, plasma_ir_source);
        
        //Set the ionisation/recombination timestep
        if (data.ion_rec_empirical) set_dt_ion_rec(data,dataNeutral);
        
        printf("dt (plasma, neutral, two-fluid)=%f %f %f \n",data.dt,dataNeutral.dt,data.two_fluid_timestep);
        
        data.dt=std::min({dataNeutral.dt,data.dt,data.two_fluid_timestep});
        dataNeutral.dt=data.dt;
    }
    
    if (data.collisions){
        using Range = portableWrapper::Range;
        portableWrapper::applyKernel(LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) {
            //Get Temperatures
            //T_dataType temperature_ion = data.gas_gamma*data.energy_ion(ix,iy,iz)*(data.gas_gamma-1.0);
            //T_dataType temperature_electron = data.gas_gamma*data.energy_electron(ix,iy,iz)*(data.gas_gamma-1.0);
            //T_dataType temperature_neutral = data.gas_gamma*dataNeutral.energy_neutral(ix,iy,iz)*(data.gas_gamma-1.0);
            
            //T_dataType ac;
            //get_ac(data.alpha0,temperature_ion,temperature_neutral);
            
            //Note that the factor of 0.5 in these is due to Strang splitting
            //Mass exchange terms
            data.rho(ix,iy,iz)+=0.5*data.dt*plasma_ir_source.source_mass(ix,iy,iz);
            dataNeutral.rho(ix,iy,iz)+=0.5*data.dt*neutral_ir_source.source_mass(ix,iy,iz);
            
            //Apply the velocity exchange terms
            data.vx(ix,iy,iz)       +=0.5*data.dt*plasma_ir_source.source_v_x(ix,iy,iz);
            dataNeutral.vx(ix,iy,iz)-=0.5*data.dt*neutral_ir_source.source_v_x(ix,iy,iz);
            
            data.vy(ix,iy,iz)       +=0.5*data.dt*plasma_ir_source.source_v_y(ix,iy,iz);
            dataNeutral.vy(ix,iy,iz)-=0.5*data.dt*neutral_ir_source.source_v_y(ix,iy,iz);
            
            data.vz(ix,iy,iz)       +=0.5*data.dt*plasma_ir_source.source_v_z(ix,iy,iz);
            dataNeutral.vz(ix,iy,iz)-=0.5*data.dt*neutral_ir_source.source_v_z(ix,iy,iz);
            
            //Energy source terms - the 3/2 here needs fixing
            data.energy_ion(ix,iy,iz)+=0.5*data.dt*plasma_ir_source.source_energy(ix,iy,iz);
            //data.energy_electron(ix,iy,iz)+=0.5*data.dt*plasma_ir_source.source_electron_energy(ix,iy,iz);
            dataNeutral.energy_neutral(ix,iy,iz)+=0.5*data.dt*neutral_ir_source.source_energy(ix,iy,iz);                 
        }, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));
    }
    
}

////////////////////////////////////////////////////////////////////////////////////////
//Get the collisional coupling coefficient
void get_ac(simulationData &data, simulationData &dataNeutral, data_two_fluid_source_ir &plasma_ir_source){

    using Range = portableWrapper::Range;
    portableWrapper::applyKernel(LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) {
        //Get Temperatures
        T_dataType temperature_electron = data.gas_gamma*data.energy_electron(ix,iy,iz)*(data.gas_gamma-1.0);
        T_dataType temperature_ion = data.gas_gamma*data.energy_ion(ix,iy,iz)*(data.gas_gamma-1.0);
        T_dataType temperature_neutral = data.gas_gamma*dataNeutral.energy_neutral(ix,iy,iz)*(data.gas_gamma-1.0);
        
        plasma_ir_source.ac(ix,iy,iz)=data.alpha0*std::sqrt(0.5*(temperature_neutral+temperature_ion));
        
//        printf("ix,iy,iz, t_i t_n ac :  %li %li %li %f %f %f \n",ix,iy,iz,dataNeutral.rho(ix,iy,iz),dataNeutral.energy_neutral(ix,iy,iz),plasma_ir_source.ac(ix,iy,iz));
        //printf("ix,iy,iz, t_i t_n ac :  %li %li %li %f %f %f \n",ix,iy,iz,temperature_ion,temperature_neutral,plasma_ir_source.ac(ix,iy,iz));
        
    	}, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));
    //return alpha0*std::sqrt(0.5*(temperature_neutral+temperature_ion));

}


////////////////////////////////////////////////////////////////////////////////////////
//Formulation from Snow+2021 paper
//Empirical estimates for the rates
//Controlled using the data.ion_rec_empirical in control.cpp
void ion_rec_rates_empirical(simulationData &data, simulationData &dataNeutral){

    //Much of this should go elsewhere
    T_dataType T0=1.0e4; //Reference temperature
    T_dataType n0=1.0e14; //Reference electron number density
    T_dataType t_ir=1.0e-5; //Reference recombination timescale (relative to collisional timescale)

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
//Formulation from Snow+2023 paper using Jeffries1968
//Controlled using the data.ion_rec_jeffries in control.cpp
//Not used yet
void ion_rec_rates_jeffries(simulationData &data, simulationData &dataNeutral){

    //Much of this should go elsewhere
    T_dataType T0=1.0e4; //Reference temperature
    T_dataType n0=1.0e14; //Reference electron number density
    T_dataType t_ir=1.0e-5; //Reference recombination timescale (relative to collisional timescale)

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
    	//data.Gm_rec(ix,iy,iz)=
    	//data.Gm_ion(ix,iy,iz)=       
    }, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));


}

////////////////////////////////////////////////////////////////////////////////////////
//Get the source terms for the IR rates
void get_ion_rec_source_terms(simulationData &data, simulationData &dataNeutral, data_two_fluid_source_ir &plasma_ir_source, data_two_fluid_source_ir &neutral_ir_source){	

    using Range = portableWrapper::Range;
    portableWrapper::applyKernel(LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) {
    
        //Mass source terms
        plasma_ir_source.source_mass(ix,iy,iz)  += data.Gm_ion(ix,iy,iz)*dataNeutral.rho(ix,iy,iz)-data.Gm_rec(ix,iy,iz)*data.rho(ix,iy,iz);
        neutral_ir_source.source_mass(ix,iy,iz) +=-data.Gm_ion(ix,iy,iz)*dataNeutral.rho(ix,iy,iz)+data.Gm_rec(ix,iy,iz)*data.rho(ix,iy,iz);
        
        //Velocity source terms
        T_dataType v_D_x  =  data.vx(ix,iy,iz) - dataNeutral.vx(ix,iy,iz); //Drift velocity in the x-direction
        T_dataType v_D_y  =  data.vy(ix,iy,iz) - dataNeutral.vy(ix,iy,iz); //Drift velocity in the y-direction
        T_dataType v_D_z  =  data.vz(ix,iy,iz) - dataNeutral.vz(ix,iy,iz); //Drift velocity in the z-direction
        plasma_ir_source.source_v_x(ix,iy,iz) += -data.Gm_ion(ix,iy,iz)*dataNeutral.rho(ix,iy,iz)*v_D_x/data.rho(ix,iy,iz);
        plasma_ir_source.source_v_y(ix,iy,iz) += -data.Gm_ion(ix,iy,iz)*dataNeutral.rho(ix,iy,iz)*v_D_y/data.rho(ix,iy,iz);
        plasma_ir_source.source_v_z(ix,iy,iz) += -data.Gm_ion(ix,iy,iz)*dataNeutral.rho(ix,iy,iz)*v_D_z/data.rho(ix,iy,iz);
        neutral_ir_source.source_v_x(ix,iy,iz) += data.Gm_rec(ix,iy,iz)*data.rho(ix,iy,iz)*v_D_x/dataNeutral.rho(ix,iy,iz);
        neutral_ir_source.source_v_y(ix,iy,iz) += data.Gm_rec(ix,iy,iz)*data.rho(ix,iy,iz)*v_D_y/dataNeutral.rho(ix,iy,iz);
        neutral_ir_source.source_v_z(ix,iy,iz) += data.Gm_rec(ix,iy,iz)*data.rho(ix,iy,iz)*v_D_z/dataNeutral.rho(ix,iy,iz);
        
        //Energy source terms
        plasma_ir_source.source_energy(ix,iy,iz) += -0.5*(data.Gm_rec(ix,iy,iz)*(data.vx(ix,iy,iz)*data.vx(ix,iy,iz)+
                                                                                data.vy(ix,iy,iz)*data.vy(ix,iy,iz)+
                                                                                data.vz(ix,iy,iz)*data.vz(ix,iy,iz))
                                                        -data.Gm_ion(ix,iy,iz)*(dataNeutral.vx(ix,iy,iz)*data.vx(ix,iy,iz)+
                                                                                dataNeutral.vy(ix,iy,iz)*data.vy(ix,iy,iz)+
                                                                                dataNeutral.vz(ix,iy,iz)*data.vz(ix,iy,iz))
                                                                              *dataNeutral.rho(ix,iy,iz)/data.rho(ix,iy,iz)                                                                                
                                                        )
                                                   -(data.Gm_rec(ix,iy,iz)*data.energy_ion(ix,iy,iz)-data.Gm_ion(ix,iy,iz)*dataNeutral.energy_neutral(ix,iy,iz)*dataNeutral.rho(ix,iy,iz)/data.rho(ix,iy,iz)); //Is this electron or ion energy (or mean energy)? is the half needed?
        
        neutral_ir_source.source_energy(ix,iy,iz) += 0.5*(data.Gm_rec(ix,iy,iz)*(data.vx(ix,iy,iz)*data.vx(ix,iy,iz)+
                                                                                data.vy(ix,iy,iz)*data.vy(ix,iy,iz)+
                                                                                data.vz(ix,iy,iz)*data.vz(ix,iy,iz))
                                                                                *data.rho(ix,iy,iz)/dataNeutral.rho(ix,iy,iz)
                                                        -data.Gm_ion(ix,iy,iz)*(dataNeutral.vx(ix,iy,iz)*data.vx(ix,iy,iz)+
                                                                                dataNeutral.vy(ix,iy,iz)*data.vy(ix,iy,iz)+
                                                                                dataNeutral.vz(ix,iy,iz)*data.vz(ix,iy,iz))
                                                        )
                                                   +(data.Gm_rec(ix,iy,iz)*data.energy_ion(ix,iy,iz)*data.rho(ix,iy,iz)/dataNeutral.rho(ix,iy,iz)-data.Gm_ion(ix,iy,iz)*dataNeutral.energy_neutral(ix,iy,iz)); //Is this electron or ion energy (or mean energy)? is the half needed?
        
        
    }, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));


}

////////////////////////////////////////////////////////////////////////////////////////
//Get the source terms for the IR rates
void get_collisional_source_terms(simulationData &data, simulationData &dataNeutral, data_two_fluid_source_ir &plasma_ir_source, data_two_fluid_source_ir &neutral_ir_source){	

    using Range = portableWrapper::Range;
    portableWrapper::applyKernel(LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) {
            
        //Get Temperatures
        T_dataType temperature_ion = data.gas_gamma*data.energy_ion(ix,iy,iz)*(data.gas_gamma-1.0);
        T_dataType temperature_electron = data.gas_gamma*data.energy_electron(ix,iy,iz)*(data.gas_gamma-1.0);
        T_dataType temperature_neutral = data.gas_gamma*dataNeutral.energy_neutral(ix,iy,iz)*(data.gas_gamma-1.0);
        
        //T_dataType ac;
        //get_ac(data.alpha0,temperature_ion,temperature_neutral);
        
        
        //Get ac and rho at the location of v (vertex)
        T_dataType ac_vertex=(plasma_ir_source.ac(ix  , iy  , iz  ) + 
                                        plasma_ir_source.ac(ix+1, iy  , iz  ) + 
                                        plasma_ir_source.ac(ix  , iy+1, iz  ) + 
                                        plasma_ir_source.ac(ix+1, iy+1, iz  ) + 
                                        plasma_ir_source.ac(ix  , iy  , iz+1) + 
                                        plasma_ir_source.ac(ix+1, iy  , iz+1) + 
                                        plasma_ir_source.ac(ix  , iy+1, iz+1) + 
                                        plasma_ir_source.ac(ix+1, iy+1, iz+1))* 
                                        0.125;
        T_dataType rho_plasma_vertex=  (data.rho(ix  , iy  , iz  ) + 
                                        data.rho(ix+1, iy  , iz  ) + 
                                        data.rho(ix  , iy+1, iz  ) + 
                                        data.rho(ix+1, iy+1, iz  ) + 
                                        data.rho(ix  , iy  , iz+1) + 
                                        data.rho(ix+1, iy  , iz+1) + 
                                        data.rho(ix  , iy+1, iz+1) + 
                                        data.rho(ix+1, iy+1, iz+1))* 
                                        0.125;
        T_dataType rho_neutral_vertex=  (dataNeutral.rho(ix  , iy  , iz  ) + 
                                         dataNeutral.rho(ix+1, iy  , iz  ) + 
                                         dataNeutral.rho(ix  , iy+1, iz  ) + 
                                         dataNeutral.rho(ix+1, iy+1, iz  ) + 
                                         dataNeutral.rho(ix  , iy  , iz+1) + 
                                         dataNeutral.rho(ix+1, iy  , iz+1) + 
                                         dataNeutral.rho(ix  , iy+1, iz+1) + 
                                         dataNeutral.rho(ix+1, iy+1, iz+1))* 
                                         0.125;
        
                
        //Apply the velocity exchange terms
        plasma_ir_source.source_v_x(ix,iy,iz)+=ac_vertex*(dataNeutral.rho(ix,iy,iz)*dataNeutral.vx(ix,iy,iz)-dataNeutral.rho(ix,iy,iz)*data.vx(ix,iy,iz));
        neutral_ir_source.source_v_x(ix,iy,iz)-=ac_vertex*(data.rho(ix,iy,iz)*dataNeutral.vx(ix,iy,iz)-data.rho(ix,iy,iz)       *data.vx(ix,iy,iz));
        
        plasma_ir_source.source_v_y(ix,iy,iz)+=ac_vertex*(dataNeutral.rho(ix,iy,iz)*dataNeutral.vy(ix,iy,iz)-dataNeutral.rho(ix,iy,iz)*data.vy(ix,iy,iz));
        neutral_ir_source.source_v_y(ix,iy,iz)-=ac_vertex*(data.rho(ix,iy,iz)*dataNeutral.vy(ix,iy,iz)-data.rho(ix,iy,iz)       *data.vy(ix,iy,iz));
        
        plasma_ir_source.source_v_z(ix,iy,iz)+=ac_vertex*(dataNeutral.rho(ix,iy,iz)*dataNeutral.vz(ix,iy,iz)-dataNeutral.rho(ix,iy,iz)*data.vz(ix,iy,iz));
        neutral_ir_source.source_v_z(ix,iy,iz)-=ac_vertex*(data.rho(ix,iy,iz)*dataNeutral.vz(ix,iy,iz)-data.rho(ix,iy,iz)       *data.vz(ix,iy,iz));
        
        //Get velocity at cell centre
        
        //Energy source terms - the 3/2 here needs fixing
        plasma_ir_source.source_energy(ix,iy,iz)=plasma_ir_source.ac(ix,iy,iz)*dataNeutral.rho(ix,iy,iz)*(0.5*(\
                        (dataNeutral.vx(ix,iy,iz)*dataNeutral.vx(ix,iy,iz)-data.vx(ix,iy,iz)*data.vx(ix,iy,iz))+\
                        (dataNeutral.vy(ix,iy,iz)*dataNeutral.vy(ix,iy,iz)-data.vy(ix,iy,iz)*data.vy(ix,iy,iz))+\
                        (dataNeutral.vz(ix,iy,iz)*dataNeutral.vz(ix,iy,iz)-data.vz(ix,iy,iz)*data.vz(ix,iy,iz)))\
                        + 3.0/data.gas_gamma/2.0*(temperature_neutral-temperature_ion));
        neutral_ir_source.source_energy(ix,iy,iz)=-plasma_ir_source.ac(ix,iy,iz)*dataNeutral.rho(ix,iy,iz)*(0.5*(\
                        (dataNeutral.vx(ix,iy,iz)*dataNeutral.vx(ix,iy,iz)-data.vx(ix,iy,iz)*data.vx(ix,iy,iz))+\
                        (dataNeutral.vy(ix,iy,iz)*dataNeutral.vy(ix,iy,iz)-data.vy(ix,iy,iz)*data.vy(ix,iy,iz))+\
                        (dataNeutral.vz(ix,iy,iz)*dataNeutral.vz(ix,iy,iz)-data.vz(ix,iy,iz)*data.vz(ix,iy,iz)))\
                        + 3.0/data.gas_gamma/2.0*(temperature_neutral-temperature_ion));  
                                   
        
    }, Range(-1,data.nx+1), Range(-1,data.ny+1), Range(-1,data.nz+1));


}

////////////////////////////////////////////////////////////////////////////////////////
//Collisional timestep calculation
//Assuming normalisation to the sound speed
void set_dt_collisional(simulationData &data,simulationData &dataNeutral, data_two_fluid_source_ir &plasma_ir_source) {

    using Range = portableWrapper::Range;

    int i0 = data.geometry == geometryType::Cartesian ? 0:1;

    //Now need to do a map and reduction
    data.two_fluid_timestep = data.dt_multiplier * 
    portableWrapper::applyReduction(LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) {
        //Get Temperatures
        //T_dataType temperature_ion = data.gas_gamma*data.energy_ion(ix,iy,iz)*(data.gas_gamma-1.0);
        //T_dataType temperature_electron = data.gas_gamma*data.energy_electron(ix,iy,iz)*(data.gas_gamma-1.0);
        //T_dataType temperature_neutral = data.gas_gamma*dataNeutral.energy_neutral(ix,iy,iz)*(data.gas_gamma-1.0);

        //This needs temeprature dependence
        //T_dataType ac=data.alpha0*std::sqrt(0.5*(temperature_neutral+temperature_ion));
        
        T_dataType collisional_timestep_plasma=0.3/(plasma_ir_source.ac(ix,iy,iz)*data.rho(ix,iy,iz));
        T_dataType collisional_timestep_neutral=0.3/(plasma_ir_source.ac(ix,iy,iz)*dataNeutral.rho(ix,iy,iz));
        
        T_dataType t1 = std::min(collisional_timestep_plasma,collisional_timestep_neutral);
        
        //printf("ix,iy,iz, t1 :  %li %li %li %f %f \n",ix,iy,iz,collisional_timestep_plasma,plasma_ir_source.ac(ix,iy,iz));
        //if (collisional_timestep_plasma < collisional_timestep_neutral){
        //    t1=collisional_timestep_plasma;
        //} else{
        //    t1=collisional_timestep_neutral;
        //}
                return t1;
    }, LAMBDA(T_dataType &a, const T_dataType &b) {
        a=portableWrapper::min(a, b);
    }, data.largest_number,
    Range(i0, data.nx), Range(0, data.ny), Range(0, data.nz));

}

////////////////////////////////////////////////////////////////////////////////////////
//Timestep associated with ionisation/recombination
void set_dt_ion_rec(simulationData &data,simulationData &dataNeutral) {

    using Range = portableWrapper::Range;

    int i0 = data.geometry == geometryType::Cartesian ? 0:1;

    //Now need to do a map and reduction
    T_dataType ir_timestep= data.dt_multiplier * 
    portableWrapper::applyReduction(LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) {        
        T_dataType t1=std::abs(1.0/(data.rho(ix,iy,iz)*data.Gm_rec(ix,iy,iz)-dataNeutral.rho(ix,iy,iz)*data.Gm_ion(ix,iy,iz)));
        return t1;
    }, LAMBDA(T_dataType &a, const T_dataType &b) {
        a=portableWrapper::min(a, b);
    }, data.largest_number,
    Range(i0, data.nx), Range(0, data.ny), Range(0, data.nz));
    
    if (ir_timestep < data.two_fluid_timestep) data.two_fluid_timestep=ir_timestep;

}
