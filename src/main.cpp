#define PRINT_PARALLELIZATION_INFO

#include "pp/range.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include "shared_data.h"
#include "include/timer.h"
#include "axisRegistry.h"
#include "variableRegistry.h"
#include "mpiManager.h"
#include "welcome.h"

int main(int argc, char *argv[]){ 


    //Initialize MPI
    SAMS::MPI::initialize(argc, argv);

    //Get the MPI manager for the default communicator
    SAMS::MPIManager<SAMS::MPI_DECOMPOSITION_RANK>& mpi = SAMS::getMPIManager<SAMS::MPI_DECOMPOSITION_RANK>();

    SAMS::printWelcomeMessage();
    //MPI auto decomposition
    //mpi.autoDecomposition({false,false,false});    
    mpi.autoDecomposition({true,true,true});
    //Initialize portable wrapper
    portableWrapper::initialize(argc, argv);
    SAMS::finishWelcomeMessage();

    //Create the simulation (LARE) and data objects
    simulation S;
    //simulation S2;
    simulationData data;
    simulationData dataNeutral;

    //Setup control variables
    S.controlvariables(data);
    data.visc2_norm=data.visc2;
    data.is_neutral=false;
    
	if (data.two_fluid) {
	    printf("Initialising two-fluid arrays \n");
	    S.controlvariables(dataNeutral);
	    dataNeutral.is_neutral=true;
		printf("Finished initialising two-fluid arrays \n");
	}

    //Register axes and attach them to MPI dimensions
    auto& axRegistry = SAMS::getaxisRegistry();
    axRegistry.registerAxis("X", SAMS::MPIAxis(0));
    axRegistry.registerAxis("Y", SAMS::MPIAxis(1));
    axRegistry.registerAxis("Z", SAMS::MPIAxis(2));
    //if (data.ion_rec_nlevel) axRegistry.registerAxis("nlevels",SAMS::MPIAxis(3));
    //Tell LARE to register its variables
    S.registerVars(data.two_fluid,data.ion_rec);
    
    //Other simulations would register their variables here too

    //Set the axis domains and decompose them
    axRegistry.setDomain("X", data.nx, data.x_min, data.x_max);
    axRegistry.setDomain("Y", data.ny, data.y_min, data.y_max);
    axRegistry.setDomain("Z", data.nz, data.z_min, data.z_max);
    //if (data.ion_rec_nlevel) axRegistry.setDomain("nlevels", 5, 1, 5);

    mpi.decomposeAllAxes();

    //Allocate all registered variables
    auto& varRegistry = SAMS::getvariableRegistry();
    auto& axisRegistry = SAMS::getaxisRegistry();
    varRegistry.allocateAll();


    //Tell LARE to grab the shared allocated variables
	S.allocate(data,dataNeutral);
    //Tell LARE to set up its grid
    S.grid(data);
    if (data.two_fluid) {
        printf("Initialising two-fluid grid \n");
        S.two_fluid_grid(data,dataNeutral);
        if (data.ion_rec_nlevel){
            printf("Getting rates \n");
            S.two_fluid_read_rates(data);
        }
    }

		portableWrapper::fence();
    S.initial_conditions(data,dataNeutral);

	portableWrapper::fence();
    S.boundary_conditions(data);
    portableWrapper::fence();
    S.boundary_conditions(dataNeutral);
    portableWrapper::fence();
    timer t;
    t.begin("Main Loop");
    data.step=0;
    //data.time=0.0;

//printf("%f %f \n",dataNeutral.rho(-1,-1,-1),dataNeutral.energy_neutral(-1,-1,-1));

    while (true)
    {

      SAMS::cout << data.step << " " << data.time << std::endl;      
      if (data.step%1==0) S.output(data,dataNeutral);
      printf("nsteps,steps,time = %ld,%ld,%f \n",data.nsteps,data.step,data.time);

      if ((data.step >= data.nsteps && data.nsteps >= 0) || (data.time >= data.t_end))
        break;
      S.set_dt(data); // timestep of fluid
      if (data.two_fluid) {
        S.set_dt(dataNeutral); //Get fluid timestep of neutrals
        S.two_fluid_source(data,dataNeutral,true); // First step of Strang-split two-fluid sources
      }
      S.lagrangian_step(data,dataNeutral);    // lagran.cpp

      S.eulerian_remap(data); // remap.cpp
      if (data.two_fluid) S.eulerian_remap(dataNeutral);
      data.step++;
      if (data.two_fluid) dataNeutral.step++;
      if (data.rke) {
        S.energy_correction(data); // diagnostics.cpp
        if (data.two_fluid) S.energy_correction(dataNeutral);
      }
      S.eta_calc(data);            // lagran.cpp
      if (data.two_fluid) {
          S.two_fluid_source(data,dataNeutral,false); // Second step of Strang-split two-fluid sources
      }
    }
    t.end();

		S.output(data,dataNeutral);

		S.manager.clear();
    axisRegistry.finalize();
    varRegistry.finalize();
    portableWrapper::finalize();
    SAMS::MPI::finalize();

}
