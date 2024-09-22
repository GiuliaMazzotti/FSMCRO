# FSMCRO - Combining the FSM2oshd model and the externalized Crocus model into one system

- FSM_SOURCE_CODE contains the source FORTRAN code for FSM.

- SHELL_SCRIPTS contains wrapper scripts

- OPTIONS_NAM contains default namelist files 

## Dependencies

This code is based on two 

## How to run the models

### Installation

Clone this repository:

```
git clone https://github.com/oshd-slf/jim_operational.git
```

### Compiling FSM
Install a 32bit version of MinGW, with gfortran >=6.3.0.
The batch file for compiling FSM is in `FSM_SOURCE_FOLDER/code/compil_FSM.bat`. It handles two *positional* arguments:
 * default behavior is with optimization level `-O3`
 * _gfortran optimization_ of computations (no argument : default behavior, poor performance but very verbose error messages, `-O3` offers the best performance but the least verbosity, `-O2` intermediate performance and verbosity).
 *  _gprof profiling_ of execution as a second positional argument after optimizatino(to evaluate code performance): (type e.g. `-O0 -pg` if you want to activate profiling combined with low optimization).

__Recommandation for operational purposes__:
1. open a command window and move to FSM_SOURCE_CODE/code/
2. `compil_FSM.bat`

### Running FSM
Run FSM using one of the start scripts in `MATLAB_SCRIPTS\SCRIPTS\` (e.g. `start_POINT.m`, `start_SPATIAL.m`).

#### Namelist
Namelists are used to set the model configuration, number of points, and list of output variables.
Default operational namelists are stored in `FSM_SOURCE_CODE\nslt\`.

#### Input data
* forcings: @TODO
* terrain/landuse: stored in `SOURCE\`
* initial snow conditions: @TODO

#### Output variables
The wrapper writes input and output variables (state variables and diagnostics) into `MODELDATA_*_FSM22.mat files`. The operational lists of output variables are written in `SOURCE\Simulation_Settings_FSM.m` and can be overwritten in any of the start scripts. Description fields for the variables can be accessed via MATLAB_SCRIPT/varname_to_desc.m (see doc in header).

## Copy forcing data to local drive

The models run faster when reading forcing data from a local drive. Below is an example script using `robocopy` for fast and incremental copying of files from network drives to a local folder. Put these commands in a batch-file.

```
rem @echo off

Robocopy I:\DATA_COSMO\PROCESSED_GRID_ANALYSIS D:\DATA_COSMO\PROCESSED_GRID_ANALYSIS /E
Robocopy K:\DATA_COSMO\PROCESSED_GRID_ANALYSIS D:\DATA_COSMO\PROCESSED_GRID_ANALYSIS /E

Robocopy I:\DATA_COSMO\PROCESSED_STAT_ANALYSIS D:\DATA_COSMO\PROCESSED_STAT_ANALYSIS /E
Robocopy K:\DATA_COSMO\PROCESSED_STAT_ANALYSIS D:\DATA_COSMO\PROCESSED_STAT_ANALYSIS /E

Robocopy I:\DATA_COSMO\PROCESSED_STAT_FARCHIVE D:\DATA_COSMO\PROCESSED_STAT_FARCHIVE /E
Robocopy K:\DATA_COSMO\PROCESSED_STAT_FARCHIVE D:\DATA_COSMO\PROCESSED_STAT_FARCHIVE /E

Robocopy I:\DATA_COSMO\PROCESSED_STAT_FORECAST D:\DATA_COSMO\PROCESSED_STAT_FORECAST /E
Robocopy K:\DATA_COSMO\PROCESSED_STAT_FORECAST D:\DATA_COSMO\PROCESSED_STAT_FORECAST /E

pause
```

# FSM-CRO research model

This branch (originally: fsmcro_sandbox branch) includes the added functionality of using the Crocus snow model instead of the FSM snow routine. Originally intended for research project on hyper-res tree- and snow-layer resolving model (SNF mobility project Giulia)

Developments started in 2022, branching off from the 2022 operational repository version and the hyper-resolution modelling capabilities Maintained for use in publication 'Exploring the potential of forest snow modelling at the tree and snowpack layer scale', in prep. 

## Technicalities: coupling FSM and Crocus codes

### Use of Externalized Crocus 

A standalone version of Crocus, intended to allow its use independent of the SURFEX platform, is provided by CEN/MeteoFrance. The code comes with the necessary dependencies and a makefile, which creates an executable aimed at testing the correct installation of the model. Running the makefile compiles the Crocus code, creating a set of fortran object files (*.o) stored in the FSM_SOURCE_CODE/obj subdirectory. The contents of this subdirectory can be directly (and manualy) copied to the equivalent /obj subdirectory in the FSM-CRO repository (i.e. the jim_operational repo) for further use by the FSM-CRO compiler (linking step creating the FSM-CRO executable). Watch out: the file prog.o needs to be deleted manually (or just not copied) because it will confuse the FSM-CRO compiler (it's a main program code, and there is another (the correct) main program in FSMCRO.F90, hence compilation will fail if prog.o is not removed)

This way of manually linking the FSM and Crocus repositories may seem cumbersome, but it's the simplest way to keep developments in the two models separate and independent. At every update of EXT-CROCUS in the Surfex repository, updates can be pulled and the standalone Crocus recompiled. Changes are ported to the FSM repository (jim_operational) by replacing contents of the ./obj subdirectory with their updated version and recompiling. Every update of the ./obj subdirectory can be pushed as commit in the jim_operational repository to keep track of these updates. Ideally, these commits should include ONLY the updates to the ./obj subdirectories, while potential required changes to the FSM-CRO code that ensure the compatibility with EXT-CROCUS updates should constitute follow-up commits.

### Compiling and running FSM-CRO

Compiling FSM-CRO requires its own compilation script, available as compil_FSMCRO.sh. This was separated from compil_FSM in the inital stages of development because it had to be structured a bit differently than the FSM compilation script Once FSM and FSM-CRO are fully integrated, this script could be obsolete. 

Consequently, the compil_FSMCRO.sh script compiles the program FSMCRO.exe rather than FSM2.exe. Also this could become obsolete at some point.

Note that the compilation of FSMCRO is done a bit differently than for FSM. in the FSM compilation script, the compilation and the linking steps are merged into one and the executable is created in the same command. In FSMCRO, since .o files are already available from EXT-CROCUS, we first create equivalent .o files from the FSM .F90 files, and then link all of them to create the FSMCRO.exe executable. 

FSMCRO is run the same way as FSM, by passing one namelist file. In a terminal: S ./FSMCRO OPTIONS.nam

### Fortran code changes

- FSMCRO.F90: added main program including call to CRO_SETUP
- CRO_MODUL: defines Crocus-specific modules
- CRO_SETUP: initialization of Crocus-specific stuff prior to timeloop

### FSM-CRO specific model switches (namelist options)

in namelist nam_modconf:
- CRO_ON: .FALSE. as default, set to .TRUE. to switch on Crocus snow
- MET1D_ON: .FALSE. as default, set to .TRUE. to apply point forcing data to entire array of modelled points (without modifications)
- TVT_ON: .TRUE. as default, reads canopy SWdir transmissivity time series, set to .FALSE. to apply constant transmissivity = sky-view fraction.

namelists nam_crooptions and nam_croparams:
include all switches required for the ESCROC ensemble runs; note that options are combinations of option tags and parameters, as defined in snowtools (see ESCROC capabilities below). namelists need to exist if CRO_ON is set; defaults corresponding to the operational Crocus settings are used if the namelists are empty (NOT the same as deterministic member in ESCROC subensembles)

Example Crocus namelists for open and forest runs are provided. 

## Additional scripts

TODO

### Wrappers

TODO

### ESCROC capabilities

The capability to run ensembles by combining alternative snow process/properties parametrizations options as in the ESCROC system is principally available for the externalized Crocus version as well. 

## Open issues - known pending developments

1. The Crocus snow model is usually run with an internal timestep that is smaller than the forcing timestep, typically 15min. This functionality has not yet been implemented here, so that the temporal resolution of FSM-CRO is currently determined by the forcing rather than being user-defined. The possibility to run shorter timesteps should be implemented (check ol_time_interp_at.f90 in SURFEX)

2. Crocus as used within Surfex is known to be sensitive to initial conditions of soil temperatures and moistures. This has not yet been tested, initializations are unchanged from FSM. 

Unloading snow is either added as rain or fresh snow because Crocus can only take that. Adding a unloading snow mass flux would be interesting and would potentially have visible impacts on snowpack evolution. Axel started doing this in Surfex, reconsider once ported to EXT-CROCUS...

Some ESCROC ensemble options dont work yet: Flanner metamorphism, Tartes...


Cosmetics: could organize the writing of Crocus-specific outputs and layer-scale outputs better. Also, there is currently no possibility to start a simulation from given initial Crocus-specific states (or existing FSM snowpack)

No fractional snow in Crocus
