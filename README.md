# FSMCRO - Combining the FSM2oshd model and the externalized Crocus model into one system

## General

This repository includes the model FSMCRO, presented in the following publication: 

Mazzotti, G., Nousu, J.-P., Vionnet, V., Jonas, T., Nheili, R., and Lafaysse, M.: Exploring the potential of forest snow modelling at the tree and snowpack layer scale, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2023-2781, 2023

The paper has been accepted for publication in The Cryosphere, an updated reference will be provided once available. Please refer to this article when using the model. 

## Contents

- FSM_SOURCE_CODE contains the source FORTRAN code for FSM.

- SHELL_SCRIPTS contains wrapper scripts created to run the code on Linux

- OPTIONS_NAM contains default namelist files 

## Dependencies

As described in the paper, the FSMCRO model is a merge of two snow models, FSM2 and Crocus. 
The FSM2 code is based on: https://github.com/oshd-slf/FSM2oshd
The Crocus standalone version is available at: https://opensource.umr-cnrm.fr/projects/surfex_git2/wiki/Install_standalone_version_of_Crocus#:~:text=Install%20standalone%20version%20of%20Crocus%C2%B6%20For

This codebase is self-sufficient, but future users should keep an eye on developments in both repositories. In particular, updates to the standalone Crocus version can be ported to FSMCRO by replacing the precompiled files in obj. 
When doing so, it the necessity for changes to the wrapper scripts should be checked. This repo branched off SLF-OSHD's operational repository in 2022, see branch 'fsmcro_sandbox'.

The wrapper scripts provided as examples here rely on functions from the snowtools repository: https://github.com/UMR-CNRM/snowtools

## How to run the models

### Installation

Clone this repository:

```
git clone https://github.com/GiuliaMazzotti/FSMCRO.git
```

### Compiling FSMCRO
Compilation requires gfortran and has been tested on Linux. 
Run the script compil_FSMCRO.sh found in 'FSM_SOURCE_FOLDER/code'.
This compilation step also links the precompiled obj files that constitute the standalone Crocus.

In case users wish to update FSMCRO with the latest version of Crocus, the instructions on the Crocus repository should be followed. 
Potentially, this will require adaptations to the coupling subroutines in the FSMCRO fortran code (and recompilation)

### Running FSM
FSMCRO can be run directly from the command line by calling ./FSMCRO OPTIONS.nam, where 'OPTIONS.nam' is the filename of the namelist to be used. 
It is recommended to create a run directory including all meteo input files, landuse files, the namelist and the program. By default, output files will be written in the same directory as the program.
The models run faster when reading forcing data from a local drive.  

#### Namelist
Namelists are used to set the model configuration, number of points, and list of output variables.
Default namelists for runs in forest and open terrains are stored in OPTIONS_CRO_for.nam

#### Input data
* forcings: provided as .bin files, one variable per file 
* terrain/landuse: provided as .bin files, one variable per file 
* initial snow conditions: the snowpack is renitialized in each simulation. reinitialization from existing snowpack pending. 

#### Output variables
The wrapper writes input and output variables (state variables and diagnostics) into .bin files

## Technicalities: coupling FSM and Crocus codes

### Use of Externalized Crocus 

A standalone version of Crocus, intended to allow its use independent of the SURFEX platform, is provided by CEN/MeteoFrance. The code comes with the necessary dependencies and a makefile, which creates an executable aimed at testing the correct installation of the model. Running the makefile compiles the Crocus code, creating a set of fortran object files (*.o) stored in the FSM_SOURCE_CODE/obj subdirectory. The contents of this subdirectory can be directly (and manualy) copied to the equivalent /obj subdirectory in the FSM-CRO repository (i.e. the jim_operational repo) for further use by the FSM-CRO compiler (linking step creating the FSM-CRO executable). Watch out: the file prog.o needs to be deleted manually (or just not copied) because it will confuse the FSM-CRO compiler (it's a main program code, and there is another (the correct) main program in FSMCRO.F90, hence compilation will fail if prog.o is not removed)

This way of manually linking the FSM and Crocus repositories may seem cumbersome, but it's the simplest way to keep developments in the two models separate and independent. At every update of EXT-CROCUS in the Surfex repository, updates can be pulled and the standalone Crocus recompiled. Changes are ported to the FSM repository (jim_operational) by replacing contents of the ./obj subdirectory with their updated version and recompiling. Every update of the ./obj subdirectory can be pushed as commit in the jim_operational repository to keep track of these updates. Ideally, these commits should include ONLY the updates to the ./obj subdirectories, while potential required changes to the FSM-CRO code that ensure the compatibility with EXT-CROCUS updates should constitute follow-up commits.

### Compiling FSM-CRO

Compiling FSM-CRO requires its own compilation script, available as compil_FSMCRO.sh. This was separated from compil_FSM in the inital stages of development because it had to be structured a bit differently than the FSM compilation script Once FSM and FSM-CRO are fully integrated, this script could be obsolete. Consequently, the compil_FSMCRO.sh script compiles the program FSMCRO.exe rather than FSM2.exe. 

Note that the compilation of FSMCRO is done a bit differently than for FSM. in the FSM compilation script, the compilation and the linking steps are merged into one and the executable is created in the same command. In FSMCRO, since .o files are already available from EXT-CROCUS, we first create equivalent .o files from the FSM .F90 files, and then link all of them to create the FSMCRO.exe executable. 


### Fortran code changes

- FSMCRO.F90: added main program including call to CRO_SETUP (used instead of the original FSM.F90)
- CRO_MODUL: defines Crocus-specific modules
- CRO_SETUP: initialization of Crocus-specific stuff prior to timeloop
- CRO_COUP: subroutine coupling the standalone Crocus to the FSM model. 

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
