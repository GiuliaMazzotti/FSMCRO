########################################################################
# Compilation script for combined FSM2-CROCUS                          #
# Giulia Mazzotti, 2022                                                #
########################################################################

# this compilation script assumes compiled .o files for all subroutines 
# / modules specific to ext-crocus exist in the ../obj subdirectory of 
# this codebase. update of the .obj files based on the latest releas of 
# EXT_CROCUS needs to be taken care of beforehand.
# this script is assumed to be located in the ./code 
# subfolder together with all FSM-CRO .F90 files. 


# 1. Compile FSM subroutines including the main programm and including files in the /obj subdir
gfortran -I../obj -g -fdefault-double-8 -fbounds-check -fcheck=all -c MODULES.F90 CRO_MODUL.F90 MODE_WRITE.F90 CANOPY.F90 CUMULATE.F90 DRIVE.F90 DUMP.F90 EBALFOR.F90 EBALSRF.F90 LUDCMP.F90 OPEN_FILES.F90 PHYSICS.F90 QSAT.F90 RADIATION.F90 SETUP.F90 SNOW.F90 SNOWCOVERFRACTION.F90 SOIL.F90 SFEXCH.F90 SUNPOS.F90 THERMAL.F90 TRIDIAG.F90 CRO_SETUP.F90 CRO_CHECKIN.F90 CRO_COUP.F90 FSMCRO.F90

# 2. Run the Linker
gfortran -o FSMCRO ../obj/*.o *.o 

# 3. Clean the code folder 
rm *.o
rm *.mod

# 4. Move the exe to the temproary run folder 
mv FSMCRO ../../TMP_RUN/FSMCRO

# reminder: to make this a runnable script (i.e. change file permissions), run: 
# chmod +x compil_FSMCRO.sh
