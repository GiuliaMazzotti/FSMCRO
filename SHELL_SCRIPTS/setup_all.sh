########################################################################
# Shell script to setup runs for FSM2-CROCUS - include update of all  #
# runs required for paper 
# Giulia Mazzotti, August 2023                                        #
########################################################################
set -x

# Usage 
#./setup_all.sh 0 0 0 0 0 0 (do nothing) or 
#./setup_all.sh 1 2 3 4 5 6 (execute all)
# $1 = 1: activate compilation
# $2 = 2: update namelists 
# $3 = 3: update open site point data
# $4 = 4: update forest points data laret
# $5 = 5: update forest points data sodankyla
# $6 = 6: update forest grids data laret 

# 1. Compile yes / no
if [ $1 -eq 1 ]; then
    cd ../FSM_SOURCE_CODE/code
    ./compil_FSMCRO.sh
    cd ../../TMP_RUN
fi 

# 2. Update namelists yes / no
if [ $2 -eq 2 ]; then
    cp ../OPTS_NAM/OPTIONS_CRO_*.nam .
fi 

# 3. Update open point data yes / no
if [ $3 -eq 3 ]; then
    folders=("FS0_sod" "FS0_lar")    
    for folder in "${folders[@]}"; do
        if [[ -d "$folder" ]]; then
            echo "$folder already exists"
        else
            mkdir $folder
        fi
    done 
    cp ../DATA_SOD/FS0/*.bin FS0_sod
    cp ../DATA_LAR/FS0/*.bin FS0_lar
fi
# 

# 4. Update forest point data laret yes / no
if [ $4 -eq 4 ]; then
    
    for ((i=1; i<=5; i++)); do
        folder="FS${i}_lar"
        if [[ -d "$folder" ]]; then
            echo "$folder already exists"
        else
            mkdir $folder
        fi
        cp ../DATA_LAR/FS0/drive*.bin $folder
        cp ../DATA_LAR/FS${i}/* $folder 
    done     

fi
# 

# 5. Update forest point data sodankyla yes / no
if [ $5 -eq 5 ]; then

    fsix=(4 5 6 7 9 10)
    for i in "${fsix[@]}"; do
        folder="FS${i}_sod"
        if [[ -d "$folder" ]]; then
            echo "$folder already exists"
        else
            mkdir $folder
        fi
        cp ../DATA_SOD/FS0/drive*.bin $folder
        cp ../DATA_SOD/FS${i}/* $folder 
    done  

fi


# 6. Update forest grid data laret yes / no
if [ $6 -eq 6 ]; then
    for ((i=1; i<=160; i++)); do
        folder="FG${i}_lar"
        if [[ -d "$folder" ]]; then
            echo "$folder already exists"
        else
            mkdir $folder
        fi
        cp ../DATA_LAR/FG0/drive*.bin $folder
        cp ../DATA_LAR/FG${i}/* $folder 
    done 
    
fi 

# TO DO: add the FMI sites at Sodankyla and the point with average forest properties at Laret -> not done for 