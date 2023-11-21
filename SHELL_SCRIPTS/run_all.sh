########################################################################
# Shell script to run FSM2-CROCUS - include option to run all types 
# required for paper 
# Giulia Mazzotti, August 2023                                        #
########################################################################
set -x

# Usage 
# to make this a runnable script: chmod +x run_all.sh
# $1 = 1: run open points 
# $2 = 2: run open ensemble 

# 1. run open points
if [ $1 -eq 1 ]; then
    folders=("FS0_sod" "FS0_lar")    
    for folder in "${folders[@]}"; do
        cp FSMCRO $folder
        cp OPTIONS_CRO_opn.nam $folder
        cd $folder
        ./FSMCRO OPTIONS_CRO_opn.nam
        rm states_out*
        cd ..
    done 
fi 

# 2. run open ensemble
if [ $2 -eq 2 ]; then
    folders=("FS0_sod" "FS0_lar")    
    for folder in "${folders[@]}"; do
        python ~/PycharmProjects/FSMCRO_scripts/setup_ensemble_point.py $folder
    done
fi

# 3. run laret and sodankyla deterministic grids
if [ $3 -eq 3 ]; then
    
    fsix=(4 5 6 7 9 10)
    fstot=11 
    for ((j=1; j<=$fstot; j++)); do
        if [[ $j -le 5 ]]; then
            i=$j
            folder="FS${i}_lar"
        else
            i=${fsix[$((j - 6))]}
            folder="FS${i}_sod"
        fi 

    cp FSMCRO $folder
    python ~/PycharmProjects/FSMCRO_scripts/edit_namelist_ptsnr.py $folder
    cd $folder
    ./FSMCRO OPTIONS_CRO_for.nam
    rm states_out*
    cd ..

        
    done

fi

# 4. run laret and sodankyla ensemble at selected grids
if [ $4 -eq 4 ]; then

    fsix=(6 7)
    fstot=2
    # quickfix to run sodankyla only, see AUGUST version for full set
    for ((j=1; j<=$fstot; j++)); do
        i=${fsix[$((j - 1))]}

        folder="FS${i}_sod"
    
        python ~/PycharmProjects/FSMCRO_scripts/setup_ensemble_point.py $folder

        
    done

fi

#5. run deterministic large grid 
if [ $5 -eq 5 ]; then

    for ((i=1; i<=160; i++)); do
        folder="FG${i}_lar"

        cp FSMCRO $folder
        python ~/PycharmProjects/FSMCRO_scripts/edit_namelist_ptsnr.py $folder
        cd $folder
        ./FSMCRO OPTIONS_CRO_for.nam

        rm states_out*

        for ((l=6; l<=50; l++)); do
            rm res_*${l}.bin
        done 

        if [[ $i -gt 1 ]]; then
            rm res_year.bin
            rm res_month.bin
            rm res_day.bin
            rm res_hour.bin
        fi
        
        rm res_lwtr.bin
        rm res_ro*.bin
        rm res_sbsc.bin
        rm res_scfe.bin
        rm res_slqt.bin
        rm res_swt*.bin
        rm res_tsfe.bin
        rm res_tsoil*.bin
        cd ..

    done
fi

# TO DO: add the FMI sites at Sodankyla and the point with average forest properties at Laret