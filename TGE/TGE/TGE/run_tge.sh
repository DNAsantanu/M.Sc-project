#!/bin/bash
#PBS -o run_tge.out
#PBS -e errorfile.err
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=5:mem=5000mb 



cd ~/scratch/TGE/TGE/

#For data
#./tge ~/scratch/drift_ph2_day00/data_cutoff/1000003720.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003720/GV_data

#./tge ~/scratch/drift_ph2_day00/data_cutoff/1000003840.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003840/GV_data

#./tge ~/scratch/drift_ph2_day00/data_cutoff/1000003960.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003960/GV_data

#./tge ~/scratch/drift_ph2_day00/data_cutoff/1000004080.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004080/GV_data

#./tge ~/scratch/drift_ph2_day00/data_cutoff/1000004200.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004200/GV_data

#./tge ~/scratch/drift_ph2_day00/data_cutoff/1000004320.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004320/GV_data

#For running TGE on .fits files after running Mg 1000003720


#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-1_UAPS_1000003720_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003720/GV-mg-1

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-2_UAPS_1000003720_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003720/GV-mg-2

./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-3_UAPS_1000003720_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003720/GV-mg-3

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-4_UAPS_1000003720_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003720/GV-mg-4

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-5_UAPS_1000003720_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003720/GV-mg-5

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-6_UAPS_1000003720_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003720/GV-mg-6

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-7_UAPS_1000003720_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003720/GV-mg-7

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-8_UAPS_1000003720_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003720/GV-mg-8

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-9_UAPS_1000003720_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003720/GV-mg-9

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-10_UAPS_1000003720_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003720/GV-mg-10

#Running TGE on 1000003840

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-1_UAPS_1000003840_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003840/GV-mg-1

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-2_UAPS_1000003840_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003840/GV-mg-2

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-3_UAPS_1000003840_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003840/GV-mg-3

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-4_UAPS_1000003840_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003840/GV-mg-4

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-5_UAPS_1000003840_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003840/GV-mg-5

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-6_UAPS_1000003840_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003840/GV-mg-6

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-7_UAPS_1000003840_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003840/GV-mg-7

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-8_UAPS_1000003840_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003840/GV-mg-8

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-9_UAPS_1000003840_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003840/GV-mg-9

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-10_UAPS_1000003840_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003840/GV-mg-10

# Running TGE on 1000003960

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-1_UAPS_1000003960_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003960/GV-mg-1

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-2_UAPS_1000003960_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003960/GV-mg-2

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-3_UAPS_1000003960_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003960/GV-mg-3

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-4_UAPS_1000003960_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003960/GV-mg-4

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-5_UAPS_1000003960_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003960/GV-mg-5

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-6_UAPS_1000003960_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003960/GV-mg-6

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-7_UAPS_1000003960_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003960/GV-mg-7

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-8_UAPS_1000003960_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003960/GV-mg-8

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-9_UAPS_1000003960_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003960/GV-mg-9

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-10_UAPS_1000003960_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000003960/GV-mg-10


# Running TGE on 1000004080

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-1_UAPS_1000004080_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004080/GV-mg-1

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-2_UAPS_1000004080_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004080/GV-mg-2

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-3_UAPS_1000004080_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004080/GV-mg-3

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-4_UAPS_1000004080_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004080/GV-mg-4

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-5_UAPS_1000004080_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004080/GV-mg-5

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-6_UAPS_1000004080_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004080/GV-mg-6

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-7_UAPS_1000004080_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004080/GV-mg-7

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-8_UAPS_1000004080_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004080/GV-mg-8

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-9_UAPS_1000004080_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004080/GV-mg-9

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-10_UAPS_1000004080_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004080/GV-mg-10


# Running TGE on 1000004200


#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-1_UAPS_1000004200_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004200/GV-mg-1

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-2_UAPS_1000004200_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004200/GV-mg-2

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-3_UAPS_1000004200_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004200/GV-mg-3

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-4_UAPS_1000004200_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004200/GV-mg-4

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-5_UAPS_1000004200_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004200/GV-mg-5

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-6_UAPS_1000004200_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004200/GV-mg-6

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-7_UAPS_1000004200_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004200/GV-mg-7

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-8_UAPS_1000004200_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004200/GV-mg-8

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-9_UAPS_1000004200_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004200/GV-mg-9

#./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-10_UAPS_1000004200_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004200/GV-mg-10


# Running TGE on 1000004320

./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-1_UAPS_1000004320_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004320/GV-mg-1

./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-2_UAPS_1000004320_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004320/GV-mg-2

./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-3_UAPS_1000004320_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004320/GV-mg-3

./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-4_UAPS_1000004320_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004320/GV-mg-4

./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-5_UAPS_1000004320_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004320/GV-mg-5

./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-6_UAPS_1000004320_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004320/GV-mg-6

./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-7_UAPS_1000004320_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004320/GV-mg-7

./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-8_UAPS_1000004320_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004320/GV-mg-8

./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-9_UAPS_1000004320_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004320/GV-mg-9

./tge ~/scratch/drift_ph2_day00_sims/data_cutoff/f90_Nrel-10_UAPS_1000004320_uvw.fits inputparameters /lfs/usrhome/phd/ph22d020/scratch/TGE_output/1000004320/GV-mg-10



