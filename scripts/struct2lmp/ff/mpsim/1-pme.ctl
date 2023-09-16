PROJECT               [p_name] 
FF                    /ul/tpascal/mpsim/AMBER_MBSCO.par
PERIODIC              T
*
STRUCTURE             [bgf_file] 
*
NB_METHOD              EWALD
COUL_SWITCH            2
EWALD_ACC              0.00000001 2.0 2
*
LEVEL                 2 
SPLINE_CUTOFF           9.0 10.0 10.0
*
NB_UPDATE_FREQ         5
CELL_REALLOC_FREQ      5
LOAD_BAL_FREQ          5
VEL_RESCALE_FREQ       0
VEL_START              NONMIN
DYN_TEMP               300 
* Time step in Femtosecond
DYN_TIME_STEP          1
DYN_STEPS              10
NOSE                   tau 0.001
RESET_GLOBAL_V     T
*
MIN_RMS_DESIRED         0.1
MIN_MAX_STEPS           50
FINAL_BGF
*
ACTION                 ONE_EF
*
BGF_TRAJ T
BGF_TRAJ_VEL T
BGF_TRAJ_FREQ 1  
FINAL_ENER
AVERAGE_FREQ 10
*
SETUP_EEX
DO
INFO
EXIT
