#!/bin/bash 

MyDIR=/wsu/home/fy/fy41/fy4125/Lattice/milc/MyFFApp
cd ${MyDIR}

for (( Trial=1; Trial<6; Trial++ )) #69
do
    echo Submitting Job \# $Trial

    #LogFile=/wsu/home/fy/fy41/fy4125/Log/ScreenOutput${Trial}.out
    ErrFile=/wsu/home/fy/fy41/fy4125/Log/Screen${Trial}.err
    LogFile=/dev/null
    #ErrFile=/dev/null
    Exec=${MyDIR}/simple_job

    ##qsub -V -q wsuq accq mwsuq  -l mem=3gb -N DoQueue -o $LogFile -e $ErrFile -- $Exec $Args $Trial
    ## qsub -V -q eamxq -l mem=2gb   -N AmitQueue  -o $LogFile -e $ErrFile --  $Exec $Trial
    qsub -V -q wsuq  -l mem=32gb -l ncpus=32 -l mpiprocs=32 -N HISQ$Trial   -o $LogFile -e $ErrFile --  $Exec $Trial
    
done

##qsub -V -I -q eamxq -- /wsu/home/fy/fy41/fy4125/RUN/PP1/simple_job 0




##      Script is submitted to this Queue: 
## To delete all the jobs given by user say fy4125
## then type in terminal " qselect -u fy4125 | xargs qdel

## To delete a process with job ID say 458393
## then type in terminal "qdel 458393" 
##   qstat -Q    "to know jobs running on different queue"
##   wsuq mwsuq eamxq mwsuq zflq ezfhq
##PBS  -q eamxq                 
##PBS -q mwsuq
                                                                             
##PBS -q zfhq 
                                                                             
##      One core and 1GB of RAM selected:
                                                  

##PBS -l select=1:ncpus=1:mem=10GB
                                                         
##PBS -l select=1:ncpus=4:mem=60GB:cpu_speed=3.3 
                                          
##PBS -l select=1:ncpus=1:mem=50GB:cpu_speed=3.0
                                           
## vmem=25GB vnode=amx3 
                                                                   

##export DISPLAY=localhost:0.0 
                                                            
##ulimit -s unlimited 
