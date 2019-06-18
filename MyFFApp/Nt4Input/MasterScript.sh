#!/bin/bash 

MyDIR=/wsu/home/fy/fy41/fy4125/Lattice/MILC/MyFFApp
cd ${MyDIR}

for (( Trial=0; Trial<2; Trial++ )) #69
do
    echo Submitting Job \# $Trial

    #LogFile=/wsu/home/fy/fy41/fy4125/Log/ScreenNewNt4Output${Trial}.out
    ErrFile=/wsu/home/fy/fy41/fy4125/Log/ScreenNt4New${Trial}.err
    LogFile=/dev/null
    #ErrFile=/dev/null
    Exec=${MyDIR}/simple_job

    ##qsub -V -q wsuq accq mwsuq  -l mem=3gb -N DoQueue -o $LogFile -e $ErrFile -- $Exec $Args $Trial
   ##qsub -V -q eamxq -l mem=32gb -l ncpus=8 -l mpiprocs=8  -N Nt4i$Trial  -o $LogFile -e $ErrFile --  $Exec $Trial

    ## wsu201-209 (40cpus each node, 1536 GBRAM)
      qsub -V -q wsuq -l vnode=wsu207 -l cpu_type=Intel -l cpu_model=E5-4627v4 -l mem=32gb -l ncpus=8 -l mpiprocs=8 -N Nt4i$Trial -o $LogFile -e $ErrFile --        $Exec $Trial

   ## wsu161-wsu184 (28cpus each node, 128GB)
   ##   qsub -V -q mwsuq  -l cpu_type=Intel -l cpu_model=E5-2697v3 -l mem=16gb -l ncpus=8 -l mpiprocs=8 -N Nt4i$Trial -o $LogFile -e $ErrFile --        $Exec $Trial
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
