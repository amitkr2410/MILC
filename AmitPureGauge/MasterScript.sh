#!/bin/bash 

MyDIR=/wsu/home/fy/fy41/fy4125/Measurement/PureGaugeColverTraceless/AmitPureGauge
cd ${MyDIR}
NT=(4  6  8  16 24 32)
NS=(16 24 32 16 24 32)
Select=$1
Mem=$2
CPU_Model=$3
i=$4
for (( m=0; m<1; m++ ))
do
    for (( Trial=0; Trial<10; Trial++ )) #69
    do
    echo Submitting Job \# $Trial

    #LogFile=/wsu/home/fy/fy41/fy4125/Log/ScreenNewNt4Output${Trial}.out
    ErrFile=/wsu/home/fy/fy41/fy4125/Log/ScreenPureGaugeNt${NT[$i]}New${Trial}.err
    LogFile=/dev/null
    #ErrFile=/dev/null
    Exec=${MyDIR}/simple_job
# ./simple_job ${NT[$i]} ${NS[$i]}  ${Trial}  16 
    ##qsub -V -q wsuq accq mwsuq  -l mem=3gb -N DoQueue -o $LogFile -e $ErrFile -- $Exec $Args $Trial
#   qsub -V  -q wsuq -l mem=8gb -l ncpus=8 -l mpiprocs=8  -N t4PG$Trial  -o $LogFile -e $ErrFile --  $Exec $Trial

   ## wsu161-wsu184 (28cpus each node, 128GB)
   ##   qsub -V -q mwsuq  -l cpu_type=Intel -l cpu_model=E5-2697v3 -l mem=16gb -l ncpus=8 -l mpiprocs=8 -N Nt4i$Trial -o $LogFile -e $ErrFile --        $Exec $Trial

   qsub -V -q wsuq  -l select=${Select}:ncpus=1:mem=${Mem}gb:mpiprocs=1:cpu_model=${CPU_Model} -N Nt${NT[$i]}i$Trial -o $LogFile -e $ErrFile -- $Exec ${NT[$i]} ${NS[$i]}  ${Trial}  ${Select}
   #qsub -V   -l select=1:ncpus=16:mem=16gb:mpiprocs=16:cpu_type=Intel  -N Nt16i$Trial -o $LogFile -e $ErrFile --  $Exec $Trial
    done
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
