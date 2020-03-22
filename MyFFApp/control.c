/************************* control.c *******************************/
/* MIMD version 7 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */
/* This file is for lattice generation with the RHMC algorithm */
#define CONTROL
#include "ks_imp_includes.h"	/* definitions files and prototypes */
#include "lattice_qdp.h"

#ifdef HAVE_QUDA
#include <quda_milc_interface.h>
#endif
#ifdef HAVE_QPHIX
#include "../include/generic_qphix.h"
#endif

#ifdef MILC_GLOBAL_DEBUG
#include "debug.h"
#endif /* MILC_GLOBAL_DEBUG */

/* For information */
#define NULL_FP -1

EXTERN gauge_header start_lat_hdr;	// Input gauge field header 
#include <stdio.h>
#include <stdlib.h>
#include<unistd.h>
#include <sys/stat.h>
int is_file_exist(const char *fileName)
{
  if(!access(fileName, F_OK )){
    printf("The File %s\t was Found\n",fileName); return 1;
  }else {printf("The File %s\t was not Found\n",fileName); return 0;}
}
/*
int is_file_exist(const char* filename)
{
  struct stat buffer;  int exist = stat(filename,&buffer);
  if(exist == 0) return 1; //File exists
  else return 0;  //File not exists
}
*/
int main( int argc, char **argv )
{
  int ComputePLoopFreeEnergy=1;
  int ComputeTraceFmunu =1;
  int SaveLattice=0; int UseSavedConfiguration=1;
  char InputDataFileDIR[100000], OutputDataFileDIR[100000], SaveLatticeDataFileDIR[100000], ReadLatticeDataFileDIR[100000];
  //sprintf(InputDataFileDIR,"%s",argv[1]);
  if(SaveLattice==1 && UseSavedConfiguration==0 ){sprintf(SaveLatticeDataFileDIR,"%s",argv[4]); }
  if(SaveLattice==0 && UseSavedConfiguration==1 ){sprintf(ReadLatticeDataFileDIR,"%s",argv[4]); }
  sprintf(OutputDataFileDIR, "%s", argv[5]);
  //sprintf(OutputDataFileDIR, "OutputTest");
  int FolderNumber=0;
  int i, MeasurementCount, traj_done, naik_index;
  int prompt;
  int s_iters=0, iters=0;
  double dtime, dclock();

  //Plaquette and Field-strength variable
  int Nc=3;//Color factor
  double SS_Plaq=0.0, ST_Plaq=0.0;
  double Current_Plaq=0.0, Sum_Plaq=0.0, Average_Plaq=0.0;
  double TadpoleFactor=0.0, TadpoleFactorNt=0.0; 
  complex CurrentPolyakovLoop, SumPolyakovLoop, AveragePolyakovLoop;
  double CurrentModPolyakovLoop=0.0, SumModPolyakovLoop=0.0, AverageModPolyakovLoop=0.0;
  double SumModPolyakovLoopTadpoleCorrected=0.0, AverageModPolyakovLoopTadpoleCorrected=0.0;
  double CurrentBareFreeEnergy=0.0,  SumBareFreeEnergy=0.0, AverageBareFreeEnergy=0.0;
  double CurrentBareFreeEnergyTadpoleCorrected=0.0,  SumBareFreeEnergyTadpoleCorrected=0.0, AverageBareFreeEnergyTadpoleCorrected=0.0;

  complex CurrentTraceF3iF3iMinusF4iF4i, SumTraceF3iF3iMinusF4iF4i, AverageTraceF3iF3iMinusF4iF4i;
  complex CurrentTraceF4iF3iPlusF3iF4i,  SumTraceF4iF3iPlusF3iF4i, AverageTraceF4iF3iPlusF3iF4i;
  complex CurrentTraceF3iDzF3iMinusF4iDzF4i, SumTraceF3iDzF3iMinusF4iDzF4i, AverageTraceF3iDzF3iMinusF4iDzF4i;
  complex CurrentTraceF4iDzF3iPlusF3iDzF4i,  SumTraceF4iDzF3iPlusF3iDzF4i,  AverageTraceF4iDzF3iPlusF3iDzF4i;  
  complex CurrentTraceF3iD2zF3iMinusF4iD2zF4i, SumTraceF3iD2zF3iMinusF4iD2zF4i, AverageTraceF3iD2zF3iMinusF4iD2zF4i;
  complex CurrentTraceF4iD2zF3iPlusF3iD2zF4i,  SumTraceF4iD2zF3iPlusF3iD2zF4i,  AverageTraceF4iD2zF3iPlusF3iD2zF4i;
  //Initialize variable to zero
  CurrentPolyakovLoop=cmplx(0.0,0.0); SumPolyakovLoop=cmplx(0.0,0.0); AveragePolyakovLoop=cmplx(0.0,0.0);
  CurrentTraceF3iF3iMinusF4iF4i =cmplx(0.0,0.0); CurrentTraceF4iF3iPlusF3iF4i =cmplx(0.0,0.0);
  CurrentTraceF3iDzF3iMinusF4iDzF4i = cmplx(0.0,0.0); CurrentTraceF4iDzF3iPlusF3iDzF4i=cmplx(0.0,0.0);
  CurrentTraceF3iD2zF3iMinusF4iD2zF4i= cmplx(0.0,0.0); CurrentTraceF4iD2zF3iPlusF3iD2zF4i= cmplx(0.0,0.0);
  SumTraceF3iF3iMinusF4iF4i =cmplx(0.0,0.0); AverageTraceF3iF3iMinusF4iF4i =cmplx(0.0,0.0);
  SumTraceF4iF3iPlusF3iF4i  =cmplx(0.0,0.0); AverageTraceF4iF3iPlusF3iF4i =cmplx(0.0,0.0);
  SumTraceF3iDzF3iMinusF4iDzF4i = cmplx(0.0,0.0); AverageTraceF3iDzF3iMinusF4iDzF4i= cmplx(0.0,0.0);
  SumTraceF4iDzF3iPlusF3iDzF4i  = cmplx(0.0,0.0); AverageTraceF4iDzF3iPlusF3iDzF4i = cmplx(0.0,0.0);
  SumTraceF3iD2zF3iMinusF4iD2zF4i = cmplx(0.0,0.0); AverageTraceF3iD2zF3iMinusF4iD2zF4i= cmplx(0.0,0.0);
  SumTraceF4iD2zF3iPlusF3iD2zF4i  = cmplx(0.0,0.0); AverageTraceF4iD2zF3iPlusF3iD2zF4i = cmplx(0.0,0.0);
  //FileName to save observables
  FILE *fploop, *ftracefmunuLO, *ftracefmunuNLO, *ftracefmunuNNLO;
  char FileNamePloop[10000], FileNameTraceFmunu[10000], FileNameTraceFmunu2[10000], FileNameTraceFmunu3[10000], SaveLatticeFileName[100000], ReadLatticeFileName[100000], FolderName[10000];

  // Initialization 
  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1); 
  g_sync();
  
  /* set up lattice parameters */
  prompt = setup();

  printf("Amit MyFFApp/control.c prompt= %d \n",prompt);
  printf("Amit MyFFApp/control.c before while(readin(prompt)==0) called \n");
  /* loop over input sets */
  while( readin(prompt) == 0)
    {
      sprintf(FileNamePloop,"%s/DataPloopNt%d_Ns%d_Beta%.4f.txt", OutputDataFileDIR, nt, nx, beta);          
      sprintf(FileNameTraceFmunu,"%s/DataTraceFmunuLO_Clover_Traceless_Nt%d_Ns%d_Beta%.4f.txt", OutputDataFileDIR, nt, nx, beta);
      sprintf(FileNameTraceFmunu2,"%s/DataTraceFmunuNLO_Clover_Traceless_Nt%d_Ns%d_Beta%.4f.txt", OutputDataFileDIR, nt, nx, beta);
      sprintf(FileNameTraceFmunu3,"%s/DataTraceFmunuNNLO_Clover_Traceless_Nt%d_Ns%d_Beta%.4f.txt", OutputDataFileDIR, nt, nx, beta);
      fploop = fopen(FileNamePloop,"w");
      ftracefmunuLO  = fopen(FileNameTraceFmunu,"w");
      ftracefmunuNLO = fopen(FileNameTraceFmunu2,"w");
      ftracefmunuNNLO = fopen(FileNameTraceFmunu3,"w");
      fprintf(fploop,"#Beta=%.4f, ml=%.6f, ms=%.6f, u0=%.3f, Nt=%d, Ns=%d^3 \n", beta, dyn_mass[0],dyn_mass[1], u0, nt, nx);
      fprintf(fploop,"#Iters \t Current_Plaq \t AvgPlaq \t TadpoleFactor \t TadpoleFactorNt \t CurrentPolyakovLoop.real \t CurrentPolyakovLoop.imag \t  CurrentModPolyakovLoop \t AverageModPolyakovLoop \t AverageModPolyakovLoopTadpoleCorrected \t  CurrentBareFreeEnergy \t AvgBareFreeEnergy  \t CurrentBareFreeEnergyTadpoleCorrected \t AvgBareFreeEnergyTadpoleCorrected \n");
      fprintf(ftracefmunuLO,"#Beta=%.4f, ml=%.6f, ms=%.6f, u0=%.3f, Nt=%d, Ns=%d^3 \n", beta, dyn_mass[0],dyn_mass[1], u0, nt, nx);
      fprintf(ftracefmunuLO,"#Iters \t TraceF3iF3iMinusF4iF4i.real \t TraceF3iF3iMinusF4iF4i.imag \t AvgTraceF3iF3iMinusF4iF4i.real \t AvgTraceF3iF3iMinusF4iF4i.imag \t TraceF4iF3iPlusF3iF4i.real \t TraceF4iF3iPlusF3iF4i.imag \t AvgTraceF4iF3iPlusF3iF4i.real \t AvgTraceF4iF3iPlusF3iF4i.imag \n");
      fprintf(ftracefmunuNLO,"#Beta=%.4f, ml=%.6f, ms=%.6f, u0=%.3f, Nt=%d, Ns=%d^3 \n", beta, dyn_mass[0],dyn_mass[1], u0, nt, nx);
      fprintf(ftracefmunuNLO,"#Iters \t TraceF3iDzF3iMinusF4iDzF4i.real \t TraceF3iDzF3iMinusF4iDzF4i.imag \t AvgTraceF3iDzF3iMinusF4iDzF4i.real \t AvgTraceF3iDzF3iMinusF4iDzF4i.imag \t TraceF4iDzF3iPlusF3iDzF4i.real \t TraceF4iDzF3iPlusF3iDzF4i.imag \t AvgTraceF4iDzF3iPlusF3iDzF4i.real \t AvgTraceF4iDzF3iPlusF3iDzF4i.imag \n");
      fprintf(ftracefmunuNNLO,"#Beta=%.4f, ml=%.6f, ms=%.6f, u0=%.3f, Nt=%d, Ns=%d^3 \n", beta, dyn_mass[0],dyn_mass[1], u0, nt, nx);
      fprintf(ftracefmunuNNLO,"#Iters \t TraceF3iD2zF3iMinusF4iD2zF4i.real \t TraceF3iD2zF3iMinusF4iD2zF4i.imag \t AvgTraceF3iD2zF3iMinusF4iD2zF4i.real \t AvgTraceF3iD2zF3iMinusF4iD2zF4i.imag \t TraceF4iD2zF3iPlusF3iD2zF4i.real \t TraceF4iD2zF3iPlusF3iD2zF4i.imag \t AvgTraceF4iD2zF3iPlusF3iD2zF4i.real \t AvgTraceF4iD2zF3iPlusF3iD2zF4i.imag \n");

      /* perform warmup trajectories */
      #ifdef MILC_GLOBAL_DEBUG
      global_current_time_step = 0;
      #endif /* MILC_GLOBAL_DEBUG */
    
      // Start application timer
      dtime = -dclock();
      printf(" Amit MyFFApp/control.c inside while(readin(prompt)==0) \n");
      for( traj_done=0; traj_done < warms; traj_done++ )
      	{
	  //rephase(OFF);
	  //SS_Plaq=0.0; ST_Plaq=0.0;
	  //d_plaquette(&SS_Plaq, &ST_Plaq);
	  //printf("Amit MyFFApp/control.c Plaquette = (%e,%e)\n",SS_Plaq, ST_Plaq);	  
	  //rephase(ON);
      	  //update();	  	  
        } 
      iters=warms;
      node0_printf("default MyFFApp/control.c WARMUPS COMPLETED\n"); fflush(stdout);
      
      /* perform measuring trajectories, reunitarizing and measuring 	*/
      MeasurementCount=0;		/* number of measurements 		*/            
      for( traj_done=0; traj_done < trajecs; traj_done++ )
	{ 
          #ifdef MILC_GLOBAL_DEBUG
          #ifdef HISQ_REUNITARIZATION_DEBUG
	  {
	    int isite, idir;
	    site *s;
	    FORALLSITES(isite,s)
	      {
		for( idir=XUP;idir<=TUP;idir++ )
		  {
		    lattice[isite].on_step_Y[idir] = 0;
		    lattice[isite].on_step_W[idir] = 0;
		    lattice[isite].on_step_V[idir] = 0;
		  }
	      }
	  }
          #endif /* HISQ_REUNITARIZATION_DEBUG */
          #endif /* MILC_GLOBAL_DEBUG */
	  
	  /* do the Measurement and then the trajectories */
	  /* measure every "propinterval" trajectories */
	  iters = 1 + iters;
	  MeasurementCount = MeasurementCount + 1;
	  if(traj_done==0){rephase(OFF);}
	  if(traj_done!=0 && UseSavedConfiguration==0){rephase(OFF);}
	  if(iters ==500) 
	    {
	      MeasurementCount = 1; 
	      Sum_Plaq=0.0;
	      SumPolyakovLoop = cmplx(0.0,0.0);
	      SumModPolyakovLoop=0.0; SumModPolyakovLoopTadpoleCorrected=0.0;
	      SumBareFreeEnergy=0.0;  SumBareFreeEnergyTadpoleCorrected=0.0;
	      SumTraceF3iF3iMinusF4iF4i = cmplx(0.0,0.0);  SumTraceF4iF3iPlusF3iF4i = cmplx(0.0,0.0);
	      SumTraceF3iDzF3iMinusF4iDzF4i = cmplx(0.0,0.0); 
	      SumTraceF4iDzF3iPlusF3iDzF4i  = cmplx(0.0,0.0);
	      SumTraceF3iD2zF3iMinusF4iD2zF4i = cmplx(0.0,0.0);
              SumTraceF4iD2zF3iPlusF3iD2zF4i  = cmplx(0.0,0.0);
	    }
	  /* call gauge_variable  measuring routines */
	  /* Compute plaquette, Polyakov loop, bare free energy and display/save in screen/file */
	  SS_Plaq=0.0; ST_Plaq=0.0;
	  d_plaquette(&SS_Plaq, &ST_Plaq);
	  Current_Plaq = ((SS_Plaq+ST_Plaq)/2.0)/Nc;
	  Sum_Plaq= Sum_Plaq + Current_Plaq;
	  Average_Plaq = Sum_Plaq/MeasurementCount;
	  TadpoleFactor = pow(Current_Plaq, 1.0/4.0);
	  TadpoleFactorNt = pow(Current_Plaq, nt/4.0);
	  printf("Amit MyFFApp/control.c beta=%.4f, u0=%.2f, Plaquette=(%e,%e), CurrentPlaq=%e, AvgPlaq=%e \n", beta, u0, SS_Plaq, ST_Plaq, Current_Plaq, Average_Plaq);
	  
	  /* Calculate trace of polyakov loop */
	  if(ComputePLoopFreeEnergy==1)
	    {
	      CurrentPolyakovLoop=cmplx(0.0,0.0); 
	      CurrentPolyakovLoop = ploop();
	      
	      CDIVREAL(CurrentPolyakovLoop,3.0,CurrentPolyakovLoop);								   
	      CADD(SumPolyakovLoop, CurrentPolyakovLoop, SumPolyakovLoop);
	      CDIVREAL(SumPolyakovLoop, MeasurementCount, AveragePolyakovLoop);
	      
	      complex *PointerPLoop = &CurrentPolyakovLoop;
	      CurrentModPolyakovLoop=cabs(PointerPLoop);
	      SumModPolyakovLoop = SumModPolyakovLoop + CurrentModPolyakovLoop;
	      SumModPolyakovLoopTadpoleCorrected = SumModPolyakovLoopTadpoleCorrected + (CurrentModPolyakovLoop*TadpoleFactorNt);
	      AverageModPolyakovLoop = SumModPolyakovLoop/MeasurementCount;
	      AverageModPolyakovLoopTadpoleCorrected = SumModPolyakovLoopTadpoleCorrected/MeasurementCount;
	      
	      CurrentBareFreeEnergy = -log(CurrentModPolyakovLoop);
	      SumBareFreeEnergy = SumBareFreeEnergy + CurrentBareFreeEnergy;
	      AverageBareFreeEnergy = SumBareFreeEnergy/MeasurementCount;
	      
	      CurrentBareFreeEnergyTadpoleCorrected = -log(TadpoleFactorNt*CurrentModPolyakovLoop);
	      SumBareFreeEnergyTadpoleCorrected = SumBareFreeEnergyTadpoleCorrected + CurrentBareFreeEnergyTadpoleCorrected;
	      AverageBareFreeEnergyTadpoleCorrected = SumBareFreeEnergyTadpoleCorrected/MeasurementCount;
	      printf("Amit MyFFApp/control.c PLoop=(%e,%e), AvgPLoop=(%e,%e), and (CurrentModPLOOP,AvgModPLOOP)=(%e,%e), and (CurrentBareFreeEnergy, AverageBareFreeEnergy)=(%e,%e), and (CurrenttBareFreeEnergyTadpoleCorr, AverageBareFreeEnergyTadpoleCorr)=(%e,%e)\n", CurrentPolyakovLoop.real, CurrentPolyakovLoop.imag, AveragePolyakovLoop.real, AveragePolyakovLoop.imag, CurrentModPolyakovLoop, AverageModPolyakovLoopTadpoleCorrected, CurrentBareFreeEnergy, AverageBareFreeEnergy, CurrentBareFreeEnergyTadpoleCorrected, AverageBareFreeEnergyTadpoleCorrected);
	      // write Plaquette, PLoop, Free Energy into a file 	      
	      fprintf(fploop,"%d \t %e \t %.4f  \t %.4f \t %.4f \t %.4f \t %.4f \t %e \t %.4f \t %.4f \t %e  \t %.4f  \t %e \t %.4f \n", iters, Current_Plaq,  Average_Plaq, TadpoleFactor, TadpoleFactorNt, CurrentPolyakovLoop.real, CurrentPolyakovLoop.imag, CurrentModPolyakovLoop, AverageModPolyakovLoop, AverageModPolyakovLoop,   CurrentBareFreeEnergy,  AverageBareFreeEnergy, CurrentBareFreeEnergyTadpoleCorrected, AverageBareFreeEnergyTadpoleCorrected);
	    } //end of ComputePLoopFreeEnergy if-condition
	  
	  if(ComputeTraceFmunu==1)
	    {
	      CurrentTraceF3iF3iMinusF4iF4i = cmplx(0.0,0.0);  CurrentTraceF4iF3iPlusF3iF4i = cmplx(0.0,0.0);
	      CurrentTraceF3iDzF3iMinusF4iDzF4i = cmplx(0.0,0.0); CurrentTraceF4iDzF3iPlusF3iDzF4i = cmplx(0.0,0.0);
	      CurrentTraceF3iD2zF3iMinusF4iD2zF4i = cmplx(0.0,0.0); CurrentTraceF4iD2zF3iPlusF3iD2zF4i = cmplx(0.0,0.0);
              fmunu_fmunu(&CurrentTraceF3iF3iMinusF4iF4i, &CurrentTraceF4iF3iPlusF3iF4i, &CurrentTraceF3iDzF3iMinusF4iDzF4i, &CurrentTraceF4iDzF3iPlusF3iDzF4i, &CurrentTraceF3iD2zF3iMinusF4iD2zF4i, &CurrentTraceF4iD2zF3iPlusF3iD2zF4i);

              CADD(SumTraceF3iF3iMinusF4iF4i, CurrentTraceF3iF3iMinusF4iF4i, SumTraceF3iF3iMinusF4iF4i);
              CADD(SumTraceF4iF3iPlusF3iF4i, CurrentTraceF4iF3iPlusF3iF4i, SumTraceF4iF3iPlusF3iF4i);
	      CADD(SumTraceF3iDzF3iMinusF4iDzF4i, CurrentTraceF3iDzF3iMinusF4iDzF4i, SumTraceF3iDzF3iMinusF4iDzF4i);
	      CADD(SumTraceF4iDzF3iPlusF3iDzF4i,  CurrentTraceF4iDzF3iPlusF3iDzF4i,  SumTraceF4iDzF3iPlusF3iDzF4i);
	      CADD(SumTraceF3iD2zF3iMinusF4iD2zF4i, CurrentTraceF3iD2zF3iMinusF4iD2zF4i, SumTraceF3iD2zF3iMinusF4iD2zF4i);
              CADD(SumTraceF4iD2zF3iPlusF3iD2zF4i,  CurrentTraceF4iD2zF3iPlusF3iD2zF4i,  SumTraceF4iD2zF3iPlusF3iD2zF4i);

              CDIVREAL(SumTraceF3iF3iMinusF4iF4i,     MeasurementCount, AverageTraceF3iF3iMinusF4iF4i);
              CDIVREAL(SumTraceF4iF3iPlusF3iF4i,      MeasurementCount, AverageTraceF4iF3iPlusF3iF4i);
	      CDIVREAL(SumTraceF3iDzF3iMinusF4iDzF4i, MeasurementCount, AverageTraceF3iDzF3iMinusF4iDzF4i);
	      CDIVREAL(SumTraceF4iDzF3iPlusF3iDzF4i,  MeasurementCount, AverageTraceF4iDzF3iPlusF3iDzF4i);
	      CDIVREAL(SumTraceF3iD2zF3iMinusF4iD2zF4i, MeasurementCount, AverageTraceF3iD2zF3iMinusF4iD2zF4i);
              CDIVREAL(SumTraceF4iD2zF3iPlusF3iD2zF4i,  MeasurementCount, AverageTraceF4iD2zF3iPlusF3iD2zF4i);

              printf("Amit MyFFApp/control.c TraceF3iF3iMinusF4iF4i=(%e,%e), AvgTrace=(%e,%e) \n", CurrentTraceF3iF3iMinusF4iF4i.real, CurrentTraceF3iF3iMinusF4iF4i.imag, AverageTraceF3iF3iMinusF4iF4i.real, AverageTraceF3iF3iMinusF4iF4i.imag);
              printf("Amit MyFFApp/control.c TraceF4iF3iPlusF3iF4i=(%e,%e), AvgTrace=(%e,%e) \n", CurrentTraceF4iF3iPlusF3iF4i.real, CurrentTraceF4iF3iPlusF3iF4i.imag, AverageTraceF4iF3iPlusF3iF4i.real, AverageTraceF4iF3iPlusF3iF4i.imag);

	      printf("Amit MyFFApp/control.c TraceF3iDzF3iMinusF4iDzF4i=(%e,%e), AvgTrace=(%e,%e) \n", CurrentTraceF3iDzF3iMinusF4iDzF4i.real, CurrentTraceF3iDzF3iMinusF4iDzF4i.imag, AverageTraceF3iDzF3iMinusF4iDzF4i.real, AverageTraceF3iDzF3iMinusF4iDzF4i.imag);
	      printf("Amit MyFFApp/control.c TraceF4iDzF3iPlusF3iDzF4i=(%e,%e), AvgTrace=(%e,%e) \n", CurrentTraceF4iDzF3iPlusF3iDzF4i.real, CurrentTraceF4iDzF3iPlusF3iDzF4i.imag, AverageTraceF4iDzF3iPlusF3iDzF4i.real, AverageTraceF4iDzF3iPlusF3iDzF4i.imag);

	      printf("Amit MyFFApp/control.c TraceF3iD2zF3iMinusF4iD2zF4i=(%e,%e), AvgTrace=(%e,%e) \n", CurrentTraceF3iD2zF3iMinusF4iD2zF4i.real, CurrentTraceF3iD2zF3iMinusF4iD2zF4i.imag, AverageTraceF3iD2zF3iMinusF4iD2zF4i.real, AverageTraceF3iD2zF3iMinusF4iD2zF4i.imag);
              printf("Amit MyFFApp/control.c TraceF4iD2zF3iPlusF3iD2zF4i=(%e,%e), AvgTrace=(%e,%e) \n", CurrentTraceF4iD2zF3iPlusF3iD2zF4i.real, CurrentTraceF4iD2zF3iPlusF3iD2zF4i.imag, AverageTraceF4iD2zF3iPlusF3iD2zF4i.real, AverageTraceF4iD2zF3iPlusF3iD2zF4i.imag);

	      fprintf(ftracefmunuLO,"%d \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n", iters, CurrentTraceF3iF3iMinusF4iF4i.real, CurrentTraceF3iF3iMinusF4iF4i.imag, AverageTraceF3iF3iMinusF4iF4i.real, AverageTraceF3iF3iMinusF4iF4i.imag, CurrentTraceF4iF3iPlusF3iF4i.real, CurrentTraceF4iF3iPlusF3iF4i.imag, AverageTraceF4iF3iPlusF3iF4i.real, AverageTraceF4iF3iPlusF3iF4i.imag );
	      fprintf(ftracefmunuNLO,"%d \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n", iters, CurrentTraceF3iDzF3iMinusF4iDzF4i.real, CurrentTraceF3iDzF3iMinusF4iDzF4i.imag, AverageTraceF3iDzF3iMinusF4iDzF4i.real, AverageTraceF3iDzF3iMinusF4iDzF4i.imag, CurrentTraceF4iDzF3iPlusF3iDzF4i.real, CurrentTraceF4iDzF3iPlusF3iDzF4i.imag, AverageTraceF4iDzF3iPlusF3iDzF4i.real, AverageTraceF4iDzF3iPlusF3iDzF4i.imag);
	      fprintf(ftracefmunuNNLO,"%d \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n", iters, CurrentTraceF3iD2zF3iMinusF4iD2zF4i.real, CurrentTraceF3iD2zF3iMinusF4iD2zF4i.imag, AverageTraceF3iD2zF3iMinusF4iD2zF4i.real, AverageTraceF3iD2zF3iMinusF4iD2zF4i.imag, CurrentTraceF4iD2zF3iPlusF3iD2zF4i.real, CurrentTraceF4iD2zF3iPlusF3iD2zF4i.imag, AverageTraceF4iD2zF3iPlusF3iD2zF4i.real, AverageTraceF4iD2zF3iPlusF3iD2zF4i.imag);
	    } // end of if-condition for trace of fmunu correlators 
	  

	  int FolderNumberIndex = (iters-1)/1000;
	  FolderNumber = 1000*(1 + FolderNumberIndex);

	  if( SaveLattice==1  )
	    { sprintf(FolderName,"mkdir %s/Nt%d_Ns%d/Beta%.4f_%d",SaveLatticeDataFileDIR, nt, nx, beta, FolderNumber);
	      if((FolderNumber-1000+1)==iters  && this_node==1)
                {system(FolderName);
                }
	      int flag=SAVE_SERIAL;
	      sprintf(SaveLatticeFileName,"%s/Nt%d_Ns%d/Beta%.4f_%d/Lattice_Nt%d_Ns%d_Beta%.4f_u0_%.3f.configuration.%d",SaveLatticeDataFileDIR, nt, nx, beta, FolderNumber, nt, nx, beta, u0, iters);
	      save_lattice( flag, SaveLatticeFileName, stringLFN );
	      //rephase( OFF );
	      // save_lattice( saveflag, savefile, stringLFN );
	      //rephase( ON );
	    }	      
	  
	  if(traj_done < trajecs - 1)
	    {
	      if(UseSavedConfiguration==0)
		{ rephase(ON);
		  printf(" Amit MyFFApp/control.c s_iters=update() called for beta= %.4f, u0=%.4f, at iters = %d, and FolderNumber=%d \n",beta, u0, iters, FolderNumber);
		  s_iters=update();
		}
	      else 
		{ int flag=RELOAD_SERIAL;		  
		  sprintf(ReadLatticeFileName,"%s/Nt%d_Ns%d/Beta%.4f_%d/Lattice_Nt%d_Ns%d_Beta%.4f_u0_%.3f.configuration.%d", ReadLatticeDataFileDIR, nt, nx, beta, FolderNumber, nt, nx, beta, u0, iters);
		  if(is_file_exist(ReadLatticeFileName)==0) 
		    {
		      while(is_file_exist(ReadLatticeFileName)==0 && traj_done < trajecs - 1)
			{
			  iters++; traj_done++;
			  FolderNumberIndex = (iters-1)/1000;
			  FolderNumber = 1000*(1 + FolderNumberIndex);
			  sprintf(ReadLatticeFileName,"%s/Nt%d_Ns%d/Beta%.4f_%d/Lattice_Nt%d_Ns%d_Beta%.4f_u0_%.3f.configuration.%d", ReadLatticeDataFileDIR, nt, nx, beta, FolderNumber, nt, nx, beta, u0, iters);
			}
		    }
		   reload_lattice( flag, ReadLatticeFileName);
		}
	    }
	}// end loop over trajectories 
      

      printf("default MyFFApp/control.c RUNNING COMPLETED, This node is %d \n",this_node); 
      fflush(stdout);	     
      dtime += dclock();	     
      printf("Default MyFFApp/control.c Time = %e seconds \n",dtime);
      printf("Default MyFFApp/control.c total_iters = %d \n",iters);      
      fflush(stdout);
      
      // save lattice if requested 
      if( saveflag != FORGET )
	{ rephase( OFF );
	  save_lattice( saveflag, savefile, stringLFN );
	  rephase( ON );
	}      

      // Destroy fermion links (created in readin() 
#if FERM_ACTION == HISQ
      destroy_fermion_links_hisq(fn_links);
#elif FERM_ACTION == HYPISQ
      destroy_fermion_links_hypisq(fn_links);
#else
      destroy_fermion_links(fn_links);
#endif
      fn_links = NULL;
    }
  
  fclose(fploop);
  fclose(ftracefmunuLO);
  fclose(ftracefmunuNLO);
  fclose(ftracefmunuNNLO);
  normal_exit(0);
  return 0;
}

