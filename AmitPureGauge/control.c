/************************ control.c ******************************/
/* MIMD version 7 */
/* Main procedure for pure gauge SU3 */

/* This version combines code for the refreshed molecular dynamics
   algorithm with the hybrid Monte Carlo algorithm, for which HMC_ALGORITHM 
   should be defined.  (Actually, the changes to control.c are minimal
   and the real differences will appear in update.c */

#define CONTROL
#include "pure_gauge_includes.h"

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
int main(int argc, char *argv[])  
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

  //int meascount,todo;
  //double dssplaq,dstplaq;
  //complex plp;

  int Nc=3.0;
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

  //Initialize variable to zero
  CurrentPolyakovLoop=cmplx(0.0,0.0); SumPolyakovLoop=cmplx(0.0,0.0); AveragePolyakovLoop=cmplx(0.0,0.0);
  CurrentTraceF3iF3iMinusF4iF4i =cmplx(0.0,0.0); CurrentTraceF4iF3iPlusF3iF4i =cmplx(0.0,0.0);
  SumTraceF3iF3iMinusF4iF4i =cmplx(0.0,0.0); AverageTraceF3iF3iMinusF4iF4i =cmplx(0.0,0.0);
  SumTraceF4iF3iPlusF3iF4i  =cmplx(0.0,0.0); AverageTraceF4iF3iPlusF3iF4i =cmplx(0.0,0.0);
  SumTraceF3iDzF3iMinusF4iDzF4i = cmplx(0.0,0.0); AverageTraceF3iDzF3iMinusF4iDzF4i= cmplx(0.0,0.0);
  SumTraceF4iDzF3iPlusF3iDzF4i  = cmplx(0.0,0.0); AverageTraceF4iDzF3iPlusF3iDzF4i = cmplx(0.0,0.0);

  //FileName to save observables 
  FILE *fploop, *ftracefmunuLO, *ftracefmunuNLO;
  char FileNamePloop[10000], FileNameTraceFmunu[10000], FileNameTraceFmunu2[10000], SaveLatticeFileName[100000], ReadLatticeFileName[100000], FolderName[10000];

  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  g_sync();
  /* set up */
  prompt = setup();

  printf("Amit MyFFApp/control.c prompt= %d \n",prompt);
  printf("Amit MyFFApp/control.c before while(readin(prompt)==0) called \n");
  /* loop over input sets */
  while( readin(prompt) == 0)
    {
      sprintf(FileNamePloop,"%s/DataPloopNt%d_Ns%d_Beta%.4f.txt", OutputDataFileDIR, nt, nx, beta);          
      sprintf(FileNameTraceFmunu,"%s/DataTraceFmunuLO_Clover_Traceless_Nt%d_Ns%d_Beta%.4f.txt", OutputDataFileDIR, nt, nx, beta);
      sprintf(FileNameTraceFmunu2,"%s/DataTraceFmunuNLO_Clover_Traceless_Nt%d_Ns%d_Beta%.4f.txt", OutputDataFileDIR, nt, nx, beta);
      fploop = fopen(FileNamePloop,"w");
      ftracefmunuLO  = fopen(FileNameTraceFmunu,"w");
      ftracefmunuNLO = fopen(FileNameTraceFmunu2,"w");
      fprintf(fploop,"#Beta=%.4f, Nt=%d, Ns=%d^3 \n", beta, nt, nx);
      fprintf(fploop,"#Iters \t Current_Plaq \t AvgPlaq \t TadpoleFactor \t TadpoleFactorNt \t CurrentPolyakovLoop.real \t CurrentPolyakovLoop.imag \t  CurrentModPolyakovLoop \t AverageModPolyakovLoop \t AverageModPolyakovLoopTadpoleCorrected \t  CurrentBareFreeEnergy \t AvgBareFreeEnergy  \t CurrentBareFreeEnergyTadpoleCorrected \t AvgBareFreeEnergyTadpoleCorrected \n");
      fprintf(ftracefmunuLO,"#Beta=%.4f, Nt=%d, Ns=%d^3 \n", beta, nt, nx);
      fprintf(ftracefmunuLO,"#Iters \t TraceF3iF3iMinusF4iF4i.real \t TraceF3iF3iMinusF4iF4i.imag \t AvgTraceF3iF3iMinusF4iF4i.real \t AvgTraceF3iF3iMinusF4iF4i.imag \t TraceF4iF3iPlusF3iF4i.real \t TraceF4iF3iPlusF3iF4i.imag \t AvgTraceF4iF3iPlusF3iF4i.real \t AvgTraceF4iF3iPlusF3iF4i.imag \n");
      fprintf(ftracefmunuNLO,"#Beta=%.4f, Nt=%d, Ns=%d^3 \n", beta, nt, nx);
      fprintf(ftracefmunuNLO,"#Iters \t TraceF3iDzF3iMinusF4iDzF4i.real \t TraceF3iDzF3iMinusF4iDzF4i.imag \t AvgTraceF3iDzF3iMinusF4iDzF4i.real \t AvgTraceF3iDzF3iMinusF4iDzF4i.imag \t TraceF4iDzF3iPlusF3iDzF4i.real \t TraceF4iDzF3iPlusF3iDzF4i.imag \t AvgTraceF4iDzF3iPlusF3iDzF4i.real \t AvgTraceF4iDzF3iPlusF3iDzF4i.imag \n");

#ifdef FUZZ
      if(this_node==0)printf("Fat Polyakov loop paramter %f\n",ALPHA_FUZZ);
#endif
      
      /* perform warmup trajectories */
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
      MeasurementCount=0;		/* number of measurements      	*/            
      for( traj_done=0; traj_done < trajecs; traj_done++ )
	{ 	  
	  /* do the Measurement and then the trajectories */
	  /* measure every "propinterval" trajectories */
	  iters = 1 + iters;
	  MeasurementCount = MeasurementCount + 1;
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
	  printf("Amit MyFFApp/control.c beta=%.4f, Plaquette=(%e,%e), CurrentPlaq=%e, AvgPlaq=%e \n", beta, SS_Plaq, ST_Plaq, Current_Plaq, Average_Plaq);
	  
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
	      
	      printf("Amit MyFFApp/control.c PLoop=(%e,%e), AvgPLoop=(%e,%e), and (CurrentModPLOOP,AvgModPLOOP)=(%e,%e), and (CurrentBareFreeEnergy, AverageBareFreeEnergy)=(%e,%e), and (CurrenttBareFreeEnergyTadpoleCorr, AverageBareFreeEnergyTadpoleCorr)=(%e,%e)\n", CurrentPolyakovLoop.real, CurrentPolyakovLoop.imag, AveragePolyakovLoop.real, AveragePolyakovLoop.imag, CurrentModPolyakovLoop, AverageModPolyakovLoopTadpoleCorrected, CurrentBareFreeEnergy, AverageBareFreeEnergy, CurrentBareFreeEnergyTadpoleCorrected, AverageBareFreeEnergyTadpoleCorrected);
	      // write Plaquette, PLoop, Free Energy into a file 	      
	      fprintf(fploop,"%d \t %e \t %.4f  \t %.4f \t %.4f \t %.4f \t %.4f \t %e \t %.4f \t %.4f \t %e  \t %.4f  \t %e \t %.4f \n", iters, Current_Plaq,  Average_Plaq, TadpoleFactor, TadpoleFactorNt, CurrentPolyakovLoop.real, CurrentPolyakovLoop.imag, CurrentModPolyakovLoop, AverageModPolyakovLoop, AverageModPolyakovLoop,   CurrentBareFreeEnergy,  AverageBareFreeEnergy, CurrentBareFreeEnergyTadpoleCorrected, AverageBareFreeEnergyTadpoleCorrected);
		} //end of ComputePLoopFreeEnergy if-condition
	  
	  if(ComputeTraceFmunu==1)
	    {
	      CurrentTraceF3iF3iMinusF4iF4i = cmplx(0.0,0.0);  CurrentTraceF4iF3iPlusF3iF4i = cmplx(0.0,0.0);
              fmunu_fmunu(&CurrentTraceF3iF3iMinusF4iF4i, &CurrentTraceF4iF3iPlusF3iF4i, &CurrentTraceF3iDzF3iMinusF4iDzF4i, &CurrentTraceF4iDzF3iPlusF3iDzF4i);

              CADD(SumTraceF3iF3iMinusF4iF4i, CurrentTraceF3iF3iMinusF4iF4i, SumTraceF3iF3iMinusF4iF4i);
              CADD(SumTraceF4iF3iPlusF3iF4i, CurrentTraceF4iF3iPlusF3iF4i, SumTraceF4iF3iPlusF3iF4i);
	      CADD(SumTraceF3iDzF3iMinusF4iDzF4i, CurrentTraceF3iDzF3iMinusF4iDzF4i, SumTraceF3iDzF3iMinusF4iDzF4i);
	      CADD(SumTraceF4iDzF3iPlusF3iDzF4i,  CurrentTraceF4iDzF3iPlusF3iDzF4i,  SumTraceF4iDzF3iPlusF3iDzF4i);

              CDIVREAL(SumTraceF3iF3iMinusF4iF4i,     MeasurementCount, AverageTraceF3iF3iMinusF4iF4i);
              CDIVREAL(SumTraceF4iF3iPlusF3iF4i,      MeasurementCount, AverageTraceF4iF3iPlusF3iF4i);
	      CDIVREAL(SumTraceF3iDzF3iMinusF4iDzF4i, MeasurementCount, AverageTraceF3iDzF3iMinusF4iDzF4i);
	      CDIVREAL(SumTraceF4iDzF3iPlusF3iDzF4i,  MeasurementCount, AverageTraceF4iDzF3iPlusF3iDzF4i);

              printf("Amit MyFFApp/control.c TraceF3iF3iMinusF4iF4i=(%e,%e), AvgTrace=(%e,%e) \n", CurrentTraceF3iF3iMinusF4iF4i.real, CurrentTraceF3iF3iMinusF4iF4i.imag, AverageTraceF3iF3iMinusF4iF4i.real, AverageTraceF3iF3iMinusF4iF4i.imag);
              printf("Amit MyFFApp/control.c TraceF4iF3iPlusF3iF4i=(%e,%e), AvgTrace=(%e,%e) \n", CurrentTraceF4iF3iPlusF3iF4i.real, CurrentTraceF4iF3iPlusF3iF4i.imag, AverageTraceF4iF3iPlusF3iF4i.real, AverageTraceF4iF3iPlusF3iF4i.imag);

	      printf("Amit MyFFApp/control.c TraceF3iDzF3iMinusF4iDzF4i=(%e,%e), AvgTrace=(%e,%e) \n", CurrentTraceF3iDzF3iMinusF4iDzF4i.real, CurrentTraceF3iDzF3iMinusF4iDzF4i.imag, AverageTraceF3iDzF3iMinusF4iDzF4i.real, AverageTraceF3iDzF3iMinusF4iDzF4i.imag);
	      printf("Amit MyFFApp/control.c TraceF4iDzF3iPlusF3iDzF4i=(%e,%e), AvgTrace=(%e,%e) \n", CurrentTraceF4iDzF3iPlusF3iDzF4i.real, CurrentTraceF4iDzF3iPlusF3iDzF4i.imag, AverageTraceF4iDzF3iPlusF3iDzF4i.real, AverageTraceF4iDzF3iPlusF3iDzF4i.imag);

	      fprintf(ftracefmunuLO,"%d \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n", iters, CurrentTraceF3iF3iMinusF4iF4i.real, CurrentTraceF3iF3iMinusF4iF4i.imag, AverageTraceF3iF3iMinusF4iF4i.real, AverageTraceF3iF3iMinusF4iF4i.imag, CurrentTraceF4iF3iPlusF3iF4i.real, CurrentTraceF4iF3iPlusF3iF4i.imag, AverageTraceF4iF3iPlusF3iF4i.real, AverageTraceF4iF3iPlusF3iF4i.imag );
	      fprintf(ftracefmunuNLO,"%d \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n", iters, CurrentTraceF3iDzF3iMinusF4iDzF4i.real, CurrentTraceF3iDzF3iMinusF4iDzF4i.imag, AverageTraceF3iDzF3iMinusF4iDzF4i.real, AverageTraceF3iDzF3iMinusF4iDzF4i.imag, CurrentTraceF4iDzF3iPlusF3iDzF4i.real, CurrentTraceF4iDzF3iPlusF3iDzF4i.imag, AverageTraceF4iDzF3iPlusF3iDzF4i.real, AverageTraceF4iDzF3iPlusF3iDzF4i.imag);
	    } /* end of if-condition for trace of fmunu correlators */
	  
	  int FolderNumberIndex = (iters-1)/1000;
          FolderNumber = 1000*(1 + FolderNumberIndex);
	
	  if( SaveLattice==1  )
	    {sprintf(FolderName,"mkdir %s/Nt%d_Ns%d/Beta%.4f_%d",SaveLatticeDataFileDIR, nt, nx, beta, FolderNumber);
	      if((FolderNumber-1000+1)==iters  && this_node==1)
                {system(FolderName);
                }
	      int flag=SAVE_SERIAL;  
	      sprintf(SaveLatticeFileName,"%s/Nt%d_Ns%d/Beta%.4f_%d/Lattice_Nt%d_Ns%d_Beta%.4f.configuration.%d", SaveLatticeDataFileDIR, nt, nx, beta, FolderNumber, nt, nx, beta, iters);
	      save_lattice( flag, SaveLatticeFileName, stringLFN );
	      //rephase( OFF );
	      // save_lattice( saveflag, savefile, stringLFN );
	      //rephase( ON );
	    }	      
	  
	  if(traj_done < trajecs - 1)
	    {
	      if(UseSavedConfiguration==0)
		{ 
		  printf(" Amit MyFFApp/control.c s_iters=update() called for beta= %.4f, at iters = %d \n",beta, iters);
		  s_iters=update();
		}
	      else 
		{int flag=RELOAD_SERIAL;
		  sprintf(ReadLatticeFileName,"%s/Nt%d_Ns%d/Beta%.4f_%d/Lattice_Nt%d_Ns%d_Beta%.4f.configuration.%d", ReadLatticeDataFileDIR, nt, nx, beta, FolderNumber, nt, nx, beta, iters);
		  if(is_file_exist(ReadLatticeFileName)==0) 
		    {
		      while(is_file_exist(ReadLatticeFileName)==0 && traj_done < trajecs - 1)
			{
			  iters++; traj_done++;
			  FolderNumberIndex = (iters-1)/1000;
			  FolderNumber = 1000*(1 + FolderNumberIndex);
			  sprintf(ReadLatticeFileName,"%s/Nt%d_Ns%d/Beta%.4f_%d/Lattice_Nt%d_Ns%d_Beta%.4f.configuration.%d", ReadLatticeDataFileDIR, nt, nx, beta, FolderNumber, nt, nx, beta, iters);
			}
		    }
		  reload_lattice( flag, ReadLatticeFileName);
		}
	    }
	}/* end loop over trajectories */
      

      printf("default MyFFApp/control.c RUNNING COMPLETED, This node is %d, exec=%s \n",this_node, argv[1]); 
      fflush(stdout);	     
      dtime += dclock();	     
      printf("Default MyFFApp/control.c Time = %e seconds \n",dtime);
      printf("Default MyFFApp/control.c total_iters = %d \n",iters);      
      fflush(stdout);
      
      /* save lattice if requested */
      if( saveflag != FORGET )
	{ 
	  save_lattice( saveflag, savefile, stringLFN );
        }
    }
  
  fclose(fploop);
  fclose(ftracefmunuLO);
  fclose(ftracefmunuNLO);
  normal_exit(0);
  return 0;
}
