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

  //Plaquette and Field-strength variable
  int Nc=3;//Color factor
  double SS_Plaq=0.0, ST_Plaq=0.0;
  double Current_Plaq=0.0, Sum_Plaq=0.0, Average_Plaq=0.0;
  double TadpoleFactor=0.0, TadpoleFactorNt=0.0; 
  complex CurrentPolyakovLoop, SumPolyakovLoop=cmplx(0.0,0.0), AveragePolyakovLoop;
  double CurrentModPolyakovLoop=0.0, SumModPolyakovLoop=0.0, AverageModPolyakovLoop=0.0;
  double SumModPolyakovLoopTadpoleCorrected=0.0, AverageModPolyakovLoopTadpoleCorrected=0.0;
  double CurrentBareFreeEnergy=0.0,  SumBareFreeEnergy=0.0, AverageBareFreeEnergy=0.0;
  double CurrentBareFreeEnergyTadpoleCorrected=0.0,  SumBareFreeEnergyTadpoleCorrected=0.0, AverageBareFreeEnergyTadpoleCorrected=0.0;

  complex CurrentSymmetricTadpole0, SumSymmetricTadpole0=cmplx(0.0,0.0), AverageSymmetricTadpole0;
  complex CurrentAntiSymmetricTadpole0, SumAntiSymmetricTadpole0=cmplx(0.0,0.0), AverageAntiSymmetricTadpole0;
  complex CurrentSymmetricTadpole2, SumSymmetricTadpole2=cmplx(0.0,0.0), AverageSymmetricTadpole2;
  complex CurrentAntiSymmetricTadpole2, SumAntiSymmetricTadpole2=cmplx(0.0,0.0), AverageAntiSymmetricTadpole2;
  complex CurrentSymmetricTadpole4, SumSymmetricTadpole4=cmplx(0.0,0.0), AverageSymmetricTadpole4;
  complex CurrentAntiSymmetricTadpole4, SumAntiSymmetricTadpole4=cmplx(0.0,0.0), AverageAntiSymmetricTadpole4;

  //Initialize variable to zero

  //FileName to save observables
  FILE *fploop, *ftracefmunuLO, *ftracefmunuNLO, *ftracefmunuNNLO;
  char FileNamePloop[10000], FileNameTraceFmunu[10000], FileNameTraceFmunu2[10000], FileNameTraceFmunu3[10000], FileNameTraceFmunu4[10000], SaveLatticeFileName[100000], ReadLatticeFileName[100000], FolderName[10000];

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
      sprintf(FileNameTraceFmunu,"%s/DataTraceFmunuTadpole0_Clover_Traceless_Nt%d_Ns%d_Beta%.4f.txt", OutputDataFileDIR, nt, nx, beta);
      sprintf(FileNameTraceFmunu2,"%s/DataTraceFmunuTadpole2_Clover_Traceless_Nt%d_Ns%d_Beta%.4f.txt", OutputDataFileDIR, nt, nx, beta);
      sprintf(FileNameTraceFmunu3,"%s/DataTraceFmunuTadpole4_Clover_Traceless_Nt%d_Ns%d_Beta%.4f.txt", OutputDataFileDIR, nt, nx, beta);

      fploop = fopen(FileNamePloop,"w");
      ftracefmunuLO  = fopen(FileNameTraceFmunu,"w");
      ftracefmunuNLO = fopen(FileNameTraceFmunu2,"w");
      ftracefmunuNNLO = fopen(FileNameTraceFmunu3,"w");

      fprintf(fploop,"#Beta=%.4f, Nt=%d, Ns=%d^3 \n", beta, nt, nx);
      fprintf(fploop,"#Iters \t Current_Plaq \t AvgPlaq \t TadpoleFactor \t TadpoleFactorNt \t CurrentPolyakovLoop.real \t CurrentPolyakovLoop.imag \t  CurrentModPolyakovLoop \t AverageModPolyakovLoop \t AverageModPolyakovLoopTadpoleCorrected \t  CurrentBareFreeEnergy \t AvgBareFreeEnergy  \t CurrentBareFreeEnergyTadpoleCorrected \t AvgBareFreeEnergyTadpoleCorrected \n");
      fprintf(ftracefmunuLO,"#Beta=%.4f, Nt=%d, Ns=%d^3 \n", beta, nt, nx);
      fprintf(ftracefmunuLO,"#Iters \t SymmetricTadpole0.real \t SymmetricTadpole0.imag \t AvgSymmetricTadpole0.real \t AvgSymmetricTadpole0.imag \t AntiSymmetricTadpole0.real \t AntiSymmetricTadpole0.imag \t AvgAntiSymmetricTadpole0.real \t AvgAntiSymmetricTadpole0.imag \n");
      fprintf(ftracefmunuNLO,"#Beta=%.4f, Nt=%d, Ns=%d^3 \n", beta, nt, nx);
      fprintf(ftracefmunuNLO,"#Iters \t SymmetricTadpole2.real \t SymmetricTadpole2.imag \t AvgSymmetricTadpole2.real \t AvgSymmetricTadpole2.imag \t AntiSymmetricTadpole2.real \t AntiSymmetricTadpole2.imag \t AvgAntiSymmetricTadpole2.real \t AvgAntiSymmetricTadpole2.imag \n");
      fprintf(ftracefmunuNNLO,"#Beta=%.4f, Nt=%d, Ns=%d^3 \n", beta, nt, nx);
      fprintf(ftracefmunuNNLO,"#Iters \t SymmetricTadpole4.real \t SymmetricTadpole4.imag \t AvgSymmetricTadpole4.real \t AvgSymmetricTadpole4.imag \t AntiSymmetricTadpole4.real \t AntiSymmetricTadpole4.imag \t AvgAntiSymmetricTadpole4.real \t AvgAntiSymmetricTadpole4.imag \n");


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
      MeasurementCount=0;		/* number of measurements 		*/            
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
	      SumSymmetricTadpole0 = cmplx(0.0,0.0);  SumAntiSymmetricTadpole0 = cmplx(0.0,0.0);
	      SumSymmetricTadpole2 = cmplx(0.0,0.0);  SumAntiSymmetricTadpole2 = cmplx(0.0,0.0);
	      SumSymmetricTadpole4 = cmplx(0.0,0.0);  SumAntiSymmetricTadpole4 = cmplx(0.0,0.0);
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
	      CurrentSymmetricTadpole0 = cmplx(0.0,0.0);  CurrentAntiSymmetricTadpole0 = cmplx(0.0,0.0);
	      CurrentSymmetricTadpole2 = cmplx(0.0,0.0);  CurrentAntiSymmetricTadpole2 = cmplx(0.0,0.0);
	      CurrentSymmetricTadpole4 = cmplx(0.0,0.0);  CurrentAntiSymmetricTadpole4 = cmplx(0.0,0.0);

              fmunu_fmunu(&CurrentSymmetricTadpole0, &CurrentSymmetricTadpole2, &CurrentSymmetricTadpole4, &CurrentAntiSymmetricTadpole0, &CurrentAntiSymmetricTadpole2, &CurrentAntiSymmetricTadpole4);

              CADD(SumSymmetricTadpole0, CurrentSymmetricTadpole0, SumSymmetricTadpole0);
	      CADD(SumSymmetricTadpole2, CurrentSymmetricTadpole2, SumSymmetricTadpole2);
	      CADD(SumSymmetricTadpole4, CurrentSymmetricTadpole4, SumSymmetricTadpole4);
	      CADD(SumAntiSymmetricTadpole0, CurrentAntiSymmetricTadpole0, SumAntiSymmetricTadpole0);
	      CADD(SumAntiSymmetricTadpole2, CurrentAntiSymmetricTadpole2, SumAntiSymmetricTadpole2);
	      CADD(SumAntiSymmetricTadpole4, CurrentAntiSymmetricTadpole4, SumAntiSymmetricTadpole4);

              CDIVREAL(SumSymmetricTadpole0,     MeasurementCount, AverageSymmetricTadpole0);
	      CDIVREAL(SumSymmetricTadpole2,     MeasurementCount, AverageSymmetricTadpole2);
	      CDIVREAL(SumSymmetricTadpole4,     MeasurementCount, AverageSymmetricTadpole4);
	      CDIVREAL(SumAntiSymmetricTadpole0,     MeasurementCount, AverageAntiSymmetricTadpole0);
	      CDIVREAL(SumAntiSymmetricTadpole2,     MeasurementCount, AverageAntiSymmetricTadpole2);
	      CDIVREAL(SumAntiSymmetricTadpole4,     MeasurementCount, AverageAntiSymmetricTadpole4);


              printf("Amit MyFFApp/control.c TraceF3iF3iMinusF4iF4i=(%e,%e), AvgTrace=(%e,%e) \n", CurrentSymmetricTadpole0.real, CurrentSymmetricTadpole0.imag, AverageSymmetricTadpole0.real, AverageSymmetricTadpole0.imag);
              printf("Amit MyFFApp/control.c TraceF4iF3iPlusF3iF4i=(%e,%e), AvgTrace=(%e,%e) \n", CurrentAntiSymmetricTadpole0.real, CurrentAntiSymmetricTadpole0.imag, AverageAntiSymmetricTadpole0.real, AverageAntiSymmetricTadpole0.imag);

	      printf("Amit MyFFApp/control.c TraceSymmetricTadpole2=(%e,%e), AvgTrace=(%e,%e) \n", CurrentSymmetricTadpole2.real, CurrentSymmetricTadpole2.imag, AverageSymmetricTadpole2.real, AverageSymmetricTadpole2.imag);
              printf("Amit MyFFApp/control.c TraceAntiSymmetricTadpole2=(%e,%e), AvgTrace=(%e,%e) \n", CurrentAntiSymmetricTadpole2.real, CurrentAntiSymmetricTadpole2.imag, AverageAntiSymmetricTadpole2.real, AverageAntiSymmetricTadpole2.imag);

	      printf("Amit MyFFApp/control.c TraceSymmetricTadpole4=(%e,%e), AvgTrace=(%e,%e) \n", CurrentSymmetricTadpole4.real, CurrentSymmetricTadpole4.imag, AverageSymmetricTadpole4.real, AverageSymmetricTadpole4.imag);
              printf("Amit MyFFApp/control.c TraceAntiSymmetricTadpole4=(%e,%e), AvgTrace=(%e,%e) \n", CurrentAntiSymmetricTadpole4.real, CurrentAntiSymmetricTadpole4.imag, AverageAntiSymmetricTadpole4.real, AverageAntiSymmetricTadpole4.imag);

	      fprintf(ftracefmunuLO,"%d \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n", iters, CurrentSymmetricTadpole0.real, CurrentSymmetricTadpole0.imag, AverageSymmetricTadpole0.real, AverageSymmetricTadpole0.imag, CurrentAntiSymmetricTadpole0.real, CurrentAntiSymmetricTadpole0.imag, AverageAntiSymmetricTadpole0.real, AverageAntiSymmetricTadpole0.imag );
	      fprintf(ftracefmunuNLO,"%d \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n", iters, CurrentSymmetricTadpole2.real, CurrentSymmetricTadpole2.imag, AverageSymmetricTadpole2.real, AverageSymmetricTadpole2.imag, CurrentAntiSymmetricTadpole2.real, CurrentAntiSymmetricTadpole2.imag, AverageAntiSymmetricTadpole2.real, AverageAntiSymmetricTadpole2.imag );
	      fprintf(ftracefmunuNNLO,"%d \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n", iters, CurrentSymmetricTadpole4.real, CurrentSymmetricTadpole4.imag, AverageSymmetricTadpole4.real, AverageSymmetricTadpole4.imag, CurrentAntiSymmetricTadpole4.real, CurrentAntiSymmetricTadpole4.imag, AverageAntiSymmetricTadpole4.real, AverageAntiSymmetricTadpole4.imag );	      
	    }
	

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
  fclose(ftracefmunuNNLO);
  normal_exit(0);
  return 0;
}
