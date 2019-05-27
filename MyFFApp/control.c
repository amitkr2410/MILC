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

EXTERN gauge_header start_lat_hdr;	/* Input gauge field header */

int main( int argc, char **argv )
{
  int i,MeasurementCount,traj_done, naik_index;
  int prompt;
  int s_iters=0, iters=0;
  double dtime, dclock();

  //Plaquette and Field-strength variable
  int Nc=3;//Color factor
  double SS_Plaq=0.0, ST_Plaq=0.0;
  double Current_Plaq=0.0, Sum_Plaq=0.0, Average_Plaq=0.0;
  double TadpoleFactor=0.0, TadpoleFactorNt=0.0; 
  complex CurrentPolyakovLoop, SumPolyakovLoop, AveragePolyakovLoop;
  double CurrentModPolyakovLoop=0.0, SumModPolyakovLoopTadpoleCorrected=0.0, AverageModPolyakovLoopTadpoleCorrected=0.0;
  double CurrentBareFreeEnergy=0.0,  SumBareFreeEnergy=0.0, AverageBareFreeEnergy=0.0;
  double CurrentBareFreeEnergyTadpoleCorrected=0.0,  SumBareFreeEnergyTadpoleCorrected=0.0, AverageBareFreeEnergyTadpoleCorrected=0.0;
  
  //Initialize variable to zero
  CurrentPolyakovLoop=cmplx(0.0,0.0); SumPolyakovLoop=cmplx(0.0,0.0); AveragePolyakovLoop=cmplx(0.0,0.0);

  //FileName to save observables
  FILE *fploop;
  char FileNamePloop[1000];

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
      sprintf(FileNamePloop,"Output/DataPloopNt%d_Beta%.2f_ml%.4f_ms%.4f.txt", nt, beta, dyn_mass[0], dyn_mass[1]);
      fploop = fopen(FileNamePloop,"w");

      fprintf(fploop,"#Beta=%.4f, ml=%.6f, ms=%.6f,  Nt=%d, Ns=%d^3 \n", beta, dyn_mass[0],dyn_mass[1], nt, nx);
      fprintf(fploop,"#Iters \t Current_Plaq \t AvgPlaq \t TadpoleFactor \t TadpoleFactorNt \t CurrentPolyakovLoop.real \t CurrentPolyakovLoop.imag \t  CurrentModPolyakovLoop \t AverageModPolyakovLoopTadpoleCorrected \t  CurrentBareFreeEnergy \t AvgBareFreeEnergy  \t CurrentBareFreeEnergyTadpoleCorrected \t AvgBareFreeEnergyTadpoleCorrected \n");

      /* perform warmup trajectories */
      #ifdef MILC_GLOBAL_DEBUG
      global_current_time_step = 0;
      #endif /* MILC_GLOBAL_DEBUG */
    
      // Start application timer
      dtime = -dclock();
      printf(" Amit MyFFApp/control.c inside while(readin(prompt)==0) \n");
      for( traj_done=0; traj_done < warms; traj_done++ )
      	{
	  rephase(OFF);
	  SS_Plaq=0.0; ST_Plaq=0.0;
	  d_plaquette(&SS_Plaq, &ST_Plaq);
	  printf("Amit MyFFApp/control.c Plaquette = (%e,%e)\n",SS_Plaq, ST_Plaq);	  
	  rephase(ON);
      	  update();	  
        }
      
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

	  rephase(OFF);
          printf(" Amit MyFFApp/control.c s_iters=update() called at iters = %d \n", iters);
	  /* measure every "propinterval" trajectories */
	  //	 if( (traj_done%propinterval)==(propinterval-1) )
	  //{
	      iters = 1 + iters;
	      MeasurementCount = MeasurementCount + 1;
	      if(iters ==200) 
		{
		  MeasurementCount = 1; 
		  Sum_Plaq=0.0;
		  SumPolyakovLoop = cmplx(0.0,0.0);
		  SumModPolyakovLoopTadpoleCorrected=0.0;
		  SumBareFreeEnergy=0.0;
		  SumBareFreeEnergyTadpoleCorrected=0.0;
		}
	      /* call gauge_variable fermion_variable measuring routines */
	      //rephase(OFF);	      
	      /* Compute plaquette and display output */
	      SS_Plaq=0.0; ST_Plaq=0.0;
	      d_plaquette(&SS_Plaq, &ST_Plaq);
	      Current_Plaq = ((SS_Plaq+ST_Plaq)/2.0)/Nc;
	      Sum_Plaq= Sum_Plaq + Current_Plaq;
	      Average_Plaq = Sum_Plaq/MeasurementCount;
	      TadpoleFactor = pow(Current_Plaq, 1.0/4.0);
	      TadpoleFactorNt = pow(Current_Plaq, nt/4.0);
	      printf("Amit MyFFApp/control.c Plaquette=(%e,%e), CurrentPlaq=%e, AvgPlaq=%e \n",SS_Plaq, ST_Plaq, Current_Plaq, Average_Plaq);

	      /* Calculate trace of polyakov loop */
	      CurrentPolyakovLoop=cmplx(0.0,0.0); 
	      CurrentPolyakovLoop = ploop();

	      CDIVREAL(CurrentPolyakovLoop,3.0,CurrentPolyakovLoop);								   
	      CADD(SumPolyakovLoop, CurrentPolyakovLoop, SumPolyakovLoop);
	      CDIVREAL(SumPolyakovLoop, MeasurementCount, AveragePolyakovLoop);

	      complex *PointerPLoop = &CurrentPolyakovLoop;
	      CurrentModPolyakovLoop=cabs(PointerPLoop);
	      SumModPolyakovLoopTadpoleCorrected = SumModPolyakovLoopTadpoleCorrected + (CurrentModPolyakovLoop*TadpoleFactorNt);
	      AverageModPolyakovLoopTadpoleCorrected = SumModPolyakovLoopTadpoleCorrected/MeasurementCount;
	      
	      CurrentBareFreeEnergy = -log(CurrentModPolyakovLoop);
              SumBareFreeEnergy = SumBareFreeEnergy + CurrentBareFreeEnergy;
              AverageBareFreeEnergy = SumBareFreeEnergy/MeasurementCount;

	      CurrentBareFreeEnergyTadpoleCorrected = -log(TadpoleFactorNt*CurrentModPolyakovLoop);
	      SumBareFreeEnergyTadpoleCorrected = SumBareFreeEnergyTadpoleCorrected + CurrentBareFreeEnergyTadpoleCorrected;
	      AverageBareFreeEnergyTadpoleCorrected = SumBareFreeEnergyTadpoleCorrected/MeasurementCount;

	      printf("Amit MyFFApp/control.c PLoop=(%e,%e), AvgPLoop=(%e,%e), and (CurrentModPLOOP,AvgModPLOOPTadpoleCorr)=(%e,%e), and (CurrentBareFreeEnergy, AverageBareFreeEnergy)=(%e,%e), and (CurrenttBareFreeEnergyTadpoleCorr, AverageBareFreeEnergyTadpoleCorr)=(%e,%e)\n", CurrentPolyakovLoop.real, CurrentPolyakovLoop.imag, AveragePolyakovLoop.real, AveragePolyakovLoop.imag, CurrentModPolyakovLoop, AverageModPolyakovLoopTadpoleCorrected, CurrentBareFreeEnergy, AverageBareFreeEnergy, CurrentBareFreeEnergyTadpoleCorrected, AverageBareFreeEnergyTadpoleCorrected);

	      /* Calculate trace of fmunu and output */	      
	      fprintf(fploop,"%d \t %e \t %.4f  \t %.4f \t %.4f \t %.4f \t %.4f \t %e \t %.4f \t %e  \t %.4f  \t %e \t %.4f \n", iters, Current_Plaq,  Average_Plaq, TadpoleFactor, TadpoleFactorNt, CurrentPolyakovLoop.real, CurrentPolyakovLoop.imag, CurrentModPolyakovLoop, AverageModPolyakovLoopTadpoleCorrected,   CurrentBareFreeEnergy,  AverageBareFreeEnergy, CurrentBareFreeEnergyTadpoleCorrected, AverageBareFreeEnergyTadpoleCorrected);
	      
	 rephase(ON);
	 /* Compute chiral condensate pbp, etc */
	 /* Make fermion links if not already done */
	 restore_fermion_links_from_site(fn_links, par_buf.prec_pbp);
	 for(i = 0; i < par_buf.num_pbp_masses; i++)
	   {
             #if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
	     naik_index = par_buf.ksp_pbp[i].naik_term_epsilon_index;
             #else
	     naik_index = 0;
             #endif
	     f_meas_imp_field( par_buf.npbp_reps, &par_buf.qic_pbp[i], par_buf.ksp_pbp[i].mass, naik_index, fn_links);
	   }
	 
	 if(traj_done < trajecs - 1)
	   {
	     s_iters=update();
	   }
	}	   	/* end loop over trajectories */
       
      if(this_node==1)
	{
	     printf("\n\n default MyFFApp/control.c RUNNING COMPLETED, This node is %d\n\n",this_node); 
	     fflush(stdout);	     
	}

      dtime += dclock();
      if(this_node==1)
	{
	  printf("\n\n Default MyFFApp/control.c Time = %e seconds\n",dtime);
	  printf("Default MyFFApp/control.c total_iters = %d \n\n",iters);
	}
      fflush(stdout);
      
      /* save lattice if requested */
      if( saveflag != FORGET )
	{
	  rephase( OFF );
	  save_lattice( saveflag, savefile, stringLFN );
	  rephase( ON );
	}
      
      /* Destroy fermion links (created in readin() */
      
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
  //  fclose(ftracefmunu);
  normal_exit(0);
  return 0;
}
