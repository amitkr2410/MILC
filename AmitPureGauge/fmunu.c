/************************* fmunu.c *******************************/
/* Calculates the temporal and spatial field strength components */
/* and the topological charge.                                   */

/* Includes */
//#include "ks_imp_includes.h"
#include "pure_gauge_includes.h"
#include "../include/field_strength.h"
//#include "../include/generic.h"
#include <assert.h>
/* Computes the real trace of the su3 matrix product: ReTr(A*B) */
Real real_trace_nn( su3_matrix *a, su3_matrix *b )
{
  register int i,j;
  register complex x;
  register Real sum;
  sum = 0.0;
  for( i=0; i<3; i++ )
    {
      for( j=0; j<3; j++ )
	{
	  CMUL( a->e[i][j], b->e[j][i], x );
	  sum += x.real;
	}
    }
  return sum;
}
/* Computes the trace of the su3 matrix product: Tr(A*B) */
complex trace_nn( su3_matrix *a, su3_matrix *b )
{
  register int i,j;
  register complex x;
  register complex sum;
  sum.real = 0.0; sum.imag=0.0;
  for( i=0; i<3; i++ )
    {
      for( j=0; j<3; j++ )
	{
	  CMUL( a->e[i][j], b->e[j][i], x );
	  CADD(x,sum,sum);
	}
    }
  return sum;
}
/* Ans= tr(a*b)*c  , c=real */
complex trace_nn_scalar_mult( su3_matrix *a, su3_matrix *b, double c)
{
  register int i,j;
  register complex x;
  register complex sum;
  sum.real = 0.0; sum.imag=0.0;
  for( i=0; i<3; i++ )
    {
      for( j=0; j<3; j++ )
	{
          CMUL( a->e[i][j], b->e[j][i], x );
          CADD(x,sum,sum);
        }
    }
  CMULREAL(sum,c,sum);// (a,b,c); c = ba with b real and a complex
  return sum;
}

void DisplayNNMatrix(su3_matrix *a )
{
  register int i,j;
  for( i=0; i<3; i++ )
    {
      for( j=0; j<3; j++ )
        {
          printf("(%e,%e) \t",a->e[i][j].real, a->e[i][j].imag);
        }
      printf("\n");
    }
}

void Display(complex a)
{
  printf("(%e,%e)", a.real, a.imag);
}


/* Computes the field strength components and topological charge */
//void fmunu_fmunu(double *time, double *space, double *charge);
//void fmunu_fmunu(complex *TraceF3iF3iMinusF4iF4i, complex *TraceF4iF3iPlusF3iF4i,  complex* TraceF3iD2zF3iMinusF4iD2zF4i,complex *TraceF4iD2zF3iPlusF3iD2zF4i, complex* TraceF3iD4zF3iMinusF4iD4zF4i, complex *TraceF4iD4zF3iPlusF3iD4zF4i)
void fmunu_fmunu(complex *SymmetricTadpole0, complex *SymmetricTadpole2, complex *SymmetricTadpole4, complex *AntiSymmetricTadpole0, complex *AntiSymmetricTadpole2, complex *AntiSymmetricTadpole4)
{
  printf("Start of fmunu_fmunu function");
  
  // Site variables 
  register int i, kplus;
  register site *s;

  // Temporary component storage 
  //su3_matrix *ft, *fs;
  // Initialize sums 
  //double *time, *space, *charge;
  //*time = *space = *charge = 0;
  *SymmetricTadpole0=cmplx(0.0,0.0);
  *SymmetricTadpole2=cmplx(0.0,0.0);
  *SymmetricTadpole4=cmplx(0.0,0.0);
  *AntiSymmetricTadpole0=cmplx(0.0,0.0);
  *AntiSymmetricTadpole2=cmplx(0.0,0.0);
  *AntiSymmetricTadpole4=cmplx(0.0,0.0);  
  // Compute 8*F_mu,nu at each site 
  make_field_strength( F_OFFSET(link), F_OFFSET(fieldstrength) );

  ////Example for phivector  at neighbouring site ( here site is displaced in ZUP direction )
  //msg_tag *tag; 
  //su3_vector phi;
  //tag = start_gather_field( phi, sizeof(su3_vector), ZUP, EVENANDODD, gen_pt[0] );
  // or
  //tag = start_gather_site( F_OFFSET(phi), sizeof(su3_vector), ZUP, EVENANDODD, gen_pt[0] );  
  //wait_gather(tag);
  // //Now gen_pt[0][i] now contains the address of the  phi vector (or a copy therof) on the neighbor of site i in the ZUP direction        for all sites i.
  //cleanup_gather(tag);

  // //Fmunu at neighbouring site
  register su3_matrix *Dummy1, *Dummy2, *Dummy3;
  su3_matrix *LinkZ_zMINUSa, *LinkZ_zMINUS2a;
  su3_matrix *LinkZ_zPLUSa;
  Dummy1 = (su3_matrix*)malloc(sizeof(su3_matrix));  assert( (Dummy1!=NULL) );
  memset ( Dummy1, 0,sizeof(su3_matrix) );
  Dummy2 = (su3_matrix*)malloc(sizeof(su3_matrix));  assert( (Dummy2!=NULL) );
  memset ( Dummy2, 0, sizeof(su3_matrix) );
  Dummy3 = (su3_matrix*)malloc(sizeof(su3_matrix));  assert( (Dummy3!=NULL) );
  memset ( Dummy3, 0, sizeof(su3_matrix) );
      
  LinkZ_zMINUSa = (su3_matrix*)malloc(sites_on_node*sizeof(su3_matrix));  assert( (LinkZ_zMINUSa!=NULL) );
  memset ( LinkZ_zMINUSa, 0, sites_on_node*sizeof(su3_matrix) );
  LinkZ_zMINUS2a = (su3_matrix*)malloc(sites_on_node*sizeof(su3_matrix));  assert( (LinkZ_zMINUS2a!=NULL) );
  memset ( LinkZ_zMINUS2a, 0, sites_on_node*sizeof(su3_matrix) );
  LinkZ_zPLUSa = (su3_matrix*)malloc(sites_on_node*sizeof(su3_matrix));  assert( (LinkZ_zPLUSa!=NULL) );
  memset ( LinkZ_zPLUSa, 0, sites_on_node*sizeof(su3_matrix) );

  msg_tag *tagF[8];//F31Plus, F32Plus, F41Plus, F42Plus, F31Minus, F32Minus, F41Minus, F42Minus;
  msg_tag *tagLinkZ_zMINUSa, *tagLinkZ_zMINUS2a, *tagLinkZ_zPLUSa; //Gather Lint at site z-a, z-2a, z+a
  //Gather Link at site=z-a in +z direction
  tagLinkZ_zMINUSa = start_gather_site( F_OFFSET(link[ZUP]), sizeof(su3_matrix), ZDOWN, EVENANDODD, gen_pt[0]);
  tagLinkZ_zPLUSa  = start_gather_site( F_OFFSET(link[ZUP]), sizeof(su3_matrix), ZUP, EVENANDODD, gen_pt[1]);
  wait_gather(tagLinkZ_zMINUSa);
  wait_gather(tagLinkZ_zPLUSa);
  FORALLSITES(i, s)
    {
      //su3_adjoint( (su3_matrix *)(gen_pt[0][i]), &(LinkZMinus[i]) ); 
      su3mat_copy( (su3_matrix *)(gen_pt[0][i]), &(LinkZ_zMINUSa[i])) ;      
      su3mat_copy( (su3_matrix *)(gen_pt[1][i]), &(LinkZ_zPLUSa[i])) ;
    }
  cleanup_gather(tagLinkZ_zMINUSa); 
  cleanup_gather(tagLinkZ_zPLUSa);

  tagLinkZ_zMINUS2a = start_gather_field( LinkZ_zMINUSa, sizeof(su3_matrix), ZDOWN, EVENANDODD, gen_pt[0]);
  wait_gather(tagLinkZ_zMINUS2a);
  FORALLSITES(i, s)
    {
      su3mat_copy( (su3_matrix *)(gen_pt[0][i]), &(LinkZ_zMINUS2a[i])) ;
    }
  cleanup_gather(tagLinkZ_zMINUS2a); 
  //tagF31Plus, tagF32Plus, tagF41Plus, tagF42Plus;    tagF31Minus, tagF32Minus, tagF41Minus, tagF42Minus
  //ZUP,        ZUP,        ZUP,        ZUP;           ZDOWN,       ZDOWN,       ZDOWN,       ZDOWN
  //gen_pt[0],  gen_pt[1],  gen_pt[2],  gen_pt[3];     gen_pt[4],   gen_pt[5],   gen_pt[6],   gen_pt[7]
  //tagF42Minus   = start_gather_site( F_OFFSET(fieldstrength[FS_YT]), sizeof(su3_matrix), ZDOWN, EVENANDODD, gen_pt[7] );

  for (int k=0;k<2;k++) 
    {
      tagF[0+4*k]   = start_gather_site( F_OFFSET(fieldstrength[FS_XZ]), sizeof(su3_matrix), ( k<1?ZUP:ZDOWN ), EVENANDODD, gen_pt[0+4*k] );
      tagF[1+4*k]   = start_gather_site( F_OFFSET(fieldstrength[FS_YZ]), sizeof(su3_matrix), ( k<1?ZUP:ZDOWN ), EVENANDODD, gen_pt[1+4*k] );
      tagF[2+4*k]   = start_gather_site( F_OFFSET(fieldstrength[FS_XT]), sizeof(su3_matrix), ( k<1?ZUP:ZDOWN ), EVENANDODD, gen_pt[2+4*k] );
      tagF[3+4*k]   = start_gather_site( F_OFFSET(fieldstrength[FS_YT]), sizeof(su3_matrix), ( k<1?ZUP:ZDOWN ), EVENANDODD, gen_pt[3+4*k] );
    }

  //for (int k=0;k<8;k++) { wait_gather(tagF[k]); } 

  su3_matrix ** Buffer=NULL;

  Buffer = (su3_matrix**)malloc(8*sizeof(su3_matrix*));  assert( (Buffer!=NULL) );
  for (int k=0;k<8;k++) 
    {
      Buffer[k] = (su3_matrix*)malloc(sites_on_node*sizeof(su3_matrix));  assert( (Buffer[k]!=NULL) );      
      memset ( Buffer[k], 0, sites_on_node*sizeof(su3_matrix) );
      wait_gather(tagF[k]);
    }

  for (int k=0;k<8;k++) { FORALLSITES(i,s) { su3mat_copy( (su3_matrix *)(gen_pt[k][i]), &(Buffer[k][i]) ); } }

  for (int k=0;k<8;k++) { cleanup_gather(tagF[k]); } 

  for (int k=0;k<8;k++) { tagF[k]  = start_gather_field( Buffer[k], sizeof(su3_matrix), ( k<4?ZUP:ZDOWN ), EVENANDODD, gen_pt[k] ); }
  for (int k=0;k<8;k++) { wait_gather(tagF[k]); } 

  // Loop over each site to sum F_mu,nu components 
  FORALLSITES(i, s) 
    {             
      //Fmunu[k] = s->fieldstrength[k+1]; //kplus=k+1=1,2,3,4 ; F31, F32, F41, F42
      //LinkZPlus  = &(s->link[ZUP]);      
      
      for(int k=0; k<8; k++)
	{
	  kplus = (k%4) +1;  //&(s->fieldstrength[kplus]), kplus=k+1=1,2,3,4 ; F31, F32, F41, F42
          
	  if(k<4) {   mult_su3_nn( &(s->link[ZUP]), &(LinkZ_zPLUSa[i]), Dummy1);}
	  else    
	    {  mult_su3_nn( &(LinkZ_zMINUS2a[i]), &(LinkZ_zMINUSa[i]), Dummy2);
	      su3_adjoint(Dummy2,Dummy1);
	    }
	  
	  mult_su3_nn(Dummy1, (su3_matrix *)(gen_pt[k][i]), Dummy2);
	  mult_su3_na(Dummy2, Dummy1, Dummy3);
	  CADD(*SymmetricTadpole4,  trace_nn_scalar_mult(&(s->fieldstrength[kplus]),Dummy3, kplus<3?1:-1), *SymmetricTadpole4);
	  CADD(*AntiSymmetricTadpole4,  trace_nn(&(s->fieldstrength[(((k+2)%4)+1)]),Dummy3), *AntiSymmetricTadpole4);
	  
	  if(k<4) {   su3mat_copy(&(s->link[ZUP]), Dummy1);}
          else    {   su3_adjoint( &(LinkZ_zMINUSa[i]), Dummy1);}
	  
	  mult_su3_nn(Dummy1, &(Buffer[k][i]), Dummy2);
	  mult_su3_na(Dummy2, Dummy1, Dummy3);
	  CADD(*SymmetricTadpole2,  trace_nn_scalar_mult(&(s->fieldstrength[kplus]),Dummy3, kplus<3?1:-1), *SymmetricTadpole2);
          CADD(*AntiSymmetricTadpole2,  trace_nn(&(s->fieldstrength[(((k+2)%4)+1)]),Dummy3), *AntiSymmetricTadpole2);
	  
	  if(k<4)
	    {
	      CADD(*SymmetricTadpole0,  trace_nn_scalar_mult(&(s->fieldstrength[kplus]), &(s->fieldstrength[kplus]), kplus<3?1:-1), *SymmetricTadpole0);
	      CADD(*AntiSymmetricTadpole0,  trace_nn(&(s->fieldstrength[(((k+2)%4)+1)]), &(s->fieldstrength[kplus]) ), *AntiSymmetricTadpole0);
	    }
	}
            
      //printf("\t TraceF4iF3iPlusF3iF4i = ");Display(*TraceF4iF3iPlusF3iF4i);
    }

  g_complexsum(SymmetricTadpole0);
  g_complexsum(SymmetricTadpole2);
  g_complexsum(SymmetricTadpole4);
  g_complexsum(AntiSymmetricTadpole0);
  g_complexsum(AntiSymmetricTadpole2);
  g_complexsum(AntiSymmetricTadpole4);
  
  
  
  for (int k=0;k<8;k++) { cleanup_gather(tagF[k]); }
  for (int k=0;k<8;k++) { free(Buffer[k]); Buffer[k]=NULL; }
  free(Buffer); Buffer=NULL;

  free(LinkZ_zMINUSa);  free(LinkZ_zMINUS2a); free(LinkZ_zPLUSa);
  free(Dummy1);  free(Dummy2);  free(Dummy3); 
  
  // Sum over all nodes 
  //g_doublesum(time);
  //g_doublesum(space);
  //g_doublesum(charge);

  // Normalizations 
  //CMULREAL(a,b,c)   c = ba with b real and a complex
  double Normalization=-1.0/(4.0*volume*16.0); //Minus-sign bc of 2i^2 in denominator of Fmunu, factor 16 for clover FF
  CMULREAL(*SymmetricTadpole0, Normalization, *SymmetricTadpole0);  
  CMULREAL(*SymmetricTadpole2, Normalization, *SymmetricTadpole2);
  CMULREAL(*SymmetricTadpole4, Normalization, *SymmetricTadpole4);
  CMULREAL(*AntiSymmetricTadpole0, Normalization, *AntiSymmetricTadpole0);
  CMULREAL(*AntiSymmetricTadpole2, Normalization, *AntiSymmetricTadpole2);
  CMULREAL(*AntiSymmetricTadpole4, Normalization, *AntiSymmetricTadpole4);
  

  printf("\n Total TraceF3iF3iMinusF4iF4i = "); Display(*SymmetricTadpole0);
  printf("\n Total TraceF4iF3iPlusF3iF4i = ");  Display(*AntiSymmetricTadpole0);
  printf("\n Total SymmetricTadpole2 = ");      Display(*SymmetricTadpole2);
  printf("\n Total AntiSymmetricTadpole2 = ");  Display(*AntiSymmetricTadpole2);
  printf("\n Total SymmetricTadpole4 = ");      Display(*SymmetricTadpole4);
  printf("\n Total AntiSymmetricTadpole4 = ");  Display(*AntiSymmetricTadpole4);
  printf("\n");
}


