/************************* fmunu.c *******************************/
/* Calculates the temporal and spatial field strength components */
/* and the topological charge.                                   */

/* Includes */
#include "ks_imp_includes.h"
#include "../include/field_strength.h"
#include "../generic/generic_includes.h"
#include <assert.h>
//#define LINK(dir) (((su3_matrix *)F_PT(s,link_src))[dir])
//#include "../include/generic.h"
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
void fmunu_fmunu(complex *TraceF3iF3iMinusF4iF4i, complex *TraceF4iF3iPlusF3iF4i, complex *TraceF3iDzF3iMinusF4iDzF4i, complex *TraceF4iDzF3iPlusF3iDzF4i, complex* TraceF3iD2zF3iMinusF4iD2zF4i,complex *TraceF4iD2zF3iPlusF3iD2zF4i, complex* TraceF3iD4zF3iMinusF4iD4zF4i, complex *TraceF4iD4zF3iPlusF3iD4zF4i)
{
  printf("Start of fmunu_fmunu function");
  
  // Site variables 
  register int i, kplus;
  register site *s;
  int *Displacement;
  Displacement = (int*)malloc(4*sizeof(int));  assert( (Displacement!=NULL) );
  memset (Displacement, 0,4*sizeof(int) );
  Displacement[0]=0; Displacement[1]=0;Displacement[2]=0;Displacement[3]=0; //X,Y,Z,T 
  // Temporary component storage 
  //su3_matrix *ft, *fs;
  // Initialize sums 
  //double *time, *space, *charge;
  //*time = *space = *charge = 0;
  *TraceF3iF3iMinusF4iF4i=cmplx(0.0,0.0);
  *TraceF4iF3iPlusF3iF4i = cmplx(0.0,0.0);
  *TraceF3iDzF3iMinusF4iDzF4i=cmplx(0.0,0.0);
  *TraceF4iDzF3iPlusF3iDzF4i = cmplx(0.0,0.0);
  *TraceF3iD2zF3iMinusF4iD2zF4i=cmplx(0.0,0.0);
  *TraceF4iD2zF3iPlusF3iD2zF4i = cmplx(0.0,0.0);
  *TraceF3iD4zF3iMinusF4iD4zF4i=cmplx(0.0,0.0);
  *TraceF4iD4zF3iPlusF3iD4zF4i = cmplx(0.0,0.0);
  // Compute 8*F_mu,nu at each site 
  make_field_strength( F_OFFSET(link), F_OFFSET(fieldstrength) );

  ////Example for phivector  at neighbouring site ( here site is displaced in ZUP direction )
  //msg_tag *tag; 
  //su3_vector phi;
  //tag = start_gather_field( F_OFFSET(phi), sizeof(su3_vector), ZUP, EVENANDODD, gen_pt[0] );
  // or
  //tag = start_gather_site( F_OFFSET(phi), sizeof(su3_vector), ZUP, EVENANDODD, gen_pt[0] );  
  //wait_gather(tag);
  // //Now gen_pt[0][i] now contains the address of the  phi vector (or a copy therof) on the neighbor of site i in the ZUP direction        for all sites i.
  //cleanup_gather(tag);

  // //Fmunu at neighbouring site
  register su3_matrix *Dummy1, *Dummy2, *Dummy3, *Dummy4;
  su3_matrix *DzFmunu, *D2zFmunu, *D4zFmunu;
  su3_matrix *LinkZ_zMINUSa, *LinkZ_zMINUS2a;
  su3_matrix *LinkZ_zPLUSa;
  Dummy1 = (su3_matrix*)malloc(sizeof(su3_matrix));  assert( (Dummy1!=NULL) );
  memset ( Dummy1, 0,sizeof(su3_matrix) );
  Dummy2 = (su3_matrix*)malloc(sizeof(su3_matrix));  assert( (Dummy2!=NULL) );
  memset ( Dummy2, 0, sizeof(su3_matrix) );
  Dummy3 = (su3_matrix*)malloc(sizeof(su3_matrix));  assert( (Dummy3!=NULL) );
  memset ( Dummy3, 0, sizeof(su3_matrix) );
  Dummy4 = (su3_matrix*)malloc(sizeof(su3_matrix));  assert( (Dummy4!=NULL) );
  memset ( Dummy4, 0, sizeof(su3_matrix) );
  
  DzFmunu = (su3_matrix*)malloc(4*sizeof(su3_matrix));  assert( (DzFmunu!=NULL) );
  memset ( DzFmunu, 0, 4*sizeof(su3_matrix) );
  D2zFmunu = (su3_matrix*)malloc(4*sizeof(su3_matrix));  assert( (D2zFmunu!=NULL) );
  memset ( D2zFmunu, 0, 4*sizeof(su3_matrix) );  
  D4zFmunu = (su3_matrix*)malloc(4*sizeof(su3_matrix));  assert( (D4zFmunu!=NULL) );
  memset ( D4zFmunu, 0, 4*sizeof(su3_matrix) );
    
  LinkZ_zMINUSa = (su3_matrix*)malloc(sites_on_node*sizeof(su3_matrix));  assert( (LinkZ_zMINUSa!=NULL) );
  memset ( LinkZ_zMINUSa, 0, sites_on_node*sizeof(su3_matrix) );
  LinkZ_zMINUS2a = (su3_matrix*)malloc(sites_on_node*sizeof(su3_matrix));  assert( (LinkZ_zMINUS2a!=NULL) );
  memset ( LinkZ_zMINUS2a, 0, sites_on_node*sizeof(su3_matrix) );
  LinkZ_zPLUSa = (su3_matrix*)malloc(sites_on_node*sizeof(su3_matrix));  assert( (LinkZ_zPLUSa!=NULL) );
  memset ( LinkZ_zPLUSa, 0, sites_on_node*sizeof(su3_matrix) );
  msg_tag *tagF31Plus, *tagF32Plus, *tagF41Plus, *tagF42Plus;
  msg_tag *tagF31Minus, *tagF32Minus, *tagF41Minus, *tagF42Minus;
  msg_tag *tagLinkZ_zMINUSa;
  //Gather Link at site=z-a in +z direction
  tagLinkZ_zMINUSa = start_gather_site( F_OFFSET(link[ZUP]), sizeof(su3_matrix), ZDOWN, EVENANDODD, gen_pt[0]);
  wait_gather(tagLinkZ_zMINUSa);
  FORALLSITES(i, s)
    {
      //su3_adjoint( (su3_matrix *)(gen_pt[0][i]), &(LinkZMinus[i]) ); 
      su3mat_copy( (su3_matrix *)(gen_pt[0][i]), &(LinkZ_zMINUSa[i])) ;      
    }
  cleanup_gather(tagLinkZ_zMINUSa);
  
  tagF31Plus   = start_gather_site( F_OFFSET(fieldstrength[FS_XZ]), sizeof(su3_matrix), ZUP, EVENANDODD, gen_pt[0] );
  tagF32Plus   = start_gather_site( F_OFFSET(fieldstrength[FS_YZ]), sizeof(su3_matrix), ZUP, EVENANDODD, gen_pt[1] );
  tagF41Plus   = start_gather_site( F_OFFSET(fieldstrength[FS_XT]), sizeof(su3_matrix), ZUP, EVENANDODD, gen_pt[2] );
  tagF42Plus   = start_gather_site( F_OFFSET(fieldstrength[FS_YT]), sizeof(su3_matrix), ZUP, EVENANDODD, gen_pt[3] );
  
  tagF31Minus   = start_gather_site( F_OFFSET(fieldstrength[FS_XZ]), sizeof(su3_matrix), ZDOWN, EVENANDODD, gen_pt[4] );
  tagF32Minus   = start_gather_site( F_OFFSET(fieldstrength[FS_YZ]), sizeof(su3_matrix), ZDOWN, EVENANDODD, gen_pt[5] );
  tagF41Minus   = start_gather_site( F_OFFSET(fieldstrength[FS_XT]), sizeof(su3_matrix), ZDOWN, EVENANDODD, gen_pt[6] );
  tagF42Minus   = start_gather_site( F_OFFSET(fieldstrength[FS_YT]), sizeof(su3_matrix), ZDOWN, EVENANDODD, gen_pt[7] );

  wait_gather(tagF31Plus);    wait_gather(tagF32Plus);   wait_gather(tagF41Plus);   wait_gather(tagF42Plus);
  wait_gather(tagF31Minus);   wait_gather(tagF32Minus);  wait_gather(tagF41Minus);  wait_gather(tagF42Minus);
  
  // Loop over each site to sum F_mu,nu components 
  FORALLSITES(i, s) 
    {             
      //Fmunu[k] = s->fieldstrength[k+1]; //kplus=k+1=1,2,3,4 ; F31, F32, F41, F42
      //LinkZPlus  = &(s->link[ZUP]);      
      
      for(int k=0; k<4; k++)
	{
	  kplus = k+1;  //&(s->fieldstrength[kplus]), kplus=k+1=1,2,3,4 ; F31, F32, F41, F42
	  mult_su3_nn( &(s->link[ZUP]), (su3_matrix *)(gen_pt[k][i]), Dummy3 );
	  mult_su3_na( Dummy3, &(s->link[ZUP]), Dummy1);
	  mult_su3_an( &(LinkZ_zMINUSa[i]), (su3_matrix *)(gen_pt[k+4][i]), Dummy3);
	  mult_su3_nn( Dummy3, &(LinkZ_zMINUSa[i]), Dummy2 );    
	  sub_su3_matrix(Dummy1, Dummy2, &DzFmunu[k]);
	  add_su3_matrix(Dummy1, Dummy2, &D2zFmunu[k]);

	  add_su3_matrix(Dummy1, Dummy2, &D4zFmunu[k]);
	  scalar_mult_su3_matrix(&D4zFmunu[k], -4.0, &D4zFmunu[k]);
	  scalar_mult_su3_matrix( &(s->fieldstrength[kplus]), 6.0, Dummy3);
	  add_su3_matrix( &D4zFmunu[k], Dummy3, &D4zFmunu[k]);

	  scalar_mult_su3_matrix( &(s->fieldstrength[kplus]), 2.0, Dummy3);
	  sub_su3_matrix( &D2zFmunu[k], Dummy3, &D2zFmunu[k]);
	  
	  
	  if(k<2)
	    {
	      CADD(*TraceF3iF3iMinusF4iF4i, trace_nn(&(s->fieldstrength[kplus]), &(s->fieldstrength[kplus])) , *TraceF3iF3iMinusF4iF4i);
	      CADD(*TraceF3iDzF3iMinusF4iDzF4i, trace_nn(&(s->fieldstrength[kplus]), &DzFmunu[k]), *TraceF3iDzF3iMinusF4iDzF4i);
	      CADD(*TraceF4iDzF3iPlusF3iDzF4i , trace_nn(&(s->fieldstrength[kplus+2]), &DzFmunu[k]), *TraceF4iDzF3iPlusF3iDzF4i);
	      CADD(*TraceF3iD2zF3iMinusF4iD2zF4i, trace_nn(&(s->fieldstrength[kplus]), &D2zFmunu[k]), *TraceF3iD2zF3iMinusF4iD2zF4i);
	      CADD(*TraceF4iD2zF3iPlusF3iD2zF4i , trace_nn(&(s->fieldstrength[kplus+2]), &D2zFmunu[k]), *TraceF4iD2zF3iPlusF3iD2zF4i);
	      CADD(*TraceF3iD4zF3iMinusF4iD4zF4i, trace_nn(&(s->fieldstrength[kplus]), &D4zFmunu[k]), *TraceF3iD4zF3iMinusF4iD4zF4i);
              CADD(*TraceF4iD4zF3iPlusF3iD4zF4i , trace_nn(&(s->fieldstrength[kplus+2]), &D4zFmunu[k]), *TraceF4iD4zF3iPlusF3iD4zF4i);
	    }
	  else
	    {
	      CSUB(*TraceF3iF3iMinusF4iF4i, trace_nn(&(s->fieldstrength[kplus]),&(s->fieldstrength[kplus])) , *TraceF3iF3iMinusF4iF4i);
	      CADD(*TraceF4iF3iPlusF3iF4i, trace_nn(&(s->fieldstrength[kplus]),&(s->fieldstrength[kplus-2]) ) , *TraceF4iF3iPlusF3iF4i);
	      CSUB(*TraceF3iDzF3iMinusF4iDzF4i, trace_nn(&(s->fieldstrength[kplus]), &DzFmunu[k]), *TraceF3iDzF3iMinusF4iDzF4i);
	      CADD(*TraceF4iDzF3iPlusF3iDzF4i , trace_nn(&(s->fieldstrength[kplus-2]), &DzFmunu[k]), *TraceF4iDzF3iPlusF3iDzF4i);
	      CSUB(*TraceF3iD2zF3iMinusF4iD2zF4i, trace_nn(&(s->fieldstrength[kplus]), &D2zFmunu[k]), *TraceF3iD2zF3iMinusF4iD2zF4i);
	      CADD(*TraceF4iD2zF3iPlusF3iD2zF4i , trace_nn(&(s->fieldstrength[kplus-2]), &D2zFmunu[k]), *TraceF4iD2zF3iPlusF3iD2zF4i);
	      CSUB(*TraceF3iD4zF3iMinusF4iD4zF4i, trace_nn(&(s->fieldstrength[kplus]), &D4zFmunu[k]), *TraceF3iD4zF3iMinusF4iD4zF4i);
              CADD(*TraceF4iD4zF3iPlusF3iD4zF4i , trace_nn(&(s->fieldstrength[kplus-2]), &D4zFmunu[k]), *TraceF4iD4zF3iPlusF3iD4zF4i);
	      
	    }
	  
	}
            
      //printf("\t TraceF4iF3iPlusF3iF4i = ");Display(*TraceF4iF3iPlusF3iF4i);
    }
  cleanup_gather(tagF31Plus); cleanup_gather(tagF32Plus);cleanup_gather(tagF41Plus);cleanup_gather(tagF42Plus);
  cleanup_gather(tagF31Minus);cleanup_gather(tagF32Minus); cleanup_gather(tagF41Minus);cleanup_gather(tagF42Minus);

  g_complexsum(TraceF3iF3iMinusF4iF4i);
  g_complexsum(TraceF4iF3iPlusF3iF4i);
  g_complexsum(TraceF3iDzF3iMinusF4iDzF4i);
  g_complexsum(TraceF4iDzF3iPlusF3iDzF4i);
  g_complexsum(TraceF3iD2zF3iMinusF4iD2zF4i);
  g_complexsum(TraceF4iD2zF3iPlusF3iD2zF4i);
  g_complexsum(TraceF3iD4zF3iMinusF4iD4zF4i);
  g_complexsum(TraceF4iD4zF3iPlusF3iD4zF4i);

  
  //Now Going to compute remaining terms of D4zFmunu i.e. at z+2a and z-2a

  msg_tag *tagLinkZ_zMINUS2a, *tagLinkZ_zPLUSa;  //Gather Link at site=z-2a and z+a in +z direction

  Displacement[0]=0; Displacement[1]=0; Displacement[2]=-2; Displacement[3]=0; //X,Y,Z,T
  tagLinkZ_zMINUS2a = start_general_gather_site( F_OFFSET(link[ZUP]), sizeof(su3_matrix), Displacement, EVENANDODD, gen_pt[0]);
  wait_general_gather(tagLinkZ_zMINUS2a);

  Displacement[0]=0; Displacement[1]=0; Displacement[2]=1; Displacement[3]=0; //X,Y,Z,T
  tagLinkZ_zPLUSa   = start_general_gather_site( F_OFFSET(link[ZUP]), sizeof(su3_matrix), Displacement, EVENANDODD, gen_pt[1]);  
  wait_general_gather(tagLinkZ_zPLUSa);
  FORALLSITES(i, s)
    {
      su3mat_copy( (su3_matrix *)(gen_pt[0][i]), &(LinkZ_zMINUS2a[i])) ;      
      su3mat_copy( (su3_matrix *)(gen_pt[1][i]), &(LinkZ_zPLUSa[i])) ;
    }
  cleanup_general_gather(tagLinkZ_zMINUS2a); cleanup_general_gather(tagLinkZ_zPLUSa); 
  
  Displacement[0]=0; Displacement[1]=0; Displacement[2]=2; Displacement[3]=0; //X,Y,Z,T
  tagF31Plus   = start_general_gather_site( F_OFFSET(fieldstrength[FS_XZ]), sizeof(su3_matrix), Displacement, EVENANDODD, gen_pt[0] );
  wait_general_gather(tagF31Plus);
  tagF32Plus   = start_general_gather_site( F_OFFSET(fieldstrength[FS_YZ]), sizeof(su3_matrix), Displacement, EVENANDODD, gen_pt[1] );
  wait_general_gather(tagF32Plus);
  tagF41Plus   = start_general_gather_site( F_OFFSET(fieldstrength[FS_XT]), sizeof(su3_matrix), Displacement, EVENANDODD, gen_pt[2] );
  wait_general_gather(tagF41Plus);
  tagF42Plus   = start_general_gather_site( F_OFFSET(fieldstrength[FS_YT]), sizeof(su3_matrix), Displacement, EVENANDODD, gen_pt[3] );
  wait_general_gather(tagF42Plus);

  Displacement[0]=0; Displacement[1]=0; Displacement[2]=-2; Displacement[3]=0; //X,Y,Z,T
  tagF31Minus   = start_general_gather_site( F_OFFSET(fieldstrength[FS_XZ]), sizeof(su3_matrix), Displacement, EVENANDODD, gen_pt[4] );
  wait_general_gather(tagF31Minus);
  tagF32Minus   = start_general_gather_site( F_OFFSET(fieldstrength[FS_YZ]), sizeof(su3_matrix), Displacement, EVENANDODD, gen_pt[5] );
  wait_general_gather(tagF32Minus);
  tagF41Minus   = start_general_gather_site( F_OFFSET(fieldstrength[FS_XT]), sizeof(su3_matrix), Displacement, EVENANDODD, gen_pt[6] );
  wait_general_gather(tagF41Minus);
  tagF42Minus   = start_general_gather_site( F_OFFSET(fieldstrength[FS_YT]), sizeof(su3_matrix), Displacement, EVENANDODD, gen_pt[7] );
  wait_general_gather(tagF42Minus);

  // Loop over each site to sum F_mu,nu components 
  FORALLSITES(i, s) 
    {
      mult_su3_nn( &(s->link[ZUP]), &(LinkZ_zPLUSa[i]), Dummy1 );
      mult_su3_nn( &(LinkZ_zMINUS2a[i]), &(LinkZ_zMINUSa[i]), Dummy2);
      for(int k=0; k<4; k++)
	{
	  kplus = k+1;  //&(s->fieldstrength[kplus]), kplus=k+1=1,2,3,4 ; F31, F32, F41, F42
	  mult_su3_nn( Dummy1, (su3_matrix *)(gen_pt[k][i]), Dummy3 );
	  mult_su3_na( Dummy3, Dummy1, &D4zFmunu[k]);

	  mult_su3_an( Dummy2, (su3_matrix *)(gen_pt[k+4][i]), Dummy3 );
          mult_su3_nn( Dummy3, Dummy2, Dummy4);	  
	  add_su3_matrix( &D4zFmunu[k], Dummy4, &D4zFmunu[k]);

	  if(k<2)
	    {
	      CADD(*TraceF3iD4zF3iMinusF4iD4zF4i, trace_nn(&(s->fieldstrength[kplus]), &D4zFmunu[k]), *TraceF3iD4zF3iMinusF4iD4zF4i);
              CADD(*TraceF4iD4zF3iPlusF3iD4zF4i , trace_nn(&(s->fieldstrength[kplus+2]), &D4zFmunu[k]), *TraceF4iD4zF3iPlusF3iD4zF4i);
	    }
	  else
	    {
	      CSUB(*TraceF3iD4zF3iMinusF4iD4zF4i, trace_nn(&(s->fieldstrength[kplus]), &D4zFmunu[k]), *TraceF3iD4zF3iMinusF4iD4zF4i);
              CADD(*TraceF4iD4zF3iPlusF3iD4zF4i , trace_nn(&(s->fieldstrength[kplus-2]), &D4zFmunu[k]), *TraceF4iD4zF3iPlusF3iD4zF4i);
	    }
	}
    }
  
  free(DzFmunu); free(D2zFmunu); free(D4zFmunu);
  free(LinkZ_zMINUSa);  free(LinkZ_zMINUS2a); free(LinkZ_zPLUSa);
  free(Dummy1);  free(Dummy2);  free(Dummy3); free(Dummy4); 
  free(Displacement);
  cleanup_general_gather(tagF31Plus); cleanup_general_gather(tagF32Plus);cleanup_general_gather(tagF41Plus);cleanup_general_gather(tagF42Plus);
  cleanup_general_gather(tagF31Minus);cleanup_general_gather(tagF32Minus); cleanup_general_gather(tagF41Minus);cleanup_general_gather(tagF42Minus);

  g_complexsum(TraceF3iD4zF3iMinusF4iD4zF4i);
  g_complexsum(TraceF4iD4zF3iPlusF3iD4zF4i);
  
  // Sum over all nodes 
  //g_doublesum(time);
  //g_doublesum(space);
  //g_doublesum(charge);

  (*TraceF4iF3iPlusF3iF4i).real = ((*TraceF4iF3iPlusF3iF4i).real)*2.0; //Minus-sign from 2i^2, Note, F4i=-Fi4, F3i=-Fi3
  (*TraceF4iF3iPlusF3iF4i).imag = ((*TraceF4iF3iPlusF3iF4i).imag)*2.0;
  // Normalizations 
  
  (*TraceF3iF3iMinusF4iF4i).real = -((*TraceF3iF3iMinusF4iF4i).real)/(4.0*volume*16.0); //Minus-sign bc of 2i^2 in denominator of Fmunu
  (*TraceF3iF3iMinusF4iF4i).imag = -((*TraceF3iF3iMinusF4iF4i).imag)/(4.0*volume*16.0);
  (*TraceF4iF3iPlusF3iF4i).real  = -((*TraceF4iF3iPlusF3iF4i).real)/(4.0*volume*16.0);
  (*TraceF4iF3iPlusF3iF4i).imag  = -((*TraceF4iF3iPlusF3iF4i).imag)/(4.0*volume*16.0);
  
  (*TraceF3iDzF3iMinusF4iDzF4i).real = -((*TraceF3iDzF3iMinusF4iDzF4i).real)/(4.0*2.0*volume*16.0);
  (*TraceF3iDzF3iMinusF4iDzF4i).imag = -((*TraceF3iDzF3iMinusF4iDzF4i).imag)/(4.0*2.0*volume*16.0);
  (*TraceF4iDzF3iPlusF3iDzF4i).real  = -((*TraceF4iDzF3iPlusF3iDzF4i).real)/(4.0*2.0*volume*16.0);
  (*TraceF4iDzF3iPlusF3iDzF4i).imag  = -((*TraceF4iDzF3iPlusF3iDzF4i).imag)/(4.0*2.0*volume*16.0);

  (*TraceF3iD2zF3iMinusF4iD2zF4i).real = -((*TraceF3iD2zF3iMinusF4iD2zF4i).real)/(4.0*volume*16.0);
  (*TraceF3iD2zF3iMinusF4iD2zF4i).imag = -((*TraceF3iD2zF3iMinusF4iD2zF4i).imag)/(4.0*volume*16.0);
  (*TraceF4iD2zF3iPlusF3iD2zF4i).real  = -((*TraceF4iD2zF3iPlusF3iD2zF4i).real)/(4.0*volume*16.0);
  (*TraceF4iD2zF3iPlusF3iD2zF4i).imag  = -((*TraceF4iD2zF3iPlusF3iD2zF4i).imag)/(4.0*volume*16.0);
    
  (*TraceF3iD4zF3iMinusF4iD4zF4i).real = -((*TraceF3iD4zF3iMinusF4iD4zF4i).real)/(4.0*volume*16.0);
  (*TraceF3iD4zF3iMinusF4iD4zF4i).imag = -((*TraceF3iD4zF3iMinusF4iD4zF4i).imag)/(4.0*volume*16.0);
  (*TraceF4iD4zF3iPlusF3iD4zF4i).real  = -((*TraceF4iD4zF3iPlusF3iD4zF4i).real)/(4.0*volume*16.0);
  (*TraceF4iD4zF3iPlusF3iD4zF4i).imag  = -((*TraceF4iD4zF3iPlusF3iD4zF4i).imag)/(4.0*volume*16.0);

  printf("\n Total TraceF3iF3iMinusF4iF4i = "); Display(*TraceF3iF3iMinusF4iF4i);
  printf("\n Total TraceF4iF3iPlusF3iF4i = ");  Display(*TraceF4iF3iPlusF3iF4i);
  printf("\n Total TraceF3iDzF3iMinusF4iDzF4i = "); Display(*TraceF3iDzF3iMinusF4iDzF4i);
  printf("\n Total TraceF4iDzF3iPlusF3iDzF4i  = "); Display(*TraceF4iDzF3iPlusF3iDzF4i);
  printf("\n Total TraceF3iD2zF3iMinusF4iD2zF4i = "); Display(*TraceF3iD2zF3iMinusF4iD2zF4i);
  printf("\n Total TraceF4iD2zF3iPlusF3iD2zF4i  = "); Display(*TraceF4iD2zF3iPlusF3iD2zF4i);
  printf("\n Total TraceF3iD4zF3iMinusF4iD4zF4i = "); Display(*TraceF3iD4zF3iMinusF4iD4zF4i);
  printf("\n Total TraceF4iD4zF3iPlusF3iD4zF4i  = "); Display(*TraceF4iD4zF3iPlusF3iD4zF4i);
  printf("\n");
}


