/************************* fmunu.c *******************************/
/* Calculates the temporal and spatial field strength components */
/* and the topological charge.                                   */

/* Includes */
#include "ks_imp_includes.h"
#include "../include/field_strength.h"
#include "../generic/generic_includes.h"
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
void fmunu_fmunu(complex *TraceF3iF3iMinusF4iF4i, complex *TraceF4iF3iPlusF3iF4i, complex *TraceF3iDzF3iMinusF4iDzF4i, complex *TraceF4iDzF3iPlusF3iDzF4i)
{
  printf("Start of fmunu_fmunu function");
  /* Site variables */
  register int i;
  register site *s;
  register site *s1, *s2;
  /* Temporary component storage */
  su3_matrix *ft, *fs;
  su3_matrix *F31, *F32, *F41, *F42;
  /* Initialize sums */
  //double *time, *space, *charge;
  //*time = *space = *charge = 0;
  *TraceF3iF3iMinusF4iF4i=cmplx(0.0,0.0);
  *TraceF4iF3iPlusF3iF4i = cmplx(0.0,0.0);
  *TraceF3iDzF3iMinusF4iDzF4i=cmplx(0.0,0.0);
  *TraceF4iDzF3iPlusF3iDzF4i = cmplx(0.0,0.0);

  /* Compute 8*F_mu,nu at each site */
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
  //su3_matrix *fieldstrength;
  su3_matrix *LinkZF31, *LinkZF32, *LinkZF41, *LinkZF42;
  su3_matrix *DzF31,    *DzF32,    *DzF41,    *DzF42;
  su3_matrix *F31DzF31, *F32DzF32, *F41DzF41, *F42DzF42,  *LinkZ;
  su3_matrix *F31DzF41, *F32DzF42, *F41DzF31, *F42DzF32;
  msg_tag *tagF31, *tagF32, *tagF41, *tagF42;
  tagF31   = start_gather_site( F_OFFSET(fieldstrength[FS_XZ]), sizeof(su3_matrix), ZUP, EVENANDODD, gen_pt[0] );
  tagF32   = start_gather_site( F_OFFSET(fieldstrength[FS_YZ]), sizeof(su3_matrix), ZUP, EVENANDODD, gen_pt[1] );
  
  tagF41   = start_gather_site( F_OFFSET(fieldstrength[FS_XT]), sizeof(su3_matrix), ZUP, EVENANDODD, gen_pt[2] );
  tagF42   = start_gather_site( F_OFFSET(fieldstrength[FS_YT]), sizeof(su3_matrix), ZUP, EVENANDODD, gen_pt[3] );
  
  wait_gather(tagF31); 
  wait_gather(tagF32);
  wait_gather(tagF41);
  wait_gather(tagF42);
  
  // Loop over each site to sum F_mu,nu components 
  FORALLSITES(i, s) {
        
    F31 = &(s->fieldstrength[FS_XZ]);
    F32 = &(s->fieldstrength[FS_YZ]);
    F41 = &(s->fieldstrength[FS_XT]);
    F42 = &(s->fieldstrength[FS_YT]);
    
    LinkZ = &(s->link[ZUP]);
    mult_su3_nn( LinkZ, F31, LinkZF31 );
    mult_su3_nn( LinkZ, F32, LinkZF32 );
    mult_su3_nn( LinkZ, F41, LinkZF41 );
    mult_su3_nn( LinkZ, F42, LinkZF42 );
    
    sub_su3_matrix( (su3_matrix *)(gen_pt[0][i]), LinkZF31, DzF31);
    sub_su3_matrix( (su3_matrix *)(gen_pt[1][i]), LinkZF32, DzF32);
    sub_su3_matrix( (su3_matrix *)(gen_pt[2][i]), LinkZF41, DzF41);
    sub_su3_matrix( (su3_matrix *)(gen_pt[3][i]), LinkZF42, DzF42);
    
    //printf("F31 is \n");DisplayNNMatrix(F31);
    //printf("F32 is \n");DisplayNNMatrix(F32);
    //printf("F41 is \n");DisplayNNMatrix(F41);
    //printf("F42 is \n");DisplayNNMatrix(F42);
    //printf("\n TraceF31F31  = ");Display(trace_nn(F31,F31));    
    CADD(*TraceF3iF3iMinusF4iF4i, trace_nn(F31,F31) , *TraceF3iF3iMinusF4iF4i);
    //printf("\t TraceF3iF3iMinusF4iF4i = ");Display(*TraceF3iF3iMinusF4iF4i);
    
    //printf("\n TraceF32F32  = ");Display(trace_nn(F32,F32));
    CADD(*TraceF3iF3iMinusF4iF4i, trace_nn(F32,F32) , *TraceF3iF3iMinusF4iF4i);
    //printf("\t TraceF3iF3iMinusF4iF4i = ");Display(*TraceF3iF3iMinusF4iF4i);
    
    //printf("\n TraceF41F41  = ");Display(trace_nn(F41,F41));
    CSUB(*TraceF3iF3iMinusF4iF4i, trace_nn(F41,F41) , *TraceF3iF3iMinusF4iF4i);
    //printf("\t TraceF3iF3iMinusF4iF4i = ");Display(*TraceF3iF3iMinusF4iF4i);
    
    //printf("\n TraceF42F42  = ");Display(trace_nn(F42,F42));
    CSUB(*TraceF3iF3iMinusF4iF4i, trace_nn(F42,F42) , *TraceF3iF3iMinusF4iF4i);
    //printf("\t TraceF3iF3iMinusF4iF4i = ");Display(*TraceF3iF3iMinusF4iF4i);

    //printf("\n TraceF41F31  = ");Display(trace_nn(F41,F31));
    CADD(*TraceF4iF3iPlusF3iF4i, trace_nn(F41,F31) , *TraceF4iF3iPlusF3iF4i);
    //printf("\t TraceF4iF3iPlusF3iF4i = ");Display(*TraceF4iF3iPlusF3iF4i);
    //printf("\n TraceF41F31  = ");Display(trace_nn(F41,F31));
    CADD(*TraceF4iF3iPlusF3iF4i, trace_nn(F42,F32) , *TraceF4iF3iPlusF3iF4i);
    //printf("\t TraceF4iF3iPlusF3iF4i = ");Display(*TraceF4iF3iPlusF3iF4i);
    
    CADD(*TraceF3iDzF3iMinusF4iDzF4i, trace_nn(F31, DzF31), *TraceF3iDzF3iMinusF4iDzF4i);
    CADD(*TraceF3iDzF3iMinusF4iDzF4i, trace_nn(F32, DzF32), *TraceF3iDzF3iMinusF4iDzF4i);
    CSUB(*TraceF3iDzF3iMinusF4iDzF4i, trace_nn(F41, DzF41), *TraceF3iDzF3iMinusF4iDzF4i);
    CSUB(*TraceF3iDzF3iMinusF4iDzF4i, trace_nn(F42, DzF42), *TraceF3iDzF3iMinusF4iDzF4i);
    
    CADD(*TraceF4iDzF3iPlusF3iDzF4i , trace_nn(F31, DzF41), *TraceF4iDzF3iPlusF3iDzF4i);
    CADD(*TraceF4iDzF3iPlusF3iDzF4i , trace_nn(F32, DzF42), *TraceF4iDzF3iPlusF3iDzF4i);
    CADD(*TraceF4iDzF3iPlusF3iDzF4i , trace_nn(F41, DzF31), *TraceF4iDzF3iPlusF3iDzF4i);
    CADD(*TraceF4iDzF3iPlusF3iDzF4i , trace_nn(F42, DzF32), *TraceF4iDzF3iPlusF3iDzF4i);
    
    }
  
  cleanup_gather(tagF42);
  cleanup_gather(tagF41);
  cleanup_gather(tagF32);   
  cleanup_gather(tagF31);

  g_complexsum(TraceF3iF3iMinusF4iF4i);    
  g_complexsum(TraceF4iF3iPlusF3iF4i);
  g_complexsum(TraceF3iDzF3iMinusF4iDzF4i);
  g_complexsum(TraceF4iDzF3iPlusF3iDzF4i);
  // Sum over all nodes 
  //g_doublesum(time);
  //g_doublesum(space);
  //g_doublesum(charge);

  (*TraceF4iF3iPlusF3iF4i).real = ((*TraceF4iF3iPlusF3iF4i).real)*2.0;
  (*TraceF4iF3iPlusF3iF4i).imag = ((*TraceF4iF3iPlusF3iF4i).imag)*2.0;
  /* Normalizations */
  
  (*TraceF3iF3iMinusF4iF4i).real = -((*TraceF3iF3iMinusF4iF4i).real)/(4.0*volume*16.0);
  (*TraceF3iF3iMinusF4iF4i).imag = -((*TraceF3iF3iMinusF4iF4i).imag)/(4.0*volume*16.0);
  (*TraceF4iF3iPlusF3iF4i).real  = -((*TraceF4iF3iPlusF3iF4i).real)/(4.0*volume*16.0);
  (*TraceF4iF3iPlusF3iF4i).imag  = -((*TraceF4iF3iPlusF3iF4i).imag)/(4.0*volume*16.0);
  
  (*TraceF3iDzF3iMinusF4iDzF4i).real = ((*TraceF3iDzF3iMinusF4iDzF4i).real)/(4.0*volume*16.0);
  (*TraceF3iDzF3iMinusF4iDzF4i).imag = ((*TraceF3iDzF3iMinusF4iDzF4i).imag)/(4.0*volume*16.0);
  (*TraceF4iDzF3iPlusF3iDzF4i).real  = ((*TraceF4iDzF3iPlusF3iDzF4i).real)/(4.0*volume*16.0);
  (*TraceF4iDzF3iPlusF3iDzF4i).imag  = ((*TraceF4iDzF3iPlusF3iDzF4i).imag)/(4.0*volume*16.0);

  printf("\n Total TraceF3iF3iMinusF4iF4i = "); Display(*TraceF3iF3iMinusF4iF4i);
  printf("\n Total TraceF4iF3iPlusF3iF4i = ");  Display(*TraceF4iF3iPlusF3iF4i);
  printf("\n Total TraceF3iDzF3iMinusF4iDzF4i = "); Display(*TraceF3iDzF3iMinusF4iDzF4i);
  printf("\n Total TraceF4iDzF3iPlusF3iDzF4i  = "); Display(*TraceF4iDzF3iPlusF3iDzF4i);
  printf("\n");
}


