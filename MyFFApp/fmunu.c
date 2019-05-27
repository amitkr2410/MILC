/************************* fmunu.c *******************************/
/* Calculates the temporal and spatial field strength components */
/* and the topological charge.                                   */

/* Includes */
#include "ks_imp_includes.h"
#include "../include/field_strength.h"
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
void fmunu_fmunu(complex *TraceF3iF3iMinusF4iF4i, complex *TraceF4iF3iPlusF3iF4i)
{
  /* Site variables */
  register int i;
  register site *s;

  /* Temporary component storage */
  su3_matrix *ft, *fs;
  su3_matrix *F31, *F32, *F41, *F42;
  /* Initialize sums */
  //double *time, *space, *charge;
  //*time = *space = *charge = 0;

  /* Compute 8*F_mu,nu at each site */
  make_field_strength( F_OFFSET(link), F_OFFSET(fieldstrength) );

  /* Loop over each site to sum F_mu,nu components */
  FORALLSITES(i, s) {
    
    //mult_su3_nn ( su3_matrix *a, su3_matrix *b, su3_matrix *c );
    F31 = &(s->fieldstrength[FS_XZ]);
    F32 = &(s->fieldstrength[FS_YZ]);
    F41 = &(s->fieldstrength[FS_XT]);
    F42 = &(s->fieldstrength[FS_YT]);

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
    
  }

  (*TraceF4iF3iPlusF3iF4i).real = ((*TraceF4iF3iPlusF3iF4i).real)*2.0;
  (*TraceF4iF3iPlusF3iF4i).imag = ((*TraceF4iF3iPlusF3iF4i).imag)*2.0;
    
  /* Sum over all nodes */
  //g_doublesum(time);
  //g_doublesum(space);
  //g_doublesum(charge);

  /* Normalizations */
  
  (*TraceF3iF3iMinusF4iF4i).real = ((*TraceF3iF3iMinusF4iF4i).real)/(volume);
  (*TraceF3iF3iMinusF4iF4i).imag = ((*TraceF3iF3iMinusF4iF4i).imag)/(volume);
  (*TraceF4iF3iPlusF3iF4i).real  = ((*TraceF4iF3iPlusF3iF4i).real)/(volume);
  (*TraceF4iF3iPlusF3iF4i).imag  = ((*TraceF4iF3iPlusF3iF4i).imag)/(volume);
  
  printf("\n Total TraceF3iF3iMinusF4iF4i = "); Display(*TraceF3iF3iMinusF4iF4i);
  printf("\n Total TraceF4iF3iPlusF3iF4i = "); Display(*TraceF4iF3iPlusF3iF4i);
  printf("\n");
}

