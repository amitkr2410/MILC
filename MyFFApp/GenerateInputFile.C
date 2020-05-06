void GenerateInputFile()
{
  int Nt=16, Ns=16;
  int TotalConfiguration=5000;

  //choose if Nt=4
  //
  double beta[18]={5.900,  6.000,  6.050,  6.125,  6.215,  6.285,  6.354,  6.423,  6.515,  6.575,  6.608,  6.664,  6.800,  6.950,  7.150,                         7.280,  7.373,  7.500};
  double   ms[18]={0.1320, 0.1138, 0.1064, 0.0966, 0.0862, 0.0790, 0.0728, 0.0670, 0.0603, 0.0564, 0.0542, 0.0514, 0.0448, 0.0386, 0.0320,                        0.0284, 0.0250, 0.0222};  
  double    T[18]={201,    221,    232,    249,    272,    291,    311,    333,    364,    386,    399,    421,    480,    554,    669,                           753,    819,    918};
  double   ml[18]={0.0};
  string   RationalFunctionFile[18]={"rat.m00660m1320", "rat.m00569m1138", "rat.m00532m1064",  "rat.m00483m0966",  "rat.m00431m0862",  "rat.m00395m0790",  "rat.m00364m0728",  "rat.m00335m0670",  "rat.m003015m0603",  "rat.m00282m0564",  "rat.m00271m0542",  "rat.m00257m0514",  "rat.m00224m0448",  "rat.m00193m0386",  "rat.m00160m0320",  "rat.m00142m0284",  "rat.m00125m0250", "rat.m00111m0222"};
  //  //end of Nt=4

  // choose if Nt=6
  /*
  double Beta[24] = {5.900,  6.000,   6.050,   6.100,   6.150,   6.195,   6.215,   6.245,   6.285,   6.341,   6.423,   6.515,   6.575,  6.608,                      6.664,  6.800,   6.950,   7.150,   7.280,   7.373,   7.500,   7.596,   7.825,   8.000};
  double ms[24]   = {0.1320, 0.1138,  0.1064,  0.0998,  0.0936,  0.0880,  0.0862,  0.0830,  0.0790,  0.0740,  0.0670,  0.0604,  0.0564, 0.0542,                     0.0514, 0.0448,  0.0386,  0.032,   0.0284,  0.0250,  0.0222,  0.0202,  0.0164,  0.0140};
  double T[24]    = {134,    147,     154,     162,     170,     178,     181,     187,     194,     205,     222,     243,     258,    266,                        281,    320,     370,     446,     502,     547,     613,     667,     815,     948};
  double ml[24]   = {0.0};
  string RationalFunctionFile[24]={"rat.m00660m1320", "rat.m00569m1138",  "rat.m00532m1064",  };
  */  //end of Nt=6

  double u0=1.0;
  ofstream fout;
  char InputFileName[10000];
  for(int i=0; i<18; i++)
    {
      ml[i]=ms[i]/20.0;
      sprintf(InputFileName, "Input/inputNt%d_Ns%d_Beta%.4f_ml%.6f_ms%.6f_u0_%.3f.txt", Nt, Ns, beta[i], ml[i], ms[i], u0);
      fout.open(InputFileName, std::ios::out);

      fout<<"prompt 0"<<endl;
      fout<<"nx "<< Ns<<endl;
      fout<<"ny "<< Ns<<endl;
      fout<<"nz "<< Ns<<endl;
      fout<<"nt "<< Nt<<endl;
      fout<<"iseed "<<86658<<endl;
      fout<<"n_pseudo "<<4<<endl;      
      fout<<"load_rhmc_params "<<"RationalFunctionTUMQCD/"<< RationalFunctionFile[i] <<endl;
      fout<<"beta "<< beta[i] <<endl;
      fout<<"n_dyn_masses "<<2<<endl;
      fout<<"dyn_mass "<< ml[i] <<"  "<< ms[i] <<endl;
      fout<<"dyn_flavors "<< 2 <<" "<< 1 <<endl;
      fout<<"u0 "<< u0<<endl;

      fout<<"warms "<< 0 <<endl;
      fout<<"trajecs "<< TotalConfiguration <<endl; 
      fout<<"traj_between_meas  1 "<<endl;
      fout<<"microcanonical_time_step  .02 "<<endl;
      fout<<"steps_per_trajectory  10 "<<endl;
      fout<<"cgresid_md_fa_gr  .0005  .0001  .0001 "<<endl;
      fout<<"max_multicg_md_fa_gr  2500  2500  2500 "<<endl;
      fout<<"cgprec_md_fa_gr 2 2 2"<<endl;
      fout<<"cgresid_md_fa_gr .000005 1e-6 1e-6 "<<endl;
      fout<<"max_multicg_md_fa_gr  2500  2500  2500 "<<endl;
      fout<<"cgprec_md_fa_gr  2 2 2 "<<endl;
      fout<<"cgresid_md_fa_gr .000005 1e-6 1e-6 "<<endl;
      fout<<"max_multicg_md_fa_gr  2500  2500  2500 "<<endl;
      fout<<"cgprec_md_fa_gr  2 2 2 "<<endl;
      fout<<"cgresid_md_fa_gr .000005 1e-6 1e-6 "<<endl;
      fout<<"max_multicg_md_fa_gr  2500  2500  2500 "<<endl;
      fout<<"cgprec_md_fa_gr  2 2 2 "<<endl;
      fout<<"prec_ff 1 "<<endl;
      fout<<"number_of_pbp_masses 2 "<<endl;
      fout<<"max_cg_prop 2000 "<<endl;
      fout<<"max_cg_prop_restarts 5 "<<endl;
      fout<<"npbp_reps 4 "<<endl;
      fout<<"prec_pbp 2 "<<endl;
      fout<<"mass "<< ml[i] <<endl;
      fout<<"naik_term_epsilon 0"<<endl;
      fout<<"error_for_propagator 1e-6 "<<endl;
      fout<<"rel_error_for_propagator 0 "<<endl;
      fout<<"mass "<< ms[i]<<endl;
      fout<<"naik_term_epsilon 0"<<endl;
      fout<<"error_for_propagator 1e-6"<<endl;
      fout<<"rel_error_for_propagator 0"<<endl;
      fout<<"fresh"<<endl;
      fout<<"forget"<<endl;
      fout.close();
    }
}
