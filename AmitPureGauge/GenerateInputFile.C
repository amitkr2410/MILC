void GenerateInputFile()
{
  int Nt=32, Ns=32;
  int TotalConfiguration=10000;

  //choose if Nt=4
  //

  //  //end of Nt=4

  // choose if Nt=6
  //

  //  //end of Nt=6

  //Choose if Nt=8
  // /* 
  double Beta[10]= {5.7,      5.95,     6,       6.1,      6.2,       6.35,     6.55,     6.7,      6.85,    6.95};

   //*/
  double u0=1.0;
  ofstream fout;
  char InputFileName[10000];
  int NumberOfFile=sizeof(Beta)/sizeof(Beta[0]);

  for(int i=0; i< NumberOfFile; i++)
    {
      sprintf(InputFileName, "Input/inputNt%d_Ns%d_Beta%.4f.txt", Nt, Ns, Beta[i]);
      fout.open(InputFileName, std::ios::out);

      fout<<"prompt 0"<<endl;
      fout<<"nx "<< Ns<<endl;
      fout<<"ny "<< Ns<<endl;
      fout<<"nz "<< Ns<<endl;
      fout<<"nt "<< Nt<<endl;
      fout<<"iseed "<<86658<<endl;
      fout<<"warms 0"<<endl;
      fout<<"trajecs "<<TotalConfiguration<<endl;
      fout<<"traj_between_meas 1"<<endl;
      fout<<"beta "<< Beta[i] <<endl;
      fout<<"steps_per_trajectory 10"<<endl;
      fout<<"qhb_steps 4"<<endl;
      fout<<"fresh"<<endl;
      fout<<"no_gauge_fix"<<endl;
      fout<<"forget"<<endl;

      fout.close();
    }
}
