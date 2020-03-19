  // Convert atomic weight to GeV/c^2 

double Bel(int z){
  double Z=(double)z;
  double B= 14.4381*pow(Z,2.39) +1.55468e-6*pow(Z,5.35);
  B=B*1.0e-9; // eV to GeV
  return B;
}

void Al_mass(){

 double  Ma=26.981538; // atomic weight 
 double  u=0.931494013; // GeV
 double Al_z=13.0;
 double me =0.510998928e-3;
 double Mp,ML;
 Mp = 0.938272046;
 ML =  1.115683;   

 
 double mEx_27Al=-17.19675e-3;
 double MAl= Ma*u + mEx_27Al -13.*me + Bel(13); 


  double Mg26=25.98259297;
  double mEx_26Mg = 16.214546e-3;

  double MMgL=Mg26*u + mEx_26Mg + ML -12.0*me + Bel(12);



  cout<<"MAl "<<MAl*1000<<" [MeV]"<<endl;
  cout<<"MMgL "<<MMgL*1000<<" [MeV]"<<endl;



}
