  // Convert atomic weight to GeV/c^2 

void mass_conv(){

 double  Ma=26.981538; // atomic weight 
 double  u=0.931494013; // MeV

  cout<<" Mass "<<Ma*u<<endl;
  double Mp,ML;
  Mp = 0.938272046;
  ML =  1.115683;   


  double Mg26=25.98259297;
  double MMgL=Mg26*u + ML;



  cout<<"MgL mass "<<MMgL<<endl;
}
