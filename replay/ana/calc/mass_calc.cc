

void mass_calc(){

  double MN,Ma,Bel;
  int Z,A,L;
  A=1; // # of Proton
  Z=3; // total neucleus
  L=1;  // # of Lambda

  const double Me=0.510998928e6;           // electron      mass (eV/c2)
  const double Mp = 0.938272046e9;         // proton        mass (eV/c2)
  const double Mn = 0.939565379e9;         // neutron       mass (eV/c2)
  const double ML = 1.115683e9;            // Lambda        mass (eV/c2)

  Ma=A*Mp + L*ML + (Z-A-L)*Mn;

  Bel=14.4381*pow(Z,2.39) + 1.55468*1.0e-6*pow(Z,5.35); //eV

  MN=Ma-Z*Me+Bel;

  

  
  cout<<"mass "<<MN*1.e-9<<" [GeV]"<<endl;
  

}
