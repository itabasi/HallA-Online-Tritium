

void mass_calc(){

  double MN,Ma,Bel;
  double Z,A,L,S0;
  Z=13; // # of Proton
  //  A=26; // total neucleus

  L=1;  // # of Lambda
  S0=0;
  const double Me=0.510998928e6;           // electron      mass (eV/c2)
  const double Mp = 0.938272046e9;         // proton        mass (eV/c2)
  const double Mn = 0.939565379e9;         // neutron       mass (eV/c2)
  const double ML = 1.115683e9;            // Lambda        mass (eV/c2)
  const double MS0 = 1.192642e9;             // Sigma         mass (eV/c2) 
  const double Mu  = 0.931494061e9;        // atomics mass

  const double Ma_27Al = 26.98153853; // atomics mass
  const double Ma_26Mg = 25.98259297; // atomics mass


  //  Ma=A*Mp + L*ML + (Z-A-L)*Mn;
  
  Ma = Ma_27Al* Mu;
  //  Ma = Ma_26Mg* Mu +ML;

  Bel=14.4381*pow(Z,2.39) + 1.55468*1.0e-6*pow(Z,5.35); //eV
  MN=Ma-Z*Me+Bel;

  

  
  cout<<"mass "<<setprecision(8)<<MN*1.e-9<<" [GeV]"<<endl;
  cout<<"Bel "<<Bel<<" [eV]"<<endl;

}
