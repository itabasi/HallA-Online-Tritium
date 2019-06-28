// conversion Energy Loss[keV] to Heat deposition [W]
//maked by Itabashi 2018/02/20

void heatdeposit_cal(){

  // One Electron Heatdeposition//

  double dE,eV,J;
  
  cout<<"Energy Loss [keV] is    "<<endl;
  cin>>dE;
  eV=dE*1e3;//[keV] ->[eV]
  J=eV/6.241e18;//[J] converse [eV] -> [J]



  // Total Heat Depositon  //

  double W,I,i,e,Ne;
  e=1.60217662e-19;//[C] elemental value of electron   
  cout<<"Current [myuA] is     "<<endl;
  cin>>I;// I [myuA] Beam Current
  i=I*1e-6;  // [myuA] ->[A]
  Ne=i/e;  // [A] -> [Counts/s]
  W=J*Ne;  // Getten Wat (total)

  cout<<endl;
  cout<<endl;
  cout<<"======RESULT OF CALCULATION=====  "<<endl;
  cout<<endl;
  cout<<endl;
  cout<<"Jule of one electron  is  "<<endl;
  cout<<J<<endl;
  cout<<endl;
  cout<<"Heat Deposition [W] is "<<endl;
  cout<<W<<endl; 
  cout<<endl;
}
