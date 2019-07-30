
void mmass_calc(){
  double Ee,Ee_,Ek,pe,pe_,pk,me,mk;
  double ang=13.2*3.14/180.;
  Ee=4.3185;
  pe_=2.1;
  pk=1.823;
  me=0.51e-3;
  mk=0.493;
  Ek=sqrt(mk*mk+pk*pk);
  Ee_=sqrt(me*me+pe_*pe_);
  pe=sqrt(Ee*Ee-me*me);
  double mt=0.938;
  double mm=sqrt((Ee+mt-Ee_-Ek)*(Ee+mt-Ee_-Ek)-
		 (pe*pe+pe_*pe_+pk*pk-2*pe*pe_*cos(ang)-2*pe*pk*cos(ang)+2*pe_*pk*cos(2*ang)));
  cout<<"mass : "<<mm<<" GeV/c^2"<<endl;
  cout<<"P : "<< (pe*pe+pe_*pe_+pk*pk-2*pe*pe_*cos(ang)-2*pe*pk*cos(ang)+2*pe_*pk*cos(2*ang))<<" E: "<<(Ee+mt-Ee_-Ek)*(Ee+mt-Ee_-Ek)<<endl;
  cout<<"pe*pe_ : "<<2*pe*pe_*cos(ang)<<endl;
}
