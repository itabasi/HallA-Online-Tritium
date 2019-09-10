// Cross Section Calculation macro
// JLab,  p(e,e'K)Lam,Sig reaction
// Not work
void cross_cal(){


  double pk,
    pk=sqrt(pow((pow(W,2)+pow(mk,2)-pow(my,2))/2W,2)-pow(mk,2));


  doble f_Q,f_W,f_res,X;

  X=2.67; // Lambda
  f_Q=1/pow(pow(Q,2)+pow(X,2));
  f_W=pk/W(pow(W,2)-pow(mp,2));
  double C1,C2,A,B;
  A=1.72; // [GeV]
  B=0.10; // [GeV]
  C1=4023.9; //[GeV^2 nb/sr]
  C2=180.0; //[GeV^2 nb/sr]

  f_res=C1*f_W+C2*pow(A*B,2)/(pow(pow(W,2)-pow(A,2))+pow(A*B,2));



}
