//////////
// Comments
// z=cos
// arc(cos)->Fill Uniform
// 
double sphere(double *x, double *par){

  double r=sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  double theta=acos(x[2]/r);
  double sign=1.0;
  if(x[1]<0)sign=1.0;
  double phi=sign*acos(x[0]/sqrt(x[0]*x[0]+ x[1]*x[1]));
  double R=1;
  if(R>r)   return r*r;
  else   return 0;

}

double sphere2(double *x, double *par){

  double X=exp( (pow(x[1]/0.2,2) )/2.0 );
  double Y=0;
  double Z=0;
  double F=X*X + Y*Y + Z*Z;
  return F;
}


double sphere3(double *x, double *par){

  double X=exp( (pow(x[1]/0.2,2) )/2.0 );
  double Y=0;
  double Z=0;
  double F=X*X + Y*Y + Z*Z;
  return F;
}

double sphere4(double *x, double *par){

  double theta = acos(-x[0]);
  double phi   = x[1];
  double r =1.0;
    return r;

}



double uniform(double *x, double* par){

  return x[0];
  
  
}
void fill3D(){

  TF3* f=new TF3("f","sphere",-1,1,-1,1,-1,1,1);
  f->SetLineColor(6);
  f->SetFillColor(6);
  f->SetFillStyle(3002);

  TF3* f2=new TF3("f2","sphere2",0,1,0,3.14/2.0,0,2*3.14,1);
  f2->SetLineColor(6);
  f2->SetFillColor(6);
  f2->SetFillStyle(3002);

  TF3* f3=new TF3("f3","sphere3",-1,1,-1,1,-1,1,1);
  f3->SetLineColor(6);
  f3->SetFillColor(6);
  f3->SetFillStyle(3002);  
  
  TH3D* h=new TH3D("h","Title;X;Y;Z",200,-1,1,200,-1,1,200,-1,1);
  h->FillRandom("f",10000);
  h->SetLineColor(4);
  h->SetFillColor(4);
  h->SetFillStyle(3002);


  TH3D* h2=new TH3D("h2","",200,-1,1,200,-1,1,200,-1,1);
  TH3D* h3=new TH3D("h3","",200,-1,1,200,-1,1,200,-1,1);
  TH2D* hh =new TH2D("hh","",200,-1,1,200,-1,1);
  TH1D* hth =new TH1D("hth","",100,0,3.14);
  TH1D* hph =new TH1D("hph","",100,0,2*3.14);
  TRandom random;
  double seed =  random.Uniform();
  double theta,phi,x,y,z,X,Y,Z;
  for(int i=0;i<10000;i++){
    theta = acos(-1.0*random.Uniform(-1.0,1.0));
    phi = random.Uniform(0.0,2.0*3.14);
    x = sin(theta)*cos(phi);
    y = sin(theta)*sin(phi);
    z = cos(theta);
    X = (x -0.2);
    Y = (y -0.2);
    if(sqrt(X*X +Y*Y)<0.5){
      h2->Fill(x,y,z);
      hh->Fill(x,y);
      hth->Fill(theta);
      hph->Fill(phi);
    }
    //    if(fabs(theta)<0.5)h3->Fill(x,y,z);
  }
  
  TCanvas*c0=new TCanvas("c0","c0");
  c0->Divide(2,1);
  c0->cd(1);
  h2->Draw();
  c0->cd(2);  
  hh->Draw("colz");
  TCanvas*c1=new TCanvas("c1","c1");
  c1->Divide(2,1);
  c1->cd(1);  
  hth->Draw();
  c1->cd(2);  
  hph->Draw();
}
