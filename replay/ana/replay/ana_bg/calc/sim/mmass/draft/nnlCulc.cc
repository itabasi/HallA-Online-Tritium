//////QF function >0////////////////////////
double f(double x,double p1,double p2,double p3)
	  {
	    double y;
	    y=p1*pow(x,2)+p2*x+p3;
	    if(x<0){ y=0;
	    }
	    return y;
 }
/////////////////////accidental function/////////////////////////
double f1(double x, double p4,double p5 )
{
  double y1;
  y1=p4*x+p5;
  return y1;
}


void nnlCulc(){
  float mm,a,o,c,o1,e,sig,b;
  for(int ii=0;ii<1000;ii++){
TFile*f1 = new TFile("test_kaon_he3_goga3.root");
  TTree*t1 = (TTree*)f1->Get("h666");
  int n;
  a=0.5;
  e=2;
  b=60;
    o=sqrt(pow(e/1000,2)+pow(a/1000,2));//o:FWHM
  sig=o/2.355;

  o1=1/(sqrt(2*3.14)*o);
  c=0.0;
  t1->SetBranchAddress("mmnuc",&mm);
  double z,q;
  q=1;
  z=q*1.73*80;


  //////////////////////////fitting funktion/////////////////////////////////////
  TF1*g = new TF1("g","[0]*exp(-0.5*pow( (x-[1]) / [2] , 2.0)) ",-1,0.05);
  TF1*fqf= new TF1("fqf","f(x,[0],[1],[2])",-1,0.05);
  TF1*fac=new TF1("fac","f1(x,[0],[1])",-1,0.05);
  fac->SetNpx(3000);
  fac->SetParameters(0,6);
  g->SetParameters(80,c,sig);
  fqf->SetParameters(12000,1000,0);
  fqf->FixParameter(2,0);
  fac->FixParameter(1,1.73*q);
  fac->FixParameter(0,0.0);
  TF1*fit =new TF1("fit","g+fqf+fac",-0.05,0.05);
  fit->FixParameter(6,0.0);
  fit->FixParameter(5,0.0);
  fit->SetParLimits(4,0.0,100000);
  ////////////////////////////hist////////////////////////////////
  TH1F*h1=new TH1F("h1","bounding energy",80,-0.01,0.03);
  TH1F*h2=new TH1F("h2","h2",80,-0.01,0.03);
  TH1F*hs=new TH1F("hs","hs",80,-0.01,0.03);
  TH1F*hf=new TH1F("hf","hf",80,-0.01,0.03);
  ///////////////////// //QF hist////////////////////////////////////////
  int ent;
  ent=t1->GetEntries();
  for(int i=0;i<2000-b;i++){
    t1->GetEntry(i);
    h1->Fill(mm-0.938272-0.939565-1.115683);
 }
  /////// /////////////accidental hist/////////////////////////////////////
  h2->FillRandom("g",b);
  hf->FillRandom("fac",floor(z));
  h1->Add(h2);//peac+QF hist
  h1->Add(hf);//+accidental hist


  h1->Fit("fit","","",-0.05,0.05);
  
 
  ////////////////Get Function parameter//////////////////////
  double mean,sigma,h,herr,qf0,qf1,qf2,ac0,qc1;
  h=fit->GetParameter(0);
  herr=fit->GetParError(0);
  mean=fit->GetParameter(1);//mean
  sigma=fit->GetParameter(2); //sigma
  qf0=fit->GetParameter(3);
  qf1=fit->GetParameter(4);
  qf2=fit->GetParameter(5);
  ac0=fit->GetParameter(6);
  ac1=fit->GetParameter(7);

  double meanerr;
  meanerr=fit->GetParError(1);//delta mean
  double sigmaerr;
  sigmaerr=fit ->GetParError(2);//delta sigma
  xm=mean*1000;
  xme=meanerr*1000;
  xw=sigma*1000;
  xwe=sigmaerr*1000;
  double sq;
  sq=2000*h*sqrt(3.14*2*pow(sigma,2));
  double dN;
  dN=sqrt(2*3.14)*sqrt(pow(sigma,2)*pow(herr,2)+pow(h,2)*pow(sigmaerr,2));
  xn=sq;
  xne=dN;
  /*
  hc1->Fill(xm);
  hc2->Fill(xme);
  hc3->Fill(xw);
  hc4->Fill(xwe);
  hc5->Fill(xn);
  hc6->Fill(xne);
  */
  }
  /*
TCanvas*cc=new TCanvas("cc","cc");
 cc->Divide(3,2);
 cc->cd(1);
 hc1->Draw();
 cc->cd(2);
 hc2->Draw();
 cc->cd(3);
 hc3->Draw();
 cc->cd(4);
 hc4->Draw();
 cc->cd(5);
 hc5->Draw();
 cc->cd(6);
 hc6->Draw();
  */

}
 
