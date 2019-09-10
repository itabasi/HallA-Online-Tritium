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


void nnL_peak_sim(){
  TCanvas *nnlambda=new TCanvas("nnlambda","nnlambda");
  //c10copy->Divide(2,1);
  // c10copy->cd(1);
  TFile*f1 = new TFile("test_kaon_he3_goga3.root");
  TTree*t1 = (TTree*)f1->Get("h666");
  float mm,a,o,c,o1,e,sig,b;


  /*
  ////// Coment out //////////////////////////////
  int n;
   cout<<"Natural wide[MeV] is  "<<endl;
   cin>>a;
   cout<<"Detecte's wide[MeV] is "<<endl;
   cin>>e;
  cout<<"nnL Counts are"<<endl;
  cin>>b;
    o=sqrt(pow(e/1000,2)+pow(a/1000,2));//o:FWHM
  sig=o/2.355;
  o1=1/(sqrt(2*3.14)*o);
 double qf;
  cout<<"Quasi-free counts are "<<endl;
  cin>>qf;
 
  double z,q;
  cout<<"Accidental B.G ratio is"<<endl;
  cin>>q;
  z=q*1.73*80; // 80 is hist counts from -0.01 to 0.03 

   cout<<"Binding Energy[MeV] is"<<endl;
   cin>>c;
   c=c/1000;//MeV->GeV
  */
  int n;
  a=0.5;
  e=2.0;
  b=120;
  double z,q;
  q=1.0;
  o=sqrt(pow(e/1000,2)+pow(a/1000,2));//o:FWHM
  sig=o/2.355;
  o1=1/(sqrt(2*3.14)*o);
  double qf=4000;

  int bin_mm;
  double min_mm,max_mm;
  min_mm=-0.05;
  max_mm=0.2;
  bin_mm=(int)((max_mm-min_mm)*2000);

  //  bin_mm=420;


   z=q*1.73*bin_mm; 
   c=1.0/1000;//MeV->GeV

  t1->SetBranchAddress("mmnuc",&mm);

  

  //////////////////////////fitting funktion/////////////////////////////////////
  TF1*g = new TF1("g","[0]*exp(-0.5*pow( (x-[1]) / [2] , 2.0)) ",-1,0.05);
  TF1*fqf= new TF1("fqf","f(x,[0],[1],[2])",-1,0.05);
  TF1*fac=new TF1("fac","f1(x,[0],[1])",-1,0.05);
  fac->SetNpx(3000);
  //fac->SetParameters(0,b/4/sig);
  fac->SetParameters(0,6);
  //TF1*fac2 =new TF1("fac2","[0]*x[1]",-0.05,0.05);
  g->SetParameters(1000,c,sig);
  fqf->SetParameters(12000,1000,0);
  fqf->FixParameter(2,0);
  //fac->FixParameter(1,b*1.25*q/sig/10000);
  fac->FixParameter(1,1.73*q);
  fac->FixParameter(0,0.0);
  // fac->FixParameter(1,0);
  //fac->SetParLimits(1,b*1.25*q/sig/10000-1,b*1.25*q/sig/10000+1);
  //fac->FixParameter(1,4);
  //fac2->SetParameters(0,b*1000/4/sig);
  g->SetLineColor(4);
  fqf->SetLineColor(2);
  fac->SetLineColor(5);
  // fqf->Draw();
   //fac->Draw("same");
  //g->Draw("same");
  //g->FixParameter(0,400);
  //fitting function
  //TF1*fit2 = new TF1("fit2","fqf+g+fac",-0.05,0.05);

  // fit2->SetNpx(1000);
  //fit2->SetLineColor(8); 
  //fit2->GetXaxis()->SetRangeUser(-0.05,0.05);
  //fit2->Draw("same");
  //fit->GetXaxis()->SetRangeUser(-0.005,0.03);
  TF1*fit =new TF1("fit","g+fqf+fac",-0.05,0.05);
//g->[0]-[2],fqf->[3]-[5],fac->[6]-[7]/////////////////
  fit->FixParameter(6,0.0);
  fit->FixParameter(5,0.0);
  //fit->FixParameter(7,0);
  fit->SetParLimits(4,0.0,100000);
  fit->SetLineColor(6);
  // fit2->Draw("same");
  // c10copy->cd(2);
  ////////////////////////////hist////////////////////////////////

  //  min_mm=-0.01;
  //  max_mm=0.03;
  //  bin_mm=80;
    //(int)(max_mm-min_mm)*2000;
  cout<<"bin : "<<bin_mm<<endl;

  TH1F*h1=new TH1F("h1","Binding Energy",bin_mm,min_mm,max_mm);
  TH1F*h2=new TH1F("h2","h2",bin_mm,min_mm,max_mm);
  TH1F*hs=new TH1F("hs","hs",bin_mm,min_mm,max_mm);
  TH1F*hf=new TH1F("hf","hf",bin_mm,min_mm,max_mm);
  //  h1->SetFillColor(3);
  ///////////////////// //QF hist////////////////////////////////////////
  int ent;
  ent=t1->GetEntries();
  for(int i=0;i<qf-b;i++){
    t1->GetEntry(i);
    h1->Fill(mm-0.938272-0.939565-1.115683);
 }
  /*int max ,qq;
  qq =h1->GetMaximumBin();
  max= h1->GetBinContent(qq);
  cout<<"maximumBin "<<max<<endl;
  double ma,ss;
  cout<<"aciidental rate is  "<<ss<<endl;
  cin>>ss;
  ma=max/10*ss;*/
  //peak hist
  /////// /////////////accidental hist/////////////////////////////////////
  h2->FillRandom("g",b);
  hf->FillRandom("fac",floor(z));
  h1->Add(h2);//peac+QF hist
  h1->Add(hf);//+accidental hist


  //  h1->Fit("fit","Rq","",-0.05,0.05);
  //h1->Fit("fac","R+","",-0.05,c-2*sig);
  
 
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
  cout<<"sigma="<<sig*1000<<endl;
  cout<<"mu="<<c*1000<<endl;
  cout<<"sigma[MeV]"<<sigma*1000<<endl;
  cout<<"sigmaerr[MeV]"<<sigmaerr*1000<<endl;
  cout<<"sigma relative[%]"<<sigmaerr*100/sigma<<endl;
  cout<<"mean [MeV]"<<mean*1000<<endl;
  cout<<"meanerr[MeV] "<<meanerr*1000<<endl;
  cout<<"mean relative[%]"<< meanerr*100/mean<<endl;
  cout<<"z= "<<z<<endl;
  cout<<"Nacc=  "<<ac1*2000*4*sig<<endl;


  g->SetParameters(h,mean,sigma);

  fqf->SetParameters(qf0,qf1,qf2);
  fac->SetParameters(ac0,ac1);
  //  h1->GetXaxis()->SetRangeUser(-0.01,0.03);
  h1->GetXaxis()->SetTitle("Binding Energy [GeV]");
  h1->GetXaxis()->CenterTitle();
  h1->GetXaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetTitle("Counts/ 0.5 MeV");
  h1->GetYaxis()->CenterTitle();
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->SetMinimum(0);
  h1->Draw();
  hf->SetLineColor(2);
  //  hf->Draw("same");

  cout<<"z= "<<z<<endl;
  double sq;
  sq=2000*h*sqrt(3.14*2*pow(sigma,2));
  cout<<"Count is "<<sq<<endl;
  double dN;
  dN=sqrt(2*3.14)*sqrt(pow(sigma,2)*pow(herr,2)+pow(h,2)*pow(sigmaerr,2));
  cout<<"Count error is "<<dN*1000<<endl;
  cout<<"Counts relative[%]"<<(dN/sq)*100000<<endl;
  cout<<"N rate "<<sq/sqrt(ac1*2000*4*sig)<<endl;
  cout<<"c-2*sig=  "<<c-2*sig<<endl;
 int sumi;
 int w;
 int r;
 int ww,www;
 r=0;
  for( int k =0; k<8000*sig; k++)
    {
      // w=c-2*sig+0.0005*k;
      ww =h1->GetXaxis()->FindBin(c-2*sig+0.0005*k);
      // cout<<"Binnumber "<<ww<<endl;
    www =h1->GetBinContent(ww);
    // cout<<"Bin counts "<<www<<endl;
     sumi=sumi+www;
     r=r+1;
    
     }
  cout<<"r=  "<<r<<endl;
  cout<<"sum    "<<sumi<<endl;
  cout<<"sumerr  "<<sq/sqrt(sumi)<<endl;



  }
