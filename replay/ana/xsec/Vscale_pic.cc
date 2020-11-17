const  int nmax =30;

void Vscale_pic(){



  TF1* fV[nmax];
  
  double va_s = -167.34; // Verma Va Singlet
  double va_t = -132.42; // Verma Va Triplet
  double b_as = 1.1;     // Verma beta Singlet
  double b_at = 1.1;     // Verma beta Triplet
  double re  = 2.91 ;   // Verma effective range
  double fac,fac_s,fac_t;
  
  double a_s = -2.29; // Verma scattering range Singlet
  double r_s =  3.15; // Verma effective  range Singlet
  double a_t = -1.77; // Verma scattering range Triplet
  double r_t =  3.25; // Verma effective  range Singlet

  double ss;
  
  for(int i=0;i<nmax;i++){

    fV[i] = new TF1(Form("fV_%d",i),"[0]*exp(-pow(x/[1],2)) + [2]*exp(-pow(x/[3],2))",0,10);
    ss = (1.0 + (double)i/(double)nmax);
  fV[i]->SetParameter(0,va_s*ss);
  fV[i]->SetParameter(1,b_as);
  fV[i]->SetParameter(2,246.8);
  fV[i]->SetParameter(3,0.82);
  fV[i]->SetLineColor(i+1);
  fV[i]->SetNpx(2000);
    
  }
  
  

  string rname="../rootfiles/fsi/ifac/vscale/verma.root";


  TFile * f = new TFile(rname.c_str());
  TGraphErrors* gV[nmax];

  for(int i=0;i<nmax;i++){

    gV[i] =(TGraphErrors*)f->Get(Form("gIs_%d",i));
    gV[i]->SetName(Form("gV_%d",i));
    gV[i]->SetMarkerStyle(7);
    gV[i]->SetMarkerColor(i+1);
    
  }



  // delta [rad]

  
  
  
  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();
  fV[0]->GetXaxis()->SetRangeUser(0,3);
  fV[0]->GetYaxis()->SetRangeUser(-150,150);
  fV[0]->Draw();
  fV[6]->Draw();
  fV[7]->Draw("same");
  fV[8]->Draw("same");
  //  for(int i=1;i<nmax;i++){
  //    fV[i]->Draw("same");
  //  }

  TCanvas* c1 = new TCanvas("c1","c1");
  c1->cd();
  gV[0]->GetXaxis()->SetRangeUser(0,500);
  gV[0]->GetYaxis()->SetRangeUser(0,10);
  gV[0]->Draw("AP");
  for(int i=1;i<nmax;i++){
    gV[i]->Draw("P");
  }

  gV[6]->Draw("AP");
  gV[7]->Draw("P");
  gV[8]->Draw("P");
}
