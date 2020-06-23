//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double expgaus2(double *x, double *par) {
  //par[0]=Total area
  //par[1]=tau of exp function
  //par[2]=Width (sigma) of convoluted Gaussian function
  //par[3]=Shift of Function Peak
  double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  double np = 500.0;      // number of convolution steps
  double sc =   8.0;      // convolution extends to +-sc Gaussian sigmas
  double xx, fland, sum = 0.0, xlow, xupp, step, i;
  double val;

// Range of convolution integral
  xlow = 0.;
  xupp = x[0] + sc * par[2];
  step = (xupp-xlow) / np;
// Convolution integral
  for(i=1.0; i<=np/2; i++){
     xx = xlow + (i-0.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);

     xx = xupp - (i-.5) * step - par[3];
     fland = TMath::Gaus(xx,x[0],par[2]);
     sum += fland * TMath::Exp(-xx/par[1]);
  }
  //  val = par[2] * step * sum * invsq2pi / par[3];
  val = par[0] * step * sum * invsq2pi / (par[2]*par[1]*exp(par[3]/par[1]));

  return val;
}

void pion_tail(){

    string ifname = "pion.root";

    TFile* ifp=new TFile(ifname.c_str());
    TH1F* hct =((TH1F*)ifp->Get("hct_peak"));
    hct->SetName("hct");
    TF1* fpi=new TF1("fpi","expgaus2",-10,10,4);

    fpi->SetParameter(0,5000);
    fpi->SetParameter(1,0.4);
    fpi->SetParameter(2,0.3);
    TF1* fk =new TF1("fk","gausn(0)",-5,5);
    hct->Fit("fk","","",-1,1);
    hct->Fit("fpi","","",-10,3.5);
    double pk[3];
    pk[0] = fk->GetParameter(0);
    pk[1] = fk->GetParameter(1);
    pk[2] = fk->GetParameter(2);

    TF1* fp =new TF1("fp","gausn(0)",-5,5);
    hct->Fit("fp","","",0,2);
    
    TF1* fpi2 =new TF1("fpi2","gausn(0)+gausn(3)",-5,5);
    fpi2->SetParameters(1650,3.,0.676,3000,3.,0.28);	
    hct->Fit("fpi2","","",2,5);

	
    TF1*fall =new TF1("fall","fpi2+fk",-10,10);
    hct->Fit("fall","","",-22,5);
    TCanvas* c0=new TCanvas("c0","c0");
    c0->cd();
    hct->SetLineColor(1);
    hct->Draw();
    fpi2->SetLineColor(2);
    //    fpi->Draw("same");
    fk->SetLineColor(4);
    fk->Draw("same");
    //    fbg->Draw("same");
    fpi2->Draw("same");
    fall->SetLineColor(3);
    fall->Draw("same");
    double Ntot = hct->Integral(hct->GetXaxis()->FindBin(-1),hct->GetXaxis()->FindBin(1));
    cout<<"=========================="<<endl;
    cout<<"Integral fk "<<fk->Integral(-1,1)*10.<<endl;
    cout<<"Integral hct"<<Ntot<<endl;
    cout<<"pion tail "<<Ntot -fk->Integral(-1,1)*10.<<endl;

}
