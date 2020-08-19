

void Vpic(){

  TF1* fV1 =new TF1("fV1","[0]*exp(-pow(x/[1],2)) + [2]*exp(-pow(x/[3],2))",0,10);
  TF1* fV2 =new TF1("fV2","[0]*exp(-pow(x/[1],2)) + [2]*exp(-pow(x/[3],2))",0,10);
  TF1* fV3 =new TF1("fV3","[0]*exp(-pow(x/[1],2)) + [2]*exp(-pow(x/[3],2))",0,10);
  //TF1* fV1 =new TF1("fV1","[0]*[1]*[1]*[1]/(8.0*3.14*3.14)*sqrt(3.14)/exp(pow(x*[1]/2.0,2)) + [2]*[3]*[3]*[3]/(8.0*3.14*3.14)*sqrt(3.14)/exp(pow(x*[3]/2.0,2))",0,100);

  TF1* fVp1 =new TF1("fVp1","[0]*[1]*[1]*[1]/(8.0*3.14*3.14)*sqrt(3.14)/exp(pow(x*[1]/197./2.0,2)) + [2]*[3]*[3]*[3]/(8.0*3.14*3.14)*sqrt(3.14)/exp(pow(x*[3]/197./2.0,2))",0,2000);

  TF1* fVp2 =new TF1("fVp2","[0]*[1]*[1]*[1]/(8.0*3.14*3.14)*sqrt(3.14)/exp(pow(x*[1]/197./2.0,2)) + [2]*[3]*[3]*[3]/(8.0*3.14*3.14)*sqrt(3.14)/exp(pow(x*[3]/197./2.0,2))",0,2000);

  TF1* fVp3 =new TF1("fVp3","[0]*[1]*[1]*[1]/(8.0*3.14*3.14)*sqrt(3.14)/exp(pow(x*[1]/197./2.0,2)) + [2]*[3]*[3]*[3]/(8.0*3.14*3.14)*sqrt(3.14)/exp(pow(x*[3]/197./2.0,2))",0,2000);  

  TF1* fVpR =new TF1("fVpR","[0]*[1]*[1]*[1]/(8.0*3.14*3.14)*sqrt(3.14)/exp(pow(x*[1]/197./2.0,2))",0,2000);
  
  fV1->SetParameter(0,-167.34);
  fV1->SetParameter(1,1.1);
  fV1->SetParameter(2,246.8);
  fV1->SetParameter(3,0.82);
  fV1->SetLineColor(1);
  fV1->SetNpx(2000);
  
  fV2->SetParameter(0,-373.94);
  fV2->SetParameter(1,0.791);
  fV2->SetParameter(2,246.8);
  fV2->SetParameter(3,0.82);
  fV2->SetLineColor(2);
  fV2->SetNpx(2000);
    
  fV3->SetParameter(0,-131.49);
  fV3->SetParameter(1,1.095);
  fV3->SetParameter(2,246.8);
  fV3->SetParameter(3,0.82);
  fV3->SetLineColor(4);
  fV3->SetNpx(2000);

  fVp1->SetParameter(0,-167.34);
  fVp1->SetParameter(1,1.1);
  fVp1->SetParameter(2,246.8);
  fVp1->SetParameter(3,0.82);
  fVp1->SetLineColor(1);
  fVp1->SetNpx(2000);


  fVp2->SetParameter(0,-373.94);
  fVp2->SetParameter(1,0.791);
  fVp2->SetParameter(2,246.8);
  fVp2->SetParameter(3,0.82);
  fVp2->SetLineColor(2);
  fVp2->SetNpx(2000);
    
  fVp3->SetParameter(0,-131.49);
  fVp3->SetParameter(1,1.095);
  fVp3->SetParameter(2,246.8);
  fVp3->SetParameter(3,0.82);
  fVp3->SetLineColor(4);
  fVp3->SetNpx(2000);

  fVpR->SetParameter(0,246.8);
  fVpR->SetParameter(1,0.82);
  fVpR->SetLineColor(8);
  fVpR->SetNpx(2000);


  TGraphErrors* grad1=new TGraphErrors();
  double ps1[11]={0.11448,0.22018,0.31094,0.38420,0.44030,0.48123,0.50953,0.5223,0.54223,0.53799,0.52338};
  TGraphErrors* grad1c=new TGraphErrors();
  double ps1c[11]={0.11656,0.22035,0.31112,0.38399,0.43910,0.47819,0.50366,0.51790,0.52093,0.50058,0.46577};
  double rad[11]={10,20,30,40,50,60,70,80,100,120,140};

  TGraphErrors* grad2=new TGraphErrors();
  double ps2[11]={0.08069,0.15906,0.2331,0.30129,0.36270,0.41696,0.46410,0.50446,0.56705,0.60955,0.63658};
  TGraphErrors* grad2c=new TGraphErrors();
  double ps2c[11]={0.08084,0.15907,0.23309,0.30124,0.36257,0.41669,0.46362,0.50368,0.56531,0.60631,0.63124};

  TGraphErrors* grad3=new TGraphErrors();
  double ps3[11]={0.02872,0.05645,0.08231,0.10565,0.12602,0.14325,0.15736,0.16851,0.18310,0.18960,0.19060};
  TGraphErrors* grad3c=new TGraphErrors();
  double ps3c[11]={0.02951,0.05728,0.08598,0.10854,0.12950,0.14601,0.16004,0.16874,0.17445,0.16700,0.14875};  

  grad1->SetName("grad1");
  grad1->SetMarkerColor(1);
  grad1->SetMarkerStyle(23);
  grad1c->SetName("grad1c");
  grad1c->SetMarkerColor(1);
  grad1c->SetMarkerStyle(8);

  grad2->SetName("grad2");
  grad2->SetMarkerColor(2);
  grad2->SetMarkerStyle(23);
  grad2c->SetName("grad2c");
  grad2c->SetMarkerColor(2);
  grad2c->SetMarkerStyle(8);  

  grad3->SetName("grad3");
  grad3->SetMarkerColor(4);
  grad3->SetMarkerStyle(23);
  grad3c->SetName("grad3c");
  grad3c->SetMarkerColor(4);
  grad3c->SetMarkerStyle(8);    
  
  for(int i=0;i<11;i++){
    grad1-> SetPoint(i+1,rad[i],ps1[i]);
    grad1c->SetPoint(i+1,rad[i],ps1c[i]);
    grad2-> SetPoint(i+1,rad[i],ps2[i]);
    grad2c->SetPoint(i+1,rad[i],ps2c[i]);
    grad3-> SetPoint(i+1,rad[i],ps3[i]);
    grad3c->SetPoint(i+1,rad[i],ps3c[i]);    
  }
  
  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();
  fV1->GetXaxis()->SetRangeUser(0,3);
  fV1->GetYaxis()->SetRangeUser(-150,150);
  
  fV1->Draw();
  fV2->Draw("same");
  fV3->Draw("same");


  TCanvas* c1 =new TCanvas("c1","c1");
  c1->cd();
  fVp1->GetYaxis()->SetRangeUser(-2,3);
  fVp1->Draw();
  fVp2->Draw("same");
  fVp3->Draw("same");
  fVpR->Draw("same");

  TFile * ofr=new TFile("../rootfiles/fsi/Vpic.root","recreate");
  fVp1->Write();
  fVp2->Write();
  fVp3->Write();
  fV1->Write();
  fV2->Write();
  fV3->Write();
  grad1->Write();
  grad1c->Write();
  grad2->Write();
  grad2c->Write();
  grad3->Write();
  grad3c->Write();

}
