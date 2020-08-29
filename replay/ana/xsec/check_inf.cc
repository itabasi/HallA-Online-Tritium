

void check_inf(){


    // k=10.
  double k =10.0;
  double delta1=0.11656;
  double delta2=0.08084;
  double delta3=0.02951;

  // k =50

   k =50.0;
   delta1=0.38399;
   delta2=0.36257;
   delta3=0.12950;
  
  
  TF1* hF1 =new TF1("hF1","pow(sin([0]+[1]*x)/sin([1]*x),2.)",-1,1);
  hF1->SetLineColor(1);
  hF1->SetNpx(2000);
  hF1->SetParameter(0,delta1);
  hF1->SetParameter(1,k);

  TF1* hF2 =new TF1("hF2","pow(sin([0]+[1]*x)/sin([1]*x),2.)",-10,10);
  hF2->SetLineColor(2);
  hF2->SetNpx(2000);
  hF2->SetParameter(0,delta2);
  hF2->SetParameter(1,k);  

  TF1* hF3 =new TF1("hF3","pow(sin([0]+[1]*x)/sin([1]*x),2.)",-10,10);
  hF3->SetLineColor(4);
  hF3->SetNpx(2000);
  hF3->SetParameter(0,delta3);
  hF3->SetParameter(1,k);    

  TCanvas* c0 =new TCanvas("c0","c0");
  c0->cd();
  hF1->Draw();
  //  hF2->Draw("same");
  //  hF3->Draw("same");

  cout<<"hF1 (r =3.15) "<<hF1->Eval(3.15)<<endl;
  cout<<"hF2 (r =3.15) "<<hF2->Eval(3.15)<<endl;
  cout<<"hF3 (r =3.15) "<<hF3->Eval(3.15)<<endl;
    
}
