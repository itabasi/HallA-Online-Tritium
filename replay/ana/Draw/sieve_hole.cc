

void sieve_hole(){

  bool rarm=true;
  //  bool rarm=false;
  TFile*f1 ;
  //  if(rarm)f1=new TFile(Form("../rootfiles/angcalib/ang_RHRS_woRas.root"));
  //  else {f1=new TFile(Form("../rootfiles/angcalib/ang_LHRS_woRas.root"));}
  f1=new TFile(Form("../rootfiles/zcalib/zt_RHRS_sieve.root"));
  //  f1=new TFile(Form("../angcalib/rootfile/ang_RHRS_5th_0830_0.root"));
  TTree* t1=(TTree*)f1->Get("T");

  TH2F* hist_a=(TH2F*)f1->Get("h3_a");
  TH2F* hist_c=(TH2F*)f1->Get("h3_c");  
  hist_a->SetName("hist_a");
  hist_c->SetName("hist_c");



  Double_t trig5;
  Double_t trig4;
  Double_t trig1;
  double ps_asum,sh_asum;
  double gs_asum;
  double Zt[100];
  double Rvz[100];
  double Lvz[100];
  double Zt_tuned[100];
  double ztR_opt[100];
  double a1,a2;
  t1->SetBranchAddress("DR.T1", &trig1);
  t1->SetBranchAddress("DR.T4", &trig4);
  t1->SetBranchAddress("R.ps.asum_c", &ps_asum);
  t1->SetBranchAddress("R.sh.asum_c", &sh_asum);  
  t1->SetBranchAddress("DR.T5", &trig5);
  t1->SetBranchAddress("R.a1.asum_c", &a1);
  t1->SetBranchAddress("R.a2.asum_c", &a2);    
  if(rarm==true){
   t1->SetBranchAddress("R.tr.vz_opt",ztR_opt);
   t1->SetBranchAddress("R.tr.vz_tuned",Zt);
   // t1->SetBranchAddress("R.tr.vz",Zt); //not tuned 
  } else{
   t1->SetBranchAddress("L.tr.vz_opt",ztR_opt);
   //   t1->SetBranchAddress("L.tr.vz",Zt);
   t1->SetBranchAddress("L.tr.vz_tuned",Zt);
  }

  int ENum=t1->GetEntries();
  cout<<"Entries : "<<ENum<<endl;
  bool tree=false;
  if(tree){

    string ofname;
    if(rarm)    ofname="sieveRHRS_drow.root";
    else     ofname="sieveLHRS_drow.root";
  TFile *ofp = new TFile(Form("%s",ofname.c_str()),"recreate");
  TTree *newtree = t1->CloneTree();


  }
  

  
const double step = 0.492 * 2.54;
const int nrow = 11; // the number of row in SS pattern
const int ncol = 7;  // the number of column in SS pattern
const int nsshole = nrow*ncol; // the number of holes to consider 
double refx[nsshole];
double refy[nsshole];
double selec_widthx = 0.60; // selection width in x (dispersive plane)
double selec_widthy = 0.45; // selection width in y   
 int nhole = 0;
   double ssy_cent_real[nrow];
   double ssx_cent_real[ncol];


  //===== Draw ======//
  TCanvas* c0=new TCanvas("c0","c0");
  TCanvas* c1=new TCanvas("c1","c1");  

  c0->cd();
  hist_a->Draw("colz");
  c1->cd();
  hist_c->Draw("colz");

  TMarker* mark[nsshole];
  

    for(int j=0; j<nrow; j++){
  for(int i=0; i<ncol ; i++){      
      ssy_cent_real[i] = -3.0*step + step*i;
      if(j%2==0)ssy_cent_real[i] = ssy_cent_real[i] - step/2.0;
      ssx_cent_real[j] = 5.0*step - step*j;
      refx[nhole] = ssx_cent_real[j];
      refy[nhole] = ssy_cent_real[i];
      //=== correction error point ====//
      if(j==10 && i==2) refy[nhole]=-1.87452;
      if(j==9  && i==1) refy[nhole]=-2.49936;
      if(j==8  && i==0) refy[nhole]=-4.377388;      
      //===============================//
      //      cout<<"j : "<<j<<" ssx "<<ssx_cent_real[j]<<" i "<<i<<" ssy "<<ssy_cent_real[i]<<endl;
      
      mark[nhole] = new TMarker(refy[nhole],refx[nhole],28);
      


      mark[nhole]->SetMarkerColor(1);
      c0->cd();
      mark[nhole]->Draw("same");
      c1->cd();
      mark[nhole]->Draw("same");      
      nhole++;

    }
  }

      mark[23]->SetMarkerColor(2);
      mark[38]->SetMarkerColor(2);

    
    if(tree){
      hist_a->Write();
      hist_c->Write();	
    }

}
