

/*
double y[5];
void chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
  double chsq=0.0;
  double func=0.0;
  for(int i=0;i<5;i++){
    func = par[0]*i+par[1];
    chisq += pow(y[i]-func/5.0,2);

    }

  f =chisq;
};
*/


int nmax=500;
void vdc_track(){
  
  //  string root_name="/data/VDC/root/LU1t0_100ns/tritium_111500_Lu1t0_100ns.root"; //100 ns shift
  //   string root_name="/data/VDC/small/initial/tritium_111500_initial.root"; // w/o shift
  string root_name="/data/VDC/root/LU1t0tuned/tritium_111500_Lu1t0tuned.root"; //Lu0 tuned    
 //  TFile* fin=new TFile(root_name.c_str(),"read");
  TChain* T=new TChain("T");
    T->Add(root_name.c_str());
  double evtype, HallA_p;  
  double Ru1_nhit[nmax],Ru2_nhit[nmax],Rv1_nhit[nmax],Rv2_nhit[nmax];
  int NRu1_nhit,NRu2_nhit,NRv1_nhit,NRv2_nhit;
  double Ru1_wire[nmax],Ru2_wire[nmax],Rv1_wire[nmax],Rv2_wire[nmax];
  int NRu1_wire,NRu2_wire,NRv1_wire,NRv2_wire;
  double Ru1_rtime[nmax],Ru2_rtime[nmax],Rv1_rtime[nmax],Rv2_rtime[nmax];
  int NRu1_rtime,NRu2_rtime,NRv1_rtime,NRv2_rtime;
  double Ru1_time[nmax],Ru2_time[nmax],Rv1_time[nmax],Rv2_time[nmax];
  int NRu1_time,NRu2_time,NRv1_time,NRv2_time;

  double Lu1_nhit[nmax],Lu2_nhit[nmax],Lv1_nhit[nmax],Lv2_nhit[nmax];
  int NLu1_nhit,NLu2_nhit,NLv1_nhit,NLv2_nhit;
  double Lu1_wire[nmax],Lu2_wire[nmax],Lv1_wire[nmax],Lv2_wire[nmax];
  int NLu1_wire,NLu2_wire,NLv1_wire,NLv2_wire;
  double Lu1_rtime[nmax],Lu2_rtime[nmax],Lv1_rtime[nmax],Lv2_rtime[nmax];
  int NLu1_rtime,NLu2_rtime,NLv1_rtime,NLv2_rtime;
 double Lu1_time[nmax],Lu2_time[nmax],Lv1_time[nmax],Lv2_time[nmax];
  int NLu1_time,NLu2_time,NLv1_time,NLv2_time;      
  double trig;
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("DR.evtypebits"          ,1);            T->SetBranchAddress("DR.evtypebits"          ,&trig); 
  T->SetBranchStatus("R.vdc.u1.nhit"          ,1);            T->SetBranchAddress("R.vdc.u1.nhit"          ,Ru1_nhit);
   T->SetBranchStatus("R.vdc.u2.nhit"          ,1);            T->SetBranchAddress("R.vdc.u2.nhit"          ,Ru2_nhit);
   T->SetBranchStatus("R.vdc.v1.nhit"          ,1);            T->SetBranchAddress("R.vdc.v1.nhit"          ,Rv1_nhit);
   T->SetBranchStatus("R.vdc.v2.nhit"          ,1);            T->SetBranchAddress("R.vdc.v2.nhit"          ,Rv2_nhit);
  T->SetBranchStatus("R.vdc.u1.wire"          ,1);            T->SetBranchAddress("R.vdc.u1.wire"          ,Ru1_wire);
  T->SetBranchStatus("R.vdc.u2.wire"          ,1);            T->SetBranchAddress("R.vdc.u2.wire"          ,Ru2_wire);
  T->SetBranchStatus("R.vdc.v1.wire"          ,1);            T->SetBranchAddress("R.vdc.v1.wire"          ,Rv1_wire);
  T->SetBranchStatus("R.vdc.v2.wire"          ,1);            T->SetBranchAddress("R.vdc.v2.wire"          ,Rv2_wire);
  T->SetBranchStatus("Ndata.R.vdc.u1.wire"          ,1);      T->SetBranchAddress("Ndata.R.vdc.u1.wire"          ,&NRu1_wire);
  T->SetBranchStatus("Ndata.R.vdc.u2.wire"          ,1);      T->SetBranchAddress("Ndata.R.vdc.u2.wire"          ,&NRu2_wire);
  T->SetBranchStatus("Ndata.R.vdc.v1.wire"          ,1);      T->SetBranchAddress("Ndata.R.vdc.v1.wire"          ,&NRv1_wire);
  T->SetBranchStatus("Ndata.R.vdc.v2.wire"          ,1);      T->SetBranchAddress("Ndata.R.vdc.v2.wire"          ,&NRv2_wire);
  T->SetBranchStatus("R.vdc.u1.rawtime"          ,1);         T->SetBranchAddress("R.vdc.u1.rawtime"          ,Ru1_rtime);         
  T->SetBranchStatus("R.vdc.u2.rawtime"          ,1);         T->SetBranchAddress("R.vdc.u2.rawtime"          ,Ru2_rtime);         
  T->SetBranchStatus("R.vdc.v1.rawtime"          ,1);         T->SetBranchAddress("R.vdc.v1.rawtime"          ,Rv1_rtime);         
  T->SetBranchStatus("R.vdc.v2.rawtime"          ,1);         T->SetBranchAddress("R.vdc.v2.rawtime"          ,Rv2_rtime);         
  T->SetBranchStatus("Ndata.R.vdc.u1.rawtime"          ,1);   T->SetBranchAddress("Ndata.R.vdc.u1.rawtime"          ,&NRu1_rtime);         
  T->SetBranchStatus("Ndata.R.vdc.u2.rawtime"          ,1);   T->SetBranchAddress("Ndata.R.vdc.u2.rawtime"          ,&NRu2_rtime);         
  T->SetBranchStatus("Ndata.R.vdc.v1.rawtime"          ,1);   T->SetBranchAddress("Ndata.R.vdc.v1.rawtime"          ,&NRv1_rtime);         
  T->SetBranchStatus("Ndata.R.vdc.v2.rawtime"          ,1);   T->SetBranchAddress("Ndata.R.vdc.v2.rawtime"          ,&NRv2_rtime);         
  T->SetBranchStatus("R.vdc.u1.time"          ,1);         T->SetBranchAddress("R.vdc.u1.time"          ,Ru1_time);         
  T->SetBranchStatus("R.vdc.u2.time"          ,1);         T->SetBranchAddress("R.vdc.u2.time"          ,Ru2_time);         
  T->SetBranchStatus("R.vdc.v1.time"          ,1);         T->SetBranchAddress("R.vdc.v1.time"          ,Rv1_time);         
  T->SetBranchStatus("R.vdc.v2.time"          ,1);         T->SetBranchAddress("R.vdc.v2.time"          ,Rv2_time);         
  //========== LHRS VDC ===============//

  T->SetBranchStatus("L.vdc.u1.nhit"          ,1);            T->SetBranchAddress("L.vdc.u1.nhit"          ,Lu1_nhit);
  T->SetBranchStatus("L.vdc.u2.nhit"          ,1);            T->SetBranchAddress("L.vdc.u2.nhit"          ,Lu2_nhit);
  T->SetBranchStatus("L.vdc.v1.nhit"          ,1);            T->SetBranchAddress("L.vdc.v1.nhit"          ,Lv1_nhit);
  T->SetBranchStatus("L.vdc.v2.nhit"          ,1);            T->SetBranchAddress("L.vdc.v2.nhit"          ,Lv2_nhit);
  T->SetBranchStatus("L.vdc.u1.wire"          ,1);            T->SetBranchAddress("L.vdc.u1.wire"          ,Lu1_wire);
  T->SetBranchStatus("L.vdc.u2.wire"          ,1);            T->SetBranchAddress("L.vdc.u2.wire"          ,Lu2_wire);
  T->SetBranchStatus("L.vdc.v1.wire"          ,1);            T->SetBranchAddress("L.vdc.v1.wire"          ,Lv1_wire);
  T->SetBranchStatus("L.vdc.v2.wire"          ,1);            T->SetBranchAddress("L.vdc.v2.wire"          ,Lv2_wire);
  T->SetBranchStatus("Ndata.L.vdc.u1.wire"          ,1);      T->SetBranchAddress("Ndata.L.vdc.u1.wire"          ,&NLu1_wire);
  T->SetBranchStatus("Ndata.L.vdc.u2.wire"          ,1);      T->SetBranchAddress("Ndata.L.vdc.u2.wire"          ,&NLu2_wire);
  T->SetBranchStatus("Ndata.L.vdc.v1.wire"          ,1);      T->SetBranchAddress("Ndata.L.vdc.v1.wire"          ,&NLv1_wire);
  T->SetBranchStatus("Ndata.L.vdc.v2.wire"          ,1);      T->SetBranchAddress("Ndata.L.vdc.v2.wire"          ,&NLv2_wire);
  T->SetBranchStatus("L.vdc.u1.rawtime"          ,1);         T->SetBranchAddress("L.vdc.u1.rawtime"          ,Lu1_rtime);         
  T->SetBranchStatus("L.vdc.u2.rawtime"          ,1);         T->SetBranchAddress("L.vdc.u2.rawtime"          ,Lu2_rtime);         
  T->SetBranchStatus("L.vdc.v1.rawtime"          ,1);         T->SetBranchAddress("L.vdc.v1.rawtime"          ,Lv1_rtime);         
  T->SetBranchStatus("L.vdc.v2.rawtime"          ,1);         T->SetBranchAddress("L.vdc.v2.rawtime"          ,Lv2_rtime);         
  T->SetBranchStatus("Ndata.L.vdc.u1.rawtime"          ,1);   T->SetBranchAddress("Ndata.L.vdc.u1.rawtime"          ,&NLu1_rtime);        
  T->SetBranchStatus("Ndata.L.vdc.u2.rawtime"          ,1);   T->SetBranchAddress("Ndata.L.vdc.u2.rawtime"          ,&NLu2_rtime);        
  T->SetBranchStatus("Ndata.L.vdc.v1.rawtime"          ,1);   T->SetBranchAddress("Ndata.L.vdc.v1.rawtime"          ,&NLv1_rtime);        
  T->SetBranchStatus("Ndata.L.vdc.v2.rawtime"          ,1);   T->SetBranchAddress("Ndata.L.vdc.v2.rawtime"          ,&NLv2_rtime);        
  T->SetBranchStatus("L.vdc.u1.time"          ,1);         T->SetBranchAddress("L.vdc.u1.time"          ,Lu1_time);         
  T->SetBranchStatus("L.vdc.u2.time"          ,1);         T->SetBranchAddress("L.vdc.u2.time"          ,Lu2_time);         
  T->SetBranchStatus("L.vdc.v1.time"          ,1);         T->SetBranchAddress("L.vdc.v1.time"          ,Lv1_time);         
  T->SetBranchStatus("L.vdc.v2.time"          ,1);         T->SetBranchAddress("L.vdc.v2.time"          ,Lv2_time);         

  


  //===== Make Hist ===============//
  double min_nhit=0.0;
  double max_nhit=5.0;
  int bin_nhit=5;
  double min_time=-1000.;
  double max_time =1000.;
  int bin_time=(int)(max_time-min_time)*2;
  
  TH2F* hLu1tdc_nhit=new TH2F("hLu1tdc_nhit","",bin_nhit,min_nhit,max_nhit,bin_time,min_time,max_time);
  TH1F* hLu1tdc[5];
  TF1* fit[5];
  TF1* ftdc=new TF1("ftdc","[0]*x+[1]",0,5);

    int ENum =T->GetEntries();
  cout<<"Entry : "<<ENum<<endl;

  int MAX=10000;
  TF1* fhit[MAX];
  TGraph* gLu1_ev[MAX];
  TH2F* hLu1_ev[MAX];
  
  for(int i=0;i<5;i++){
    hLu1tdc[i]=new TH1F(Form("hLu1tdc_%d",i),"",bin_time,min_time,max_time);
    hLu1tdc[i]->SetTitle("LU1 TDC [ns] ; TDC [ns] ; Counts ");
    fit[i]=new TF1(Form("fit_%d",i),"gausn(0)",min_time,max_time);


  }

  TGraphErrors* gLu1=new TGraphErrors();
  gLu1->SetTitle("LHRS U1 VDC TDC peak ; #hits ; TDC peak [ns] ");
  
  double time[MAX][5];
  double ns=1.0e9;


  int nn=0;
  
  for(int k=0;k<ENum;k++){
    T->GetEntry(k);
    if(NLu1_wire==5 && trig==32){
      hLu1tdc[0]->Fill(-Lu1_time[0]*ns);
      hLu1tdc[1]->Fill(-Lu1_time[1]*ns);
      hLu1tdc[2]->Fill(Lu1_time[2]*ns);
      hLu1tdc[3]->Fill(Lu1_time[3]*ns);
      hLu1tdc[4]->Fill(Lu1_time[4]*ns);      
      hLu1tdc_nhit->Fill(0.0,-Lu1_time[0]*ns);
      hLu1tdc_nhit->Fill(1.0,-Lu1_time[1]*ns);
      hLu1tdc_nhit->Fill(2.0,Lu1_time[2]*ns);
      hLu1tdc_nhit->Fill(3.0,Lu1_time[3]*ns);
      hLu1tdc_nhit->Fill(4.0,Lu1_time[4]*ns);


    if(nn<MAX){

      gLu1_ev[nn]=new TGraph();
      gLu1_ev[nn]->SetTitle("LU1 5hists TDC [ns] ; #hits ; TDC [ns]");
      gLu1_ev[nn]->SetPoint(0,0.0,-Lu1_time[0]*ns);
      gLu1_ev[nn]->SetPoint(1,1.0,-Lu1_time[1]*ns);
      gLu1_ev[nn]->SetPoint(2,2.0,Lu1_time[2]*ns);
      gLu1_ev[nn]->SetPoint(3,3.0, Lu1_time[3]*ns);
      gLu1_ev[nn]->SetPoint(4,4.0, Lu1_time[4]*ns);
      time[nn][0]=-Lu1_time[0]*ns;
      time[nn][1]=-Lu1_time[1]*ns;      
      time[nn][2]= Lu1_time[2]*ns;
      time[nn][3]= Lu1_time[3]*ns;
      time[nn][4]= Lu1_time[4]*ns;      
    }
    
    nn++;

      
    }
  }



  double dist[MAX],a[MAX],b[MAX];
  double conv=4.24; //[mm]
  TGraph* gdist=new TGraph();
  TH1F* hchi2=new TH1F("hchi2","",1000,0.0,5000);
  TH1F* hdist=new TH1F("hdist","",400,-2.0,4.0);
  hchi2->SetTitle("Tracking #chi^{2} ; #chi^{2};Counts");  
  hdist->SetTitle("Tracking distance ; dist [mm];Counts");

  /*
  TMinuit *min=new TMinuit(2);
  for(int i=0;i<MAX;i++){
    
    for(int j=0;j<5;j++)y[j]=time[i][j];
    
    min->SetPrintLevel(1);
    min->SetFCN(chi2);

    // ~~~ Chi-square ~~~~
    arglist[0] = 1;
    int ierflag=0;    
    minuit -> mnexcm("SET ERR",arglist,1,ierflg);
    minuit -> SetPrintLevel(-1);
    

    double vstart[3];
    vstart[0] = 1.1;
    vstart[1] = 0.1;
    vstart[2] = 0.1;
    double step[3];
    step[0] = 0.01;
    step[1] = 0.01;
    step[2] = 0.01;
    min->mnparm(0, "p0", vstart[0], step[0], 0, 0, ierflg);
    min->mnparm(1, "p1", vstart[1], step[1], 0, 0, ierflg);

    arglist[0] = 1000;//引数1:maxcalls
    arglist[1] = 1;//引数2:tolerance
    min->mnexcm("MIGRAD", arglist, 2, ierflg);

  }
  */
  double chi2[MAX];
  
  for(int i=0;i<MAX;i++){

    fhit[i]=new TF1(Form("fhit_%d",i),"[0]*x+[1]",0,5);  
    gLu1_ev[i]->Fit(Form("fhit_%d",i),"QR","",0.0,4.0);
    chi2[i]= fhit[i]->GetChisquare();

    
    a[i]=fhit[i]->GetParameter(0);
    b[i]=fhit[i]->GetParameter(1);
    dist[i]=b[i]/a[i];
    gdist->SetPoint(i,i,dist[i]);
    hchi2->Fill(chi2[i]);
    hdist->Fill(-(dist[i]+2.)*conv);    
    //    cout<<"dist : "<<dist[i]<<endl;
  }
  

  double max[5],max_err[5];
  for(int i=0;i<5;i++){
    max[i]=hLu1tdc[i]->GetBinCenter(hLu1tdc[i]->GetMaximumBin());

    fit[i]->SetParameter(1,max[i]);
    fit[i]->SetParameter(2,40);
    hLu1tdc[i]->Fit(Form("fit_%d",i),"QR","",max[i]-200.,max[i]+200);
    max[i]=fit[i]->GetParameter(1);
 //    max_err[i]=fit[i]->GetParError(1);    
    cout<<"Get mean parameter "<<i<<" : "<<max[i]<<endl;
    gLu1->SetPoint(i,i,max[i]);
    //    gLu1->SetPointError(i,0,max_err[i]);    
  }




  

  TCanvas* c0=new TCanvas("c0","c0");
  c0->cd();


  hLu1tdc_nhit->SetTitle("TDC vs nhits in LU1; #hits ; TDC [ns]");
  hLu1tdc_nhit->Draw("colz");
  
  TCanvas* c1=new TCanvas("c1","c1");
  c1->cd();
  hchi2->Draw();

  TCanvas* c2=new TCanvas("c2","c2");
  c2->cd();
  gLu1_ev[100]->SetMarkerSize(1.0);
  gLu1_ev[100]->SetMarkerColor(2);
  gLu1_ev[100]->SetMarkerStyle(20);  
  gLu1_ev[100]->Draw("AP");
  //  gLu1->Fit("ftdc","RQ","",0,5); 
  fhit[100]->Draw("same");


  TCanvas* c3=new TCanvas("c3","c3");
  c3->cd();
  gdist->SetMarkerSize(1.0);
  gdist->SetMarkerColor(2);
  gdist->SetMarkerStyle(20);  
  //  gdist->Draw("AP");
  hdist->Draw();
}
