using namespace std;
#include "t0corr.h"
#include "Param.h"
#include <TMinuit.h>


double coin_offset;


///////////////////////////////////////////////////////////////////////

void t0corr::SetRoot(string ifname){
  add_tree(ifname);
  pack_tree();
  readtreeHRSR();
  readtreeHRSL();

}
//////////////////////////////////////////////////////////////////////////

void t0corr::SetRunList(string ifname){

  ifstream ifp(Form("%s",ifname.c_str()),ios::in);
  if(!ifp){ cout<<"no input file "<<ifname<<endl; exit(1); }
  string buf, runname;
  while(1){
    getline(ifp,buf);
    if( buf[0]=='#' ){ continue; }
    if( ifp.eof() ) break;
    istringstream sbuf(buf);
    sbuf >> runname;
    add_tree(runname);
  }

  pack_tree();
  readtreeHRSR();
  readtreeHRSL();
}


////////////////////////////////////////////////////////////////////////////

void t0corr::ReadParam(string name){

  param = new ParamMan(name.c_str());
  cout<<"param name : "<<name<<endl;
  if(param -> SetVal())cout<<"F1TDC parameter setted"<<endl; 
  tdc_time=param->F1Res();
  coin_offset=param->GetF1CoinOffset();
  cout<<" coin_offset : "<<coin_offset<<endl; 
  for(int i=0;i<ns2seg;i++){
    LS2T_off[i]=param->GetTdcOffset(1, i , 0, 0);
    LS2B_off[i]=param->GetTdcOffset(1, i , 0, 1);
    RS2T_off[i]=param->GetTdcOffset(1, i , 1, 0);
    RS2B_off[i]=param->GetTdcOffset(1, i , 1, 1);
  }
}

/////////////////////////////////////////////////////////////////////////////


void t0corr::NewRoot(string ofname){

  ofr = new TFile(Form("%s",ofname.c_str()),"recreate");
  tnew =new TTree("T","Coincalib matrix tuning");

  tnew=tree->CloneTree(0);
  tnew->Branch("ct",&coint);
  tnew->Branch("ct_c",&coint_c);
  //  tnew->Branch("RS2T",RS2_F1time);
  //  tnew->Branch("LS2T",LS2_F1time);
  tnew->Branch("RS2T_ref",&RF1Ref[0]);
  tnew->Branch("LS2T_ref",&LF1Ref[0]);
  tnew->Branch("rtof",&Rtof);
  tnew->Branch("ltof",&Ltof);
  tnew->Branch("rtof_c",&Rtof_c);
  tnew->Branch("ltof_c",&Ltof_c);
  tnew->Branch("rpathl",&R_pathl);
  tnew->Branch("lpathl",&L_pathl);  
  tnew->Branch("rpathl_c",&R_pathl_c);
  tnew->Branch("lpathl_c",&L_pathl_c);
  tnew->Branch("Rp",&Rp);
  tnew->Branch("Beta_R",&Beta_R);
  tnew->Branch("Beta_L",&Beta_L);
  
  
}



////////////////////////////////////////////////////////////////////////////




void t0corr::MakeHist(){

  min_ct=-20;
  max_ct= 20.;
  bin_ct=1000.;



  //======= Definition of Hist ============//


  for(int i=0;i<ns2seg;i++){
    hRs2t[i]=new TH1D(Form("hRs2t_%d",i),Form("hRs2t_%d; coin time [ns] ;Counts",i),bin_ct,min_ct,max_ct);
    hRs2b[i]=new TH1D(Form("hRs2b_%d",i),Form("hRs2b_%d; coin time [ns] ;Counts",i),bin_ct,min_ct,max_ct);
    hLs2t[i]=new TH1D(Form("hLs2t_%d",i),Form("hLs2t_%d; coin time [ns] ;Counts",i),bin_ct,min_ct,max_ct);
    hLs2b[i]=new TH1D(Form("hLs2b_%d",i),Form("hLs2b_%d; coin time [ns] ;Counts",i),bin_ct,min_ct,max_ct);



    fRs2t[i]=new TF1(Form("fRs2t_%d",i),"gausn(0)",min_ct,max_ct);
    fRs2b[i]=new TF1(Form("fRs2b_%d",i),"gausn(0)",min_ct,max_ct);
    fLs2t[i]=new TF1(Form("fLs2t_%d",i),"gausn(0)",min_ct,max_ct);
    fLs2b[i]=new TF1(Form("fLs2b_%d",i),"gausn(0)",min_ct,max_ct);
   

    fRs2t[i]->SetLineColor(2);
    fRs2b[i]->SetLineColor(2);
    fLs2t[i]->SetLineColor(2);
    fLs2b[i]->SetLineColor(2);

 

  }


  hcoin=new TH1D("hcoin","",bin_ct,min_ct,max_ct);
  hcoin_k=new TH1D("hcoin_k","",bin_ct,min_ct,max_ct);
  hcoin_r=new TH2D("hcoin_r","",bin_ct,min_ct,max_ct,ns2seg,0,ns2seg);
  hcoin_l=new TH2D("hcoin_l","",bin_ct,min_ct,max_ct,ns2seg,0,ns2seg);
  fcoin= new TF1("fcoin","gausn(0)",min_ct,max_ct);
  fcoin->SetLineColor(2);
  fcoin->SetNpx(2000);
  fcoin_k= new TF1("fcoin","gausn(0)",min_ct,max_ct);
  fcoin_k->SetLineColor(2);
  fcoin_k->SetNpx(2000);
  gRct=new TGraphErrors();
  gLct=new TGraphErrors();
  gRct->SetName("gRct");
  gLct->SetName("gLct");
  gRct->SetMarkerStyle(21);
  gRct->SetMarkerColor(kRed);
  gRct->SetMarkerSize(0.5);
  gLct->SetMarkerStyle(21);
  gLct->SetMarkerColor(kBlue);
  gLct->SetMarkerSize(0.5);

}


/////////////////////////////////////////////////////////////  




void t0corr::Fit(){


  cout<<"================================="<<endl;
  cout<<"======= Fitting ================="<<endl;
  cout<<"================================="<<endl;



  for(int i=0;i<ns2seg;i++){


  //===== Get Peak Position ======//

    Rs2t_posi[i]=hRs2t[i]->GetBinCenter(hRs2t[i]->GetMaximumBin());
    Rs2t_posi[i]=hRs2b[i]->GetBinCenter(hRs2b[i]->GetMaximumBin());
    Ls2t_posi[i]=hLs2t[i]->GetBinCenter(hLs2t[i]->GetMaximumBin());
    Ls2b_posi[i]=hLs2b[i]->GetBinCenter(hLs2b[i]->GetMaximumBin());
  
    hRs2t[i]->Fit(Form("fRs2t_%d",i),"QR","RQ",Rs2t_posi[i]-1.0, Rs2t_posi[i]+1.0);
    hRs2b[i]->Fit(Form("fRs2b_%d",i),"QR","RQ",Rs2b_posi[i]-1.0, Rs2b_posi[i]+1.0);
    hLs2t[i]->Fit(Form("fLs2t_%d",i),"QR","RQ",Ls2t_posi[i]-1.0, Ls2t_posi[i]+1.0);
    hLs2b[i]->Fit(Form("fLs2b_%d",i),"QR","RQ",Ls2b_posi[i]-1.0, Ls2b_posi[i]+1.0);

    Rs2t_posi[i]=fRs2t[i]->GetParameter(1);
    Rs2b_posi[i]=fRs2b[i]->GetParameter(1);
    Ls2b_posi[i]=fLs2b[i]->GetParameter(1);
    Ls2t_posi[i]=fLs2t[i]->GetParameter(1);

    gRct->SetPoint(i,i,Rs2t_posi[i]);
    gLct->SetPoint(i,i,Ls2t_posi[i]);


  }

  double mean_ct=hcoin->GetBinCenter(hcoin->GetMaximumBin());
  hcoin->Fit("fcoin","RQ","RQ",mean_ct-1.,mean_ct+1.0);
  double mean_ct_k =hcoin_k->GetBinCenter(hcoin_k->GetMaximumBin());
  //  double n_ct_k =hcoin_k->GetBinContent
  fcoin_k->SetParameter(1,mean_ct_k);
  hcoin_k->Fit("fcoin_k","RQ","RQ",mean_ct_k-1.,mean_ct_k+1.0);




}

//////////////////////////////////////////////////////////////////////////

void t0corr::PathCalc(int rhit, int lhit){


  R_pathl= R_tr_pathl[rhit] + R_s2_trpath[rhit];
  L_pathl= L_tr_pathl[lhit] + L_s2_trpath[lhit];
  R_pathl_c= R_tr_pathl[rhit] + R_s2_trpath[rhit];
  L_pathl_c= L_tr_pathl[lhit] + L_s2_trpath[lhit];
  
  
  R_tr_x[rhit]  = (R_tr_x[rhit]-XFPm)/XFPr;
  R_tr_th[rhit] = (R_tr_th[rhit]-XpFPm)/XpFPr;
  R_tr_y[rhit]  = (R_tr_y[rhit]-YFPm)/YFPr;
  R_tr_ph[rhit] = (R_tr_ph[rhit]-YpFPm)/YpFPr;
  R_tr_vz[rhit] = (R_tr_vz[rhit] - Ztm)/Ztr;
  
  L_tr_x[lhit]  = (L_tr_x[lhit]-XFPm)/XFPr;
  L_tr_th[lhit] = (L_tr_th[lhit]-XpFPm)/XpFPr;
  L_tr_y[lhit]  = (L_tr_y[lhit]-YFPm)/YFPr;
  L_tr_ph[lhit] = (L_tr_ph[lhit]-YpFPm)/YpFPr;
  L_tr_vz[lhit] = (L_tr_vz[lhit] - Ztm)/Ztr;      
  
  
  
  //==== Calc Path Length =========//
  
      
  
      
  R_tr_x[rhit]  = R_tr_x[rhit] * XFPr + XFPm;
  R_tr_th[rhit] = R_tr_th[rhit] * XpFPr + XpFPm;
  R_tr_y[rhit]  = R_tr_y[rhit] * YFPr + YFPm;
  R_tr_ph[rhit] = R_tr_ph[rhit] * YpFPr   + YpFPm;
  R_tr_vz[rhit] = R_tr_vz[rhit] * Ztr + Ztm;
  
  L_tr_x[lhit]  = L_tr_x[lhit] * XFPr + XFPm;
  L_tr_th[lhit] = L_tr_th[lhit] * XpFPr + XpFPm;
  L_tr_y[lhit]  = L_tr_y[lhit] * YFPr + YFPm;
  L_tr_ph[lhit] = L_tr_ph[lhit] * YpFPr   + YpFPm;
  L_tr_vz[lhit] = L_tr_vz[lhit] * Ztr + Ztm;
  
  R_pathl_c = R_pathl_c * PaRr + PaRm;
  L_pathl_c = L_pathl_c * PaLr + PaLm;



  R_pathl_c = R_pathl_c + R_s2_trpath[rhit];
  L_pathl_c = L_pathl_c + L_s2_trpath[lhit];
      
}





///////////////////////////////////////////////////////////////////////////

void t0corr::Fill(){



  cout<<"=================================="<<endl;
  cout<<"======== Fill Events ============="<<endl;
  cout<<"=================================="<<endl;

  ENum=100000;
  cout<<"Events : "<<ENum<<endl;

  for(int k=0;k<ENum;k++){


    // Initialization //
    
    for(int i=0;i<MAX;i++){
    R_tr_x[i]  = -2000.0;
    R_tr_th[i] = -2000.0;
    R_tr_y[i]  = -2000.0;
    R_tr_ph[i] = -2000.0;
    R_tr_p[i]  = -2000.0;
    L_tr_x[i]  = -2000.0;
    L_tr_th[i] = -2000.0;
    L_tr_y[i]  = -2000.0;
    L_tr_ph[i] = -2000.0;
    L_tr_p[i]  = -2000.0;
    R_s2_trpath[i]=-2000.;
    L_s2_trpath[i]=-2000.;
    }

    
    
    tree->GetEntry(k);


    int NLtr = (int)L_tr_n;  if(NLtr>MAX) NLtr = MAX;
    int NRtr = (int)R_tr_n;  if(NRtr>MAX) NRtr = MAX;

    for(int lt=0;lt<NLtr;lt++){
      for(int rt=0;rt<NRtr;rt++){
	
	Rtof=-1000.0; Rtof_c=-1000.;Ltof=-1000., Ltof_c=-1000.;
	coint=-1000.0; coint_c=-1000.0;
	R_pathl=-1000.0; R_pathl_c=-1000.; L_pathl=-1000.; L_pathl_c=-1000.;
	Rp=-1000.;
	pid_flag=false;
	z_flag=false;

	if(fabs(R_tr_vz[rt] + L_tr_vz[lt])/2.0<0.1 && fabs(R_tr_vz[rt]-L_tr_vz[lt])<0.025)z_flag=true;
	  if(R_a1_asum_p<50. && R_a2_asum_p>2000. && L_cer_asum_c>2000.)pid_flag=true;

	int R_s2pad=(int)R_s2_t_pads[rt];
	int L_s2pad=(int)L_s2_t_pads[lt]; 

   
	Rp=R_tr_p[rt];	
	CoinCalc(R_s2pad, L_s2pad, rt, lt);


	//====== Fill Hist ======//



	hRs2t[R_s2pad]->Fill(coint);
	hRs2b[R_s2pad]->Fill(coint);
	hLs2t[L_s2pad]->Fill(coint);
	hLs2b[L_s2pad]->Fill(coint);
	hcoin->Fill(coint);
	hcoin_r->Fill(coint, R_s2pad);
	hcoin_l->Fill(coint, L_s2pad);

	if(z_flag && pid_flag){
	  hcoin_k->Fill(coin_k);
	}

	tnew->Fill();
      }
    }
    
    if(k%100000==0)cout<<"Event Fill : "<<k<<" / "<<ENum<<endl;

 
  }// End Fill
  


}


////////////////////////////////////////////////////////////////////////////





void t0corr::GetParam(string ofname){


  ofp= new ofstream(ofname.c_str());

  mean_pi=0.0;
  ref_ch=mean_pi/tdc_time;

  string  name[4];
  name[0] = "R.s2.R.off_F1 = " ; 
  name[1] = "R.s2.L.off_F1 = " ; 
  name[2] = "L.s2.R.off_F1 = " ; 
  name[3] = "L.s2.L.off_F1 = " ; 



  for(int i=0;i<ns2seg;i++){


    

    Rs2t_posi_ch[i]= RS2T_off[i] - Rs2t_posi[i]/tdc_time - ref_ch;
    Rs2b_posi_ch[i]= RS2B_off[i] - Rs2b_posi[i]/tdc_time - ref_ch;
    Ls2t_posi_ch[i]= LS2T_off[i] + Ls2t_posi[i]/tdc_time - ref_ch;
    Ls2b_posi_ch[i]= LS2B_off[i] + Ls2b_posi[i]/tdc_time - ref_ch;

    Offset[0][i]=Rs2t_posi_ch[i];
    Offset[1][i]=Rs2b_posi_ch[i];
    Offset[2][i]=Ls2t_posi_ch[i];
    Offset[3][i]=Ls2b_posi_ch[i];



  }



  for(int j=0;j<4;j++){
    *ofp << name[j];
    for(int i=0;i<ns2seg;i++)     
      *ofp << Offset[j][i]<<" ";
   
     *ofp <<endl;

  }
  


}


///////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////

void t0corr::CoinCalc(int RS2_seg, int LS2_seg, int rhit, int lhit){

  // ==== Initialization =====//
  
  Rs2_t=-100.;
  Ls2_t=-100.;
  Beta_R=-100.;
  Beta_L=-100.;
  Beta_K=-100.;
  Rtof=-100.;
  Ltof=-100.;
  Rtof_c=-100.;
  Ltof_c=-100.;
  coint=-1000.;
  coint_c=-1000.;
  coin_k=-1000.;
  
  convertF1TDCR(param);
  convertF1TDCL(param);
  
  PathCalc(rhit,lhit);
  
  Beta_K = R_tr_p[rhit]/sqrt(R_tr_p[rhit]*R_tr_p[rhit]+MK*MK);
  Beta_R = R_tr_p[rhit]/sqrt(R_tr_p[rhit]*R_tr_p[rhit]+Mpi*Mpi);
  Beta_L = L_tr_p[lhit]/sqrt(L_tr_p[lhit]*L_tr_p[lhit]+Me*Me);


  
  //====== w/o Path Calibration =========//  
  Rtof = RS2_F1time[RS2_seg] - R_pathl/(Beta_R*LightVelocity);
  Ltof = LS2_F1time[LS2_seg] - L_pathl/(Beta_L*LightVelocity);
  Rtof_K = RS2_F1time[RS2_seg] - R_pathl/(Beta_K*LightVelocity);

  Rtof_c = RS2_F1time[RS2_seg] - R_pathl/(Beta_R*LightVelocity);
  Ltof_c = LS2_F1time[LS2_seg] - L_pathl/(Beta_L*LightVelocity); 
  
  Rs2_t=RS2_F1time[RS2_seg];
  Ls2_t=LS2_F1time[LS2_seg];


  if(RS2_F1time[RS2_seg]!=-9999. &&LS2_F1time[LS2_seg]!=-9999.){
    coint = - Rtof + Ltof - coin_offset;
    coint_c = - Rtof_c + Ltof_c - coin_offset;
    coin_k = - Rtof_K + Ltof - coin_offset;
  }
  else{
    coint=-1000;

  }



  
}

//////////////////////////////////////////////////////////////////////


void t0corr::Write(){

  for(int i=0; i<ns2seg;i++){

    hRs2t[i]->Write();
    hRs2b[i]->Write();
    hLs2t[i]->Write();
    hLs2b[i]->Write();
    fRs2t[i]->Write();
    fRs2b[i]->Write();
    fLs2t[i]->Write();
    fLs2b[i]->Write();

}

  hcoin->Write();
  hcoin_k->Write();
  hcoin_r->Write();
  hcoin_l->Write();
  fcoin->Write();
  fcoin_k->Write();
  gRct->Write();
  gLct->Write();
  tnew->Write();
}

/////////////////////////////////////////////////////////////////////////

void t0corr::Draw(){

  for(int i=0;i<ncanvas;i++){
    c[i]=new TCanvas(Form("c%d",i),Form("c%d",i)); 
    if(i<ns2seg){
    c[i]->Divide(2,2);
    c[i]->cd(1);
    hRs2t[i]->Draw();
    fRs2t[i]->Draw("same");
    c[i]->cd(2);
    hRs2b[i]->Draw();
    fRs2b[i]->Draw("same");
    c[i]->cd(3);
    hLs2t[i]->Draw();
    fLs2t[i]->Draw("same");
    c[i]->cd(4);
    hLs2b[i]->Draw();
    fLs2b[i]->Draw("same");
    }

  }
  
  c[16]->Divide(2,1);
  c[16]->cd(1);
  hcoin->Draw();
  fcoin->Draw("same"); 
  c[16]->cd(2);
  hcoin_k->Draw();
  fcoin_k->Draw("same"); 
  c[17]->Divide(2,1);
  c[17]->cd(1);
  hcoin_r->Draw("colz");
  c[17]->cd(2);
  hcoin_l->Draw("colz");
  c[18]->cd();
  gRct->Draw("AP");
  gLct->Draw("P");

  


}


////////////////////////////////////////////////////////////////////////



void t0corr::Print(string ofname){

  cout<<endl;
  cout<<"============================="<<endl;
  cout<<"======== Print =============="<<endl;
  cout<<"============================="<<endl;

  c[0]->Print(Form("%s[",ofname.c_str()));
  for(int i=0;i<ncanvas;i++)
    c[i]->Print(Form("%s",ofname.c_str()));

  c[ncanvas-1]->Print(Form("%s]",ofname.c_str()));

  cout<<"print is done "<<endl;

}


////////////////////////////////////////////////////////////////////////

void t0corr::Close(){
  ofr->Close();
  ofp->close();
  
}

/////////////////////////////////////////////////////////////////////////

//=======================================================================//
//===============   Main   ==============================================//
//=======================================================================//


int main(int argc, char** argv){



  gStyle->SetOptFit(111111111);
  int ch;
  //  double tdc_time=58.0e-3;//[ns]
  string ifname;
  string ofname;

  string opt_file="./scale_offset_20190210.dat";  

  bool output_flag = false;
  bool output_tree_flag = false;
  bool draw_flag = false;
  bool root_flag=false;
  bool matrix_flag=false;
  bool RHRS_flag=false; 
  bool tuning_flag=false;
  string pngname;
  extern char *optarg;
  string F1tdc="1";
  int f1tdc=1;
  string ofparam;
  string root_init="../rootfiles/";
  string root_end=".root";
  string dat_init="../matrix/";
  string dat_end=".dat";
  string matrix="matrix/zt_RHRS_2.dat";
  string pname;
  string print_name;
  bool print=false;

  while((ch=getopt(argc,argv,"h:P:s:w:t:p:f:r:z:L:R:m:o:O:i:bcoC"))!=-1){
    switch(ch){
      
      
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;

    case 's':
      root_flag = true;
      draw_flag = false;
      single=true;
      ifname = optarg;
      cout<<"output root filename : "<<ofname<<endl;      
      break;      

      
    case 'r':
      root_flag = true;
      draw_flag = false;
      ofname = optarg;
      cout<<"output root filename : "<<ofname<<endl;      
      break;


    case 'p':
      pname = optarg;
      break;
       
    case 'P':
      print_name =optarg;
      print=true;
      break;

    case 'w':
      ofparam = optarg;
      break;

    case 'o':
      root_flag = true;
      draw_flag = false;
      matrix_flag=true;
      ofname = optarg;
      ofname =root_init+ ofname+ root_end;
      cout<<"output root filename : "<<ofname<<endl;      

      break;
      
    case 'b':
      draw_flag = false;
      cout<<"BACH MODE!"<<endl;
      break;
      


    case 'h':
      cout<<"-f : input root  filename"<<endl;
      cout<<"-p : input matrix filename"<<endl;      
      cout<<"-w : output matrix filename"<<endl;
      cout<<"-r : output root filename"<<endl;
      cout<<"-t : target mode H:hydrogen T:trititum"<<endl;
      cout<<"-m : tuning mode C:R&L-HRS, L:L-HRS, R:R-HRS "<<endl;      
      cout<<"-o : output root & tuning file name "<<endl;
      return 0;
      break;
    case '?':
      cout<<"unknown option...."<<endl;
      return 0;
      break;
    default:
      cout<<"type -h to see help!!"<<endl;
      return 0;

    }
  }




  
  TApplication *theApp =new TApplication("App",&argc,argv);
  gSystem->Load("libMinuit");
  
  t0corr* T0Corr = new t0corr();
  if(single)T0Corr->SetRoot(ifname);
  else T0Corr->SetRunList(ifname);
  T0Corr->NewRoot(ofname);
  T0Corr->MakeHist();
  T0Corr->ReadParam(pname);
  T0Corr->Fill();
  T0Corr->Fit();
  T0Corr->GetParam(ofparam);
  T0Corr->Write();
  T0Corr->Draw();
  if(print)T0Corr->Print(print_name);
  T0Corr->Close();

  cout<<"=====================================" <<endl;
  cout<<"========== Comment Out ==============" <<endl;
  cout<<"=====================================" <<endl;
  cout<<endl;
  cout<<"new root : "<<ofname.c_str()<<endl;
  cout<<"new pdf  : "<<print_name.c_str()<<endl;

  
  gSystem->Exit(1);
  theApp->Run();

  return 0;
  
}//end main



