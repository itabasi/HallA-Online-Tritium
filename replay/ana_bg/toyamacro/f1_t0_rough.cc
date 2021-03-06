//////////////////////////////////////////////////////////////
//t0 calibration of S0&S2 F1TDC(rough parameter as 1st step)//
//                                    by Y. Toyama Nov. 2018//
//////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <iomanip>
#include <csignal>
#include <stdlib.h>
#include <climits>
using namespace std;

#include "TApplication.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TArrow.h"

#include "Tree.h"
#include "Setting.h"
#include "ParamMan.h"

#define Calibration

const int NCanvas = 5;//num of canvas

class f1tdc_t0_rough : public Tree
{
 public:
         f1tdc_t0_rough();
        ~f1tdc_t0_rough();
  void makehist();
  void loop();
  void fit();
  void draw(); 
  void savecanvas(string ofname); 
  void SetMaxEvent( int N )  { ENumMax = N; }
  void SetRoot(string ifname);
  void WriteParam(string ofname);

  private:
    int GetMaxEvent() { return ENumMax; }
    int ENumMax;

    TH1F *h_Rs2t[16], *h_Rs2b[16];
    TH1F *h_Ls2t[16], *h_Ls2b[16];
    TH1F *h_Rs0t, *h_Rs0b;
    TH1F *h_Ls0t, *h_Ls0b;
    Setting *set;
    ParamMan *param;
    TCanvas *c[NCanvas];
    TArrow *arrow;

};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
f1tdc_t0_rough::f1tdc_t0_rough()
{

  gErrorIgnoreLevel = kError;
  gROOT->SetStyle("Plain");
  gROOT->SetBatch(1);

  gStyle->SetOptDate(0);
  gStyle->SetOptFit(1);
  gStyle->SetHistFillStyle(3002);
  gStyle->SetHistFillColor(0);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetOptDate(0);
//  gStyle->SetStatW(0.15);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleTextColor(1);
  gStyle->SetGridWidth(1);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetNdivisions(510); // tertiary*10000 + secondary*100 + first
  gStyle->SetOptStat("iMenr");
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.13);

  for(int i=0;i<NCanvas;i++){
    c[i]= new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1400,800 );
  }

  set = new Setting();
  param = new ParamMan("param/default.param");

  arrow = new TArrow(0,0,0,0,0.01,"");
  arrow ->SetAngle(45);
  arrow ->SetFillColor(2);
  arrow ->SetLineColor(2);
  arrow ->SetLineWidth(1);
}
////////////////////////////////////////////////////////////////////////////
f1tdc_t0_rough::~f1tdc_t0_rough(){
}
////////////////////////////////////////////////////////////////////////////
void f1tdc_t0_rough::SetRoot(string ifname){
  add_tree(ifname);
  pack_tree();
  readtreeF1TDCL();
  readtreeF1TDCR();

}
////////////////////////////////////////////////////////////////////////////
void f1tdc_t0_rough::WriteParam(string ofname){
  param->WriteToFile(ofname.c_str());
}
////////////////////////////////////////////////////////////////////////////
void f1tdc_t0_rough::makehist(){

  for(int i=0;i<16;i++){
    h_Rs2t[i]     = new TH1F(Form("h_Rs2t%d",i)    , Form("h_Rs2t%d",i)     , 400,-15000,    0);
    h_Rs2b[i]     = new TH1F(Form("h_Rs2b%d",i)    , Form("h_Rs2b%d",i)     , 400,-15000,    0);
    h_Ls2t[i]     = new TH1F(Form("h_Ls2t%d",i)    , Form("h_Ls2t%d",i)     , 400,-20000,-5000);
    h_Ls2b[i]     = new TH1F(Form("h_Ls2b%d",i)    , Form("h_Ls2b%d",i)     , 400,-20000,-5000);
    set->SetTH1(h_Rs2t[i]      ,Form("F1TDC S2R%d top"              ,i),"TDC [ch]"    ,"counts");
    set->SetTH1(h_Rs2b[i]      ,Form("F1TDC S2R%d bottom"           ,i),"TDC [ch]"    ,"counts");
    set->SetTH1(h_Ls2t[i]      ,Form("F1TDC S2L%d top"              ,i),"TDC [ch]"    ,"counts");
    set->SetTH1(h_Ls2b[i]      ,Form("F1TDC S2L%d bottom"           ,i),"TDC [ch]"    ,"counts");
  }

  h_Rs0t     = new TH1F("h_Rs2t%d"    , "h_Rs2t%d"     , 400, -15000,    0);
  h_Rs0b     = new TH1F("h_Rs2b%d"    , "h_Rs2b%d"     , 400, -15000,    0);
  h_Ls0t     = new TH1F("h_Ls2t%d"    , "h_Ls2t%d"     , 400, -20000,-5000);
  h_Ls0b     = new TH1F("h_Ls2b%d"    , "h_Ls2b%d"     , 400, -20000,-5000);
  set->SetTH1(h_Rs0t      ,"F1TDC S0R top"     ,"TDC [ch]"    ,"counts");
  set->SetTH1(h_Rs0b      ,"F1TDC S0R bottom"  ,"TDC [ch]"    ,"counts");
  set->SetTH1(h_Ls0t      ,"F1TDC S0L top"     ,"TDC [ch]"    ,"counts");
  set->SetTH1(h_Ls0b      ,"F1TDC S0L bottom"  ,"TDC [ch]"    ,"counts");

}
////////////////////////////////////////////////////////////////////////////
void f1tdc_t0_rough::loop(){

  if( GetMaxEvent()>0 && GetMaxEvent()<ENum) ENum = GetMaxEvent();
  for(int n=0;n<ENum;n++){
    if(n%10000==0)cout<<n <<" / "<<ENum<<endl;
    convertF1TDCR(param);
    convertF1TDCL(param);
    tree->GetEntry(n);

    if(LF1Ref[0]>1.&&LF1Ref[1]>1.&&RF1Ref[0]>1.&&RF1Ref[1]>1.){
      //S0
      if(RS0T_F1TDC[0]>1. && DR_T5>1.)h_Rs0t ->Fill(RS0T_F1TDC[0]-RF1Ref[1]);
      if(RS0B_F1TDC[0]>1. && DR_T5>1.)h_Rs0b ->Fill(RS0B_F1TDC[0]-RF1Ref[1]);
      if(LS0T_F1TDC[0]>1. && DR_T5>1.)h_Ls0t ->Fill(LS0T_F1TDC[0]-LF1Ref[1]);
      if(LS0B_F1TDC[0]>1. && DR_T5>1.)h_Ls0b ->Fill(LS0B_F1TDC[0]-LF1Ref[1]);

      //S2
      for(int i=0;i<16;i++){
        if(RS2T_F1TDC[i]>1. && DR_T5>1.)h_Rs2t[i] ->Fill(RS2T_F1TDC[i]-RF1Ref[0]);
        if(RS2B_F1TDC[i]>1. && DR_T5>1.)h_Rs2b[i] ->Fill(RS2B_F1TDC[i]-RF1Ref[1]);
        if(LS2T_F1TDC[i]>1. && DR_T5>1.)h_Ls2t[i] ->Fill(LS2T_F1TDC[i]-LF1Ref[0]);
        if(LS2B_F1TDC[i]>1. && DR_T5>1.)h_Ls2b[i] ->Fill(LS2B_F1TDC[i]-LF1Ref[1]);
      }
    }
  }

}
////////////////////////////////////////////////////////////////////////////
void f1tdc_t0_rough::fit(){
  param->SetTdcOffset(CID_F1S0,0,1,0,h_Rs0t ->GetXaxis()->GetBinCenter(h_Rs0t ->GetMaximumBin()));
  param->SetTdcOffset(CID_F1S0,0,1,1,h_Rs0b ->GetXaxis()->GetBinCenter(h_Rs0b ->GetMaximumBin()));
  param->SetTdcOffset(CID_F1S0,0,0,0,h_Ls0t ->GetXaxis()->GetBinCenter(h_Ls0t ->GetMaximumBin()));
  param->SetTdcOffset(CID_F1S0,0,0,1,h_Ls0b ->GetXaxis()->GetBinCenter(h_Ls0b ->GetMaximumBin()));

  for(int i=0;i<16;i++){
    param->SetTdcOffset(CID_F1S2,i,1,0,h_Rs2t[i] ->GetXaxis()->GetBinCenter(h_Rs2t[i] ->GetMaximumBin()));
    param->SetTdcOffset(CID_F1S2,i,1,1,h_Rs2b[i] ->GetXaxis()->GetBinCenter(h_Rs2b[i] ->GetMaximumBin()));
    param->SetTdcOffset(CID_F1S2,i,0,0,h_Ls2t[i] ->GetXaxis()->GetBinCenter(h_Ls2t[i] ->GetMaximumBin()));
    param->SetTdcOffset(CID_F1S2,i,0,1,h_Ls2b[i] ->GetXaxis()->GetBinCenter(h_Ls2b[i] ->GetMaximumBin()));
  }

}
////////////////////////////////////////////////////////////////////////////
void f1tdc_t0_rough::draw(){

  c[0]->Clear();c[0]->Divide(2,2);
  c[0]->cd(1);h_Rs0t ->Draw();arrow ->DrawArrow(param->GetTdcOffset(CID_F1S0,0,1,0),1,param->GetTdcOffset(CID_F1S0,0,1,0),10,0.01,"<|");
  c[0]->cd(2);h_Rs0b ->Draw();arrow ->DrawArrow(param->GetTdcOffset(CID_F1S0,0,1,1),1,param->GetTdcOffset(CID_F1S0,0,1,1),10,0.01,"<|");

  c[0]->cd(3);h_Ls0t ->Draw();arrow ->DrawArrow(param->GetTdcOffset(CID_F1S0,0,0,0),1,param->GetTdcOffset(CID_F1S0,0,0,0),10,0.01,"<|");

  c[0]->cd(4);h_Ls0b ->Draw();arrow ->DrawArrow(param->GetTdcOffset(CID_F1S0,0,0,1),1,param->GetTdcOffset(CID_F1S0,0,0,1),10,0.01,"<|");


  c[1]->Clear();c[1]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[1]->cd(i+1);gPad->SetLogy(1);h_Rs2t[i]->Draw();
    arrow ->DrawArrow(param->GetTdcOffset(CID_F1S2,i,1,0),1,param->GetTdcOffset(CID_F1S2,i,1,0),10,0.01,"<|");
  }

  c[2]->Clear();c[2]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[2]->cd(i+1);gPad->SetLogy(1);h_Rs2b[i]->Draw();
    arrow ->DrawArrow(param->GetTdcOffset(CID_F1S2,i,1,1),1,param->GetTdcOffset(CID_F1S2,i,1,1),10,0.01,"<|");
  }

  c[3]->Clear();c[3]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[3]->cd(i+1);gPad->SetLogy(1);h_Ls2t[i]->Draw();
    arrow ->DrawArrow(param->GetTdcOffset(CID_F1S2,i,0,0),1,param->GetTdcOffset(CID_F1S2,i,0,0),10,0.01,"<|");
  }

  c[4]->Clear();c[4]->Divide(4,4);
  for(int i=0;i<16;i++){
    c[4]->cd(i+1);gPad->SetLogy(1);h_Ls2b[i]->Draw();
    arrow ->DrawArrow(param->GetTdcOffset(CID_F1S2,i,0,1),1,param->GetTdcOffset(CID_F1S2,i,0,1),10,0.01,"<|");
  }

  
}
////////////////////////////////////////////////////////////////////////////
void f1tdc_t0_rough::savecanvas(string ofname){
  c[0]->Print(Form("%s[",ofname.c_str()) );
  for(int i=0;i<NCanvas;i++){
    c[i]->Print(Form("%s" ,ofname.c_str()) );
  }
  c[NCanvas-1]->Print(Form("%s]",ofname.c_str()) );
  cout<<ofname<<" saved"<<endl;
}
////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ifname = "rootfiles/cosmic1020.root";
  string ofname = "toyamacro/pdf/f1tdc_t0_rough1020.pdf";
  string paramname = "param/default.param";
  int ch;
  int MaxNum = 0;
  bool output_flag = false;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hf:w:n:bp:"))!=-1){
    switch(ch){
    case 'f':
      ifname = optarg;
      cout<<"input filename : "<<ifname<<endl;
      break;
    case 'w':
      output_flag = true;
      ofname = optarg;
      cout<<"output pdf filename : "<<ofname<<endl;
      break;
    case 'n':
      MaxNum = atoi(optarg);
      break;
    case 'b':
      cout<<"BACH MODE!"<<endl;
      break;
    case 'p':
      paramname = optarg;
      cout<<"input param name : "<<paramname<<endl;
      break;
    case 'h':
      cout<<"-f : input root filename"<<endl;
      cout<<"-w : output pdf filename"<<endl;
      cout<<"-n : maximum number of analysed events"<<endl;
      cout<<"-p : output param file"<<endl;
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

  TApplication *theApp = new TApplication("App", &argc, argv);
  f1tdc_t0_rough *calib = new f1tdc_t0_rough();

  calib->SetMaxEvent(MaxNum);
  calib->SetRoot(ifname);
  calib->makehist();
  calib->loop();
  calib->fit();
  calib->draw();
  calib->WriteParam(paramname);
  if(output_flag)calib->savecanvas(ofname);
  delete calib;

  gSystem->Exit(1);
  theApp->Run();
  return 0;
}

