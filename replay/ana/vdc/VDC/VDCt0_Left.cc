#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <sstream>
using namespace std;
#include <TApplication.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>
#include <TCut.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLatex.h>
#include <TText.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TColor.h>
#include <TPaveText.h>
#include <TRandom.h>
#include <TMinuit.h>
#include <TLegend.h>
#include "GmpVDCt0Calib.h"

string ofname("output.pdf");
string ofroot("output.root");

//===================================================================//
//============================= Main ================================//
//==================================================================//

int main(int argc, char** argv){


 TApplication *theApp =new TApplication("App",&argc,argv);

    gStyle->SetTitleX(0.5);
    gStyle->SetTitleY(1.05);
    gStyle->SetTitleH(0.20);
    gStyle->SetTitleW(0.50);
    gStyle->SetLabelSize(0.04,"X");
    gStyle->SetLabelSize(0.04,"Y");
    gStyle->SetTitleSize(0.04,"X");
    gStyle->SetTitleSize(0.04,"Y");

    //const TString ROOTfilePath = "/work/halla/gmp12/longwu/gmp_analysis/rootfiles";
    // const TString ROOTfilePath = "/work/halla/gmp12/longwu/gmp_analysis/rootfiles/Fall2016/K3-8";
   // const TString ROOTfilePath = "/chafs1/work1/tritium/rootfiles_1.6_root6";
    
    //    const TString ROOTfilePath = "/chafs1/work1/tritium/Rootfiles";
    const TString ROOTfilePath = "/data1/root";    
    TChain* T = new TChain("T");
    int run=111200;
    std::ostringstream str;
    // str << ROOTfilePath << "/gmp_" << run;
    str << ROOTfilePath << "/tritium_" << run;

    TString basename = str.str().c_str();
    TString rootfile = basename + ".root";
    cout<<rootfile<<endl;
    Long_t i=0;
    while ( !gSystem->AccessPathName(rootfile.Data()) ) {
        T->Add(rootfile.Data());
        cout << "ROOT file " << rootfile << " added to TChain." << endl;
        i++;
        //rootfile = basename + Form("_%d",i) + ".root";
        rootfile = basename + "_" + i + ".root";
    }

    if (i==0) {
        cerr << "The specified run does not exist or has not been replayed." << endl;
        return 0;
    }

    //gROOT->LoadMacro("GmpVDCt0Calib.C+");

    TString plotname = Form("L_vdct0_%d.pdf",run);
    GmpVDCt0Calib* vdct0 = new GmpVDCt0Calib;
    //    vdct0->Calibrate(T,1,"(DL.evtypebits>>2)&1",plotname);
    vdct0->Calibrate(T,1,"DL.evtypebits==14",plotname);
    vdct0->SaveDatabase(Form("./db_L.vdc.%d.dat",run));

    Double_t min = 10000, max = 0;
    for (UInt_t ii=0; ii<GmpVDCt0Calib::GetNumberOfPlanes(); ii++) {
        for (UInt_t jj=0; jj<GmpVDCt0Calib::GetNumberOfWiresPerPlane(); jj++) {
            Double_t offset = vdct0->GetOffset(ii,jj);
            if (offset<min) min=offset;
            if (offset>max) max=offset;
        }
    }

    TCanvas* c1 = new TCanvas;
    //TH1F* hFrame = c1->DrawFrame(0,1480,GmpVDCt0Calib::GetNumberOfWiresPerPlane(),1560,"VDC t0 calibration reuslt ( (DL.evtypebits&1<<3) == 1<<3 );Wire number;Time offset (TDC channel)");

    TH1F* hFrame = c1->DrawFrame(0,0,GmpVDCt0Calib::GetNumberOfWiresPerPlane(),3000,"VDC t0 calibration result;Wire number;Time offset (TDC channel)");
    hFrame->GetXaxis()->CenterTitle();
    hFrame->GetYaxis()->CenterTitle();
    hFrame->GetYaxis()->SetTitleOffset(1.05);

    TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    for (UInt_t ii=0; ii<GmpVDCt0Calib::GetNumberOfPlanes(); ii++) {
        TGraph* graph = new TGraph;
        graph->SetName(vdct0->GetPlaneName(ii));
        UInt_t npoint = 0;
        for (UInt_t j=0; j<GmpVDCt0Calib::GetNumberOfWiresPerPlane(); j++) {
            Double_t offset = vdct0->GetOffset(ii,j);
            if (offset>0) graph->SetPoint(npoint++,j,offset);
        }
        graph->SetMarkerColor(ii+1);
        graph->SetMarkerStyle(ii+20);
        graph->Draw("P");
        leg->AddEntry(vdct0->GetPlaneName(ii),vdct0->GetPlaneName(ii),"p");
    }

    leg->Draw();
    c1->Update();

    c1->SaveAs(Form("plot_L/vdct0plot_%d.pdf",run));
 
  gSystem->Exit(1);
  theApp->Run();



  return 0;

}//end main




