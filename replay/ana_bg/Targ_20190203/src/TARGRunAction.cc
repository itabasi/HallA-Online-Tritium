/*
  TARG RunAction.cc
  5 Oct 2010
  T.Gogami
*/

#include "TARGRunAction.hh"
#include "TARGParamManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TH1F.h>
#include <stdlib.h>

TARGRunAction::TARGRunAction()
{}

TARGRunAction::~TARGRunAction()
{}

TARGRunAction::TARGRunAction(TARGParamManager* PManDragon)
  :paramMan(PManDragon),RfileName("DRAGON.root")
{
  RfileName = paramMan->GetROOTFileName();
  if( paramMan->GetRFnum()!=0 ){
    char tmp[80];
    sprintf(tmp,"%d",paramMan->GetRFnum());
    RfileName = RfileName+"_"+tmp;

  }
}

void TARGRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << " " << G4endl;
  G4cout << "// Simulation Start //" << G4endl;
  G4cout << " " << G4endl;
  file = new TFile(RfileName,"recreate"); //Open ROOT file (Recreate)
  //file = new TFile("../analysis/root/test.root","recreate");
  TTree *tree = new TTree("tree"," TARGET Simulation ");//Tree
  Tree = (TTree*)file->Get("tree");
}

void TARGRunAction::EndOfRunAction(const G4Run* aRun)
{
  Tree->Write();  //Write tree into ROOT file
  file -> Close();//Close ROOT file
  G4cout << "  END... bye ! " << G4endl;
  
}
