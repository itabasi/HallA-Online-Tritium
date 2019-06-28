/*
  TARG RunAction.hh
  5 Oct 2010
  T.Gogami
*/

#ifndef TARGRunAction_h
#define TARGRunAction_h 1

#include "TARGParamManager.hh"

#include "G4UserRunAction.hh"
#include "globals.hh"

#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TH1F.h>

class G4Run;

class TARGRunAction : public G4UserRunAction
{
public:
  TARGRunAction();
  virtual ~TARGRunAction();
  TARGRunAction(TARGParamManager*);
  
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);
  
private:
  TARGParamManager* paramMan;
  
  //File
private:
  G4String RfileName;
  TFile *file;
  TTree *Tree;
  
public:
  TTree* GetTree() { return Tree; };
  TFile* GetFile() { return file; };

};

#endif
