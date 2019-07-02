/*
  TARGvdSD.hh
  5 Oct 2010
  T.Gogami
*/

#ifndef TARGvdSD_h
#define TARGvdSD_h 1

#include "TARGvdHit.hh"

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"

class TARGvdSD : public G4VSensitiveDetector
{
  
public:
  TARGvdSD(G4String name, G4int nCells, G4String colName);
  ~TARGvdSD();
  
  void Initialize(G4HCofThisEvent*HCE);
  G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
  void EndOfEvent(G4HCofThisEvent*HCE);
  void clear();
  void DrawAll();
  void PrintAll();
  
private:
  TARGvdHitsCollection *VDCollection;
  
  int* CellID;
  int numberOfCells;
  int HCID;
  int procID;
};




#endif

