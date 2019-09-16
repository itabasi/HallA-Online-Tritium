/*
  "TARGPhysicsList.hh"
  
  Toshi Gogami , 15July2015
*/

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
//#include "globals.hh"
#include "G4Decay.hh"
#include "TARGParamManager.hh"


class TARGPhysicsList: public G4VUserPhysicsList
{
public:
  TARGPhysicsList();
  TARGPhysicsList(TARGParamManager*);
  ~TARGPhysicsList();
  
protected:
  // Construct particle and physics
  void ConstructParticle();
  void ConstructProcess();
  void SetCuts();
  
protected:
  // these methods Construct particles 
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBaryons();
  void ConstructHeavyIon();
  void ConstructOp();
  
private:
  // these methods Construct physics processes and register them
  void AddTransportation();
  void ConstructDecay();
  void ConstructEM();
  void ConstructHadronic();
  void SetNbOfPhotonsCerenkov(G4int);
  
private:
  TARGParamManager* paramMan;
  G4int emflag;
  G4int decayflag;
  G4int hadronflag;
  G4double defaultCutValue;
  G4int opflag;
  
};

#endif

 
