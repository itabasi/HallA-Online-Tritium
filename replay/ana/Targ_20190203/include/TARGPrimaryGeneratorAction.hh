/*
  TARGPrimaryGeneratorAction.cc
  
  1 Oct 2010
  T.Gogami
*/

#ifndef TARGPrimaryGeneratorAction_h
#define TARGPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "TARGParamManager.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4String.hh"

class G4ParticleGun;
class G4Event;

class TARGPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  TARGPrimaryGeneratorAction();
    ~TARGPrimaryGeneratorAction();
  TARGPrimaryGeneratorAction(TARGParamManager*);
  
  void GeneratePrimaries( G4Event* anEvent );
  void GenPosRandinTarget( G4Event* anEvent );
  
private:
  G4ParticleGun* particleGun;
  TARGParamManager *paramMan;
  G4int GenFlag;
  G4int SeedFlag;
  G4String SeedFile;
  FILE* sfile;
  //char sfile[500];
  G4double rasterX,rasterY,rasterZ;
  G4double bpc,bpw;  // Beam momentum (cetner,width)
  G4double bthc,bthw;// Beam theta    (cetner,width)
  G4double bphc,bphw;// Beam phi      (cetner,width)
  //G4double BeamOffx,BeamOffy,BeamOffz;
  G4ThreeVector BeamOffset;
  G4double Beamx,Beamy,Beamz;
  G4double BeamMom;
  G4int ParticleFlag;
  G4ParticleDefinition* particle;
  void GenParPoint(G4Event* );
  void GenParUni(G4Event* );
  void SetMom(G4ThreeVector , G4ThreeVector , G4Event* , G4int );
  //void GenParSeed(G4Event* );

};
#endif
