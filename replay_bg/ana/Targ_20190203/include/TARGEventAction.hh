/*
  TARGEventAction.hh
  1 Oct 2010
  T.Gogami
*/

#ifndef TARGEventAction_h
#define TARGEventAction_h 1

#include "G4UserEventAction.hh"
#include <TTree.h>
#include <TFile.h>
#include <TObjArray.h>
#include "TARGParamManager.hh"
#include <time.h>
class G4Event;

class TARGEventAction : public G4UserEventAction
{
public:
  TARGEventAction();
  ~TARGEventAction();
  TARGEventAction(TARGParamManager*);
  
public:
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
private:
  G4int primaryEventID;
  TARGParamManager * paramMan;
  time_t start,end;
  
  //ROOT FILE
public:
  void definetree(const G4Event*);
  void filldata(const G4Event*);
private:
  const int maxhit = 20;
  G4int nhit;
  TTree* tree;      //Tree
  TFile* file;      //ROOT file
  G4int evID;       //event ID
  G4int evID2;      //event ID
  G4int trackID;    //track ID
  G4int SeedFlag;
  G4double beamx,beamy,beamz; // Beam x,y,z position 
  G4double beamp;             // Beam Momentum
  G4double beampx;            // Beam Momentum
  G4double beampy;            // Beam Momentum
  G4double beampz;            // Beam Momentum
  G4double beam_theta;        // Beam Momentum
  G4double beam_phi;          // Beam Momentum
  //G4double x,y,z;   //position on virtual detector
  G4double x[20],y[20],z[20];   //position on virtual detector
  //G4double px,py,pz;//momentum on virtual detector
  //G4double p;       //momentum on virtual detector
  G4double px[20],py[20],pz[20];//momentum on virtual detector
  G4double p[20];       //momentum on virtual detector

  // ---- LHRS ------ //
  G4double x2[20],y2[20],z2[20];   //position on virtual detector
  G4double px2[20],py2[20],pz2[20];//momentum on virtual detector
  G4double p2[20];       //momentum on virtual detector
  G4int nhit2;
  G4double theta2[20];   //theta on virtual detector
  G4double phi2[20];     //phi on virtual detector
  G4double charge2[20];     //phi on virtual detector
  G4double pid2[20];     // Particle ID
  G4int ntrig2;
  
  // ---- RHRS ------ //
  G4double x3[20],y3[20],z3[20];   //position on virtual detector
  G4double px3[20],py3[20],pz3[20];//momentum on virtual detector
  G4double p3[20];       //momentum on virtual detector
  G4int nhit3;
  G4double theta3[20];   //theta on virtual detector
  G4double phi3[20];     //phi on virtual detector
  G4double charge3[20];     //phi on virtual detector
  G4double pid3[20];     // Particle ID
  G4int ntrig3;

  // ---- Beam ------ //
  G4double x4[20],y4[20],z4[20];   //position on virtual detector
  G4double px4[20],py4[20],pz4[20];//momentum on virtual detector
  G4double p4[20];       //momentum on virtual detector
  G4int nhit4;
  G4double theta4[20];   //theta on virtual detector
  G4double phi4[20];     //phi on virtual detector
  G4double charge4[20];     //phi on virtual detector
  G4double pid4[20];     // Particle ID
  G4int ntrig4;
  
  //G4double theta;   //theta on virtual detector
  //G4double phi;     //phi on virtual detector
  G4double theta[20];   //theta on virtual detector
  G4double phi[20];     //phi on virtual detector
  //nG4double charge;  // charge of a detected particle
  G4double charge[20];  // charge of a detected particle
  G4String pname;   // particle name
  G4int pioflag;    // pion flag
  G4int kaoflag;    // kaon flag
  G4int proflag;    // proton flag
  G4int eflag;      // ?
  G4int phoflag;    // photon flag
  G4int posflag;    // positron flag
  //G4int eleflag[20];    // electron flag
  G4int eleflag;    // electron flag
  G4int muoflag;    // muon flag
  G4int PosiAnalFlag;
  G4int ProcID;     //Process ID
  G4int passth;     //beam passed through the target
  G4int eIoni;      //e Ionization
  G4int eBrem;      //e Brems
  G4int compt;      //compton scatterting
  G4int conv;       //photon conversion
  G4int phoele;     //photoelectric effect
  G4double minMom;  //min Momentum
  G4double maxMom;  //max Momentum
  G4double minTheta;//min Theta
  G4double maxTheta;//max Theta
  G4int btrig;      // trigger flag (test , 12Dec2013)

};
#endif
