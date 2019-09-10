/*
  TARG ParamManager.hh
  8 Oct 2010
  T.Gogami
*/

#ifndef TARGParamManager_h
#define TARGParamManager_h 1

//#include<TFile.h>
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include <fstream>
using namespace std;

class TARGParamManager
{
public:
  TARGParamManager(G4String fname);
  ~TARGParamManager();
  
public:  
  G4bool readparam();
private:
  G4String Filename;    //Input file name
  G4String ROOTFilename;//Output ROOT file name
  G4int RFnum;          //Number of Output ROOT file name
  G4double BeamEnergy;  //Incident beam momentum (input)
  G4double BeamP_w;     //Incident beam momentum width (input)
  G4double BeamTheta;   //Incident beam theta (input)
  G4double BeamTheta_w; //Incident beam theta width (input)
  G4double BeamPhi;     //Incident beam phi (input)
  G4double BeamPhi_w;   //Incident beam phi width (input)
  G4double GenP;        // Generated particle momentum
  G4double GenPx;// Generated particle momentum
  G4double GenPy;// Generated particle momentum
  G4double GenPz;// Generated particle momentum
  G4double GenTheta;    // Generated particle theta
  G4double GenPhi;      // Generated particle phi
  G4int PosiAnalFlag;   //POSITRON analysis flag
  G4String Target;      //Target
  G4String SeedFile;    //Seed file
  G4double TThickness;  //Target thickness [cm]
  G4double minMom;      //Mimimum momentum of analyzed particle
  G4double maxMom;      //Maximum momentum of analyzed particle
  G4double mintheta;    //Mimimum theta of analyzed particle
  G4double maxtheta;    //Maximum theta of analyzed particle
  G4int EMFlag;         //ElectroMagnetic Flag
  G4int DecayFlag;      //Decay Flag
  G4int HadronFlag;     //Hadron Flag
  G4int pGenerationFlag;//Particle Generation Flag
  G4double rasterx,rastery,rasterz;   //raster[cm]
  G4double beamoffx,beamoffy,beamoffz;//Beam Position Offset [cm]
  G4double beamx , beamy , beamz; 
  G4int ParticleFlag;   // Particle Flag
  G4int SeedFlag;       // Seed Flag
  G4int evid;        // Event ID 
  G4int evnum;       // Number of event

public:
  G4String GetFileName()     { return Filename;       };
  G4String GetROOTFileName() { return ROOTFilename;   };
  G4int    GetRFnum()        { return RFnum;     };
  //G4int    SetRFnum(G4int rfn){ RFnum=rfn;       };
  void     SetRFnum(G4int rfn){ RFnum=rfn;       };
  void     SetBeamPos(G4ThreeVector bpos){
    beamx = bpos.x();
    beamy = bpos.y();
    beamz = bpos.z();
  }
  void SetGenP(G4double bp)        { GenP = bp; }
  void SetGenPvec(G4double a , G4double b, G4double c){ 
    GenPx = a;
    GenPy = b;
    GenPz = c;
  }
  void SetGenTheta(G4double btheta){ GenTheta = btheta; }
  void SetGenPhi(G4double bphi)    { GenPhi = bphi; }
  G4double GetGenP()         { return GenP; }
  G4ThreeVector GetGenPvec() { return G4ThreeVector(GenPx,GenPy,GenPz); }
  G4double GetGenTheta()     { return GenTheta; }
  G4double GetGenPhi()       { return GenPhi; }
  G4int GetEVID()            { return evid; }
  G4int SetEVID(G4int eventid){ evid = eventid; }
  
  G4ThreeVector GetBeamPos() { return G4ThreeVector(beamx,beamy,beamz);}
  G4double GetBeamEnergy()   { return BeamEnergy;     };//momentum
  G4double GetBeamP_w()      { return BeamP_w;        };//momentum
  G4double GetBeamTheta()    { return BeamTheta;      };//Theta 
  G4double GetBeamTheta_w()  { return BeamTheta_w;    };//Theta
  G4double GetBeamPhi()      { return BeamPhi;        };//Phi
  G4double GetBeamPhi_w()    { return BeamPhi_w;      };//Phi
  G4int    GetPosiAnalFlag() { return PosiAnalFlag;   };
  G4String GetTarget()       { return Target;         };
  G4double GetTThickness()   { return TThickness;     };
  G4double GetMinMom()       { return minMom;         };
  G4double GetMaxMom()       { return maxMom;         };
  G4double GetMinTheta()     { return mintheta;       };
  G4double GetMaxTheta()     { return maxtheta;       };
  G4int    GetEMFlag()       { return EMFlag;         };
  G4int    GetDecayFlag()    { return DecayFlag;      };
  G4int    GetHadronFlag()   { return HadronFlag;     };
  G4int    GetPGenFlag()     { return pGenerationFlag;};
  G4double GetRasterX()      { return rasterx;        };
  G4double GetRasterY()      { return rastery;        };
  G4double GetRasterZ()      { return rasterz;        };
  G4String GetSeedFile()     { return SeedFile;       };
  G4int    GetSeedFlag()     { return SeedFlag;       };
  
  G4ThreeVector GetBeamOffset()
  { return G4ThreeVector(beamoffx,beamoffy,beamoffz);};
  G4int GetParticleFlag()    { return ParticleFlag;   };
  
};

#endif
