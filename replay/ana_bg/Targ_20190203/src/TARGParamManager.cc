/*
  TARG ParamManger.cc
  8 Oct 2010
  T.Gogami
*/

#include "TARGParamManager.hh"
//#include "fstream"
#include <fstream>
using namespace std;

TARGParamManager::TARGParamManager(G4String inputfilename)
  :Filename(inputfilename),ROOTFilename("DRAGON.root")
  ,RFnum(0)
  ,BeamEnergy(0),PosiAnalFlag(0),Target("DRAGON")
  ,minMom(0),maxMom(0),mintheta(0),maxtheta(0), ParticleFlag(0)
{
  readparam();
}
TARGParamManager::~TARGParamManager()
{}

G4bool TARGParamManager::readparam()
{
  G4String dummy;
  G4String Logfile;
  ifstream *file  = new ifstream(Filename);
  *file >> dummy >> ROOTFilename;                //ROOT File name
  Logfile = ROOTFilename + "_Log";               //Logfile name
  ofstream *ofile = new ofstream(Logfile);       //Open Logfile
  // ~~~ (1) ~~~~~
  *ofile<< dummy << " "<< ROOTFilename << G4endl;//ROOT File name
  // ~~~ (2) ~~~~~
  *file >> dummy >> SeedFile >> SeedFlag;        // Seed file
  *ofile<< dummy 
	<< " " << SeedFile 
	<< " " << SeedFlag << G4endl;            // Seed file
  // ~~~ (3) ~~~~~
  *file >> dummy >> BeamEnergy >> BeamP_w;       //Incident Beam Momentum [GeV]
  *ofile<< dummy 
	<< " " << BeamEnergy
	<< " " << BeamP_w << G4endl;             //Incident Beam Momentum [GeV]
  // ~~~ (4) ~~~~~
  *file >> dummy >> BeamTheta >> BeamTheta_w;    // Theta [rad]
  *ofile<< dummy 
	<< " " << BeamTheta 
	<< " " << BeamTheta_w << G4endl;         // Theta [rad]
  // ~~~ (5) ~~~~~
  *file >> dummy >> BeamPhi >> BeamPhi_w;        // Phi [rad]
  *ofile<< dummy 
	<< " " << BeamPhi 
	<< " " << BeamPhi_w << G4endl;           // Phi [rad]
  // ~~~ (6) ~~~~~
  *file >> dummy;   *file >> Target >> TThickness;  //Target thickness [cm]
  *ofile<< dummy 
	<< " " << Target 
	<< " " << TThickness << " cm" << G4endl;    //Target and Target thickness
  // ~~~ (7) ~~~~~
  *file >> dummy;   *file >> minMom;             //Minimum momentum of analyzed particle
  *ofile<< dummy << " "<< minMom << G4endl;      //Minimum momentum of analyzed particle
  // ~~~ (8) ~~~~~
  *file >> dummy;   *file >> maxMom;             //Maximum momentum of analyzed particle
  *ofile<< dummy << " "<< maxMom << G4endl;      //Maximum momentum of analyzed particle
  // ~~~ (9) ~~~~~
  *file >> dummy;   *file >> mintheta;           //Minimum theta of analyzed particle
  *ofile<< dummy << " "<< mintheta << G4endl;    //Minimum theta of analyzed particle
  // ~~~ (10) ~~~~~
  *file >> dummy;   *file >> maxtheta;           //Maximum theta of analyzed particle
  *ofile<< dummy << " "<< maxtheta << G4endl;    //Maximum theta of analyzed particle
  // ~~~ (11) ~~~~~
  *file >> dummy;   *file >> PosiAnalFlag;       //Positron Analysis flag
  *ofile<< dummy << " "<< PosiAnalFlag << G4endl;//Positron Analysis flag
  // ~~~ (12) ~~~~~
  *file >> dummy;   *file >> EMFlag;             //ElectroMagnetic flag
  *ofile<< dummy << " "<<  EMFlag << G4endl;     //ElectroMagnetic flag
  // ~~~ (13) ~~~~~
  *file >> dummy;   *file >> DecayFlag;          //Decay flag
  *ofile<< dummy << " "<< DecayFlag << G4endl;   //Decay flag
  // ~~~ (14) ~~~~~
  *file >> dummy;   *file >> HadronFlag;         //hadron flag
  *ofile<< dummy << " "<< HadronFlag << G4endl;  //hadron flag
  // ~~~ (15) ~~~~~
  *file >> dummy;   *file >> pGenerationFlag;         //Particle Generation Flag
  *ofile<< dummy << " "<< pGenerationFlag << G4endl;  //Particle Generation Flag
  // ~~~ (16) ~~~~~
  *file >> dummy;// raster
  *file >> rasterx; *file >> rastery ; *file >> rasterz; // raster x,y,z
  *ofile << dummy << " " 
	 << rasterx << "," 
	 << rastery << ","
	 << rasterz << endl;
  // ~~~ (17) ~~~~~
  *file >> dummy;
  *file >> beamoffx; *file >> beamoffy ; *file >> beamoffz; // Beam position Offset [cm]
  *ofile << dummy << " " 
	 << beamoffx << "," 
	 << beamoffy << ","
	 << beamoffz << G4endl;
  // ~~~ (18) ~~~~~
  *file >> dummy;  *file >> ParticleFlag;           //Particle Flag
  *ofile<< dummy << " " << ParticleFlag << G4endl;  //Particle Flag
  
  //for(int i=0 ; i<1000000;i++){
  //  G4cout << ROOTFilename << G4endl;
  //  G4cout << BeamEnergy << G4endl;
  //G4cout << "Logfile =====" << Logfile << G4endl;
  //}
  
  ofile -> close();
  file  -> close();
  
  return true;
}
