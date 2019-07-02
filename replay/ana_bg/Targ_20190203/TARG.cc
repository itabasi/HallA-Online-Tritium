/*
  TARG simulation
  1 Oct 2010
  T.Gogami
*/

#include "TARGDetectorConstruction.hh"
#include "TARGPhysicsList.hh"
#include "TARGPrimaryGeneratorAction.hh"
#include "TARGSteppingAction.hh"
#include "TARGRunAction.hh"
#include "TARGEventAction.hh"
#include "TARGvdHit.hh"
#include "TARGvdSD.hh"
#include "TARGParamManager.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"

#include "GetNumberFromKernelEntropyPool.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif



int main(int argc,char** argv)
{
  
  //  ====== Choose the Random Engine ========= //
  int initSeed=GetIntFromKernelEntropyPool()&0x7FFFFFFF;
  CLHEP::HepRandom::setTheSeed(initSeed);
  int startSeed=CLHEP::HepRandom::getTheSeed();
  G4cout << "*** Initial Seed = " << startSeed << G4endl;
  CLHEP::HepRandom::showEngineStatus();
  
  
  G4String macFile   = "./macro/run.mac" ;//Default Mac file
  G4String inputFile = "./input/param.in";//Default input file
  int rfnum = 0;
  //G4int rfnum = 0;
  if(argc == 1){
    macFile   = "./macro/run.mac" ;//Default Mac file  
    inputFile = "./input/param.in";//Default input file
  }
  else if(argc == 2){
    macFile = argv[1];
  }
  else if(argc == 3){
    macFile   = argv[1];
    inputFile = argv[2];
    //rfnum = atoi(argv[3]);
  }
  else if(argc == 4){
    macFile   = argv[1];
    inputFile = argv[2];
    rfnum = atoi(argv[3]);
  }
  else {
    macFile   = argv[1];
    inputFile = argv[2];
    rfnum     = atoi(argv[3]);
  }

  //Parameter Manager
  TARGParamManager * paramMan = new TARGParamManager(inputFile);
  if(rfnum!=0)paramMan->SetRFnum(rfnum);
  
  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // User Initialization classes (mandatory)
  TARGDetectorConstruction* detector = new TARGDetectorConstruction(paramMan);
  runManager->SetUserInitialization(detector);
  //G4VUserPhysicsList* physics = new TARGPhysicsList;
  G4VUserPhysicsList* physics = new TARGPhysicsList(paramMan);
  runManager->SetUserInitialization(physics);
  
  // User Action classes
  G4VUserPrimaryGeneratorAction* gen_action = new TARGPrimaryGeneratorAction(paramMan);
  runManager->SetUserAction(gen_action);
  G4UserRunAction* run_action = new TARGRunAction(paramMan);
  runManager->SetUserAction(run_action);
  G4UserEventAction* event_action 
    = new TARGEventAction(paramMan);
  runManager->SetUserAction(event_action);
  G4UserSteppingAction* stepping_action = new TARGSteppingAction();
  runManager->SetUserAction(stepping_action);

  // Initialize G4 kernel
  runManager->Initialize();
  
  // Get the pointer to the User Interface manager
  G4UImanager * UI = G4UImanager::GetUIpointer();  

#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif    
  
  G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
  session = new G4UIterminal(new G4UItcsh);      
#else
  session = new G4UIterminal();
#endif
  G4String command = "/control/execute " + macFile;
  //UI->ApplyCommand("/control/execute vis.mac");     
  //UI->ApplyCommand("/control/execute vrml.mac");     
  //UI->ApplyCommand("/control/execute run.mac");     
  UI->ApplyCommand( command );
  //session->SessionStart();
  delete session;
  
#ifdef G4VIS_USE
  delete visManager;
#endif     
  
  //cout << endl;
  G4cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << G4endl;
  G4cout << " Macro file        :  " << macFile   << G4endl;
  G4cout << " Input file        :  " << inputFile << G4endl;
  G4cout << " Output ROOT file  :  " << paramMan->GetROOTFileName() 
	 << " - " << rfnum << G4endl;
  G4cout << " Beam Momentum     :  " << paramMan->GetBeamEnergy()   
	 << " GeV/c " << G4endl;
  G4cout << " Target            :  " << paramMan->GetTarget()
	 << " " << paramMan->GetTThickness() << " cm" << G4endl;
  //G4cout << " POSI ANALYSYS FLAG:  " << paramMan->GetPosiAnalFlag() << G4endl;
  G4cout << " EM Flag           :  " << paramMan->GetEMFlag()       << G4endl;
  G4cout << " Decay Flag        :  " << paramMan->GetDecayFlag()    << G4endl;
  G4cout << " Hadron Flag       :  " << paramMan->GetHadronFlag()   << G4endl;
  G4cout << " Generation Flag   :  " << paramMan->GetPGenFlag()     << G4endl;
  G4cout << " Raster[cm]        :  " 
	 << paramMan->GetRasterX() << " , " 
	 << paramMan->GetRasterY() << " , " 
	 << paramMan->GetRasterZ() << G4endl;
  G4cout << " Beam Particle     :  " << paramMan->GetParticleFlag() << G4endl;
  G4cout << " ( 0:e- , 1:e+ , 2:pi- , 3:pi+ , 4:mu-, 5:p, 6:kaon+ )" << G4endl;
  //G4cout << "File Number       :  " << rfnum <<G4endl;
  G4cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << G4endl;
  
  delete paramMan;
  delete runManager;
  
  return 0;
}


