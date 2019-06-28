/*
  TARGvdSD.cc
  5 Oct 2010
  T.Gogami
*/

#include "TARGvdSD.hh"
#include "TARGvdHit.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"

TARGvdSD::TARGvdSD( G4String name,
		    G4int nCells,
		    G4String colName )
  : G4VSensitiveDetector(name),
    numberOfCells(nCells),//interger
    HCID(-1),//integer
    procID(0)//Process ID
{
  G4String HCname;
  collectionName.insert(HCname=colName);
  CellID = new G4int[numberOfCells];//interger
}

TARGvdSD::~TARGvdSD()
{
  delete [] CellID;//interger
}

void TARGvdSD::Initialize(G4HCofThisEvent*)
{
  VDCollection = new TARGvdHitsCollection
    (SensitiveDetectorName,collectionName[0]); 
  for(G4int j=0;j<numberOfCells;j++)
  {
    CellID[j] = -1;//integer
  }
}

G4bool TARGvdSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  //========TouchableHistory=======///
  G4TouchableHistory* hist = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  const G4VPhysicalVolume* physVol = hist->GetVolume();
  G4int copyID = hist->GetReplicaNumber();

  //=========Get Track ===========///
  G4Track *aTrack = aStep->GetTrack();
  
  //==========Get Process=========///
  const G4VProcess*  vpro = aTrack -> GetCreatorProcess() ;//Get creator process
  const G4String& proname = vpro->GetProcessName(); //Process name
  if(vpro!=0){
    //G4cout << " " << proname << G4endl;
    if (proname=="eIoni")procID = 1 ;   //electron ionization
    else if(proname=="eBrem")procID = 2;//electron brems
    else if(proname=="compt")procID = 3;//compton scattering
    else if(proname=="phot" )procID = 4;//photoelectric effect
    else if(proname=="conv" )procID = 5;//photon conversion
  }
  else procID = 0;
  //G4cout << procID <<G4endl;
  //if(procID!=0 && procID!=2 && procID!=1){
  //  G4cout << procID << ":::: " <<proname <<G4endl;
  //}
  
  //==========Get Particle info========///
  G4ParticleDefinition *particle =
    aTrack -> GetDynamicParticle() -> GetDefinition();
  
  //==========Virtual Detector hit====///
  TARGvdHit* vdHit =
    new TARGvdHit(physVol->GetLogicalVolume());
  vdHit->SetEdep( edep );
  G4AffineTransform aTrans = hist->GetHistory()->GetTopTransform();
  aTrans.Invert();
  vdHit->SetPos( aStep->GetPreStepPoint()->GetPosition() );
  vdHit->SetMom( aTrack->GetMomentum() );
  vdHit->SetPname( particle->GetParticleName() );
  vdHit->SetCharge( particle->GetPDGCharge() );
  vdHit->SetProcID( procID );
  G4int icell = VDCollection->insert( vdHit );
  
  CellID[copyID] = icell - 1;
//    if(verboseLevel>0)
//      { G4cout << " New vd Hit on CellID " << copyID << G4endl; }

  return true;
}

void TARGvdSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, VDCollection );
}

void TARGvdSD::clear()
{
} 

void TARGvdSD::DrawAll()
{
} 

void TARGvdSD::PrintAll()
{
} 
