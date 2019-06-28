/*
  "TARGPrimaryGeneratorAction.cc"
  
  T.Gogami, November 23, 2018 
*/

#include "TARGPrimaryGeneratorAction.hh"
#include "TARGParamManager.hh"
#include "G4ParticleTable.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ThreeVector.hh"
#include "G4Geantino.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4String.hh"
#include <TMath.h>

class G4ParticleGun;
class G4Event;


TARGPrimaryGeneratorAction::TARGPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  particleGun -> SetParticleDefinition( G4Geantino::GeantinoDefinition() );
  particleGun -> SetParticleEnergy( 2.344*CLHEP::GeV );
  particleGun -> SetParticlePosition( G4ThreeVector(0*CLHEP::m , 0*CLHEP::m , -2.0*CLHEP::cm) );
}
TARGPrimaryGeneratorAction::TARGPrimaryGeneratorAction(TARGParamManager* PManDragon)
  :paramMan(PManDragon),
   Beamx(0.0),Beamy(0.0),Beamz(0.0),
   BeamOffset( G4ThreeVector( 0.0 , 0.0 , 0.0 ) ),
   GenFlag(0),rasterX(0.0),rasterY(0.0),rasterZ(0.0),
   ParticleFlag(0),
   bpc(0.0),bpw(0.0),
   bthc(0.0),bthw(0.0),
   bphc(0.0),bphw(0.0),
   SeedFlag(0)
{
  //G4double beamenergy = paramMan->GetBeamEnergy();//Momentum
  BeamMom = paramMan->GetBeamEnergy();//Momentum
  GenFlag = paramMan->GetPGenFlag();
  rasterX = paramMan->GetRasterX();
  rasterY = paramMan->GetRasterY();
  rasterZ = paramMan->GetRasterZ();
  BeamOffset = paramMan->GetBeamOffset();
  ParticleFlag = paramMan->GetParticleFlag();
  bpc  = BeamMom*CLHEP::GeV;
  bpw  = paramMan->GetBeamP_w()*CLHEP::GeV;
  bthc = paramMan->GetBeamTheta();
  bthw = paramMan->GetBeamTheta_w();
  bphc = paramMan->GetBeamPhi();
  bphw = paramMan->GetBeamPhi_w();
  SeedFlag = paramMan->GetSeedFlag();
  //SeedFile = paramMan->GetSeedFile();
  sfile=fopen( paramMan->GetSeedFile(), "r");
  
  //for(int i=0 ; i<10000 ; i++){
  //  cout << GenFlag << endl;
  //}
  //particleGun = new G4ParticleGun(n_particle);
  //particleGun = new G4ParticleGun();
  //particleGun -> SetParticleDefinition( G4Geantino::GeantinoDefinition() );
  G4ParticleTable* particleTable  = G4ParticleTable::GetParticleTable();
  //particleGun -> SetParticleEnergy( beamenergy * CLHEP::GeV );
  //particleGun -> SetParticleMomentum( BeamMom * CLHEP::GeV );//Momentum
  //particleGun -> SetParticlePosition( G4ThreeVector(0*m , 0*m , -2.0*cm) );
  //for(int i=0 ; i<900000 ; i++){
  //    G4cout << ParticleFlag  << G4endl;
  //}
  if(ParticleFlag==0){
    particle = particleTable->FindParticle("e-");
  }
  else if(ParticleFlag==1){
    particle = particleTable->FindParticle("e+");
  }
  else if(ParticleFlag==2){ 
    particle = particleTable->FindParticle("pi-");
  }
  else if(ParticleFlag==3){
    particle = particleTable->FindParticle("pi+");
  }
  else if(ParticleFlag==4){
    particle = particleTable->FindParticle("mu-");
  }
  else if(ParticleFlag==5){
    particle = particleTable->FindParticle("proton");
  }
  else if(ParticleFlag==6){
    particle = particleTable->FindParticle("kaon+");
  }
  // ----------------------------------------------------
  else{
    particle = particleTable->FindParticle("e-"); // comment out 17Oct2013
    //particle = particleTable->FindParticle("pi+");  // for test
  }
  
}

TARGPrimaryGeneratorAction::~TARGPrimaryGeneratorAction()
{
  delete particleGun;
}

void TARGPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // ~~~~~~~~ Original ~~~~~~~~~~~
  //GenParPoint(anEvent);
  
  // ~~~~~~~~ Uniform generation ~~~~~
  //GenParUni(anEvent);
  if( SeedFlag==0 || SeedFlag==1 ){
    GenParUni(anEvent);
  }
  //else{
  //  GenParSeed();
  //}
}


//############################################ 
//   GenParPoint                        ######
//############################################ 
void TARGPrimaryGeneratorAction::GenParPoint(G4Event* anEvent){
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  particleGun -> SetParticleEnergy(0.0*CLHEP::GeV);// reset
  paramMan->SetBeamPos( G4ThreeVector(0.0*CLHEP::cm,0.0*CLHEP::cm,0.0*CLHEP::cm) );// reset
  //particleGun -> SetParticleDefinition( G4Geantino::GeantinoDefinition() );
  //BeamMom = BeamMom + ( 50.0e-3 * ( 2.0*G4UniformRand()-1.0 ) ) ;
  particleGun -> SetParticleMomentum( BeamMom * CLHEP::GeV );// Momentum in GeV
  particleGun -> SetParticleDefinition( particle );   // Particle type
  particleGun -> SetParticleMomentumDirection( G4ThreeVector(0.0,0.0,1.0) );
  G4double boffx = BeamOffset.x();
  G4double boffy = BeamOffset.y();
  G4double boffz = BeamOffset.z();
  if(GenFlag==2){ // Rastering
    //Beamx = rasterX*(G4UniformRand()-0.5)*2./2.0*cm;//toshi
    //Beamy = rasterY*(G4UniformRand()-0.5)*2./2.0*cm;//toshi
    //Beamz = rasterZ*(G4UniformRand()-0.5)*2./2.0*cm;//toshi
    Beamx = ( rasterX*(G4UniformRand()-0.5)*2./2.0 + boffx )*CLHEP::cm;//toshi
    Beamy = ( rasterY*(G4UniformRand()-0.5)*2./2.0 + boffy )*CLHEP::cm;//toshi
    Beamz = ( rasterZ*(G4UniformRand()-0.5)*2./2.0 + boffz )*CLHEP::cm;//toshi
    particleGun -> SetParticlePosition( G4ThreeVector( Beamx , Beamy , Beamz ) );
    //G4cout << G4UniformRand() << G4endl;
    //G4cout << Beamx/cm << " , " << Beamy/cm << " , " << Beamz/cm << G4endl;
  }
  else{ // No Rastering
    Beamx  =  (0.0 + boffx )*CLHEP::cm;
    Beamy  =  (0.0 + boffy )*CLHEP::cm;
    Beamz  =  (-2.0 + boffz)*CLHEP::cm;
    particleGun -> SetParticlePosition( G4ThreeVector( Beamx , Beamy , Beamz ) );
  }
  paramMan->SetBeamPos( G4ThreeVector(Beamx/CLHEP::cm,Beamy/CLHEP::cm,Beamz/CLHEP::cm) );
  //G4int i = anEvent -> GetEventID() % 3;
  //switch(i)
  // {
  //  case 0:
  //    particleGun -> SetParticleMomentumDirection(G4ThreeVector(1.0,0.0,0.0));
  //  case 1:
  //    particleGun -> SetParticleMomentumDirection(G4ThreeVector(0.0,1.0,0.0));
  //  case 2:
  //    particleGun -> SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,1.0));
  //  }
  particleGun -> GeneratePrimaryVertex(anEvent);
}

//############################################ 
//   GenParUni                          ######
//############################################ 
void TARGPrimaryGeneratorAction::GenParUni( G4Event* anEvent )
{
  //HKSParamMan *paramMan = HKSParamMan::GetParamMan();

  G4double px,py,pz; 
  //G4double CentMom = 0.1;//GeV/c
  G4double CentMom = bpc;
  G4double AcptMom = bpw;
  G4double CentTheta = bthc;
  G4double AcptTheta = bthw;
  G4double CentPhi = bphc;
  G4double AcptPhi = bphw;
  //G4double ThetaLimit = paramMan->GetThetaLimit();
  G4double xR = rasterX;
  G4double yR = rasterY;
  G4double zR = rasterZ;
  G4double Mom = (CentMom + (G4UniformRand()-0.5)*2.0*AcptMom);
  G4double xp = 0.;
  G4double yp = 0.;
  //G4int RanSign=rand()%2;
  G4double phi=0;
  G4double min = cos(CentTheta+AcptTheta);
  G4double max = cos(CentTheta-AcptTheta);
  G4double da = (max-min)/2. ;
  G4double a = min + da ;

  G4double theta=acos( a + ( 2.*da*(G4UniformRand()-0.5) ) );
  //if(RanSign==1)phi=-1.*(CentPhi+( 2.0*AcptPhi*( G4UniformRand()-0.5 ) ));
  //else if(RanSign==0)phi=1.*(CentPhi+( 2.0*AcptPhi*( G4UniformRand()-0.5 ) ));
  phi=1.*(CentPhi+( 2.0*AcptPhi*( G4UniformRand()-0.5 ) ));
  xp=sin(theta)*cos(phi)/cos(theta);
  yp=sin(theta)*sin(phi)/cos(theta);
  pz = Mom/sqrt(1+xp*xp+yp*yp);
  px = pz*xp;
  py = pz*yp;

  G4double boffx = BeamOffset.x();
  G4double boffy = BeamOffset.y();
  G4double boffz = BeamOffset.z();
  const int aa=1000;
  char str[aa];
  G4double temp;
  G4int evID;
  if(SeedFlag==0){
    // ~~~~~~~~~~ Rastering ~~~~~~~~~~~~~~~
    if(GenFlag==2){ 
      Beamx = ( xR*(G4UniformRand()-0.5)*2./2.0 + boffx )*CLHEP::cm;
      Beamy = ( yR*(G4UniformRand()-0.5)*2./2.0 + boffy )*CLHEP::cm;
      Beamz = ( zR*(G4UniformRand()-0.5)*2./2.0 + boffz )*CLHEP::cm;
    }
    // ~~~~~~~~~~ No rastering ~~~~~~~~~~~~~~~
    else{ 
      Beamx  =  (0.0 + boffx )*CLHEP::cm;
      Beamy  =  (0.0 + boffy )*CLHEP::cm;
      ///Beamz  =  (-2.0 + boffz)*CLHEP::cm;
      Beamz  =  (0.0 + boffz)*CLHEP::cm;
    }
    
  }
  else{ // ~~~~~~~~~~~~~ From Seed file  ~~~~~~~~~~~~~~~~
    fgets( str, aa, sfile );
    sscanf( str, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	    &evID,  &Beamx, &temp, &Beamy, &temp, &temp, &Beamz, &temp, &temp, &temp);
    paramMan->SetEVID(evID);
    //Double_t tang = 17.0;   // Target angle (HES-HKS)
    Double_t tang = 0.0;   // Target angle
    tang = tang * 3.14159 / 180.0; // degree --> rad
    Beamz = -1.0 * (Beamz + 2.0 * Beamx*tan(tang) ); // 11Dec2013 , Toshi
    Beamx = Beamx * CLHEP::cm;
    Beamy = Beamy * CLHEP::cm;
    Beamz = Beamz * CLHEP::cm;
  }
  //paramMan->SetBeamPos( G4ThreeVector(Beamx/CLHEP::cm,Beamy/CLHEP::cm,Beamz/CLHEP::cm) );
  paramMan->SetBeamPos( G4ThreeVector(Beamx,Beamy,Beamz) );
  paramMan->SetGenP(Mom); // [GeV]
  paramMan->SetGenPvec(px,py,pz);// [GeV]
  paramMan->SetGenTheta(theta);
  paramMan->SetGenPhi(phi);
  
  //G4ThreeVector gPos(x, y, z);
  G4ThreeVector gPos(Beamx, Beamy, Beamz);
  //G4ThreeVector gMom(px*CLHEP::GeV, py*CLHEP::GeV, pz*CLHEP::GeV);
  G4ThreeVector gMom(px, py, pz);
    
  SetMom(gMom,gPos,anEvent,0);
}


////############################################ 
////   GenParSeed                         ######
////############################################
//void TARGPrimaryGeneratorAction::GenParSeed( G4Event *anEvent )
//{
//  //HKSParamMan *paramMan = HKSParamMan::GetParamMan();
//  const int aa=500;
//  char str[aa];
//  fgets( str, aa, SeedFile );
//  
//  G4double px, py, pz, Mom;
//  G4double x0, y0, z0;
//  G4double pe;
//  G4double xp, yp;
//  G4double xp0, yp0;
//  G4double xpk,ypk,pk;
//  G4double xR = paramMan->GetRasterX();
//  G4double yR = paramMan->GetRasterY();
//  G4double zR = paramMan->GetRasterZ();
//  G4int evID;
//  
//  sscanf( str, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
//          &evID,  &x0, &xp0, &y0, &yp0, &pe, &z0, &xpk, &ypk, &pk);
//  Mom = pk;
//  
//  
//  xp = xpk/1000.;
//  yp = ypk/1000.;
//  pz = (Mom/sqrt(1+xp*xp+yp*yp));
//  px = pz*xp;
//  py = pz*yp;
//  G4ThreeVector gMom( px*GeV, py*GeV,pz*GeV);  
//  G4ThreeVector gPos( xR*x0*cm, yR*y0*cm, zR*z0*cm);
//  //G4ThreeVector gPos( 0.*cm, 0.*cm, 0.*cm);
//  
//  SetMom(gMom,gPos,anEvent,evID);
//}

//############################################ 
//   SetMom                             ######
//############################################
void TARGPrimaryGeneratorAction::SetMom
(G4ThreeVector gMom, G4ThreeVector gPos, G4Event* anEvent, G4int evID)
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  
  particleGun -> SetParticleEnergy(0.0*CLHEP::GeV); //reset
  particleGun -> SetParticleDefinition(particle); 
  particleGun -> SetParticleMomentum(gMom);
  particleGun -> SetParticlePosition(gPos);
  particleGun -> GeneratePrimaryVertex(anEvent);
  //  G4cout  << gPos.x()/cm << " " << gPos.y()/cm << " "  << gPos.z()/cm << G4endl;

  //if( anaMan ){
  //  anaMan->PrimaryGeneration( gPos, gMom, evID );
  //}
  //delete particleGun;

}
