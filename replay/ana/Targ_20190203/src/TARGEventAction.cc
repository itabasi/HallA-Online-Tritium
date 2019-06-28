/*
  "TARGEventAction.cc"
  
  T. Gogami, November 23, 2018
*/

#include "TARGRunAction.hh"
#include "TARGEventAction.hh"
#include "TARGvdHit.hh"
#include "TARGvdSD.hh"
#include "TARGParamManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"

#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <math.h>
#include <TObjArray.h>
#include <time.h>

TARGEventAction::TARGEventAction()//Constructor 
  :primaryEventID(0)
  ,beamx(0.0),beamy(0.0),beamz(0.0),beamp(0.0)
  ,evID(0),trackID(0)
   //,x(0),y(0),z(0),px(0),py(0),pz(0),p(0)
   //,theta(0),phi(0),charge(0)
  ,pname("DRAGON")
  ,pioflag(0),kaoflag(0),proflag(0),eflag(0),phoflag(0),posflag(0)
  ,paramMan(0),PosiAnalFlag(0)
  ,ProcID(0),passth(0),eIoni(0),eBrem(0),compt(0),phoele(0),conv(0)
  ,minMom(0),maxMom(0),minTheta(0),maxTheta(0)
{}

TARGEventAction::~TARGEventAction()//Destructor
{}

TARGEventAction::TARGEventAction(TARGParamManager* PmanDragon)//Constructor 
  :primaryEventID(0)
  ,beamx(0.0),beamy(0.0),beamz(0.0),beamp(0.0)
  ,evID(0),evID2(0),trackID(0),
   ntrig2(0),ntrig3(0),ntrig4(0)
   //,x(0),y(0),z(0),px(0),py(0),pz(0),p(0)
   //,theta(0),phi(0),charge(0)
  , pname("DRAGON")
  ,pioflag(0),kaoflag(0),proflag(0),eflag(0),phoflag(0),posflag(0)
  ,muoflag(0),paramMan(PmanDragon),PosiAnalFlag(0)
  ,ProcID(0),passth(0),eIoni(0),eBrem(0),compt(0),phoele(0),conv(0)
  ,minMom(0),maxMom(0),minTheta(0),maxTheta(0)
  ,beampx(0.0),beampy(0.0),beampz(0.0)
  ,beam_theta(0.0),beam_phi(0.0)
  ,SeedFlag(0),btrig(0)
{
  PosiAnalFlag = paramMan->GetPosiAnalFlag();
  minMom   = paramMan->GetMinMom();
  maxMom   = paramMan->GetMaxMom();
  minTheta = paramMan->GetMinTheta();
  maxTheta = paramMan->GetMaxTheta();
  SeedFlag = paramMan->GetSeedFlag();
  start = time(NULL);
  time(&start);
  //for(int i =0 ; i<10000000 ; i++){
  //  G4cout << "PmanDragon ===" << PmanDragon << G4endl;
  //}
}

void TARGEventAction::BeginOfEventAction(const G4Event*evt)
{
  primaryEventID++;
  //if( (evt->GetEventID() ) == 1){
  if( (evt->GetEventID() ) == 0){ // fixed (Toshi , 13Dec2013)
    definetree(evt); //Create Branches for ROOT file
  }
  //G4cout << "Event ID = " << evt->GetEventID() << G4endl;
}

void TARGEventAction::EndOfEventAction(const G4Event* evt)
{
  G4int event_id = evt->GetEventID();
  ///  Int_t event_id = evt->GetEventID();
  char anatime[100];
  
  //Periodic printing
  if((int)event_id % 10000 == 0 ){
    end = time(NULL);
    time(&end);
    sprintf( anatime,"%.0f Sec",difftime(end,start) );
    G4cout << "EvID: " << evt->GetEventID()
	   << ", L:" << ntrig2 << ", R:" << ntrig3 << ", B:" << ntrig4 
	   << " ( " << anatime << " ) " << G4endl;
   }
  //G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  //G4int n_trajectories = 0;
  //if(trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  //G4cout << "trajectory ontainer = " << n_trajectories << G4endl;
  //Analysis and Filling data to ROOT file
  filldata(evt);
}

//=============================================TOP//
//Define Branches for ROOT file
void TARGEventAction::definetree(const G4Event*){
  //TH1F* h;
  tree = dynamic_cast<TTree*>(gFile->Get("tree"));
  //tree->Branch("evID" ,  &evID  , "evID/I"   );
  tree->Branch("eventid",&evID  , "eventid/I"   );
  //tree->Branch("evID2",  &evID2 , "evID2/I"  );
  //tree->Branch("btrig",  &btrig , "btrig/I"  );
  //tree->Branch("trackID",&trackID,"trackID/I");
  tree->Branch("xBeam",  &beamx , "beamx/D");
  tree->Branch("yBeam",  &beamy , "beamy/D");
  tree->Branch("zBeam",  &beamz , "beamz/D");
  tree->Branch("pxBeam", &beampx ,"beampx/D");
  tree->Branch("pyBeam", &beampy ,"beampy/D");
  tree->Branch("pzBeam", &beampz ,"beampz/D");
  tree->Branch("pBeam",  &beamp , "beamp/D");
  tree->Branch("thetaBeam", &beam_theta, "thetaBeam/D");
  tree->Branch("phiBeam", &beam_phi  , "phiBeam/D");
  //tree->Branch("xVD"  ,  &x     , "x[20]/D"    );
  //tree->Branch("yVD"  ,  &y     , "y[20]/D"    );
  //tree->Branch("zVD"  ,  &z     , "z[20]/D"    );
  //tree->Branch("pxVD" ,  &px    , "px[20]/D"   );
  //tree->Branch("pyVD" ,  &py    , "py[20]/D"   );
  //tree->Branch("pzVD" ,  &pz    , "pz[20]/D"   );
  //tree->Branch("pVD"  ,  &p     , "p[20]/D"    );
  //tree->Branch("thetaVD",  &theta , "thetaVD[20]/D");
  //tree->Branch("phiVD"  ,  &phi   , "phiVD[20]/D"  );
  ////tree->Branch("pname",  &pname , "pname/C"  );
  //tree->Branch("charge" , &charge , "charge[20]/D"  );
  ////tree->Branch("posflag", &posflag, "posflag/I" );
  //tree->Branch("eleflag", &eleflag, "eleflag[20]/I" );
  //
  //tree->Branch("eleflag", &eleflag, "eleflag[20]/I" );
  
  //tree->Branch("phoflag", &phoflag, "phoflag/I" );
  //tree->Branch("proflag", &proflag, "proflag/I" );
  //tree->Branch("kaoflag", &kaoflag, "kaoflag/I" );
  //tree->Branch("pioflag", &pioflag, "pioflag/I" );
  //tree->Branch("muoflag", &muoflag, "muoflag/I" );
  //tree->Branch("passth" , &passth , "passth/I"  );
  //tree->Branch("eIoni"  , &eIoni  , "eIoni/I"   );
  //tree->Branch("eBrem"  , &eBrem  , "eBrem/I"   );
  //tree->Branch("compt"  , &compt  , "compt/I"   );
  //tree->Branch("phoele" , &phoele , "phoele/I"  );
  //tree->Branch("conv"   , &conv   , "conv/I"    );
  tree->Branch("nhit", &nhit,  "nhit/I" );
  tree->Branch("lnhit", &nhit2, "lnhit/I" );
  tree->Branch("lpx", &px2, "lpx[20]/D" );
  tree->Branch("lpy", &py2, "lpy[20]/D" );
  tree->Branch("lpz", &pz2, "lpz[20]/D" );
  tree->Branch("lp", &p2, "lp[20]/D" );
  tree->Branch("lx", &x2, "lx[20]/D" );
  tree->Branch("ly", &y2, "ly[20]/D" );
  //tree->Branch("lz", &z2, "lz[20]/D" );
  tree->Branch("ltheta", &theta2, "ltheta[20]/D" );
  tree->Branch("lphi"  , &phi2,   "lphi[20]/D" );
  tree->Branch("lpid"  , &pid2,   "lpid[20]/D" );
  tree->Branch("lcharge"  , &charge2,   "lcharge[20]/D" );

  tree->Branch("rnhit", &nhit3, "rnhit/I" );
  tree->Branch("rpx", &px3, "rpx[20]/D" );
  tree->Branch("rpy", &py3, "rpy[20]/D" );
  tree->Branch("rpz", &pz3, "rpz[20]/D" );
  tree->Branch("rp", &p3, "rp[20]/D" );
  tree->Branch("rx", &x3, "rx[20]/D" );
  tree->Branch("ry", &y3, "ry[20]/D" );
  //tree->Branch("rz", &z3, "rz[20]/D" );
  tree->Branch("rtheta", &theta3, "rtheta[20]/D" );
  tree->Branch("rphi"  , &phi3,   "rphi[20]/D" );
  tree->Branch("rpid"  , &pid3,   "rpid[20]/D" );
  tree->Branch("rcharge"  , &charge3,   "rcharge[20]/D" );

  tree->Branch("bnhit", &nhit4, "bnhit/I" );
  tree->Branch("bpx", &px4, "bpx[20]/D" );
  tree->Branch("bpy", &py4, "bpy[20]/D" );
  tree->Branch("bpz", &pz4, "bpz[20]/D" );
  tree->Branch("bp", &p4, "bp[20]/D" );
  tree->Branch("bx", &x4, "bx[20]/D" );
  tree->Branch("by", &y4, "by[20]/D" );
  //tree->Branch("bz", &z4, "bz[20]/D" );
  tree->Branch("btheta", &theta4, "btheta[20]/D" );
  tree->Branch("bphi"  , &phi4,   "bphi[20]/D" );
  tree->Branch("bpid"  , &pid4,   "bpid[20]/D" );
  tree->Branch("bcharge"  , &charge4,   "bcharge[20]/D" );
}
//==========================================BOTTOM//


//=============================================TOP//
//Analysis of data and Filling data to ROOT file
void TARGEventAction::filldata(const G4Event*evt){
  //const int nn;
  
  G4int event_id = evt->GetEventID();
  tree = dynamic_cast<TTree*>(gFile->Get("tree"));
  
  G4int colidVD;
  G4HCofThisEvent *HCE = evt -> GetHCofThisEvent();//Hits collection of this event
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();//SD manager
  //
  G4double xBeam = -2222.0 , yBeam= -2222.0 , zBeam= -2222.0;// Beam Position
  G4double pBeam = -2222.0 ; // Beam Momentum
  G4double pxBeam= -2222.0 ; // Beam Momentum
  G4double pyBeam= -2222.0 ; // Beam Momentum
  G4double pzBeam= -2222.0 ; // Beam Momentum
  G4double thBeam= -2222.0 ; // Beam theta
  G4double phBeam= -2222.0 ; // Beam phi
  //G4double xBeam[nn] , yBeam[nn], zBeam[nn];// Beam Position
  //G4double pBeam[nn]; // Beam Momentum
  //G4double pxBeam[nn]; // Beam Momentum
  //G4double pyBeam[nn]; // Beam Momentum
  //G4double pzBeam[nn]; // Beam Momentum
  //G4double thBeam[nn]; // Beam theta
  //G4double phBeam[nn]; // Beam phi
  G4ThreeVector pvecBeam(0.0,0.0,0.0);
  G4double xVD  = -2222.0 , yVD  = -2222.0 , zVD  = -2222.0;//position
  G4double pxVD = -2222.0 , pyVD = -2222.0 , pzVD = -2222.0;//momentum(x,y,z)
  G4double pVD  = -2222.0 ;   // momentum
  G4double thetaVD= -2222.0 ; // theta
  G4double phiVD  = -2222.0 ; // phi
  G4double chargeVD=-2222.0 ; //Charge
  G4int piflagVD  = 0; //pion flag
  G4int posflagVD = 0; //positron flag
  G4int eleflagVD = 0; //electron flag
  G4int phoflagVD = 0; //photon flag
  G4int proflagVD = 0; //proton flag
  G4int kaoflagVD = 0; //proton flag
  G4int pioflagVD = 0; //pion flag
  G4int muoflagVD = 0; //muon flag
  G4int passthVD  = 0; //beam passing through the target
  G4int eIoniVD   = 0; //electron ionization
  G4int eBremVD   = 0; //electron brems
  G4int comptVD   = 0; //compton scattering
  G4int phoeleVD  = 0; //photoelectric effect
  G4int convVD    = 0; //photon conversion
  
  G4ThreeVector bpos(0.0,0.0,0.0); //Beam position
  G4ThreeVector pos;//Hit position on VD
  G4ThreeVector mom;//Momentum vector of the particle
  bpos   = paramMan -> GetBeamPos();
  xBeam  = bpos.x();
  yBeam  = bpos.y();
  zBeam  = bpos.z();
  beamx = xBeam/CLHEP::cm; 
  beamy = yBeam/CLHEP::cm; 
  beamz = zBeam/CLHEP::cm; 
  pBeam  = paramMan -> GetGenP();
  //beamp= pBeam*1000.0;//GeV--> MeV
  beamp = pBeam;
  beamp = beamp / CLHEP::MeV;
  thBeam = paramMan -> GetGenTheta();
  beam_theta = thBeam;
  phBeam = paramMan -> GetGenPhi();
  beam_phi   = phBeam;
  pvecBeam = paramMan ->GetGenPvec();
  beampx = pvecBeam.x()/CLHEP::MeV;
  beampy = pvecBeam.y()/CLHEP::MeV;
  beampz = pvecBeam.z()/CLHEP::MeV;
  evID2 = paramMan->GetEVID();  //Event ID
  
  //mom = 0.0;
  pxVD = -2222.0;     //Momentum x [MeV]
  pyVD = -2222.0;     //Momentum y [MeV]
  pzVD = -2222.0;     //Momentum z [MeV]
  //G4cout << evID2 << " " << xBeam << " " << yBeam << " " << zBeam << G4endl;
  btrig = 1;
  
  colidVD = SDMan -> GetCollectionID("VDCollection1");

  TARGvdHitsCollection *vdHC;
  vdHC  = HCE ? (TARGvdHitsCollection*)( HCE->GetHC(colidVD) ) : 0 ;
 // 
 // for(int i=0 ; i<nn ; i++){
 //   xBeam[i] = -2222.0;
 //   yBeam[i] = -2222.0;
 //   zBeam[i] = -2222.0; // Beam Position
 //   pBeam[i] = -2222.0; // Beam Momentum
 //   pxBeam[i]= -2222.0; // Beam Momentum
 //   pyBeam[i]= -2222.0; // Beam Momentum
 //   pzBeam[i]= -2222.0; // Beam Momentum
 //   thBeam[i]= -2222.0; // Beam theta
 //   phBeam[i]= -2222.0; // Beam phi
 // }
 //
  for(int i=0 ; i<maxhit ; i++){
    x[i]      = -2222.0;
    y[i]      = -2222.0;
    z[i]      = -2222.0;
    px[i]     = -2222.0;
    py[i]     = -2222.0;
    pz[i]     = -2222.0;
    p[i]      = -2222.0;
    theta[i]  = -2222.0;
    phi[i]    = -2222.0;
    charge[i] = -10;

    x2[i]    = -2222.0; // LHRS
    y2[i]    = -2222.0;
    z2[i]    = -2222.0;
    px2[i]   = -2222.0;
    py2[i]   = -2222.0;
    pz2[i]   = -2222.0;
    p2[i]    = -2222.0;
    theta2[i]= -2222.0;
    charge2[i]= -2222.0;
    phi2[i]  = -2222.0;

    x3[i]    = -2222.0; // RHRS
    y3[i]    = -2222.0;
    z3[i]    = -2222.0;
    px3[i]   = -2222.0;
    py3[i]   = -2222.0;
    pz3[i]   = -2222.0;
    p3[i]    = -2222.0;
    theta3[i]= -2222.0;
    charge3[i]=-2222.0;
    phi3[i]  = -2222.0;

    x4[i]    = -2222.0; // Beam
    y4[i]    = -2222.0;
    z4[i]    = -2222.0;
    px4[i]   = -2222.0;
    py4[i]   = -2222.0;
    pz4[i]   = -2222.0;
    p4[i]    = -2222.0;
    theta4[i]= -2222.0;
    charge4[i]=-2222.0;
    phi4[i]  = -2222.0;
  }
  
  G4int entVD = vdHC->entries();
  nhit = 0;
  if(vdHC && entVD<maxhit){
    //G4ThreeVector bpos(0.0,0.0,0.0); //Beam position
    //G4ThreeVector pos;//Hit position on VD
    //G4ThreeVector mom;//Momentum vector of the particle
    //G4cout << entVD << G4endl;
    //G4cout << "entries() " <<entVD <<G4endl;
    //xBeam  =  -2222.0 ; yBeam= -2222.0 ; zBeam= -2222.0 ;
    //pBeam  =  -2222.0 ;
    bpos   = paramMan -> GetBeamPos();
    xBeam  = bpos.x();
    yBeam  = bpos.y();
    zBeam  = bpos.z();
    //pBeam  = paramMan -> GetBeamEnergy();
    pBeam  = paramMan -> GetGenP();
    //pBeam  = pBeam*1000.0; // GeV --> MeV
    //cout << entVD << endl;
    
    // oOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoO
    // Only one evnent will be stored when SeedFlag==1
    if(SeedFlag==1) entVD=1;
    // oOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoO
    
    for(G4int j=0 ; j<entVD ; j++){
      //Initializing
      beamx = -2222.0;
      beamy = -2222.0 ;
      beamz = -2222.0 ;
      beamp =  -2222.0 ;
      beampx=  -2222.0 ;
      beampy=  -2222.0 ;
      beampz=  -2222.0 ;
      xVD   =  -2222.0 ; yVD  = -2222.0 ; zVD  = -2222.0 ;
      pxVD  =  -2222.0 ; pyVD = -2222.0 ; pzVD = -2222.0 ;
      pVD   =  -2222.0 ;
      thetaVD=  -2222.0 ;
      phiVD  =  -2222.0 ;
      chargeVD= -2222.0 ;
      piflagVD  = 0; 
      posflagVD = 0;
      eleflagVD = 0;
      phoflagVD = 0;
      proflagVD = 0;
      kaoflagVD = 0;
      pioflagVD = 0;
      muoflagVD = 0;
      passthVD  = 0;
      eIoniVD   = 0;
      eBremVD   = 0;
      comptVD   = 0;
      phoeleVD  = 0;
      convVD    = 0;
      
      TARGvdHit *vdHit = (*vdHC)[j];
      //beamx= xBeam;
      //beamy= yBeam;
      //beamz= zBeam;
      //beamp= pBeam*1000.0;//GeV--> MeV
      //beamp = pBeam/CLHEP::MeV;
      //beam_theta = thBeam;
      //beam_phi   = phBeam;
      pos  = vdHit->GetPos(); //Get hit position on Virtual detector
      xVD  = pos.x()/CLHEP::cm;      //x on VD [cm]
      yVD  = pos.y()/CLHEP::cm;      //y on VD [cm]
      zVD  = pos.z()/CLHEP::cm;      //z on VD [cm]
      mom  = vdHit->GetMom(); //Get momentum on Virtual detector
      pxVD = mom.x()/CLHEP::MeV;     //Momentum x [MeV]
      pyVD = mom.y()/CLHEP::MeV;     //Momentum y [MeV]
      pzVD = mom.z()/CLHEP::MeV;     //Momentum z [MeV]
      pVD  = sqrt( pxVD*pxVD + pyVD*pyVD + pzVD*pzVD );//Momentum [MeV]
      thetaVD = acos( pzVD/pVD  ); // Theta
      //phiVD   = atan( pyVD/pxVD ); // Phi
      phiVD   = atan2( pyVD,pxVD ); // Phi
      chargeVD= vdHit->GetCharge();// Charge
      pname = vdHit->GetPname();       //Get particle name
      if(pname=="e+")    posflagVD = 1;//positron  flag
      else posflagVD = 0;
      if(pname=="e-")    eleflagVD = 1;//electron flag
      else eleflagVD = 0;
      if(pname=="gamma") phoflagVD = 1;//photon flag
      else phoflagVD = 0;
      if(pname=="proton")proflagVD = 1;//proton flag
      else proflagVD = 0;
      if(pname=="kaon+" || pname=="kaon-" || pname =="kaon0") kaoflagVD = 1;//kaon flag
      else kaoflagVD = 0;
      if(pname=="pi+" || pname=="pi-" || pname =="pi0") pioflagVD = 1;//pion flag
      else pioflagVD = 0;
      if(pname=="mu-" || pname=="mu+") muoflagVD = 1;//muon flag
      else muoflagVD = 0;
      ProcID = vdHit -> GetProcID();   //Get Process ID
      if(ProcID == 0)passthVD= 1;      //Beam passing through the target
      else passthVD= 0;
      if(ProcID == 1)eIoniVD = 1;      //electron ionization
      else eIoniVD = 0;
      if(ProcID == 2)eBremVD = 1;      //electron brems
      else eBremVD = 0;
      if(ProcID == 3)comptVD = 1;      //compton scattering
      else comptVD = 0;
      if(ProcID == 4)phoeleVD= 1;      //photoelectric effect
      else phoeleVD= 0;
      if(ProcID == 5)convVD  = 1;      //photon conversion
      else convVD  = 0;
      
      // Fill data to tree
      //evID = primaryEventID;  //Event ID
      //trackID = event_id;     //Track ID
      evID  = event_id;  //Event ID
      //evID2 = paramMan->GetEVID();  //Event ID      
      
      x[j]   =   xVD ;     //x on VD [cm]
      y[j]   =   yVD ;     //y on VD [cm]
      z[j]   =   zVD ;     //z on VD [cm]
      px[j]  =   pxVD;     //Momentum x [GeV/c]
      py[j]  =   pyVD;     //Momentum y [GeV/c]
      pz[j]  =   pzVD;     //Momentum z [GeV/c]
      p[j]   =   pVD ;     //Momentum   [GeV/c]
      theta[j] = thetaVD;  //Theta
      if(phiVD<3.14159/2.){
	phiVD = phiVD;
      }
      else if (phiVD<3.14159*3.0/2.0){
	phiVD = phiVD + 3.14159;
      }
      else{
	phiVD = phiVD + 2.0*3.14159;
      }
      phi[j]   = phiVD;    //Phi
      charge[j] = chargeVD;//Charge
      //posflag = posflagVD;//positron flag
      //eleflag = eleflagVD;//electron flag
      //phoflag = phoflagVD;//photon flag
      //proflag = proflagVD;//proton flag
      //kaoflag = kaoflagVD;//proton flag
      //pioflag = pioflagVD;//pion flag
      //muoflag = muoflagVD;//muon flag
      //passth  = passthVD; //Beam passing through the target
      //eIoni   = eIoniVD;  //electron ionization
      //eBrem   = eBremVD;  //electron brems
      //compt   = comptVD;  //compton scattering
      //phoele  = phoeleVD; //photoelectric effect
      //conv    = convVD;   //photon conversion
      nhit++;
      
      //if( minTheta/1000.0 < theta && theta < maxTheta/1000.0 ){ //Theta cut
      //	if( minMom < p && p < maxMom ){                         //Momentum cut
      //	  if( PosiAnalFlag==0 ){
      //	    tree->Fill();  //Filling data for all particles
      //	  }// if 
      //	  else{
      //	    if(posflag == 1)tree->Fill(); //Filling data for only positrons
      //	  }
      //	}//if
      //}//if
    }//for
  }//if

  // ----------------------------------- //
  // ---- Virtual detector for LHRS ---- //
  // ----------------------------------- //
  colidVD = SDMan -> GetCollectionID("VDCollection2");
  TARGvdHitsCollection *vdHC2;
  vdHC2  = HCE ? (TARGvdHitsCollection*)( HCE->GetHC(colidVD) ) : 0 ;
  G4int entVD2 = vdHC2->entries();
  nhit2 = entVD2;
  
  if(vdHC2 && entVD2<maxhit){
    
    // oOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoO
    // Only one evnent will be stored when SeedFlag==1
    if(SeedFlag==1) entVD2=1;
    // oOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoO
    
    for(G4int j=0 ; j<entVD2 ; j++){
      TARGvdHit *vdHit = (*vdHC2)[j];
      pos  = vdHit->GetPos(); //Get hit position on Virtual detector
      x2[j]  = pos.x()/CLHEP::cm;      //x on VD [cm]
      y2[j]  = pos.y()/CLHEP::cm;      //y on VD [cm]
      z2[j]  = pos.z()/CLHEP::cm;      //z on VD [cm]
      mom  = vdHit->GetMom(); //Get momentum on Virtual detector
      px2[j] = mom.x()/CLHEP::MeV;     //Momentum x [MeV]
      py2[j] = mom.y()/CLHEP::MeV;     //Momentum y [MeV]
      pz2[j] = mom.z()/CLHEP::MeV;     //Momentum z [MeV]
      p2[j]  = sqrt( px2[j]*px2[j] + py2[j]*py2[j] + pz2[j]*pz2[j] );
      theta2[j] = acos( pz2[j]/p2[j]  ); // Theta
      phi2[j]   = atan2( py2[j],px2[j] ); // Phi
      charge2[j]= vdHit->GetCharge();  // Charge
      pname = vdHit->GetPname();       //Get particle name
      if(pname=="e-" || pname=="e+" || pname=="gamma") {
	pid2[j] = 0;
	if(charge2[j]<0) ntrig2++;
      }
      else if(pname=="proton"){
	pid2[j] = 1;
      }
      else if(pname=="kaon+" || pname=="kaon-" || pname =="kaon0"){
	pid2[j] = 2;
      }
      else if(pname=="pi+" || pname=="pi-" || pname =="pi0") {
	pid2[j] = 3;//pion flag
      }
      else pid2[j] = -1;
    }
  }

  // ----------------------------------- //
  // ---- Virtual detector for RHRS ---- //
  // ----------------------------------- //
  colidVD = SDMan -> GetCollectionID("VDCollection3");
  TARGvdHitsCollection *vdHC3;
  vdHC3  = HCE ? (TARGvdHitsCollection*)( HCE->GetHC(colidVD) ) : 0 ;
  G4int entVD3 = vdHC3->entries();
  nhit3 = entVD3;

  if(vdHC3 && entVD3<maxhit){
    
    // oOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoO
    // Only one evnent will be stored when SeedFlag==1
    if(SeedFlag==1) entVD3=1;
    // oOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoO
    
    for(G4int j=0 ; j<entVD3 ; j++){
      TARGvdHit *vdHit = (*vdHC3)[j];
      pos  = vdHit->GetPos(); //Getting hit position at the virtual detector
      x3[j]  = pos.x()/CLHEP::cm;      //x at VD [cm]
      y3[j]  = pos.y()/CLHEP::cm;      //y at VD [cm]
      z3[j]  = pos.z()/CLHEP::cm;      //z at VD [cm]
      mom  = vdHit->GetMom(); //Getting momentum at the virtual detector
      px3[j] = mom.x()/CLHEP::MeV;     //Momentum x [MeV]
      py3[j] = mom.y()/CLHEP::MeV;     //Momentum y [MeV]
      pz3[j] = mom.z()/CLHEP::MeV;     //Momentum z [MeV]
      p3[j]  = sqrt( px3[j]*px3[j] + py3[j]*py3[j] + pz3[j]*pz3[j] );
      theta3[j] = acos( pz3[j]/p3[j]  ); // Theta
      phi3[j]   = atan2( py3[j],px3[j] ); // Phi
      charge3[j]= vdHit->GetCharge();  // Charge
      pname = vdHit->GetPname();       //Get particle name
      
      if(pname=="e-" || pname=="e+" || pname=="gamma") {
	pid3[j] = 0;
      }
      else if(pname=="proton"){
	pid3[j] = 1;
      }
      else if(pname=="kaon+" || pname=="kaon-" || pname =="kaon0"){
	pid3[j] = 2;
	if(charge3[j]>0) ntrig3++;
      }
      else if(pname=="pi+" || pname=="pi-" || pname =="pi0") {
	pid3[j] = 3;//pion flag
      }
      else pid3[j] = -1;

      //cout << pname << " " << pid3[j] << endl;
    }
  }

  // ----------------------------------- //
  // ---- Virtual detector for Beam ---- //
  // ----------------------------------- //
  colidVD = SDMan -> GetCollectionID("VDCollection4");
  TARGvdHitsCollection *vdHC4;
  vdHC4  = HCE ? (TARGvdHitsCollection*)( HCE->GetHC(colidVD) ) : 0 ;
  G4int entVD4 = vdHC4->entries();
  nhit4 = entVD4;

  if(vdHC4 && entVD4<maxhit){
    
    // oOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoO
    // Only one evnent will be stored when SeedFlag==1
    if(SeedFlag==1) entVD4=1;
    // oOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoO
    
    for(G4int j=0 ; j<entVD4 ; j++){
      TARGvdHit *vdHit = (*vdHC4)[j];
      pos  = vdHit->GetPos(); //Getting hit position at the virtual detector
      x4[j]  = pos.x()/CLHEP::cm;      //x at VD [cm]
      y4[j]  = pos.y()/CLHEP::cm;      //y at VD [cm]
      z4[j]  = pos.z()/CLHEP::cm;      //z at VD [cm]
      mom  = vdHit->GetMom(); //Getting momentum at the virtual detector
      px4[j] = mom.x()/CLHEP::MeV;     //Momentum x [MeV]
      py4[j] = mom.y()/CLHEP::MeV;     //Momentum y [MeV]
      pz4[j] = mom.z()/CLHEP::MeV;     //Momentum z [MeV]
      p4[j]  = sqrt( px4[j]*px4[j] + py4[j]*py4[j] + pz4[j]*pz4[j] );
      theta4[j] = acos( pz4[j]/p4[j]  ); // Theta
      phi4[j]   = atan2( py4[j],px4[j] ); // Phi
      charge4[j]= vdHit->GetCharge();  // Charge
      pname = vdHit->GetPname();       //Get particle name
      if(pname=="e-" || pname=="e+" || pname=="gamma") {
	pid4[j] = 0;
	if(charge4[j]<0) ntrig4++;
      }
      else if(pname=="proton"){
	pid4[j] = 1;
      }
      else if(pname=="kaon+" || pname=="kaon-" || pname =="kaon0"){
	pid4[j] = 2;
      }
      else if(pname=="pi+" || pname=="pi-" || pname =="pi0") {
	pid4[j] = 3;//pion flag
      }
      else pid4[j] = -1;
    }
  }
  
  

  //G4cout << evID2 << " " << xBeam << " " << yBeam << " " << zBeam << G4endl;
  if(nhit>0 || nhit2>0 || nhit3>0 || nhit4>0){
    tree->Fill();  //Filling data for all particles
  }
  btrig = 0;
  
}
//==========================================BOTTOM//

