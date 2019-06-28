/*
  TARGvdHit.cc
  5 Oct 2010
  T.Gogami
*/

#include "TARGvdHit.hh"

#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4ios.hh"

G4Allocator<TARGvdHit> TARGvdHitAllocator;

TARGvdHit::TARGvdHit()
{pLogV=NULL;}

TARGvdHit::TARGvdHit(G4LogicalVolume* logVol)
:pLogV(logVol)
{;}

TARGvdHit::~TARGvdHit()
{;}

TARGvdHit::TARGvdHit(const TARGvdHit &right)
  : G4VHit()
{
  edep = right.edep;
  pos = right.pos;
  rot = right.rot;
  mom = right.mom;
  pname = right.pname;
  procID= right.procID;
  pLogV = right.pLogV;
}

const TARGvdHit& TARGvdHit::operator=(const TARGvdHit &right)
{
  edep = right.edep;
  pos = right.pos;
  rot = right.rot;
  mom = right.mom;
  pname = right.pname;
  procID= right.procID;
  pLogV = right.pLogV;
  return *this;
}

G4int TARGvdHit::operator==(const TARGvdHit &right) const
{
  return (this==&right) ? 1 : 0;
}

void TARGvdHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Transform3D trans(rot,pos);
    G4VisAttributes attribs;
    const G4VisAttributes* pVA = pLogV->GetVisAttributes();
    if(pVA) attribs = *pVA;
    G4Colour colour(1.,0.,0.);
    attribs.SetColour(colour);
    attribs.SetForceWireframe(false);
    attribs.SetForceSolid(true);
    pVVisManager->Draw(*pLogV,attribs,trans);
  }
}

void TARGvdHit::Print()
{
}


