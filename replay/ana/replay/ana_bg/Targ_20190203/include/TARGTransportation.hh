/*
  TARGTransportation.hh

  T. Gogami
*/

#ifndef TARGTransportation_h

#define TARGTransportation_h 1

#include "G4Transportation.hh"

class TARGTransportation : public G4Transportation
{
public:
  TARGTransportation()
    : G4Transportation() {}
  ~TARGTransportation() {}
  
  G4double AlongStepGetPhysicalInteractionLength( const G4Track & track, 
						  G4double previousStepSize,
						  G4double currentMinimumStep, 
						  G4double & currentSafety,
						  G4GPILSelection* selection );
  
};

#endif
