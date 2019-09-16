/*
  Stepping Action
  2 Oct 2010
  T.Gogami
*/

#ifndef TARGSteppingAction_h
#define TARGSteppingAction_h 1

#include "G4UserSteppingAction.hh"


class TARGSteppingAction : public G4UserSteppingAction
{
  public:
    TARGSteppingAction();
   ~TARGSteppingAction(){};

    void UserSteppingAction(const G4Step*);
};


#endif
