/*
  Stepping Action
  2 Oct 2010
  T.Gogami
*/

#include "TARGSteppingAction.hh"
#include "G4SteppingManager.hh"

#include "TARGEventAction.hh"
#include "TARGRunAction.hh"


TARGSteppingAction::TARGSteppingAction()
{ }

void TARGSteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
  //G4StepPoint* endPoint = aStep->GetPostStepPoint();
  //G4String procName = endPoint->GetProcessDefinedStep()->GetProcessName();
  //G4cout << procName <<G4endl;

}
