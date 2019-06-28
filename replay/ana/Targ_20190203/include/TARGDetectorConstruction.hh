/*
  TARGDetectorConstruction.hh
  1 Oct 2010
  T.Gogami
*/

#ifndef TARGDetectorConstruction_h

#define TARGDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "TARGParamManager.hh"
#include "G4LogicalVolume.hh"

class MaterialList;
//class TARGField;

class TARGDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  TARGDetectorConstruction();
  ~TARGDetectorConstruction();
  TARGDetectorConstruction(TARGParamManager*);

  G4VPhysicalVolume* Construct();

private:
  G4double dist_t2window = 169.0*CLHEP::cm;
  G4double lhrs_ang = 13.2*CLHEP::deg;
  G4double rhrs_ang = 13.2*CLHEP::deg;
  MaterialList *mList_;
  MaterialList *DefineMaterials( void );
  G4String TargetMaterial;
  G4double TThickness;    // Added by Toshi , 18Oct2013
  G4LogicalVolume *targetLV;
  G4LogicalVolume *target_cellLV;
  G4LogicalVolume *scat_chambLV;
  G4LogicalVolume *scat_chambwinLV;
  G4LogicalVolume *be_isolatorLV;
  G4LogicalVolume *hrs_window1LV;
  G4LogicalVolume *hrs_window2LV;  

  //TARGField *MakeDetectorField( void );
  //TARGField *EMField_;

  G4VPhysicalVolume* ConstructPayload();
  TARGParamManager* paramMan;
};

#endif

