/*
  "TARGDetectorConstruction.cc"
  1 Oct 2010, T.Gogami
*/

#include "TARGDetectorConstruction.hh"
#include "TARGvdHit.hh"
#include "TARGvdSD.hh"
#include "TARGParamManager.hh"

//=== Material Define ===//
#include "MaterialList.hh"

//=== Geometry Define ===//
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4DisplacedSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4AffineTransform.hh"
#include "G4SDManager.hh"
#include "G4Colour.hh"

//=== ElectoMagnetic Field ===//
//#include "TARGField.hh"

#include <sstream>
#include <iomanip>
#include <string>

//#include "TARGParamMan.hh"
//#include "TARGPrimaryGeneratorAction.hh"

//=== Program Start ===//
///////////////////////////////////////////////////
TARGDetectorConstruction::TARGDetectorConstruction()
  : mList_(0)
///////////////////////////////////////////////////
{
}
TARGDetectorConstruction::TARGDetectorConstruction(TARGParamManager*PManDragon)
  : mList_(0),paramMan(PManDragon)
///////////////////////////////////////////////////
{
  TargetMaterial = paramMan->GetTarget();     // Target material
  TThickness     = paramMan->GetTThickness(); // Target thickness [cm]
  //G4cout << paramMan->GetBeamEnergy() << G4endl;
  
  // ---- Check the target information -------------------------------
  /* 
  for(int i=0 ; i<10000 ; i++){
    G4cout << TargetMaterial << " " << TThickness << " cm "  << G4endl;
  }
  */
  // ------------------------------------------------------------------

}

////////////////////////////////////////////////////
TARGDetectorConstruction::~TARGDetectorConstruction()
///////////////////////////////////////////////////
{
  delete mList_;
}

/////////////////////////////////////////////////////////
MaterialList *TARGDetectorConstruction::DefineMaterials()
/////////////////////////////////////////////////////////
{
  MaterialList *ml = new MaterialList;

  return ml;
}

////////////////////////////////////////////////////////
G4VPhysicalVolume *TARGDetectorConstruction::Construct()
///////////////////////////////////////////////////////
{
  mList_ = DefineMaterials();

  G4VPhysicalVolume *world = ConstructPayload();

  return world;
}

/////////////////////////////////////////////////////////
//TARGField *TARGDetectorConstruction::MakeDetectorField()
/////////////////////////////////////////////////////////
//{
//  return new TARGField();
//}

///////////////////////////////////////////////////////////////
G4VPhysicalVolume *TARGDetectorConstruction::ConstructPayload()
//////////////////////////////////////////////////////////////
{ 
  //=== World ===//
  G4double WorldSizeX = 300.0*CLHEP::cm;
  G4double WorldSizeY = 100.0*CLHEP::cm;
  G4double WorldSizeZ = 400.0*CLHEP::cm;
  G4Box *worldSolid = new G4Box( "World",WorldSizeX/2.0,
				 WorldSizeY/2.0,WorldSizeZ/2.0 );
  
  G4LogicalVolume *worldLV =
    new G4LogicalVolume( worldSolid, mList_->Vacuum, "World LV");
    //new G4LogicalVolume( worldSolid, mList_->Air, "World LV");
  
  G4VPhysicalVolume *world = 
    new G4PVPlacement( 0, 
		       G4ThreeVector( 0.*CLHEP::m, 0.*CLHEP::m, 0.*CLHEP::m ),
		       worldLV,
		       "World",
		       0,
		       false, 
		       0);
  
  G4double zpos_for_VD = 0.5*CLHEP::cm;
  //=== Target ===//
  //G4double TargetThickness = 100.0;
  //std::string TargetMaterial = "Carbon";
  G4Material *Target =mList_->Vacuum;
  G4double TargetSizeX = 1.27 * CLHEP::cm ;
  G4double TargetSizeY = 1.27 * CLHEP::cm ;
  G4double TargetSizeZ = 0.1 * 0.1 * CLHEP::cm ;
  TargetSizeZ = TThickness * CLHEP::cm ;
  
  G4RotationMatrix* rotTarget = new G4RotationMatrix();//Rotaion of target
  //rotTarget->rotateY(-17.0 * CLHEP::deg);    //Rotation of the target for HKS-HES
  //rotTarget->rotateY(-54.0 * deg);    //Rotation of the target for KaoS
  rotTarget->rotateY( 0.0 * CLHEP::deg);      // For energy loss study, 18Oct2013
  G4Tubs* tcell_tube; // target cell tube
  G4Sphere* tcell_end;
  G4Tubs* tcell_front; // target cell tube
  G4Sphere* tcell_front_sphere; //
  G4Tubs* gastar_main;
  G4Sphere* gastar_end;
  G4Sphere* gastar_front;
  G4RotationMatrix* noRot = new G4RotationMatrix();
  G4double scat_chamber_rad = 1037.0/2.0*CLHEP::mm;
  G4double scat_chamber_window_thick = 0.279*CLHEP::mm;

  G4Tubs* scat_chamber_tube
    = new G4Tubs("scat_chamber_tube",
		 0.0, scat_chamber_rad,
		 600.0/2.0 * CLHEP::mm,
		 0.0, 360.0*CLHEP::deg);
  G4Material * sc_material =mList_->Vacuum;
  scat_chambLV = new G4LogicalVolume(scat_chamber_tube,
				     sc_material, 
				     "Scattering Chamber LV");

  G4Tubs* scat_chamber_window_tube
    = new G4Tubs("scat_chamber_tube",
		 scat_chamber_rad + 0.01*CLHEP::mm,
		 scat_chamber_rad + scat_chamber_window_thick + 0.01*CLHEP::mm,
		 200.0/2.0 * CLHEP::mm,
		 0.0, 180.0*CLHEP::deg);
  G4Material * scwin_material =mList_->Al;
  scat_chambwinLV
    = new G4LogicalVolume(scat_chamber_window_tube,
			  scwin_material, 
			  "Scattering Chamber Window LV");

  G4Tubs* be_isolator_tube
    = new G4Tubs("be_isolator_tube",
		 0.0, 2.0*CLHEP::cm,
		 0.216/2.0 * CLHEP::mm, // ~40 mg/cm2
		 0.0, 360.0*CLHEP::deg);
  G4Material * beiso_material =mList_->Be9;
  be_isolatorLV = new G4LogicalVolume(be_isolator_tube,
				      beiso_material,
				      "Scattering Chamber LV");

  G4Box *hrs_window1 = new G4Box("hrs_window1",
				 70.0*CLHEP::cm/2.0,
				 70.0*CLHEP::cm/2.0,
				 0.127*CLHEP::mm/2.0);
  G4Material * hrs_window1_material =mList_->Ti48;
  //G4Material * hrs_window1_material =mList_->Vacuum; // for test
  hrs_window1LV = new G4LogicalVolume(hrs_window1,
				      hrs_window1_material,
				      "HRS window material 1 LV");
  
  G4Box *hrs_window2 = new G4Box("hrs_window2",
				 73.0*CLHEP::cm/2.0,
				 73.0*CLHEP::cm/2.0,
				 0.127*CLHEP::mm/2.0);
  G4Material * hrs_window2_material =mList_->Ti48;
  //G4Material * hrs_window2_material =mList_->Vacuum; // for test
  hrs_window2LV = new G4LogicalVolume(hrs_window2,
				      hrs_window2_material,
				      "HRS window material 2 LV");
  
  
  
  if( TargetMaterial == "Carbon"){
    Target = mList_->C12;
  }
  //else if( TargetMaterial == "Silicon"){
  //  Target = mList_->Si28;
  //}
  else if( TargetMaterial == "ch2"){
    Target = mList_->CH2;
  }
  else if( TargetMaterial == "Lithium"){
    Target = mList_->Li7;
  }
  else if( TargetMaterial == "Berylium"){
    Target = mList_->Be9;
  }
  else if( TargetMaterial == "Boron"){
    Target = mList_->B10;
  }
  else if( TargetMaterial == "Li3N"){
    Target = mList_->Li3N;
  }
  else if (TargetMaterial == "LHe4" ||
	   TargetMaterial == "LHe3" ){
    if(TargetMaterial == "LHe4") Target = mList_->LHe4;
    else Target = mList_->LHe3;
    G4double target_zlen = TargetSizeZ;
    G4double tcell_tube_zlen = TargetSizeZ + 0.2*CLHEP::mm;
    G4double cellcover_thick = 0.2*CLHEP::mm;
    zpos_for_VD = tcell_tube_zlen/2.0 + cellcover_thick + 0.5*CLHEP::cm; // virtual detector position
    G4Material* Target_cell = mList_->Al;
    G4double innrad = 10.0/2.0 * CLHEP::mm; // target inner radius
    G4double outrad = innrad + 0.5*CLHEP::mm;
    tcell_tube =
      new G4Tubs("tcell_tube",
		 innrad+0.2*CLHEP::mm, outrad+0.2*CLHEP::mm,
		 tcell_tube_zlen/2.0,
		 0.0, 360.0*CLHEP::deg);
    tcell_front =
      new G4Tubs("tcell_front",
		  0.0, outrad+0.2*CLHEP::mm,
		  cellcover_thick/2.0,
		  0.0, 360.0*CLHEP::deg);
    G4UnionSolid* tcell_front_tube =
      new G4UnionSolid("tcell_front_tube",
		       tcell_tube,
		       tcell_front,
		       noRot,
		       G4ThreeVector(0.0,
				     0.0,
				     -(cellcover_thick/2.0 + tcell_tube_zlen/2.0 + 0.002*CLHEP::mm))
		       );
    G4UnionSolid* tcell_front_tube_end =
      new G4UnionSolid("tcell_front_tube_end",
		       tcell_front_tube,
		       tcell_front,
		       noRot,
		       G4ThreeVector(0.0,
				     0.0,
				     +(cellcover_thick/2.0 + tcell_tube_zlen/2.0 + 0.002*CLHEP::mm))
		       );
    target_cellLV = new G4LogicalVolume(tcell_front_tube_end,Target_cell,"Target cell LV");
    G4Colour colourTCell(0.8, 0.5, 1.0); //
    G4VisAttributes *TCellVisAtt = new G4VisAttributes(true, colourTCell);
    target_cellLV->SetVisAttributes(TCellVisAtt);

    gastar_main = // (this is not a gas target, though)
      new G4Tubs("gastar_main",
		 0.0, innrad,
		 target_zlen/2.0,
		 0.0, 360.0*CLHEP::deg);
    targetLV = new G4LogicalVolume(gastar_main, Target, "Target LV");
    rotTarget->rotateY( 0.0 * CLHEP::deg );

//    new G4PVPlacement( noRot
//		       , G4ThreeVector(0*CLHEP::cm,
//				       0*CLHEP::cm,
//				       0*CLHEP::cm)
//		       , targetLV
//		       , "Target in target cell" 
//		       , target_cellLV
//		       , false
//		       , 0 );
//    new G4PVPlacement( rotTarget
//		       , G4ThreeVector(0*CLHEP::cm,
//				       0*CLHEP::cm,
//				       //-(cellcover_thick + tcell_tube_zlen/2.0 + 0.02*CLHEP::mm))
//				       0*CLHEP::cm)
//		       , target_cellLV
//		       , "Target in world " 
//		       , worldLV 
//		       , false
//		       , 0 );
//
    new G4PVPlacement( rotTarget
		       , G4ThreeVector(0*CLHEP::cm,
				       0*CLHEP::cm,
				       0*CLHEP::cm)
		       , targetLV
		       , "Target in world" 
		       , worldLV
		       , false
		       , 0 );
    new G4PVPlacement( rotTarget
		       , G4ThreeVector(0*CLHEP::cm,
				       0*CLHEP::cm,
				       //-(cellcover_thick + tcell_tube_zlen/2.0 + 0.02*CLHEP::mm))
				       0*CLHEP::cm)
		       , target_cellLV
		       , "Target cell in world " 
		       , worldLV 
		       , false , 0 );
    
  }
  else if( TargetMaterial == "Tritium"){
    // --------------------------------------------- //
    // ------ Cell volume is 33.38 cm^{3}     ------ //
    // ------ H2 = 77.00 +/- 0.01 mg/cm^{2}   ------ //
    // --------------------------------------------- //
    // The cell and target geometries are not incorporated yet //
    // T. Gogami, November 20, 2018                            //
    // 
    G4Material *Target_cell =mList_->Al7075;
    //G4Material *Target_cell =mList_->Al;
    Target = mList_->H3Gas;
    //Target = mList_->Vacuum;
    zpos_for_VD = 13.0*CLHEP::cm; // virtual detector position
    G4double innrad = 12.7/2.0 * CLHEP::mm;
    G4double outrad1 = innrad + 0.5472*CLHEP::mm;
    //G4double outrad2 = innrad + 0.254 *CLHEP::mm;
    G4double front_thickness =  0.254 *CLHEP::mm;
    G4double tcell_tube_zlen = 250.0*CLHEP::mm - innrad;
    tcell_tube =
      new G4Tubs("tcell_tube",
		 //innrad, outrad1,
		 innrad+0.2*CLHEP::mm, outrad1+0.2*CLHEP::mm,
		 tcell_tube_zlen/2.0,
		 0.0, 360.0*CLHEP::deg);
    tcell_end =
      new G4Sphere("tcell_end",
		   //innrad, outrad1,
		   innrad+0.2*CLHEP::mm, outrad1+0.2*CLHEP::mm,
		   0.0*CLHEP::deg, 180.0*CLHEP::deg,
		   0.0*CLHEP::deg, 180.0*CLHEP::deg);
    G4RotationMatrix* rottmp1 = new G4RotationMatrix();
    rottmp1->rotateX(270.0*CLHEP::deg);
    G4UnionSolid* tcell_tube_end
      = new G4UnionSolid("tcell_tube_end",
			 tcell_tube,
			 tcell_end,
			 rottmp1,
			 G4ThreeVector( 0.0 , 0.0 , tcell_tube_zlen/2.0 )
			 );
    tcell_front
      = new G4Tubs("tcell_front",
		   0.0, outrad1,
		   front_thickness/2.0,
		   0.0, 360.0*CLHEP::deg);
    G4UnionSolid* tcell_front_tube_end
      = new G4UnionSolid("tcell_front_tube_end",
			 tcell_tube_end,
			 tcell_front,
			 noRot,
			 G4ThreeVector( 0.0,
					0.0,
					-(tcell_tube_zlen+front_thickness)/2.0)
			 );
    G4UnionSolid* tcell_front_tube
      = new G4UnionSolid("tcell_front_tube",
			 tcell_tube,
			 tcell_front,
			 noRot,
			 G4ThreeVector( 0.0,
					0.0,
					-(tcell_tube_zlen+front_thickness)/2.0)
			 );

    target_cellLV =  new G4LogicalVolume( tcell_front_tube, Target_cell, "Target cell LV"); // removed the front cell wall
    //target_cellLV =  new G4LogicalVolume( tcell_tube_end, Target_cell, "Target cell LV"); // removed the front cell wall
    //target_cellLV =  new G4LogicalVolume( tcell_front_tube_end, Target_cell, "Target cell LV");
    G4Colour colourTCell(0.8, 0.5, 1.0); //
    G4VisAttributes *TCellVisAtt = new G4VisAttributes(true, colourTCell);
    target_cellLV->SetVisAttributes(TCellVisAtt);
    
    gastar_main =
      new G4Tubs("gastar_main",
		 0.0, innrad, 
		 tcell_tube_zlen/2.0,
		 0.0, 360.0*CLHEP::deg);
    gastar_end =
      new G4Sphere("gastar_end",
		   0.0, innrad,
		   0.0*CLHEP::deg, 180.0*CLHEP::deg,
		   0.0*CLHEP::deg, 180.0*CLHEP::deg);
    G4UnionSolid* gastar_main_end
      = new G4UnionSolid("gastar_main_end",
			 gastar_main,
			 gastar_end,
			 rottmp1,
			 G4ThreeVector( 0.0 , 0.0 , tcell_tube_zlen/2.0 )
			 );
    
    targetLV =  new G4LogicalVolume( gastar_main_end, Target, "Target LV");
//    new G4PVPlacement( noRot 
//		       , G4ThreeVector(0.0*CLHEP::cm, 0.0*CLHEP::cm, 0.0*CLHEP::cm)
//		       , target_cellLV
//		       , "GasTarget" 
//		       , targetLV
//		       , false
//		       , 0 );
    
    //rotTarget->rotateZ( 90.0 * CLHEP::deg );
    //rotTarget->rotateZ( 270.0 * CLHEP::deg );
    rotTarget->rotateY( 180.0 * CLHEP::deg ); // front side back, Toshi, 29 April 2018
  }
  else if( TargetMaterial == "H2"){
    // --------------------------------------------- //
    // ------ Cell volume is 33.38 cm^{3}     ------ //
    // ------ H2 = 70.8 +/- 0.4 mg/cm^{2}     ------ //
    // ------ T. Gogami, Nov 20, 2018         ------ //
    // --------------------------------------------- //
    G4Material *Target_cell =mList_->Al7075;
    //G4Material *Target_cell =mList_->Vacuum;
    //G4Material *Target_cell =mList_->Al;
    Target = mList_->HGas;
    //Target = mList_->Vacuum;
    zpos_for_VD = 13.0*CLHEP::cm; // virtual detector position
    G4double innrad = 12.7/2.0 * CLHEP::mm;
    //G4double outrad1 = innrad + 0.5472*CLHEP::mm;
    //G4double outrad1 = innrad + 0.4*CLHEP::mm;
    G4double outrad1 = innrad + 0.3*CLHEP::mm;
    //G4double outrad2 = innrad + 0.254 *CLHEP::mm;
    G4double front_thickness =  0.311 *CLHEP::mm; // measured
    G4double exit_thickness =   0.330 *CLHEP::mm; // measured
    G4double outrad_front = innrad + front_thickness;
    G4double outrad_exit  = innrad + exit_thickness;
    G4double tcell_tube_zlen = 250.0*CLHEP::mm - innrad;
    tcell_tube =
      new G4Tubs("tcell_tube",
		 //innrad, outrad1,
		 innrad+0.001*CLHEP::mm, outrad1+0.001*CLHEP::mm,
		 tcell_tube_zlen/2.0,
		 0.0, 360.0*CLHEP::deg);
    tcell_end =
      new G4Sphere("tcell_end",
		   //innrad, outrad1,
		   innrad+0.001*CLHEP::mm, outrad_exit+0.001*CLHEP::mm,
		   0.0*CLHEP::deg, 180.0*CLHEP::deg,
		   0.0*CLHEP::deg, 180.0*CLHEP::deg);
    G4RotationMatrix* rottmp1 = new G4RotationMatrix();
    rottmp1->rotateX(270.0*CLHEP::deg);
    G4UnionSolid* tcell_tube_end
      = new G4UnionSolid("tcell_tube_end",
			 tcell_tube,
			 tcell_end,
			 rottmp1,
			 G4ThreeVector( 0.0 , 0.0 , tcell_tube_zlen/2.0 )
			 );
    tcell_front_sphere =
      new G4Sphere("tcell_front_sphere",
		   //innrad+0.05*CLHEP::mm, outrad_front+0.05*CLHEP::mm + 10*CLHEP::cm, // for check
		   innrad+0.001*CLHEP::mm, outrad_front+0.001*CLHEP::mm,
		   0.0*CLHEP::deg, 180.0*CLHEP::deg,
		   0.0*CLHEP::deg, 180.0*CLHEP::deg);
		   //180.0*CLHEP::deg, 360.0*CLHEP::deg,
		   //180.0*CLHEP::deg, 360.0*CLHEP::deg);
    G4RotationMatrix* rottmp2 = new G4RotationMatrix();
    //rottmp2->rotateX(90.0*CLHEP::deg);
    rottmp2->rotateX(90.0*CLHEP::deg);
//    tcell_front
//      = new G4Tubs("tcell_front",
//		   0.0, outrad1,
//		   front_thickness/2.0,
//		   0.0, 360.0*CLHEP::deg);
    G4UnionSolid* tcell_front_tube_end
      = new G4UnionSolid("tcell_front_tube_end",
			 tcell_tube_end,
			 //tcell_tube,
			 tcell_front_sphere,
			 rottmp2,
			 G4ThreeVector( 0.0,
					0.0,
					-(tcell_tube_zlen+front_thickness)/2.0)
			 );
//    G4UnionSolid* tcell_front_tube
//      = new G4UnionSolid("tcell_front_tube",
//			 tcell_tube,
//			 tcell_front,
//			 noRot,
//			 G4ThreeVector( 0.0,
//					0.0,
//					-(tcell_tube_zlen+front_thickness)/2.0)
//			 );

    //target_cellLV =  new G4LogicalVolume( tcell_front_tube, Target_cell, "Target cell LV"); // removed the front cell wall
    //target_cellLV =  new G4LogicalVolume( tcell_tube_end, Target_cell, "Target cell LV"); // removed the front cell wall
    target_cellLV =  new G4LogicalVolume( tcell_front_tube_end, Target_cell, "Target cell LV");
    G4Colour colourTCell(0.8, 0.5, 1.0); //
    G4VisAttributes *TCellVisAtt = new G4VisAttributes(true, colourTCell);
    target_cellLV->SetVisAttributes(TCellVisAtt);
    
    gastar_main =
      new G4Tubs("gastar_main",
		 0.0, innrad, 
		 tcell_tube_zlen/2.0,
		 0.0, 360.0*CLHEP::deg);
    gastar_end =
      new G4Sphere("gastar_end",
		   0.0, innrad,
		   0.0*CLHEP::deg, 180.0*CLHEP::deg,
		   0.0*CLHEP::deg, 180.0*CLHEP::deg);
    G4UnionSolid* gastar_main_end
      = new G4UnionSolid("gastar_main_end",
			 gastar_main,
			 gastar_end,
			 rottmp1,
			 G4ThreeVector( 0.0 , 0.0 , tcell_tube_zlen/2.0 )
			 );
    gastar_front =
      new G4Sphere("gastar_front",
		   0.0, innrad*CLHEP::mm,
		   0.0*CLHEP::deg, 180.0*CLHEP::deg,
		   0.0*CLHEP::deg, 180.0*CLHEP::deg);
		   //180.0*CLHEP::deg, 360.0*CLHEP::deg,
    //180.0*CLHEP::deg, 270.0*CLHEP::deg);
    G4UnionSolid* gastar_main_end_front
      = new G4UnionSolid("gastar_main_end_front",
			 gastar_main_end,
			 gastar_front,
			 rottmp2,
			 G4ThreeVector( 0.0 , 0.0 , -tcell_tube_zlen/2.0 )
			 //G4ThreeVector( 0.0 , 0.0 , -tcell_tube_zlen/1.5 )
			 );
    
    targetLV =  new G4LogicalVolume( gastar_main_end_front, Target, "Target LV");
//    new G4PVPlacement( noRot 
//		       , G4ThreeVector(0.0*CLHEP::cm, 0.0*CLHEP::cm, 0.0*CLHEP::cm)
//		       , target_cellLV
//		       , "GasTarget" 
//		       , targetLV
//		       , false
//		       , 0 );
    
    //rotTarget->rotateZ( 90.0 * CLHEP::deg );
    //rotTarget->rotateZ( 270.0 * CLHEP::deg );
    //rotTarget->rotateX( 180.0 * CLHEP::deg );
    //rotTarget->rotateY( 180.0 * CLHEP::deg ); // front side back, Toshi, 29 April 2018

//    tcell_front
//      = new G4Tubs("tcell_front",
//		   0.0, outrad1,
//		   front_thickness/2.0,
//		   0.0, 360.0*CLHEP::deg);
    
  }
  
  else if( TargetMaterial == "Calcium40"){
    Target = mList_->Ca40;
  }
  else if( TargetMaterial == "Calcium48"){
    Target = mList_->Ca48;
  }
  else if( TargetMaterial == "Copper"){
    Target = mList_->Cu;
  }
  else if( TargetMaterial == "Water"){
    Target = mList_->Water;
    //Foil for water target
    //G4Material *Foil =mList_->Vacuum;
    G4Material *Foil =mList_->Harver;
    G4double foilz = 0.025 * 0.1 * CLHEP::cm ; // 25 [um]
    //G4double foilz = 0.25 * 0.1 * CLHEP::cm ; // 25 [um]
    G4Box *foilSolid 
      = new G4Box( "foilSolid_box",
		   TargetSizeX/2+2.0/2.0*CLHEP::mm,TargetSizeY/2,foilz/2 );
    //G4RotationMatrix* noRot = new G4RotationMatrix();
    G4DisplacedSolid* targetf1 
      = new G4DisplacedSolid("targetf1",
			     foilSolid,
			     noRot ,
			     G4ThreeVector( 0.0 , 0.0 , -TargetSizeZ/2.0 - foilz/2.0 )
			     );
    G4UnionSolid* targetf2 
      = new G4UnionSolid("targetf2",
			 targetf1,
			 foilSolid,
			 noRot,
			 G4ThreeVector( 0.0 , 0.0 , TargetSizeZ/2.0 + foilz/2.0 )
			 );
    G4LogicalVolume *foilLV = 
      //new G4LogicalVolume( foilSolid, Foil , "Foil LV");
      new G4LogicalVolume( targetf2, Foil , "Foil LV");
    //    G4double poszFoil_1  = -TargetSizeZ/2.0 - foilz/2.0 ;
    //    G4double poszFoil_2  = TargetSizeZ/2.0 + foilz/2.0 ;
    //    new G4PVPlacement( rotTarget
    //		       , G4ThreeVector(0*CLHEP::cm, 0*CLHEP::cm, poszFoil_1)
    //		       , foilLV
    //		       , "Foil" 
    //		       , worldLV 
    //		       , false
    //		       , 0 );
    //    new G4PVPlacement( rotTarget
    //		       , G4ThreeVector(0*CLHEP::cm, 0*CLHEP::cm, poszFoil_2)
    //		       , foilLV
    //		       , "Foil" 
    //		       , worldLV 
    //		       , false
    //		       , 1 );
    new G4PVPlacement( rotTarget
		       , G4ThreeVector(0.0*CLHEP::cm, 0.0*CLHEP::cm, 0.0*CLHEP::cm)
		       , foilLV
		       , "Foil" 
		       , worldLV 
		       , false
		       , 0 );
    //=== Visualization for foil of water target ===//
    //G4Colour colourFoil(1.0, 1.0, 1.0); // White
    G4Colour colourFoil(1.0, 0.8, 1.0); //
    G4VisAttributes *foilVisAtt = new G4VisAttributes(true, colourFoil);
    foilLV->SetVisAttributes(foilVisAtt);
  }
  else if( TargetMaterial == "Chromium"){
    Target = mList_->Cr52;
  }
  else if( TargetMaterial == "Lead"){
    Target = mList_->Pb208;
  }
  else if( TargetMaterial == "Tungsten"){
    Target = mList_->W;
  }
  else {
    G4cout << "You have chosen a target " << TargetMaterial 
	   << " which is not registered in this simulation." 
	   << G4endl;
    G4cout << " ..... So, vacuum is set to be used in this simulation." 
	   <<G4endl;
    Target = mList_->Vacuum;
  }
  //G4double TargetDensity = Target->GetDensity();
  //G4double TargetSizeX = 1.*cm;
  //G4double TargetSizeY = 1.*cm;
  //G4double TargetSizeZ = (TargetThickness)/(TargetDensity/(g/cm3)*1000)*cm;
  //G4RotationMatrix* rotTarget = new G4RotationMatrix();//Rotaion of target
  //rotTarget->rotateY(-17.0 * deg);//Rotation of target
  //G4LogicalVolume *targetLV; 

  if(TargetMaterial == "Tritium" || TargetMaterial == "H2" ){
    
    // ------- Scattering chamber ---------- //
    G4RotationMatrix* scat_chamb_rot = new G4RotationMatrix();
    scat_chamb_rot->rotateX(90.0*CLHEP::deg);
    new G4PVPlacement( scat_chamb_rot 
		       , G4ThreeVector(0*CLHEP::cm, 0*CLHEP::cm, 0.0*CLHEP::cm)
		       //, targetLV
		       , scat_chambLV
		       , "Scattering chamber" 
		       , worldLV 
		       , false
		       , 0 );
    // ------- Scattering chamber window ---------- //
    G4RotationMatrix* scat_chamb_window_rot = new G4RotationMatrix();
    scat_chamb_window_rot->rotateX(270.0*CLHEP::deg);
    new G4PVPlacement( scat_chamb_window_rot 
		       , G4ThreeVector(0*CLHEP::cm, 0*CLHEP::cm, 0.0*CLHEP::cm)
		       //, targetLV
		       , scat_chambwinLV
		       , "Scattering chamber" 
		       , worldLV 
		       , false
		       , 0 );

    // ------- Gass cell ---------- //
    G4RotationMatrix* gastar_rot = new G4RotationMatrix();
    //gastar_rot->rotateX(90.0*CLHEP::deg); // backside front 
    gastar_rot->rotateX(270.0*CLHEP::deg);
    new G4PVPlacement( gastar_rot
		       , G4ThreeVector(0*CLHEP::cm, 0*CLHEP::cm, 0.0* CLHEP::cm)
		       //, targetLV
		       , target_cellLV
		       , "Target cell" 
		       , scat_chambLV 
		       , false
		       , 0 );

    // ------- Gaseous target in the gass cell ---------- //
    new G4PVPlacement( 0
		       , G4ThreeVector(0*CLHEP::cm, 0*CLHEP::cm, 0.0 * CLHEP::cm)
		       //, targetLV
		       , targetLV
		       , "Target" 
		       , target_cellLV
		       , false
		       , 0 );

    // ------- Be isolator ---------- //
    new G4PVPlacement( gastar_rot
		       //, G4ThreeVector(0*CLHEP::cm, 0*CLHEP::cm, (-12.5-5.0)* CLHEP::cm)
		       , G4ThreeVector(0*CLHEP::cm, 40.0*CLHEP::cm, 0.0* CLHEP::cm)
		       //, targetLV
		       , be_isolatorLV 
		       , "Be isolator" 
		       , scat_chambLV
		       , false
		       , 0 );

    // ------- HRS windows ---------- //
    G4RotationMatrix* lhrs_window1_rot = new G4RotationMatrix();
    //G4double dist_t2window = 169.0*CLHEP::cm;
    //G4double lhrs_ang = 13.2*CLHEP::deg;
    lhrs_window1_rot->rotateY(-lhrs_ang);
    new G4PVPlacement( lhrs_window1_rot
		       , G4ThreeVector(dist_t2window*sin(lhrs_ang),
				       0.0*CLHEP::cm,
				       dist_t2window*cos(lhrs_ang))
		       , hrs_window1LV
		       , "LHRS window 1" 
		       , worldLV 
		       , false
		       , 0 );

    G4RotationMatrix* rhrs_window1_rot = new G4RotationMatrix();
    //G4double rhrs_ang = 13.2*CLHEP::deg;
    rhrs_window1_rot->rotateY(rhrs_ang);
    new G4PVPlacement( rhrs_window1_rot
		       , G4ThreeVector(-dist_t2window*sin(rhrs_ang),
				       0.0*CLHEP::cm,
				       dist_t2window*cos(rhrs_ang))
		       , hrs_window1LV
		       , "RHRS window 1" 
		       , worldLV 
		       , false
		       , 1 );

    // ------- HRS windows (before VDC) ---------- //
    G4double dist_t2window2 = dist_t2window + 0.5*CLHEP::cm;
    new G4PVPlacement( lhrs_window1_rot
		       , G4ThreeVector(dist_t2window2*sin(lhrs_ang),
				       0.0*CLHEP::cm,
				       dist_t2window2*cos(lhrs_ang))
		       , hrs_window2LV
		       , "LHRS window 1" 
		       , worldLV 
		       , false
		       , 0 );

    new G4PVPlacement( rhrs_window1_rot
		       , G4ThreeVector(-dist_t2window2*sin(rhrs_ang),
				       0.0*CLHEP::cm,
				       dist_t2window2*cos(rhrs_ang))
		       , hrs_window2LV
		       , "RHRS window 1" 
		       , worldLV 
		       , false
		       , 1 );
    
//    new G4PVPlacement( rotTarget
//		       , G4ThreeVector(0*CLHEP::cm, 0*CLHEP::cm, -12.7/2.0/2.0 * CLHEP::mm)
//		       //, targetLV
//		       , target_cellLV
//		       , "Target cell" 
//		       , worldLV 
//		       , false
//		       , 0 );
//    new G4PVPlacement( rotTarget
//		       , G4ThreeVector(0*CLHEP::cm, 0*CLHEP::cm, -12.7/2.0/2.0 * CLHEP::mm)
//		       //, targetLV
//		       , targetLV
//		       , "Target" 
//		       , worldLV 
//		       , false
//		       , 0 );
  }
  else if (TargetMaterial=="LHe4" || TargetMaterial=="LHe3"){
  }
  else if (TargetMaterial != "Vacuum"){
    G4Box *targetSolid = new G4Box( "Target",
				    TargetSizeX/2.,TargetSizeY/2.,TargetSizeZ/2. );
    targetLV =  new G4LogicalVolume( targetSolid, Target, "Target LV");
    new G4PVPlacement( rotTarget
		       , G4ThreeVector(0.0*CLHEP::cm, 0.0*CLHEP::cm, 0.0*CLHEP::mm)
		       , targetLV
		       , "Target" 
		       , worldLV 
		       , false
		       , 0 );
  }
  
///  if (TargetMaterial != "Vacuum"){
///    new G4PVPlacement( rotTarget
///		       , G4ThreeVector(0*CLHEP::cm, 0*CLHEP::cm, 0*CLHEP::cm)
///		       ,targetLV
///		       , "Target" 
///		       , worldLV 
///		       , false
///		       , 0 );
///  }
  
  //Virtual Detector
  G4ThreeVector vd1pos(0.0*CLHEP::cm , 0.0*CLHEP::cm , zpos_for_VD);
  //G4ThreeVector vd1pos(0.0*CLHEP::cm , 0.0*CLHEP::cm , 0.5*CLHEP::cm);
  //G4double vd1_x = 4.0*CLHEP::cm;   // [CLHEP::cm]
  //G4double vd1_y = 4.0*CLHEP::cm;   // [CLHEP::cm]
  //G4ThreeVector vd1pos(0.0*CLHEP::cm , 0.0*CLHEP::cm , 13.0*CLHEP::cm); // Tritium target
  G4double vd1_x = 10.0*CLHEP::cm;   // [cm]
  G4double vd1_y = 10.0*CLHEP::cm;   // [cm]
  G4double vd1_z = 0.001*CLHEP::mm; // [mm]
  //G4double vd1_z = 1.0*CLHEP::cm; // [mm]
  G4Box *vd1_box = new G4Box( "VD1_BOX",vd1_x/2.0 , vd1_y/2.0 , vd1_z/2.0 );
  G4LogicalVolume *vd1LV
    = new G4LogicalVolume( vd1_box , mList_->Vacuum , "VD1_LV");
  G4PVPlacement *vd1 
    = new G4PVPlacement( 0 ,      //rotaion
			 vd1pos , //position
			 vd1LV ,  //its logical volume
			 "VD1" ,  //its name
			 worldLV ,//its mother
			 false ,  //no boolean operation
			 0 );     //copy number
  G4SDManager *SDman = G4SDManager::GetSDMpointer();
  G4String vdSDname1 = "/TARG/VD1";
  TARGvdSD* vdSD1 = new TARGvdSD( vdSDname1 , 1 , "VDCollection1" );
  SDman -> AddNewDetector( vdSD1 );
  vd1LV -> SetSensitiveDetector( vdSD1 );


  // ----- Virtual detector for Left HRS ------- //
  G4double vd_x = 300.0*CLHEP::mm;
  G4double vd_y = 200.0*CLHEP::mm;
  G4double vd_z =  0.01*CLHEP::mm;
  G4double dist_t2vdleft = dist_t2window + 1.0*CLHEP::cm;
  G4double dist_t2vdright = dist_t2vdleft;

  G4RotationMatrix* lhrs_vd_rot = new G4RotationMatrix();
  lhrs_vd_rot->rotateY(-lhrs_ang);
  //G4Box *vd_left_box = new G4Box( "VD_LEFT_BOX", vd_x/2.0 , vd_y/2.0 , vd_z/2.0 );
  G4Tubs* vd_left_tube =  new G4Tubs("VD LEFT TUBS",
				      0.0, 15.0*CLHEP::cm,
				      0.01*CLHEP::mm,
				      0.0, 360.0*CLHEP::deg);
  G4LogicalVolume *vd_left_LV
    //= new G4LogicalVolume( vd_left_box , mList_->Vacuum , "VD_LEFT_LV");
    = new G4LogicalVolume( vd_left_tube , mList_->Vacuum , "VD_LEFT_LV");
  G4PVPlacement *vd_left 
    = new G4PVPlacement( lhrs_vd_rot,
			 G4ThreeVector(dist_t2vdleft*sin(lhrs_ang),
				       0.0*CLHEP::cm,
				       dist_t2vdleft*cos(lhrs_ang)),
			 vd_left_LV ,
			 "VD LEFT LV" ,
			 worldLV ,
			 false ,  
			 0 );
  //G4SDManager *SDman = G4SDManager::GetSDMpointer();
  G4String vdSDname2 = "/TARG/VD_LEFT";
  TARGvdSD* vdSD2 = new TARGvdSD( vdSDname2 , 2 , "VDCollection2" );
  SDman -> AddNewDetector( vdSD2 );
  vd_left_LV -> SetSensitiveDetector( vdSD2 );

  
  // ----- Virtual detector for Right HRS ------- //
  G4RotationMatrix* rhrs_vd_rot = new G4RotationMatrix();
  rhrs_vd_rot->rotateY(rhrs_ang);
  //G4Box *vd_right_box = new G4Box( "VD_RIGHT_BOX", vd_x/2.0 , vd_y/2.0 , vd_z/2.0 );
  G4Tubs* vd_right_tube =  new G4Tubs("VD RIGHT TUBS",
				      0.0, 15.0*CLHEP::cm,
				      0.01*CLHEP::mm,
				      0.0, 360.0*CLHEP::deg);
  G4LogicalVolume *vd_right_LV
    //= new G4LogicalVolume( vd_right_box , mList_->Vacuum , "VD_RIGHT_LV");
    = new G4LogicalVolume( vd_right_tube , mList_->Vacuum , "VD_RIGHT_LV");
  G4PVPlacement *vd_right
    = new G4PVPlacement( rhrs_vd_rot,
			 G4ThreeVector(-dist_t2vdright*sin(rhrs_ang),
				       0.0*CLHEP::cm,
				       dist_t2vdright*cos(rhrs_ang)),
			 vd_right_LV ,
			 "VD RIGHT LV" ,  //its name
			 worldLV ,//its mother
			 false ,  //no boolean operation
			 0 );     //copy number
  //G4SDManager *SDman = G4SDManager::GetSDMpointer();
  G4String vdSDname3 = "/TARG/VD_RIGHT";
  TARGvdSD* vdSD3 = new TARGvdSD( vdSDname3 , 3 , "VDCollection3" );
  SDman -> AddNewDetector( vdSD3 );
  vd_right_LV -> SetSensitiveDetector( vdSD3 );

  // --------------------------------------------- //
  // ----- Virtual detector for beam       ------- //
  // --------------------------------------------- //
  G4RotationMatrix* beam_vd_rot = new G4RotationMatrix();
  beam_vd_rot->rotateX(90*CLHEP::deg);
  //G4Box *vd_right_box = new G4Box( "VD_RIGHT_BOX", vd_x/2.0 , vd_y/2.0 , vd_z/2.0 );
  G4Tubs* vd_beam_tube =  new G4Tubs("VD BEAM TUBS",
				      0.0, 5.0*CLHEP::cm,
				      0.01*CLHEP::mm,
				      0.0, 360.0*CLHEP::deg);
  G4LogicalVolume *vd_beam_LV
    = new G4LogicalVolume( vd_beam_tube , mList_->Vacuum , "VD_BEAM_LV");
  G4PVPlacement *vd_beam
    = new G4PVPlacement( beam_vd_rot,
			 G4ThreeVector(0.0*CLHEP::cm,
				       42.0*CLHEP::cm,
				       0.0*CLHEP::cm),
			 vd_beam_LV ,
			 "VD BEAM LV" ,  //its name
			 //worldLV ,//its mother
			 scat_chambLV ,//its mother
			 false ,  //no boolean operation
			 0 );     //copy number
  //G4SDManager *SDman = G4SDManager::GetSDMpointer();
  G4String vdSDname4 = "/TARG/VD_BEAM";
  TARGvdSD* vdSD4 = new TARGvdSD( vdSDname4 , 4 , "VDCollection4" );
  SDman -> AddNewDetector( vdSD4 );
  vd_beam_LV -> SetSensitiveDetector( vdSD4 );
  
  //VD2
  //G4ThreeVector vd2pos(0.0*CLHEP::cm , 0.0*CLHEP::cm , 2.5*CLHEP::cm);
  //G4LogicalVolume *vd2LV
  //  = new G4LogicalVolume( vd1_box , mList_->Vacuum , "VD1_LV");
  //G4PVPlacement *vd2 
  //  = new G4PVPlacement( 0 ,      //rotaion
  //			 vd2pos , //position
  //			 vd2LV ,  //its logical volume
  //			 "VD2" ,  //its name
  //			 worldLV ,//its mother
  //			 false ,  //no boolean operation
  //			 0 );     //copy number
  ////G4SDManager *SDman = G4SDManager::GetSDMpointer();
  //G4String vdSDname2 = "/TARG/VD2";

  //TARGvdSD* vdSD2 = new TARGvdSD( vdSDname2 , 1 , "VDCollection2" );
  //SDman -> AddNewDetector( vdSD2 );
  //vd2LV -> SetSensitiveDetector( vdSD2 );
  
  ///////////////////
  // Visualisation //
  ///////////////////
  
  //=== World ===//
  worldLV->SetVisAttributes(G4VisAttributes::Invisible);
  
  //=== Target ===//
  G4Colour colourTarget(0., 1., 1.); // Cyan
  G4VisAttributes *targetVisAtt = new G4VisAttributes(true, colourTarget);
  targetLV->SetVisAttributes(targetVisAtt);

  G4Colour colourBeIso(0.3, 1.0, 0.5);
  G4VisAttributes *BeIsoVisAtt = new G4VisAttributes(true, colourBeIso);
  be_isolatorLV->SetVisAttributes(BeIsoVisAtt);

  G4Colour colourWindow1(0.8, 0.2, 0.2);
  G4VisAttributes *Window1VisAtt = new G4VisAttributes(true, colourWindow1);
  hrs_window1LV->SetVisAttributes(Window1VisAtt);

  G4Colour colourWindow2(0.2, 0.2, 0.8);
  G4VisAttributes *Window2VisAtt = new G4VisAttributes(true, colourWindow2);
  hrs_window2LV->SetVisAttributes(Window2VisAtt);

  
  //=== Virtual Detector ===//
  vd1LV->SetVisAttributes(G4VisAttributes::Invisible);
  G4Colour colourvd1(0.6, 0.5, 0.6); // 
  G4VisAttributes *vd1VisAtt = new G4VisAttributes(true, colourvd1);
  //vd1LV->SetVisAttributes(vd1VisAtt);

  G4Colour colourvd2(1.0, 0.5, 0.3);
  G4VisAttributes *vd2VisAtt = new G4VisAttributes(true, colourvd2);
  vd_left_LV->SetVisAttributes(vd2VisAtt);
  //vd_left_LV->SetVisAttributes(G4VisAttributes::Invisible);

  G4Colour colourvd3(1.0, 0.5, 0.3);
  G4VisAttributes *vd3VisAtt = new G4VisAttributes(true, colourvd3);
  vd_right_LV->SetVisAttributes(vd3VisAtt);
  //vd_right_LV->SetVisAttributes(G4VisAttributes::Invisible);

  G4Colour colourvd4(1.0, 0.5, 0.3);
  G4VisAttributes *vd4VisAtt = new G4VisAttributes(true, colourvd4);
  vd_beam_LV->SetVisAttributes(vd4VisAtt);
  //vd_beam_LV->SetVisAttributes(G4VisAttributes::Invisible);
  
  //vd1LV->SetVisAttributes(vd1VisAtt);

  return world;
}

