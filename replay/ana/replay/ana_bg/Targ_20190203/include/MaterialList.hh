/*
  "MaterialList.hh"
*/

#ifndef MaterialList_h

#define MaterialList_h 1

#include "G4Element.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"

class MaterialList
{
public:
  MaterialList();
  ~MaterialList();

  G4Element *elH;
  G4Element *el3H;
  G4Element *elLi;
  G4Element *elBe;
  G4Element *elB;
  G4Element *elC;
  G4Element *elN;
  G4Element *elO;
  G4Element *elMg;
  G4Element *elAl;
  G4Element *elSi;
  G4Element *elAr;
  G4Element *elCa;
  G4Element *elCr;
  G4Element *elMn;
  G4Element *elFe;
  G4Element *elCo;
  G4Element *elNi;
  G4Element *elCu;
  G4Element *elZn;
  G4Element *elMo;
  G4Element *elW;

  G4Material *LHe3;
  G4Material *LHe4;
  G4Material *HeGas;
  G4Material *H3Gas;
  G4Material *HGas;
  G4Material *ArGas;
  G4Material *Mylar;
  G4Material *Li3N;
  G4Material *Kevlar;
  G4Material *Al;
  G4Material *W;
  G4Material *Fe;
  G4Material *Cu;
  G4Material *C12;
  G4Material *Si28;
  G4Material *Ti48;
  G4Material *V51;
  G4Material *Cr52;
  G4Material *Y89;
  G4Material *Li7;
  G4Material *Be9;
  G4Material *B10;
  G4Material *Ca40;
  G4Material *Ca48;
  G4Material *Pb208;
  G4Material *N2Gas;
  G4Material *O2Gas;

  G4Material *Scinti;
  G4Material *AC;
  G4Material *Water;
  G4Material *C2H6;
  G4Material *CH2;
  G4Material *Al7075;

  G4Material *Vacuum;
  G4Material *Air;
  G4Material *DCGas;
  G4Material *Heavymet;
  G4Material *Harver;

private:
  MaterialList( const MaterialList & );
  MaterialList & operator = ( const MaterialList & );
};

#endif
