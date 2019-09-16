/*
  "MaterialList.cc"
*/

#include "MaterialList.hh"

MaterialList::MaterialList()
{
  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel, nAtoms;
  G4double fractionMass;

  // =========================== //
  // ======== Elements ========= //
  // =========================== //
  name = "Hydrogen"; symbol = "H";
  a = 1.00794*CLHEP::g/CLHEP::mole;
  elH = new G4Element( name, symbol, iz=1., a);
  
  name = "Tritium"; symbol = "3H";
  a = 3.0*CLHEP::g/CLHEP::mole;
  el3H = new G4Element( name, symbol, iz=1., a);

  name = "Lithium"; symbol = "Li";
  a = 6.941*CLHEP::g/CLHEP::mole;
  elLi = new G4Element( name, symbol, iz=3., a);
  
  name = "Beryllium"; symbol = "Be";
  a = 9.012182*CLHEP::g/CLHEP::mole;
  elBe = new G4Element( name, symbol, iz=4., a);
  
  name = "Boron"; symbol = "B";
  a = 10.811*CLHEP::g/CLHEP::mole;
  elB = new G4Element( name, symbol, iz=5., a);

  name = "Carbon"; symbol = "C";
  a = 12.011*CLHEP::g/CLHEP::mole;
  elC = new G4Element( name, symbol, iz=6., a);

  name = "Nitrogen"; symbol = "N";
  a = 14.00674*CLHEP::g/CLHEP::mole;
  elN = new G4Element( name, symbol, iz=7., a);

  name = "Oxygen"; symbol = "O";
  a = 15.9994*CLHEP::g/CLHEP::mole;
  elO = new G4Element( name, symbol, iz=8., a);

  name = "Magnesium"; symbol = "Mg";
  a = 24.304*CLHEP::g/CLHEP::mole;
  elMg = new G4Element( name, symbol, iz=12., a); // new

  name = "Alminium"; symbol = "Al";
  a = 26.9815385*CLHEP::g/CLHEP::mole;
  elAl = new G4Element( name, symbol, iz=13., a); // new
  
  name = "Silicon"; symbol = "Si";
  a = 28.0855*CLHEP::g/CLHEP::mole;
  elSi = new G4Element( name, symbol, iz=14., a);

  name = "Argon"; symbol = "Ar";
  a = 39.948*CLHEP::g/CLHEP::mole;
  elAr = new G4Element( name, symbol, iz=18., a);
  
  name = "Calcium"; symbol = "Ca";
  a = 40.078*CLHEP::g/CLHEP::mole;
  elCa = new G4Element( name, symbol, iz=20., a);
  
  name = "Chromium"; symbol = "Cr";
  a = 51.996*CLHEP::g/CLHEP::mole;
  elCr = new G4Element( name, symbol, iz=24., a);

  name = "Manganese"; symbol = "Mn";
  a = 54.938*CLHEP::g/CLHEP::mole;
  elMn = new G4Element( name, symbol, iz=25., a);

  name = "Iron"; symbol = "Fe";
  a = 55.845*CLHEP::g/CLHEP::mole;
  elFe = new G4Element( name, symbol, iz=26., a);

  name = "Cobalt"; symbol = "Co";
  a = 58.933*CLHEP::g/CLHEP::mole;
  elCo = new G4Element( name, symbol, iz=27., a);

  name = "Nickel"; symbol = "Ni";
  a = 58.69*CLHEP::g/CLHEP::mole;
  elNi = new G4Element( name, symbol, iz=28., a);

  name = "Copper"; symbol = "Cu";
  a = 63.546*CLHEP::g/CLHEP::mole;
  elCu= new G4Element( name, symbol, iz=29., a);

  name = "Zinc"; symbol = "Zn";
  a = 65.38*CLHEP::g/CLHEP::mole;
  elZn = new G4Element( name, symbol, iz=30., a); // new

  name = "Molybd"; symbol = "Mo";
  a = 95.94*CLHEP::g/CLHEP::mole;
  elMo = new G4Element( name, symbol, iz=42., a);

  name = "Tungsten"; symbol = "W";
  a = 183.84*CLHEP::g/CLHEP::mole;
  elW = new G4Element( name, symbol, iz=74., a );

  //=== Simple Materials ==//
  name = "Tungsten";
  a = 183.84*CLHEP::g/CLHEP::mole;
  density = 19.3*CLHEP::g/CLHEP::cm3;
  W = new G4Material( name, z=74., a, density);
  
  name = "Aluminuium";
  a = 26.98*CLHEP::g/CLHEP::mole;
  density = 2.70*CLHEP::g/CLHEP::cm3;
  Al = new G4Material( name, z=13., a, density);

  name = "Iron";
  a = 55.843*CLHEP::g/CLHEP::mole;
  density = 7.86*CLHEP::g/CLHEP::cm3;
  Fe = new G4Material( name, z=26., a, density);

  name = "Copper";
  a = 63.546*CLHEP::g/CLHEP::mole;
  density = 8.93*CLHEP::g/CLHEP::cm3;
  Cu = new G4Material( name, z=29., a, density);

  name = "Carbon";
  a = 12.0*CLHEP::g/CLHEP::mole;
  //density = 2.25*CLHEP::g/CLHEP::cm3; // original 
  density = 1.75*CLHEP::g/CLHEP::cm3; // real measurement in E05-115, Toshi(16Dec2013)
  C12 = new G4Material( name, z=6., a, density);

  name = "Silicon";
  a = 28.0*CLHEP::g/CLHEP::mole;
  density = 2.33*CLHEP::g/CLHEP::cm3;
  Si28 = new G4Material( name, z=14., a, density);
  
  name = "Titanium";
  a = 47.867*CLHEP::g/CLHEP::mole;
  density = 4.54*CLHEP::g/CLHEP::cm3;
  Ti48 = new G4Material( name, z=22., a, density);

  name = "Vanadium";
  a = 51.0*CLHEP::g/CLHEP::mole;
  density = 6.11*CLHEP::g/CLHEP::cm3;
  V51 = new G4Material( name, z=23., a, density);

  name = "Cromium";
  a = 52.0*CLHEP::g/CLHEP::mole;
  density = 7.15*CLHEP::g/CLHEP::cm3;
  Cr52 = new G4Material( name, z=24., a, density);

  name = "Yttrium";
  a = 89.0*CLHEP::g/CLHEP::mole;
  density = 4.469*CLHEP::g/CLHEP::cm3;
  Y89 = new G4Material( name, z=39., a, density);
  
  name = "Lithium";
  a = 7.0*CLHEP::g/CLHEP::mole;
  density = 0.534*CLHEP::g/CLHEP::cm3;
  Li7 = new G4Material( name, z=3., a, density);
  
  name = "Berylium";
  a = 9.012*CLHEP::g/CLHEP::mole;
  density = 1.848*CLHEP::g/CLHEP::cm3;
  Be9 = new G4Material( name, z=4., a, density);
  
  name = "Boron";
  a = 10.0*CLHEP::g/CLHEP::mole;
  density = 2.340*CLHEP::g/CLHEP::cm3;
  B10 = new G4Material( name, z=5., a, density);
  
  name = "Calcium40";
  a = 40.0*CLHEP::g/CLHEP::mole;
  density = 1.55*CLHEP::g/CLHEP::cm3;
  Ca40 = new G4Material( name, z=20., a, density);
  
  name = "Calcium48";
  a = 48.0*CLHEP::g/CLHEP::mole;
  density = 1.55*CLHEP::g/CLHEP::cm3;
  Ca48 = new G4Material( name, z=20., a, density);

  name = "Lead";
  a = 208.0*CLHEP::g/CLHEP::mole;
  density = 11.342*CLHEP::g/CLHEP::cm3;
  Pb208 = new G4Material( name, z=82., a, density);

  name = "LHe3";
  a = 3.0*CLHEP::g/CLHEP::mole;
  density = 70.0*CLHEP::mg/CLHEP::cm3;
  LHe3 = new G4Material( name, z=2., a, density );
  
  name = "LHe4";
  a = 4.002602*CLHEP::g/CLHEP::mole;
  density = 130.0*CLHEP::mg/CLHEP::cm3;
  LHe4 = new G4Material( name, z=2., a, density );
  
  name = "HeliumGas";
  a = 4.002602*CLHEP::g/CLHEP::mole;
  density = 0.1786*CLHEP::mg/CLHEP::cm3;
  HeGas = new G4Material( name, z=2., a, density );
  
  name = "TritiumGas";
  a = 3.0*CLHEP::g/CLHEP::mole;
  density = 3.16*CLHEP::mg/CLHEP::cm3;
  H3Gas = new G4Material( name, z=1., a, density );

  name = "HydrogenGas";
  a = 1.0*CLHEP::g/CLHEP::mole;
  density = 2.832*CLHEP::mg/CLHEP::cm3;
  HGas = new G4Material( name, z=1., a, density );

  name = "Argon Gas";
  a = 39.948*CLHEP::g/CLHEP::mole;
  density = 0.001782*CLHEP::g/CLHEP::cm3;
  ArGas = new G4Material( name, z=18., a, density );

  name = "Nitrogen Gas";
  a = 14.01*CLHEP::g/CLHEP::mole;
  N2Gas = new G4Material( name, z=7., a, density );

  name = "Oxygen Gas";
  a = 16.00*CLHEP::g/CLHEP::mole;
  density = 1.429*CLHEP::mg/CLHEP::cm3;
  O2Gas = new G4Material( name, z=8., a, density );

  //=== Componds Materials ===//
  name = "Mylar";
  density = 1.39*CLHEP::g/CLHEP::cm3;
  Mylar = new G4Material( name, density, nel=3 );
  Mylar->AddElement( elH, nAtoms=8 );
  Mylar->AddElement( elC, nAtoms=10 );
  Mylar->AddElement( elO, nAtoms=4 );
  
  //=== Componds Materials ===//
  name = "Li3N";
  density = 1.27*CLHEP::g/CLHEP::cm3;
  Li3N = new G4Material( name, density, nel=2 );
  Li3N->AddElement( elLi, nAtoms=3 );
  Li3N->AddElement( elN,  nAtoms=1 );
  
  name = "Kevlar";
  density = 0.74*CLHEP::g/CLHEP::cm3;
  Kevlar = new G4Material( name, density, nel=4 );
  Kevlar->AddElement( elH, nAtoms=10 );
  Kevlar->AddElement( elC, nAtoms=14 );
  Kevlar->AddElement( elO, nAtoms=2 );
  Kevlar->AddElement( elN, nAtoms=2 );

  name = "Scintillator";
  density = 1.032*CLHEP::g/CLHEP::cm3;
  Scinti = new G4Material( name, density, nel=2);
  Scinti->AddElement(elH, nAtoms=8);
  Scinti->AddElement(elC, nAtoms=8);

  name = "Aerogel";
  density = 0.06*CLHEP::g/CLHEP::cm3;
  AC = new G4Material( name, density, nel=2);
  AC->AddElement(elSi, nAtoms=1);
  AC->AddElement(elO,  nAtoms=2);

  name = "Water";
  density = 1.*CLHEP::g/CLHEP::cm3;
  Water = new G4Material( name, density, nel=2);
  Water->AddElement(elH, nAtoms=1);
  Water->AddElement(elO, nAtoms=2);

  name = "Etane";
  density = 1.356*CLHEP::mg/CLHEP::cm3;
  C2H6 = new G4Material( name, density, nel=2);
  C2H6->AddElement(elC, nAtoms=2);
  C2H6->AddElement(elH, nAtoms=6);
  
  name = "Polyethylene";
  density = 0.93*CLHEP::g/CLHEP::cm3;
  CH2 = new G4Material( name, density, nel=2);
  CH2->AddElement(elC, nAtoms=1);
  CH2->AddElement(elH, nAtoms=2);

  name = "Al7075";
  density = 2.81*CLHEP::g/CLHEP::cm3;
  Al7075 = new G4Material( name, density, nel=5);
  Al7075->AddElement(elAl, fractionMass=0.900);
  Al7075->AddElement(elZn, fractionMass=0.055);
  Al7075->AddElement(elMg, fractionMass=0.025);
  Al7075->AddElement(elCu, fractionMass=0.018);
  Al7075->AddElement(elCr, fractionMass=0.002);
  
  density= CLHEP::universe_mean_density;
  Vacuum= new G4Material(name="Vacuum", z=1, a=1.01*CLHEP::g/CLHEP::mole,
			 density, kStateGas, 300*CLHEP::kelvin, 3.e-18*CLHEP::pascal);
  //Vacuum= new G4Material(name="Vacuum", density, nel=2);
  //   Vacuum-> AddElement(elN, .7);
  //   Vacuum-> AddElement(elO, .3);
  //Vacuum->AddMaterial( N2Gas, .78 );
  //Vacuum->AddMaterial( O2Gas, .22 );
  
  name = "Air";
  density = 1.293*CLHEP::mg/CLHEP::cm3;
  Air = new G4Material( name, density, nel=2);
  //   Air->AddElement(elN, .7);
  //   Air->AddElement(elO, .3);
  Air->AddMaterial( N2Gas, .78 );
  Air->AddMaterial( O2Gas, .22 );
  
  name = "Ar+C2H6";
  density = (0.5*39.948+0.5*(2*12.011+6*1.00794))/22.3*CLHEP::mg/CLHEP::cm3;
  DCGas = new G4Material( name, density, nel=2);
  fractionMass = 39.948/(39.948+2*12.011+6*1.00794);
  DCGas->AddMaterial(ArGas, fractionMass);
  fractionMass = (2*12.011+6*1.00794)/(39.948+2*12.011+6*1.00794);
  DCGas->AddMaterial(C2H6, fractionMass);

  /*name = "EDCLayer";
  density = 1.39*CLHEP::g/CLHEP::cm3;
  EDCLayer = new G4Material( name, density, nel=3 );
  EDCLayer->AddElement( elH, nAtoms=8 );
  EDCLayer->AddElement( elC, nAtoms=10 );
  EDCLayer->AddElement( elO, nAtoms=4 );*/

  name = "Heavy Metal";
  density = 16.9*CLHEP::g/CLHEP::cm3;
  Heavymet = new G4Material( name, density, nel=3 );
  Heavymet->AddElement( elW, .9 );
  Heavymet->AddElement(elNi, .06 );
  Heavymet->AddElement(elCu, .04 );

  name = "Harver";
  density = 8.3*CLHEP::g/CLHEP::cm3;
  Harver = new G4Material( name, density, nel=7 );
  Harver->AddElement(elCr, .200 );
  Harver->AddElement(elMn, .016 );
  Harver->AddElement(elFe, .181 );
  Harver->AddElement(elCo, .425 );
  Harver->AddElement(elNi, .130 );
  Harver->AddElement(elMo, .020 );
  Harver->AddElement(elW , .028 );
}

MaterialList::~MaterialList()
{
  delete Harver;
  delete Heavymet;
  delete DCGas;
  delete Air;
  delete Vacuum;

  delete C2H6;
  delete CH2;
  delete Water;
  delete AC;
  delete Scinti;

  delete Li7;
  delete Be9;
  delete B10;
  delete C12;
  delete Si28;
  delete Ti48;
  delete V51;
  delete Ca40;
  delete Ca48;
  delete Cr52;
  delete Cu;
  delete Fe;
  delete Al;
  delete Y89;
  delete Pb208;
  delete W;
  delete Mylar;
  delete HeGas;
  delete H3Gas;
  delete HGas;
  delete LHe3;
  delete LHe4;
  delete Kevlar;
  delete ArGas;
  delete N2Gas;
  delete O2Gas;
  delete Al7075;

  delete elH;
  delete el3H;
  delete elLi;
  delete elBe;
  delete elB;
  delete elC;
  delete elN;
  delete elO;
  delete elSi;
  delete elAr;
  delete elCa;
  delete elCr;
  delete elMn;
  delete elFe;
  delete elCo;
  delete elNi;
  delete elCu;
  delete elMo;
  delete elW;
}
