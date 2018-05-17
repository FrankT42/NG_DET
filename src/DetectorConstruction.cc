

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Polycone.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction()
{
  fExpHall_x = fExpHall_y = fExpHall_z = 10*m;
  fGuide_x    = fGuide_y    = fGuide_z    =  5.0*m;
  fScint_x  = fScint_y  = fScint_z  =  0.5*mm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{

// ------------- Materials -------------

  G4double a, z, density;
  G4int nelements;
  G4NistManager* man = G4NistManager::Instance();

// Air
//
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  G4Element* C = new G4Element("Carbon"  , "C",  z=6 , a=12.011*g/mole);
  G4Element* H = new G4Element("Hydrogen" , "H", z=1, a=1.008*g/mole);
  G4Element* Cs = new G4Element("Caesium",  "Cs", z=55, a=132.905*g/mole);
  G4Element* Li = new G4Element("Lithium",  "Li", z=3,  a=6.94*g/mole);
  G4Element* Y  = new G4Element("Yttrium",  "Y",  z=39, a=88.905*g/mole);
  G4Element* Cl = new G4Element("Chlorine", "Cl", z=17, a=35.45*g/mole);
  G4Element* Ce = new G4Element("Cerium",   "Ce", z=58, a=140.116*g/mole);
  G4Element* Li_6 = new G4Element("Lithium", "Li", z=3, a=6.*g/mole); 

  G4Material* air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  air->AddElement(N, 70.*perCent);
  air->AddElement(O, 30.*perCent);

// Water
//
  //G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);

  G4Material* water = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
  water->AddElement(H, 2);
  water->AddElement(O, 1);
// perspex
//
  //G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);

  G4Material* perspex = new G4Material("Perspex", density= 1.18*g/cm3, nelements=3);
  perspex->AddElement(H, 8);
  perspex->AddElement(O, 2);
  perspex->AddElement(C, 5);

  G4Material* Ej228 = new G4Material("Ej228", density = 1.023*g/cm3, nelements = 2);
  Ej228->AddElement(H,10);
  Ej228->AddElement(C,4);

    G4Material* SiPM_mat = new G4Material("SiPM_mat", density= 1.18*g/cm3, nelements=3);
  SiPM_mat->AddElement(H, 8);
  SiPM_mat->AddElement(O, 2);
  SiPM_mat->AddElement(C, 5);
//wrapping, density about perspex for now
  G4Material* wrapping = new G4Material("wrapping", density= 1.18*g/cm3, nelements=3);
  wrapping->AddElement(H, 8);
  wrapping->AddElement(O, 2);
  wrapping->AddElement(C, 5);
// CLYC
//

  G4Material* CLYC = new G4Material("CLYC", density= 3.31*g/cm3, nelements=5);
  CLYC->AddElement(Cs, 2);
  CLYC->AddElement(Li_6, 1);
  CLYC->AddElement(Y, 1);
  CLYC->AddElement(Cl,6);
  CLYC->AddElement(Ce,1);

  
  G4Material* quartz = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4Material* Al = man->FindOrBuildMaterial("G4_Al");
//
// ------------ Generate & Add Material Properties Table ------------
//
  G4double photonEnergy[] =
            { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
              2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
              2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
              2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
              2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
              3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
              3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
              3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };

  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);

//
// Water
//
  G4double refractiveIndex1[] =
            { 1.3435, 1.344,  1.3445, 1.345,  1.3455,
              1.346,  1.3465, 1.347,  1.3475, 1.348,
              1.3485, 1.3492, 1.35,   1.3505, 1.351,
              1.3518, 1.3522, 1.3530, 1.3535, 1.354,
              1.3545, 1.355,  1.3555, 1.356,  1.3568,
              1.3572, 1.358,  1.3585, 1.359,  1.3595,
              1.36,   1.3608};

  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));

  G4double absorption[] =
           {3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m,
           15.152*m, 17.241*m, 18.868*m, 20.000*m, 26.316*m, 35.714*m,
           45.455*m, 47.619*m, 52.632*m, 52.632*m, 55.556*m, 52.632*m,
           52.632*m, 47.619*m, 45.455*m, 41.667*m, 37.037*m, 33.333*m,
           30.000*m, 28.500*m, 27.000*m, 24.500*m, 22.000*m, 19.500*m,
           17.500*m, 14.500*m };

  assert(sizeof(absorption) == sizeof(photonEnergy));

  G4double scintilFast[] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };

  assert(sizeof(scintilFast) == sizeof(photonEnergy));

  G4double scintilSlow[] =
            { 0.01, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00,
              7.00, 8.00, 9.00, 8.00, 7.00, 6.00, 4.00,
              3.00, 2.00, 1.00, 0.01, 1.00, 2.00, 3.00,
              4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 8.00,
              7.00, 6.00, 5.00, 4.00 };

  
  
assert(sizeof(scintilSlow) == sizeof(photonEnergy));

  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX",       photonEnergy, refractiveIndex1,nEntries)
        ->SetSpline(true);
  myMPT1->AddProperty("ABSLENGTH",    photonEnergy, absorption,     nEntries)
        ->SetSpline(true);
  myMPT1->AddProperty("FASTCOMPONENT",photonEnergy, scintilFast,     nEntries)
        ->SetSpline(true);
  myMPT1->AddProperty("SLOWCOMPONENT",photonEnergy, scintilSlow,     nEntries)
        ->SetSpline(true);

  myMPT1->AddConstProperty("SCINTILLATIONYIELD",50./MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT1->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  myMPT1->AddConstProperty("SLOWTIMECONSTANT",1000.*ns);
  myMPT1->AddConstProperty("YIELDRATIO",0.91);

  G4double energy_water[] = {
     1.56962*eV, 1.58974*eV, 1.61039*eV, 1.63157*eV,
     1.65333*eV, 1.67567*eV, 1.69863*eV, 1.72222*eV,
     1.74647*eV, 1.77142*eV, 1.7971 *eV, 1.82352*eV,
     1.85074*eV, 1.87878*eV, 1.90769*eV, 1.93749*eV,
     1.96825*eV, 1.99999*eV, 2.03278*eV, 2.06666*eV,
     2.10169*eV, 2.13793*eV, 2.17543*eV, 2.21428*eV,
     2.25454*eV, 2.29629*eV, 2.33962*eV, 2.38461*eV,
     2.43137*eV, 2.47999*eV, 2.53061*eV, 2.58333*eV,
     2.63829*eV, 2.69565*eV, 2.75555*eV, 2.81817*eV,
     2.88371*eV, 2.95237*eV, 3.02438*eV, 3.09999*eV,
     3.17948*eV, 3.26315*eV, 3.35134*eV, 3.44444*eV,
     3.54285*eV, 3.64705*eV, 3.75757*eV, 3.87499*eV,
     3.99999*eV, 4.13332*eV, 4.27585*eV, 4.42856*eV,
     4.59258*eV, 4.76922*eV, 4.95999*eV, 5.16665*eV,
     5.39129*eV, 5.63635*eV, 5.90475*eV, 6.19998*eV
  };

  const G4int numentries_water = sizeof(energy_water)/sizeof(G4double);

  //assume 100 times larger than the rayleigh scattering for now.
  G4double mie_water[] = {
     167024.4*m, 158726.7*m, 150742  *m,
     143062.5*m, 135680.2*m, 128587.4*m,
     121776.3*m, 115239.5*m, 108969.5*m,
     102958.8*m, 97200.35*m, 91686.86*m,
     86411.33*m, 81366.79*m, 76546.42*m,
     71943.46*m, 67551.29*m, 63363.36*m,
     59373.25*m, 55574.61*m, 51961.24*m,
     48527.00*m, 45265.87*m, 42171.94*m,
     39239.39*m, 36462.50*m, 33835.68*m,
     31353.41*m, 29010.30*m, 26801.03*m,
     24720.42*m, 22763.36*m, 20924.88*m,
     19200.07*m, 17584.16*m, 16072.45*m,
     14660.38*m, 13343.46*m, 12117.33*m,
     10977.70*m, 9920.416*m, 8941.407*m,
     8036.711*m, 7202.470*m, 6434.927*m,
     5730.429*m, 5085.425*m, 4496.467*m,
     3960.210*m, 3473.413*m, 3032.937*m,
     2635.746*m, 2278.907*m, 1959.588*m,
     1675.064*m, 1422.710*m, 1200.004*m,
     1004.528*m, 833.9666*m, 686.1063*m
  };

  assert(sizeof(mie_water) == sizeof(energy_water));

  // gforward, gbackward, forward backward ratio
  G4double mie_water_const[3]={0.99,0.99,0.8};

  myMPT1->AddProperty("MIEHG",energy_water,mie_water,numentries_water)
        ->SetSpline(true);
  myMPT1->AddConstProperty("MIEHG_FORWARD",mie_water_const[0]);
  myMPT1->AddConstProperty("MIEHG_BACKWARD",mie_water_const[1]);
  myMPT1->AddConstProperty("MIEHG_FORWARD_RATIO",mie_water_const[2]);

  G4cout << "Water G4MaterialPropertiesTable" << G4endl;
  myMPT1->DumpTable();

  water->SetMaterialPropertiesTable(myMPT1);

  // Set the Birks Constant for the Water scintillator

  water->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

//
// Air
//
  G4double refractiveIndex2[] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };

  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", photonEnergy, refractiveIndex2, nEntries);

  G4cout << "Air G4MaterialPropertiesTable" << G4endl;
  myMPT2->DumpTable();

  air->SetMaterialPropertiesTable(myMPT2);

//
// Perspex
//
  G4double refractiveIndex3[] =
            { 1.486, 1.486, 1.486, 1.486, 1.486, 1.486, 1.486,
              1.486, 1.486, 1.486, 1.486, 1.486, 1.486, 1.486,
              1.486, 1.486, 1.486, 1.486, 1.486, 1.486, 1.486,
              1.486, 1.486, 1.486, 1.486, 1.486, 1.486, 1.486,
              1.486, 1.486, 1.486, 1.486 };

  G4MaterialPropertiesTable* myMPT3 = new G4MaterialPropertiesTable();
  myMPT3->AddProperty("RINDEX", photonEnergy, refractiveIndex3, nEntries);

  G4cout << "PERRSPEX G4MaterialPropertiesTable" << G4endl;
  myMPT3->DumpTable();

  perspex->SetMaterialPropertiesTable(myMPT3);
  SiPM_mat->SetMaterialPropertiesTable(myMPT3);


  // Ej228
  //
  G4MaterialPropertiesTable* MPT_Ej228 = new G4MaterialPropertiesTable();
  G4double refractiveIndexEj[sizeof(photonEnergy)/sizeof(double)];
   G4double ScintSlowEj[sizeof(photonEnergy)/sizeof(double)];
    G4double ScintFastEj[sizeof(photonEnergy)/sizeof(double)];
  for(int i = 0; i<sizeof(photonEnergy)/sizeof(double);i++){
    refractiveIndexEj[i] = 1.58;
    ScintSlowEj[i] = 1;
    ScintFastEj[i] = 1;
  }

  
  MPT_Ej228->AddProperty("RINDEX",       photonEnergy, refractiveIndexEj,nEntries)
        ->SetSpline(true);
  MPT_Ej228->AddProperty("FASTCOMPONENT",photonEnergy, ScintFastEj,     nEntries)
        ->SetSpline(true);
  MPT_Ej228->AddProperty("SLOWCOMPONENT",photonEnergy, ScintSlowEj,     nEntries)
        ->SetSpline(true);

  MPT_Ej228->AddConstProperty("SCINTILLATIONYIELD",10200./MeV);
  MPT_Ej228->AddConstProperty("RESOLUTIONSCALE",1.0);
  MPT_Ej228->AddConstProperty("FASTTIMECONSTANT", 0.5*ns);
  MPT_Ej228->AddConstProperty("SLOWTIMECONSTANT",1.4*ns);
  MPT_Ej228->AddConstProperty("YIELDRATIO",1);
  Ej228->SetMaterialPropertiesTable(MPT_Ej228);

  
  
//
// CLYC
//
    G4MaterialPropertiesTable* myMPT4 = new G4MaterialPropertiesTable();
  G4double refractiveIndex4[] =
            { 1.81, 1.81,  1.81, 1.81,  1.81,
              1.81,  1.81, 1.81,  1.81, 1.81,
              1.81, 1.81, 1.81,   1.81, 1.81,
              1.81, 1.81, 1.81, 1.81, 1.81,
              1.81, 1.81,  1.81, 1.81,  1.81,
              1.81, 1.81,  1.81, 1.81,  1.81,
              1.81,   1.81};

  assert(sizeof(refractiveIndex4) == sizeof(photonEnergy));


  G4double scintilFast4[] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };

  assert(sizeof(scintilFast4) == sizeof(photonEnergy));

  /*  G4double scintilSlow4[] =
            { 0.01, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00,
              7.00, 8.00, 9.00, 8.00, 7.00, 6.00, 4.00,
              3.00, 2.00, 1.00, 0.01, 1.00, 2.00, 3.00,
              4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 8.00,
              7.00, 6.00, 5.00, 4.00 };
  */
  
  G4double scintilSlow4[] =
            { 1., 1., 1., 1., 1., 1., 1.,
              1., 1., 1., 1., 1., 1., 1.,
              1., 1., 1., 1., 1., 1., 1.,
              1., 1., 1., 1., 1., 1., 1.,
              1., 1., 1., 1. };
  assert(sizeof(scintilSlow4) == sizeof(photonEnergy));

  //G4MaterialPropertiesTable* myMPT4 = new G4MaterialPropertiesTable();

  myMPT4->AddProperty("RINDEX",       photonEnergy, refractiveIndex4,nEntries)
        ->SetSpline(true);
  myMPT4->AddProperty("FASTCOMPONENT",photonEnergy, scintilFast4,     nEntries)
        ->SetSpline(true);
  myMPT4->AddProperty("SLOWCOMPONENT",photonEnergy, scintilSlow4,     nEntries)
        ->SetSpline(true);

  myMPT4->AddConstProperty("SCINTILLATIONYIELD",20000./MeV);
  myMPT4->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT4->AddConstProperty("FASTTIMECONSTANT", 50.*ns);
  myMPT4->AddConstProperty("SLOWTIMECONSTANT",1000.*ns);
  myMPT4->AddConstProperty("YIELDRATIO",0.91);
  CLYC->SetMaterialPropertiesTable(myMPT4);
  //remember to change back to CLYC times 1ns and 1000ns and 20000/Mev
  
  //quartz
  G4double refractiveIndex5[] =
            { 1.46, 1.46,  1.46, 1.46,  1.46,
              1.46,  1.46, 1.46,  1.46, 1.46,
              1.46, 1.46, 1.46,   1.46, 1.46,
              1.46, 1.46, 1.46, 1.46, 1.46,
              1.46, 1.46,  1.46, 1.46,  1.46,
              1.46, 1.46,  1.46, 1.46,  1.46,
              1.46,   1.46};

  assert(sizeof(refractiveIndex4) == sizeof(photonEnergy));
  G4MaterialPropertiesTable* myMPT5 = new G4MaterialPropertiesTable();
  myMPT5->AddProperty("RINDEX",photonEnergy, refractiveIndex5,nEntries);
  quartz->SetMaterialPropertiesTable(myMPT5);
//
// ------------- Volumes --------------


//==============================================
// Volume Overhall
//==============================================
// Geometry has been overhalled to match discussed geometry
//
//===============================================
// The experimental Hall
//
  G4Box* expHall_box = new G4Box("World",fExpHall_x,fExpHall_y,fExpHall_z);

  G4LogicalVolume* expHall_log
    = new G4LogicalVolume(expHall_box,air,"World",0,0,0);

  //G4VPhysicalVolume* expHall_phys
  expHall_phys
  = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);
    expHall_log->SetVisAttributes(G4VisAttributes::Invisible);
    



// The Scintillator
//
    G4Box* Scint_box = new G4Box("Scint",50*mm,0.5*mm,0.5*mm);

  G4LogicalVolume* Scint_log
    = new G4LogicalVolume(Scint_box,Ej228,"Scint",0,0,0);
//G4VPhysicalVolume* Scint_phys =
      fScint = new G4PVPlacement(0,G4ThreeVector(0,0,0*mm),Scint_log,"Scint",
                        expHall_log,false,0);
//SiPM
 G4Box* SiPM_box = new G4Box("SiPM", 0.025*mm,0.5*mm, 0.5*mm);
 G4LogicalVolume* SiPM_log = new G4LogicalVolume(SiPM_box, SiPM_mat, "SiPM",0,0,0);
 fSiPM = new G4PVPlacement(0,G4ThreeVector(-50.025*mm,0,0), SiPM_log,"SiPM", expHall_log,false,0);//make 2 copies
  fSiPM = new G4PVPlacement(0,G4ThreeVector(50.025*mm,0,0), SiPM_log,"SiPM", expHall_log,false,1);//make 2 copies
 

// ------------- Surfaces --------------
//
//Wrapping
//
 
  G4double pp[] = {2.0*eV, 3.5*eV};
  const G4int num2 = sizeof(pp)/sizeof(G4double);
  G4double reflectivity3[] = {1.,1.};
  assert(sizeof(reflectivity2) == sizeof(pp));
  G4double efficiency2[] = {0.0,0.0};
  assert(sizeof(effieciency2) == sizeof(pp));
//SiPM
//
  G4OpticalSurface* opSiPMSurface = new G4OpticalSurface("SiPMSurface");
  opSiPMSurface->SetType(dielectric_dielectric);
  opSiPMSurface->SetFinish(polished);
  opSiPMSurface->SetModel(unified);

  // new G4LogicalBorderSurface("SiPMSurface", fSiPM,Guide_phys,opSiPMSurface);
// Scintillator
//
  G4OpticalSurface* opScintSurface = new G4OpticalSurface("ScintSurface");
  opScintSurface->SetType(dielectric_dielectric);
  opScintSurface->SetFinish(polished);
  opScintSurface->SetModel(unified);

  //G4LogicalSkinSurface* airSurface =
  //      new G4LogicalSkinSurface("ScintSurface", Scint_log, opScintSurface);

  //G4OpticalSurface* opticalSurface = dynamic_cast <G4OpticalSurface*>
  //    (airSurface->GetSurface(Scint_log)->GetSurfaceProperty());

  //if (opticalSurface) opticalSurface->DumpInfo();

//
// Generate & Add Material Properties Table attached to the optical surfaces
//
  const G4int num = 2;
  G4double ephoton[num] = {2.034*eV, 4.136*eV};

  //OpticalGuideSurface
  G4double refractiveIndex[num] = {1.48, 1.48};
  G4double specularLobe[num]    = {0.3, 0.3};
  G4double specularSpike[num]   = {0.2, 0.2};
  G4double backScatter[num]     = {0.2, 0.2};

  G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();

  myST1->AddProperty("RINDEX",                ephoton, refractiveIndex, num);
  myST1->AddProperty("SPECULARLOBECONSTANT",  ephoton, specularLobe,    num);
  myST1->AddProperty("SPECULARSPIKECONSTANT", ephoton, specularSpike,   num);
  myST1->AddProperty("BACKSCATTERCONSTANT",   ephoton, backScatter,     num);

  G4cout << "Water Surface G4MaterialPropertiesTable" << G4endl;
  myST1->DumpTable();

 
  opSiPMSurface->SetMaterialPropertiesTable(myST1);

  //OpticalScintSurface
  G4double reflectivity[num] = {0.3, 0.5};
  G4double efficiency[num]   = {0.8, 1.0};
  G4double reflectivity2[num] = {1.0, 1.0};

  G4MaterialPropertiesTable *myST2 = new G4MaterialPropertiesTable();

  myST2->AddProperty("REFLECTIVITY", ephoton, reflectivity, num);
  myST2->AddProperty("EFFICIENCY",   ephoton, efficiency,   num);

  G4cout << "Air Surface G4MaterialPropertiesTable" << G4endl;
  myST2->DumpTable();

  opScintSurface->SetMaterialPropertiesTable(myST2);
  G4MaterialPropertiesTable *myST3 = new G4MaterialPropertiesTable();

  myST3->AddProperty("REFLECTIVITY", ephoton, reflectivity2, num);
  myST3->AddProperty("EFFICIENCY",   ephoton, efficiency,   num);

  G4cout << "Air Surface G4MaterialPropertiesTable" << G4endl;
  myST3->DumpTable();

  G4MaterialPropertiesTable *myST4 = new G4MaterialPropertiesTable();

  myST4->AddProperty("REFLECTIVITY", ephoton, reflectivity2, num);
  myST4->AddProperty("EFFICIENCY",   ephoton, efficiency,   num);

  G4cout << "Air Surface G4MaterialPropertiesTable" << G4endl;
  myST4->DumpTable();

 



//always return the physical World
  return expHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
