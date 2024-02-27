//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class
#include <vector>
#include "B1DetectorConstruction.hh"
#include "B1CalorimeterSD.hh"

#include "G4SDParticleFilter.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"

#include "G4MaterialTable.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4VPhysicalVolume.hh"

using std::vector;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction(RootManager *rootMng)
    : G4VUserDetectorConstruction(),
      fScoringVolume(0)
{
    fRootMng = rootMng;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *B1DetectorConstruction::Construct()
{
    // Get nist material manager
    G4NistManager *nist = G4NistManager::Instance();
    DefineMaterials();

    // Option to switch on/off checking of volumes overlaps
    G4bool checkOverlaps = true;
    
    // World
    G4double world_sizeX = 200 * cm;
    G4double world_sizeY = 200 * cm;
    G4double world_sizeZ = 400 * cm;
    G4Material *world_mat = nist->FindOrBuildMaterial("G4_AIR");

    G4Box *solidWorld =
        new G4Box("World",                                                  // its name
                  0.5 * world_sizeX, 0.5 * world_sizeY, 0.5 * world_sizeZ); // its size

    G4LogicalVolume *logicWorld =
        new G4LogicalVolume(solidWorld, // its solid
                            world_mat,  // its material
                            "World");   // its name

    G4VPhysicalVolume *physWorld =
        new G4PVPlacement(0,               // no rotation
                          G4ThreeVector(), // at (0,0,0)
                          logicWorld,      // its logical volume
                          "World",         // its name
                          0,               // its mother  volume
                          false,           // no boolean operation
                          0,               // copy number
                          checkOverlaps);  // overlaps checking

    // Geometry Parameter
    G4double Scint_X = 5 * cm, Scint_Y = 75 * cm, Scint_Z = 1 * cm; // Scintillator geometry parameter

    G4double groove_R = 0.5 * mm, groove_Y = Scint_Y; // Fiber geometry parameter
    G4double fiber_extrude = 0. * mm; // this is to simulate the fiber is lightly longer than the scintalor, thus extrude a little bit.

    G4double clad_R = groove_R, clad_Y = (Scint_Y + fiber_extrude * 2) / 2; // Fiber clad geometry parameter
    G4double core_R = clad_R * 0.98, core_Y = clad_Y;
    
    G4double ESR_thickness = 0.08 * mm; // ESR geometry parameter
    G4double ESR_X = Scint_X + 2 * ESR_thickness, ESR_Y = Scint_Y + 2 * ESR_thickness, ESR_Z = Scint_Z + 2 * ESR_thickness;

    G4double SD_X = 3. * mm, SD_Y = 1. * mm, SD_Z = 3. * mm; // SiPM geometry parameter

    G4double Cell_X = ESR_X, Cell_Y = ESR_Y, Cell_Z = ESR_Z; // Cell geometry parameter

    G4double Absorber_X = Cell_Y, Absorber_Y = Cell_Y, Absorber_Z = 1 * cm; // Absorber geometry parameter

    auto* alongY = new G4RotationMatrix();
    alongY->rotateX(-90 * deg);

    auto* alongX = new G4RotationMatrix();
    alongX->rotateZ(-90 * deg);

    // Element
    auto* fPMMA = new G4Material("PMMA", density = 1190 * kg / m3, 3);
    auto* fH = new G4Element("H", "H", z = 1., a = 1.01 * g / mole);
    auto* fC = new G4Element("C", "C", z = 6., a = 12.01 * g / mole);
    auto* fN = new G4Element("N", "N", z = 7., a = 14.01 * g / mole);
    auto* fO = new G4Element("O", "O", z = 8., a = 16.00 * g / mole);

    // Cell
    G4Material* Cell_mat = nist->FindOrBuildMaterial("G4_AIR"); // air
    G4Box* solidCell =
        new G4Box("Cell",  // its name
            0.5 * Cell_X, 0.5 * Cell_Y, 0.5 * Cell_Z); // its size

    G4LogicalVolume* logicCell =
        new G4LogicalVolume(solidCell, // its solid
            Cell_mat,  // its material
            "Cell");   // its name

    // Absorber
    G4Material* Absorber_mat = nist->FindOrBuildMaterial("G4_Fe");
    G4Box* solidAbsorber =
        new G4Box("Absorber", 0.5 * Absorber_X, 0.5 * Absorber_Y, 0.5 * Absorber_Z);

    G4LogicalVolume* logicAbsorber =
        new G4LogicalVolume(solidAbsorber, // its solid
                            Absorber_mat,  // its material
                            "Absorber");   // its name

    // ESR surface
    ESR_mat = new G4Material("ESR", density = 1.38 * g / cm3, ncomponents = 3);
    ESR_mat->AddElement(C, 10);
    ESR_mat->AddElement(H, 8);
    ESR_mat->AddElement(O, 4);

    G4Box* solidESR_out =
        new G4Box("ESR", 0.5 * ESR_X, 0.5 * ESR_Y, 0.5 * ESR_Z);

    G4Box* solidESR_inner =
        new G4Box("ESR", 0.5 * ESR_X - ESR_thickness, 0.5 * ESR_Y - ESR_thickness, 0.5 * ESR_Z - ESR_thickness);

    G4SubtractionSolid* solidESR =
        new G4SubtractionSolid("ESR_Without_Groove", solidESR_out, solidESR_inner, alongY, G4ThreeVector(0, 0, 0));

    G4double ESR_groove_Y = ESR_Y;

    G4Tubs* solidESR_Groove =
        new G4Tubs("Groove", 0 * mm, groove_R, ESR_groove_Y, 0 * deg, 360 * deg);

    G4SubtractionSolid* solidESR_sub_Groove1 =
        new G4SubtractionSolid("ESRWithGroove", solidESR, solidESR_Groove, alongY, G4ThreeVector(0, 0, (Scint_Z / 2 - groove_R)));

    G4SubtractionSolid* solidESR_sub_Groove2 =
        new G4SubtractionSolid("ESRWithGroove", solidESR_sub_Groove1, solidESR_Groove, alongY, G4ThreeVector((Scint_X / 3), 0, ((Scint_Z / 2 - groove_R) * (-1))));

    G4SubtractionSolid* solidESR_sub_Groove =
        new G4SubtractionSolid("ESRWithGroove", solidESR_sub_Groove2, solidESR_Groove, alongY, G4ThreeVector((Scint_X / (-3)), 0, ((Scint_Z / 2 - groove_R) * (-1))));

    G4LogicalVolume* logicESR =
        new G4LogicalVolume(solidESR_sub_Groove, // its solid
                            ESR_mat,             // its material
                            "ESRWithGroove_LV");

    new G4PVPlacement(0,               // no rotation
                      G4ThreeVector(), // no translation
                      logicESR,      // its logical volume;
                      "ESR",         // its name
                      logicCell,        // its mother volume
                      false,           // no boolean operation
                      0,               // copy number
                      checkOverlaps);  // overlaps checking

    const G4int cNum = 2;
    G4double ephoton[cNum] = {1 * eV, 7 * eV}; // overflow and underflow will take first/last bin
    G4double reflectivity[cNum] = {0.95, 0.95};  // VIP
    G4double efficiency[cNum] = {0.0, 0.0};
    G4double transmittance[cNum] = {0.0, 0.0};
    auto Wrap_Surface_Mat = new G4MaterialPropertiesTable();
    Wrap_Surface_Mat->AddProperty("REFLECTIVITY", ephoton, reflectivity, cNum); // reflect fraction, default=1
    Wrap_Surface_Mat->AddProperty("EFFICIENCY", ephoton, efficiency,
                                  cNum);                                          // detection  fraction (abs=1-reflet-trans, then at efficiency, detectition(invoke post-step SD)),default=0
    Wrap_Surface_Mat->AddProperty("TRANSMITTANCE", ephoton, transmittance, cNum); // transmission fraction, default=0
    auto Wrap_Surface = new G4OpticalSurface("WrapSurfaceOptical");
    Wrap_Surface->SetType(dielectric_LUT);
    Wrap_Surface->SetModel(LUT);
    Wrap_Surface->SetFinish(polishedtyvekair);
    Wrap_Surface->SetMaterialPropertiesTable(Wrap_Surface_Mat);

    new G4LogicalSkinSurface("ESR_surface", logicESR, Wrap_Surface); // here just use the inner skin surface, which is just between ESR and the scintallator

    // scintalator sub groove
    // ############ material start
    //  material defination #################################################################################
    //  Scintillation Box, PS, 10000/MeV, enable Birks
    G4Material *scintillatorMat = nist->FindOrBuildMaterial("G4_POLYSTYRENE"); // PS
    G4double wls_Energy[] = {2.00 * eV, 2.87 * eV, 2.90 * eV, 3.47 * eV};
    const G4int wlsnum = sizeof(wls_Energy) / sizeof(G4double);

    G4double rIndexPstyrene[] = {1.58, 1.58, 1.58, 1.58};
    assert(sizeof(rIndexPstyrene) == sizeof(wls_Energy));
    G4double absorption1[] = {50. * cm, 50. * cm, 50. * cm, 50. * cm};
    assert(sizeof(absorption1) == sizeof(wls_Energy));
    G4double scintilFast[] = {0.00, 0.00, 1.00, 1.00};
    assert(sizeof(scintilFast) == sizeof(wls_Energy));
    auto *MPTPStyrene = new G4MaterialPropertiesTable();
    MPTPStyrene->AddProperty("RINDEX", wls_Energy, rIndexPstyrene, wlsnum);
    MPTPStyrene->AddProperty("ABSLENGTH", wls_Energy, absorption1, wlsnum);
    MPTPStyrene->AddProperty("FASTCOMPONENT", wls_Energy, scintilFast, wlsnum);
    MPTPStyrene->AddConstProperty("SCINTILLATIONYIELD", 7700. / MeV); // control of photon yield
    MPTPStyrene->AddConstProperty("RESOLUTIONSCALE", 1.0);
    MPTPStyrene->AddConstProperty("FASTTIMECONSTANT", 2.4 * ns);
    scintillatorMat->SetMaterialPropertiesTable(MPTPStyrene);
    // Set the Birks Constant for the Polystyrene scintillator
    scintillatorMat->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
    // ############ material definition end
    
    G4Box *solidScintilator =
        new G4Box("Scintillator_NoGroove",
                  0.5 * Scint_X, 0.5 * Scint_Y, 0.5 * Scint_Z);
    
    //打孔
    G4Tubs *solidGroove =
        new G4Tubs("Groove", 0 * mm, groove_R, groove_Y, 0 * deg, 360 * deg);

    G4SubtractionSolid *solidScint_sub_Groove1 =
        new G4SubtractionSolid("ScintWithGroove", solidScintilator, solidGroove, alongY, G4ThreeVector(0, 0, (Scint_Z/2-groove_R))); 

    G4SubtractionSolid *solidScint_sub_Groove2 =
        new G4SubtractionSolid("ScintWithGroove", solidScint_sub_Groove1, solidGroove, alongY, G4ThreeVector((Scint_X/3), 0, ((Scint_Z/2-groove_R)*(-1))));

    G4SubtractionSolid *solidScint_sub_Groove =
        new G4SubtractionSolid("ScintWithGroove", solidScint_sub_Groove2, solidGroove, alongY, G4ThreeVector((Scint_X/(-3)), 0, ((Scint_Z/2-groove_R)*(-1))));

    G4LogicalVolume *logicScint =
        new G4LogicalVolume(solidScint_sub_Groove, // its solid
                            scintillatorMat,       // its material
                            "ScintWithGroove_LV");

    new G4PVPlacement(0,               // no rotation
                      G4ThreeVector(), // no translation
                      logicScint,      // its logical volume;
                      "Scint",         // its name
                      logicCell,        // its mother volume
                      false,           // no boolean operation
                      0,               // copy number
                      checkOverlaps);  // overlaps checking

    // clad
    // ############## materiali start, PE, n1.49
    G4Material *cladMat = nist->FindOrBuildMaterial("G4_POLYETHYLENE"); // PE

    G4double RefractiveIndexClad1[] = {1.49, 1.49, 1.49, 1.49};
    G4double AbsFiber[] = {3.5 * m, 3.5 * m, 0.1 * mm, 0.1 * mm};
    assert(sizeof(RefractiveIndexClad1) == sizeof(wls_Energy));
    G4MaterialPropertiesTable *clad1Property = new G4MaterialPropertiesTable();
    clad1Property->AddProperty("RINDEX", wls_Energy, RefractiveIndexClad1, wlsnum);
    clad1Property->AddProperty("ABSLENGTH", wls_Energy, AbsFiber, wlsnum);
    cladMat->SetMaterialPropertiesTable(clad1Property);
    // ############# material ends

    G4Tubs *solidClad =
        new G4Tubs("Clad_Solid", core_R, clad_R, clad_Y, 0 * deg, 360 * deg);

    G4LogicalVolume *logicClad =
        new G4LogicalVolume(solidClad,  // its solid
                            cladMat,    // its material
                            "Clad_LV"); // its name
    new G4PVPlacement(alongY,           // no rotation
                      G4ThreeVector(0, 0, (Scint_Z/2-groove_R)),  // no translation
                      logicClad,        // its logical volume;
                      "Clad",           // its name
                      logicCell,         // its mother volume
                      false,            // no boolean operation
                      0,                // copy number
                      checkOverlaps);   // overlaps checking

    new G4PVPlacement(alongY,           // no rotation
                      G4ThreeVector((Scint_X/3), 0, ((Scint_Z/2-groove_R)*(-1))),  // no translation
                      logicClad,        // its logical volume;
                      "Clad",           // its name
                      logicCell,         // its mother volume
                      false,            // no boolean operation
                      1,                // copy number
                      checkOverlaps);   // overlaps checking

    new G4PVPlacement(alongY,           // no rotation
                      G4ThreeVector((Scint_X/(-3)), 0, ((Scint_Z/2-groove_R)*(-1))),  // no translation
                      logicClad,        // its logical volume;
                      "Clad",           // its name
                      logicCell,         // its mother volume
                      false,            // no boolean operation
                      2,                // copy number
                      checkOverlaps);   // overlaps checking

    // wls fiber core
    // ############## material start, PMMA, n1.60
    G4double density, a, z; // atmoic mass ,atomic number
    G4int polyPMMA = 1;
    G4int nC_PMMA = 3 + 2 * polyPMMA;
    G4int nH_PMMA = 6 + 2 * polyPMMA;
    fPMMA->AddElement(fH, 8);
    fPMMA->AddElement(fC, 5);
    fPMMA->AddElement(fO, 2);
    G4Material *FiberMat = new G4Material("PMMA", 1190 * kg / m3, 3); // PMMA
    FiberMat->AddElement(fH, nH_PMMA);
    FiberMat->AddElement(fC, nC_PMMA);
    FiberMat->AddElement(fO, 2);

    G4double RefractiveIndexFiber[] = {1.60, 1.60, 1.60, 1.60};
    assert(sizeof(RefractiveIndexFiber) == sizeof(wls_Energy));
    // G4double AbsFiber[]={9.00*m,9.00*m,0.1*mm,0.1*mm};
    assert(sizeof(AbsFiber) == sizeof(wls_Energy));
    G4double EmissionFib[] = {1.0, 1.0, 0.0, 0.0};
    assert(sizeof(EmissionFib) == sizeof(wls_Energy));
    G4MaterialPropertiesTable *fiberProperty = new G4MaterialPropertiesTable();
    fiberProperty->AddProperty("RINDEX", wls_Energy, RefractiveIndexFiber, wlsnum);
    fiberProperty->AddProperty("WLSABSLENGTH", wls_Energy, AbsFiber, wlsnum);
    fiberProperty->AddProperty("WLSCOMPONENT", wls_Energy, EmissionFib, wlsnum);
    fiberProperty->AddConstProperty("WLSTIMECONSTANT", 8.8 * ns);
    FiberMat->SetMaterialPropertiesTable(fiberProperty);
    // ############# material ends

    G4Tubs *solidFiber =
        new G4Tubs("Fiber_Solid", 0 * mm, core_R, core_Y, 0 * deg, 360 * deg);

    G4LogicalVolume *logicFiber =
        new G4LogicalVolume(solidFiber,               // its solid
                            FiberMat,                 // its material
                            "Fiber_LV");              // its name
    
    auto *wls_PV = new G4PVPlacement(alongY,          // no rotation
                                     G4ThreeVector(0, 0, (Scint_Z/2-groove_R)), // no translation
                                     logicFiber,      // its logical volume;
                                     "Fiber",         // its name
                                     logicCell,       // its mother volume
                                     false,           // no boolean operation
                                     0,               // copy number
                                     checkOverlaps);  // overlaps checking

    auto *wls_PV1 = new G4PVPlacement(alongY,          // no rotation
                                      G4ThreeVector((Scint_X/3), 0, ((Scint_Z/2-groove_R)*(-1))), // no translation
                                      logicFiber,      // its logical volume;
                                      "Fiber",         // its name
                                      logicCell,       // its mother volume
                                      false,           // no boolean operation
                                      1,               // copy number
                                      checkOverlaps);  // overlaps checking


    auto *wls_PV2 = new G4PVPlacement(alongY,          // no rotation
                                      G4ThreeVector((Scint_X/(-3)), 0, ((Scint_Z/2-groove_R)*(-1))), // no translation
                                      logicFiber,      // its logical volume;
                                      "Fiber",         // its name
                                      logicCell,       // its mother volume
                                      false,           // no boolean operation
                                      2,               // copy number
                                      checkOverlaps);  // overlaps checking

    // sipm
    G4Material *SD_mat = nist->FindOrBuildMaterial("G4_Si");

    G4Box *solidSIPM =
        new G4Box("SIPM_Solid",
                  0.5 * SD_X, 0.5 * SD_Y, 0.5 * SD_Z);

    G4LogicalVolume *logicSIPM =
        new G4LogicalVolume(solidSIPM,  // its solid
                            SD_mat,     // its material
                            "SIPM_LV"); // its name
                            
    auto *sipm_PV = new G4PVPlacement(0,
                                      G4ThreeVector(0, 0.5 * Scint_Y + fiber_extrude + 0.5 * SD_Y, (Scint_Z/2-groove_R)), // ensure by yourself not exceed the logicBar
                                      logicSIPM,
                                      "SIPM",
                                      logicCell,
                                      false,
                                      0,
                                      checkOverlaps);


    auto *sipm_PV1 = new G4PVPlacement(0,
                                      G4ThreeVector((Scint_X/3), 0.5 * Scint_Y + fiber_extrude + 0.5 * SD_Y, ((Scint_Z/2-groove_R)*(-1))), // ensure by yourself not exceed the logicBar
                                      logicSIPM,
                                      "SIPM",
                                      logicCell,
                                      false,
                                      1,
                                      checkOverlaps); 
    
    auto *sipm_PV2 = new G4PVPlacement(0,
                                      G4ThreeVector((Scint_X/(-3)), 0.5 * Scint_Y + fiber_extrude + 0.5 * SD_Y, ((Scint_Z/2-groove_R)*(-1))), // ensure by yourself not exceed the logicBar
                                      logicSIPM,
                                      "SIPM",
                                      logicCell,
                                      false,
                                      2,
                                      checkOverlaps);

    // SiPM surface
    G4double reflectivityAPD[cNum] = {0.0, 0.0};
    G4double efficiencyAPD[cNum] = {0.2, 0.2}; // VIP
    G4double transmittanceAPD[cNum] = {0.0, 0.0};
    auto APD_Surface_Mat = new G4MaterialPropertiesTable();
    APD_Surface_Mat->AddProperty("REFLECTIVITY", ephoton, reflectivityAPD, cNum); // reflect fraction, default=1
    APD_Surface_Mat->AddProperty("EFFICIENCY", ephoton, efficiencyAPD,
                                 cNum);                                             // detection  fraction (abs=1-reflet-trans, then at efficiency, detectition(invoke post-step SD)),default=0
    APD_Surface_Mat->AddProperty("TRANSMITTANCE", ephoton, transmittanceAPD, cNum); // transmission fraction, default=0

    auto APD_Surface = new G4OpticalSurface("APDSurfaceOptical");
    APD_Surface->SetType(dielectric_LUT);
    APD_Surface->SetModel(LUT);
    APD_Surface->SetFinish(polishedvm2000glue);
    APD_Surface->SetMaterialPropertiesTable(APD_Surface_Mat);

    new G4LogicalBorderSurface("SIPM_surface", wls_PV, sipm_PV, APD_Surface);
    new G4LogicalBorderSurface("SIPM_surface", wls_PV1, sipm_PV1, APD_Surface);
    new G4LogicalBorderSurface("SIPM_surface", wls_PV2, sipm_PV2, APD_Surface);

    // Put cell
    G4int nLayer = 10, nCell = 15;
    
    // each layer of the prototype   
    for (G4int i = 1; i < nLayer; i++)
    {
        new G4PVPlacement(0, // no rotation
                          G4ThreeVector(0, 0, i * (Absorber_Z + Cell_Z)),
                          logicAbsorber,          // its logical volume;
                          "absorber",        // its name
                          logicWorld,        // its mother volume
                          false,             // no boolean operation
                          i, // copy number
                          checkOverlaps);    // overlaps checking

        if ((i%2) == 0)
            for (G4int j = 0; j < nCell; j++)
            {
                new G4PVPlacement(0, // no rotation
                                  G4ThreeVector((j - 7) * Cell_X, 0, (i + 0.5) * (Absorber_Z + Cell_Z)),
                                  logicCell,          // its logical volume;
                                  "scin_bar",        // its name
                                  logicWorld,        // its mother volume
                                  false,             // no boolean operation
                                  i * nCell + j, // copy number
                                  checkOverlaps);    // overlaps checking
            }
        else
            for (G4int j = 0; j < nCell; j++)
            {
                new G4PVPlacement(alongX, // no rotation
                                  G4ThreeVector(0, (j - 7) * Cell_X, (i + 0.5) * (Absorber_Z + Cell_Z)),
                                  logicCell,          // its logical volume;
                                  "scin_bar",        // its name
                                  logicWorld,        // its mother volume
                                  false,             // no boolean operation
                                  i * nCell + j, // copy number
                                  checkOverlaps);    // overlaps checking
            }
    }
    
                      
    /*loop part for building stand alnoe HCAL

    G4double looplength_x = 5.2 * cm, looplength_z = 4.71 * cm;
    G4double box_hx = 100 * cm, box_hy = 100 * cm, box_hz = 2 * cm; // size of absorber ,absorbersize

    G4double absorber_x, absorber_y, absorber_z = 0.5 * box_hz + 0.5 * env_sizeZ;
    G4Box *absorber_box =
        new G4Box("absorber_box", 0.5 * box_hx, 0.5 * box_hy, 0.5 * box_hz);
    G4LogicalVolume *absorber =
        new G4LogicalVolume(absorber_box,
                            absorber_mat,
                            "absorber_Fe");

    int copy_num = 0, looptime_scin = 20, looptime_double = 83; // looptime_double is the loop time of

    G4Box *solidScin_layer =
        new G4Box("Scintillator_layer",                                                      // its name
                  0.5 * env_sizeX * looptime_scin, 0.5 * env_sizeY + SD_Y, 0.5 * env_sizeZ); // its size
    G4LogicalVolume *logiclayer =                                                            // logicvolume of scintillator
        new G4LogicalVolume(solidScin_layer,                                                 // its solid
                            world_mat,                                                       // its material
                            "layer_envlope");                                                // its name
    int copynum_scinlayer = 0;
    G4double bar_x = -49.4 * cm;
    for (int i = 0; i < looptime_scin; i++)
    {
        new G4PVPlacement(0, // no rotation
                          G4ThreeVector(bar_x, 0, 0),
                          logicBar,          // its logical volume;
                          "scin_bar",        // its name
                          logiclayer,        // its mother volume
                          false,             // no boolean operation
                          copynum_scinlayer, // copy number
                          checkOverlaps);    // overlaps checking
        bar_x += env_sizeX;
        copynum_scinlayer++;
    }

    G4Box *solid_envelope =                                                                           // solid of the envelope include scintillator and the absorber(box)
        new G4Box("scin_env",                                                                         // its name
                  0.5 * env_sizeX * looptime_scin, 0.5 * env_sizeY + SD_Y, env_sizeZ + 0.5 * box_hz); // its size
    G4LogicalVolume *logic_envelope =
        new G4LogicalVolume(solid_envelope, // its solid
                            world_mat,      // its material
                            "scin_env");    // its name

    // here we create a whole hcal layer
    G4RotationMatrix *rotm = new G4RotationMatrix();
    rotm->rotateZ(90 * deg);
    int copynum_envelope = 0;
    new G4PVPlacement(rotm, // rotation 90 degrees
                      G4ThreeVector(0, 0, -0.5 * box_hz),
                      logiclayer,       // its logical volume
                      "scin_layer1",    // its name
                      logic_envelope,   // its mother volume
                      false,            // no boolean operation
                      copynum_envelope, // copy number
                      checkOverlaps);
    copynum_envelope++;
    new G4PVPlacement(0, // rotation
                      G4ThreeVector(0, 0, -0.5 * box_hz + env_sizeZ),
                      logiclayer,       // its logical volume
                      "scin_layer2",    // its name
                      logic_envelope,   // its mother volume
                      false,            // no boolean operation
                      copynum_envelope, // copy number
                      checkOverlaps);
    copynum_envelope++;
    new G4PVPlacement(0, // rotation
                      G4ThreeVector(0, 0, 1.5 * env_sizeZ),
                      absorber,         // its logical volume
                      "Absorber_Fe",    // its name
                      logic_envelope,   // its mother volume
                      false,            // no boolean operation
                      copynum_envelope, // copy number
                      checkOverlaps);

    G4double logenv_z = 0 * cm;
    // int copynum_whole = 0;
    for (int i = 0; i < looptime_double; i++)
    {
        new G4PVPlacement(0, // rotation
                          G4ThreeVector(0, 0, logenv_z),
                          logic_envelope,     // its logical volume
                          "scin_layer_whole", // its name
                          logicWorld,         // its mother volume
                          false,              // no boolean operation
                          i,                  // copy number
                          checkOverlaps);

        logenv_z += (box_hz + 2 * env_sizeZ);
    }
    */
    //  G4Colour ESR_colour(1.,1.,0.);
    //     logicESR->SetVisAttributes(G4VisAttributes(ESR_colour));
    G4Colour ESR_colour(1., 1., 0.); // yellow
    // logiclayer->SetVisAttributes(G4VisAttributes(ESR_colour));

    G4Colour Steel_colour(0.45, 0.25, 0.0); // brown
    // absorber->SetVisAttributes(G4VisAttributes(Steel_colour));

    // logic_envelope->SetVisAttributes(G4VisAttributes(false));

    // for (int i=0 ; i<looptime_scin*looptime_double ; i++){
    // G4LogicalBorderSurface* WLSSurface1 = new G4LogicalBorderSurface("WLS_surface1",fiber_PV1,scint_PV,Wrap_Surface);
    //     G4LogicalBorderSurface* WLSSurface1 = new G4LogicalBorderSurface("WLS_surface1",fiber_PV1vector.at(i),scint_PVvector.at(i),Wrap_Surface);
    //     WLSSurface1_vector.push_back(WLSSurface1);
    // }
    //     for (int i=0 ; i<looptime_scin*looptime_double ; i++)
    //     G4LogicalBorderSurface* WLSSurface1 = new G4LogicalBorderSurface("WLS_surface1",fiber_PV1vector_x.at(i),scint_PVvector_x.at(i),Wrap_Surface);
    //     WLSSurface1_vector_x.push_back(WLSSurface1);
    // }
    G4Colour core_colour(0., 0., 1.); // Blue
    logicFiber->SetVisAttributes(G4VisAttributes(core_colour));
    G4Colour SD_colour(1., 4., 2);
    logicSIPM->SetVisAttributes(G4VisAttributes(SD_colour));

    // set sipm to sd
    auto sipmSD = new B1CalorimeterSD("/SiliconPMSD", fRootMng);
    G4SDManager::GetSDMpointer()->AddNewDetector(sipmSD);
    logicSIPM->SetSensitiveDetector(sipmSD);
    // Set scintillator to sensitive detector
    // auto scinSD = new B1CalorimeterSD_scin("/ScintillatorSD", fRootMng);
    // G4SDManager::GetSDMpointer()->AddNewDetector(scinSD);
    // logicScint->SetSensitiveDetector(scinSD);

    // auto scinSD_x = new B1CalorimeterSD_scin("/ScintillatorSD_x" , fRootMng);
    // G4SDManager::GetSDMpointer()->AddNewDetector(scinSD_x);
    // logicShape->SetSensitiveDetector(scinSD_x);

    // sipm surface
    //  APD related
    //  G4double reflectivityAPD[cNum] = {0.0, 0.0};
    //  G4double efficiencyAPD[cNum] = {1.0, 1.0};
    //  G4double transmittanceAPD[cNum] = {0.0, 0.0};
    //  auto APD_Surface_Mat = new G4MaterialPropertiesTable();
    //  APD_Surface_Mat->AddProperty("REFLECTIVITY", ephoton, reflectivityAPD, cNum); //reflect fraction, default=1
    //  APD_Surface_Mat->AddProperty("EFFICIENCY", ephoton, efficiencyAPD,
    //                               cNum); //detection  fraction (abs=1-reflet-trans, then at efficiency, detectition(invoke post-step SD)),default=0
    //  APD_Surface_Mat->AddProperty("TRANSMITTANCE", ephoton, transmittanceAPD, cNum); //transmission fraction, default=0

    // auto APD_Surface = new G4OpticalSurface("APDSurfaceOptical");
    // APD_Surface->SetType(dielectric_LUT);
    // APD_Surface->SetModel(LUT);
    // APD_Surface->SetFinish(polishedvm2000glue);
    // APD_Surface->SetMaterialPropertiesTable(APD_Surface_Mat);

    // auto *SIPMSurface1 = new G4LogicalSkinSurface("SIPM_surface1",logicSIPM,APD_Surface);
    // auto *SIPMSurface2 = new G4LogicalSkinSurface("SIPM_surface2",logicSD2,APD_Surface);

    //
    // always return the physical World
    //
    return physWorld;
}

void B1DetectorConstruction::DefineMaterials()
{
    /* Define Optical Properties */
    /*
     */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
