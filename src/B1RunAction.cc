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
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
#include "B1CalorimeterSD.hh"
// #include "B1Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction(RootManager *rootMng)
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.)
{ 
    fRootMng = rootMng;  
    fRootMng->book();
  //// add new units for dose
   
  //const G4double milligray = 1.e-3*gray;
  //const G4double microgray = 1.e-6*gray;
  //const G4double nanogray  = 1.e-9*gray;  
  //const G4double picogray  = 1.e-12*gray;
   
  //new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  //new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  //new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  //new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray); 

  //// Register accumulable to the accumulable manager
  //G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  //accumulableManager->RegisterAccumulable(fEdep);
  //accumulableManager->RegisterAccumulable(fEdep2); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run* run)
{ ////about clock
  // start = clock();
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  //// reset accumulables to their initial values

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* run)
{

  G4int nofEvents = run->GetNumberOfEvent();
    G4cout << " ================= Number of  Event : "<< nofEvents << G4endl;
  if (nofEvents == 0) return;

  //  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }

  fRootMng->saveTree();

  ////About clock
  // end = clock();
  // dur = (double)(end-start);
  // std::cout<<"Run using time: "  << (dur/CLOCKS_PER_SEC) << "s" <<"\n" ;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::AddEdep(G4double edep)
{
    fEdep  += edep;
}
void B1RunAction::AddSipmEdep(G4double edep)
{
    SipmEdep += edep;
}

void B1RunAction::AddWLSEdep(G4double WLSEdep)
{
    fWLSEdep  += WLSEdep;
}

void B1RunAction::AddNumber(G4int num)
{
    fNoPhotons += num;
}
void B1RunAction::AddSipmNumber(G4int num)
{
    SipmPhotons += num;
}

void B1RunAction::GetPID(G4int pid)
{
    fParticleID = pid;
}
void B1RunAction::GetEventID(G4int id)
{
    SipmEid = id;
}


void B1RunAction::GetTime(G4double Partime)
{
    fParticleTime.push_back(Partime);
}
void B1RunAction::GetSipmTime(G4double Partime)
{
    SipmEtime.push_back(Partime);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

