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
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1DetectorConstruction.hh"
#include "B1CalorimeterSD.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction,RootManager *rootMng)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
    
    if(step->GetTrack()->GetTrackID()==1 && step->GetTrack()->GetCurrentStepNumber()==1)
    {
       G4double preEnergy = step->GetTrack()->GetKineticEnergy();
    //    G4cout << "preEnergy is "<<preEnergy<< G4endl;
    }
    G4int particleID = step->GetTrack()->GetDefinition()->GetPDGEncoding();
    G4String name = step->GetTrack()->GetDefinition()->GetParticleName();
    G4String processN = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    G4String Volume = step->GetTrack()->GetVolume()->GetName();

    G4int OPinFiber = 0, OPinScint = 0;
    G4double ab_edep = 0;
    if ( Volume == "Scintillator" && name == "opticalphoton" )
    {
        OPinScint = 1;
    }
    if( Volume == "Fiber" && name == "opticalphoton" )
    {
        OPinFiber = 1;
    }
    if (Volume == "Absorber_Fe"){
        ab_edep = step->GetTotalEnergyDeposit();
        // std::cout << "AbsorberEdep: "<<ab_edep <<" MeV" << std::endl; 
        
    }
    G4double WLSEdep;
    if (processN == "OpWLS")
    {
        WLSEdep = step->GetTrack()->GetKineticEnergy();
        fEventAction->AddWLSEdep(WLSEdep);
    }

    G4double ParTime = step->GetPreStepPoint()->GetGlobalTime();
    G4double edep = step->GetTotalEnergyDeposit();
    fEventAction->AddEdep(edep);
    fEventAction->GetPID(particleID);
    fEventAction->GetTime(ParTime);
    fEventAction->AddFiberOP(OPinFiber);
    fEventAction->AddScintOP(OPinScint);
    fEventAction->GetAbsoEdep(ab_edep);
    

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

