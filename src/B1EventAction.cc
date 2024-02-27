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
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "B1EventAction.hh"
#include "B1RunAction.hh"
#include "B1CalorimeterSD.hh"
#include "B1CalorHit.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction(B1RunAction* runAction,RootManager *rootMng)
: G4UserEventAction(),
  fRunAction(runAction),
    fEdep(0.),
    fsipmHCID (-1)
{
    fRootMng = rootMng;
    fRootMng->initialize();
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
*/
void B1EventAction::BeginOfEventAction(const G4Event*)
{   ////about clock
//   start = clock();
  fStart = 0;
  fEdep = fWLSEdep = 0.;
  fFiberOP = fScintOP = 0;
  fParticleID=0;
  fParticleTime.clear();
  fAbEdep = 0 ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::EndOfEventAction(const G4Event* evt)
{ /*  
    */
    // accumulate statistics in run action
    G4int evtNb = evt->GetEventID();
    G4cout << " ================= Event : "<< evtNb << G4endl;

    // fRootMng->FillSim(fEdep, fWLSEdep, evtNb, fFiberOP, fScintOP, fParticleID, fParticleTime, fAbEdep);
    fRootMng->Fill();

    // //about clock
    // end = clock();
    // dur = (double)(end-start);
    // std::cout<< "Event Use Time: " << (dur/CLOCKS_PER_SEC) << "s" << "\n" ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
