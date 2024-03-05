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
/// \file B1CalorimeterSD.hh
/// \brief Definition of the B1CalorimeterSD class

#ifndef B1CalorimeterSD_h
#define B1CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"

#include "B1CalorHit.hh"
//now we use rootmgr or IO, skip calohit vector
#include "RootManager.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

/// Calorimeter sensitive detector class
///
/// In Initialize(), it creates one hit for each calorimeter layer and one more
/// hit for accounting the total quantities in all layers.
///
/// The values are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step.

class B1CalorimeterSD : public G4VSensitiveDetector
{
  public:
    B1CalorimeterSD(const G4String& name, RootManager *rootMng, G4int nofLayers, G4int nofCells);
    virtual ~B1CalorimeterSD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory*ROhist);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);
    // void clean();

  private:
    //B1EventAction* fEventAction;
    B1CalorHitsCollection* fHitsCollection;
    G4int fNofLayers;
    G4int fNofCells;
    RootManager * fRootMgr;

  
    // int eID;
    // double eEnergy;
    // double eTime;
    int nPhoton;
    std::vector<int> sipm_photons;
    int sipm_photon[5] = {0};
    int sipm_photon_num = 0;
    int photon_layer_cell[200][200] = {0};
    // int Layer_num;
    int Cell_num;
    std::vector<int> photon_num;
    // int layer_id;
    // int bar_id;
    // int xy_id;
    // int layer_num = 83 ;
    // int  perbar_pho_x[300][25]={0};
    // int  perbar_pho_y[300][25]={0};
    // std::vector<int> pho_bar_x;
    // std::vector<int> pho_bar_y;
    // std::vector<int> photon_num = std::vector<int>(layer_num,0);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

