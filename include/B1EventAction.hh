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
/// \file B1EventAction.hh
/// \brief Definition of the B1EventAction class

#ifndef B1EventAction_h
#define B1EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "RootManager.hh"
#include "B1CalorimeterSD.hh"
#include "B1ScintHit.hh"
#include "B1ScintillatorSD.hh"
#include <time.h>

class B1RunAction;

/// Event action class
///

class B1EventAction : public G4UserEventAction
{
  public:
    B1EventAction(B1RunAction* runAction, RootManager * rootMng);
    virtual ~B1EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    void AddEdep(G4double edep) { fEdep += edep; }
    void AddWLSEdep(G4double WLSEdep) { fWLSEdep += WLSEdep; }
    void AddFiberOP(G4int num) { fFiberOP += num; }
    void AddScintOP(G4int num) { fScintOP += num; }
    void GetPID(G4int pid)      { fParticleID = pid; }
    void GetTime(G4double time)  { fParticleTime.push_back(time);}

    void AddSipmEdep(G4double edep) { SipmEdep += edep; }
    void AddSipmNumber(G4int num) { SipmPhotons += num; }
    void GetEventID(G4int id) { SipmEid = id; }
    void GetSipmTime(G4double time) { SipmEtime.push_back(time); }

    void GetAbsoEdep(G4double edep){fAbEdep += edep;}
    clock_t start,end;
    G4double dur;

  private:
    B1RunAction* fRunAction;
    G4double     fEdep;
    G4double     fEvtNb;
    G4double     fStart;
    G4double     fWLSEdep;
    G4double     fFiberOP;
    G4double     fScintOP;
    G4int        fParticleID;
    G4double     fTime;
    std::vector<double>     fParticleTime;
    RootManager * fRootMng;

    G4double SipmEdep;
    G4int SipmPhotons;
    G4int SipmEid;
    std::vector<double> SipmEtime;
    G4double     fAbEdep;

    B1ScintHitsCollection* GetHitsCollection(G4int hcID,
                                             const G4Event* event) const;
    void PrintEventStatistics(G4double absoEdep, G4double absoTrackLength) const;

    G4int fsipmHCID;
    G4int fScinHCID = -1;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
