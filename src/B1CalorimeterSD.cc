#include "B1CalorimeterSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "G4Event.hh"


B1CalorimeterSD::B1CalorimeterSD(const G4String& name, const G4String& hitsCollectionName, RootManager *rootMng, G4int nofLayers, G4int nofCells) : G4VSensitiveDetector(name),
fNofLayers(nofLayers),
fNofCells(nofCells),
fRootMgr(rootMng)
{
    collectionName.insert(hitsCollectionName);
    G4cout << "SensitiveDetector Processed Successfully " << G4endl;
}

B1CalorimeterSD::~B1CalorimeterSD()
{
}

void B1CalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  // eID=-1;
  // eEnergy=0.;
  // eTime=-1.;
  nPhoton=0;
  sipm_photon_num = 0;
  // Layer_num = 15;
  // Cell_num = 15;
}

G4bool B1CalorimeterSD::ProcessHits(G4Step* step, G4TouchableHistory*ROhist)
{  //since I am a SIPM sensor, no deposition here is considered!
/*
  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();
  
  // step length
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
  }
  if ( edep==0. && stepLength == 0. ) return false;

  // Get hit accounting data for this cell
*/
  // Get hit for total accounting

  //ok, here we directly call rootmgr ro record the optical photon 
  // G4cout << " [SD] ==> Catch one photon \n" << G4endl;
  // const G4Event* event =  G4RunManager::GetRunManager()->GetCurrentEvent();
  // if(event)
  // {
  //   eID=event->GetEventID();
  // }
  // eEnergy=step->GetTrack()->GetTotalEnergy();
  // eTime=step->GetPreStepPoint()->GetGlobalTime();
  G4String name = step->GetTrack()->GetDefinition()->GetParticleName();
  // G4cout<< "name is "<< name << G4endl;
  if (name=="opticalphoton")
  {
      auto Layer_id = step->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(2);
      auto Cell_id = step->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(1);
      // auto Sipm_id =  step->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber();
      // std:: cout << "sipm id: " << sipmid <<std::endl;
      // sipm_photon[sipmid]++;
      photon_layer_cell[Layer_id][Cell_id] += 1;
  //     // nPhoton+=1;
  //     layer_id = step->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(3);
  //     // std::cout << "Photon layer id: " << layer_id << std::endl;
  //     photon_num.at(layer_id)+=1;
  //     xy_id = step->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(2);
  //     bar_id = step->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(1);
  // // bar_edep.push_back(eEnergy_scin);
  // //eEnergy_scin = step->GetTrack()->GetTotalEnergy();
  // // edep_perbar.at(copy_num) += eEnergy_scin;
  //   if(xy_id == 0){
  //       perbar_pho_x[layer_id][bar_id] += 1;
  //   }
  //   else if(xy_id==1){
  //      perbar_pho_y[layer_id][bar_id] += 1;
  //   }
  }

  //fEventAction->AddSipmEdep(eEnergy);
  //fEventAction->GetSipmTime(eTime);

  //fRootMgr->FillSipmPhoton(step->GetTrack()->GetTotalEnergy(),  step->GetPreStepPoint()->GetGlobalTime(), eID);
  return true;
}

void B1CalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
  // for(int i = 0 ; i < 20 ; i++){
  //       for (int j = 0 ; j < layer_num ; j++){
  //           pho_bar_x.push_back(perbar_pho_x[j][i]);
  //           // if(perbar_edep_x[j][i] < 0 )std::cout<<"barE0less exists" << "\n";
  //           perbar_pho_x[j][i] = 0;
  //       }
  //   }
  //    for(int i = 0 ; i < 20 ; i++){
  //       for (int j = 0 ; j < layer_num ; j++){
  //           pho_bar_y.push_back(perbar_pho_y[j][i]);
  //           perbar_pho_y[j][i] = 0 ;
  //       }
    // }
  // fRootMgr->FillSipmPhoton(step->GetTrack()->GetTotalEnergy(),  step->GetPreStepPoint()->GetGlobalTime(), eID);
  // for(int idx = 0 ; idx < 2; idx++){
  //   sipm_photon_num+=sipm_photon[idx];
  // }
  // sipm_photons.emplace_back(sipm_photon_num);
  // sipm_photon_num = 0;
  // fRootMgr->FillSipmPhoton(sipm_photons);
  
  // sipm_photons.clear();
  // photon_num.clear();
  // pho_bar_x.clear();
  // pho_bar_y.clear();
  // std::fill(photon_num.begin(), photon_num.end(), 0);
  // std::fill(pho_bar_x.begin(), pho_bar_x.end(), 0);
  // std::fill(pho_bar_y.begin(), pho_bar_y.end(), 0);
  // std::fill(sipm_photon, sipm_photon+4, 0);
  // eEnergy=0;
  // eTime=0;
  // eID=0;
    for(G4int i = 0; i < fNofLayers; i++)
    {
        for (G4int j = 0; j < fNofCells; j++)
        {
            photon_num.emplace_back(photon_layer_cell[i][j]);
            photon_layer_cell[i][j] = 0;
        }
    }
    fRootMgr->FillSipmPhoton(photon_num);
    photon_num.clear();
    std::fill(photon_num.begin(), photon_num.end(), 0);
}



