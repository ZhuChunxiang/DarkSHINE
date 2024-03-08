#include "B1ScintillatorSD.hh"
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


B1ScintillatorSD::B1ScintillatorSD(const G4String& name, const G4String& hitsCollectionName, RootManager *rootMng, G4int NofLayers, G4int NofCells) : G4VSensitiveDetector(name),
FNofLayers(NofLayers),
FNofCells(NofCells),
fRootMgr(rootMng)
{
    collectionName.insert(hitsCollectionName);
    G4cout << "SensitiveDetector Processed Successfully " << G4endl;
}

B1ScintillatorSD::~B1ScintillatorSD()
{
}

void B1ScintillatorSD::Initialize(G4HCofThisEvent* hce)
{
  eID=-1;
  eEnergy=0.;
  eTime=-1.;
  // nPhoton=0;
  // sipm_photon_num = 0;
  // Layer_n = 15;
  // Cell_n = 15;

  // Create hits collection
  fHitsCollection = new B1ScintHitsCollection(SensitiveDetectorName, collectionName[0]);
  
  // Add this collection in hce
  G4int fhcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(fhcID, fHitsCollection);
  
  // Create hits
  // fill calorimeter hits with zero energy deposition
  for (G4int i = 0; i < (FNofLayers * FNofCells + 1); i++)
  {
      fHitsCollection->insert(new B1ScintHit());
  }
}

G4bool B1ScintillatorSD::ProcessHits(G4Step* step, G4TouchableHistory*ROhist)
{ 
  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();
  
  // step length
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
  }
  if ( edep==0. && stepLength == 0. ) return false;

  // Get hit accounting data for this cell
  auto iLayer = step->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(2);
  auto iCell = step->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(1);
  // Energy_dep_layer_cell[iLayer][iCell] += edep;
  G4int hit_id = iLayer * FNofCells + iCell;
  auto hit_ = (*fHitsCollection)[hit_id]; // 解指针，判断放入哪个hit_id空间
  if (!hit) {
      G4ExceptionDescription msg;
      msg << "Cannot access hit " << copyNo;
      G4Exception("CalorimeterSD::ProcessHits()",
          "MyCode0004", FatalException, msg);
  }
  
  // Get hit for total accounting
  auto hit_ = (*fHitsCollection)[fHitsCollection->entries()-1];
  
  // Add values
  hit->Add(edep, stepLength);
  hitTotal->Add(edep, stepLength);

  // ok, here we directly call rootmgr ro record the optical photon 
  // G4cout << " [SD] ==> Catch one photon \n" << G4endl;
  // const G4Event* event =  G4RunManager::GetRunManager()->GetCurrentEvent();
  // if(event)
  // {
  //   eID=event->GetEventID();
  // }
  // eEnergy=step->GetTrack()->GetTotalEnergy();
  // eTime=step->GetPreStepPoint()->GetGlobalTime();
  // G4String name = step->GetTrack()->GetDefinition()->GetParticleName();
  // G4cout<< "name is "<< name << G4endl;
  // if (name=="opticalphoton")
  // {
      // auto Layer_id = step->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(2);
      // auto Cell_id = step->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(1);
      // auto Sipm_id =  step->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber();
      // std:: cout << "sipm id: " << sipmid <<std::endl;
      // sipm_photon[sipmid]++;
      // photon_layer_cell[Layer_id][Cell_id] += 1;
  //     // nPhoton+=1;
  //     layer_id = step->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(3);
  //     // std::cout << "Photon layer id: " << layer_id << std::endl;
  //     photon_num.at(layer_id)+=1;
  //     xy_id = step->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(2);
  //     bar_id = step->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(1);
  // // bar_edep.push_back(eEnergy_scin);
  // // eEnergy_scin = step->GetTrack()->GetTotalEnergy();
  // // edep_perbar.at(copy_num) += eEnergy_scin;
  //   if(xy_id == 0){
  //       perbar_pho_x[layer_id][bar_id] += 1;
  //   }
  //   else if(xy_id==1){
  //      perbar_pho_y[layer_id][bar_id] += 1;
  //   }
  // }

  //fEventAction->AddSipmEdep(eEnergy);
  //fEventAction->GetSipmTime(eTime);

  //fRootMgr->FillSipmPhoton(step->GetTrack()->GetTotalEnergy(),  step->GetPreStepPoint()->GetGlobalTime(), eID);
  return true;
}

void B1ScintillatorSD::EndOfEvent(G4HCofThisEvent*)
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
    for(G4int i = 0; i < FNofLayers; i++)
    {
        for (G4int j = 0; j < FNofCells; j++)
        {
            Energy_dep_per_cell.emplace_back(Energy_dep_layer_cell[i][j]);
            Energy_dep_layer_cell[i][j] = 0;
        }
    }
    fRootMgr->FillScinEdep(Energy_dep_per_cell);
    Energy_dep_per_cell.clear();
    std::fill(Energy_dep_per_cell.begin(), Energy_dep_per_cell.end(), 0);
}



