#include "RootManager.hh"

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TString.h"

#include <map>
#include <vector>


//....ooo0000ooo........ooo0000ooo.......ooo0000ooo.........ooo0000ooo.......
//....ooo0000ooo........ooo0000ooo.......ooo0000ooo.........ooo0000ooo.......

RootManager::RootManager()
{
    outfilename = "output.root";
    outTreename = "darkPhoton";

}

void RootManager::initialize()
{
    Edep = 0;
    WLSEdep = 0;
    NFiberOP = 0;
    NScintOP = 0;
    fStart = 0;
    fEvtNb = 0;
    fStepTime.clear();
}

RootManager::~RootManager()
{
     //if(rootFile) rootFile->Close();
     //delete rootFile;
}

void RootManager::book()
{
     G4String filename = outfilename;
     G4String treename = outTreename;
     G4cout<<filename<<G4endl;
     rootFile = new TFile(filename, "RECREATE");

     if (!rootFile)
     {
         G4cout << " RootManager::book :"
             << " problem creating the ROOT TFile "
             << G4endl;
         return;
     }

     tr = new TTree(treename, "Dark_Photon");
     if (!tr)
     {
         G4cout << "RootManager::book :"
                << "problem creating the TTree Object "
                << G4endl;
         return;
     }
     else {
        //  tr->Branch("depEnergy", &Edep, "depEnergy/D");
        //  tr->Branch("WLSdepEnergy", &WLSEdep, "WLSdepEnergy/D");
        //  tr->Branch("EventNumber", &fEvtNb, "fEvtNb/I");
        //  tr->Branch("NFiberOP", &NFiberOP, "NFiberOP/I");
        //  tr->Branch("NScintOP", &NScintOP, "NScintOP/I");
        //  tr->Branch("ParticleID",&PID, "ParticleID/I");
        //  tr->Branch("ParticleTime", &fStepTime);
        //  tr->Branch("Ab_Edep",&Ab_edep,"Ab_Edep/D");

        //  tr->Branch("EnergyDep_per_bar",&edep_per_scintbar);
        //  tr->Branch("TotalEdep_scin",&TotalEdep_scin,"TotalEdep_scin/D");
        
        // tr->Branch("layer_Edep",&layer_Edep);
        // tr->Branch("scin_Edep",&scin_Edep);
        // tr->Branch("edep_Bar_x",&edep_Bar_x);
        // tr->Branch("edep_Bar_y",&edep_Bar_y);
        

        //  record in sipm SD
        //  tr->Branch("SipmPhoton_energy", &SipmPhoton_E, "SipmPhoton_energy/D");
        //  tr->Branch("SipmPhoton_time", &SipmPhoton_T, "SipmPhoton_time/D");
        //  tr->Branch("SipmPhoton_eventID", &SipmPhoton_eID, "SipmPhoton_eventID/I");
        //  tr->Branch("SipmPhoton_Number", &SipmPhoton_No, "SipmPhoton_Number/I");
        //  tr->Branch("Photon_num",&photon_);
        //  tr->Branch("Photon_num_x",&pho_num_x);
        //  tr->Branch("Photon_num_y",&pho_num_y);
        //  SipmPhoton_No = 0;
         for(Int_t i = 0; i < layer_num; i++)
         {
             for (Int_t j = 0; j < cell_num; j++)
             {
                 s_name = s1 + std::to_string(i+1) + s2 + std::to_string(j+1);
                 tr->Branch(s_name, &Cell_photon[i * cell_num + j]);
             }
         }
         
     }

}

void RootManager::FillSim(double Energy, double WLSEnergy, int EventNb, int OPinFiber, int OPinScint, 
                           int pid, std::vector<double> Partime , double ab_edep )
{
    //this is for deposition-type information record: need accumulate in all the steps
    Edep     = Energy;
    WLSEdep  = WLSEnergy;
    fEvtNb   = EventNb;
    NFiberOP = OPinFiber;
    NScintOP = OPinScint;
    PID      = pid;
    fStepTime = Partime;
    Ab_edep  = ab_edep;
}

// void RootManager::FillScinEdep(std::vector<double>&layer_edep , std::vector<double>&edep_bar_x , std::vector<double>&edep_bar_y){
//     // edep_per_scintbar = edep_perbar;
//     // TotalEdep_scin = 0;
//     // for(it = edep_perbar.begin(); it != edep_perbar.end() ; it++){
//     //     TotalEdep_scin += *it;
//     //     if( *it != 0)
//     //     G4cout << "edep_perbar:"<<*it <<G4endl;
//     // }
//     // G4cout <<"TotalEdep_scin: " <<TotalEdep_scin <<G4endl;
//     // G4cout << "[Root Manager] ==> Deposited energy from per scintillator bar has been recorded \n "<<G4endl;
   
//     scin_Edep  = layer_edep;
//     edep_Bar_x  = edep_bar_x;
//     edep_Bar_y  = edep_bar_y;
 
// }


void RootManager::FillSipmPhoton( std::vector<int>& sipm_photons)
{
    //this is for deposition-type information record: need accumulate in all the steps
   
    photon_ = sipm_photons;

    // for(int idx = 0; idx < sipm_photons.size(); idx++)
    // G4cout << " [Root Manager] Sipm: "<<idx<<" ==> Catch " << sipm_photons[idx]<<  " photon \n" << G4endl;
    for(Int_t i = 0; i < Num; i++)
    {
        Cell_photon[i].emplace_back(photon_.at(i));
    }
}

void RootManager::Fill()
{
    tr->Fill();
}

void RootManager::saveTree() 
{
     if (rootFile)
     {
         rootFile->cd();
         rootFile->Write("", TObject::kOverwrite);
         G4cout << " [Root Manager] ==> Simulation Tree is saved \n" << G4endl;
     }
     rootFile->Close();
}
