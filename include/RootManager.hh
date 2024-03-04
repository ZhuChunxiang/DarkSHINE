#ifndef ROOTMANAGER_H
#define ROOTMANAGER_H

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TArrayD.h"

#include "TH3D.h"
#include "TH3I.h"
#include "TH2D.h"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <map>
#include <type_traits>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TFile;

class TTree;

class TRandom3;

class RootManager {
    public:
        RootManager();

        ~RootManager();

        void book();

        void bookCollection(const G4String&);

        void saveTree();

        void initialize();

        /* set methods */
        void SetOutFileName(const G4String& in) { outfilename = in; };

        void SetStartID(int id) { fStart = id; };

        void SetNbEvent(int id) { fEvtNb = id; };

        /* fill methods */
        void FillSim(Double_t Edep, Double_t WLSEdep, int fEvtNb, int OPinFiber, int OPinScint, int pid, std::vector<double> Partime,Double_t ab_edep); //fill from endOfEvent
        void FillScinEdep(std::vector<double>& cell_edep); // fill from scintillator
        void FillSipmPhoton(std::vector<int>& sipm_photons); //fill from SD
      //   void FillScinEdep(std::vector<double>&layer_edep , std::vector<double>&edep_bar_x , std::vector<double>&edep_bar_y);
        void Fill();

     public:
        G4String outfilename;
        G4String outTreename;
        TFile *rootFile;
        TTree *tr;
        

     private:
        Int_t fStart; // Run Number, Initialized to 0. Set method: RootManager::SetStartID(int id)
        Int_t fEvtNb; // Event Number, Initialized to 100000. Set method: RootManager::SetNbEvent(int id)
        Int_t PID;    //Particle IDentification Number 
        Int_t NFiberOP;
        Int_t NScintOP;
        Double_t Edep;
        Double_t WLSEdep;
        Double_t TotalEdep_scin; //Total energy for the whole scintillator
        Double_t Ab_edep;//Energy of every single layer absorber
        std::vector<double> edep_per_scintbar; //Deposited energy for per scintillator bar
        std::vector<double> fStepTime;
        std::vector<double>::iterator it;

        //Sipm SD needed variable
        Double_t SipmPhoton_E;
        Double_t SipmPhoton_T;
        Int_t SipmPhoton_eID;
        Int_t SipmPhoton_No; 
        Int_t layer_num = 5;
        Int_t cell_num = 15;
        Int_t Num = layer_num * cell_num;
        std::string s1 = "Photon_num_Layer_";
        std::string s2 = "_Cell_";
        std::string s_name;
        std::vector<double> energy_dep_;
        std::vector<int> photon_;
        std::vector<std::vector<double> > Cell_energy_dep;
        std::vector<std::vector<int> > Cell_photon;
      //   std::vector<int> pho_num_x;
      //   std::vector<int> pho_num_y;
      //   std::vector<int> layer_id;
      //   std::vector<double> layer_Edep;
      //   std::vector<double> scin_Edep;
      //   std::vector<double> edep_Bar_x;
      //   std::vector<double> edep_Bar_y;
        
};

#endif
