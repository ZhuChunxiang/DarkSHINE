
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "rnd.h"
// #include <iostream>
// #include <cmath>
#include <TROOT.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <string.h>
#include <stdio.h>
#include "TSystem.h"
#include <vector>

void rnd::Loop(string filename)
{
Int_t           nbytes=0 , nb=0; 



initialize(filename);
// TFile *f = TFile::Open(filename,"READ");
// f->GetObject("darkPhoton" , fChain);

// fChain->SetMakeClass(1);

// fChain->SetBranchAddress("layer_Edep", &layer_Edep, &b_layer_Edep);
// fChain->SetBranchAddress("Photon_num", &Photon_num, &b_Photon_num);


// TTree *fChain=(TTree*)large_file->Get("darkPhoton") ;//so here,should I input the tree name of my root file?

//  DAOD_tree large_DAOD(fChain_DAOD);//  lage_DAOD.Loop();  //134,what this line work?

 if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();  
   }

Long64_t nentries = fChain->GetEntriesFast();
for(Long64_t jentry=0 ; jentry<nentries ; jentry++){
   Long64_t ientry = LoadTree(jentry);
   if(ientry < 0)
      break;
   nb = fChain->GetEntry(jentry);
   nbytes += nb;

   execute();
}
   finalize();

}
void exectute()
{

}