#define dijet_cxx
#define vbswy_cxx
#define DAOD_tree_cxx

#include <iostream>
#include <DAOD_tree.h>
#include <dijet.h>
#include <vbswy.h>
#include "TROOT.h"
#include "TH1F.h"
#include "THStack.h"
#include "TRatioPlot.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include <TStyle.h>
#include <TTree.h>
#include <TText.h>
#include <dirent.h>
#include <string>
#include <boost/algorithm/string.hpp>
#include <vector>
#include "TMath.h"
#include "TString.h"
#include <ctime>
#include <TLorentzVector.h>
#include <sstream>
#include <math.h>


#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

using namespace std;

double mpi_pi(double angle){

  while (angle >= TMath::Pi()) angle -= 2*TMath::Pi();
  while (angle < -1*TMath::Pi()) angle += 2*TMath::Pi();

  return angle;

}

TString *testGetOpt(int argc, char *argv[]){
  int opt;
  const char *optstring = "w:m:j:r:p:n:c:f:";
  static TString outc[8]={"20","20","30","0.4","0", "10","0","./"};

  while ((opt = getopt(argc, argv ,optstring))!=-1) {
    printf("opt = %c\n", opt);
    printf("optarg = %s\n", optarg);
    printf("optind = %d\n", optind);
    printf("argv[optind - 1] = %s\n", argv[optind - 1]);
    switch(opt){
      case 'w':
        printf("wmt is %s\n\n", optarg);
        outc[0] = optarg;
        break;
      case 'm':
        printf("met is %s\n\n", optarg);
        outc[1] = optarg;
        break;
      case 'j':
        printf("jet is %s\n\n", optarg);
        outc[2] = optarg;
        break;
      case 'r':
        printf("dr is %s\n\n", optarg);
        outc[3] = optarg;
        break;
      case 'p':
        printf("dphi is %s\n\n", optarg);
        outc[4] = optarg;
        break;
      case 'n':
        printf("njet is %s\n\n", optarg);
        outc[5] = optarg;
        break;
      case 'c':
        printf("data is %s\n\n", optarg);
        outc[6] = optarg;
        break;
      case 'f':
        printf("dir_out is %s\n\n", optarg);
        outc[7] = optarg;
        break;
      default:
        printf("You should look for help! ");
        exit(1);
        break; 
    }
  }
  return outc;
}


bool SplitString(const TString& theOpt, const char separator,
                 int &id)
{
    TString splitOpt(theOpt);
    splitOpt.ReplaceAll("\n", " ");
    splitOpt = splitOpt.Strip(TString::kBoth, separator); //Split splitOpt into TSubStrings by seperator

    while(splitOpt.Length() > 0)
    {
        if(!splitOpt.Contains(separator))
        {
            if(splitOpt.IsFloat()){
                id = splitOpt.Atoi();
                return true;
            }
            break;
        }
        //If the TSubString contains the seperator, fill v from 1st char to 1st seperator
        else
        {
            TString toSave = splitOpt(0, splitOpt.First(separator));
            //cout << toSave << endl;
            if(toSave.IsFloat()){
                id = toSave.Atoi();
                return true;
            }
            splitOpt = splitOpt(splitOpt.First(separator), splitOpt.Length()); //Discard the chars added to the vector
        }
        splitOpt = splitOpt.Strip(TString::kLeading, separator);
    }

    return false;
}

Double_t ScaleFactor(TString f_name, bool is_dijet ){  
  TFile *large_file = TFile::Open(f_name,"READ"); //TFile* myfile = new TFile("","OPEN");

  TTree *fChain_DAOD=(TTree*)large_file->Get("DAOD_tree");
  DAOD_tree large_DAOD(fChain_DAOD);//  lage_DAOD.Loop(); 

  Long64_t nentries_DAOD = fChain_DAOD->GetEntriesFast();
  Double_t lumi=0, xsec=0, geff=0, kfac=0, sw=0;
  for (Long64_t jentry=0; jentry<nentries_DAOD;jentry++) {
    fChain_DAOD->GetEntry(jentry);
    if (jentry==0) {
      lumi = large_DAOD.luminosity;
      xsec = large_DAOD.xsection;
      geff = large_DAOD.geneff;
      kfac = large_DAOD.kfactor;
    }
    sw += large_DAOD.sum_of_weights;
  }
  if (is_dijet)
    return lumi*xsec*geff*kfac/16000000/1.3;
  else
    return lumi*xsec*geff*kfac/sw;
}


Double_t FakeFactor(double met_modify, double wmt_modify, double jet_modify, double dr_modify, double dphi_modify, double njet_modify, int run_data){
  vector<TString> dir_in;
  dir_in.push_back("/lustre/collider/wangxi/data_dijet/mc16de_data/"); // end with /
  dir_in.push_back("/lustre/collider/wangxi/data_dijet/mc16a_data/"); // end with /
  //ihep dir_in.push_back("/publicfs/atlas/atlasnew/SM/VBS/xiw/mini_files/dijet/data_mc1516/"); // end with /
  //ihep dir_in.push_back("/publicfs/atlas/atlasnew/SM/VBS/xiw/mini_files/dijet/dijet_v7_mc1718/"); // end with /
  //ihep dir_in.push_back("/publicfs/atlas/atlasnew/SM/VBS/xiw/mini_files/dijet/data1718/"); // end with /
  //ihep dir_in.push_back("/publicfs/atlas/atlasnew/SM/VBS/xiw/download_merge/dijet_new/"); // end with /

  double n_dijet = 0;
  double n_bkg = 0;
  double n_dijet_tight = 0;
  double n_bkg_tight = 0;
  TH1F * dr_dijet = new TH1F("dr_dijet", "", 50 ,0 , 5);
  TH1F * dr_bkg = new TH1F("dr_bkg", "", 50 ,0 , 5);
  TH1F * dphi_dijet = new TH1F("dphi_dijet", "", 40 ,0 , 4);
  TH1F * dphi_bkg = new TH1F("dphi_bkg", "", 40 ,0 , 4);

  for (int di = 0; di< dir_in.size(); di++){
    if (auto dir = opendir(dir_in[di])) {
      while (auto f = readdir(dir)) {
        if (!f->d_name || f->d_name[0] == '.')
            continue; // Skip everything that starts with a dot

        TString f_name = f->d_name;
        printf("File: %s\n", f->d_name);

        std::string fileName_str = f->d_name;
        Int_t isData =0;
        if (boost::algorithm::contains(fileName_str, "period")) 
          isData = 1;

        Double_t sf = ScaleFactor(dir_in[di]+f_name, false);
        Double_t dijet_sf = ScaleFactor(dir_in[di]+f_name, true);

        TFile *large_file = TFile::Open(dir_in[di]+f_name,"READ"); //TFile* myfile = new TFile("","OPEN");
                
        TTree *fChain_dijet;
        if (run_data==0){
          int id;
          if(SplitString(f_name, '.', id)){
            if(!(id >=364701 && id <= 364712))
              continue; 
            else
              cout<<id<<endl;
          }
          else
            continue;
        }
                                                           //change 2
        fChain_dijet=(TTree*)large_file->Get("dijet_tree");
        dijet dijet(fChain_dijet); // large_dijet.Loop();


        //Long64_t nentries_dijet = fChain_dijet->GetEntriesFast();
        Long64_t nentries_dijet = fChain_dijet->GetEntries();         //change 1
        cout << "dijet tree nentries: " << nentries_dijet << endl;

        for (Long64_t jentry=0; jentry<nentries_dijet;jentry++) {
          fChain_dijet->GetEntry(jentry);
          Double_t mc_weight = sf * dijet.weight_all*dijet.weight_m *dijet.weight_e;    //        change 3
          Double_t mc_dijet_weight = dijet_sf * dijet.weight_all*dijet.weight_m *dijet.weight_e;

          //int baseline_mu = dijet.mu_pt>20 && fabs(dijet.mu_eta)<2.5 && fabs(dijet.mu_z0)<0.5 && dijet.mu_idTight && dijet.met<20 && dijet.w_mt<20 && ( (dijet.j0_pt>30 && dijet.j0lep_dr2>3) || (dijet.j1_pt>30 && dijet.j_n >= 2 && dijet.j1lep_dr2>3));
          TLorentzVector j0, j1, mu;
          mu.SetPtEtaPhiE(dijet.mu_pt, dijet.mu_eta, dijet.mu_phi, dijet.mu_e);
          j0.SetPtEtaPhiE(dijet.j0_pt, dijet.j0_eta, dijet.j0_phi, dijet.j0_e);
          double j0mu_dphi = j0.DeltaPhi(mu);
          double muj0_dphi = mu.DeltaPhi(j0);
          j1.SetPtEtaPhiE(dijet.j1_pt, dijet.j1_eta, dijet.j1_phi, dijet.j1_e);
          double j1mu_dphi = j1.DeltaPhi(mu);

          //int baseline_mu = dijet.mu_pt>20 && fabs(dijet.mu_eta)<2.5 && fabs(dijet.mu_z0)<0.5 && dijet.mu_idTight && dijet.met<20 && dijet.w_mt<20 && ( (dijet.j0_pt>30 && dijet.j0lep_dr2>1) || (dijet.j1_pt>30 && dijet.j_n == 2 && dijet.j1lep_dr2>1));
          //int baseline_mu = dijet.mu_pt>20 && fabs(dijet.mu_eta)<2.5 && fabs(dijet.mu_z0)<0.5 && dijet.mu_idTight && dijet.met<20 && dijet.w_mt<20 && dijet.j0_pt>30 && dijet.j0lep_dr2>0.4 ; //to be used

          //int baseline_mu = dijet.mu_pt>20 && fabs(dijet.mu_eta)<2.5 && fabs(dijet.mu_z0)<0.5 && dijet.mu_idTight && dijet.met<20 && dijet.w_mt<20 && dijet.j0_pt>30 && dijet.j0lep_dr2>0.4; //to be used

          double j0met_dphi = fabs(mpi_pi(dijet.j0_phi - dijet.met_phi)); 
          double j1met_dphi = fabs(mpi_pi(dijet.j1_phi - dijet.met_phi)); 
          int pass_dphi = ( fabs(j0met_dphi > 0.4) && fabs(j1met_dphi > 0.4) );

          //int baseline_mu = dijet.mu_n==1 && dijet.mu_pt>20 && fabs(dijet.mu_eta)<2.5 && dijet.met<met_modify && dijet.w_mt<wmt_modify && dijet.j0_pt>jet_modify && fabs(dijet.j0lep_dr)> 3.2 && fabs(j0mu_dphi)>3;// cut best
          //int baseline_mu = dijet.mu_n==1 && dijet.pass_triggermatch_mu && dijet.mu_pt>20 && fabs(dijet.mu_eta)<2.5 && dijet.met<met_modify && dijet.w_mt<wmt_modify && dijet.j0_pt>jet_modify && dijet.gam_n==0 && dijet.mu_passIP && dijet.mu_isoFCLoose && dijet.mu_idLoose && dijet.j1_pt>jet_modify;
          //int optimize_mu = ((fabs(dijet.j0lep_dr)>dr_modify && fabs(j0mu_dphi)>dphi_modify) || (fabs(dijet.j1lep_dr)>dr_modify && fabs(j1mu_dphi)>dphi_modify)) && pass_dphi;
          //int baseline_mu = dijet.mu_n==1 && dijet.pass_triggermatch_mu && dijet.mu_pt>20 && fabs(dijet.mu_eta)<2.5 && dijet.met<met_modify && dijet.w_mt<wmt_modify && dijet.j0_pt>jet_modify && dijet.gam_n==0 && dijet.mu_passIP && dijet.mu_isoFCLoose && dijet.mu_idLoose;
          int baseline_mu = dijet.mu_n==1 && dijet.pass_triggermatch_mu && dijet.mu_pt>20 && fabs(dijet.mu_eta)<2.5 && dijet.met<met_modify && dijet.w_mt<wmt_modify && dijet.j0_pt>jet_modify && dijet.gam_n==0 && dijet.mu_passIP;
          int optimize_mu = fabs(dijet.j0lep_dr)>dr_modify && fabs(j0mu_dphi)>dphi_modify;
              
          int is_dijet_mc = 0, isData = 0, is_sherpa = 0, is_bkg_mc = 0;

          if ((dijet.dsid == 364250 || dijet.dsid == 364253) || (dijet.dsid >= 364100 && dijet.dsid <= 364141) || (dijet.dsid >= 364156 && dijet.dsid <= 364183) || (dijet.dsid >= 308092 && dijet.dsid <= 308094))
            is_sherpa =1;
          else is_sherpa = 0;
          //if (is_sherpa == 1 && (fabs(dijet.weight_all*dijet.weight_m *dijet.weight_e)>100)) continue;
          //if (is_sherpa == 1 && (fabs(dijet.weight_all*dijet.weight_m *dijet.weight_e)>100)) continue;
          //if (is_sherpa == 1 && (fabs(mc_weight)>50)) continue;

          if (!baseline_mu) continue;
          if (dijet.dsid >= 364702 && dijet.dsid <= 364712) {
              dr_dijet->Fill(fabs(dijet.j0lep_dr), mc_weight); 
              dphi_dijet->Fill(fabs(j0mu_dphi), mc_weight); 
          }

          if (dijet.dsid == 364250 || dijet.dsid == 364253) {dr_bkg->Fill(fabs(dijet.j0lep_dr), mc_weight); dphi_bkg->Fill(fabs(j0mu_dphi),mc_weight);}
          if (dijet.dsid >= 364100 && dijet.dsid <= 364141) {dr_bkg->Fill(fabs(dijet.j0lep_dr), mc_weight); dphi_bkg->Fill(fabs(j0mu_dphi),mc_weight);}
          if (dijet.dsid >= 364156 && dijet.dsid <= 364183) {dr_bkg->Fill(fabs(dijet.j0lep_dr), mc_weight); dphi_bkg->Fill(fabs(j0mu_dphi),mc_weight);}
          if (dijet.dsid >= 410400 && dijet.dsid <= 410700) {dr_bkg->Fill(fabs(dijet.j0lep_dr), mc_weight); dphi_bkg->Fill(fabs(j0mu_dphi),mc_weight);}
          if (dijet.dsid >= 308092 && dijet.dsid <= 308094) {dr_bkg->Fill(fabs(dijet.j0lep_dr), mc_weight); dphi_bkg->Fill(fabs(j0mu_dphi),mc_weight);}

          if (!optimize_mu) continue;

          //if (dijet.dsid == 4294967295) {h_data->Fill(dijet.jj_m, 1); isData=1;}  //change 4
          if (dijet.dsid >1e7) {h_data->Fill(dijet.jj_m, 1); isData=1;}  //change 4
          else isData=0;

          if (dijet.dsid >= 364701 && dijet.dsid <= 364712) {
            is_dijet_mc=1;
            if(dijet.dsid == 364701){
              h_dijet701->Fill(dijet.jj_m, mc_weight);
            }
            else{
              dijet_origin->Fill(dijet.mu_truth_origin, mc_weight);
              h_dijet->Fill(dijet.jj_m, mc_weight); dijet_event_yield->Fill(mc_weight);
              n_dijet+=mc_weight;
            }
          }
          else is_dijet_mc = 0;

          if      (dijet.dsid == 364250 || dijet.dsid == 364253) {is_bkg_mc = 1; h_diboson->Fill(dijet.jj_m, mc_weight);n_bkg+=mc_weight;}
          else if (dijet.dsid >= 364100 && dijet.dsid <= 364141) {is_bkg_mc = 1; h_zjets->Fill(dijet.jj_m, mc_weight);n_bkg+=mc_weight;}
          else if (dijet.dsid >= 364156 && dijet.dsid <= 364183) {is_bkg_mc = 1; h_wjets->Fill(dijet.jj_m, mc_weight);n_bkg+=mc_weight;}
          else if (dijet.dsid >= 410400 && dijet.dsid <= 410700) {is_bkg_mc = 1; h_top->Fill(dijet.jj_m, mc_weight);n_bkg+=mc_weight;}
          else if (dijet.dsid >= 308092 && dijet.dsid <= 308094) {is_bkg_mc = 1; h_zjets_vbf->Fill(dijet.jj_m, mc_weight);n_bkg+=mc_weight;}
          else is_bkg_mc = 0;
          if (dijet.mu_isoFCTight && dijet.mu_passIP && dijet.mu_idTight ){
          //if (dijet.mu_isoFCTight && dijet.mu_idTight){
          //if (dijet.mu_isoFCTight && pass_dphi){
          //if (dijet.mu_isoFCTight){
            //isData? h_pt_tight->Fill(dijet.mu_pt, 1): h_pt_tight->Fill(dijet.mu_pt, -mc_weight);
            //isData? n_data_tight +=1: n_mc_tight += mc_weight;
            if(isData){
              n_data_tight +=1;
              h_pt_tight->Fill(dijet.mu_pt, 1);
              h_eta_tight->Fill(fabs(dijet.mu_eta), 1);
              h2_pt_tight->Fill(dijet.mu_pt, 1);
              h2_eta_tight->Fill(fabs(dijet.mu_eta), 1);
              t_h2->Fill(dijet.mu_pt, fabs(dijet.mu_eta), 1);
              t_h2_bin2->Fill(dijet.mu_pt, fabs(dijet.mu_eta), 1);
              h_mupt_data->Fill(dijet.mu_pt, 1);
              if (fabs(dijet.mu_eta)<1.05)
                h_mupt_loweta_data->Fill(dijet.mu_pt,1);
            }
            else if(is_dijet_mc){
              if(dijet.dsid==364701){
                h_dijet701_tight->Fill(dijet.mu_pt, mc_weight);
                h_dijet701_tight_jet->Fill(dijet.j0_pt, mc_weight);
              }
              else{
                h_dijet_pt_tight->Fill(dijet.mu_pt, mc_weight);
                h_dijet_eta_tight->Fill(fabs(dijet.mu_eta), mc_weight);
                t_h2_dijet->Fill(dijet.mu_pt, fabs(dijet.mu_eta), mc_weight);
                h2_dijet_pt_tight->Fill(dijet.mu_pt, mc_weight);
                n_dijet_tight += mc_weight;
                h_mupt_dijet->Fill(dijet.mu_pt, mc_weight);
                if (fabs(dijet.mu_eta)<1.05)
                   h_mupt_loweta_dijet->Fill(dijet.mu_pt, mc_weight);
              }
            }
            else if(is_bkg_mc){

              if      (dijet.dsid == 364250 || dijet.dsid == 364253) {h_mupt_diboson->Fill(dijet.mu_pt, mc_weight);}
              else if (dijet.dsid >= 364100 && dijet.dsid <= 364141) {
                h_mupt_zjets->Fill(dijet.mu_pt, mc_weight);
                zjets_event_yield_tight->Fill(mc_weight);
              }
              else if (dijet.dsid >= 364156 && dijet.dsid <= 364183) {
                h_mupt_wjets->Fill(dijet.mu_pt, mc_weight);
                if (dijet.mu_pt>80 && (dijet.mu_eta <1.05))
                wjets_event_yield_tight->Fill(mc_weight);
              }
              else if (dijet.dsid >= 410400 && dijet.dsid <= 410700) {h_mupt_top->Fill(dijet.mu_pt, mc_weight);}
              else if (dijet.dsid >= 308092 && dijet.dsid <= 308094) {h_mupt_zjets_vbf->Fill(dijet.mu_pt, mc_weight);}

              if(fabs(dijet.mu_eta)<1.05 ){
                if      (dijet.dsid == 364250 || dijet.dsid == 364253) {h_mupt_loweta_diboson->Fill(dijet.mu_pt, mc_weight);}
                else if (dijet.dsid >= 364100 && dijet.dsid <= 364141) {
                  h_mupt_loweta_zjets->Fill(dijet.mu_pt, mc_weight);
                }
                else if (dijet.dsid >= 364156 && dijet.dsid <= 364183) {
                  h_mupt_loweta_wjets->Fill(dijet.mu_pt, mc_weight);
                  if (dijet.mu_pt>100) cout << "wjets dsid "<< dijet.dsid << ", pt " <<  dijet.mu_pt << ", eta " << dijet.mu_eta << "sf " << mc_weight<<endl;
                }
                else if (dijet.dsid >= 410400 && dijet.dsid <= 410700) {h_mupt_loweta_top->Fill(dijet.mu_pt, mc_weight);}
                else if (dijet.dsid >= 308092 && dijet.dsid <= 308094) {h_mupt_loweta_zjets_vbf->Fill(dijet.mu_pt, mc_weight);}
              }

              n_mc_tight += mc_weight;   
              h_pt_tight->Fill(dijet.mu_pt, -mc_weight);
              h_eta_tight->Fill(fabs(dijet.mu_eta), -mc_weight);
              h2_pt_tight->Fill(dijet.mu_pt, -mc_weight);
              h2_eta_tight->Fill(fabs(dijet.mu_eta), -mc_weight);
              t_h2->Fill(dijet.mu_pt, fabs(dijet.mu_eta), -mc_weight);
              t_h2_bin2->Fill(dijet.mu_pt, fabs(dijet.mu_eta), -mc_weight);
              n_bkg_tight += mc_weight;
            }
            else cout<< "unknown dsid: " << dijet.dsid <<endl;
          }
	  //else if (!dijet.mu_isoFCTight || (fabs(dijet.mu_d0sig)>3)){
          else{                               //change 4
            //isData? h_pt_antiTight->Fill(dijet.mu_pt, 1): h_pt_antiTight->Fill(dijet.mu_pt, -mc_weight);
            //isData? n_data_antiTight +=1: n_mc_antiTight += mc_weight;
            if(isData){
              n_data_antiTight +=1;
              h_pt_antiTight->Fill(dijet.mu_pt, 1);
              h_eta_antiTight->Fill(fabs(dijet.mu_eta), 1);
              h2_pt_antiTight->Fill(dijet.mu_pt, 1);
              h2_eta_antiTight->Fill(fabs(dijet.mu_eta), 1);
              anti_h2->Fill(dijet.mu_pt, fabs(dijet.mu_eta), 1);
              anti_h2_bin2->Fill(dijet.mu_pt, fabs(dijet.mu_eta), 1);
            }
            else if(is_dijet_mc){
              if(dijet.dsid==364701){
                h_dijet701_antiTight->Fill(dijet.mu_pt, mc_weight);
                h_dijet701_antiTight_jet->Fill(dijet.j0_pt, mc_weight);
              }
              else{
                h_dijet_pt_antiTight->Fill(dijet.mu_pt, mc_weight);
                h_dijet_eta_antiTight->Fill(fabs(dijet.mu_eta), mc_weight);
                anti_h2_dijet->Fill(dijet.mu_pt, fabs(dijet.mu_eta), mc_weight);
                h2_dijet_pt_antiTight->Fill(dijet.mu_pt, mc_weight);
              }
            }
            else if(is_bkg_mc){
              n_mc_antiTight += mc_weight;   
              h_pt_antiTight->Fill(dijet.mu_pt, -mc_weight);
              h_eta_antiTight->Fill(fabs(dijet.mu_eta), -mc_weight);
              h2_pt_antiTight->Fill(dijet.mu_pt, -mc_weight);
              h2_eta_antiTight->Fill(fabs(dijet.mu_eta), -mc_weight);
              anti_h2->Fill(dijet.mu_pt, fabs(dijet.mu_eta), -mc_weight);
              anti_h2_bin2->Fill(dijet.mu_pt, fabs(dijet.mu_eta), -mc_weight);

              if (dijet.dsid >= 364100 && dijet.dsid <= 364141) {
                zjets_event_yield_antiTight->Fill(mc_weight);
              }
              else if (dijet.dsid >= 364156 && dijet.dsid <= 364183) {
                if (dijet.mu_pt>80 && (dijet.mu_eta <1.05))
                wjets_event_yield_antiTight->Fill(mc_weight);
              }
            }
            else cout<< "unknown dsid: " << dijet.dsid <<endl;
          }
        }  
      }
    }
  }
  /* 
  cout << "tight muon bin 1: " << h_pt_tight->GetBinContent(1) << endl;
  cout << "antiTight muon bin 1: " << h_pt_antiTight->GetBinContent(1) << endl;
  cout << "antiTight muon: " << h_pt_antiTight->GetSum() << endl;
  cout << "tight data, mc: " << n_data_tight << ", " << n_mc_tight << endl;
  cout << "antiTight data, mc: " << n_data_antiTight << ", " << n_mc_antiTight <<endl;
  */
  
  //h_dijet_pt_antiTight->Add(h_dijet701_antiTight);
  //h_dijet_pt_tight->Add(h_dijet701_tight);
  //h_dijet->Add(h_dijet701);
//----------------------------------------------Draw fake factor-----------------------------------------------
  gStyle->SetOptStat(0);
  TH2F *ff2d =  (TH2F*)t_h2->Clone();
  TCanvas *c1=new TCanvas("c1","",1300,350);
  c1->Divide(3,1);
  c1->cd(3);
  gPad->SetLogx();
  c1->SetLogx();
  anti_h2->Draw("COLZTEXT e0");
  anti_h2->SetTitle("Number of Events with anti-tight muons");
  anti_h2->GetXaxis()->SetTitle("muon p_{T}[GeV]");
  anti_h2->GetYaxis()->SetTitle("muon #eta");
  anti_h2->GetXaxis()->SetRangeUser(20,300);
  //anti_h2->Write();

  c1->cd(2);
  c1->SetLogx();
  gPad->SetLogx();
  t_h2->Draw("COLZTEXT e0");
  t_h2->GetXaxis()->SetTitle("muon p_{T}[GeV]");
  t_h2->GetYaxis()->SetTitle("muon #eta");
  t_h2->SetTitle("Number of Events with tight muons");
  t_h2->GetXaxis()->SetRangeUser(20,300);
  //t_h2->Write();
//
  c1->cd();
  c1->SetLogx();
  gPad->SetLogx();
  ff2d->Divide(anti_h2);
  ff2d->Draw("COLZTEXT e0");
  ff2d->SetName("ff_2D");
  ff2d->SetTitle("fake factor");
  ff2d->GetXaxis()->SetRangeUser(20,300);
  ff2d->GetXaxis()->SetTitle("muon p_{T}[GeV]");
  ff2d->GetYaxis()->SetTitle("muon #eta");
  //ff2d->Write();
  c1->SaveAs(dir_out+"/2D_ff.pdf");
//---------------------------------------bin check---------------
  TH2F *ff2d_bin2 =  (TH2F*)t_h2_bin2->Clone();
  TCanvas *c12=new TCanvas("c1","",1300,350);

  c12->Divide(3,1);
  c12->cd(3);
  gPad->SetLogx();
  c12->SetLogx();
  anti_h2_bin2->Draw("COLZTEXT e0");
  anti_h2_bin2->SetTitle("Number of Events with anti-tight muons");
  anti_h2_bin2->GetXaxis()->SetTitle("muon p_{T}[GeV]");
  anti_h2_bin2->GetYaxis()->SetTitle("muon #eta");
  anti_h2_bin2->GetXaxis()->SetRangeUser(20,300);
  //anti_h2->Write();

  c12->cd(2);
  c12->SetLogx();
  gPad->SetLogx();
  t_h2_bin2->Draw("COLZTEXT e0");
  t_h2_bin2->GetXaxis()->SetTitle("muon p_{T}[GeV]");
  t_h2_bin2->GetYaxis()->SetTitle("muon #eta");
  t_h2_bin2->SetTitle("Number of Events with tight muons");
  t_h2_bin2->GetXaxis()->SetRangeUser(20,300);
  //t_h2->Write();

  c12->cd();
  c12->SetLogx();
  gPad->SetLogx();
  ff2d_bin2->Divide(anti_h2_bin2);
  ff2d_bin2->Draw("COLZTEXT e0");
  ff2d_bin2->SetName("ff_2D");
  ff2d_bin2->SetTitle("fake factor");
  ff2d_bin2->GetXaxis()->SetRangeUser(20,300);
  ff2d_bin2->GetXaxis()->SetTitle("muon p_{T}[GeV]");
  ff2d_bin2->GetYaxis()->SetTitle("muon #eta");
  //ff2d->Write();
  c12->SaveAs(dir_out+"/2D_ff_bin2.pdf");



  TCanvas *c5=new TCanvas("c5","",800,600);
  c5->cd();
  c5->SetLogx();
  //c5->SetLogy();
  h2_pt_antiTight->Draw("e0");
  h2_pt_antiTight->SetLineColor(4);
  h2_pt_antiTight->GetXaxis()->SetTitle("muon p_{T}[GeV]");
  h2_pt_antiTight->GetYaxis()->SetTitle("events/bin");
  //h2_pt_antiTight->GetYaxis()->SetRangeUser(0.1,20000000.0);
  h2_pt_antiTight->GetYaxis()->SetRangeUser(-1000,1000);
  h2_pt_tight->Draw("e0 same");
  h2_pt_tight->SetLineColor(2);

  TLegend *leg1=new TLegend(0.60,0.8,0.9,0.9);
  leg1->AddEntry(h2_pt_tight,"Tight");
  leg1->AddEntry(h2_pt_antiTight,"anti-Tight");
  leg1->Draw("same");
  c5->SaveAs(dir_out+"/mu_pt_distribution.pdf");

  TCanvas *c6=new TCanvas("c6","",800,600);
  c6->cd();
  c6->SetLogy();
  h2_eta_antiTight->Draw("e0");
  h2_eta_antiTight->SetLineColor(4);
  h2_eta_antiTight->GetXaxis()->SetTitle("muon |#eta|");
  h2_eta_antiTight->GetYaxis()->SetTitle("events/bin");
  h2_eta_antiTight->GetYaxis()->SetRangeUser(10.0,20000000.0);
  h2_eta_tight->Draw("e0 same");
  h2_eta_tight->SetLineColor(2);

  TLegend *leg2=new TLegend(0.60,0.8,0.9,0.9);
  leg2->AddEntry(h2_eta_tight,"Tight");
  leg2->AddEntry(h2_eta_antiTight,"anti-Tight");
  leg2->Draw("same");
  c6->SaveAs(dir_out+"/mu_eta_distribution.pdf");

  TCanvas *c66=new TCanvas("c66","",800,600);
  c66->cd();
  c66->SetLogy();
  dijet_event_yield->Draw();
  c66->SaveAs(dir_out+"/dijet_event_yield.png");

  TCanvas *c77 = new TCanvas("c77","",800,600);
  c77->cd();
  dr_dijet->Draw();
  dr_dijet->SetLineColor(2);
  dr_bkg->Draw("same");
  dr_bkg->SetLineColor(4);
  TLegend *leg77=new TLegend(0.60,0.8,0.9,0.9);
  leg77->AddEntry(dr_dijet,"dijet");
  leg77->AddEntry(dr_bkg,"other background");
  leg77->Draw("same");
  c77->SaveAs(dir_out+"/dr_distribution.pdf");

  TCanvas *c88 = new TCanvas("c88","",800,600);
  c88->cd();
  dphi_dijet->Draw();
  dphi_dijet->SetLineColor(2);
  dphi_bkg->Draw("same");
  dphi_bkg->SetLineColor(4);
  TLegend *leg88=new TLegend(0.60,0.8,0.9,0.9);
  leg88->AddEntry(dphi_dijet,"dijet");
  leg88->AddEntry(dphi_bkg,"other background");
  leg88->Draw("same");
  c88->SaveAs(dir_out+"/dphi_distribution.pdf");


  /*
  //TCanvas *c2=new TCanvas("c2","",800,600);   
  c2->cd();
  c2->SetLogx();
  h_pt_tight->Divide(h_pt_antiTight);
  h_dijet_pt_tight->Divide(h_dijet_pt_antiTight);
  h_pt_tight->Draw("e0");
  h_pt_tight->SetLineColor(2);
  h_dijet_pt_tight->Draw("e0 same");
  h_dijet_pt_tight->SetLineColor(4);
  h_pt_tight->GetXaxis()->SetTitle("muon p_{T}[GeV]");
  h_pt_tight->GetYaxis()->SetRangeUser(0,2);
  h_pt_tight->GetYaxis()->SetTitle("fake factor");
  leg_ff->AddEntry(h_pt_tight,"Dijet DD");
  leg_ff->AddEntry(h_dijet_pt_tight,"Dijet MC");
  c2->SaveAs(dir_out+"/mu_pt_ff.pdf");
  */
  
  TCanvas *c3=new TCanvas("c3","",800,600);   
  c3->cd();
  h_eta_tight->Divide(h_eta_antiTight);
  h_dijet_eta_tight->Divide(h_dijet_eta_antiTight);
  h_eta_tight->Draw("e0");
  h_dijet_eta_tight->Draw("e0 same");
  h_dijet_eta_tight->SetLineColor(4);
  h_eta_tight->GetXaxis()->SetTitle("muon |#eta|");
  h_eta_tight->GetYaxis()->SetRangeUser(0,1);
  h_eta_tight->GetYaxis()->SetTitle("fake factor");
  c3->SaveAs(dir_out+"/mu_eta_ff.pdf");

  TCanvas *c37 = new TCanvas("c7", "", 800, 600);
  c37->cd();
  h2_dijet_pt_antiTight->Draw("e0");
  h2_dijet_pt_antiTight->GetYaxis()->SetRangeUser(-100000,600000);
  h2_dijet_pt_antiTight->SetLineColor(4);
  h2_dijet_pt_tight->Draw("e0 same");
  h2_dijet_pt_tight->SetLineColor(2);
  TLegend *leg_d=new TLegend(0.60,0.8,0.9,0.9); 
  leg_d->AddEntry(h2_dijet_pt_antiTight, "antiTight");
  leg_d->AddEntry(h2_dijet_pt_tight, "tight");
  leg_d->Draw("same");
  c37->SaveAs(dir_out+"/dijet_distribution.png");

  TCanvas *c137 = new TCanvas("c7", "", 800, 600);
  c137->cd();
  h_dijet701_antiTight->Draw("e0");
  h_dijet701_antiTight->GetYaxis()->SetRangeUser(-100000,5000000);
  h_dijet701_antiTight->SetLineColor(4);
  h_dijet701_tight->Draw("e0 same");
  h_dijet701_tight->SetLineColor(2);
  TLegend *leg_d2=new TLegend(0.60,0.8,0.9,0.9); 
  leg_d2->AddEntry(h_dijet701_antiTight, "antiTight");
  leg_d2->AddEntry(h_dijet701_tight, "tight");
  leg_d2->Draw("same");
  c137->SaveAs(dir_out+"/dijet701_distribution.png");

  TCanvas *c138 = new TCanvas("c7", "", 800, 600);
  c138->cd();
  h_dijet701_antiTight_jet->Draw("e0");
  h_dijet701_antiTight_jet->GetYaxis()->SetRangeUser(-100000,5000000);
  h_dijet701_antiTight_jet->SetLineColor(4);
  h_dijet701_tight_jet->Draw("e0 same");
  h_dijet701_tight_jet->SetLineColor(2);
  TLegend *leg_d22=new TLegend(0.60,0.8,0.9,0.9); 
  leg_d22->AddEntry(h_dijet701_antiTight_jet, "antiTight");
  leg_d22->AddEntry(h_dijet701_tight_jet, "tight");
  leg_d22->Draw("same");
  c138->SaveAs(dir_out+"/dijet701_jet_distribution.png");

  TCanvas *c107 = new TCanvas("c107","", 800,600);
  c107->cd();
  c107->SetLogy();
  wjets_event_yield_tight->Draw();
  c107->SaveAs(dir_out+"/wjets_event_yield_tight.png");

  TCanvas *c108 = new TCanvas("c108","", 800,600);
  c108->cd();
  c108->SetLogy();
  wjets_event_yield_antiTight->Draw();
  c108->SaveAs(dir_out+"/wjets_event_yield_antiTight.png");

  TCanvas *c109 = new TCanvas("c109","", 800,600);
  c109->cd();
  c109->SetLogy();
  zjets_event_yield_tight->Draw();
  c109->SaveAs(dir_out+"/zjets_event_yield_tight.png");

  TCanvas *c110 = new TCanvas("c110","", 800,600);
  c110->cd();
  c110->SetLogy();
  zjets_event_yield_antiTight->Draw();
  c110->SaveAs(dir_out+"/zjets_event_yield_antiTight.png");

 /* 
  char nm[]="/lustre/collider/wangxi/dijet_fake/run15/optimization.txt";
  ofstream outfile(nm,ios::app|ios::out);
  outfile << dr_modify << " " << dphi_modify << " " << n_dijet/n_bkg << " " << n_dijet/sqrt(n_bkg) <<" "<< n_dijet/sqrt(n_dijet+n_bkg)<<" "<< n_dijet << " "<<n_bkg<< endl;
  outfile.close();
  
  char nm2[]="/lustre/collider/wangxi/dijet_fake/run15/optimization_tight.txt";
  ofstream outfile2(nm2,ios::app|ios::out);
  outfile2 << dr_modify << " " << dphi_modify << " " << n_dijet_tight/n_bkg_tight << " " << n_dijet_tight/sqrt(n_bkg_tight) <<" "<< n_dijet_tight/sqrt(n_dijet_tight+n_bkg_tight)<<" "<< n_dijet_tight << " "<<n_bkg_tight<< endl;
  outfile2.close();
 */ 
}
void initialize(){
  h_pt_tight->Sumw2();
  h_pt_antiTight->Sumw2();
  h_dijet_pt_antiTight->Sumw2();
  h_dijet_pt_tight->Sumw2();
  h_sp_pt_antiTight->Sumw2();
  h_sp_pt_tight->Sumw2();
  h_eta_tight->Sumw2();
  h_eta_antiTight->Sumw2();
  t_h2->Sumw2();
  anti_h2->Sumw2();  
}
void ratioplot1() {
  TCanvas *c4 = new TCanvas("c4", "" , 800, 800);
  c4->cd();
  //c4->SetLogy();
  THStack *hs = new THStack("hs","");    
  //hs->Add(h_zjets_vbf); h_zjets_vbf->SetFillColor(kRed);
  //hs->Add(h_diboson); h_diboson->SetFillColor(kBlue);
  //hs->Add(h_top); h_top->SetFillColor(kGreen);
  //hs->Add(h_zjets); h_top->SetFillColor(kCyan);
  //hs->Add(h_wjets); h_wjets->SetFillColor(kViolet);
  hs->Add(h_zjets_vbf); h_zjets_vbf->SetFillColor(2);
  hs->Add(h_diboson); h_diboson->SetFillColor(3);
  hs->Add(h_top); h_top->SetFillColor(4);
  hs->Add(h_zjets); h_zjets->SetFillColor(5);
  hs->Add(h_wjets); h_wjets->SetFillColor(6);
  hs->Add(h_dijet); h_dijet->SetFillColor(7);
  //hs->Draw();
  //hs->GetXaxis()->SetTitle("jj_m[GeV]");
  //hs->GetYaxis()->SetTitle("evetns");
  //h_data->Draw("same e0");
  //hs->SetMaximum(1e7);
  //hs->SetMinimum(0.1);
  
  TLegend *leg=new TLegend(0.60,0.7,0.9,0.9);
  leg->AddEntry(h_data,"Data");
  leg->AddEntry(h_dijet,"Dijet");
  leg->AddEntry(h_wjets,"Wjets");
  leg->AddEntry(h_zjets,"Zjets");
  leg->AddEntry(h_top,"Top quark");
  leg->AddEntry(h_diboson,"Diboson");
  leg->AddEntry(h_zjets_vbf,"VBF Zjets");

  auto rp = new TRatioPlot(hs, h_data , "pois");
  //auto rp = new TRatioPlot(h_data, hs, "pois");

  c4->SetTicks(0, 1);
  rp->Draw();
  rp->GetLowerRefYaxis()->SetTitle("Ratio");
  rp->GetUpperRefYaxis()->SetTitle("Events / bin");
  rp->GetUpperRefXaxis()->SetTitle("jj_{m}[GeV]");
  leg->Draw("same");
  c4->Update();
  c4->SaveAs(dir_out+"/ratio.pdf");

  TCanvas *c5 = new TCanvas("c5", "" , 800, 800);
  c5->cd();
  c5->SetLogy();
  THStack *hs_mupt = new THStack("hs_mupt","");    
  //hs_mupt->Add(h_zjets_vbf); h_zjets_vbf->SetFillColor(kRed);
  //hs_mupt->Add(h_diboson); h_diboson->SetFillColor(kBlue);
  //hs_mupt->Add(h_top); h_top->SetFillColor(kGreen);
  //hs_mupt->Add(h_zjets); h_top->SetFillColor(kCyan);
  //hs_mupt->Add(h_wjets); h_wjets->SetFillColor(kViolet);
  hs_mupt->Add(h_mupt_zjets_vbf);  h_mupt_zjets_vbf->SetFillColor(2);
  hs_mupt->Add(h_mupt_diboson);    h_mupt_diboson->SetFillColor(3);
  hs_mupt->Add(h_mupt_top);        h_mupt_top->SetFillColor(4);
  hs_mupt->Add(h_mupt_zjets);      h_mupt_zjets->SetFillColor(5);
  hs_mupt->Add(h_mupt_wjets);      h_mupt_wjets->SetFillColor(6);
  hs_mupt->Add(h_mupt_dijet);      h_mupt_dijet->SetFillColor(7);
  //hs_mupt->Draw();
  //hs_mupt->GetXaxis()->SetTitle("jj_m[GeV]");
  //hs_mupt->GetYaxis()->SetTitle("evetns");
  //h_data->Draw("same e0");
  hs_mupt->SetMaximum(1e7);
  hs_mupt->SetMinimum(0.1);
  
  TLegend *leg_mupt=new TLegend(0.60,0.7,0.9,0.9);
  leg_mupt->AddEntry(h_mupt_data,"Data");
  leg_mupt->AddEntry(h_mupt_dijet,"Dijet");
  leg_mupt->AddEntry(h_mupt_wjets,"Wjets");
  leg_mupt->AddEntry(h_mupt_zjets,"Zjets");
  leg_mupt->AddEntry(h_mupt_top,"Top quark");
  leg_mupt->AddEntry(h_mupt_diboson,"Diboson");
  leg_mupt->AddEntry(h_mupt_zjets_vbf,"VBF Zjets");

  auto rp_mupt = new TRatioPlot(hs_mupt, h_mupt_data , "pois");
  //auto rp_mupt = new TRatioPlot(h_data, hs_mupt, "pois");

  c5->SetTicks(0, 1);
  rp_mupt->Draw();
  rp_mupt->GetLowerRefYaxis()->SetTitle("Ratio");
  rp_mupt->GetUpperRefYaxis()->SetTitle("Events / bin");
  rp_mupt->GetUpperRefXaxis()->SetTitle("mu_{pt}[GeV]");
  leg_mupt->Draw("same");
  c5->Update();
  c5->SaveAs(dir_out+"/ratio_mupt.pdf");


  TCanvas *c6 = new TCanvas("c6", "" , 800, 800);
  c6->cd();
  c6->SetLogy();
  THStack *hs_mupt_loweta = new THStack("hs_mupt_loweta","");    
  //hs_mupt_loweta->Add(h_zjets_vbf); h_zjets_vbf->SetFillColor(kRed);
  //hs_mupt_loweta->Add(h_diboson); h_diboson->SetFillColor(kBlue);
  //hs_mupt_loweta->Add(h_top); h_top->SetFillColor(kGreen);
  //hs_mupt_loweta->Add(h_zjets); h_top->SetFillColor(kCyan);
  //hs_mupt_loweta->Add(h_wjets); h_wjets->SetFillColor(kViolet);
  hs_mupt_loweta->Add(h_mupt_loweta_zjets_vbf);  h_mupt_loweta_zjets_vbf->SetFillColor(2);
  hs_mupt_loweta->Add(h_mupt_loweta_diboson);    h_mupt_loweta_diboson->SetFillColor(3);
  hs_mupt_loweta->Add(h_mupt_loweta_top);        h_mupt_loweta_top->SetFillColor(4);
  hs_mupt_loweta->Add(h_mupt_loweta_zjets);      h_mupt_loweta_zjets->SetFillColor(5);
  hs_mupt_loweta->Add(h_mupt_loweta_wjets);      h_mupt_loweta_wjets->SetFillColor(6);
  hs_mupt_loweta->Add(h_mupt_loweta_dijet);      h_mupt_loweta_dijet->SetFillColor(7);
  //hs_mupt_loweta->Draw();
  //hs_mupt_loweta->GetXaxis()->SetTitle("jj_m[GeV]");
  //hs_mupt_loweta->GetYaxis()->SetTitle("evetns");
  //h_data->Draw("same e0");
  hs_mupt_loweta->SetMaximum(1e7);
  hs_mupt_loweta->SetMinimum(0.1);
  
  TLegend *leg_mupt_loweta=new TLegend(0.60,0.7,0.9,0.9);
  leg_mupt_loweta->AddEntry(h_mupt_loweta_data,"Data");
  leg_mupt_loweta->AddEntry(h_mupt_loweta_dijet,"Dijet");
  leg_mupt_loweta->AddEntry(h_mupt_loweta_wjets,"Wjets");
  leg_mupt_loweta->AddEntry(h_mupt_loweta_zjets,"Zjets");
  leg_mupt_loweta->AddEntry(h_mupt_loweta_top,"Top quark");
  leg_mupt_loweta->AddEntry(h_mupt_loweta_diboson,"Diboson");
  leg_mupt_loweta->AddEntry(h_mupt_loweta_zjets_vbf,"VBF Zjets");

  auto rp_mupt_loweta = new TRatioPlot(hs_mupt_loweta, h_mupt_loweta_data , "pois");
  //auto rp_mupt_loweta = new TRatioPlot(h_data, hs_mupt_loweta, "pois");

  c6->SetTicks(0, 1);
  rp_mupt_loweta->Draw();
  rp_mupt_loweta->GetLowerRefYaxis()->SetTitle("Ratio");
  rp_mupt_loweta->GetUpperRefYaxis()->SetTitle("Events / bin");
  rp_mupt_loweta->GetUpperRefXaxis()->SetTitle("mu_{pt}[GeV]");
  leg_mupt_loweta->Draw("same");
  c6->Update();
  c6->SaveAs(dir_out+"/ratio_mupt_loweta_loweta.pdf");
}
void singlephoton(){
  vector<TString> dir_in;
  //dir_in.push_back("/lustre/collider/wangxi/data_dijet/common_ntupes_v5/mc16a/"); // end with /
  //dir_in.push_back("/lustre/collider/wangxi/data_dijet/common_ntupes_v5/mc16d/"); // end with /
  //dir_in.push_back("/lustre/collider/wangxi/data_dijet/common_ntupes_v5/mc16e/"); // end with /
  dir_in.push_back("/lustre/collider/wangxi/data_dijet/common_ntuple_v6/sp_convert/");
  //ihep // dir_in.push_back("/publicfs/atlas/atlasnew/SM/VBS/xiw/root_files/vbswy_tree5/convert/mc16a/"); // end with /
  //ihep // dir_in.push_back("/publicfs/atlas/atlasnew/SM/VBS/xiw/root_files/vbswy_tree5/convert/mc16d/"); // end with /
  //ihep // dir_in.push_back("/publicfs/atlas/atlasnew/SM/VBS/xiw/root_files/vbswy_tree5/convert/mc16e/"); // end with /
  //ihep dir_in.push_back("/publicfs/atlas/atlasnew/SM/VBS/zhenw/result/v05/ntuple/"); // end with /
        TH1F *mc_weight1 = new TH1F("mc_weight1","",50,-600,600);
        TH1F *mc_weight2 = new TH1F("mc_weight2","",50,-300,300);
   
  for (int di = 0; di< dir_in.size(); di++){
    if (auto dir = opendir(dir_in[di])) {
      while (auto f = readdir(dir)) {
        if (!f->d_name || f->d_name[0] == '.')
          continue; // Skip everything that starts with a dot

        TString f_name = f->d_name;
        printf("File: %s\n", f->d_name);

        std::string fileName_str = f->d_name;

        Double_t sf = ScaleFactor(dir_in[di]+f_name, false);

        int id;
        TTree *fChain_sp;
        if(SplitString(f_name, '.', id)){
          if(id >=364541 && id <= 364547){
            TFile *large_file = TFile::Open(dir_in[di]+f_name,"READ"); //TFile* myfile = new TFile("","OPEN");
            fChain_sp=(TTree*)large_file->Get("vbswy");
          }
          else
            continue;
        }
        else
          continue;
        
        vbswy vbswy(fChain_sp); // large_dijet.Loop();


        Long64_t nentries_sp = fChain_sp->GetEntriesFast();
        cout << "vbswy tree nentries: " << nentries_sp << endl;

        for (Long64_t jentry=0; jentry<nentries_sp;jentry++) {
          fChain_sp->GetEntry(jentry);
          Double_t mc_weight = sf * vbswy.weight_all*vbswy.weight_m *vbswy.weight_e;
          if (vbswy.dsid==364542) mc_weight1->Fill(vbswy.weight_all*vbswy.weight_m);
          if (vbswy.dsid==364543) mc_weight2->Fill(vbswy.weight_all*vbswy.weight_m);
          if (!vbswy.pass_init) continue;
          if (!vbswy.pass_trigger) continue;
          if (!vbswy.pass_dq) continue;
          if (!vbswy.pass_primary_vertex) continue;
          if (!(vbswy.pass_twojet && vbswy.pass_jetclean && vbswy.j0_pt>30 && vbswy.j1_pt>30)) continue;
          if (!vbswy.pass_vy_OR) continue;
          if (!vbswy.pass_duplicate) continue;
          if (!(vbswy.pass_onegam && vbswy.gam_idTight && vbswy.gam_isoTight)) continue;
          if (!(vbswy.pass_onelep && vbswy.pass_second_lepton_veto)) continue;
          if (!(vbswy.n_bjets_85 == 0)) continue;
          if (!(vbswy.is_Wmunu)) continue;
          if (!(vbswy.pass_dr)) continue;
          //if (!(vbswy.pass_ly_Z_veto && vbswy.pass_dphi)) continue;
          //if (!(vbswy.met>30 && vbswy.w_mt>30)) continue;
          
          //int baseline_mu = vbswy.gam_pt >22 && fabs(vbswy.gam_eta)<2.37 && vbswy.mu_pt>20 && fabs(vbswy.mu_eta)<2.5  && vbswy.met>30 && vbswy.w_mt>30  && vbswy.mu_passIP;
          //int baseline_mu = vbswy.gam_pt >22 && fabs(vbswy.gam_eta)<2.37 && vbswy.mu_pt>20 && fabs(vbswy.mu_eta)<2.5 && vbswy.mu_passIP && vbswy.mu_idLoose && vbswy.mu_isoFCLoose;
          int baseline_mu = vbswy.gam_pt >22 && fabs(vbswy.gam_eta)<2.37 && vbswy.mu_pt>20 && fabs(vbswy.mu_eta)<2.5 && vbswy.mu_passIP;
          if (!baseline_mu) continue;
          sp_event_yield->Fill(mc_weight);

          if (!(fabs(vbswy.weight_all*vbswy.weight_m *vbswy.weight_e)<100)) continue;
          sp2_event_yield->Fill(mc_weight);
          if (fabs(mc_weight)>300)
            continue;
          /*
          if (fabs(mc_weight)>100){
              cout<< "mc_weight: " << mc_weight << ", mu_pt: " << vbswy.mu_pt << ", dsid: " << vbswy.dsid << ", mu isoTight: "<< vbswy.mu_isoFCTight<< ", mu_idTight: " << vbswy.mu_idTight << ", mu_passIP: " << vbswy.mu_passIP << endl;
          } 
          */
          gamjets_origin->Fill(vbswy.mu_truth_origin, mc_weight);
	  if (vbswy.mu_isoFCTight && vbswy.mu_passIP && vbswy.mu_idTight){
          //if (vbswy.mu_isoFCTight && fabs(vbswy.mu_d0sig)<3 && vbswy.pass_dphi){
            //if(fabs(mc_weight)>20) continue;
            sp_event_yield_tight->Fill(mc_weight);
          /*  
            if (abs(mc_weight)>20 && vbswy.mu_pt<60){
              cout<< "mc_weight: " << mc_weight << ", mu_pt: " << vbswy.mu_pt << ", dsid: " << vbswy.dsid << ", mu isoTight: "<< vbswy.mu_isoFCTight<< ", mu_idTight: " << vbswy.mu_idTight << ", mu_passIP: " << vbswy.mu_passIP << endl;
              if(mc_weight<-20) continue;
            }
          */  
            h_sp_pt_tight->Fill(vbswy.mu_pt, mc_weight);
            h_sp_eta_tight->Fill(fabs(vbswy.mu_eta), mc_weight);
            h2_sp_pt_tight->Fill(vbswy.mu_pt, mc_weight);
            h2_sp_eta_tight->Fill(fabs(vbswy.mu_eta), mc_weight);
            t_h2_sp->Fill(vbswy.mu_pt, fabs(vbswy.mu_eta), mc_weight);
          } 
          //else if (!vbswy.mu_isoFCTight || (fabs(vbswy.mu_d0sig)>3 && fabs(vbswy.mu_d0sig)<10)){
          else{
            sp_event_yield_antiTight->Fill(mc_weight);
            h_sp_pt_antiTight->Fill(vbswy.mu_pt, mc_weight);
            h_sp_eta_antiTight->Fill(fabs(vbswy.mu_eta), mc_weight);
            h2_sp_pt_antiTight->Fill(vbswy.mu_pt, mc_weight);
            h2_sp_eta_antiTight->Fill(fabs(vbswy.mu_eta), mc_weight);
            anti_h2_sp->Fill(vbswy.mu_pt, fabs(vbswy.mu_eta), mc_weight);
          }
        } 
      }
    }
  }
  for (int i=1;i<5;i++){
    cout<< "gamjets mc bin anti-tight "   << i << ":	"  <<  h_sp_pt_antiTight->GetBinContent(i)     << endl;
    cout<< "gamjets mc bin tight "   << i << ":	"  <<  h_sp_pt_tight->GetBinContent(i)     << endl;
  }
  
  TCanvas *c7 = new TCanvas("c7", "", 800, 600);
  c7->cd();
  h2_sp_pt_antiTight->Draw("e0");
  h2_sp_pt_antiTight->GetYaxis()->SetRangeUser(-200,3000);
  h2_sp_pt_antiTight->SetLineColor(4);
  h2_sp_pt_tight->Draw("e0 same");
  h2_sp_pt_tight->SetLineColor(2);
  TLegend *leg_d=new TLegend(0.60,0.8,0.9,0.9); 
  leg_d->AddEntry(h2_sp_pt_antiTight, "antiTight");
  leg_d->AddEntry(h2_sp_pt_tight, "tight");
  leg_d->Draw("same");
  c7->SaveAs(dir_out+"/gamjet_distribution.png");

  TCanvas *c8 = new TCanvas("c8", "", 800, 600);
  c8->cd();
  h_sp_pt_antiTight->Draw("e0");
  h_sp_pt_antiTight->SetLineColor(2);
  h_sp_pt_tight->Draw("e0 same");
  h_sp_pt_tight->SetLineColor(4);
  TLegend *leg_d22=new TLegend(0.60,0.8,0.9,0.9); 
  leg_d22->AddEntry(h_sp_pt_antiTight, "antiTight");
  leg_d22->AddEntry(h_sp_pt_tight, "tight");
  leg_d22->Draw("same");
  c8->SaveAs(dir_out+"/gamjet_distribution_binning.png");

  TCanvas *c9 = new TCanvas("c9", "", 800, 600);
  c9->cd();
  mc_weight1->Draw();
  c9->SaveAs(dir_out+"/check_364542.png");

  TCanvas *c92 = new TCanvas("c92", "", 800, 600);
  c92->cd();
  mc_weight2->Draw();
  c92->SaveAs(dir_out+"/check_364543.png");

  TCanvas *c102 = new TCanvas("c102", "", 800, 600);
  c102->cd();
  c102->SetLogy();
  sp_event_yield->Draw();
  sp_event_yield->SetLineColor(4);
  sp2_event_yield->Draw("same");
  sp2_event_yield->SetLineColor(2);
  sp2_event_yield->SetLineWidth(2);
  c102->SaveAs(dir_out+"/sp_event_yield.png");

  TCanvas *c103 = new TCanvas("c103","", 800,600);
  c103->cd();
  c103->SetLogy();
  sp_event_yield_tight->Draw();
  c103->SaveAs(dir_out+"/sp_event_yield_tight.png");

  TCanvas *c104 = new TCanvas("c104","", 800,600);
  c104->cd();
  c104->SetLogy();
  sp_event_yield_antiTight->Draw();
  c104->SaveAs(dir_out+"/sp_event_yield_antiTight.png");
   
 /*  
  c2->cd();
  h_sp_pt_tight->Divide(h_sp_pt_antiTight);
  h_sp_pt_tight->Draw("e0 same");
  h_sp_pt_tight->SetLineColor(6);
  leg_ff->AddEntry(h_sp_pt_tight,"Gam+jets MC");
  leg_ff->Draw("same");
  c2->SaveAs(dir_out+"/mu_pt_ff_sp.pdf"); 
  for(int i=1; i<5;i++){
    cout<< "dijet mc bin "     << i << ":	"  <<  h_dijet_pt_tight->GetBinContent(i)  << endl;
    cout<< "gamjets mc bin "   << i << ":	"  <<  h_sp_pt_tight->GetBinContent(i)     << endl;
  }
  */
}

void plot(){
  c2->cd();
  c2->SetLogx();
  h_pt_tight->Divide(h_pt_antiTight);
  h_dijet_pt_tight->Divide(h_dijet_pt_antiTight);
  h_pt_tight->Draw("e0");
  h_pt_tight->SetLineColor(2);
  h_dijet_pt_tight->Draw("e0 same");
  h_dijet_pt_tight->SetLineColor(4);
  h_pt_tight->GetXaxis()->SetTitle("muon p_{T}[GeV]");
  h_pt_tight->GetYaxis()->SetRangeUser(-1,3);
  h_pt_tight->GetYaxis()->SetTitle("fake factor");
  leg_ff->AddEntry(h_pt_tight,"Dijet DD");
  leg_ff->AddEntry(h_dijet_pt_tight,"Dijet MC");

  h_sp_pt_tight->Divide(h_sp_pt_antiTight);
  h_sp_pt_tight->Draw("e0 same");
  h_sp_pt_tight->SetLineColor(6);
  leg_ff->AddEntry(h_sp_pt_tight,"Gam+jets MC");
  leg_ff->Draw("same");
  c2->SaveAs(dir_out+"/mu_pt_ff_sp.pdf");


  double bQuark_dj=0.,cQuark_dj=0.,other_dj=0.;
  double bQuark_yj=0.,cQuark_yj=0.,other_yj=0.;
  for(int i=2;i<=dijet_origin->GetNbinsX();++i)
  {
          if(i==27||i==34)
          {
                  bQuark_dj+=dijet_origin->GetBinContent(i);
                  bQuark_yj+=gamjets_origin->GetBinContent(i);
          }
          else if(i==26||i==33)
          {
                  cQuark_dj+=dijet_origin->GetBinContent(i);
                  cQuark_yj+=gamjets_origin->GetBinContent(i);
          }
          else
          {
                  other_dj+=dijet_origin->GetBinContent(i);
                  other_yj+=gamjets_origin->GetBinContent(i);
          }
  }
  double norm_dj=bQuark_dj+cQuark_dj+other_dj,norm_yj=bQuark_yj+cQuark_yj+other_yj;
  bQuark_dj=bQuark_dj/norm_dj;
  cQuark_dj=cQuark_dj/norm_dj;
  other_dj=other_dj/norm_dj;
 
  bQuark_yj=bQuark_yj/norm_yj;
  cQuark_yj=cQuark_yj/norm_yj;
  other_yj=other_yj/norm_yj;
  TH1F *hm=new TH1F();
  hm->Fill("GamJet",bQuark_yj);hm->Fill("DiJet",bQuark_dj);
  TH1F *hb=new TH1F();
  hb->Fill("GamJet",cQuark_yj);hb->Fill("DiJet",cQuark_dj);
  TH1F *ho=new TH1F();
  ho->Fill("GamJet",other_yj);ho->Fill("DiJet",other_dj);
  hm->SetFillColor(kAzure+1);
  hb->SetFillColor(kOrange+1);
  ho->SetFillColor(kBlue-5);
  THStack *horigin=new THStack("horigin","Jet Composition");
  horigin->Add(hm);
  horigin->Add(hb);
  horigin->Add(ho);
  TLegend *leg3=new TLegend(0.7,0.7,0.9,0.9);
  leg3->AddEntry(hm,"bQuark");
  leg3->AddEntry(hb,"cQuark");
  leg3->AddEntry(ho,"Other");
  TText *tyj[3],*tdj[3];
  double value_dj[3]={bQuark_dj,cQuark_dj,other_dj};
  double value_yj[3]={bQuark_yj,cQuark_yj,other_yj};
  for(int i=0;i<3;++i)
  {
          char cdj[100],cyj[100];
          double hdj=0,hyj=0;
          for(int j=0;j<i;++j)hdj+=value_dj[j];
          for(int j=0;j<i;++j)hyj+=value_yj[j];
          hdj+=value_dj[i]/2.;
          hyj+=value_yj[i]/2.;
          sprintf(cdj,"%.2e",value_dj[i]);
          sprintf(cyj,"%.2e",value_yj[i]);
          tyj[i]=new TText(hm->GetBinCenter(1),hyj,cyj);
          cout<<hm->GetBinCenter(1)<<" "<<hyj<<" "<<cyj<<endl;
          tdj[i]=new TText(hm->GetBinCenter(2),hdj,cdj);
          tyj[i]->SetTextSize(0.02);
          tdj[i]->SetTextSize(0.02);
  }
  TCanvas *c34 = new TCanvas("c34", "", 800, 600);
  c34->cd();
  gStyle->SetOptStat("");
  horigin->SetMaximum(1.5);
  horigin->Draw("HIST");
  horigin->GetYaxis()->SetTitle("Ratio");
  leg3->Draw();
  for(int i=0;i<3;++i)
  {
          tyj[i]->Draw();
          tdj[i]->Draw();
  }
  gPad->SaveAs(dir_out+"/Composition.pdf");
  //c3->SaveAs(dir_out+"/Composition.eps");



  for(int i=1; i<5;i++){
    cout<< "dijet data bin "   << i << ":	"  <<  h_pt_tight->GetBinContent(i) << " +-: "<< h_pt_tight->GetBinError(i)    << endl;
    cout<< "dijet mc bin "     << i << ":	"  <<  h_dijet_pt_tight->GetBinContent(i) << " +-: "<< h_dijet_pt_tight->GetBinError(i) << endl;
    cout<< "gamjets mc bin "   << i << ":	"  <<  h_sp_pt_tight->GetBinContent(i)  << " +-: "<< h_sp_pt_tight->GetBinError(i)     << endl;
  }
  for(int i=1; i<5;i++){
    cout<< "ratio bin " <<  i   << ":              " <<   fabs(h_sp_pt_tight->GetBinContent(i) - h_dijet_pt_tight->GetBinContent(i))/h_dijet_pt_tight->GetBinContent(i) << endl;
  }
}


int main(int argc,char *argv[]){
  cout<<"hello mlq"<<endl;

  TString *out = testGetOpt(argc, argv);
  cout << "wmt2 is " << out[0] << endl;
  wmt_modify = atoi(out[0]);
  cout << "met2 is " << out[1] << endl;
  met_modify = atoi(out[1]);
  cout << "jet2 is " << out[2] << endl;
  jet_modify = atoi(out[2]);
  cout << "dr2 is " << out[3] << endl;
  cout << "dphi2 is " << out[4] << endl;
  stringstream ss,cc;
  ss << out[3];
  cc << out[4];
  ss >> dr_modify;
  cc >> dphi_modify;
  cout << "dr3 is " << dr_modify << endl;
  cout << "dphi3 is " << dphi_modify << endl;
  njet_modify = atoi(out[5]);
  run_data = atoi(out[6]);
  cout << "rundata2 is " << run_data << endl;
  dir_out = out[7];

  if(system("mkdir -p "+dir_out)) return 0;  
  initialize();
  FakeFactor(met_modify, wmt_modify, jet_modify, dr_modify, dphi_modify , njet_modify, run_data);
  ratioplot1();
  singlephoton();
  plot(); 
  cout << "The run time is:" << (double)clock() /CLOCKS_PER_SEC<< "s" << endl;
}
