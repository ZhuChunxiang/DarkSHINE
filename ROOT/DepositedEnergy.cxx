#include<iostream>
#include<fstream>
#include<TFile.h>
#include<TGraph.h>
#include<TROOT.h>
#include<TCanvas.h>
#include<TH1.h>

using namespace std;

void DepositedEnergy()
{
    string filename;
    const char * Name;
    double depEnergy;

    vector<int> InitialEnerge;
    vector<double> MeanEnergy;
    vector<double> RMSEnergy;

    ifstream datalist("/lustre/collider/liudanning/DarkPhoton/StandAlone/Test/ROOT/datalist.txt");

    TH1D * hist = new TH1D("hist","hist",1000,0,1000);

    for(int i = 0 ; !datalist.eof() ; i++)
    {
        if(datalist.eof()) break;


        datalist>>filename;
        Name=filename.c_str();
        cout<<Name<<endl;

        int pos1 = filename.find("photon");
        int pos2 = filename.find("MeV");

        int energy;
        string Temp;
        double MeanValue = 0;
        double RMSValue = 0;

        Temp = filename.substr(pos1+7,pos2-pos1-7);
        std::stringstream stream(Temp);
        stream>>energy;
        cout<<"i = "<<i<<" energy is "<<energy<<endl;

        TFile * file = TFile::Open(Name,"READ");
        if(!file) cout<<"Can not Read ROOT File"<<endl;
        
        TTree * tree = (TTree*)file->Get("TTree");
        tree->SetBranchAddress("depEnergy", &depEnergy);

        int Entries = tree->GetEntries();
        for (int j = 0 ; j < Entries ; j++)
        {
            tree->GetEntry(i);
            hist->Fill(depEnergy);
        }

        MeanValue = hist->GetMean();
        RMSValue = hist->GetStdDev();
        MeanEnergy.push_back(MeanValue);
        RMSEnergy.push_back(RMSValue);
        InitialEnerge.push_back(energy);
        cout<<"Mean Energy value is "<<MeanValue<<endl;
        cout<<"RMS Value is "<<RMSValue<<endl;

        hist->Reset();

    }

    int n = MeanEnergy.size();

    TGraph * graph = new TGraph();
    TGraphErrors * graError = new TGraphErrors();
    double Xvalue, Yvalue, XRMS, YRMS = 0;

    for(int i = 0 ; i < 10 ; i++)
    {
        Xvalue = InitialEnerge[i];
        Yvalue = MeanEnergy[i];
        XRMS = 0;
        YRMS = RMSEnergy[i];
        graph->SetPoint(i,Xvalue,Yvalue/Xvalue);
        graError->SetPoint(i, Xvalue, Yvalue/Xvalue);
        graError->SetPointError(i, XRMS, YRMS*10);
    }

    TCanvas * c1 = new TCanvas();
    TPad * thePad1 = new TPad(); 
    thePad1->Draw();
    thePad1->cd();
    thePad1->SetLogx();

    graph->SetMarkerStyle(21);
    graph->SetMarkerSize(2.);
    graph->GetXaxis()->SetTitle("Incident energy[MeV]");
    graph->GetYaxis()->SetTitle("Deposited energy[MeV]");
    graph->Draw("AP");

    TCanvas * c2 = new TCanvas();
    TPad * thePad2 = new TPad(); 
    thePad2->Draw();
    thePad2->cd();
    thePad2->SetLogx();

    graError->SetMarkerStyle(20);
    graError->SetMarkerSize(2.);
    graError->SetMarkerColor(kRed);
    graError->GetXaxis()->SetTitle("Incident energy[MeV]");
    graError->GetYaxis()->SetTitle("Deposited energy[MeV]");
    graError->Draw("AP");
}
