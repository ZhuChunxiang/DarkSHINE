#include <iostream>
#include <vector>
#include <string>
#include <TGraph.h>
#include <TH1F.h>
#include <TCanvas.h>

using namespace std;

string prefix="./time-consuming-image/";
vector<string> graphname;
vector<double> totaltime;
vector<double> cputime;
vector<double> transtime;
vector<double> kerneltime;
int Colorset[3] = {kBlack, kRed, kBlue};

void ManageConfigure(){
    graphname.push_back("photon1e6");    graphname.push_back("photon1e7");
    graphname.push_back("photon1e8");    graphname.push_back("photon1e9");
    graphname.push_back("totalcurve");

    // totaltime.push_back(0.4291106);   totaltime.push_back(1.9083494);   totaltime.push_back(16.8815507);   //photonnum=1E6: 1sim*100sets; 10sim*10sets; 100sim*1sets
    // totaltime.push_back(1.48885);   totaltime.push_back(2.39683);   totaltime.push_back(11.5043);   //photonnum=1E7: 1sim*100sets; 10sim*10sets; 100sim*1sets
    // totaltime.push_back(12.2612);   totaltime.push_back(13.4029);   totaltime.push_back(25.0729);   //photonnum=1E8: 1sim*100sets; 10sim*10sets; 100sim*1sets
    // totaltime.push_back(118.637635);   totaltime.push_back(120.343359);   totaltime.push_back(134.76841);   //photonnum=1E9: 1sim*100sets; 10sim*10sets; 100sim*1sets

    // transtime.push_back(0.1176615);   transtime.push_back(0.1156945);   transtime.push_back(0.1704366);   //photonnum=1E6: 1sim*100sets; 10sim*10sets; 100sim*1sets
    // transtime.push_back(0.1260228);   transtime.push_back(0.129776);   transtime.push_back(0.1723077);   //photonnum=1E7: 1sim*100sets; 10sim*10sets; 100sim*1sets
    // transtime.push_back(0.1721112);   transtime.push_back(0.1751934);   transtime.push_back(0.1909129);   //photonnum=1E8: 1sim*100sets; 10sim*10sets; 100sim*1sets
    // transtime.push_back(0.145021);   transtime.push_back(0.152324);   transtime.push_back(0.170392);   //photonnum=1E8: 1sim*100sets; 10sim*10sets; 100sim*1sets


    // kerneltime.push_back(0.3097437);   kerneltime.push_back(1.7899917);   kerneltime.push_back(16.6959744);   //photonnum=1E6: 1sim*100sets; 10sim*10sets; 100sim*1sets
    // kerneltime.push_back(1.3611513);   kerneltime.push_back(2.2643877);   kerneltime.push_back(11.3167448);   //photonnum=1E7: 1sim*100sets; 10sim*10sets; 100sim*1sets
    // kerneltime.push_back(12.0865692);   kerneltime.push_back(13.2241795);   kerneltime.push_back(24.8660006);   //photonnum=1E8: 1sim*100sets; 10sim*10sets; 100sim*1sets
    // kerneltime.push_back(118.490673);   kerneltime.push_back(120.187931);   kerneltime.push_back(134.582863);   //photonnum=1E8: 1sim*100sets; 10sim*10sets; 100sim*1sets

    totaltime.push_back( 1.960442);   totaltime.push_back(2.869855);   totaltime.push_back(18.071138);   //photonnum=1E6: 1sim*100sets; 10sim*10sets; 100sim*1sets
    totaltime.push_back(2.440206);   totaltime.push_back(3.352606);   totaltime.push_back(12.989590);   //photonnum=1E7: 1sim*100sets; 10sim*10sets; 100sim*1sets
    totaltime.push_back(13.004854);   totaltime.push_back(14.152341);   totaltime.push_back(25.919428);   //photonnum=1E8: 1sim*100sets; 10sim*10sets; 100sim*1sets
    totaltime.push_back(119.352382);   totaltime.push_back(119.966427);   totaltime.push_back(133.907592);   //photonnum=1E9: 1sim*100sets; 10sim*10sets; 100sim*1sets

    transtime.push_back( 1.609312);   transtime.push_back( 1.056949);   transtime.push_back(1.633788);   //photonnum=1E6: 1sim*100sets; 10sim*10sets; 100sim*1sets
    transtime.push_back(1.051096);   transtime.push_back(1.059043);   transtime.push_back(1.652039);   //photonnum=1E7: 1sim*100sets; 10sim*10sets; 100sim*1sets
    transtime.push_back(1.053089);   transtime.push_back(1.070379);   transtime.push_back( 1.642617);   //photonnum=1E8: 1sim*100sets; 10sim*10sets; 100sim*1sets
    transtime.push_back(1.624635);   transtime.push_back(1.609734);   transtime.push_back( 1.639190);   //photonnum=1E8: 1sim*100sets; 10sim*10sets; 100sim*1sets


    kerneltime.push_back(0.348224);   kerneltime.push_back( 1.804054);   kerneltime.push_back( 16.343340);   //photonnum=1E6: 1sim*100sets; 10sim*10sets; 100sim*1sets
    kerneltime.push_back(1.386277);   kerneltime.push_back(2.281943);   kerneltime.push_back(11.225342);   //photonnum=1E7: 1sim*100sets; 10sim*10sets; 100sim*1sets
    kerneltime.push_back(11.948542);   kerneltime.push_back(13.068029);   kerneltime.push_back(24.148891);   //photonnum=1E8: 1sim*100sets; 10sim*10sets; 100sim*1sets
    kerneltime.push_back(117.724327);   kerneltime.push_back( 118.343519);   kerneltime.push_back( 132.139942);   //photonnum=1E8: 1sim*100sets; 10sim*10sets; 100sim*1sets
   
    for(int i=0; i<totaltime.size(); i++){
        cputime.push_back( totaltime[i]-transtime[i]-kerneltime[i]);
    }


}


void DrawTimeConsuming(){
    gStyle->SetOptStat(0);
    ManageConfigure();
    vector<TH1F*> vhist;
    vhist.push_back( new TH1F("sim1set100","sim1set100",4,1,4) );   
    vhist.push_back( new TH1F("sim10set10","sim10set10",4,1,4) );
    vhist.push_back( new TH1F("sim100set1","sim100set1",4,1,4) );
    
    TLegend leg(0.1,0.7,0.4,0.9,"event configure");

    for(int i=0; i<vhist.size(); i++){
        vhist[i]->SetLineColor(Colorset[i]);
        vhist[i]->SetLineWidth(2);
        vhist[i]->Fill("CPU time",0);
        vhist[i]->Fill("transfer time",0);
        vhist[i]->Fill("Kernel time", 0);
        vhist[i]->Fill("Total time", 0);
        leg.AddEntry(vhist[i],vhist[i]->GetName());
    }



    for(int evtid=0; evtid<graphname.size()-1; evtid++){
        TCanvas* c=new TCanvas();
        for(int i=vhist.size()-1; i>=0;i--){
            vhist[i]->SetBinContent(1, cputime[evtid*3+i]);
            vhist[i]->SetBinContent(2, transtime[evtid*3+i]);
            vhist[i]->SetBinContent(3, kerneltime[evtid*3+i]);
            vhist[i]->SetBinContent(4, totaltime[evtid*3+i]);
            vhist[i]->Draw("hist l * same");
        }
        leg.Draw("same");
        c->SaveAs( (prefix+graphname[evtid]+".png").c_str() );
        c->Clear();
    }


    // TLegend leg2(0.1,0.7,0.4,0.9,"");
    // TCanvas* c=new TCanvas();
    // c->SetLogy();
    // TH1F* h=new TH1F("evttime","evttime",4,0,4);
    // h->Fill("1e6",totaltime[0]);
    // h->Fill("1e7",totaltime[3]);
    // h->Fill("1e8",totaltime[6]);
    // h->Fill("1e9",totaltime[9]);
    // h->SetLineWidth(2);     h->SetLineColor(kBlack);
    // h->Draw("hist l *");
    // leg2.AddEntry(h,"Minimal Time");

    // TH1F* assist_h=new TH1F("assevttime","assevttime",4,0,4);
    // // assist_h->Fill("1e6",totaltime[0]);
    // // assist_h->Fill("1e7",10*totaltime[0]);
    // // assist_h->Fill("1e8",100*totaltime[0]);
    // // assist_h->Fill("1e9",1000*totaltime[0]);
    // assist_h->Fill("1e6",totaltime[6]/100.);
    // assist_h->Fill("1e7",totaltime[6]/10.);
    // assist_h->Fill("1e8",1*totaltime[6]);
    // assist_h->Fill("1e9",10*totaltime[6]);
    // assist_h->SetLineWidth(2);     assist_h->SetLineColor(kRed);   assist_h->SetLineStyle(9);
    // assist_h->Draw("hist l same");
    // leg2.AddEntry(assist_h,"Assist_line");
    // leg2.Draw("same");
    // c->SaveAs( (prefix+"eventstime.png").c_str() );
}