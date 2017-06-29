#include "TLatex.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TFile.h"
#include <iostream>
#include <TPad.h>
#include "TH1F.h"
#include "TString.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TMath.h"

#include <cmath>
#include <iomanip>
using namespace std;
double n_denom;

class TriggerThresholdDetermination {

    private:
        TFile * f_dijet;
        bool firstHist;
        TCanvas *c1, *c11, *c2;
        TGraph *graph,*graph2;

    public:
    TriggerThresholdDetermination();
    ~TriggerThresholdDetermination();

    void SetFile(TString fileName);

    void BuildRatio(TString nominator, TString denominator);

    void GetOfflineThreshold(double triggerThreshold);

};

TriggerThresholdDetermination::TriggerThresholdDetermination()
{
    firstHist = true;
    gStyle -> SetFrameBorderMode(0);

    gStyle -> SetNdivisions(10);
    gStyle -> SetCanvasBorderMode(0);
    gStyle -> SetPadBorderMode(1);
    gStyle -> SetOptTitle(1);
    gStyle -> SetStatFont(42);
    gStyle -> SetCanvasColor(10);
    gStyle -> SetPadColor(0);
    gStyle -> SetTitleFont(62,"xy");
    gStyle -> SetLabelFont(62,"xy");
    gStyle -> SetTitleFontSize(0.06);
    gStyle -> SetTitleSize(0.08,"xy");
    gStyle -> SetLabelSize(0.06,"xy");
    gStyle -> SetHistFillStyle(3001);
    gStyle -> SetHistFillColor(0);
    gStyle -> SetHistLineStyle(1);
    gStyle -> SetHistLineWidth(2);
    gStyle -> SetHistLineColor(2);
    gStyle -> SetOptStat(1110);
    gStyle -> SetOptStat(kFALSE);
    gStyle -> SetOptFit(0111);
    gStyle -> SetStatH(0.1);
    gStyle -> SetErrorX(0.001);
    gStyle -> SetLabelOffset(0.02,"xy");
    gStyle -> SetTitleOffset(1.0,"y");
    gStyle -> SetTitleOffset(0.99,"x");

    c1 = new TCanvas("c1","c1",0,0,800,600);
    gStyle->SetOptStat(0);
    c1->Divide(1,1,0,0);
    c1->SetFrameFillColor(0);

    c11 = new TCanvas("c11","c11",0,0,800,600);
    gStyle->SetOptStat(0);
    c11->Divide(1,1,0,0);
    c11->SetFrameFillColor(0);


    c2 = new TCanvas("c2","c2",0,0,800,600);
    c2->Divide(1,1);
    c2->SetFrameFillColor(0);

    graph = new TGraph();
    graph->SetMarkerStyle(20);

    graph2 = new TGraph();
    graph2->SetMarkerStyle(20);
    graph2->SetMarkerColor(4);
    float oldThreshold_run2[] = {60, 140, 260, 320}; 
    float Threshold_run2[] = {80, 135, 290, 365}; 


    graph2 = new TGraphAsymmErrors(4, oldThreshold_run2,Threshold_run2);
}

TriggerThresholdDetermination::~TriggerThresholdDetermination()
{
}

void TriggerThresholdDetermination::SetFile(TString fileName){
    f_dijet = new TFile (fileName, "read");
}

Double_t SmoothFit(Double_t *x, Double_t *par) {
   // cout << "n_denom = " << n_denom << endl;
    if (x[0] < n_denom) {
        TF1::RejectPoint();
        //  return 0;
    }
    Double_t p0 = par[0];
    Double_t p1 = par[1];
    Double_t N  = par[2];

    //Double_t fitval = N * (0.5 * (TMath::Erf(p0 * (x[0] - p1)) +1));
    Double_t fitval = 0.5 * N * (1. + TMath::Erf((x[0] - p0)/(pow(2,0.5)*p1)));

    return fitval;
}

void TriggerThresholdDetermination::BuildRatio(TString nominator, TString denominator)
{
    TH1F * nom   = (TH1F*) f_dijet -> Get(nominator);
    TH1F * denom = (TH1F*) f_dijet -> Get(denominator);

   // TH1F * eff_norm = new TH1F("eff_norm","eff_norm", 500, 0, 500);

    TH1F * eff = (TH1F*) nom -> Clone("eff");
    eff -> Divide(denom);

    double offlineThreshold = eff->GetBinLowEdge(eff->FindFirstBinAbove(0.99,1));

    TString sThreshold = nominator;
    int startIdx = sThreshold.Index("DiPFJetAve", 10, 0, TString::kExact);
    sThreshold  = sThreshold(startIdx+10,4);
    double triggerThreshold = sThreshold.Atof();
    std::cout << triggerThreshold << " -> " << offlineThreshold << std::endl;

    int oldSize = graph->GetN();


    TString sThreshold_fit = denominator;
    int startIdx_fit = sThreshold_fit.Index("DiPFJetAve", 10, 0, TString::kExact);
    sThreshold_fit  = sThreshold_fit(startIdx_fit+10,4);
    double triggerThreshold_fit = sThreshold_fit.Atof();
    n_denom = triggerThreshold_fit;

    TString sThreshold_tf1 = nominator;
    //TF1 * fEff_fi = sThreshold_tf1.Atof();

    TF1 * fEff_fit = new TF1(sThreshold_tf1, SmoothFit, n_denom, 500., 3);
   //Double_t Parameters[3] = {214.1630, 11.26654, 4.912097};
///    Double_t Parameters[3] = {84.1630, 7.26654, 27.912097};
   // Double_t Parameters[3] = {426.247, 11.26654, 0.912097};

    Double_t Parameters[3] = {426.247, 16.6424, 2.54782};






    fEff_fit -> SetParameters(Parameters);
    fEff_fit -> SetLineColor(4);
    fEff_fit -> SetParNames ("p0","p1","N");

    c1->cd();
    eff -> SetTitle(0);
    eff->GetYaxis()->SetTitle("Efficiency");
    eff->GetXaxis()->SetTitle("p_{T}^{ave} (GeV)");
    eff->GetXaxis()->SetTitleOffset(1.25);
    eff->GetYaxis()->SetTitleOffset(1.25);
    eff->SetMarkerStyle(20);
    eff->SetMarkerColor(kRed+oldSize-1);

    TLegend *leg = new TLegend(0.15,0.35+(oldSize*0.1),0.3,0.40+(oldSize*0.1));
    leg -> SetBorderSize(0);
    leg -> AddEntry(eff,sThreshold,"p");
    leg -> SetTextSize(0.035);
    leg -> SetFillColor(0);


    if(firstHist) {
        eff -> Draw("");
        eff -> Fit(fEff_fit,"R");
        leg -> Draw("same");
    }
    else {
        eff -> Draw("same");
        eff -> Fit(fEff_fit,"R","same");
        leg -> Draw("same");
    }
    double Norm = fEff_fit -> GetParameter(2);
    cout << Norm<<endl;
    TH1F * eff_n = (TH1F*) eff -> Clone("eff_n");
    eff_n -> Scale(1/Norm);
    double offlineThreshold_n = eff_n->GetBinLowEdge(eff_n->FindFirstBinAbove(0.985,1));
    std::cout <<"NORMALIZED: "<< triggerThreshold << " -> " << offlineThreshold_n << std::endl;
    graph->SetPoint(oldSize, triggerThreshold, offlineThreshold_n);

    c11->cd();
    if(firstHist) {
        eff_n->Draw("");
        eff_n -> Fit(fEff_fit,"R");
        leg -> Draw("same");
    }
    else {
        eff_n->Draw("same");
        eff_n -> Fit(fEff_fit,"R","same");
        leg -> Draw("same");
    }
    firstHist = false;

    c2->cd();
    graph->GetXaxis()->SetLabelOffset(0.01);
    graph->GetYaxis()->SetLabelOffset(0.01);
    graph->GetXaxis()->SetLabelSize(0.035);
    graph->GetYaxis()->SetLabelSize(0.035);
    graph->GetXaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetTitle("99% efficiency threshold (GeV)");
    graph->GetXaxis()->SetTitle("Nominal trigger threshold (GeV)");
    graph->Draw("AP");


    graph2->Draw("same P");
}

void TriggerThresholdDetermination::GetOfflineThreshold(double triggerThreshold)
{
    TF1* pol1 = new TF1("pol1","[0]+[1]*x",0,500);
    graph->Fit(pol1);


    std::cout << "extrapolated: " << 500 << " -> " << pol1->Eval(500) << std::endl;
}

void trigger_study_data()
{
    TriggerThresholdDetermination trig;
    trig.SetFile("/nfs/dust/cms/user/kovalch/sFrame/JEC/uhh2.AnalysisModuleRunner.DATA.data_8TeV_bacon_data_new_bins.root");
    //trig.SetFile("/nfs/dust/cms/user/kovalch/sFrame/JEC/uhh2.AnalysisModuleRunner.DATA.data_8TeV_bacon_test_new.root");

    //Selection  or noCuts
    trig.BuildRatio("noCuts/pt_ave_hltDiPFJetAve80",    "noCuts/pt_ave_hltDiPFJetAve40");
    trig.BuildRatio("noCuts/pt_ave_hltDiPFJetAve140",   "noCuts/pt_ave_hltDiPFJetAve80");
    trig.BuildRatio("noCuts/pt_ave_hltDiPFJetAve200",   "noCuts/pt_ave_hltDiPFJetAve140");
    trig.BuildRatio("noCuts/pt_ave_hltDiPFJetAve260",   "noCuts/pt_ave_hltDiPFJetAve200");
    trig.BuildRatio("noCuts/pt_ave_hltDiPFJetAve320",   "noCuts/pt_ave_hltDiPFJetAve260");
    trig.BuildRatio("noCuts/pt_ave_hltDiPFJetAve400",   "noCuts/pt_ave_hltDiPFJetAve320");

    trig.GetOfflineThreshold(40);
}