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
        TGraph *graph;

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
    gStyle -> SetOptFit(0);
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

    TF1 * fEff_fit = new TF1(sThreshold_tf1, SmoothFit, n_denom, 600., 3);
   //Double_t Parameters[3] = {214.1630, 11.26654, 4.912097};
///    Double_t Parameters[3] = {84.1630, 7.26654, 27.912097};
   // Double_t Parameters[3] = {426.247, 11.26654, 0.912097};

    ////Double_t Parameters[3] = {426.247, 16.6424, 2.54782};

    Double_t Parameters[3] = {triggerThreshold, 6.4, 2.8};





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
    eff->SetMarkerColor(kRed+oldSize-3);
//     eff->SetMarkerColor(n_denom+kRed);

    TLegend *leg = new TLegend(0.15,0.5+(oldSize*0.05),0.3,0.5+(oldSize*0.05));
    leg -> SetBorderSize(0);
    leg -> AddEntry(eff,"HLT_DiPFJetAve"+sThreshold+"_HFJEC","p");
    leg -> SetTextSize(0.035);
    leg -> SetFillColor(0);


    if(firstHist) {
        eff -> Draw("");
        eff -> Fit(fEff_fit,"RLI");
        leg -> Draw("same");
    }
    else {
        eff -> Draw("same");
        eff -> Fit(fEff_fit,"RLI","same");
        leg -> Draw("same");
    }
    double x0 = fEff_fit -> GetParameter(0);
     eff -> Fit(fEff_fit,"RLI","same",0.5*x0,3*x0);
    double Norm = fEff_fit -> GetParameter(2);
    cout << "Norm = "<< Norm<<endl;
    TH1F * eff_n = (TH1F*) eff -> Clone("eff_n");
    eff_n -> Scale(1/Norm);


    c11->cd();
    if(firstHist) {
        eff_n->Draw("");
        eff_n -> Fit(fEff_fit,"RLI","",0.5*x0,3*x0);
        leg -> Draw("same");
    }
    else {
        eff_n->Draw("same");
        eff_n -> Fit(fEff_fit,"RLI","same",0.5*x0,3*x0);
        leg -> Draw("same");
    }
    firstHist = false;
//     if (fEff_fit->FitResult() == 0.99) cout <<"fit = "<< p0 <<endl;
    TLatex *mylatex = new TLatex(40,7.9, "CMS, 2.11 fb^{-1} (13 TeV)");
    mylatex->SetTextSize(0.04);
    mylatex->Draw();
    double offlineThreshold_n_h = eff_n->GetBinLowEdge(eff_n->FindFirstBinAbove(0.985,1));
    double offlineThreshold_n = fEff_fit ->GetX(0.99, 0, 600);

    std::cout <<"NORMALIZED: "<< triggerThreshold << " -> " << offlineThreshold_n_h << std::endl;
    cout <<"Thresholds from the fit for: "<< triggerThreshold << " -> " <<offlineThreshold_n<<endl;

    graph->SetPoint(oldSize, triggerThreshold, offlineThreshold_n);


   // cout <<"eff_n = "<<fEff_fit ->GetX(0.99, 0, 600)<<endl;
/*   
    TF1 * fEff_fit = new TF1(sThreshold_tf1, SmoothFit, n_denom, 600., 3);

nu dyvysj tobi treba vytjahnuty funkciju (rezuljtat fita), TF1 *
a potim
function -> GetX(0.99, ptmin, ptmax)
de ptmin i ptmax zadajutj mezhi de shukajetjsa rozvjazok
napr. ptmin =0, ptmax=600*/


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
}

void TriggerThresholdDetermination::GetOfflineThreshold(double triggerThreshold)
{
    TF1* pol1 = new TF1("pol1","[0]+[1]*x",0,600);
    graph->Fit(pol1);
    std::cout << "extrapolated: " << 40 << " -> " << pol1->Eval(40) << std::endl;
}

void trigger_study_data_run2()
{
    TriggerThresholdDetermination trig;
///    trig.SetFile("/nfs/dust/cms/user/kovalch/sFrame/JEC/uhh2.AnalysisModuleRunner.DATA.all_data_pu_dist_only_for_251721_goldenjson_40pb.root");
    //trig.SetFile("/nfs/dust/cms/user/kovalch/sFrame/JEC/uhh2.AnalysisModuleRunner.DATA.data_8TeV_bacon_test_new.root");
// trig.SetFile(" /nfs/dust/cms/user/kovalch/sFrame/JEC/uhh2.AnalysisModuleRunner.DATA.data_no_L1data_L2L3_new_triggers_convert.root");
    trig.SetFile("/nfs/dust/cms/user/kovalch/sFrame/JEC/V6/uhh2.AnalysisModuleRunner.DATA.DATAdata_1200pt_ave_V6_CDv3Dv4.root");

    //Selection  or noCuts
//    trig.BuildRatio("noCuts/pt_ave_hltDiPFJetAve60",    "noCuts/pt_ave_hltDiPFJetAve40");
//      trig.BuildRatio("noCuts/pt_ave_hltDiPFJetAve80",    "noCuts/pt_ave_hltDiPFJetAve60");
//       trig.BuildRatio("noCuts/pt_ave_hltDiPFJetAve140",   "noCuts/pt_ave_hltDiPFJetAve80");
//     trig.BuildRatio("noCuts/pt_ave_hltDiPFJetAve200",   "noCuts/pt_ave_hltDiPFJetAve140");
//  trig.BuildRatio("noCuts/pt_ave_hltDiPFJetAve260",   "noCuts/pt_ave_hltDiPFJetAve200");
//     trig.BuildRatio("noCuts/pt_ave_hltDiPFJetAve320",   "noCuts/pt_ave_hltDiPFJetAve260");
//       trig.BuildRatio("noCuts/pt_ave_hltDiPFJetAve400",   "noCuts/pt_ave_hltDiPFJetAve320");
//       trig.BuildRatio("noCuts/pt_ave_hltDiPFJetAve500",   "noCuts/pt_ave_hltDiPFJetAve400");

    trig.BuildRatio("noCuts/HLT_DiPFJetAve80_HFJEC",   "noCuts/HLT_DiPFJetAve60_HFJEC");
    trig.BuildRatio("noCuts/HLT_DiPFJetAve100_HFJEC",   "noCuts/HLT_DiPFJetAve80_HFJEC");
   trig.BuildRatio("noCuts/HLT_DiPFJetAve160_HFJEC",   "noCuts/HLT_DiPFJetAve100_HFJEC");
    trig.BuildRatio("noCuts/HLT_DiPFJetAve220_HFJEC",   "noCuts/HLT_DiPFJetAve160_HFJEC");
    trig.BuildRatio("noCuts/HLT_DiPFJetAve300_HFJEC",   "noCuts/HLT_DiPFJetAve220_HFJEC");
    trig.GetOfflineThreshold(40);
}