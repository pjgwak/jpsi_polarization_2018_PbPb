#include <iostream>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include <TMath.h>
#include "../headers/commonUtility.h"
#include "../headers/cutsAndBin.h"
#include "../headers/HiEvtPlaneList.h"
#include "../headers/Style.h"
#include "../headers/tdrstyle.C"
#include "../headers/CMS_lumi_v2mass.C"
#include "../headers/rootFitHeaders.h"



static const long MAXTREESIZE = 1000000000000;
double getCorrection(TH2D *h, double pt, double angle);
double getAccWeight(TH1D *h = 0, double pt = 0);
double getEffWeight(TH1D *h = 0, double pt = 0);

void skim_to_roodata(
    string data_label_ = "250221",
    int nEvt = -1,
     int cLow = 0, int cHigh = 180,
    float massLow = 2.6, float massHigh = 3.5,
    // float massLow = 3.3, float massHigh =4.1, // 2S mass range
    bool dimusign = true, 
    bool isMC = false, int MCtype = 1,
    bool fAccW = true, bool fEffW = true,
    bool isTnP = true, bool isPtW = true,
    int hiHFBinEdge = 0)
{
    TStopwatch *t = new TStopwatch;
    t->Start();
    cout << "\n=================================\n";
    cout << "\n Convert FlowSkim to RooDataset\n";
    cout << "\n=================================\n";

    // ===== macro configure ===== //
    // namespaces
    using namespace std;
    using namespace RooFit;
    using namespace hi;

    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(0);
    gStyle->SetOptStat(000000000);
    gROOT->ForceStyle();
    // setTDRStyle();
    writeExtraText = true;

    // labeling
    TString DATE = data_label_;
    TString dimusignString;
    if (dimusign)
        dimusignString = "OS";
    else if (!dimusign)
        dimusignString = "SS";

    TString fCentSelHF = "Nom";
    if (hiHFBinEdge == 1)
        fCentSelHF = "Up";
    else if (hiHFBinEdge == -1)
        fCentSelHF = "Down";

    TString MCtype_tag = "";
    if (MCtype == 1)
        MCtype_tag = "PR";
    else if (MCtype == 2)
        MCtype_tag = "NP";

    // int iPeriod = 2; // Not used
    // int iPos = 33; // Nout used

    // TString kineLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, SiMuPtCut, cLow, cHigh, 1) ;
    // TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0, cLow, cHigh) ;

    // TH1::SetDefaultSumw2(); // Make histograms have errors
    // TH2::SetDefaultSumw2(); // 


    // ===== import FlowSkim files ===== //
    TChain *muon_chain = new TChain("mmepevt");
    if (!isMC)
    {
        // Use merged sample (DoubleMuon + Peripheral)
        TString f1 = "../files_skim/OniaFlowSkim_JpsiTrig_DBAll_AOD_isMC0_HFNom_250221.root";
        muon_chain->Add(f1.Data());
    }
    else if (MCtype == 1)
    {
        TString f1 = "../files_skim/OniaFlowSkim_Jpsi_AOD_isMC1_Signal_HFNom_250217.root";
        muon_chain->Add(f1.Data());
    }
    else if (MCtype == 2)
    {
        TString f1 = "../files_skim/OniaFlowSkim_Jpsi_AOD_isMC1_NPOnly_HFNom_250217.root";
        muon_chain->Add(f1.Data());
    }

    // ===== import Acc x Eff files ===== //
    // TFile *fin_weight = new TFile(Form("eff_acc_inputs/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_Jpsi_PtW%d_tnp%d_210913.root", isPtW, isTnP), "read");
    TFile *fin_correct = new TFile("../eff_acc/roots/mc_eff_vs_pt_cent_rap_prompt_pbpb_Jpsi_PtW1_tnp1_250221.root");

    TH2D *h_correct[6];
    // h_correct[0] = (TH1D *)fEff->Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_180_absy0_1p2", isTnP, isPtW));
    
    h_correct[0] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent0to20_absy0_1p6");
    h_correct[1] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent0to20_absy1p6_2p4");
    h_correct[2] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent20_60_absy0_1p6");
    h_correct[3] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent20_60_absy1p6_2p4");
    h_correct[4] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent60to180_absy0_1p6");
    h_correct[5] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent60to180_absy1p6_2p4");

    // ===== connect to input tree ===== //
    // SetBranchAddress
    Int_t event;
    Int_t cBin;
    Int_t nDimu;
    const int nMaxDimu = 1000;
    int recoQQsign[nMaxDimu];
    float vz;
    float mass[nMaxDimu];
    float pt[nMaxDimu];
    float y[nMaxDimu];
    float pt1[nMaxDimu];
    float pt2[nMaxDimu];
    float eta[nMaxDimu];
    float eta1[nMaxDimu];
    float eta2[nMaxDimu];
    float ctau3D[nMaxDimu];
    float ctau3DErr[nMaxDimu];
    float ctau3DTrue[nMaxDimu];
    double weight;

    float cos_theta[nMaxDimu];
    float cos_theta1[nMaxDimu];
    float phi[nMaxDimu];
    float phi1[nMaxDimu];
    float cos_cs[nMaxDimu];
    float phi_cs[nMaxDimu];

    float cos_hx[nMaxDimu];
    float phi_hx[nMaxDimu];
    float cos_ep[nMaxDimu];
    float phi_ep[nMaxDimu];

    TBranch *b_event;
    TBranch *b_cBin;
    TBranch *b_nDimu;
    TBranch *b_recoQQsign;
    TBranch *b_vz;
    TBranch *b_mass;
    TBranch *b_pt;
    TBranch *b_y;
    TBranch *b_eta;
    TBranch *b_eta1;
    TBranch *b_eta2;
    TBranch *b_ctau3D;
    TBranch *b_ctau3DErr;
    TBranch *b_ctau3DTrue;
    TBranch *b_pt1;
    TBranch *b_pt2;
    TBranch *b_weight;

    TBranch *b_cos_theta;
    TBranch *b_cos_theta1;
    TBranch *b_phi;
    TBranch *b_phi1;
    TBranch *b_cos_cs;
    TBranch *b_phi_cs;

    TBranch *b_cos_hx;
    TBranch *b_phi_hx;
    TBranch *b_cos_ep;
    TBranch *b_phi_ep;

    muon_chain->SetBranchAddress("event", &event, &b_event);
    muon_chain->SetBranchAddress("cBin", &cBin, &b_cBin);
    muon_chain->SetBranchAddress("nDimu", &nDimu, &b_nDimu);
    muon_chain->SetBranchAddress("recoQQsign", recoQQsign, &b_recoQQsign);
    muon_chain->SetBranchAddress("vz", &vz, &b_vz);
    muon_chain->SetBranchAddress("mass", mass, &b_mass);
    muon_chain->SetBranchAddress("y", y, &b_y);
    muon_chain->SetBranchAddress("pt", pt, &b_pt);
    muon_chain->SetBranchAddress("pt1", pt1, &b_pt1);
    muon_chain->SetBranchAddress("pt2", pt2, &b_pt2);
    muon_chain->SetBranchAddress("eta", eta, &b_eta);
    muon_chain->SetBranchAddress("eta1", eta1, &b_eta1);
    muon_chain->SetBranchAddress("eta2", eta2, &b_eta2);
    muon_chain->SetBranchAddress("ctau3D", ctau3D, &b_ctau3D);
    muon_chain->SetBranchAddress("ctau3DTrue", ctau3DTrue, &b_ctau3DTrue);
    muon_chain->SetBranchAddress("ctau3DErr", ctau3DErr, &b_ctau3DErr);
    // muon_chain->SetBranchAddress("cos_theta_hx", cos_theta_hx, &b_cos_theta_hx);
    muon_chain->SetBranchAddress("weight", &weight, &b_weight);
    muon_chain->SetBranchAddress("cos_theta", cos_theta, &b_cos_theta);
    muon_chain->SetBranchAddress("cos_theta1", cos_theta1, &b_cos_theta1);
    muon_chain->SetBranchAddress("phi", phi, &b_phi);
    muon_chain->SetBranchAddress("phi1", phi1, &b_phi1);
    muon_chain->SetBranchAddress("cos_cs", cos_cs, &b_cos_cs);
    muon_chain->SetBranchAddress("phi_cs", phi_cs, &b_phi_cs);
    muon_chain->SetBranchAddress("cos_hx", cos_hx, &b_cos_hx);
    muon_chain->SetBranchAddress("phi_hx", phi_hx, &b_phi_hx);
    muon_chain->SetBranchAddress("cos_ep", cos_ep, &b_cos_ep);
    muon_chain->SetBranchAddress("phi_ep", phi_ep, &b_phi_ep);


    // ===== Outfile ===== //
    TFile *out_file = nullptr;
    if (isMC)
    {
        out_file = new TFile(Form("../files_roodata/RooDataSet_miniAOD_isMC%d_%s_Jpsi_cent%i_%i_Effw%d_Accw%d_PtW%d_TnP%d_HF%s_%s.root",
                            isMC, MCtype_tag.Data(), cLow, cHigh, fEffW, fAccW, isPtW, isTnP, fCentSelHF.Data(), DATE.Data()),
                       "recreate");
    }
    else
    {
        out_file = new TFile(Form("../files_roodata/RooDataSet_miniAOD_isMC%d_Jpsi_cent%i_%i_Effw%d_Accw%d_PtW%d_TnP%d_HF%s_%s.root",
                            isMC, cLow, cHigh, fEffW, fAccW, isPtW, isTnP, fCentSelHF.Data(), DATE.Data()),
                       "recreate");
    }
    out_file->cd();


    // ==== build output RooDataset ===== //
    RooRealVar *massVar = new RooRealVar("mass", "mass variable", 1.0, 6.0, "GeV/c^{2}");
    RooRealVar *ptVar = new RooRealVar("pt", "pt variable", 0, 100, "GeV/c");
    RooRealVar *yVar = new RooRealVar("y", "rapidity of the dimuon pair", -5, 5, "");
    RooRealVar *pt1Var = new RooRealVar("pt1", "pt of muon+", 0, 500, "GeV/c");
    RooRealVar *eta1Var = new RooRealVar("eta1", "eta of muon+", -4, 4, "");
    RooRealVar *pt2Var = (RooRealVar *)pt1Var->Clone("pt2");
    RooRealVar *eta2Var = (RooRealVar *)eta1Var->Clone("eta2");
    RooRealVar *cBinVar = new RooRealVar("cBin", "Centrality bin", 0, 200, "");
    RooRealVar *ep2Var = new RooRealVar("ep2", "2nd order event plane", -100, 100, "");
    RooRealVar *evtWeight = new RooRealVar("weight", "corr weight", 0, 10000, "");
    // local name은 evtWeight인데 root 이름연 weight이라 output에 weight으로 저장된다.
    RooRealVar *recoQQ = new RooRealVar("recoQQsign", "qq sign", -1, 3, "");
    RooRealVar *ctau3DVar = new RooRealVar("ctau3D", "c_{#tau}", -10.0, 10.0, "mm");
    RooRealVar *ctau3DErrVar = new RooRealVar("ctau3DErr", "#sigma_{c#tau}", 0, 10.0, "mm");
    RooRealVar *ctau3DResVar = new RooRealVar("ctau3DRes", "c_{#tau} res", -100.0, 100.0, "");
    RooRealVar *ctau3DTrueVar = new RooRealVar("ctau3DTrue", "c_{#tau} true", -10.0, 10.0, "");
    RooRealVar *NumDimu = new RooRealVar("NumDimu", "number of dimuon", 0, 1000, "");
    
    // polarization variables
    RooRealVar *cos_thetaVar = new RooRealVar("cos_theta", "", -1.0, 1.0, "");
    RooRealVar *cos_theta1Var = new RooRealVar("cos_theta1", "", -1.0, 1.0, "");
    RooRealVar *phiVar = new RooRealVar("phi", "", -5.0, 5.0, "");
    RooRealVar *phi1Var = new RooRealVar("phi1", "", -5.0, 5.0, "");

    RooRealVar *cos_csVar = new RooRealVar("cos_cs", "", -1.0, 1.0, "");
    RooRealVar *phi_csVar = new RooRealVar("phi_cs", "", -5.0, 5.0, "");
    RooRealVar *cos_hxVar = new RooRealVar("cos_hx", "", -1.0, 1.0, "");
    RooRealVar *phi_hxVar = new RooRealVar("phi_hx", "", -5.0, 5.0, "");
    RooRealVar *cos_epVar = new RooRealVar("cos_ep", "", -1.0, 1.0, "");
    RooRealVar *phi_epVar = new RooRealVar("phi_ep", "", -5.0, 5.0, "");

    RooArgSet *argSet = new RooArgSet(
        *massVar, *ptVar, *yVar, *pt1Var, *pt2Var, *eta1Var, *eta2Var, *evtWeight,
        *cBinVar, *recoQQ, *NumDimu, *ctau3DVar, *ctau3DErrVar, *ctau3DResVar,
        *cos_thetaVar, *cos_theta1Var, *phiVar, *phi1Var, *cos_csVar, *phi_csVar,
        *cos_hxVar, *phi_hxVar, *cos_epVar, *phi_epVar);
    if (isMC)
        argSet->add(*ctau3DTrueVar);
    RooDataSet *ds = new RooDataSet("dataset", "", *argSet);

    
    // ===== weight and counters ===== //
    double weight_acc_eff = 1;
    double weight_eff = 1;

    int nDimu_all = 0;
    int nDimuPass = 0;
    int nDimu_one = 0;
    int nDimu_more = 0;
    int count_dimu = 0;

    // set a number of events
    if (nEvt == -1)
        nEvt = muon_chain->GetEntries();
    cout << "Total events = " << nEvt << " / " << muon_chain->GetEntries() << endl;


    // ===== begin main loop ===== //
    for (int i = 0; i < nEvt; i++)
    {
        muon_chain->GetEntry(i);
        
        if (fabs(vz) > 15)
            continue;
        
        nDimuPass = 0;

        if ( !(cBin >= cLow && cBin < cHigh))
            continue;
        

        // ===== start dimuon loop ===== //
        for (int j = 0; j < nDimu; j++)
        {
            // check cut once more with mass cut
            if ( !((double)pt[j] < 50 && recoQQsign[j] == 0 && mass[j] > massLow && mass[j] < massHigh && abs(y[j]) < 2.4 && IsAcceptanceQQ(pt1[j], eta1[j]) && IsAcceptanceQQ(pt2[j], eta2[j])))
                continue;

            // ===== apply acc, eff weight ===== //
            weight_acc_eff = 1;
            // weight_eff = 1;

            if (fAccW || fEffW) {
                if(cBin<20) {
                    if (abs((double)y[j]) < 1.6)
                        weight_acc_eff = getCorrection(h_correct[0], pt[j], cos_ep[j]);
                    else if (abs((double)y[j]) > 1.6 && abs((double)y[j]) < 2.4)
                        weight_acc_eff = getCorrection(h_correct[1], pt[j], cos_ep[j]);
                } else if (cBin > 20 && cBin < 60) {
                    if (abs((double)y[j]) < 1.6)
                        weight_acc_eff = getCorrection(h_correct[2], pt[j], cos_ep[j]);
                    else if (abs((double)y[j]) > 1.6 && abs((double)y[j]) < 2.4)
                        weight_acc_eff = getCorrection(h_correct[3], pt[j], cos_ep[j]);
                } else if (cBin > 60 && cBin < 180) {
                    if (abs((double)y[j]) < 1.6)
                        weight_acc_eff = getCorrection(h_correct[4], pt[j], cos_ep[j]);
                    else if (abs((double)y[j]) > 1.6 && abs((double)y[j]) < 2.4)
                        weight_acc_eff = getCorrection(h_correct[5], pt[j], cos_ep[j]);
                }
            }

            double weight_ = weight * weight_acc_eff;

            recoQQ->setVal((int)recoQQsign[j]);
            massVar->setVal((double)mass[j]);
            ptVar->setVal((double)pt[j]);
            yVar->setVal((double)y[j]);
            pt1Var->setVal((double)pt1[j]);
            eta1Var->setVal((double)eta1[j]);
            pt2Var->setVal((double)pt2[j]);
            eta2Var->setVal((double)eta2[j]);
            cBinVar->setVal((double)cBin);
            ctau3DVar->setVal((double)ctau3D[j]);
            ctau3DErrVar->setVal((double)ctau3DErr[j]);
            ctau3DResVar->setVal((double)ctau3D[j] / ctau3DErr[j]);
            if (isMC)
                ctau3DTrueVar->setVal((double)ctau3DTrue[j]);
            // cos_theta_hxVar->setVal((double)cos_theta_hx[j]);
            evtWeight->setVal((double)weight_);
            NumDimu->setVal((int)nDimu);
            // cout<<"Evt: "<<j<<", Cent: "<<cBin<<", mass: "<<mass[j]<<", pt: "<<pt[j]<<", pt1:"<<pt1[j]<<", pt2: "<<pt2[j]<<", eta1: "<<eta1[j]<<", eta2: "<<eta2[j]<<", vz: "<<vz<<endl;
            cos_thetaVar->setVal((float)cos_theta[j]);
            cos_theta1Var->setVal((float)cos_theta1[j]);
            phiVar->setVal((float)phi[j]);
            phi1Var->setVal((float)phi1[j]);
            
            cos_csVar->setVal((float)cos_cs[j]);
            phi_csVar->setVal((float)phi_cs[j]);

            cos_hxVar->setVal((float)cos_hx[j]);
            phi_hxVar->setVal((float)phi_hx[j]);

            cos_epVar->setVal((float)cos_ep[j]);
            phi_epVar->setVal((float)phi_ep[j]);

            ds->add(*argSet);

            count_dimu++;
        } // end of dimuon loop
    }
    cout << "How many Jpsi??: " << count_dimu << endl;


    // ==== save output ===== //
    ds->Write();
    out_file->Close();

    cout << "\n=================================\n";
    cout << "\n Finish converting\n";
    cout << "\n=================================\n";
    t->Stop();
    printf("RealTime=%f seconds, CpuTime=%f seconds\n", t->RealTime(), t->CpuTime());
}

// ===== Other functions ===== //
double getCorrection(TH2D *h, double pt, double angle)
{
    int bin_angle = h->GetXaxis()->FindBin(angle);
    int bin_pt = h->GetYaxis()->FindBin(pt);
    double weight_ = 1. / (h->GetBinContent(bin_angle, bin_pt));
    return weight_;
}

double getAccWeight(TH1D *h, double pt)
{
    // cout<<"pt: "<<pt<<endl;
    double binN = h->FindBin(pt);
    double weight_ = 1. / (h->GetBinContent(binN));
    return weight_;
}

double getEffWeight(TH1D *h, double pt)
{
    double binN = h->FindBin(pt);
    double weight_ = 1. / (h->GetBinContent(binN));
    return weight_;
}
// ===== Other functions ===== //