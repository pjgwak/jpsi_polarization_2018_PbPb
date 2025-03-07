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

using namespace std;
using namespace RooFit;
using namespace hi;

static const long MAXTREESIZE = 1000000000000;

double getAccWeight(TH1D *h = 0, double pt = 0);
double getEffWeight(TH1D *h = 0, double pt = 0);
void GetHistSqrt(TH1D *h1 = 0, TH1D *h2 = 0);
void Get2DHistSqrt(TH2D *h1 = 0, TH2D *h2 = 0);
double GetHistBkg(TH1D *h1 = 0, TH1D *h2 = 0);

void skim_to_roodata_gen_only(
    string data_label_ = "250307",
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
    cout << "\n Start conversion from FlowSkim to RooDataset\n";
    cout << "\n=================================\n";

    // Basic Setting
    gStyle->SetOptStat(0);

    TString DATE = data_label_;
    gStyle->SetEndErrorSize(0);
    gStyle->SetOptStat(000000000);
    gROOT->ForceStyle();
    setTDRStyle();
    writeExtraText = true;

    // int iPeriod = 2; // Not used
    // int iPos = 33; // Nout used

    // TString kineLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, SiMuPtCut, cLow, cHigh, 1) ;
    // TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0, cLow, cHigh) ;

    TH1::SetDefaultSumw2(); // Make hist have errors

    // labels
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

    // ===== import sample files ===== //
    TChain *tree = new TChain("mmgentree");
    if (!isMC)
    {
        // Use merged sample (DoubleMuon + Peripheral)
        TString f1 = "../files_skim/OniaFlowSkim_JpsiTrig_DBAll_miniAOD_isMC0_HFNom_250201.root";
        tree->Add(f1.Data());
    }
    else if (MCtype == 1)
    {
        TString f1 = "../files_skim/OniaFlowSkim_Jpsi_miniAOD_isMC1_Signal_HFNom_250211.root";
        tree->Add(f1.Data());
    }
    else if (MCtype == 2)
    {
        TString f1 = "../files_skim/OniaFlowSkim_Jpsi_AOD_isMC1_NPOnly_HFNom_250217.root";
        tree->Add(f1.Data());
    }

    // ===== import Acc and Eff files ===== //
    // // TFile *fin_weight = new TFile(Form("eff_acc_inputs/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_Jpsi_PtW%d_tnp%d_210913.root", isPtW, isTnP), "read");
    // TFile *fin_correct = new TFile("../eff_acc/roots/mc_eff_vs_pt_cent_rap_prompt_pbpb_Jpsi_PtW1_tnp1_250221.root");

    // TH2D *h_correct[6];
    // // h_correct[0] = (TH1D *)fEff->Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_180_absy0_1p2", isTnP, isPtW));

    // h_correct[0] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent0to20_absy0_1p6");
    // h_correct[1] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent0to20_absy1p6_2p4");
    // h_correct[2] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent20_60_absy0_1p6");
    // h_correct[3] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent20_60_absy1p6_2p4");
    // h_correct[4] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent60to180_absy0_1p6");
    // h_correct[5] = (TH2D *)fin_correct->Get("mc_eff_vs_pt_TnP1_PtW1_cent60to180_absy1p6_2p4");


    // ===== Tree connection ===== //
    const int nMaxDimu = 1000;

    Int_t nDimuGen;
    float mass[nMaxDimu];
    float y[nMaxDimu];

    float pt[nMaxDimu];
    float pt1[nMaxDimu];
    float pt2[nMaxDimu];
    
    float eta[nMaxDimu];
    float eta1[nMaxDimu];
    float eta2[nMaxDimu];

    float phi[nMaxDimu];
    float phi1[nMaxDimu];
    float phi2[nMaxDimu];

    float ctau3D[nMaxDimu]; // = ctau3DTrue
    
    TBranch *b_nDimuGen;
    TBranch *b_mass;
    TBranch *b_y;
    
    TBranch *b_pt;
    TBranch *b_pt1;
    TBranch *b_pt2;

    TBranch *b_eta;
    TBranch *b_eta1;
    TBranch *b_eta2;

    TBranch *b_phi;
    TBranch *b_phi1;
    TBranch *b_phi2;

    TBranch *b_ctau3D;
    

    tree->SetBranchAddress("nDimuGen", &nDimuGen, &b_nDimuGen);
    tree->SetBranchAddress("mass", mass, &b_mass);
    tree->SetBranchAddress("y", y, &b_y);
    tree->SetBranchAddress("pt", pt, &b_pt);
    tree->SetBranchAddress("pt1", pt1, &b_pt1);
    tree->SetBranchAddress("pt2", pt2, &b_pt2);
    tree->SetBranchAddress("eta", eta, &b_eta);
    tree->SetBranchAddress("eta1", eta1, &b_eta1);
    tree->SetBranchAddress("eta2", eta2, &b_eta2);
    tree->SetBranchAddress("phi", phi, &b_phi);
    tree->SetBranchAddress("phi1", phi1, &b_phi1);
    tree->SetBranchAddress("phi2", phi2, &b_phi2);
    tree->SetBranchAddress("ctau3D", ctau3D, &b_ctau3D);
    

    // ===== Outfile ===== //
    TFile *wf = nullptr;
    if (isMC)
    {
        wf = new TFile(Form("../files_roodata/RooDataSet_miniAOD_isMC%d_%s_Jpsi_GenOnly_cent%i_%i_Effw%d_Accw%d_PtW%d_TnP%d_HF%s_%s.root",
                            isMC, MCtype_tag.Data(), cLow, cHigh, fEffW, fAccW, isPtW, isTnP, fCentSelHF.Data(), DATE.Data()),
                       "recreate");
    }
    else
    {
        cout << "Gen Only code needs to use MC sample" << endl;
        cout << "Check your input file" << endl;
        exit(1);
    }
    wf->cd();

    ////////////////////////////////////////////////////////////////////////
    //////////////////  build RooDataSet
    ////////////////////////////////////////////////////////////////////////
    RooRealVar *NumDimu = new RooRealVar("NumDimu", "number of dimuon", 0, 100, "");
    RooRealVar *massVar = new RooRealVar("mass", "mass variable", 1.0, 6.0, "GeV/c^{2}");
    RooRealVar *yVar = new RooRealVar("y", "rapidity of the dimuon pair", -5, 5, "");
    RooRealVar *ptVar = new RooRealVar("pt", "pt variable", 0, 100, "GeV/c");
    RooRealVar *pt1Var = new RooRealVar("pt1", "pt of muon+", 0, 500, "GeV/c");
    RooRealVar *pt2Var = (RooRealVar *)pt1Var->Clone("pt2");
    RooRealVar *etaVar = new RooRealVar("eta", "eta of muon+", -5, 5, "");
    RooRealVar *eta1Var = new RooRealVar("eta1", "eta of mupl", -5, 5, "");
    RooRealVar *eta2Var = new RooRealVar("eta2", "eta of mumi", -5, 5, "");
    // RooRealVar *eta2Var = (RooRealVar *)eta1Var->Clone("eta2");
    RooRealVar *phiVar = new RooRealVar("phi", "phi of muon+", -4, 4, "");
    RooRealVar *phi1Var = new RooRealVar("phi1", "phi of mupl", -4, 4, "");
    RooRealVar *phi2Var = new RooRealVar("phi2", "phi of mumi", -4, 4, "");
    RooRealVar *ctau3DVar = new RooRealVar("ctau3D", "c_{#tau}", 0, 20.0, "mm");

    RooArgSet *argSet = new RooArgSet(
        *NumDimu, *massVar, *yVar, *ptVar, *pt1Var, *pt2Var, *etaVar, *eta1Var, *eta2Var,
        *phiVar, *phi1Var, *phi2Var, *ctau3DVar);
    RooDataSet *dataSet = new RooDataSet("dataset", "", *argSet);


    // ===== Prepare event loop ===== //
    double weight_acc = 1;
    double weight_eff = 1;

    if (nEvt == -1)
        nEvt = tree->GetEntries();
    cout << "nEvt : " << nEvt << endl;

    int nDimuGen_all = 0;
    int nDimuGenPass = 0;
    int nDimuGen_one = 0;
    int nDimuGen_more = 0;
    int count_dimu = 0;

    // Begin Loop
    for (int i = 0; i < nEvt; i++)
    {
        tree->GetEntry(i);
        
        nDimuGenPass = 0;
        // ===== start dimuon loop ===== //

        for (int j = 0; j < nDimuGen; j++)
        {
            // check cut once more with mass cut
            if ( !((double)pt[j] < 50 && mass[j] > massLow && mass[j] < massHigh && abs(y[j]) < 2.4 && IsAcceptanceQQ(pt1[j], eta1[j]) && IsAcceptanceQQ(pt2[j], eta2[j])))
                continue;
            // IsAcceptanceQQ()에 pt, eta cut만 들어가 있다. Gen에 넣어야 하나? dimuon이 맞기는 해야하니까 일단 넣어보자.
            // if (!((double)pt[j] < 50 && mass[j] > massLow && mass[j] < massHigh && abs(y[j]) < 2.4 && IsAcceptanceQQ(pt1[j], eta1[j]) && IsAcceptanceQQ(pt2[j], eta2[j])))

            // ===== apply acc, eff weight ===== //
            weight_acc = 1;
            weight_eff = 1;
                        
            // if (fAccW)
            // {
            //     if (abs((double)y[j]) < 1.6)
            //     {
            //         weight_acc = getAccWeight(hAccPt[1], pt[j]);
            //     }
            //     else if (abs((double)y[j]) > 1.6 && abs((double)y[j]) < 2.4)
            //     {
            //         weight_acc = getAccWeight(hAccPt[2], pt[j]);
            //     }
            // }
            // if (fEffW)
            // {
            //     if (abs((double)y[j]) >= 0.0 && abs((double)y[j]) < 1.2)
            //     {
            //         weight_eff = getEffWeight(hEffPt[0], pt[j]);
            //     }
            //     else if (abs((double)y[j]) >= 1.2 && abs((double)y[j]) < 1.6)
            //     {
            //         weight_eff = getEffWeight(hEffPt[1], pt[j]);
            //     }
            //     else if (abs((double)y[j]) >= 1.6 && abs((double)y[j]) < 1.8)
            //     {
            //         weight_eff = getEffWeight(hEffPt[2], pt[j]);
            //     }
            //     else if (abs((double)y[j]) >= 1.8 && abs((double)y[j]) < 2.4)
            //     {
            //         weight_eff = getEffWeight(hEffPt[3], pt[j]);
            //     }
            // }

            // double weight_ = weight * weight_eff * weight_acc;

            NumDimu->setVal((int)nDimuGen);
            massVar->setVal((double)mass[j]);
            yVar->setVal((double)y[j]);
            ptVar->setVal((double)pt[j]);
            pt1Var->setVal((double)pt1[j]);
            pt2Var->setVal((double)pt2[j]);
            etaVar->setVal((double)eta[j]);
            eta1Var->setVal((double)eta1[j]);
            eta2Var->setVal((double)eta2[j]);
            phiVar->setVal((double)phi[j]);
            phi1Var->setVal((double)phi1[j]);
            phi2Var->setVal((double)phi2[j]);
            ctau3DVar->setVal((double)ctau3D[j]); // = ctau3DTrue
            // evtWeight->setVal((double)weight_);
            
            dataSet->add(*argSet);

            count_dimu++;
        } // end of dimuon loop
    }
    dataSet->Write();
    wf->Close();

    cout << "How many Jpsi??: " << count_dimu << endl;

    cout << "\n=================================\n";
    cout << "\n Finish Skim to Roodataset conversion\n";
    cout << "\n=================================\n";
    t->Stop();
    printf("RealTime=%f seconds, CpuTime=%f seconds\n", t->RealTime(), t->CpuTime());
}

// ===== Other functions ===== //
void GetHistSqrt(TH1D *h1, TH1D *h2)
{
    if (h1->GetNbinsX() != h2->GetNbinsX())
    {
        cout << "Inconsistent # of bins b/w histograms !! " << endl;
    }
    double content;
    double err;
    for (int i = 1; i <= h1->GetNbinsX(); i++)
    {
        content = 0;
        err = 0;
        content = h1->GetBinContent(i);
        err = h1->GetBinError(i);
        err = 0.5 * err * TMath::Power(content, -0.5);
        h2->SetBinContent(i, TMath::Sqrt(content));
        h2->SetBinError(i, err);
    }
}

void Get2DHistSqrt(TH2D *h1, TH2D *h2)
{
    if (h1->GetNbinsX() != h2->GetNbinsX())
    {
        cout << "Inconsistent # of bins b/w histograms !! " << endl;
    }
    double content;
    double err;
    for (int i = 1; i <= h1->GetNbinsX(); i++)
    {
        content = 0;
        err = 0;
        content = h1->GetBinContent(i);
        err = h1->GetBinError(i);
        err = 0.5 * err * TMath::Power(content, -0.5);
        h2->SetBinContent(i, TMath::Sqrt(content));
        h2->SetBinError(i, err);
    }
}

double GetHistBkg(TH1D *h1, TH1D *h2)
{
    double Nh1;
    double Nh2;
    double Nh12;
    Nh1 = h1->GetEntries();
    Nh2 = h2->GetEntries();
    Nh12 = Nh1 + Nh2;
    return Nh1 / (Nh1 + Nh2);
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