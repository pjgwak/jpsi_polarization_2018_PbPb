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
#include "../cms_headers/commonUtility.h"
#include "../cms_headers/cutsAndBin.h"
#include "../cms_headers/HiEvtPlaneList.h"
#include "../cms_headers/Style.h"
#include "../cms_headers/tdrstyle.C"
#include "../cms_headers/CMS_lumi_v2mass.C"
#include "../cms_headers/rootFitHeaders.h"

using namespace std;
using namespace RooFit;
using namespace hi;

static const long MAXTREESIZE = 1000000000000;

double getAccWeight(TH1D *h = 0, double pt = 0); // 1D case
double getEffWeight(TH1D *h = 0, double pt = 0);
double getWeight(TH2D *h, double angle, double pt);
// void GetHistSqrt(TH1D *h1 = 0, TH1D *h2 = 0); // Legacies. Not used.
// void Get2DHistSqrt(TH2D *h1 = 0, TH2D *h2 = 0);
// void GetHistBkg(TH1D *h1 = 0, TH1D *h2 = 0);

void cos_skim_to_rooDataset_data(
    int cLow = 0, int cHigh = 180,
    float massLow = 2.6, float massHigh = 3.5,
    // float massLow = 3.3, float massHigh =4.1,
    bool dimusign = true, bool isMC = true,
    bool fAccW = true, bool fEffW = true,
    bool isPtW = true, bool isTnP = true,
    int hiHFBinEdge = 0,
    int MCtype = 1) // 1: PR, 2: NP
{
  // Basic Setting
  gStyle->SetOptStat(0);

  string folder_label = "241014_temporal_Aall_Weightings_on_test";
  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(000000000);
  gROOT->ForceStyle();
  setTDRStyle();
  writeExtraText = true;
  int iPeriod = 2;
  int iPos = 33;
  TH1::SetDefaultSumw2();
  // TString kineLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, SiMuPtCut, cLow, cHigh, 1) ;
  // TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0, cLow, cHigh) ;
  TString dimusignString;
  if (dimusign)
    dimusignString = "OS";
  else if (!dimusign)
    dimusignString = "SS";

  TString sample;
  TString wName;
  TString cutName;

  // Not used
  if (MCtype == 1)
    wName = "PR";
  else if (MCtype == 2)
    wName = "NP";

  TChain *tree = new TChain("mmepevt");
  if (!isMC)
  {
    TString f1 = "skimmed_files/OniaFlowSkim_JpsiTrig_DBAll_miniAOD_isMC0_HFNom_240820.root";
    tree->Add(f1.Data());
    sample = "RD"; // Means RealData (?)
  }
  else if (MCtype == 1)
  {
    TString f1 = "skimmed_files/OniaFlowSkim_Jpsi_miniAOD_isMC1_Signal_HFNom_240910.root";
    tree->Add(f1.Data());
    sample = "MC_PR";
  }
  else if (MCtype == 2)
  {
    TString f1 = "skimmed_files/OniaFlowSkim_Jpsi_miniAOD_isMC1_NPOnly_HFNom_240822.root";
    tree->Add(f1.Data());
    sample = "MC_NP";
  }

  auto fEff = new TFile("/work/pjgwak/work/pub_2024_jpsi_pol_2018PbPb/data_selection/local_cos_eff_acc_study/weight_outputs/eff_pt_theta_hx_cent_0_to_180_Signal_pbpb_Jpsi_PtW1_tnp1_angle.root");
  TH2D *hEffPt_fwd[9];
  hEffPt_fwd[0] = (TH2D *)fEff->Get("heff_cos_pT_y1p6_2p4_cent0_20");
  hEffPt_fwd[1] = (TH2D *)fEff->Get("heff_cos_pT_y1p6_2p4_cent20_40");
  hEffPt_fwd[2] = (TH2D *)fEff->Get("heff_cos_pT_y1p6_2p4_cent40_60");
  hEffPt_fwd[3] = (TH2D *)fEff->Get("heff_cos_pT_y1p6_2p4_cent60_80");
  hEffPt_fwd[4] = (TH2D *)fEff->Get("heff_cos_pT_y1p6_2p4_cent80_100");
  hEffPt_fwd[5] = (TH2D *)fEff->Get("heff_cos_pT_y1p6_2p4_cent100_120");
  hEffPt_fwd[6] = (TH2D *)fEff->Get("heff_cos_pT_y1p6_2p4_cent120_140");
  hEffPt_fwd[7] = (TH2D *)fEff->Get("heff_cos_pT_y1p6_2p4_cent120_140");
  hEffPt_fwd[8] = (TH2D *)fEff->Get("heff_cos_pT_y1p6_2p4_cent160_180");

  TH2D *hEffPt_mid[9];
  hEffPt_mid[0] = (TH2D *)fEff->Get("heff_cos_pT_y0_1p6_cent0_20");
  hEffPt_mid[1] = (TH2D *)fEff->Get("heff_cos_pT_y0_1p6_cent20_40");
  hEffPt_mid[2] = (TH2D *)fEff->Get("heff_cos_pT_y0_1p6_cent40_60");
  hEffPt_mid[3] = (TH2D *)fEff->Get("heff_cos_pT_y0_1p6_cent60_80");
  hEffPt_mid[4] = (TH2D *)fEff->Get("heff_cos_pT_y0_1p6_cent80_100");
  hEffPt_mid[5] = (TH2D *)fEff->Get("heff_cos_pT_y0_1p6_cent100_120");
  hEffPt_mid[6] = (TH2D *)fEff->Get("heff_cos_pT_y0_1p6_cent120_140");
  hEffPt_mid[7] = (TH2D *)fEff->Get("heff_cos_pT_y0_1p6_cent120_140");
  hEffPt_mid[8] = (TH2D *)fEff->Get("heff_cos_pT_y0_1p6_cent160_180");

  TFile *fAcc = new TFile("/work/pjgwak/work/pub_2024_jpsi_pol_2018PbPb/data_selection/local_cos_eff_acc_study/weight_outputs/acceptance_Prompt_Jpsi_GenOnly_wgt1_PbPb_SysUp0_cos_hx.root");
  TH2D *hAccPt[2];
  hAccPt[0] = (TH2D *)fAcc->Get("hAcc_pT_theta_hx_mid");
  hAccPt[1] = (TH2D *)fAcc->Get("hAcc_pT_theta_hx_fwd");

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
  double cos_theta_hx[nMaxDimu];
  double weight;

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
  TBranch *b_pt1;
  TBranch *b_pt2;
  TBranch *b_weight;
  TBranch *b_cos_theta_hx;

  tree->SetBranchAddress("event", &event, &b_event);
  tree->SetBranchAddress("cBin", &cBin, &b_cBin);
  tree->SetBranchAddress("nDimu", &nDimu, &b_nDimu);
  tree->SetBranchAddress("recoQQsign", recoQQsign, &b_recoQQsign);
  tree->SetBranchAddress("vz", &vz, &b_vz);
  tree->SetBranchAddress("mass", mass, &b_mass);
  tree->SetBranchAddress("y", y, &b_y);
  tree->SetBranchAddress("pt", pt, &b_pt);
  tree->SetBranchAddress("pt1", pt1, &b_pt1);
  tree->SetBranchAddress("pt2", pt2, &b_pt2);
  tree->SetBranchAddress("eta", eta, &b_eta);
  tree->SetBranchAddress("eta1", eta1, &b_eta1);
  tree->SetBranchAddress("eta2", eta2, &b_eta2);
  tree->SetBranchAddress("ctau3D", ctau3D, &b_ctau3D);
  tree->SetBranchAddress("ctau3DErr", ctau3DErr, &b_ctau3DErr);
  tree->SetBranchAddress("cos_theta_hx", cos_theta_hx, &b_cos_theta_hx);
  tree->SetBranchAddress("weight", &weight, &b_weight);

  const int nV2Bin = 6;

  float massBinDiff[nV2Bin + 1] = {-3, -2, -1, 0, 1, 2, 3};
  float massBin_[nV2Bin + 1];

  // kineLabel = kineLabel + Form("_m%.1f-%.1f",massLow,massHigh) + "_" + dimusignString;
  // cout<<kineLabel<<endl;

  // Define drawing histogram
  TH1D *h_ljpsi;

  double nMassBin = 45;
  double nQBin = 100;
  double Q_avg_low = -10000;
  double Q_avg_high = 40000;
  double Q_avg_low_dimu = -200;
  double Q_avg_high_dimu = 200;
  h_ljpsi = new TH1D("h_ljBkg", ";l_{J/#psi(mm)};Counts/(0.03 mm)", 1050, -1.5, 2);

  const static int countMax = 1000000;

  ////////////////////////////////////////////////////////////////////////
  //////////////////  RooDataSet
  ////////////////////////////////////////////////////////////////////////
  RooRealVar *massVar = new RooRealVar("mass", "mass variable", 1.0, 6.0, "GeV/c^{2}");
  RooRealVar *ptVar = new RooRealVar("pt", "pt variable", 0, 100, "GeV/c");
  RooRealVar *yVar = new RooRealVar("y", "rapidity of the dimuon pair", -5, 5, "");
  RooRealVar *pt1Var = new RooRealVar("pt1", "pt of muon+", 0, 500, "GeV/c");
  RooRealVar *eta1Var = new RooRealVar("eta1", "eta of muon+", -4, 4, "");
  RooRealVar *pt2Var = (RooRealVar *)pt1Var->Clone("pt2");
  RooRealVar *eta2Var = (RooRealVar *)eta1Var->Clone("eta2");
  RooRealVar *cBinVar = new RooRealVar("cBin", "Centrality bin", -100, 500, "");
  RooRealVar *ep2Var = new RooRealVar("ep2", "2nd order event plane", -100, 100, "");
  RooRealVar *evtWeight = new RooRealVar("weight", "corr weight", 0, 10000, "");
  RooRealVar *recoQQ = new RooRealVar("recoQQsign", "qq sign", -1, 3, "");
  RooRealVar *ctau3DVar = new RooRealVar("ctau3D", "c_{#tau}", -100000.0, 100000.0, "mm");
  RooRealVar *ctau3DErrVar = new RooRealVar("ctau3DErr", "#sigma_{c#tau}", -100000.0, 100000.0, "mm");
  RooRealVar *ctau3DResVar = new RooRealVar("ctau3DRes", "c_{#tau}", -100000.0, 100000.0, "");
  RooRealVar *cos_theta_hxVar = new RooRealVar("cos_theta_hx", "cosTheta in HX", -1.0, 1.0, "");
  RooRealVar *NumDimu = new RooRealVar("NumDimu", "number of dimuon", 0, 100, "");
  RooArgSet *argSet = new RooArgSet(*massVar, *ptVar, *yVar, *pt1Var, *pt2Var, *eta1Var, *eta2Var, *evtWeight, *cos_theta_hxVar);
  argSet->add(*cBinVar);
  argSet->add(*recoQQ);
  argSet->add(*NumDimu);
  argSet->add(*ctau3DVar);
  argSet->add(*ctau3DErrVar);
  argSet->add(*ctau3DResVar);
  RooDataSet *dataSet = new RooDataSet("dataset", "", *argSet);

  ////////////////////////////////////////////////////////////////////////
  //////////////////  TreeDataSet
  ////////////////////////////////////////////////////////////////////////
  TString fCentSelHF = "HFNom";
  if (hiHFBinEdge == 1)
    fCentSelHF = "HFUp";
  else if (hiHFBinEdge == -1)
    fCentSelHF = "HFDown";

  TFile *wf = nullptr;

  if (isMC == true)
  {
    wf = new TFile(Form("roo_dataset/RooDataSet_miniAOD_isMC%d_Jpsi_%s_cent%i_%i_Effw%d_Accw%d_PtW%d_TnP%d_cos_weight_%s.root", isMC, wName.Data(), cLow, cHigh, fEffW, fAccW, isPtW, isTnP, folder_label.c_str()), "recreate");
  }
  else
  {
    wf = new TFile(Form("roo_dataset/RooDataSet_miniAOD_isMC%d_Jpsi_cent%i_%i_Effw%d_Accw%d_PtW%d_TnP%d_cos_weight_%s.root", isMC, cLow, cHigh, fEffW, fAccW, isPtW, isTnP, folder_label.c_str()), "recreate");
  }
  // TFile *newfile;
  // newfile = new TFile(Form("roo_dataset/HFFlowSkim_isMC%d_%s_cos_weight.root", isMC, fCentSelHF.Data()), "recreate");

  const static int nMaxDimu_ = 1000;
  int nDimu_;
  float mass_[nMaxDimu_];
  float ctau3D_[nMaxDimu_];
  float ctau3DErr_[nMaxDimu_];

  TTree *mmevttree = new TTree("mmepevt", "dimuonAndEventPlanes in event based");
  mmevttree->SetMaxTreeSize(MAXTREESIZE);
  mmevttree->Branch("nDimu_", &nDimu_, "nDimu_/I");
  mmevttree->Branch("mass_", mass_, "mass_[nDimu_]/F");
  mmevttree->Branch("ctau3D_", ctau3D_, "ctau3D_[nDimu_]/F");
  mmevttree->Branch("ctau3DErr_", ctau3DErr_, "ctau3DErr_[nDimu_]/F");

  double weight_acc = 1;
  double weight_eff = 1;

  Int_t nEvt = tree->GetEntries();
  cout << "nEvt : " << nEvt << endl;

  int nDimu_all = 0;
  int nDimuPass = 0;
  int nDimu_one = 0;
  int nDimu_more = 0;

  int counts[10] = {0};

  // Begin Loop
  for (int i = 0; i < nEvt; i++)
  {
    counts[0]++;
    tree->GetEntry(i);
    if (fabs(vz) > 15)
      continue;
    counts[1]++;
    nDimuPass = 0;
    // if(i==1000000) break;
    if (cBin >= cLow && cBin < cHigh)
    {
      counts[2]++;
      for (int j = 0; j < nDimu; j++)
      {
        // cout<<"Evt: "<<i<<", Cent: "<<cBin<<", mass: "<<mass[j]<<", pt: "<<pt[j]<<", pt1:"<<pt1[j]<<", pt2: "<<pt2[j]<<", eta1: "<<eta1[j]<<", eta2: "<<eta2[j]<<", vz: "<<vz<<endl;
        if (!((double)pt[j] < 50 && recoQQsign[j] == 0 && abs((double)y[j]) < 2.4 && IsAcceptanceQQ(pt1[j], eta1[j]) && IsAcceptanceQQ(pt2[j], eta2[j])))
          continue;
        nDimuPass++;
      }
      counts[3]++;
      nDimu_all++;
      if (nDimuPass > 1)
      {
        nDimu_more++;
        continue;
      }
      counts[4]++;
      if (nDimuPass == 1)
        nDimu_one++;
      counts[5]++;
      // Fill Dimuon Loop
      nDimu_ = 0;
      for (int j = 0; j < nDimu; j++)
      {
        if ((double)pt[j] < 50 && recoQQsign[j] == 0 && mass[j] > massLow && mass[j] < massHigh && abs(y[j]) < 2.4 && IsAcceptanceQQ(pt1[j], eta1[j]) && IsAcceptanceQQ(pt2[j], eta2[j]))
        {
          counts[6]++;
          weight_acc = 1;
          weight_eff = 1;
          //	cout << "pt : " << pt[j] << " y : " << y[j] << endl;
          if (fAccW)
          {
            if (abs((double)y[j]) >= 1.6 && abs((double)y[j]) < 2.4)
            {
              weight_acc = getWeight(hAccPt[1], cos_theta_hx[j], pt[j]);
            }
            else
            {
              weight_acc = getWeight(hAccPt[0], cos_theta_hx[j], pt[j]);
            }
          }
          if (fEffW)
          {
            if (abs((double)y[j]) < 1.6)
            {
              if (cBin > 0 && cBin < 20)
                weight_eff = getWeight(hEffPt_mid[0], cos_theta_hx[j], pt[j]);
              else if (cBin > 20 && cBin < 40)
                weight_eff = getWeight(hEffPt_mid[1], cos_theta_hx[j], pt[j]);
              else if (cBin > 40 && cBin < 60)
                weight_eff = getWeight(hEffPt_mid[2], cos_theta_hx[j], pt[j]);
              else if (cBin > 60 && cBin < 80)
                weight_eff = getWeight(hEffPt_mid[3], cos_theta_hx[j], pt[j]);
              else if (cBin > 80 && cBin < 100)
                weight_eff = getWeight(hEffPt_mid[4], cos_theta_hx[j], pt[j]);
              else if (cBin > 100 && cBin < 120)
                weight_eff = getWeight(hEffPt_mid[5], cos_theta_hx[j], pt[j]);
              else if (cBin > 120 && cBin < 140)
                weight_eff = getWeight(hEffPt_mid[6], cos_theta_hx[j], pt[j]);
              else if (cBin > 140 && cBin < 160)
                weight_eff = getWeight(hEffPt_mid[7], cos_theta_hx[j], pt[j]);
              else if (cBin > 160 && cBin < 180)
                weight_eff = getWeight(hEffPt_mid[8], cos_theta_hx[j], pt[j]);
            }
            else
            {
              if (cBin > 0 && cBin < 20)
                weight_eff = getWeight(hEffPt_fwd[0], cos_theta_hx[j], pt[j]);
              else if (cBin > 20 && cBin < 40)
                weight_eff = getWeight(hEffPt_fwd[1], cos_theta_hx[j], pt[j]);
              else if (cBin > 40 && cBin < 60)
                weight_eff = getWeight(hEffPt_fwd[2], cos_theta_hx[j], pt[j]);
              else if (cBin > 60 && cBin < 80)
                weight_eff = getWeight(hEffPt_fwd[3], cos_theta_hx[j], pt[j]);
              else if (cBin > 80 && cBin < 100)
                weight_eff = getWeight(hEffPt_fwd[4], cos_theta_hx[j], pt[j]);
              else if (cBin > 100 && cBin < 120)
                weight_eff = getWeight(hEffPt_fwd[5], cos_theta_hx[j], pt[j]);
              else if (cBin > 120 && cBin < 140)
                weight_eff = getWeight(hEffPt_fwd[6], cos_theta_hx[j], pt[j]);
              else if (cBin > 140 && cBin < 160)
                weight_eff = getWeight(hEffPt_fwd[7], cos_theta_hx[j], pt[j]);
              else if (cBin > 160 && cBin < 180)
                weight_eff = getWeight(hEffPt_fwd[8], cos_theta_hx[j], pt[j]);
            }
          }

          double weight_ = weight * weight_eff * weight_acc;

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
          cos_theta_hxVar->setVal((double)cos_theta_hx[j]);
          evtWeight->setVal((double)weight_);
          NumDimu->setVal((int)nDimu);
          // cout<<"Evt: "<<j<<", Cent: "<<cBin<<", mass: "<<mass[j]<<", pt: "<<pt[j]<<", pt1:"<<pt1[j]<<", pt2: "<<pt2[j]<<", eta1: "<<eta1[j]<<", eta2: "<<eta2[j]<<", vz: "<<vz<<endl;
          dataSet->add(*argSet);

          mass_[j] = mass[j];
          ctau3D_[j] = ctau3D[j];
          ctau3DErr_[j] = ctau3DErr[j];
          nDimu_++;
        }
        if (nDimu_ > 0)
          mmevttree->Fill();
      }
    }
  }

  // newfile->cd();
  // mmevttree->Write();
  // newfile->Close();

  wf->cd();
  dataSet->Write();
  wf->Close();

  cout << "How many Jpsi??: " << nDimu_one << endl;
  cout << "0 : " << counts[0] << endl;
  cout << "1 : " << counts[1] << endl;
  cout << "2 : " << counts[2] << endl;
  cout << "3 : " << counts[3] << endl;
  cout << "4 : " << counts[4] << endl;
  cout << "5 : " << counts[5] << endl;
  cout << "6 : " << counts[6] << endl;
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

void GetHistBkg(TH1D *h1, TH1D *h2)
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

double getWeight(TH2D *h, double angle, double pt)
{
  // cout<<"pt: "<<pt<<endl;
  // double binN = h->FindBin(fabs(angle), pt);
  Int_t x = h->GetXaxis()->FindBin(fabs(angle));
  Int_t y = h->GetYaxis()->FindBin(pt);
  double weight_ = 1. / (h->GetBinContent(x, y));
  return weight_;
}
// ===== Other functions ===== //