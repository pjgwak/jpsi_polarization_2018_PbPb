#include "RooMakerRun3PbPb.h"
#include <iostream>
#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "TMath.h"
#include "../headers/HiEvtPlaneList.h"

using namespace std;
// using namespace RooFit;
using namespace hi;

static const long MAXTREESIZE = 1000000000000;

RooMakerRun3PbPb::RooMakerRun3PbPb(int cLow, int cHigh, float massLow, float massHigh,
                                   bool isMC, bool isAccW, bool isEffW,
                                   bool isTnP, bool isPtW, int hiHFBinEdge, int mcType,
                                   std::string fInputPath, std::string inputChainName, std::string userTag)
    : cLow(cLow), cHigh(cHigh), massLow(massLow), massHigh(massHigh),
      isMC(isMC), isAccW(isAccW), isEffW(isEffW),
      isTnP(isTnP), isPtW(isPtW), hiHFBinEdge(hiHFBinEdge), mcType(mcType), 
      fInputPath(fInputPath), inputChainName(inputChainName), userTag(userTag)
{
  // TH1::SetDefaultSumw2();
  // gStyle->SetOptStat(0); // no histogram stat box
  // gStyle->SetEndErrorSize(0); // RooPlot errbar size
  // setTDRStyle(); // use TDRStyle(). need header
  // gROOT->ForceStyle();
}

RooMakerRun3PbPb::~RooMakerRun3PbPb()
{
  delete fEff;
  delete fAcc;
  delete inChain;
  delete outfile;

  // RooFit variables
  delete massVar;
  delete ptVar;
  delete yVar;
  delete pt1Var;
  delete pt2Var;
  delete eta1Var;
  delete eta2Var;
  delete cBinVar;
  delete evtWeight;
  delete recoQQ;
  delete ctau3DVar;
  delete ctau3DErrVar;
  delete ctau3DResVar;
  // delete ctau3D2SVar;
  // delete ctau3DErr2SVar;
  // delete ctau3DRes2SVar;
  delete NumDimu;

  delete argSet;
  delete dataSet;

  delete cosHXVar;
  delete phiHXVar;
  delete cosCSVar;
  delete phiCSVar;
  delete cosEPVar;
  delete phiEPVar;
}

void RooMakerRun3PbPb::init()
{
  std::cout << "===== init() =====\n\n";
  setLabels();
  openInputFile();
  setOutfile();
  // loadEffHistograms();
  // loadAccHistograms();
  setupTree();
  defineRooVariables();
  setupRooDataSet();
}

void RooMakerRun3PbPb::run()
{
  std::cout << "===== run() =====\n\n";
  loopAndFillData();
  saveResults();
}

void RooMakerRun3PbPb::setLabels()
{
  std::cout << "===== setLabels() =====\n\n";
  sampleTag = (isMC) ? "MC" : "Data";
  mcTag = (mcType == 0) ? "PR" : "NP";

  gSystem->mkdir(Form("../roots"));

  if (hiHFBinEdge == 1)
    centHfTag = "HFUp";
  else if (hiHFBinEdge == -1)
    centHfTag = "HFDown";
  else centHfTag = "HFNom";

  // we always use Opposit Sign
  // dimusignString;
  // if (dimusign)
  //   dimusignString = "OS";
  // else if (!dimusign)
  //   dimusignString = "SS";
}

void RooMakerRun3PbPb::openInputFile()
{
  std::cout << "===== openInputFile() =====\n\n";
  inChain = new TChain(inputChainName.c_str());
  inChain->Add(fInputPath.c_str());
  if (inChain->GetEntries() == 0)
  {
    std::cerr << "[ERROR]: Failed to open input file: " << fInputPath << "\n";
    std::cerr << "       : Tree name: " << inputChainName << "\n";
    exit(1);
  }
  // inChain->Print();

  // Todo: make a naming function. Use fEffPath and fAccPath
  if (isEffW)
  {
    fEff = new TFile(Form("./Eff_Acc/mc_eff_vs_pt_cent_0_to_200_rap_prompt_pbpb_psi2s_PtW%d_tnp%d_20221116.root", isPtW, isTnP), "read");

    // if (!fEff || fEff->IsZombie())
    // {
    //   std::cerr << "[Error] Failed to open efficiency file: " << effName << std::endl;
    //   exit(1);
    // }
  }

  if (isAccW) {
    fAcc = new TFile(Form("./Eff_Acc/acceptance_Prompt_Psi2s_GenOnly_wgt%d_20220125.root", isPtW), "read");
    // if (!fAcc || fAcc->IsZombie())
    // {
    //   std::cerr << "[Error] Failed to open acceptance file: " << accName << std::endl;
    //   exit(1);
    // }
  }
}

void RooMakerRun3PbPb::setupTree()
{
  std::cout << "===== setupTree() =====\n\n";
  inChain->SetBranchAddress("event", &event);
  inChain->SetBranchAddress("cBin", &cBin);
  inChain->SetBranchAddress("nDimu", &nDimu);
  inChain->SetBranchAddress("recoQQsign", recoQQsign);
  inChain->SetBranchAddress("vz", &vz);
  inChain->SetBranchAddress("mass", mass);
  inChain->SetBranchAddress("pt", pt);
  inChain->SetBranchAddress("y", y);
  inChain->SetBranchAddress("pt1", pt1);
  inChain->SetBranchAddress("pt2", pt2);
  inChain->SetBranchAddress("eta1", eta1);
  inChain->SetBranchAddress("eta2", eta2);
  inChain->SetBranchAddress("ctau3D", ctau3D);
  inChain->SetBranchAddress("ctau3DErr", ctau3DErr);
  // inChain->SetBranchAddress("ctau3D2S", ctau3D2S);
  // inChain->SetBranchAddress("ctau3DErr2S", ctau3DErr2S);
  // inChain->SetBranchAddress("ctau3DRes2S", ctau3DRes2S);
  inChain->SetBranchAddress("weight", &weight);

  inChain->SetBranchAddress("cosHX", cosHX);
  inChain->SetBranchAddress("phiHX", phiHX);
  inChain->SetBranchAddress("cosCS", cosCS);
  inChain->SetBranchAddress("phiCS", phiCS);
  inChain->SetBranchAddress("cosEP", cosEP);
  inChain->SetBranchAddress("phiEP", phiEP);
}

void RooMakerRun3PbPb::defineRooVariables()
{
  std::cout << "===== defineRooVariables() =====\n\n";
  massVar = new RooRealVar("mass", "mass", 1.0, 6.0);
  ptVar = new RooRealVar("pt", "pt", 0, 100);
  yVar = new RooRealVar("y", "y", -5, 5);
  pt1Var = new RooRealVar("pt1", "pt1", 0, 500);
  pt2Var = new RooRealVar("pt2", "pt2", 0, 500);
  eta1Var = new RooRealVar("eta1", "eta1", -4, 4);
  eta2Var = new RooRealVar("eta2", "eta2", -4, 4);
  cBinVar = new RooRealVar("cBin", "cBin", -100, 500);
  evtWeight = new RooRealVar("weight", "weight", 0, 10000);
  recoQQ = new RooRealVar("recoQQsign", "recoQQsign", -1, 3);
  ctau3DVar = new RooRealVar("ctau3D", "ctau3D", -1e5, 1e5);
  ctau3DErrVar = new RooRealVar("ctau3DErr", "ctau3DErr", -1e5, 1e5);
  ctau3DResVar = new RooRealVar("ctau3DRes", "ctau3DRes", -1e5, 1e5);
  // ctau3D2SVar = new RooRealVar("ctau3D2S", "ctau3D2S", -1e5, 1e5);
  // ctau3DErr2SVar = new RooRealVar("ctau3DErr2S", "ctau3DErr2S", -1e5, 1e5);
  // ctau3DRes2SVar = new RooRealVar("ctau3DRes2S", "ctau3DRes2S", -1e5, 1e5);
  NumDimu = new RooRealVar("NumDimu", "NumDimu", 0, 100);

  cosHXVar = new RooRealVar("cosHX", "cosHX", -1.5, 1.5);
  phiHXVar = new RooRealVar("phiHX", "phiHX", -TMath::Pi(), TMath::Pi());
  cosCSVar = new RooRealVar("cosCS", "cosCS", -1.5, 1.5);
  phiCSVar = new RooRealVar("phiCS", "phiCS", -TMath::Pi(), TMath::Pi());
  cosEPVar = new RooRealVar("cosEP", "cosEP", -1.5, 1.5);
  phiEPVar = new RooRealVar("phiEP", "phiEP", -TMath::Pi(), TMath::Pi());

  argSet = new RooArgSet(*massVar, *ptVar, *yVar, *pt1Var, *pt2Var, *eta1Var, *eta2Var, *evtWeight);
  argSet->add(*cBinVar);
  argSet->add(*recoQQ);
  argSet->add(*NumDimu);
  argSet->add(*ctau3DVar);
  argSet->add(*ctau3DErrVar);
  argSet->add(*ctau3DResVar);
  // argSet->add(*ctau3D2SVar);
  // argSet->add(*ctau3DErr2SVar);
  // argSet->add(*ctau3DRes2SVar);

  argSet->add(*cosHXVar);
  argSet->add(*phiHXVar);
  argSet->add(*cosCSVar);
  argSet->add(*phiCSVar);
  argSet->add(*cosEPVar);
  argSet->add(*phiEPVar);
}

void RooMakerRun3PbPb::setupRooDataSet()
{
  std::cout << "===== setupRooDataSet() =====\n\n";
  dataSet = new RooDataSet("dataset", " a dataset", *argSet);

  // weight should be turned on by fitting step!
  // dataSet = new RooDataSet("dataset", " a dataset", *argSet, RooFit::WeightVar(*evtWeight));
}

void RooMakerRun3PbPb::setOutfile()
{
  std::cout << "===== setOutfile() =====\n\n";
  // outfile = new TFile(Form("roots/OniaFlowSkim_isMC%d_%s_230119.root", isMC, centHfTag.c_str()), "recreate");

  // centHfTag.c_str()
  outfile= new TFile(Form("../roots/OniaRooDataSet_miniAOD_isMC%d_Jpsi_cent%i_%i_Effw%d_Accw%d_PtW%d_TnP%d_%s.root", isMC, cLow, cHigh, isEffW, isAccW, isPtW, isTnP, userTag.c_str()), "recreate");
}

void RooMakerRun3PbPb::loopAndFillData()
{
  std::cout << "===== loopAndFillData() =====\n\n";
  long nEvt = inChain->GetEntries();
  for (long i = 0; i < nEvt; i++)
  {
    inChain->GetEntry(i);
    if (fabs(vz) > 15)
      continue;

    for (int j = 0; j < nDimu; j++)
    {
      // we apply these cut in fitting step too.
      // if (pt[j] >= 50 || fabs(y[j]) >= 2.4)
      //   continue;

      massVar->setVal(mass[j]);
      ptVar->setVal(pt[j]);
      yVar->setVal(y[j]);
      pt1Var->setVal(pt1[j]);
      pt2Var->setVal(pt2[j]);
      eta1Var->setVal(eta1[j]);
      eta2Var->setVal(eta2[j]);
      cBinVar->setVal(cBin);
      recoQQ->setVal(recoQQsign[j]);
      ctau3DVar->setVal(ctau3D[j]);
      ctau3DErrVar->setVal(ctau3DErr[j]);
      ctau3DResVar->setVal(ctau3D[j] / (ctau3DErr[j] > 0 ? ctau3DErr[j] : 1));
      // ctau3D2SVar->setVal(ctau3D2S[j]);
      // ctau3DErr2SVar->setVal(ctau3DErr2S[j]);
      // ctau3DRes2SVar->setVal(ctau3DRes2S[j]);
      evtWeight->setVal(weight);
      NumDimu->setVal(nDimu);

      cosHXVar->setVal(cosHX[j]);
      phiHXVar->setVal(phiHX[j]);
      cosCSVar->setVal(cosCS[j]);
      phiCSVar->setVal(phiCS[j]);
      cosEPVar->setVal(cosEP[j]);
      phiEPVar->setVal(phiEP[j]);

      dataSet->add(*argSet);
    }
  }
}

// Todo: make this function - include weighitng
// void RooMakerRun3PbPb::loopAndFillData()
// {
//   std::cout << "===== loopAndFillData() =====\n\n";
//   // double weight_acc = 1;
//   // double weight_eff = 1;

//   long nEvt = inChain->GetEntries();
//   cout << "nEvt : " << nEvt << endl;

//   int nDimu_all = 0;
//   int nDimuPass = 0;
//   int nDimu_one = 0;
//   int nDimu_more = 0;

//   // Begin Loop
//   for (long i = 0; i < nEvt; i++) {
//     inChain->GetEntry(i);
    
//     if (fabs(vz) > 15)
//       continue;
//     nDimuPass = 0;
//     // if(i==1000000) break;

//     if (cBin >= cLow && cBin < cHigh)
//     {

//       // default Algorithm already got this cut options - keep it now because it's a test
//       for (int j = 0; j < nDimu; j++)
//       {
//         // cout<<"Evt: "<<i<<", Cent: "<<cBin<<", mass: "<<mass[j]<<", pt: "<<pt[j]<<", pt1:"<<pt1[j]<<", pt2: "<<pt2[j]<<", eta1: "<<eta1[j]<<", eta2: "<<eta2[j]<<", vz: "<<vz<<endl;
//         if (!((double)pt[j] < 50 && recoQQsign[j] == 0 && abs((double)y[j]) < 2.4 && IsAcceptanceQQ(pt1[j], eta1[j]) && IsAcceptanceQQ(pt2[j], eta2[j])))
//           continue;
//         nDimuPass++;
//       }

//       // counting
//       nDimu_all++;

//       // counting
//       if (nDimuPass > 1)
//       {
//         nDimu_more++;
//         continue;
//       }

//       // counting
//       if (nDimuPass == 1) nDimu_one++;
      
//       // --- Fill Dimuon Loop ---
//       nDimu_ = 0;
//       for (int j = 0; j < nDimu; j++)
//       {
//         if ((double)pt[j] < 50 && recoQQsign[j] == 0 && mass[j] > massLow && mass[j] < massHigh && abs(y[j]) < 2.4 && IsAcceptanceQQ(pt1[j], eta1[j]) && IsAcceptanceQQ(pt2[j], eta2[j]))
//         {
//           weight_acc = 1;
//           weight_eff = 1;
//           //	cout << "pt : " << pt[j] << " y : " << y[j] << endl;
//           if (isAccW)
//           {
//             if (abs((double)y[j]) < 1.6)
//             {
//               weight_acc = getAccWeight(hAccPt[1], pt[j]);
//             }
//             else if (abs((double)y[j]) > 1.6 && abs((double)y[j]) < 2.4)
//             {
//               weight_acc = getAccWeight(hAccPt[2], pt[j]);
//             }
//           }
//           if (isEffW)
//           {
//             if (abs((double)y[j]) >= 0.0 && abs((double)y[j]) < 1.2)
//             {
//               weight_eff = getEffWeight(hEffPt[0], pt[j]);
//             }
//             else if (abs((double)y[j]) >= 1.2 && abs((double)y[j]) < 1.6)
//             {
//               weight_eff = getEffWeight(hEffPt[1], pt[j]);
//             }
//             else if (abs((double)y[j]) >= 1.6 && abs((double)y[j]) < 1.8)
//             {
//               weight_eff = getEffWeight(hEffPt[2], pt[j]);
//             }
//             else if (abs((double)y[j]) >= 1.8 && abs((double)y[j]) < 2.4)
//             {
//               weight_eff = getEffWeight(hEffPt[3], pt[j]);
//             }
//           }

//           double weight_ = weight * weight_eff * weight_acc;

//           recoQQ->setVal((int)recoQQsign[j]);
//           massVar->setVal((double)mass[j]);
//           ptVar->setVal((double)pt[j]);
//           yVar->setVal((double)y[j]);
//           pt1Var->setVal((double)pt1[j]);
//           eta1Var->setVal((double)eta1[j]);
//           pt2Var->setVal((double)pt2[j]);
//           eta2Var->setVal((double)eta2[j]);
//           cBinVar->setVal((double)cBin);
//           ctau3DVar->setVal((double)ctau3D[j]);
//           ctau3DErrVar->setVal((double)ctau3DErr[j]);
//           ctau3DResVar->setVal((double)ctau3D[j] / ctau3DErr[j]);
//           ctau3D2SVar->setVal((double)ctau3D[j]);
//           ctau3DErr2SVar->setVal((double)ctau3DErr[j]);
//           ctau3DRes2SVar->setVal((double)ctau3D2S[j]);
//           evtWeight->setVal((double)weight_);
//           NumDimu->setVal((int)nDimu);
//           // cout<<"Evt: "<<j<<", Cent: "<<cBin<<", mass: "<<mass[j]<<", pt: "<<pt[j]<<", pt1:"<<pt1[j]<<", pt2: "<<pt2[j]<<", eta1: "<<eta1[j]<<", eta2: "<<eta2[j]<<", vz: "<<vz<<endl;
//           dataSet->add(*argSet);

//           mass_[j] = mass[j];
//           ctau3D_[j] = ctau3D[j];
//           ctau3DErr_[j] = ctau3DErr[j];
//           nDimu_++;
//         }
//       }
//     }
//   }
//   // cout << "How many Jpsi??: " << nDimu_one << endl;
// }

void RooMakerRun3PbPb::saveResults()
{
  std::cout << "===== saveResults() =====\n\n";
  outfile->cd();
  
  dataSet->Write();
  outfile->Close();
}

// --- Make later ---
void RooMakerRun3PbPb::loadEffHist()
{
  std::cout << "===== loadEffHist() =====\n\n";
  // TH1D *hEffPt[4];
  // hEffPt[0] = (TH1D *)fEff->Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_200_absy0_1p2", isTnP, isPtW));
  // hEffPt[1] = (TH1D *)fEff->Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_200_absy1p2_1p6", isTnP, isPtW));
  // hEffPt[2] = (TH1D *)fEff->Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_200_absy1p6_1p8", isTnP, isPtW));
  // hEffPt[3] = (TH1D *)fEff->Get(Form("mc_eff_vs_pt_TnP%d_PtW%d_cent_0_to_200_absy1p8_2p4", isTnP, isPtW));
}

void RooMakerRun3PbPb::loadAccHist()
{
  std::cout << "===== loadAccHist() =====\n\n";
  // TH1D *hAccPt[3];
  // hAccPt[0] = (TH1D *)fAcc->Get("hAccPt_2021_ally");
  // hAccPt[1] = (TH1D *)fAcc->Get("hAccPt_2021_midy");
  // hAccPt[2] = (TH1D *)fAcc->Get("hAccPt_2021_Fory");
}

void RooMakerRun3PbPb::getEffWeight()
{
  std::cout << "===== getEffWeight() =====\n\n";
  // double binN = h->FindBin(pt);
  // double weight_ = 1. / (h->GetBinContent(binN));
  // return weight_;
}

void RooMakerRun3PbPb::getAccWeight()
{
  std::cout << "===== getAccWeight() =====\n\n";
  // double binN = h->FindBin(pt);
  // double weight_ = 1. / (h->GetBinContent(binN));
  // return weight_;
}