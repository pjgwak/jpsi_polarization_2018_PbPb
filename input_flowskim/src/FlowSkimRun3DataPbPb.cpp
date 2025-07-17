#include "FlowSkimRun3DataPbPb.h"
#include <iostream>
#include <algorithm>
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "../../headers/polFlowSkimUtility.h"

// 
// 
// 

// #include "commonUtility.h"
// #include "HiEvtPlaneList.h"
// #include "cutsAndBin.h"
// #include "tnp_weight_lowptPbPb.h"

// using namespace std;
// using namespace hi;

// Todo: check variables where is it? meanings?
// void SkimTree_Event_Psi2S_5p36TeV(int nevt = -1, bool isMC = false, int kTrigSel = kTrigJpsi, int hiHFBinEdge = 0, int PDtype = 1);

static const long MAXTREESIZE = 1000000000000; // 1 TB

FlowSkimRun3DataPbPb::FlowSkimRun3DataPbPb(
  bool isMC, int kTrigSel, int hiHFBinEdge, int PDtype,
  std::string inputFilePath, std::string inputTreeName, std::string outputFilePath,
  std::string outputTreeName, std::string userTag)
    : isMC(isMC), kTrigSel(kTrigSel), hiHFBinEdge(hiHFBinEdge), PDtype(PDtype), inputFilePath(inputFilePath),
    outputFilePath(outputFilePath), outputTreeName(outputTreeName), inputTreeName(inputTreeName), userTag(userTag)
{
  // std::cout << " Index of " << EPNames[HFm2] << " = " << HFm2 << "\n";
  // std::cout << " Index of " << EPNames[HFp2] << " = " << HFp2 << "\n";
  // std::cout << " Index of " << EPNames[trackmid2] << " = " << trackmid2 << "\n";
}

FlowSkimRun3DataPbPb::~FlowSkimRun3DataPbPb()
{
  delete inputTree;
  delete outputTree;
  delete outputFile;
  delete JP_Reco;
  delete mupl_Reco;
  delete mumi_Reco;
  delete Reco_QQ_4mom;
  delete Reco_mu_4mom;
}

void FlowSkimRun3DataPbPb::run(long nevt)
{
  std::cout << "===== run() =====\n\n";
  setLabels();
  openInputFile();
  setInputBranches();
  setOutputFile();
  setOutputBranches();

  // event loop
  long nEntries = inputTree->GetEntries();
  totalEvents = (nevt == -1) ? nEntries : std::min(nevt, nEntries);
  std::cout << "Total events: " << totalEvents << "\n";

  std::cout << "===== processEvent() =====\n\n";
  for (long iev = 0; iev < totalEvents; ++iev)
    processEvent(iev);
  
  saveResults();

  // output summary
  // std::cout << "count " << count << "\n";
  // std::cout << "count soft " << count_soft << "\n";
  // std::cout << "count vtx " << count_vtx << "\n";
  // std::cout << "counttnp " << counttnp << "\n";
}

void FlowSkimRun3DataPbPb::setLabels()
{
  std::cout << "===== setLabels() =====\n\n";
  sampleType = (PDtype == 1) ? "DB" : "DBPeri";

  if (hiHFBinEdge == 1)
    centHfTag = "HFUp";
  else if (hiHFBinEdge == -1)
    centHfTag = "HFDown";
  else centHfTag = "HFNom";

  if (kTrigSel == 12) // Todo: use cutsAndBins
    trigIndx = 0;
  
  // not used
  // else if (kTrigSel == kTrigUps)
  //   trigIndx = 1;
  // else if (kTrigSel == kTrigL1DBOS40100)
  //   trigIndx = 2;
  // else if (kTrigSel == kTrigL1DB50100)
  //   trigIndx = 3;
}

void FlowSkimRun3DataPbPb::openInputFile()
{
  std::cout << "===== openInputFile() =====\n\n";

  inputTree = new TChain(inputTreeName.c_str());
  inputTree->Add(inputFilePath.c_str());
  if (inputTree->GetEntries() == 0)
  {
    std::cerr << "[ERROR]: Failed to open input file: " << inputFilePath << "\n";
    std::cerr << "       : Tree name: " << inputTreeName << "\n";
    exit(1);
  }
}

void FlowSkimRun3DataPbPb::setInputBranches()
{
  std::cout << "===== setInputBranches() =====\n\n";
  // TClonesArray declarion - must do it befor connecting
  Reco_QQ_4mom = new TClonesArray("TLorentzVector");
  Reco_mu_4mom = new TClonesArray("TLorentzVector");


  inputTree->SetBranchAddress("runNb", &runNb);
  inputTree->SetBranchAddress("eventNb", &eventNb);
  inputTree->SetBranchAddress("LS", &LS);
  inputTree->SetBranchAddress("zVtx", &zVtx);
  inputTree->SetBranchAddress("Centrality", &Centrality);
  inputTree->SetBranchAddress("HLTriggers", &HLTriggers);
  inputTree->SetBranchAddress("SumET_HF", &SumET_HF);

  // Reco_QQ 
  inputTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
  inputTree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
  inputTree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);
  inputTree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign);
  inputTree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig);
  inputTree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb);
  inputTree->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D);
  inputTree->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D);
  inputTree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom);

  // Reco_mu
  inputTree->SetBranchAddress("Reco_mu_size", &Reco_mu_size);
  inputTree->SetBranchAddress("Reco_mu_charge", Reco_mu_charge);
  inputTree->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity);
  inputTree->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight);
  inputTree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits);
  inputTree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits);
  inputTree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched);
  inputTree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global);
  inputTree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy);
  inputTree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr);
  inputTree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz);
  inputTree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr);
  inputTree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea);
  inputTree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea);
  inputTree->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits);
  inputTree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType);
  inputTree->SetBranchAddress("Reco_mu_type", Reco_mu_type);
  inputTree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig);
  inputTree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom);

  // MC 
  if (isMC)
  {
    inputTree->SetBranchAddress("Gen_weight", &Gen_weight);
    inputTree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen);
  }
}

void FlowSkimRun3DataPbPb::setOutputFile()
{
  std::cout << "===== setOutputFile() =====\n\n";
  // Todo: change naming rule
  outputFile = new TFile(Form("OniaFlowSkim_isMC%d_%s.root", isMC, userTag.c_str()), "recreate");
  // outputFile = new TFile(Form("OniaFlowSkim_%sTrig_5p36TeV_%sPD_miniAOD_isMC%d_%s_231019.root", fTrigName[trigIndx].Data(), sampleType.Data(), isMC, centHfTag.Data()), "recreate");
}

void FlowSkimRun3DataPbPb::setOutputBranches()
{
  std::cout << "===== setOutputBranches() =====\n\n";
  outputTree = new TTree(outputTreeName.c_str(), "");
  outputTree->SetMaxTreeSize(MAXTREESIZE);

  outputTree->Branch("event", &evt, "event/I");
  outputTree->Branch("runN", &runN, "runN/I");
  outputTree->Branch("lumi", &lumi, "lumi/I");
  outputTree->Branch("cBin", &cBin, "cBin/I");
  outputTree->Branch("vz", &vz, "vz/F");
  outputTree->Branch("nDimu", &nDimu, "nDimu/I");

  outputTree->Branch("mass", mass, "mass[nDimu]/F");
  outputTree->Branch("pt", pt, "pt[nDimu]/F");
  outputTree->Branch("y", y, "y[nDimu]/F");
  outputTree->Branch("phi", phi, "phi[nDimu]/F");
  outputTree->Branch("eta", eta, "eta[nDimu]/F");

  outputTree->Branch("pt1", pt1, "pt1[nDimu]/F");
  outputTree->Branch("pt2", pt2, "pt2[nDimu]/F");
  outputTree->Branch("eta1", eta1, "eta1[nDimu]/F");
  outputTree->Branch("eta2", eta2, "eta2[nDimu]/F");
  outputTree->Branch("phi1", phi1, "phi1[nDimu]/F");
  outputTree->Branch("phi2", phi2, "phi2[nDimu]/F");

  outputTree->Branch("recoQQsign", recoQQsign, "recoQQsign[nDimu]/I");
  outputTree->Branch("ctau3D", ctau3D, "ctau3D[nDimu]/F");
  outputTree->Branch("ctau3DErr", ctau3DErr, "ctau3DErr[nDimu]/F");
  outputTree->Branch("ctau3DRes", ctau3DRes, "ctau3DRes[nDimu]/F");
  // outputTree->Branch("ctau3D2S", ctau3D2S, "ctau3D2S[nDimu]/F");
  // outputTree->Branch("ctau3DErr2S", ctau3DErr2S, "ctau3DErr2S[nDimu]/F");
  // outputTree->Branch("ctau3DRes2S", ctau3DRes2S, "ctau3DRes2S[nDimu]/F");

  outputTree->Branch("weight", &weight, "weight/D");
  if (isMC)
    outputTree->Branch("TnPweight", TnPweight, "TnPweight[nDimu]/D");
}


void FlowSkimRun3DataPbPb::processEvent(long iev)
{
  inputTree->GetEntry(iev);

  // print progress
  if (iev % 100000 == 0)
    std::cout << ">>>>> EVENT " << iev << " / " << inputTree->GetEntries() << " ("
              << (int)(100. * iev / inputTree->GetEntries()) << "%) <<<<<\n";

  // event value
  evt = eventNb;
  runN = runNb;
  lumi = LS;
  vz = zVtx;

  // decide centrality bins 
  if (hiHFBinEdge == 0)
    cBin = getHiBinFromhiHF(SumET_HF);
  else if (hiHFBinEdge == 1)
    cBin = getHiBinFromhiHF_Up(SumET_HF);
  else if (hiHFBinEdge == -1)
    cBin = getHiBinFromhiHF_Down(SumET_HF);
  if (cBin == -999)
  {
    std::cout << "ERROR!!! No HF Centrality Matching!!" << "\n";
    return;
  }
  
  nDimu = 0; // number of dimuon in this event

  // apply high level trigger
  // if((  (HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ){ std::cout << "iev : " << iev << " HLTriggers : " << HLTriggers << " kTrigSel : " << kTrigSel << "\n";	}
  if (!((HLTriggers & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel))))
    return;

  for (Int_t irqq = 0; irqq < Reco_QQ_size; ++irqq)
  {
    // --- trigger --- 
    // if( ( (Reco_QQ_trig[irqq]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ){ std::cout << "Reco_QQ_trig : " << Reco_QQ_trig[irqq] << "\n";}
    if (!((Reco_QQ_trig[irqq] & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel))))
      continue;

    // --- soft muon cut ---
    bool passMuonTypePl = true;
    passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]] & ((int)pow(2, 1)));
    passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]] & ((int)pow(2, 3)));

    bool passMuonTypeMi = true;
    passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]] & ((int)pow(2, 1)));
    passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]] & ((int)pow(2, 3)));

    bool muplSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]]==true) &&
        (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]] > 5) &&
        (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]] > 0) &&
        (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]]) < 0.3) &&
        (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]]) < 20.) &&
        passMuonTypePl //                       &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true)
    );

    bool mumiSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]]==true) &&
        (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5) &&
        (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0) &&
        (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]]) < 0.3) &&
        (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]]) < 20.) &&
        passMuonTypeMi //                        &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true)
    );

    if ((muplSoft && mumiSoft))
      count_soft++;
    if (!(muplSoft && mumiSoft))
      continue;

    // --- vertex probability ---
    if (!(Reco_QQ_VtxProb[irqq] < 0.01)) // counter
      count_vtx++;
    if (Reco_QQ_VtxProb[irqq] < 0.01)
      continue;
 
    // --- legacy - array style ---
    // float mupl_pt, mumi_pt;
    // float mupl_phi, mumi_phi;
    // float mupl_eta, mumi_eta;
    // float mupl_m, mumi_m;

    // if (Reco_mu_charge[irqq] == 1)
    // {
    //   mupl_m = Reco_mu_4mom_m->at(Reco_QQ_mupl_idx[irqq]);
    //   mupl_pt = Reco_mu_4mom_pt->at(Reco_QQ_mupl_idx[irqq]);
    //   mupl_phi = Reco_mu_4mom_phi->at(Reco_QQ_mupl_idx[irqq]);
    //   mupl_eta = Reco_mu_4mom_eta->at(Reco_QQ_mupl_idx[irqq]);
    // }
    // else if (Reco_mu_charge[irqq] == -1)
    // {
    //   mumi_m = Reco_mu_4mom_m->at(Reco_QQ_mumi_idx[irqq]);
    //   mumi_pt = Reco_mu_4mom_pt->at(Reco_QQ_mumi_idx[irqq]);
    //   mumi_phi = Reco_mu_4mom_phi->at(Reco_QQ_mumi_idx[irqq]);
    //   mumi_eta = Reco_mu_4mom_eta->at(Reco_QQ_mumi_idx[irqq]);
    // }
    
    // get 4-vectors
    JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
    mupl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
    mumi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);

    weight = 1.0; // reset to 1
    if (isMC)
    {
      if (Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1)
        continue;
      if (Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1)
        continue;
      weight = findNcoll(Centrality) * Gen_weight;
    }

    count++;
    // Todo:
    //  1. Which filter should we use?
    //  [Done] 2. where is tnp realted functions?
    if (isMC)
    {
      tnp_weight = 1; // reset
      tnp_trig_weight_mupl = -1; // negative one
      tnp_trig_weight_mumi = -1; // negative one
      tnp_weight = tnp_weight * tnp_weight_muid_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 0) * tnp_weight_muid_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 0); // mu id
      tnp_weight = tnp_weight * tnp_weight_trk_pbpb(mupl_Reco->Eta(), 0) * tnp_weight_trk_pbpb(mumi_Reco->Eta(), 0);                                     // inner tracker

      // Trigger part
      if (!((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) && (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))))
      {
        //         std::cout << "irqq : " << irqq << " - iev : " << iev << "\n";
        //         std::cout << "TnP ERROR !!!! ::: No matched L2 filter1 " << "\n";
        continue;
      }
      bool mupl_L2Filter = ((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))) ? true : false;
      bool mupl_L3Filter = ((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter))) ? true : false;
      bool mumi_L2Filter = ((Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))) ? true : false;
      bool mumi_L3Filter = ((Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter))) ? true : false;
      if (mupl_L2Filter == false || mumi_L2Filter == false)
      {
        std::cout << "TnP ERROR !!!! ::: No matched L2 filter2 " << "\n";
        std::cout << "\n";
      }

      bool mupl_isL2 = (mupl_L2Filter && !mupl_L3Filter) ? true : false;
      bool mupl_isL3 = (mupl_L2Filter && mupl_L3Filter) ? true : false;
      bool mumi_isL2 = (mumi_L2Filter && !mumi_L3Filter) ? true : false;
      bool mumi_isL3 = (mumi_L2Filter && mumi_L3Filter) ? true : false;
      bool SelDone = false;

      if (mupl_isL2 && mumi_isL3)
      {
        tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 2, 0);
        tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 3, 0);
        SelDone = true;
      }
      else if (mupl_isL3 && mumi_isL2)
      {
        tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 3, 0);
        tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 2, 0);
        SelDone = true;
      }
      else if (mupl_isL3 && mumi_isL3)
      {
        int t[2] = {-1, 1}; // mupl, mumi
        int l = rand() % (2);
        // pick up what will be L2
        if (t[l] == -1)
        {
          tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 2, 0);
          tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 3, 0);
        }
        else if (t[l] == 1)
        {
          tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 3, 0);
          tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 2, 0);
        }
        else
        {
          std::cout << "ERROR :: No random selection done !!!!" << "\n";
          continue;
        }
        SelDone = true;
      }
      if (SelDone == false || (tnp_trig_weight_mupl == -1 || tnp_trig_weight_mumi == -1))
      {
        std::cout << "ERROR :: No muon filter combination selected !!!!" << "\n";
        continue;
      }
      tnp_weight = tnp_weight * tnp_trig_weight_mupl * tnp_trig_weight_mumi;
      counttnp++;
    }

    // --- fill the output tree ---
    mass[nDimu] = JP_Reco->M();
    pt[nDimu]   = JP_Reco->Pt();
    phi[nDimu]  = JP_Reco->Phi();
    eta[nDimu]  = JP_Reco->Eta();
    y[nDimu]    = JP_Reco->Rapidity();

    pt1[nDimu]  = mupl_Reco->Pt();
    eta1[nDimu] = mupl_Reco->Eta();
    phi1[nDimu] = mupl_Reco->Phi();

    pt2[nDimu]  = mumi_Reco->Pt();
    eta2[nDimu] = mumi_Reco->Eta();
    phi2[nDimu] = mumi_Reco->Phi();

    ctau3D[nDimu]     = Reco_QQ_ctau3D[irqq];
    ctau3DErr[nDimu]  = Reco_QQ_ctauErr3D[irqq];
    ctau3DRes[nDimu]  = (Reco_QQ_ctauErr3D[irqq] > 0) ? Reco_QQ_ctau3D[irqq] / Reco_QQ_ctauErr3D[irqq] : -999;
    // ctau3D2S[nDimu]   = ctau3D[nDimu] * (pdgMass.Psi2S / pdgMass.JPsi);
    // ctau3DErr2S[nDimu] = ctau3DErr[nDimu] * (pdgMass.Psi2S / pdgMass.JPsi);
    // ctau3DRes2S[nDimu] = ctau3DRes[nDimu] * (pdgMass.Psi2S / pdgMass.JPsi);

    recoQQsign[nDimu] = Reco_QQ_sign[irqq];
    if (isMC)
      TnPweight[nDimu] = tnp_weight;

    nDimu++;
  } // end of dimuon loop

  if (nDimu > 0)
    outputTree->Fill();
}

void FlowSkimRun3DataPbPb::saveResults()
{
  std::cout << "===== saveResults() =====\n\n";
  outputFile->cd();
  outputTree->Write();
  outputFile->Close();
}

// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
//

