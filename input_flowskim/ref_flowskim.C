#include "FlowSkimmer.h"
#include <cmath>
#include <iostream>
#include <fstream>

FlowSkimmer::FlowSkimmer(bool isMC_, int kTrigSel_, int hiHFBinEdge_, int PDtype_)
    : isMC(isMC_), kTrigSel(kTrigSel_), hiHFBinEdge(hiHFBinEdge_), PDtype(PDtype_) {}

FlowSkimmer::~FlowSkimmer()
{
  if (outputFile)
  {
    outputFile->Write();
    outputFile->Close();
  }
  if (exceptionLog.is_open())
    exceptionLog.close();
}

void FlowSkimmer::run(int nevt)
{
  openInputFile();
  setupInputBranches();
  setupOutputFile();
  setupOutputBranches();

  totalEvents = (nevt < 0) ? inputTree->GetEntries() : nevt;

  for (int iev = 0; iev < totalEvents; ++iev)
  {
    processEvent(iev);
  }

  // --- Selection 요약 출력 ---
  std::cout << "\n========== Selection Summary ==========\n";
  std::cout << "Total Events       : " << totalEvents << "\n";
  std::cout << "Trigger Passed     : " << passedTrigger << "\n";
  std::cout << "MuonSel Passed     : " << passedMuonSel << "\n";
  std::cout << "Vertex Passed      : " << passedVertex << "\n";
  std::cout << "Fully Selected     : " << passedAll << "\n";
}

void FlowSkimmer::openInputFile()
{
  inputTree = new TChain("mytree");
  inputTree->Add("your_input_file.root");
}

void FlowSkimmer::setupInputBranches()
{
  inputTree->SetBranchAddress("runNb", &runNb);
  inputTree->SetBranchAddress("eventNb", &eventNb);
  inputTree->SetBranchAddress("LS", &LS);
  inputTree->SetBranchAddress("zVtx", &zVtx);
  inputTree->SetBranchAddress("Centrality", &Centrality);
  inputTree->SetBranchAddress("HLTriggers", &HLTriggers);
  inputTree->SetBranchAddress("SumET_HF", &SumET_HF);
  inputTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
  inputTree->SetBranchAddress("Reco_mu_size", &Reco_mu_size);
  inputTree->SetBranchAddress("Reco_mu_charge", Reco_mu_charge);
  inputTree->SetBranchAddress("Reco_QQ_4mom_m", &Reco_QQ_4mom_m);
  inputTree->SetBranchAddress("Reco_QQ_4mom_eta", &Reco_QQ_4mom_eta);
  inputTree->SetBranchAddress("Reco_QQ_4mom_phi", &Reco_QQ_4mom_phi);
  inputTree->SetBranchAddress("Reco_QQ_4mom_pt", &Reco_QQ_4mom_pt);
  inputTree->SetBranchAddress("Reco_mu_4mom_m", &Reco_mu_4mom_m);
  inputTree->SetBranchAddress("Reco_mu_4mom_eta", &Reco_mu_4mom_eta);
  inputTree->SetBranchAddress("Reco_mu_4mom_phi", &Reco_mu_4mom_phi);
  inputTree->SetBranchAddress("Reco_mu_4mom_pt", &Reco_mu_4mom_pt);
  inputTree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig);
  inputTree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig);
  inputTree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb);
  inputTree->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity);
  inputTree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
  inputTree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);
  inputTree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits);
  inputTree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global);
  inputTree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits);
  inputTree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched);
  inputTree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy);
  inputTree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr);
  inputTree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz);
  inputTree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr);
  inputTree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea);
  inputTree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea);
  inputTree->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight);
  inputTree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign);
  inputTree->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits);
  inputTree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType);
  inputTree->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D);
  inputTree->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D);
  if (isMC)
  {
    inputTree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen);
    inputTree->SetBranchAddress("Gen_weight", &Gen_weight);
  }
}

void FlowSkimmer::setupOutputFile()
{
  outputFile = new TFile("flowSkim_output.root", "RECREATE");
  outputTree = new TTree("mmepevt", "output tree");
}

void FlowSkimmer::setupOutputBranches()
{
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
  outputTree->Branch("ctau3D", ctau3D, "ctau3D[nDimu]/F");
  outputTree->Branch("ctau3DErr", ctau3DErr, "ctau3DErr[nDimu]/F");
  outputTree->Branch("ctau3DRes", ctau3DRes, "ctau3DRes[nDimu]/F");
  outputTree->Branch("ctau3D2S", ctau3D2S, "ctau3D2S[nDimu]/F");
  outputTree->Branch("ctau3DErr2S", ctau3DErr2S, "ctau3DErr2S[nDimu]/F");
  outputTree->Branch("ctau3DRes2S", ctau3DRes2S, "ctau3DRes2S[nDimu]/F");
  outputTree->Branch("recoQQsign", recoQQsign, "recoQQsign[nDimu]/I");
  outputTree->Branch("weight", &weight, "weight/D");
  outputTree->Branch("TnPweight", TnPweight, "TnPweight[nDimu]/D");
}

void FlowSkimmer::processEvent(int iev)
{
  inputTree->GetEntry(iev);

  // 기본 이벤트 정보 저장
  evt = eventNb;
  runN = runNb;
  lumi = LS;
  vz = zVtx;

  // centrality bin 계산
  if (hiHFBinEdge == 0)
    cBin = getHiBinFromhiHF(SumET_HF); // 일반 HF 기준
  else if (hiHFBinEdge == 1)
    cBin = getHiBinFromhiHF_Up(SumET_HF); // HF 업쉬프트
  else if (hiHFBinEdge == -1)
    cBin = getHiBinFromhiHF_Down(SumET_HF); // HF 다운쉬프트
  else
    cBin = -999;

  if (cBin == -999)
    return; // 에러 시 무시

  nDimu = 0; // 이번 이벤트의 디뮤온 수

  // 트리거 조건 (HLT 비트)
  if ((HLTriggers & (1ULL << kTrigSel)) == 0)
    return;

  for (int iQQ = 0; iQQ < Reco_QQ_size; ++iQQ)
  {
    // 디뮤온에 대한 트리거 패턴 검사
    if ((Reco_QQ_trig[iQQ] & (1ULL << kTrigSel)) == 0)
      continue;

    int idx_mupl = Reco_QQ_mupl_idx[iQQ];
    int idx_mumi = Reco_QQ_mumi_idx[iQQ];
    if (idx_mupl < 0 || idx_mumi < 0)
      continue;

    // Muon Soft Selection 조건
    auto isSoftMuon = [&](int idx) -> bool
    {
      return (Reco_mu_nTrkWMea[idx] > 5) &&
             (Reco_mu_nPixWMea[idx] > 0) &&
             (std::fabs(Reco_mu_dxy[idx]) < 0.3) &&
             (std::fabs(Reco_mu_dz[idx]) < 20.) &&
             (Reco_mu_SelectionType[idx] & (1 << 1)) && // tracker muon
             (Reco_mu_SelectionType[idx] & (1 << 3));   // global muon
    };

    if (!(isSoftMuon(idx_mupl) && isSoftMuon(idx_mumi)))
      continue;

    // Vertex 조건
    if (Reco_QQ_VtxProb[iQQ] < 0.01)
      continue;

    // TLorentzVector 재구성
    TLorentzVector mupl, mumi, dimu;
    mupl.SetPtEtaPhiM(Reco_mu_4mom_pt->at(idx_mupl), Reco_mu_4mom_eta->at(idx_mupl),
                      Reco_mu_4mom_phi->at(idx_mupl), Reco_mu_4mom_m->at(idx_mupl));
    mumi.SetPtEtaPhiM(Reco_mu_4mom_pt->at(idx_mumi), Reco_mu_4mom_eta->at(idx_mumi),
                      Reco_mu_4mom_phi->at(idx_mumi), Reco_mu_4mom_m->at(idx_mumi));
    dimu.SetPtEtaPhiM(Reco_QQ_4mom_pt->at(iQQ), Reco_QQ_4mom_eta->at(iQQ),
                      Reco_QQ_4mom_phi->at(iQQ), Reco_QQ_4mom_m->at(iQQ));

    // MC weight 적용
    weight = 1.0;
    if (isMC)
    {
      if (Reco_mu_whichGen[idx_mupl] < 0 || Reco_mu_whichGen[idx_mumi] < 0)
        continue;
      weight = findNcoll(Centrality) * Gen_weight;
    }

    // TnP weight 적용
    TnPweight[nDimu] = 1.0;
    if (isMC)
    {
      double w_mupl_id = tnp_weight_muid_pbpb(mupl.Pt(), mupl.Eta(), 0);
      double w_mumi_id = tnp_weight_muid_pbpb(mumi.Pt(), mumi.Eta(), 0);
      double w_mupl_trk = tnp_weight_trk_pbpb(mupl.Eta(), 0);
      double w_mumi_trk = tnp_weight_trk_pbpb(mumi.Eta(), 0);

      // 트리거 필터 weight은 아래와 같이 예시 처리
      double w_mupl_trg = 1.0;
      double w_mumi_trg = 1.0;
      TnPweight[nDimu] = w_mupl_id * w_mumi_id * w_mupl_trk * w_mumi_trk * w_mupl_trg * w_mumi_trg;
    }

    // 출력 배열에 채우기
    mass[nDimu] = dimu.M();
    pt[nDimu] = dimu.Pt();
    eta[nDimu] = dimu.Eta();
    phi[nDimu] = dimu.Phi();
    y[nDimu] = 0.5 * log((dimu.E() + dimu.Pz()) / (dimu.E() - dimu.Pz()));

    pt1[nDimu] = mupl.Pt();
    pt2[nDimu] = mumi.Pt();
    eta1[nDimu] = mupl.Eta();
    eta2[nDimu] = mumi.Eta();
    phi1[nDimu] = mupl.Phi();
    phi2[nDimu] = mumi.Phi();

    recoQQsign[nDimu] = Reco_QQ_sign[iQQ];

    ctau3D[nDimu] = Reco_QQ_ctau3D[iQQ];
    ctau3DErr[nDimu] = Reco_QQ_ctauErr3D[iQQ];
    ctau3DRes[nDimu] = ctau3D[nDimu] / ctau3DErr[nDimu];

    constexpr double mJPsi = 3.0969; // PDG mass
    constexpr double mPsi2S = 3.6861;

    ctau3D2S[nDimu] = ctau3D[nDimu] * (mPsi2S / mJPsi);
    ctau3DErr2S[nDimu] = ctau3DErr[nDimu] * (mPsi2S / mJPsi);
    ctau3DRes2S[nDimu] = ctau3DRes[nDimu] * (mPsi2S / mJPsi);

    ++nDimu;
  }

  if (nDimu > 0)
    outputTree->Fill();
}