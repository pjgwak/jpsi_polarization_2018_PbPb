#pragma once
#include <string>
#include <vector>
#include "Rtypes.h" // UInt_t, etc

class TTree;
class TChain;
class TFile;
class TLorentzVector;
class TClonesArray;

class FlowSkimRun3DataPbPb
{
  public:
    FlowSkimRun3DataPbPb(bool isMC, int kTrigSel, int hiHFBinEdge, int PDtype,
      std::string inputFilePath, std::string inputTreeName, std::string outputFilePath,
      std::string outputTreeName, std::string userTag);
    ~FlowSkimRun3DataPbPb();

    void run(long nevt = -1);

    // in(out)put paths
    std::string inputFilePath, inputTreeName;
    std::string outputFilePath, outputTreeName;
    std::string userTag;

  private:
    void setLabels();
    void openInputFile();
    void setInputBranches();
    void setOutputFile();
    void setOutputBranches();
    void processEvent(long iev);
    void saveResults();

    bool isMC;
    int kTrigSel, hiHFBinEdge, PDtype;
    std::string sampleType, centHfTag;

    // triger label
    int trigIndx = 0;

    // triger and filters
    int kTrigJpsi = 12;
    int kTrigJpsipp = 3;
    int kL2filter = 16; // Todo: should be checked. -> PbPb Jpsi filter. Only for MC
    int kL3filter = 17;

    TChain *inputTree = nullptr;
    TTree *outputTree = nullptr;
    TFile *outputFile = nullptr;

    // event counter
    long totalEvents = 0;
    int count = 0;
    int count_soft = 0;
    int count_vtx = 0;
    int counttnp = 0;
    // int passedTrigger = 0;
    // int passedMuonSel = 0;
    // int passedVertex = 0;
    // int passedAll = 0;

    // tnp weight
    double tnp_weight = 1;
    double tnp_trig_weight_mupl = -1; // negative one
    double tnp_trig_weight_mumi = -1; // negative one

    static constexpr int maxBranchSize = 10000; // max input items
    static constexpr int nMaxDimu = 10000;      // max output items

    // --- input tree variabels ---
    // scalar values
    UInt_t runNb, eventNb, LS;
    Float_t zVtx;
    Int_t Centrality;
    ULong64_t HLTriggers;
    Float_t SumET_HF;
    Short_t Reco_QQ_size, Reco_mu_size;
    Float_t Gen_weight;

    // Reco_QQ 
    Short_t Reco_QQ_mupl_idx[maxBranchSize], Reco_QQ_mumi_idx[maxBranchSize], Reco_QQ_sign[maxBranchSize];
    ULong64_t Reco_QQ_trig[maxBranchSize];
    Float_t Reco_QQ_VtxProb[maxBranchSize], Reco_QQ_ctau3D[maxBranchSize], Reco_QQ_ctauErr3D[maxBranchSize];

    // Reco_mu 
    Int_t Reco_mu_whichGen[maxBranchSize];
    Bool_t Reco_mu_highPurity[maxBranchSize], Reco_mu_TMOneStaTight[maxBranchSize];
    Int_t Reco_mu_nTrkHits[maxBranchSize], Reco_mu_nMuValHits[maxBranchSize], Reco_mu_StationsMatched[maxBranchSize];
    Float_t Reco_mu_normChi2_global[maxBranchSize], Reco_mu_dxy[maxBranchSize], Reco_mu_dxyErr[maxBranchSize];
    Float_t Reco_mu_dz[maxBranchSize], Reco_mu_dzErr[maxBranchSize];
    Int_t Reco_mu_nTrkWMea[maxBranchSize], Reco_mu_nPixWMea[maxBranchSize], Reco_mu_nPixValHits[maxBranchSize];
    Int_t Reco_mu_SelectionType[maxBranchSize];
    Short_t Reco_mu_type[maxBranchSize], Reco_mu_charge[maxBranchSize];
    ULong64_t Reco_mu_trig[maxBranchSize];

    // TLorentzVector - must be nullptr
    TClonesArray *Reco_QQ_4mom = nullptr, *Reco_mu_4mom = nullptr;

    // --- output tree variables ---
    int evt, runN, lumi, cBin, nDimu;
    float vz;
    float mass[nMaxDimu], pt[nMaxDimu], y[nMaxDimu], phi[nMaxDimu], eta[nMaxDimu];
    float eta1[nMaxDimu], eta2[nMaxDimu], phi1[nMaxDimu], phi2[nMaxDimu], pt1[nMaxDimu], pt2[nMaxDimu];
    float ctau3D[nMaxDimu], ctau3DErr[nMaxDimu], ctau3DRes[nMaxDimu];
    float ctau3D2S[nMaxDimu], ctau3DErr2S[nMaxDimu], ctau3DRes2S[nMaxDimu];
    int recoQQsign[nMaxDimu];
    double TnPweight[nMaxDimu], weight = 1.0;

    TLorentzVector *JP_Reco = nullptr;
    TLorentzVector *mupl_Reco = nullptr;
    TLorentzVector *mumi_Reco = nullptr;
};