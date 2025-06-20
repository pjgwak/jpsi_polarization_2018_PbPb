# pragma once

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include <iostream>

class Data {
  public:
    /**
     * @brief constructor and destructor
     */
    Data(TTree *tree, bool isMC, bool isGenOnly, bool isEP, bool isRun2);
    virtual ~Data() = default;

    // common flags
    bool isMC_flag;
    bool isGenOnly_flag;
    bool isEP_flag;
    bool isRun2_flag;

    // Reco common branches - All samples have it.
    UInt_t runNb, eventNb, LS;
    Int_t Centrality;
    ULong64_t HLTriggers;
    Float_t zVtx, Gen_weight, SumET_HF;
    Float_t ep_flat;
    TClonesArray *Reco_QQ_4mom = nullptr;
    TClonesArray *Reco_mu_4mom = nullptr;
    Float_t Reco_QQ_VtxProb[1500];
    Float_t Reco_mu_dxy[1500];
    Float_t Reco_mu_dz[1500];
    Int_t Reco_mu_nTrkWMea[1500];
    Int_t Reco_mu_nPixWMea[1500];
    Int_t Reco_mu_SelectionType[1500];
    Float_t Reco_QQ_ctau3D[1500];
    Float_t Reco_QQ_ctauErr3D[1500];
    ULong64_t Reco_mu_trig[1500];
    ULong64_t Reco_QQ_trig[1500];

    // Gen common variables
    TClonesArray *Gen_QQ_4mom = nullptr;
    TClonesArray *Gen_mu_4mom = nullptr;
    Float_t Gen_QQ_ctau3D[1500];

    // Run2GenOnly branches
    TClonesArray *Gen_QQ_mupl_4mom = nullptr;
    TClonesArray *Gen_QQ_mumi_4mom = nullptr;
    Float_t Gen_QQ_ctau[1500];

    /**
     * @brief variables to use in algorithm. Not stored in TTree.
     * TLoerntzVector pointers are non-owning pointers to objects within TClonesArrays
     * We must not delete these pointers manually. (automatic handling now)
     */
    TLorentzVector *JP_Reco = nullptr;
    TLorentzVector *mupl_Reco = nullptr;
    TLorentzVector *mumi_Reco = nullptr;
    TLorentzVector *JP_Gen = nullptr;
    TLorentzVector *mupl_Gen = nullptr;
    TLorentzVector *mumi_Gen = nullptr;

    // virtual getter functions - Samples have different branch types
    // = 0 is also the grammer.
    virtual Long64_t getGenQQSize() const = 0;
    virtual Int_t getGenQQMuplIdx(Long64_t idx) const = 0;
    virtual Int_t getGenQQMumiIdx(Long64_t idx) const = 0;
    virtual Int_t getGenMuCharge(Long64_t idx) const = 0;
    
    virtual Long64_t getRecoQQSize() const = 0;
    virtual Int_t getRecoQQSign(Long64_t idx) const = 0;
    virtual Int_t getRecoQQMuplIdx(Long64_t idx) const = 0;
    virtual Int_t getRecoQQMumiIdx(Long64_t idx) const = 0;
    virtual Int_t getRecoMuWhichGen(Long64_t idx) const = 0;

  protected:
    // TTree (or TChain) class
    TTree* m_tree = nullptr;
};
