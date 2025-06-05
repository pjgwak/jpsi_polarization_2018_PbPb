# pragma once

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

class Data {
  public:
    /**
     * @brief constructor and destructor
     */
    Data(TTree* tree, bool isMC);
    virtual ~Data() = default;

    /**
     * @brief is sample MC?
     */
    bool isMC;

    /**
     * @brief tree variables
     */
    UInt_t runNb, eventNb, LS;
    Int_t Centrality;
    Float_t zVtx;
    TClonesArray *Reco_QQ_4mom = nullptr;
    TClonesArray *Reco_mu_4mom = nullptr;
    Float_t Reco_QQ_VtxProb[1000];
    Float_t Reco_mu_dxy[1000];
    Float_t Reco_mu_dz[1000];
    Int_t Reco_mu_nTrkWMea[1000];
    Int_t Reco_mu_nPixWMea[1000];
    Int_t Reco_mu_SelectionType[1000];
    Float_t Reco_QQ_ctau3D[1000];
    Float_t Reco_QQ_ctauErr3D[1000];

    /**
     * @brief variables to use in algorithm. Not stored in TTree.
     * TLoerntzVector pointers are non-owning pointers to objects within TClonesArrays
     * We must not delete these pointers manually. (automatic handling now)
     */
    TLorentzVector *JP_Reco = nullptr;
    TLorentzVector *mupl_Reco = nullptr;
    TLorentzVector *mumi_Reco = nullptr;

  protected:
    /**
     * @brief pointer to the TTree (or TChain) class
     */
    TTree* m_tree = nullptr;
};


// ====== Childern classes ===== //
// I wrote them here due to their small sizes
class DataRun2 : public Data
{
public:
  /**
   * @brief constructor and destructor
   */
  DataRun2(TTree *tree, bool isMC);
  ~DataRun2() override = default;

  /**
   * @brief run2 data use Int_t (Run3 uses Short_t)
   */
  Int_t Reco_QQ_size;
  Int_t Reco_QQ_mupl_idx[1000];
  Int_t Reco_QQ_mumi_idx[1000];
  Int_t Reco_QQ_sign[1000];
};

class DataRun3 : public Data
{
public:
  /**
   * @brief constructor and destructor
   */
  DataRun3(TTree *tree, bool isMC);
  ~DataRun3() override = default;

  /**
   * @brief run3 data use Short_t
   */
  Short_t Reco_QQ_size;
  Short_t Reco_QQ_mupl_idx[1000];
  Short_t Reco_QQ_mumi_idx[1000];
  Short_t Reco_QQ_sign[1000];
};