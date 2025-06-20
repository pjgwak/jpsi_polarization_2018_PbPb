# pragma once

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

class Data {
  public:
    /**
     * @brief construct a new Data object
     */
    Data(TTree* tree);

    /**
     * @brief tree variables
     */
    UInt_t runNb;
    Int_t Centrality;
    Short_t Reco_QQ_size;
    TClonesArray *Reco_QQ_4mom = nullptr;
    TClonesArray *Reco_mu_4mom = nullptr;
    Float_t Reco_QQ_VtxProb[1000];
    Short_t Reco_QQ_mupl_idx[1000];
    Short_t Reco_QQ_mumi_idx[1000];
    Float_t Reco_mu_dxy[1000];
    Float_t Reco_mu_dz[1000];
    Int_t Reco_mu_nTrkWMea[1000];
    Int_t Reco_mu_nPixWMea[1000];
    Short_t Reco_QQ_sign[1000];
    Int_t Reco_mu_SelectionType[1000];

    /**
     * variables to use in algorithm. Not stored in TTree.
     */
    TLorentzVector *JP_Reco = nullptr;
    TLorentzVector *mupl_Reco = nullptr;
    TLorentzVector *mumi_Reco = nullptr;

  protected:
    /**
     * @brief pointer to the TTree (or TChain) class
     */
    TTree* m_tree = nullptr;

    

    /**
     * @brief
     */

    /**
     * @brief
     */

    /**
     * @brief
     */


};