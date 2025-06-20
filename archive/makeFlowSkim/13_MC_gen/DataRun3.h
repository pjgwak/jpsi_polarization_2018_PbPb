#pragma once
#include "Data.h"

class DataRun3 : public Data
{
public:
  DataRun3(TTree *tree, bool isMC, bool isGenOnly, bool isEP, bool isRun2);
  ~DataRun3() override = default;

  // run3 specific branches
  Short_t Gen_QQ_size_run3;
  Short_t Gen_QQ_mupl_idx_run3[1500];
  Short_t Gen_QQ_mumi_idx_run3[1500];
  Short_t Gen_mu_charge_run3[1500];

  Short_t Reco_QQ_size_run3;
  Short_t Reco_mu_whichGen_run3[1500];
  Short_t Reco_QQ_mupl_idx_run3[1500];
  Short_t Reco_QQ_mumi_idx_run3[1500];
  Short_t Reco_QQ_sign_run3[1500];

  // override getters - explict type casting from Short_t to Int_t
  Long64_t getGenQQSize() const override { return Gen_QQ_size_run3; }
  Int_t getGenQQMuplIdx(Long64_t idx) const override { return static_cast<Int_t>(Gen_QQ_mupl_idx_run3[idx]); }
  Int_t getGenQQMumiIdx(Long64_t idx) const override { return static_cast<Int_t>(Gen_QQ_mumi_idx_run3[idx]); }
  Int_t getGenMuCharge(Long64_t idx) const override { return static_cast<Int_t>(Gen_mu_charge_run3[idx]); }

  Long64_t getRecoQQSize() const override { return Reco_QQ_size_run3; }
  Int_t getRecoQQSign(Long64_t idx) const override { return static_cast<Int_t>(Reco_QQ_sign_run3[idx]); }
  Int_t getRecoQQMuplIdx(Long64_t idx) const override { return static_cast<Int_t>(Reco_QQ_mupl_idx_run3[idx]); }
  Int_t getRecoQQMumiIdx(Long64_t idx) const override { return static_cast<Int_t>(Reco_QQ_mumi_idx_run3[idx]); }
  Int_t getRecoMuWhichGen(Long64_t idx) const override { return static_cast<Int_t>(Reco_mu_whichGen_run3[idx]); }
};