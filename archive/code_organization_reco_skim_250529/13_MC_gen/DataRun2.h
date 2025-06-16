#pragma once
#include "Data.h"

class DataRun2 : public Data
{
public:
  DataRun2(TTree *tree, bool isMC, bool isGenOnly, bool isEP, bool isRun2);
  ~DataRun2() override = default;

  // Run2 Reco specific variables - they're named: orignal_branch_name + _run2
  Int_t Gen_QQ_size_run2;
  Int_t Gen_QQ_mupl_idx_run2[1500];
  Int_t Gen_QQ_mumi_idx_run2[1500];
  Int_t Gen_mu_charge_run2[1500];
  
  Int_t Reco_QQ_size_run2;
  Int_t Reco_mu_whichGen_run2[1500];
  Int_t Reco_QQ_mupl_idx_run2[1500];
  Int_t Reco_QQ_mumi_idx_run2[1500];
  Int_t Reco_QQ_sign_run2[1500];

  // override getters
  // be careful - we never return the whole array or its pointer, instead we return the value!
  Long64_t getGenQQSize() const override { return Gen_QQ_size_run2; }
  Int_t getGenQQMuplIdx(Long64_t idx) const override { return Gen_QQ_mupl_idx_run2[idx]; }
  Int_t getGenQQMumiIdx(Long64_t idx) const override { return Gen_QQ_mumi_idx_run2[idx]; }
  Int_t getGenMuCharge(Long64_t idx) const override { return Gen_mu_charge_run2[idx]; }

  Long64_t getRecoQQSize() const override { return Reco_QQ_size_run2; }
  Int_t getRecoQQSign(Long64_t idx) const override { return Reco_QQ_sign_run2[idx]; }
  Int_t getRecoQQMuplIdx(Long64_t idx) const override { return Reco_QQ_mupl_idx_run2[idx]; }
  Int_t getRecoQQMumiIdx(Long64_t idx) const override { return Reco_QQ_mumi_idx_run2[idx]; }
  Int_t getRecoMuWhichGen(Long64_t idx) const override { return Reco_mu_whichGen_run2[idx]; }
};
