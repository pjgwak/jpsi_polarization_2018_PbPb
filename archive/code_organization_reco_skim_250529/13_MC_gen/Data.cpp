#include "Data.h"

Data::Data(TTree *tree, bool isMC, bool isGenOnly, bool isEP, bool isRun2)
    : isMC_flag(isMC), isGenOnly_flag(isGenOnly), isEP_flag(isEP), isRun2_flag(isRun2), m_tree(tree) {
  // Basic connecting grammer - m_tree->SetBranchAddress("XX", &XX);
  // TClonesArray pointers - m_tree->SetBranchAddress("XX", &XX);
  // fixed-size array (no &) - m_tree->SetBranchAddress("XX", XX);

  // common branches - check existences of branches first.
  if (!isMC_flag) {
    // data only branches
    if (m_tree->GetBranch("runNb")) m_tree->SetBranchAddress("runNb", &runNb);
    if (m_tree->GetBranch("LS")) m_tree->SetBranchAddress("LS", &LS);
  }

  if (m_tree->GetBranch("eventNb")) m_tree->SetBranchAddress("eventNb", &eventNb);
  if (m_tree->GetBranch("Centrality")) m_tree->SetBranchAddress("Centrality", &Centrality);
  if (m_tree->GetBranch("zVtx")) m_tree->SetBranchAddress("zVtx", &zVtx);
  if (m_tree->GetBranch("HLTriggers")) m_tree->SetBranchAddress("HLTriggers", &HLTriggers);
  if (m_tree->GetBranch("SumET_HF")) m_tree->SetBranchAddress("SumET_HF", &SumET_HF);
  
  if (m_tree->GetBranch("Reco_QQ_4mom")) m_tree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom);
  if (m_tree->GetBranch("Reco_mu_4mom")) m_tree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom);
  if (m_tree->GetBranch("Reco_QQ_VtxProb")) m_tree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb);
  if (m_tree->GetBranch("Reco_mu_dxy")) m_tree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy);
  if (m_tree->GetBranch("Reco_mu_dz")) m_tree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz);
  if (m_tree->GetBranch("Reco_mu_nTrkWMea")) m_tree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea);
  if (m_tree->GetBranch("Reco_mu_nPixWMea")) m_tree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea);
  if (m_tree->GetBranch("Reco_mu_SelectionType")) m_tree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType);
  if (m_tree->GetBranch("Reco_QQ_ctau3D")) m_tree->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D);
  if (m_tree->GetBranch("Reco_QQ_ctauErr3D")) m_tree->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D);
  if (m_tree->GetBranch("Reco_mu_trig")) m_tree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig);
  if (m_tree->GetBranch("Reco_QQ_trig")) m_tree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig);

  // MC only branches
  if (isMC_flag) {
    if (m_tree->GetBranch("Gen_weight")) m_tree->SetBranchAddress("Gen_weight", &Gen_weight);
    if (m_tree->GetBranch("Gen_QQ_4mom")) m_tree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);
    if (m_tree->GetBranch("Gen_mu_4mom")) m_tree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom);
    if (m_tree->GetBranch("Gen_QQ_ctau3D")) m_tree->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D);

    // GenOnly branches
    if (isGenOnly_flag) {
      if (m_tree->GetBranch("Gen_QQ_mupl_4mom")) m_tree->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom);
      if (m_tree->GetBranch("Gen_QQ_mumi_4mom")) m_tree->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom);
      if (m_tree->GetBranch("Gen_QQ_ctau")) m_tree->SetBranchAddress("Gen_QQ_ctau", Gen_QQ_ctau);
    }
  }

  // HI EP angles
  if (isEP_flag)
    m_tree->SetBranchAddress("ep_flat", &ep_flat);
}
