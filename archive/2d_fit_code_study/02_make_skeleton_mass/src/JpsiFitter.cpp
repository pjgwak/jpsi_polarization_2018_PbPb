#include "JpsiFitter.h"
#include "TFile.h"
#include "TTree.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgSet.h"

// RDataFrame
#include <ROOT/RDataFrame.hxx>
#include <RooAbsDataHelper.h>
#include <ROOT/RSnapshotOptions.hxx>
#include <iostream>

JpsiFitter::JpsiFitter(RooWorkspace &ws) : m_ws(ws) {}
JpsiFitter::~JpsiFitter() {}

void JpsiFitter::processTree(const std::string &filePath, const std::string &treeName, const std::string &obs, const std::string &cut){
  auto f_pr_mc = std::unique_ptr<TFile>(TFile::Open("/work/pjgwak/pol24/files_roodata/RooDataSet_miniAOD_isMC1_NP_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root"));
  f_pr_mc->SetName("f_pr_mc");
  // f_pr_mc.Print("V");
  RooDataSet *ds_mc_temp = (RooDataSet *)f_pr_mc->Get("dataset");
  ds_mc_temp->SetName("ds_mc_temp");
  m_ws.import(*ds_mc_temp);
}

void JpsiFitter::print()
{
  // show the contents of workspace
  m_ws.Print("V");
}