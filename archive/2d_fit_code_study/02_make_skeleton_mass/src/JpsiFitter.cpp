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

void JpsiFitter::processTree(const std::string &filePath, const std::string &dsname, const std::string &obs, const std::string &cut){
  auto inputFile = std::unique_ptr<TFile>(TFile::Open(filePath.c_str()));
  inputFile->SetName("inputFile");
  // inputFile.Print("V");

  RooDataSet *ds_temp = (RooDataSet *)inputFile->Get("dataset");
  ds_temp->SetName(dsname.c_str());
  m_ws.import(*ds_temp);
}

void JpsiFitter::print()
{
  // show the contents of workspace
  m_ws.Print("V");
}