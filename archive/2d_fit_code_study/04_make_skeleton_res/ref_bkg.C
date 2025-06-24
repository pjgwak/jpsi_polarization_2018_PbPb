#include <iostream>
#include <utility> // std::pair
#include "RooMsgService.h"


// #include <RooFormulaVar.h>
// #include <RooCBShape.h>
// #include <RooWorkspace.h>
// #include <RooChebychev.h>
// #include <RooPolynomial.h>
// #include "RooPlot.h"
// #include "TText.h"
// #include "TArrow.h"
// #include "TFile.h"
// #include "RooDataHist.h"
// #include "RooCategory.h"
// #include "RooSimultaneous.h"
// #include "RooStats/SPlot.h"
// #include "../headers/rootFitHeaders.h"
// #include "../headers/commonUtility.h"
// #include "../headers/JpsiUtility.h"
// #include "../headers/cutsAndBin.h"
// #include "../headers/CMS_lumi_v2mass.C"
// #include "../headers/tdrstyle.C"

// global variable - will be memober variable
const double nCtauBins = 200;
const float ctauLow = -4.;
const float ctauHigh = 7.;

// void setConfig();
void setupWorkspace(RooWorkspace *ws, const std::string &mass_path, const std::string &err_path);
std::pair<double, double> processDataset(RooWorkspace *ws, bool isForceMin);
// void buildBkgModel(RooWorkspace *ws, int nGauss);
// RooFitResult *perfromBkgFit(RooWorkspace *ws);
// void drawBkgPlot(RooWorkspace *ws, double ctauResMin, double ctauResMax, int nGauss, const std::string &outpath); // Todo: use outPath
// void saveOutput(RooWorkspace *ws, const std::string &outpath);


//
//
//

using namespace std;
using namespace RooFit;

std::pair<double, double> processDataset(RooWorkspace *ws, bool isForceMin)
{
  std::cout << "===== Start processDataset() =====\n";
  
  // ===== set ctau3DErr range =====
  double ctauErrMin, ctauErrMax;
  ws->data("ds_errBkg")->getRange(*ws->var("ctau3DErr"), ctauErrMin, ctauErrMax);
  ws->var("ctau3D")->setRange(ctauLow, ctauHigh);
  ws->var("ctau3D")->setRange("ctauRange", ctauLow, ctauHigh);
  ws->var("ctau3DErr")->setRange(ctauErrMin, ctauErrMax);
  ws->var("ctau3DErr")->setRange("ctauErrRange", ctauErrMin, ctauErrMax);
  // ws->var("ctau3D")->Print();
  // ws->var("ctau3DErr")->Print();

  // ===== et ctau3D range =====
  TH1D *h_tmp = (TH1D *)ws->data("ds_errBkg")->createHistogram(("h_tmp"), *ws->var("ctau3D"), Binning(nCtauBins, ctauLow, ctauHigh));
  double ctauMin = h_tmp->GetBinLowEdge(h_tmp->FindFirstBinAbove(1, 1));
  double ctauMax = h_tmp->GetBinLowEdge(h_tmp->FindLastBinAbove(2, 1)) + h_tmp->GetBinWidth(h_tmp->FindLastBinAbove(2, 1));

  if (isForceMin) {
    // ctauLow, High: default, ctauMin, Max -> user custom
    // skip: user custom ctauMin - make later
    // if (cLow == 80 && cHigh == 100)
    //     ctauMin = -1.0;
    // else if (cLow == 100 && cHigh == 180)
    //     ctauMin = -0.6;
    ctauMin = -1.0;
  }
  ws->var("ctau3D")->setRange("ctauFitRange", ctauLow, ctauHigh);

  // choose one of them
  RooDataSet* dataToFit = (RooDataSet*)ds_errBkg->reduce(Form("ctau3D>=%.f&&ctau3D<=%.f",ctauMin, ctauMax))->Clone("ds_errBkg");
  // RooDataSet *dataToFit = (RooDataSet *)ds_errBkg->reduce(Form("ctau3D>=%.f&&ctau3D<=%.f", ctauLow, ctauHigh))->Clone("dataToFit");
  ws->import(*dataToFit);

  std::cout << "===== Finish processDataset() =====\n";
}



void ref_bkg()
{
  std::cout << "===== Start ref_bkg.C =====\n";
  TStopwatch *t = new TStopwatch;
  t->Start();

  // get silence
  RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling); // ingnore import success log
  // RooMsgService::instance().setGlobalKillBelow(WARNING);

  gStyle->SetEndErrorSize(0); // error bar will be a long bar shape -> to avoid overlapping

  // void setConfig(); -> skip: will be  member method

  // setup workspace
  RooWorkspace *ws_bkg = new RooWorkspace("ws_bkg", "");
  setupWorkspace(ws_bkg, "roots/mass.root", "roots/splot.root", "roots/res.root");

  // process dataset
  processDataset(ws_bkg, false);

  // build ctauBkg model
  // void buildBkgModel(RooWorkspace *ws, int nGauss);

  // do fit
  // RooFitResult *perfromBkgFit(RooWorkspace *ws);

  // draw plots
  // void drawBkgPlot(RooWorkspace *ws, double ctauResMin, double ctauResMax, int nGauss, const std::string &outpath);

  // draw ratio plot
  // void drawBkgRatioPlot(RooWorkspace *ws, double ctauResMin, double ctauResMax, int nGauss, const std::string &outpath)

  // save output
  // void saveOutput(RooWorkspace *ws, const std::string &outpath);

  t->Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t->RealTime(), t->CpuTime());
  std::cout << "===== Start ref_bkg.C =====\n";
}


void setupWorkspace(RooWorkspace *ws, const std::string &mass_path, const std::string &err_path, const std::string &res_path)
{
  std::cout << "===== start setupWorkspace() =====\n";
  // get 3 input file paths
  TFile *f_mass = new TFile(mass_path.c_str());
  TFile *f_err = new TFile(err_path.c_str());
  TFile *f_res = new TFile(res_path.c_str());

  if (!f_mass || !f_err || !f_res)
  {
    std::cout << "Error: Could not open one of the files.\n";
    return;
  }

  // 하나씩 import 테스트하기!

  // import mass dataset and model
  RooWorkspace *ws_mass = (RooWorkspace *)f_mass->Get("wsMy");
  RooDataSet *ds_mass = (RooDataSet *)ws_mass->data("ds_mass");
  RooAddPdf *massModel = (RooAddPdf *)ws_mass->pdf("massModel");
  ws->import(*ds_mass);
  ws->import(*massModel);

  // splot: sWeighted bkg dataset and model
  RooWorkspace *ws_err = (RooWorkspace *)f_err->Get("wsMy");
  RooDataSet *ds_errBkg = (RooDataSet *)ws_err->data("ds_errBkg");
  RooAddPdf *errBkgModel = (RooAddPdf *)ws_err->pdf("errBkgModel");
  ws->import(*ds_errBkg);
  ws->import(*errBkgModel);

  // import ctauRes fit model
  RooWorkspace *ws_res = (RooWorkspace *)f_res->Get("ws_res");
  RooAddPdf *gaus_res_model = (RooAddPdf *)ws_res->pdf("gaus_res_model");
  ws->import(*gaus_res_model);

  // ws->Print("V");
  std::cout << "===== finish setupWorkspace() =====\n";
}

void setConfig()
{
    // // input parameters
    // float cos_low = 0.0, cos_high = 0.1;
    // float ptLow = 3;
    // float ptHigh = 6.5;
    // float yLow = 1.6;
    // float yHigh = 2.4;
    // int cLow = 60;
    // int cHigh = 180;

    // // Usually not used
    // int PR = 0; // 0=PR, 1=NP, 2=Inc.
    // int PRw = 1;
    // bool fEffW = false;
    // bool fAccW = false;
    // bool isPtW = false;
    // bool isTnP = false;
    // double massLow = 2.6;
    // double massHigh = 3.5; // Jpsi mass range

    // // be quiet please
    // // RooMsgService::instance().getStream(0).removeTopic(Caching); // in RooFit::
    // // RooMsgService::instance().getStream(1).removeTopic(Caching);
    // RooMsgService::instance().getStream(0).removeTopic(Plotting);
    // RooMsgService::instance().getStream(1).removeTopic(Plotting);
    // // RooMsgService::instance().getStream(0).removeTopic(Integration);
    // // RooMsgService::instance().getStream(1).removeTopic(Integration);
    // // RooMsgService::instance().getStream(0).removeTopic(InputArguments);
    // // RooMsgService::instance().getStream(1).removeTopic(InputArguments);
    // // RooMsgService::instance().getStream(0).removeTopic(NumIntegration);
    // // RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
    // // RooMsgService::instance().getStream(0).removeTopic(Minimization);
    // // RooMsgService::instance().getStream(1).removeTopic(Minimization);
    // RooMsgService::instance().setGlobalKillBelow(WARNING);

    // // ===== labeling ===== //
    // std::ostringstream out_file_path;
    // out_file_path << "bkg/bkg_pt" << ptLow << "-" << ptHigh
    //               << "_cent" << cLow << "-" << cHigh << "_absY" << yLow << "-" << yHigh
    //               << "_cos" << cos_low << "_" << cos_high;
    // string out_ss = out_file_path.str();

    // std::ostringstream input_file;
    // input_file << "pt" << ptLow << "-" << ptHigh
    //            << "_cent" << cLow << "-" << cHigh << "_absY" << yLow << "-" << yHigh
    //            << "_cos" << cos_low << "_" << cos_high
    //            << ".root";
    // string in_ss = input_file.str();

    // // ===== make output folders ===== //
    // gSystem->mkdir(Form("roots/bkg/"), kTRUE);
    // gSystem->mkdir(Form("figs/bkg"), kTRUE);

    // // ===== kinematic cut ===== //
    // // cut is applied in f_err dataset
}



void buildBkgModel(RooWorkspace *ws, int nGauss)
{
    // // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // // *-* Build model *-*//
    // // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    // // ===== make res model ===== //
    // cout << "\nn_bkg: " << ws->var("n_bkg")->getVal() << "+/-" << ws->var("n_bkg")->getError() << endl;

    // auto n_bkg_ctau = (RooRealVar *)ws->var("n_bkg")->Clone("n_bkg_ctau");
    // ws->import(*n_bkg_ctau);

    // cout << "Bias of gaus1: " << ws->var("ctau1_CtauRes")->getVal() << endl;
    // cout << "Fraction of gaus1: " << ws->var("f_CtauRes")->getVal() << "+/-" << ws->var("f_CtauRes")->getError() << endl;
    // cout << "Sigma of gaus1: " << ws->var("s1_CtauRes")->getVal() << endl;

    // ws->factory("zeroMean[0.0]");
    // ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes1", "ctau3D",
    //                      "ctau1_CtauRes", //"ctau1_CtauRes",
    //                      "s1_CtauRes",
    //                      "One",
    //                      "ctau3DErr"));
    // ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes2", "ctau3D",
    //                      "ctau2_CtauRes", //"ctau2_CtauRes",
    //                      "s2_CtauRes",
    //                      "One",
    //                      "ctau3DErr"));
    // // ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes3", "ctau3D",
    // //       "ctau3_CtauRes", //"ctau3_CtauRes",
    // //       "s3_CtauRes",
    // //       "zeroMean",
    // //       "ctau3DErr"
    // //       ));
    // // ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "ctauRes32", "ctauRes3", "ctauRes2", "f2_CtauRes"));
    // ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "pdfCTAURES", "ctauRes1", "ctauRes2", "f_CtauRes"));

    // // ===== buld bkg decay model ===== //
    // ws->factory("b_Bkg[0.1, 1e-3, 1]"); // NP fraction for bkg - Not final b fraciton
    // ws->factory("fDFSS[0.1, 1e-3, 1.]");
    // ws->factory("fDLIV[0.1, 1e-3, 1.]");
    // ws->factory("lambdaDDS_Bkg[0.5, 0, 10]");
    // ws->factory("lambdaDF_Bkg[ 0.5, 0, 10]");
    // // ws->factory("lambdaDF_Bkg2[ 0.5, 0, 10]");
    // ws->factory("lambdaDSS_Bkg[0.5, 0, 10]");
    // // ws->factory("lambdaDSS_Bkg2[0.5, 0, 10]");
    // // ws->factory("fSS[0.1, 1e-3, 1.]");
    // // ws->factory("fFF[0.1, 1e-3, 1.]");

    // // make 3 exp
    // ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUDSS1", "ctau3D", "lambdaDSS_Bkg", "pdfCTAURES"));
    // ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUDSS2", "ctau3D", "lambdaDSS_Bkg2[0.1,0.01,1]", "pdfCTAURES"));
    // ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUDSS", "fSS[0.5,0.01,1]", "pdfCTAUDSS1", "pdfCTAUDSS2"));

    // ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", "pdfCTAUDF1", "ctau3D", "lambdaDF_Bkg", "pdfCTAURES"));
    // ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", "pdfCTAUDF2", "ctau3D", "lambdaDF_Bkg2[0.1,0.01,1]", "pdfCTAURES"));
    // ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUDF", "fFF[0.5,0.01,1]", "pdfCTAUDF1", "pdfCTAUDF2"));
    // // ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", "pdfCTAUDF1", "ctau3D", "lambdaDF_Bkg2[0.1,0.01,1]", "pdfCTAURES"));
    // // ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUDF", "fFS[0.5,0.01,1]", "pdfCTAUDF1", "pdfCTAUDF1"));

    // ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", "pdfCTAUDDS", "ctau3D", "lambdaDDS_Bkg", "pdfCTAURES"));
    // // ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", "pdfCTAUDDS2", "ctau3D", "lambdaDDS_Bkg2[0.1,0.01,1]", "pdfCTAURES"));
    // // ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUDDS", "fDDS[0.5,0.01,1]", "pdfCTAUDDS1", "pdfCTAUDDS2"));

    // ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAU1", "fDFSS", "pdfCTAUDSS", "pdfCTAUDF"));
    // ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUCOND_BkgNoPR", "fDLIV", "pdfCTAU1", "pdfCTAUDDS")); // NP

    // // ctauBkgPR
    // ws->factory(Form("SUM::%s(%s)", "pdfCTAUCOND_BkgPR", "pdfCTAURES"));
    // ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUCOND_Bkg", "b_Bkg", "pdfCTAUCOND_BkgNoPR", "pdfCTAUCOND_BkgPR"));

    // // total
    // ws->factory(Form("RooExtendPdf::%s(%s,%s)", "pdfTot_Bkg", "pdfCTAUCOND_Bkg", "n_bkg_ctau")); // n_bkg_ctau is number of bkg from ds_errBkg

    // RooProdPdf pdfPR("pdfCTAU_BkgPR", "", *ws->pdf("errBkgModel"), Conditional(*ws->pdf("pdfCTAUCOND_BkgPR"), RooArgList(*ws->var("ctau3D"))));
    // ws->import(pdfPR);
    // RooProdPdf pdfNoPR("pdfCTAU_BkgNoPR", "", *ws->pdf("errBkgModel"), Conditional(*ws->pdf("pdfCTAUCOND_BkgNoPR"), RooArgList(*ws->var("ctau3D"))));
    // ws->import(pdfNoPR);
    // ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAU_Bkg", "b_Bkg", "pdfCTAU_BkgNoPR", "pdfCTAU_BkgPR"));
    // RooAbsPdf *ctauBkgModel = new RooAddPdf("ctauBkgModel", "ctauBkgModel", RooArgList(*ws->pdf("pdfCTAU_Bkg")),  // originalName bkg_fit_pdfRooArgList(*ws->var("n_bkg_ctau")));
    // ws->import(*ctauBkgModel); // originalName bkg_fit_pdf
    // // RooAbsPdf* ctauBkgModel = ctauBkgModel = new RooAddPdf("pdfTot_Bkg","pdfTot_Bkg",*ws->pdf("pdfCTAUCOND_Bkg"), *ws->var("n_bkg_ctau"));
    // // ws->import(*pdfCTAUCOND_Bkg);
}

void perfromBkgFit(RooWorkspace *ws)
{
    // // ===== fit here ===== //
    // // buld dataset for fit

    // ws->pdf("ctauBkgModel")->setNormRange("ctauWindow"); // originalName bkg_fit_pdf

    // // to avoid ternimal breaking due to too long Tracing Error message
    // // It's okay as long as minimizer find good fitting point
    // // Always check fit results.
    // RooMsgService::instance().getStream(0).removeTopic(Tracing);

    // bool isWeighted = ws->data("ds_errBkg")->isWeighted();
    // // RooFitResult* fit_bkg = ws->pdf("ctauBkgModel")->fitTo(*ds_errBkg, Save(), Range("ctauRange"), Extended(kTRUE),  // originalName bkg_fit_pdfNumCPU(12));

    // RooFitResult *fit_bkg = ws->pdf("ctauBkgModel")->fitTo(*dataToFit, Save(), Extended(1), NumCPU(12), SumW2Error // originalName bkg_fit_pdf(isWeighted), Strategy(2), RecoverFromUndefinedRegions(1.), PrintEvalErrors(-1), PrintLevel(0));

    // ws->import(*fit_bkg, "fit_bkg");

    // // fis parameters
    // ws->var("b_Bkg")->setConstant(kTRUE);
    // ws->var("fDFSS")->setConstant(kTRUE);
    // ws->var("fDLIV")->setConstant(kTRUE);
    // ws->var("lambdaDDS_Bkg")->setConstant(kTRUE);
    // ws->var("lambdaDF_Bkg")->setConstant(kTRUE);
    // ws->var("lambdaDF_Bkg2")->setConstant(kTRUE);
    // ws->var("lambdaDSS_Bkg")->setConstant(kTRUE);
    // ws->var("lambdaDSS_Bkg2")->setConstant(kTRUE);
    // ws->var("fSS")->setConstant(kTRUE);
    // ws->var("fFF")->setConstant(kTRUE);
}

void drawBkgPlot(RooWorkspace *ws, double ctauResMin, double ctauResMax, int nGauss, const std::string &outpath)
{
  // // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  // // *-* Plotting - very long but simple *-*//
  // // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  // normalization
  // double normDSTot = 1.0;
  // if (ws->data("dataToFit"))
  // {
  //   normDSTot = ws->data("dataToFit")->sumEntries() / ws->data("dataToFit")->sumEntries();
  // }
  // double normBkg = 1.0;
  // if (ws->data("ds_errBkg"))
  // {
  //   normBkg = ws->data("dataToFit")->sumEntries() * normDSTot / ws->data("ds_errBkg")->sumEntries();
  // }

    // TCanvas *c_bkg = new TCanvas("c_bkg", "My plots", 800, 800);
    // c_bkg->cd();
    // TPad *pad_bkg1 = new TPad("pad_bkg1", "pad_bkg1", 0, 0.16, 1, 1.0);
    // pad_bkg1->SetTicks(1, 1);
    // pad_bkg1->Draw();
    // pad_bkg1->cd();
    // RooPlot *ctau_frame1 = ws->var("ctau3D")->frame(Bins(nCtauBins), Range(ctauLow, ctauHigh)); // bins
    // ctau_frame1->SetTitle("");
    // pad_bkg1->cd();
    // gPad->SetLogy();
    // RooPlot *ctau_frame2 = (RooPlot *)ctau_frame1->Clone("ctau_frame2");

    // ctau_frame2->updateNormVars(RooArgSet(*ws->var("ctau3D")));
    // // ctau_frame2->updateNormVars(RooArgSet(*ws->var("ctau3D")));

    // // ws->data("ds_errBkg")->plotOn(ctau_frame2, Name("data_ctauBkg"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.7)); // to show full ranges

    // RooDataSet *expDataY = (RooDataSet *)dataToFit->reduce(*ws->var("ctau3DErr"));
    // RooAbsData *binnedDataY = expDataY->binnedClone();

    // // draw data point to provide normalization info to PDFs
    // ws->data("dataToFit")->plotOn(ctau_frame2, Name("data_ctauBkg"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.7)); // data points inside range

    // // ws->pdf("ctauBkgModel")->plotOn(ctau_frame2, ProjWData(*binnedDataY), LineColor(kCyan), LineStyle(kDotted)); // originalName bkg_fit_pdf
    // // ws->pdf("ctauBkgModel")->plotOn(ctau_frame2, LineColor(kCyan), LineStyle(kDotted)); // originalName bkg_fit_pdf
    // // ws->pdf("ctauBkgModel")->forceNumInt(kTRUE); // originalName bkg_fit_pdf
    // ws->pdf("ctauBkgModel")->plotOn(ctau_frame2, Name("ctauBkg_Tot"), NumCPU(12), LineColor(kBlack)); // originalName bkg_fit_pdf

    // // draw data point again - to bring point forward
    // // ws->data("dataToFit")->plotOn(ctau_frame2, Name("data_ctauBkg"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.7)); // data points inside range

    // if (ws->pdf("pdfCTAU_BkgPR"))
    // {
    //     ws->pdf("ctauBkgModel")->plotOn(ctau_frame2, Name("BKGPR"), Components(RooArgSet(*ws->pdf("pdfCTAU_BkgPR"))),  // originalName bkg_fit_pdfLineColor(44), NumCPU(12), Range("ctauRange"), LineStyle(kDashed), LineWidth(2));
    // }
    // if (ws->pdf("pdfCTAU_BkgNoPR"))
    // {
    //     ws->pdf("ctauBkgModel")->plotOn(ctau_frame2, Name("BKGNoPR"), Components(RooArgSet(*ws->pdf // originalName bkg_fit_pdf("pdfCTAU_BkgNoPR"))), LineColor(8), NumCPU(12), Range("ctauRange"), LineStyle(kDashed), LineWidth(2));
    // }
    // // ws->pdf("ctauBkgModel")->plotOn(ctau_frame2, Name("PDF"), , LineColor(kBlack), Precision(1e-4), Range // originalName bkg_fit_pdf("ctauRange"));
    // // Normalization(normBkg, RooAbsReal::NumEvent), NumCPU(12), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("ds_errBkg"), kTRUE)

    // // set proper y range
    // ctau_frame2->GetYaxis()->SetRangeUser(10e-2, 10e7);
    // TH1 *h = ws->data("ds_errBkg")->createHistogram("hist", *ws->var("ctau3D"), Binning(ctau_frame1->GetNbinsX(), ctau_frame1->GetXaxis()->GetXmin(), ctau_frame1->GetXaxis()->GetXmax()));
    // Double_t YMax = h->GetBinContent(h->GetMaximumBin());
    // Double_t YMin = 1e99;
    // for (int i = 1; i <= h->GetNbinsX(); i++)
    //     if (h->GetBinContent(i) > 0)
    //         YMin = min(YMin, h->GetBinContent(i));
    // Double_t Yup(0.), Ydown(0.);
    // // Yup = YMax*TMath::Power((YMax/0.1), 0.5);
    // Yup = YMax * TMath::Power((YMax / 0.01), 0.5);
    // Ydown = 0.01;
    // ctau_frame2->GetYaxis()->SetRangeUser(Ydown, Yup);
    // ctau_frame2->GetXaxis()->SetRangeUser(-4, 7);
    // ctau_frame2->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
    // ctau_frame2->SetFillStyle(4000);
    // ctau_frame2->GetYaxis()->SetTitleOffset(1.43);
    // ctau_frame2->GetXaxis()->SetLabelSize(0);
    // ctau_frame2->GetXaxis()->SetTitleSize(0);

    // // ===== draw lines to present Res rantge ===== //
    // TLine *minline = new TLine(ctauMin, 0.0, ctauMin, Ydown * TMath::Power((Yup / Ydown), 0.4));
    // minline->SetLineStyle(2);
    // minline->SetLineColor(1);
    // minline->SetLineWidth(3);
    // ctau_frame2->addObject(minline);
    // TLine *maxline = new TLine(ctauMax, 0.0, ctauMax, Ydown * TMath::Power((Yup / Ydown), 0.4));
    // maxline->SetLineStyle(2);
    // maxline->SetLineColor(1);
    // maxline->SetLineWidth(3);
    // ctau_frame2->addObject(maxline);

    // // ===== legend ===== //
    // ctau_frame2->Draw();
    // TLegend *leg_bkg = new TLegend(text_x + 0.25, text_y + 0.04, text_x + 0.38, text_y - 0.15);
    // leg_bkg->SetTextSize(text_size);
    // leg_bkg->SetTextFont(43);
    // leg_bkg->SetBorderSize(0);
    // leg_bkg->AddEntry(ctau_frame2->findObject("data_ctauBkg"), "Data_Bkg", "pe");
    // // leg_bkg->AddEntry(ctau_frame2->findObject("BKGPR"), "Bkg PR", "fl");
    // // leg_bkg->AddEntry(ctau_frame2->findObject("BKGNoPR"), "Bkg NonPR", "fl");
    // leg_bkg->Draw("same");

    // // ===== draw latex ===== //
    // drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x, text_y, text_color, text_size);
    // if (yLow == 0)
    //     drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x, text_y - y_diff, text_color, text_size);
    // else if (yLow != 0)
    //     drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x, text_y - y_diff, text_color, text_size);
    // drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x, text_y - y_diff * 2, text_color, text_size);
    // drawText(Form("N_{Bkg} = %.f #pm %.f", ws->var("n_bkg")->getVal(), ws->var("n_bkg")->getError()), text_x + 0.5, text_y, text_color, text_size);
    // drawText(Form("b_{Bkg} = %.4f #pm %.4f", ws->var("b_Bkg")->getVal(), ws->var("b_Bkg")->getError()), text_x + 0.5, text_y - y_diff * 1, text_color, text_size);
    // drawText(Form("fDFSS = %.4f #pm %.4f", ws->var("fDFSS")->getVal(), ws->var("fDFSS")->getError()), text_x + 0.5, text_y - y_diff * 2, text_color, text_size);
    // drawText(Form("fDLIV = %.4f #pm %.4f", ws->var("fDLIV")->getVal(), ws->var("fDLIV")->getError()), text_x + 0.5, text_y - y_diff * 3, text_color, text_size);
    // drawText(Form("#lambdaDDS_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDDS_Bkg")->getVal(), ws->var("lambdaDDS_Bkg")->getError()), text_x + 0.5, text_y - y_diff * 4, text_color, text_size);
    // drawText(Form("#lambdaDF_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDF_Bkg")->getVal(), ws->var("lambdaDF_Bkg")->getError()), text_x + 0.5, text_y - y_diff * 5, text_color, text_size);
    // drawText(Form("#lambdaDSS_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDSS_Bkg")->getVal(), ws->var("lambdaDSS_Bkg")->getError()), text_x + 0.5, text_y - y_diff * 6, text_color, text_size);

    // // ===== draw pull ===== //
    // TPad *pull_pad = new TPad("pull_pad", "pull_pad", 0, 0.006, 1, 0.227);
    // c_bkg->cd();
    // pull_pad->Draw();
    // pull_pad->cd();
    // pull_pad->SetTopMargin(0); // Upper and lower plot are joined
    // pull_pad->SetBottomMargin(0.67);
    // pull_pad->SetBottomMargin(0.4);
    // pull_pad->SetFillStyle(4000);
    // pull_pad->SetFrameFillStyle(4000);
    // pull_pad->SetTicks(1, 1);

    // RooPlot *tmp_pull_frame = (RooPlot *)ctau_frame2->Clone("tmp_pull_frame");
    // RooHist *pull_hist = tmp_pull_frame->pullHist("data_ctauBkg", "ctauBkg_Tot", true);
    // // RooHist *pull_hist = tmp_pull_frame->pullHist("data_ctauBkg", "ctauBkg_Tot", true); // flow code
    // pull_hist->SetMarkerSize(0.8);
    // RooPlot *pull_frame = ws->var("ctau3D")->frame(Title("Pull Distribution"), Bins(nCtauBins), Range(ctauLow, ctauHigh));
    // pull_frame->addPlotable(pull_hist, "PX");
    // pull_frame->SetTitle("");
    // pull_frame->SetTitleSize(0);
    // pull_frame->GetYaxis()->SetTitleOffset(0.3);
    // pull_frame->GetYaxis()->SetTitle("Pull");
    // pull_frame->GetYaxis()->SetTitleSize(0.08);
    // pull_frame->GetYaxis()->SetLabelSize(0.08);
    // pull_frame->GetYaxis()->SetRangeUser(-50, 50);
    // pull_frame->GetYaxis()->CenterTitle();
    // pull_frame->GetYaxis()->SetTickSize(0.04);
    // pull_frame->GetYaxis()->SetNdivisions(404);
    // pull_frame->GetXaxis()->SetTickSize(0.03);

    // pull_frame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
    // pull_frame->GetXaxis()->SetTitleOffset(1.55);
    // pull_frame->GetXaxis()->SetLabelOffset(0.04);
    // pull_frame->GetXaxis()->SetLabelSize(0.08);
    // pull_frame->GetXaxis()->SetTitleSize(0.08);
    // pull_frame->GetXaxis()->CenterTitle();
    // pull_frame->Draw();

    // TLine *lD = new TLine(ctauLow, 0, ctauHigh, 0);
    // lD->SetLineStyle(1);
    // lD->Draw("same");

    // // ===== horizontoal lines on pull hist ===== //
    // // maybe we can make function for them...
    // TLine *line_at_2 = new TLine(pull_frame->GetXaxis()->GetXmin(), 2, pull_frame->GetXaxis()->GetXmax(), 2);
    // line_at_2->SetLineColor(kBlack);
    // line_at_2->SetLineStyle(3);
    // line_at_2->SetLineWidth(1);
    // line_at_2->Draw();

    // TLine *line_at_4 = new TLine(pull_frame->GetXaxis()->GetXmin(), 4, pull_frame->GetXaxis()->GetXmax(), 4);
    // line_at_4->SetLineColor(kBlack);
    // line_at_4->SetLineStyle(3);
    // line_at_4->SetLineWidth(2);
    // line_at_4->Draw();

    // TLine *line_at_8 = new TLine(pull_frame->GetXaxis()->GetXmin(), 8, pull_frame->GetXaxis()->GetXmax(), 8);
    // line_at_8->SetLineColor(kBlack);
    // line_at_8->SetLineStyle(3);
    // line_at_8->SetLineWidth(1);
    // line_at_8->Draw();

    // // negative
    // TLine *line_at_n2 = new TLine(pull_frame->GetXaxis()->GetXmin(), -2, pull_frame->GetXaxis()->GetXmax(), -2);
    // line_at_n2->SetLineColor(kBlack);
    // line_at_n2->SetLineStyle(3);
    // line_at_n2->SetLineWidth(1);
    // line_at_n2->Draw();

    // TLine *line_at_n4 = new TLine(pull_frame->GetXaxis()->GetXmin(), -4, pull_frame->GetXaxis()->GetXmax(), -4);
    // line_at_n4->SetLineColor(kBlack);
    // line_at_n4->SetLineStyle(3);
    // line_at_n4->SetLineWidth(2);
    // line_at_n4->Draw();

    // TLine *line_at_n8 = new TLine(pull_frame->GetXaxis()->GetXmin(), -8, pull_frame->GetXaxis()->GetXmax(), -8);
    // line_at_n8->SetLineColor(kBlack);
    // line_at_n8->SetLineStyle(3);
    // line_at_n8->SetLineWidth(1);
    // line_at_n8->Draw();

    // // ===== chi2 ===== //
    // printChi2(ws, pull_pad, tmp_pull_frame, fit_bkg, "ctau3D", "data_ctauBkg", "ctauBkg_Tot", nCtauBins, false);

    // pull_pad->Update();

    // // ===== draw canvas ===== //
    // c_bkg->Update();
    // c_bkg->SaveAs(("figs/" + out_ss + ".png").c_str());
    // c_bkg->SaveAs(("figs/" + out_ss + ".pdf").c_str());
}

void drawBkgRatioPlot(RooWorkspace *ws, double ctauResMin, double ctauResMax, int nGauss, const std::string &outpath)
{
    // // ===== ratio test cnavas ===== //
    // auto c = new TCanvas();
    // c->cd();
    // TH1 *h2 = dataToFit->createHistogram("h2", *ws->var("ctau3D"), Binning(ctau_frame1->GetNbinsX(), ctau_frame1->GetXaxis()->GetXmin(), ctau_frame1->GetXaxis()->GetXmax()));
    // // auto h1 = data->createHistogram("x");

    // auto f1 = ws->pdf("ctauBkgModel")->asTF(*ws->var("ctau3D"), RooArgSet(), *ws->var("ctau3D")); // originalName bkg_fit_pdf
    // auto fcor = new TF1("fcor", [&](double *x, double *p)
    //                     { return f1->EvalPar(x, p) * h2->GetSumOfWeights() * h2->GetBinWidth(1); }, ctau_frame1->GetXaxis()->GetXmin(), ctau_frame1->GetXaxis()->GetXmax(), 0);
    // auto f2 = (TF1 *)fcor->Clone("f2"); // for drawing
    // double f2_max = f2->Eval(0);
    // cout << "Eval f2 at ctau3D = 0: " << f2->Eval(0) << endl;
    // h2->SetMaximum(f2_max * 5); // usually, f2's maximum is bigger than a hist's
    // gPad->SetLogy();
    // h2->Draw("");
    // f2->Draw("same");

    // int n_tmp = h->GetNbinsX();
    // for (int bin = 1; bin <= n_tmp; ++bin)
    // {
    //     Double_t binCenter = h2->GetBinCenter(bin);
    //     Double_t result = f2->Eval(binCenter);
    //     Double_t binContent = h2->GetBinContent(bin);

    //     cout << "Bin center: " << binCenter
    //          << ", Bin content: " << binContent
    //          << ", f(" << binCenter << ") = " << result << std::endl;
    // }

    // Double_t chi2 = h2->Chisquare(f2, "L");
    // cout << "Chi2: " << chi2 << endl;
    // c->Draw();
    // c->SaveAs(("figs/" + out_ss + "_TF1.png").c_str());
}

void saveOutput(RooWorkspace *ws, const std::string &outpath)
{
    // // ===== Export results ===== //
    // auto out_file = new TFile(("roots/" + out_ss + ".root").c_str(), "recreate");

    // fit_bkg->Write();
    // ws->Write();
    // out_file->Close();

    // fit_bkg->Print("V");
}