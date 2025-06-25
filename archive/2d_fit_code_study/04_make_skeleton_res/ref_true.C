#include <iostream>
#include <string>

#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "TH1.h"
#include "TString.h"
#include "TSystem.h" // gSystem->mkdir

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFormulaVar.h"
#include "RooTruthModel.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooMsgService.h"


// global variable - will be memober variable
const float ctauLow = 0; // is it good for ctau3Dtrue too? -> try anyway
const float ctauHigh = 10.;
const int nCtauTrueBins = 50;

using namespace std;
using namespace RooFit;

void setupWorkspace(RooWorkspace *ws, const std::string &gen_path, const std::string &cut);
void processDataset(RooWorkspace *ws);
void buildTrueModel(RooWorkspace *ws, int nExp);
RooFitResult *perfromTrueFit(RooWorkspace *ws);
void drawTruePlot(RooWorkspace *ws, int nExp, const std::string &outpath);
void saveOutput(RooWorkspace *ws, RooFitResult *fitRes, const std::string &outpath);

void ref_true()
{
  /**
   * @brief import RecoGen RooDataset, use ctau3DTrue (ctau3D from Gen) to make NP decay model
   */
  std::cout << "===== Start ref_true.C =====\n";
  TStopwatch *t = new TStopwatch;
  t->Start();

  // make a cut -> member function and variables
  float ptLow = 9, ptHigh = 12;
  float yLow = 0, yHigh = 1.6;
  double massLow = 2.6, massHigh = 3.5;
  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f", ptLow, ptHigh, yLow, yHigh, massLow, massHigh);
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2) && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1) && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&"; // 2018 acceptance cut
  TString final_cut = accCut + kineCut; // Gen: no opposite sign cut.

  // setup workspace
  RooWorkspace *ws_true = new RooWorkspace("ws_true", "");
  setupWorkspace(ws_true, "../../../files_roodata/RooDataSet_miniAOD_isMC1_NP_Jpsi_GenOnly_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250307.root", final_cut.Data());

  // process dataset
  processDataset(ws_true);

  // build ctautrue model;
  int nExp = 2;
  buildTrueModel(ws_true, nExp);

  // do fit
  RooFitResult *fit_result = perfromTrueFit(ws_true);
  // fit_result->correlationMatrix().Print();

  // draw plots
  drawTruePlot(ws_true, nExp, "figs/true.png");

  // save output
  saveOutput(ws_true, fit_result, "roots/true.root");

  t->Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t->RealTime(), t->CpuTime());
  std::cout << "===== Start ref_true.C =====\n";
}

void setupWorkspace(RooWorkspace *ws, const std::string &gen_path, const std::string &cut)
{
  std::cout << "===== start setupWorkspace() =====\n";

  TFile *f_recoGen = new TFile(gen_path.c_str());

  if (!f_recoGen)
  {
    std::cout << "Error: Could not open one of the files.\n";
    return;
  }

  // import recoGenOnly - FYI: RecoGenOnly doesn't have weight - (Check: Or applied making RooDataset??)
  RooDataSet *ds_trueTmp = (RooDataSet *)f_recoGen->Get("dataset");
  RooDataSet *ds_reduced = (RooDataSet *)ds_trueTmp->reduce(cut.c_str());
  ds_reduced->SetName("ds_reduced");
  ws->import(*ds_reduced);

  // ws->Print("V");
  delete f_recoGen;
  std::cout << "===== finish setupWorkspace() =====\n";
}

void processDataset(RooWorkspace *ws)
{
  std::cout << "===== Start processDataset() =====\n";

  // set fit range -> User config?
  const double fitRangeMin = 0.0001; // should be bigger than 0.00001 -> Make a distance from the edge (0)
  const double fitRangeMax = 10;

  // apply fit range cut - dataset for fitting and plotting
  RooDataSet *ds_true = (RooDataSet *)ws->data("ds_reduced")->reduce(Form("ctau3D>=%.4f&&ctau3D<%.4f", fitRangeMin, fitRangeMax));
  ds_true->SetName("ds_true");
  ws->import(*ds_true);

  // set ranges
  ws->var("ctau3D")->setRange("plotRange", ctauLow, ctauHigh);
  ws->var("ctau3D")->setRange("fitRange", fitRangeMin, fitRangeMax);
  std::cout << "===== Finish processDataset() =====\n";
}

void buildTrueModel(RooWorkspace *ws, int nExp)
{
  std::cout << "===== Start buildTrueModel() =====\n";

  // parameters
  double entryTrue = ws->data("ds_true")->numEntries();
  ws->factory(Form("nSig_true[%.f, 0., %.f]", entryTrue, entryTrue * 2));
  ws->factory("lambdaDSS_true[0.33, 1e-3, 10.]");
  if (nExp >= 2)
  {
    ws->factory("delta_lambda2[0.2, 1e-3, 10.]");
    ws->factory("RooFormulaVar::lambdaDSS2_true('@0+@1', {lambdaDSS_true, delta_lambda2})");
    ws->factory("fDSS_true[0.45, 0.01, 1.]");
  }
  if (nExp >= 3)
  {
    ws->factory("delta_lambda3[0.1, 1e-3, 10.]");
    ws->factory("RooFormulaVar::lambdaDSS3_true('@0+@1', {lambdaDSS2_true, delta_lambda3})");
    ws->factory("fDSS23_true[0.4, 0.01, 1.]");
  }

  // res model for Gen
  ws->factory("TruthModel::pdfRES_true(ctau3D)"); // truth
  // ws->factory("sigmaRes_true[0.01, 0.00001, 0.1]");
  // ws->factory("meanRes_true[0]");
  // ws->var("meanRes_true")->setConstant(kTRUE);
  // ws->factory("GaussModel::pdfRES_true(ctau3D, meanRes_true, sigmaRes_true)"); // consider gauss

  // decay model
  ws->factory("Decay::pdfDSS1_true(ctau3D, lambdaDSS_true, pdfRES_true, RooDecay::SingleSided)");
  if (nExp >= 2)
    ws->factory("Decay::pdfDSS2_true(ctau3D, lambdaDSS2_true, pdfRES_true, RooDecay::SingleSided)");
  if (nExp >= 3)
    ws->factory("Decay::pdfDSS3_true(ctau3D, lambdaDSS3_true, pdfRES_true, RooDecay::SingleSided)");

  // combine decay models
  if (nExp == 1)
  {
    ws->factory("SUM::pdfDECAY_true(pdfDSS1_true)");
  }
  else if (nExp == 2)
  {
    ws->factory("SUM::pdfDECAY_true(fDSS_true * pdfDSS1_true, pdfDSS2_true)");
  }
  else if (nExp == 3)
  {
    ws->factory("SUM::pdfDSS23_true(fDSS23_true * pdfDSS2_true, pdfDSS3_true)");
    ws->factory("SUM::pdfDECAY_true(fDSS_true * pdfDSS1_true, pdfDSS23_true)");
  }

  // true fit model
  // ws->factory("SUM::trueModel(nSig_true * pdfDECAY_true)");
  auto trueModel = new RooExtendPdf("trueModel", "", *ws->pdf("pdfDECAY_true"), *ws->var("nSig_true"));
  ws->import(*trueModel);

  std::cout << "===== Finish buildTrueModel() =====\n";
}

RooFitResult *perfromTrueFit(RooWorkspace *ws)
{
  std::cout << "===== start perfromTrueFit() =====\n";

  bool isWeighted = ws->data("ds_true")->isWeighted();
  auto fit_true = ws->pdf("trueModel")->fitTo(*ws->data("ds_true"), Save(1), Extended(1), NumCPU(24), Range("fitRange"), Strategy(2), RecoverFromUndefinedRegions(1.), PrintEvalErrors(-1));

  fit_true->Print("V");
  ws->import(*fit_true, "fit_true");

  // fix parameters for next steps
  std::cout << "fixing floating parameters to constant (after fitting)\n";
  const RooArgList &floatted_params = fit_true->floatParsFinal();
  for (auto arg : floatted_params)
  {
    RooRealVar *param = dynamic_cast<RooRealVar *>(arg);
    if (param)
    {
      std::cout << "Fixing parameter: " << param->GetName() << "\n";
      ws->var(param->GetName())->setConstant(kTRUE);
    }
  }

  std::cout << "===== finish perfromTrueFit() =====\n";
  return fit_true;
}

void drawTruePlot(RooWorkspace *ws, int nExp, const std::string &outpath)
{
  std::cout << "===== Start drawTruePlot() =====\n";

  // make canvas
  TCanvas *c_true = new TCanvas("c_true", "", 800, 600);
  c_true->cd();

  // make true pad
  TPad *truePad = new TPad("truePad", "", 0, 0.25, 0.98, 1.0);
  truePad->SetTicks(1, 1);
  truePad->Draw();

  truePad->cd();

  // true frame
  RooPlot *trueFrame = ws->var("ctau3D")->frame(Bins(nCtauTrueBins), Range("plotRange"));

  // draw dataset and pdfs
  ws->data("ds_true")->plotOn(trueFrame, Name("ds_true")); // DataError(RooAbsData::SumW2)
  ws->pdf("trueModel")->plotOn(trueFrame, Range("plotRange"), NormRange("fitRange"), Name("trueModel"));
  if (nExp >= 2)
  {
    ws->pdf("trueModel")->plotOn(trueFrame, Name("comp1"), Components(*ws->pdf("pdfDSS1_true")), LineStyle(kDashed), LineColor(kGreen + 2), Range("plotRange"), NormRange("fitRange"));
    ws->pdf("trueModel")->plotOn(trueFrame, Name("comp2"), Components(*ws->pdf("pdfDSS2_true")), LineStyle(kDashed), LineColor(kAzure - 4), Range("plotRange"), NormRange("fitRange"));
  }
  if (nExp >= 3)
  {
    ws->pdf("trueModel")->plotOn(trueFrame, Name("comp3"), Components(*ws->pdf("pdfDSS3_true")), LineStyle(kDotted), LineColor(kOrange + 1), Range("plotRange"), NormRange("fitRange"));
  }

  // set y-axis range
  gPad->SetLogy();
  TH1 *h_tmp = ws->data("ds_true")->createHistogram("h_tmp", *ws->var("ctau3D"), Binning(trueFrame->GetNbinsX(), trueFrame->GetXaxis()->GetXmin(), trueFrame->GetXaxis()->GetXmax()));
  Double_t YMax = h_tmp->GetBinContent(h_tmp->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= h_tmp->GetNbinsX(); i++)
    if (h_tmp->GetBinContent(i) > 0)
      YMin = min(YMin, h_tmp->GetBinContent(i));
  double Ydown;
  double Yup;
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.1 / (1.0 - 0.1 - 0.4))));
  if (Ydown <= 0)
    Ydown = 0.1;
  Yup = YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.1 - 0.4)));
  trueFrame->GetYaxis()->SetRangeUser(Ydown, Yup);
  delete h_tmp;

  // cosmetics
  trueFrame->SetFillStyle(4000);
  trueFrame->GetYaxis()->SetTitleOffset(1.43);
  trueFrame->GetXaxis()->SetLabelSize(0); // to hide x-axis label
  trueFrame->GetXaxis()->SetTitleSize(0);
  trueFrame->Draw();

  // legend
  TLegend *leg = new TLegend(0.6, 0.65, 0.93, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(trueFrame->findObject("ds_true"), "Gen ctau3D", "pe");
  leg->AddEntry(trueFrame->findObject("trueModel"), "Total fit", "l");
  if (nExp >= 2)
  {
    leg->AddEntry(trueFrame->findObject("comp1"), "Decay Exp1", "l");
    leg->AddEntry(trueFrame->findObject("comp2"), "Decay Exp2", "l");
  }
  if (nExp >= 3)
  {
    leg->AddEntry(trueFrame->findObject("comp3"), "Decay Exp3", "l");
  }
  leg->Draw();

  // pull part
  TPad *pullPad = new TPad("pullPad", "", 0, 0.001, 0.98, 0.32);
  c_true->cd();
  pullPad->Draw();
  pullPad->cd();
  pullPad->SetTopMargin(0); // Upper and lower plot are joined
  pullPad->SetBottomMargin(0.67);
  pullPad->SetBottomMargin(0.4);
  pullPad->SetFillStyle(4000);
  pullPad->SetFrameFillStyle(4000);
  pullPad->SetTicks(1, 1);

  // to ignore empty residual bin warning - it's not a matter.
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);

  RooPlot *pullTmp = (RooPlot *)trueFrame->Clone("pullTmp");
  RooHist *pullHist = pullTmp->pullHist("ds_true", "trueModel", true);
  pullHist->SetMarkerSize(0.8);
  RooPlot *pullFrame = ws->var("ctau3D")->frame(Title(""), Range(ctauLow, ctauHigh));
  pullFrame->addPlotable(pullHist, "PX");

  pullFrame->SetTitle("");
  pullFrame->SetTitleSize(0);
  pullFrame->GetYaxis()->SetTitleOffset(0.3);
  pullFrame->GetYaxis()->SetTitle("Pull");
  pullFrame->GetYaxis()->SetTitleSize(0.08);
  pullFrame->GetYaxis()->SetLabelSize(0.08);
  pullFrame->GetYaxis()->SetRangeUser(-9.5, 9.5);
  pullFrame->GetYaxis()->CenterTitle();
  pullFrame->GetYaxis()->SetTickSize(0.04);
  pullFrame->GetYaxis()->SetNdivisions(404);
  pullFrame->GetXaxis()->SetTickSize(0.03);

  pullFrame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  pullFrame->GetXaxis()->SetTitleOffset(1.55);
  pullFrame->GetXaxis()->SetLabelOffset(0.04);
  pullFrame->GetXaxis()->SetLabelSize(0.08);
  pullFrame->GetXaxis()->SetTitleSize(0.08);
  pullFrame->GetXaxis()->CenterTitle();
  pullFrame->Draw();

  // ===== horizontoal lines on pull hist ===== //
  TLine *line_at_2 = new TLine(pullFrame->GetXaxis()->GetXmin(), 2, pullFrame->GetXaxis()->GetXmax(), 2);
  line_at_2->SetLineColor(kBlack);
  line_at_2->SetLineStyle(3);
  line_at_2->SetLineWidth(1);
  line_at_2->Draw();

  TLine *line_at_4 = new TLine(pullFrame->GetXaxis()->GetXmin(), 4, pullFrame->GetXaxis()->GetXmax(), 4);
  line_at_4->SetLineColor(kBlack);
  line_at_4->SetLineStyle(3);
  line_at_4->SetLineWidth(2);
  line_at_4->Draw();

  TLine *line_at_8 = new TLine(pullFrame->GetXaxis()->GetXmin(), 8, pullFrame->GetXaxis()->GetXmax(), 8);
  line_at_8->SetLineColor(kBlack);
  line_at_8->SetLineStyle(3);
  line_at_8->SetLineWidth(1);
  line_at_8->Draw();

  // negative lines
  TLine *line_at_n2 = new TLine(pullFrame->GetXaxis()->GetXmin(), -2, pullFrame->GetXaxis()->GetXmax(), -2);
  line_at_n2->SetLineColor(kBlack);
  line_at_n2->SetLineStyle(3);
  line_at_n2->SetLineWidth(1);
  line_at_n2->Draw();

  TLine *line_at_n4 = new TLine(pullFrame->GetXaxis()->GetXmin(), -4, pullFrame->GetXaxis()->GetXmax(), -4);
  line_at_n4->SetLineColor(kBlack);
  line_at_n4->SetLineStyle(3);
  line_at_n4->SetLineWidth(2);
  line_at_n4->Draw();

  TLine *line_at_n8 = new TLine(pullFrame->GetXaxis()->GetXmin(), -8, pullFrame->GetXaxis()->GetXmax(), -8);
  line_at_n8->SetLineColor(kBlack);
  line_at_n8->SetLineStyle(3);
  line_at_n8->SetLineWidth(1);
  line_at_n8->Draw();

  // ===== chi2 ===== //
  // printChi2(ws, pullPad, pullTmp, fit_true, "ctau3D", "ds_true", "trueModel", nCtauTrueBins, false);
  pullPad->Update();

  // ===== draw main canvas ===== //
  c_true->Draw();
  c_true->SaveAs(outpath.c_str());

  delete c_true;
  std::cout << "===== Finish drawTruePlot() =====\n";
}

void saveOutput(RooWorkspace *ws, RooFitResult *fitRes, const std::string &outpath)
{
  std::cout << "===== start saveOutput() =====\n";
  ws->SetName("old_ws_true");

  RooWorkspace ws_to_save("ws_to_save", "Workspace with model and fit results");
  ws_to_save.SetName("ws_true");
  ws_to_save.import(*ws->pdf("trueModel"));
  ws_to_save.import(*fitRes);

  auto outfile = new TFile(outpath.c_str(), "recreate");
  ws_to_save.Write();
  outfile->Close();

  delete outfile;

  std::cout << "output saved to: " << outpath << "\n";
  std::cout << "===== finish saveOutput() =====\n";
}

// void setConfig() {
//-> will be member function later
// input parameters
// float cos_low = 0.0, cos_high = 0.1;
// float ptLow = 3;
// float ptHigh = 6.5;
// float yLow = 1.6;
// float yHigh = 2.4;
// int cLow = 0;
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
// // RooMsgService::instance().getStream(0).removeTopic(Caching); // in RooFit:: namespace
// // RooMsgService::instance().getStream(1).removeTopic(Caching);
// // RooMsgService::instance().getStream(0).removeTopic(Plotting);
// // RooMsgService::instance().getStream(1).removeTopic(Plotting);
// // RooMsgService::instance().getStream(0).removeTopic(InputArguments);
// // RooMsgService::instance().getStream(1).removeTopic(InputArguments);
// // RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
// // RooMsgService::instance().getStream(1).removeTopic(Minimization);
// // RooMsgService::instance().getStream(1).removeTopic(Caching);
// RooMsgService::instance().setGlobalKillBelow(WARNING);

// // ===== kinematic cut ===== //
// TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f", ptLow, ptHigh, yLow, yHigh, massLow, massHigh);
// TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&"; // 2018 acceptance cut

// TString final_cut = accCut + kineCut;

// // ===== labeling ===== //
// std::ostringstream out_file_path;
// out_file_path << "true/true_pt" << ptLow << "-" << ptHigh
//               << "_absY" << yLow << "-" << yHigh;
// // << "_cos" << cos_low << "_" << cos_high;
// string out_ss = out_file_path.str();

// // ===== make output folders ===== //
// gSystem->mkdir(Form("roots/true/"), kTRUE);
// gSystem->mkdir(Form("figs/true"), kTRUE);
// }