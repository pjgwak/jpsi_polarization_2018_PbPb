#include "SPlotter.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include "TStopwatch.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TH1D.h"
#include "TMath.h"
#include "TStyle.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooMsgService.h"
#include "RooStats/SPlot.h"
#include "../../../headers/polarizationUtilities.h"

using namespace RooFit;

SPlotter::SPlotter(float ptLow_, float ptHigh_, float yLow_, float yHigh_, int cLow_, int cHigh_, float cosLow, float cosHigh)
    : ptLow_(ptLow_), ptHigh_(ptHigh_), yLow_(yLow_), yHigh_(yHigh_), cLow_(cLow_), cHigh_(cHigh_), cosLow_(cosLow), cosHigh_(cosHigh)
{
  // RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
  // RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  gStyle->SetEndErrorSize(0);
  std::cout << "SPlotter object is created\n";
}

SPlotter::~SPlotter()
{
  delete dataFile_;
  delete massFitFile_;
  delete ws_err_;
  delete sPlotData_;
  std::cout << "SPlotter object is destroyed.\n";
}

void SPlotter::run()
{
  TStopwatch t;
  t.Start();
  std::cout << "\n===== Start ctau3DErr sPlot work =====\n";
  setupInputOutput();
  prepareDataset();
  performSPlot();
  buildWeighedObjects();
  plotResults();
  saveResults();
  std::cout << "\n ===== Finish ctau3DErr sPlot work =====\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}

// ===== setup, prepare, sPlot =====
void SPlotter::setupInputOutput()
{
  std::ostringstream out_path;
  out_path << "err/err_pt" << ptLow_ << "-" << ptHigh_ << "_cent" << cLow_ << "-" << cHigh_ << "_absY" << yLow_ << "-" << yHigh_ << "_cos" << cosLow_ << "_" << cosHigh_;
  outputPath_ = out_path.str();

  std::ostringstream mass_in_path;
  mass_in_path << "roots/mass/mass_pt" << ptLow_ << "-" << ptHigh_ << "_cent" << cLow_ << "-" << cHigh_ << "_absY" << yLow_ << "-" << yHigh_ << "_cos" << cosLow_ << "_" << cosHigh_ << ".root";
  massFitInputPath_ = mass_in_path.str();

  gSystem->mkdir("roots/err/", kTRUE);
  gSystem->mkdir("figs/err", kTRUE);
}

void SPlotter::prepareDataset()
{
  massFitFile_ = new TFile(massFitInputPath_.c_str());
  RooWorkspace *ws_mass = (RooWorkspace *)massFitFile_->Get("ws_mass_");
  RooAbsPdf *mass_pdf = (RooAbsPdf *)ws_mass->pdf("mass_pdf_");
  mass_pdf->SetName("mass_pdf");

  dataFile_ = new TFile(dataInputPath_.c_str());
  RooDataSet *ds_tmp = (RooDataSet *)dataFile_->Get("dataset");
  RooDataSet ds_tmpW("ds_tmpW", "weighted dataset", ds_tmp, *ds_tmp->get(), 0, "weight");

  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>2.6 && mass<3.5 && cBin>=%d && cBin<%d", ptLow_, ptHigh_, yLow_, yHigh_, cLow_, cHigh_);
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&"; // 2018 acceptance cut
  TString OS = "recoQQsign==0 &&";
  TString angle_cut = Form("&& cos_ep>%.2f && cos_ep<%.2f", cosLow_, cosHigh_);
  TString final_cut = OS + accCut + kineCut + angle_cut;

  RooDataSet *ds_err = nullptr;
  if (true) // weight on
    ds_err = (RooDataSet *)ds_tmpW.reduce(final_cut.Data());
  else
    ds_err = (RooDataSet *)ds_tmp->reduce(final_cut.Data());

  ds_err->SetName("ds_err");

  ws_err_ = new RooWorkspace("ws_err");
  ws_err_->import(*mass_pdf);
  ws_err_->import(*ds_err);

  ctau3DErr_ = ws_err_->var("ctau3DErr");
  ctau3DErr_->setRange(ctauErrLow, ctauErrHigh);
  nBins_ = std::min(int(std::round((ws_err_->var("ctau3DErr")->getMax() - ws_err_->var("ctau3DErr")->getMin()) / 0.0025)), 100);
}

void SPlotter::performSPlot()
{
  RooRealVar *sigYield = ws_err_->var("nSig_mass_");
  sigYield->SetName("nSig_mass");
  RooRealVar *bkgYield = ws_err_->var("nBkg_mass_");
  bkgYield->SetName("nBkg_mass");
  sigYield->setMin(0);
  bkgYield->setMin(0);

  // nSig_mass, nBkg_mass will be used for sPlot
  RooArgList yieldList;
  yieldList.add(*ws_err_->var("nSig_mass"));
  yieldList.add(*ws_err_->var("nBkg_mass"));

  std::cout << "Sig Yield: " << sigYield->getVal() << " +/- " << sigYield->getError() << "\n";
  std::cout << "Bkg Yield: " << bkgYield->getVal() << " +/- " << bkgYield->getError() << "\n";

  // Signal and Bkg PDFs of mass fit will be used
  // code below is messy, I keep it because it's legacy
  RooDataSet *ds_splot = (RooDataSet *)ws_err_->data("ds_err")->Clone("ds_splot");
  RooArgSet *cloneSet = (RooArgSet *)RooArgSet(*ws_err_->pdf("mass_pdf"), "mass_pdf").snapshot(kTRUE);
  auto clone_mass_pdf = (RooAbsPdf *)cloneSet->find("mass_pdf");
  clone_mass_pdf->setOperMode(RooAbsArg::ADirty, kTRUE);

  // ===== sPlot fit ===== //
  // sPlot: fds_splotit variable O1 to get the weightings of sig and bkg. We will aplly the weightings to variable O2 later.
  RooStats::SPlot sData = RooStats::SPlot("sData", "An SPlot", *ds_splot, clone_mass_pdf, yieldList);

  // the total yield and total weight should be same
  std::cout << "[INFO] Jpsi yield -> Mass Fit:" << ws_err_->var("nSig_mass")->getVal() << ", sWeights :" << sData.GetYieldFromSWeight("nSig_mass") << "\n";
  std::cout << "[INFO] Bkg  yield -> Mass Fit:" << ws_err_->var("nBkg_mass")->getVal() << ", sWeights :" << sData.GetYieldFromSWeight("nBkg_mass") << "\n";

  // individual events have different weights
  for (int i = 0; i < 10; i++)
  {
    // check weights of some events
    std::cout << "nSig_mass weight: " << i << std::right << std::setw(12) << sData.GetSWeight(i, "nSig_mass") << " /  nBkg_mass weight: "
              << std::setw(12) << sData.GetSWeight(i, "nBkg_mass") << " / Total Weight" << std::setw(12) << sData.GetSumOfEventSWeight(i)
              << "\n";
  }

  ws_err_->import(*ds_splot);

  // ===== build RooPlot and histograms ===== //
  // a lot of legacy codes
  RooPlot *err_frame1 = ws_err_->var("ctau3DErr")->frame(Range(ctauErrLow, ctauErrHigh));

  TH1D *hTot = (TH1D *)ws_err_->data("ds_err")->createHistogram(("hTot"), *ws_err_->var("ctau3DErr"), Binning(err_frame1->GetNbinsX(), err_frame1->GetXaxis()->GetXmin(), err_frame1->GetXaxis()->GetXmax()));
  TH1D *hTot_M = (TH1D *)ws_err_->data("ds_err")->createHistogram(("hTot_M"), *ws_err_->var("ctau3DErr"), Binning(nBins_, ctauErrLow, ctauErrHigh));

  // ===== build temporal sPlot weighted datasets ===== //
  RooDataSet *dataw_Sig_b = new RooDataSet("dataw_Sig_b", "TMP_SIG_DATA", (RooDataSet *)ws_err_->data("ds_splot"),
                                           RooArgSet(*ws_err_->var("ctau3DErr"), *ws_err_->var("nSig_mass_sw"), *ws_err_->var("ctau3DRes"), *ws_err_->var("ctau3D"), *ws_err_->var("mass")), 0, "nSig_mass_sw");
  TH1D *hSig = (TH1D *)dataw_Sig_b->createHistogram(("hSig"), *ws_err_->var("ctau3DErr"), Binning(nBins_, ctauErrLow, ctauErrHigh));

  // Bkg
  RooDataSet *dataw_Bkg_b = new RooDataSet("dataw_Bkg_b", "TMP_BKG_DATA", (RooDataSet *)ws_err_->data("ds_splot"),
                                           RooArgSet(*ws_err_->var("ctau3DErr"), *ws_err_->var("nBkg_mass_sw"), *ws_err_->var("ctau3DRes"), *ws_err_->var("ctau3D"), *ws_err_->var("mass")), 0, "nBkg_mass_sw");
  TH1D *hBkg = (TH1D *)dataw_Bkg_b->createHistogram(("hBkg"), *ws_err_->var("ctau3DErr"), Binning(nBins_, ctauErrLow, ctauErrHigh));
  // ===== build sPlot weighted datasets and RooHistPdfs ===== //
  // total ds and pdf
  TH1D *h_tot_test = (TH1D *)ws_err_->data("ds_err")->createHistogram(("h_tot_test"), *ws_err_->var("ctau3DErr"), Binning(newBins_, ctauErrMin_, ctauErrMax_));
  // TH1D* h_tot_test = (TH1D*)ws_err_->data("ds_err")->createHistogram(("h_tot_test"), *ws_err_->var("ctau3DErr"),Binning(nBins_,ctauErrMin_,ctauErrMax_));
  RooDataHist *totHist = new RooDataHist("ds_err", "", *ws_err_->var("ctau3DErr"), h_tot_test);
  RooHistPdf *err_tot_pdf = new RooHistPdf("err_tot_pdf", "Total RooHistPdf of ctauErr", *ws_err_->var("ctau3DErr"), *totHist);

  // bkg ds and pdf
  TH1D *hBkg_w = (TH1D *)dataw_Bkg_b->createHistogram(("hBkg_w"), *ws_err_->var("ctau3DErr"), Binning(newBins_, ctauErrMin_, ctauErrMax_));
  RooDataHist *bkgHist = new RooDataHist("ds_err", "", *ws_err_->var("ctau3DErr"), hBkg_w);
  RooHistPdf *err_bkg_pdf = new RooHistPdf("err_bkg_pdf", "hist pdf", *ws_err_->var("ctau3DErr"), *bkgHist);

  // sig ds and pdf
  TH1D *hSig_w = (TH1D *)dataw_Sig_b->createHistogram(("hSig_w"), *ws_err_->var("ctau3DErr"), Binning(newBins_, ctauErrMin_, ctauErrMax_));
  RooDataHist *sigHist = new RooDataHist("ds_err", "", *ws_err_->var("ctau3DErr"), hSig_w);
  RooHistPdf *err_sig_pdf = new RooHistPdf("err_sig_pdf", "hist pdf", *ws_err_->var("ctau3DErr"), *sigHist);

  // Re-build datasets
  RooDataSet *ds_bkg_err = new RooDataSet("ds_bkg_err", "", (RooDataSet *)ws_err_->data("ds_splot"),
                                          RooArgSet(*ws_err_->var("ctau3DErr"), *ws_err_->var("nBkg_mass_sw"), *ws_err_->var("ctau3DRes"), *ws_err_->var("ctau3D"), *ws_err_->var("mass")), 0, "nBkg_mass_sw");
  RooDataSet *ds_sig_err = new RooDataSet("ds_sig_err", "", (RooDataSet *)ws_err_->data("ds_splot"),
                                          RooArgSet(*ws_err_->var("ctau3DErr"), *ws_err_->var("nSig_mass_sw"), *ws_err_->var("ctau3DRes"), *ws_err_->var("ctau3D"), *ws_err_->var("mass")), 0, "nSig_mass_sw");

  // ===== import sPlot results ===== //
  ws_err_->import(*ds_sig_err);
  ws_err_->import(*ds_bkg_err);
  ws_err_->import(*err_tot_pdf);
  ws_err_->import(*err_sig_pdf);
  ws_err_->import(*err_bkg_pdf);
}

// ===== build, plot, save 함수 =====
void SPlotter::buildWeighedObjects()
{

}

void SPlotter::plotResults()
{
  double minRange = (double)(floor(ctauErrMin_ * 100.) / 100.);
  double maxRange = (double)(ceil(ctauErrMax_ * 100.) / 100.);
  ws_err_->var("ctau3DErr")->setRange("ctauErrWindow", ctauErrMin_, ctauErrMax_);
  RooPlot *err_frame1 = ws_err_->var("ctau3DErr")->frame(Range(ctauErrLow, ctauErrHigh));

  TH1D *hTot = (TH1D *)ws_err_->data("ds_err")->createHistogram(("hTot"), *ws_err_->var("ctau3DErr"), Binning(err_frame1->GetNbinsX(), err_frame1->GetXaxis()->GetXmin(), err_frame1->GetXaxis()->GetXmax()));
  TH1D *hTot_M = (TH1D *)ws_err_->data("ds_err")->createHistogram(("hTot_M"), *ws_err_->var("ctau3DErr"), Binning(nBins_, ctauErrLow, ctauErrHigh));

  // ===== set new ranges of ctau3DErr ===== //
  double ctauErrMin_calc;
  double ctauErrMax_calc;
  // new ErrMin
  for (int i = 0; i < hTot_M->GetNbinsX() / 2; i++)
  {
    if (hTot_M->GetBinContent(i) > 1)
    {
      ctauErrMin_calc = hTot_M->GetBinLowEdge(i + 1);
      break;
    }
  }
  // new ErrMax
  for (int i = 0; i < hTot_M->GetNbinsX(); i++)
  {
    if (hTot_M->GetBinContent(i) >= 1 && hTot_M->GetBinContent(i + 1) < 1 && hTot_M->GetBinContent(i + 2) < 1)
    {
      // ctauErrMax_calc = hSig->GetBinLowEdge(i)+hSig->GetBinWidth(i);
      ctauErrMax_calc = hTot_M->GetBinLowEdge(i) + hTot_M->GetBinWidth(i);
      break;
    }
    // else { ctauErrMax_calc = hSig->GetBinLowEdge(i); }
    else
    {
      ctauErrMax_calc = hTot_M->GetBinLowEdge(i);
    }
  }
  ctauErrMin_ = ctauErrMin_calc;
  ctauErrMax_ = ctauErrMax_calc;

  std::cout << "\nSet new bin width\n";
  double BinWidth = (ctauErrHigh - ctauErrLow) / nBins_;
  std::cout << "BinWidth : " << BinWidth << "\n";
  newBins_ = (ctauErrMax_ - ctauErrMin_) / BinWidth;
  std::cout << "newBins_ : " << newBins_ << "\n\n";

  std::cout << "Check number of events\n";
  std::cout << "ds_err: " << ws_err_->data("ds_err")->numEntries() << "\n";
  std::cout << "ds_splot: " << ws_err_->data("ds_splot")->numEntries() << "\n\n";

  TCanvas *c_err = new TCanvas("c_err", "ctau3D Err plot", 554, 4, 550, 520);
  c_err->cd();

  // ===== ctau3DErr dist. ===== //
  // buld pad1
  TPad *pad_err1 = new TPad("pad_err1", "pad_err1", 0, 0.16, 0.98, 1.0);
  pad_err1->SetTicks(1, 1);
  pad_err1->Draw();
  pad_err1->cd();
  err_frame1->SetTitle("");

  c_err->cd();
  c_err->SetLogy();

  pad_err1->cd();
  gPad->SetLogy();

  // draw dist. on pad1
  // keep the legacy Name("legacy") -> Need to track in legend and pull.
  RooPlot *err_frame2 = (RooPlot *)err_frame1->Clone("err_frame2");
  ws_err_->data("ds_err")->plotOn(err_frame2, Name("dataCTAUERR_Tot"), MarkerSize(.7), Binning(newBins_)); // Normalization(ws_errmc->ds_splot("reducedDS_MC")->sumEntries()
  ws_err_->pdf("err_tot_pdf")->plotOn(err_frame2, Name("err_tot_pdf"), LineColor(kGreen + 1), Range(ctauErrMin_, ctauErrMax_), LineWidth(2), Normalization(ws_err_->data("ds_err")->sumEntries(), RooAbsReal::NumEvent));
  ws_err_->data("ds_sig_err")->plotOn(err_frame2, Name("dataHist_Sig"), MarkerSize(.7), LineColor(kRed + 2), MarkerColor(kRed + 2), Binning(newBins_));
  ws_err_->pdf("err_sig_pdf")->plotOn(err_frame2, Name("err_sig_pdf"), LineColor(kRed + 2), LineWidth(2), Range(ctauErrMin_, ctauErrMax_));
  ws_err_->data("ds_bkg_err")->plotOn(err_frame2, Name("dataHist_Bkg"), MarkerSize(.7), LineColor(kBlue + 2), MarkerColor(kBlue + 2), Binning(newBins_));
  ws_err_->pdf("err_bkg_pdf")->plotOn(err_frame2, Name("err_bkg_pdf"), LineColor(kBlue + 2), LineWidth(2), Range(ctauErrMin_, ctauErrMax_));

  // set range of y axis
  Double_t YMax = hTot->GetBinContent(hTot->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= hTot->GetNbinsX(); i++)
    if (hTot->GetBinContent(i) > 0)
      YMin = std::min(YMin, hTot->GetBinContent(i));
  Double_t Yup(0.), Ydown(0.);
  Yup = YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.4 - 0.3)));
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.3 / (1.0 - 0.4 - 0.3))));
  err_frame2->GetYaxis()->SetRangeUser(Ydown, Yup);

  // ===== print lost events due to new ctau3DErr range ===== //
  std::cout << "final ctau3DErr range\n";
  std::cout << ws_err_->var("ctau3DErr")->getMin() << ", " << ws_err_->var("ctau3DErr")->getMax() << "\n\n";

  RooDataSet *ctauResCutDS = (RooDataSet *)ws_err_->data("ds_sig_err")->reduce(RooArgSet(*(ws_err_->var("ctau3DRes")), *(ws_err_->var("ctau3D")), *(ws_err_->var("ctau3DErr"))), Form("ctau3DErr>%f&&ctau3DErr<%f", ctauErrMin_, ctauErrMax_));
  ctauResCutDS->SetName("ctauResCutDS"); // used only for comparison
  ws_err_->import(*ctauResCutDS);
  Double_t outTot = ws_err_->data("ds_err")->numEntries();
  Double_t outRes = ws_err_->data("ctauResCutDS")->numEntries();
  std::cout << "Check how many events we lost due to new range cut\n";
  std::cout << "Total evt: " << outTot << "" << "\n";
  std::cout << "Residual evt: " << outRes << "" << "\n";
  std::cout << "lost evt: " << ((outTot - outRes) * 100) / outTot << " %\n\n";

  // ===== draw vertical lines to show ctau3DErr range ===== //
  TLine *minline = new TLine(ctauErrMin_, 0.0, ctauErrMin_, (Ydown * TMath::Power((Yup / Ydown), 0.5)));
  minline->SetLineStyle(2);
  minline->SetLineColor(1);
  minline->SetLineWidth(3);
  err_frame2->addObject(minline);
  TLine *maxline = new TLine(ctauErrMax_, 0.0, ctauErrMax_, (Ydown * TMath::Power((Yup / Ydown), 0.5)));
  maxline->SetLineStyle(2);
  maxline->SetLineColor(1);
  maxline->SetLineWidth(3);
  err_frame2->addObject(maxline);

  err_frame2->GetXaxis()->CenterTitle();
  err_frame2->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} Error (mm)");
  err_frame2->SetFillStyle(4000);
  err_frame2->GetYaxis()->SetTitleOffset(1.43);
  err_frame2->GetXaxis()->SetLabelSize(0);
  err_frame2->GetXaxis()->SetTitleSize(0);
  err_frame2->Draw();

  // ===== draw legends ===== //
  TLegend *leg_err = new TLegend(text_x + 0.5, text_y - 0.2, text_x + 0.7, text_y);
  leg_err->SetTextSize(text_size);
  leg_err->SetTextFont(43);
  leg_err->SetBorderSize(0);
  leg_err->AddEntry(err_frame2->findObject("dataCTAUERR_Tot"), "Data", "pe");
  leg_err->AddEntry(err_frame2->findObject("err_tot_pdf"), "Total PDF", "l");
  leg_err->AddEntry(err_frame2->findObject("err_sig_pdf"), "Signal", "l");
  leg_err->AddEntry(err_frame2->findObject("err_bkg_pdf"), "Background", "l");
  leg_err->Draw("same");
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow_, ptHigh_), text_x, text_y, text_color, text_size);
  if (yLow_ == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh_), text_x, text_y - y_diff, text_color, text_size);
  else if (yLow_ != 0)
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow_, yHigh_), text_x, text_y - y_diff, text_color, text_size);
  drawText(Form("Cent. %d - %d%s", cLow_ / 2, cHigh_ / 2, "%"), text_x, text_y - y_diff * 2, text_color, text_size);
  drawText(Form("%.2f < cos#theta_{EP} < %.2f", cosLow_, cosHigh_), text_x, text_y - y_diff * 3, text_color, text_size);
  drawText(Form("Loss: (%.4f%s) %.f evts", (outTot - outRes) * 100 / outTot, "%", outTot - outRes), text_x, text_y - y_diff * 4, text_color, text_size);

  // ===== pull distribution ===== //
  TPad *pull_pad = new TPad("pull_pad", "pull_pad", 0, 0.006, 0.98, 0.227);
  c_err->cd();
  pull_pad->Draw();
  pull_pad->cd();
  pull_pad->SetTopMargin(0); // Upper and lower plot are joined
  pull_pad->SetBottomMargin(0.67);
  pull_pad->SetBottomMargin(0.4);
  pull_pad->SetFillStyle(4000);
  pull_pad->SetFrameFillStyle(4000);
  pull_pad->SetTicks(1, 1);

  RooPlot *pull_tmp = (RooPlot *)err_frame2->Clone("pull_tmp");
  RooHist *pull_hist = pull_tmp->pullHist("dataCTAUERR_Tot", "err_tot_pdf");
  pull_hist->SetMarkerSize(0.8);
  RooPlot *pull_frame = ws_err_->var("ctau3DErr")->frame(Bins(nBins_), Range(minRange - 0.01, maxRange + 0.01));
  pull_frame->addPlotable(pull_hist, "PX");
  pull_frame->SetTitle("");
  pull_frame->SetTitleSize(0);
  pull_frame->GetYaxis()->SetTitleOffset(0.3);
  pull_frame->GetYaxis()->SetTitle("Pull");
  pull_frame->GetYaxis()->SetTitleSize(0.15);
  pull_frame->GetYaxis()->SetLabelSize(0.15);
  pull_frame->GetYaxis()->SetRangeUser(-3.8, 3.8);
  pull_frame->GetYaxis()->CenterTitle();

  pull_frame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} Error (mm)");
  pull_frame->GetXaxis()->SetTitleOffset(1.05);
  pull_frame->GetXaxis()->SetLabelOffset(0.04);
  pull_frame->GetXaxis()->SetLabelSize(0.15);
  pull_frame->GetXaxis()->SetTitleSize(0.15);
  pull_frame->GetXaxis()->CenterTitle();

  pull_frame->GetYaxis()->SetTickSize(0.04);
  pull_frame->GetYaxis()->SetNdivisions(404);
  pull_frame->GetXaxis()->SetTickSize(0.03);
  pull_frame->Draw();

  // draw horizontal line at center -> It's not a result of fit so pull points should be 0
  TLine *l_pull = new TLine(minRange - 0.01, 0, maxRange + 0.01, 0);
  l_pull->SetLineStyle(1);
  l_pull->Draw("same");
  pull_pad->Update();

  // ===== draw plots ===== //
  c_err->Update();

  c_err->SaveAs(("figs/" + outputPath_ + ".png").c_str());
  delete c_err;
}

void SPlotter::saveResults()
{
  TFile *outFile = new TFile(("roots/" + outputPath_ + ".root").c_str(), "RECREATE");

  ws_err_->Write();
  // ds_bkg_err->Write(); // inside ws_err
  // ds_sig_err->Write();
  // err_tot_pdf->Write();
  // err_bkg_pdf->Write();
  // err_sig_pdf->Write();

  outFile->Close();
  delete outFile;
}