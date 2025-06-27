#include "MassFitter.h"

#include <iostream>
#include <algorithm>
#include <sstream>
#include "TStopwatch.h"
#include "TSystem.h"
#include "TFile.h"
#include "TROOT.h" // gROOT
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TH1.h"
#include "TMath.h"

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooChebychev.h"
#include "RooHist.h"
#include "RooMsgService.h"

#include "../../../headers/polarizationUtilities.h"
using namespace RooFit;


MassFitter::MassFitter(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh, float cosLow, float cosHigh)
    : ptLow_(ptLow), ptHigh_(ptHigh), yLow_(yLow), yHigh_(yHigh), cLow_(cLow), cHigh_(cHigh), cosLow_(cosLow), cosHigh_(cosHigh)
{
  // RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
  // RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  std::cout << "MassFitter object is created." << std::endl;
}

MassFitter::~MassFitter()
{
  delete dataFile_;
  delete mcFitFile_;
  delete ws_mass_;
  delete fit_mass_;

  // model Parameters
  delete mean_mass_;
  delete sigma1_mass_;
  delete sigma2_mass_;
  delete alpha1_mass_;
  delete power1_mass_;
  delete fracG1_mass_;
  delete sl1_mass_;
  delete sl2_mass_;
  delete nSig_mass_;
  delete nBkg_mass_;

  std::cout << "MassFitter object is destroyed." << std::endl;
}

void MassFitter::run()
{
  TStopwatch t;
  t.Start();
  std::cout << "\n ===== Start Data mass fit =====\n";

  setupInputOutput();
  prepareDataset();
  loadMcFitResults();
  buildModel();
  performFit();
  plotResults();
  saveResults();

  std::cout << "\n ===== Finish Data mass fit =====\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}

void MassFitter::setupInputOutput()
{
  // output path
  std::ostringstream out_path;
  out_path << "mass/mass_pt" << ptLow_ << "-" << ptHigh_
           << "_cent" << cLow_ << "-" << cHigh_ << "_absY" << yLow_ << "-" << yHigh_
           << "_cos" << cosLow_ << "_" << cosHigh_;
  outputPath_ = out_path.str();

  // build mc fit input path
  std::ostringstream mc_in_path;
  mc_in_path << "roots/mc_mass/mc_mass_pt" << ptLow_ << "-" << ptHigh_
             << "_cent" << cLow_ << "-" << cHigh_ << "_absY" << yLow_ << "-" << yHigh_
             << "_cos" << cosLow_ << "_" << cosHigh_
             << ".root";
  mcFitInputPath_ = mc_in_path.str();

  // output folder
  gSystem->mkdir("roots/mass/", kTRUE);
  gSystem->mkdir("figs/mass", kTRUE);
}

void MassFitter::prepareDataset()
{
  dataFile_ = new TFile(dataInputPath_.c_str());
  RooDataSet *ds_tmp = (RooDataSet *)dataFile_->Get("dataset");
  RooDataSet *ds_tmp_weight = new RooDataSet("ds_tmp_weight", "", *ds_tmp->get(), Import(*ds_tmp), WeightVar("weight"));

  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f && cBin>=%d && cBin<%d", ptLow_, ptHigh_, yLow_, yHigh_, massLow_, massHigh_, cLow_, cHigh_);
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";
  TString OS = "recoQQsign==0 &&";
  TString angle_cut = Form("&& cos_ep>%.2f && cos_ep<%.2f", cosLow_, cosHigh_);
  TString final_cut = OS + accCut + kineCut + angle_cut;

  ds_mass_ = (RooDataSet *)ds_tmp_weight->reduce(final_cut.Data());
  ds_mass_->SetName("ds_mass_");

  mass_ = (RooRealVar *)ds_mass_->get(0)->find("mass");
  mass_->setRange(massLow_, massHigh_);
}

void MassFitter::loadMcFitResults()
{
  mcFitFile_ = new TFile(mcFitInputPath_.c_str());
  RooWorkspace *ws_mc = (RooWorkspace *)mcFitFile_->Get("ws_mc");

  // deep-cloning MC parameters
  mean_mass_ = (RooRealVar *)ws_mc->var("mean_mc_")->Clone("mean_mass_");
  sigma1_mass_ = (RooRealVar *)ws_mc->var("sigma1_mc_")->Clone("sigma1_mass_");
  sigma2_mass_ = (RooRealVar *)ws_mc->var("sigma2_mc_")->Clone("sigma2_mass_");
  alpha1_mass_ = (RooRealVar *)ws_mc->var("alpha1_mc_")->Clone("alpha1_mass_");
  power1_mass_ = (RooRealVar *)ws_mc->var("power1_mc_")->Clone("power1_mass_");
  fracG1_mass_ = (RooRealVar *)ws_mc->var("fracG1_mc_")->Clone("fracG1_mass_");

  alpha1_mass_->setConstant(kTRUE);
  power1_mass_->setConstant(kTRUE);
  fracG1_mass_->setConstant(kTRUE);
}

void MassFitter::buildModel()
{
  ws_mass_ = new RooWorkspace("ws_mass_");

  // signal model
  auto G1Sig = new RooGaussian("G1Sig", "", *mass_, *mean_mass_, *sigma1_mass_);
  auto CB1Sig = new RooCBShape("CB1Sig", "", *mass_, *mean_mass_, *sigma2_mass_, *alpha1_mass_, *power1_mass_);
  auto mass_sig_ = new RooAddPdf("mass_sig_", "", RooArgList(*G1Sig, *CB1Sig), RooArgList(*fracG1_mass_));

  // bkg model
  sl1_mass_ = new RooRealVar("sl1_mass_", "sl1", 0.1, -1., 1.);
  sl2_mass_ = new RooRealVar("sl2_mass_", "sl2", 0.01, -1., 1.);
  auto mass_bkg_ = new RooChebychev("mass_bkg_", "Background", *mass_, RooArgList(*sl1_mass_, *sl2_mass_));

  // total model
  nSig_mass_ = new RooRealVar("nSig_mass_", "", 5000, 1, 250000);
  nBkg_mass_ = new RooRealVar("nBkg_mass_", "", 20000, 1, 1000000);
  auto mass_pdf_ = new RooAddPdf("mass_pdf_", "", RooArgList(*mass_sig_, *mass_bkg_), RooArgList(*nSig_mass_, *nBkg_mass_));

  ws_mass_->import(*mass_pdf_);
}

void MassFitter::performFit()
{
  std::cout << "Is weighted: " << ds_mass_->isWeighted() << std::endl;
  fit_mass_ = ws_mass_->pdf("mass_pdf_")->fitTo(*ds_mass_, Extended(kTRUE), Save(kTRUE), AsymptoticError(ds_mass_->isWeighted()), NumCPU(12), Strategy(2), PrintLevel(-1));
  fit_mass_->Print("v");
}

void MassFitter::saveResults()
{
  ws_mass_->var("mean_mass_")->setConstant(kTRUE);
  ws_mass_->var("sigma1_mass_")->setConstant(kTRUE);
  ws_mass_->var("sigma2_mass_")->setConstant(kTRUE);
  ws_mass_->var("sl1_mass_")->setConstant(kTRUE);
  ws_mass_->var("sl2_mass_")->setConstant(kTRUE);

  TFile *outFile = new TFile(("roots/" + outputPath_ + ".root").c_str(), "RECREATE");

  TCanvas *canvas = (TCanvas *)gROOT->FindObject("c_mass");
  if (canvas)
    canvas->Write();
  if (fit_mass_)
    fit_mass_->Write();
  if (ws_mass_) {
    // nSig_mass_->SetName("nSig_mass");
    // nBkg_mass_->SetName("nBkg_mass");
    ws_mass_->import(*ds_mass_);
    ws_mass_->Write();
  }
  
  outFile->Close();
  delete outFile;
}

void MassFitter::plotResults()
{
  TCanvas *c_mass = new TCanvas("c_mass", "", 800, 800);
  c_mass->cd();
  TPad *massPad = new TPad("massPad", "massPad", 0, 0.25, 0.98, 1.0);
  massPad->SetTicks(1, 1);
  massPad->Draw();
  massPad->cd();
  gPad->SetLogy();

  // mass dist.
  // RooPlot *massFrame = (RooPlot *)mass_->frame(Range(2.6, 3.5));
  RooPlot *massFrame = (RooPlot *)mass_->frame();
  ds_mass_->plotOn(massFrame, Name("ds_mass_"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7));
  ws_mass_->pdf("mass_pdf_")->plotOn(massFrame, Name("mass_pdf_"), LineColor(kBlack));

  ws_mass_->pdf("mass_pdf_")->plotOn(massFrame, LineStyle(kDashed), Components("mass_bkg_,CB1Sig"), Name("CB1Sig"), LineColor(44), LineWidth(2));

  ws_mass_->pdf("mass_pdf_")->plotOn(massFrame, LineStyle(kDashed), Components("mass_bkg_,G1Sig"), Name("G1Sig"), LineColor(8), LineWidth(2));

  ws_mass_->pdf("mass_pdf_")->plotOn(massFrame, LineStyle(kDashed), Components("mass_bkg_"), Name("mass_bkg"), LineColor(kBlue + 2), LineWidth(2));


  // set y plot range
  // RooPlot *massFrame = mass_->frame(nMassBin); // bins
  
  TH1 *h_tmp = ds_mass_->createHistogram("h_tmp", *mass_, Binning(massFrame->GetNbinsX(), massFrame->GetXaxis()->GetXmin(), massFrame->GetXaxis()->GetXmax()));

  Double_t YMax = h_tmp->GetBinContent(h_tmp->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= h_tmp->GetNbinsX(); i++)
    if (h_tmp->GetBinContent(i) > 0)
      YMin = std::min(YMin, h_tmp->GetBinContent(i));
  double Ydown;
  double Yup;
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.1 / (1.0 - 0.1 - 0.4))));
  Yup = YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.1 - 0.4)));
  // massFrame->SetMaximum(1e6);
  // massFrame->SetMinimum(20); // yLog can't get 0 as minimum
  massFrame->GetYaxis()->SetRangeUser(Ydown, Yup);

  massFrame->SetFillStyle(4000);
  massFrame->GetYaxis()->SetTitleOffset(1.43);
  massFrame->GetXaxis()->SetLabelSize(0); // to hide x-axis label
  massFrame->GetXaxis()->SetTitleSize(0);
  // massFrame->GetXaxis()->CenterTitle();
  // massFrame->GetXaxis()->SetRangeUser(massLow, massHigh);
  // massFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  massFrame->Draw();
  
  // ===== draw legends for PDFs and dataset ===== //
  TLegend *leg_pdfs = new TLegend(text_x + 0.45, text_y - 0.2, text_x + 0.65, text_y);
  leg_pdfs->SetTextSize(text_size);
  leg_pdfs->SetTextFont(43);
  leg_pdfs->SetBorderSize(0);
  leg_pdfs->AddEntry(massFrame->findObject("ds_mass_"), "Data", "pe");
  leg_pdfs->AddEntry(massFrame->findObject("mass_pdf_"), "Total", "l");
  leg_pdfs->AddEntry(massFrame->findObject("CB1Sig"), "CB1+Bkg", "l");
  leg_pdfs->AddEntry(massFrame->findObject("G1Sig"), "Gauss1+Bkg", "l");
  leg_pdfs->AddEntry(massFrame->findObject("mass_bkg"), "2nd Chebychev", "l");
  leg_pdfs->Draw("same");

  // ws_mass_->var(;

  // ===== draw legends ===== //
  // drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c ; Cent. %d - %d%s; %.2f < cos#theta_{EP} < %.2f", ptLow_, ptHigh_, cLow_ / 2, cHigh_ / 2, "%", cosLow_, cosHigh_), text_x, text_y, text_color, text_size);
  // if (yLow_ == 0)
  //   drawText(Form("|y^{#mu#mu}| < %.1f", yHigh_), text_x, text_y - y_diff, text_color, text_size);
  // else if (yLow_ != 0)
  //   drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow_, yHigh_), text_x, text_y - y_diff, text_color, text_size);
  // // drawText(Form("Cent. %d - %d%s", cLow/2, cHigh/2, "%"),text_x,text_y-y_diff*2,text_color,text_size);
  // drawText(Form("N_{J/#psi} = %.f #pm %.f; N_{Bkg} = %.f #pm %.f", ws_mass_->var("n_sig")->getVal(), ws_mass_->var("n_sig")->getError(), ws_mass_->var("n_bkg")->getVal(), ws_mass_->var("n_bkg")->getError()), text_x, text_y - y_diff * 2, text_color, text_size);
  // drawText(Form("m_{J/#psi} = %.4f #pm %.4f", ws_mass_->var("mass_mean")->getVal(), ws_mass_->var("mass_mean")->getError()), text_x, text_y - y_diff * 3, text_color, text_size);
  // drawText(Form("#alpha_{J/#psi} = %.4f (fixed); n_{J/#psi} = %.4f (fixed)", ws_mass_->var("mass_alpha1")->getVal(), ws_mass_->var("mass_power1")->getVal()), text_x, text_y - y_diff * 4, text_color, text_size);
  // drawText(Form("f_{G1} = %.4f (fixed)", fracG1_mass_->getVal()), text_x, text_y - y_diff * 5, text_color, text_size);
  // drawText(Form("#sigma1_{J/#psi} = %.2f #pm %.2f MeV/c^{2}", ws_mass_->var("mass_sigma1")->getVal() * 1000, ws_mass_->var("mass_sigma1")->getError() * 1000), text_x, text_y - y_diff * 6, text_color, text_size);
  // drawText(Form("#sigma2_{J/#psi} = %.2f #pm %.2f MeV/c^{2}", ws_mass_->var("mass_sigma2")->getVal() * 1000, ws_mass_->var("mass_sigma2")->getError() * 1000), text_x, text_y - y_diff * 7, text_color, text_size);
  // drawText(Form("slope1_{Bkg} = %.2f #pm %.2f MeV/c^{2}", (ws_mass_->var("mass_sl1")->getVal()), (ws_mass_->var("mass_sl1")->getError())), text_x, text_y - y_diff * 8, text_color, text_size);
  // drawText(Form("slope2_{Bkg} = %.2f #pm %.2f MeV/c^{2}", (ws_mass_->var("mass_sl2")->getVal()), (ws_mass_->var("mass_sl2")->getError())), text_x, text_y - y_diff * 9, text_color, text_size);

  // ===== draw pull dist. ===== //
  TPad *pullPad = new TPad("pullPad", "", 0, 0.001, 0.98, 0.32);
  c_mass->cd();
  pullPad->Draw();
  pullPad->cd();
  pullPad->SetTopMargin(0); // Upper and lower plot are joined
  pullPad->SetBottomMargin(0.67);
  pullPad->SetBottomMargin(0.4);
  pullPad->SetFillStyle(4000);
  pullPad->SetFrameFillStyle(4000);
  pullPad->SetTicks(1, 1);

  auto pullHist = massFrame->pullHist("ds_mass_", "mass_pdf_", true); // true: isExtended?
  auto pullFrame = mass_->frame(Title(";;Pull Distribution"), Range(2.6, 3.5));
  pullHist->SetMarkerSize(0.8);
  pullFrame->addPlotable(pullHist, "P");

  // cosmetics
  pullFrame->SetTitle("");
  pullFrame->SetTitleSize(0);
  pullFrame->GetYaxis()->SetTitleOffset(0.3);
  pullFrame->GetYaxis()->SetTitle("Pull");
  pullFrame->GetYaxis()->SetTitleSize(0.08);
  pullFrame->GetYaxis()->SetLabelSize(0.08);
  pullFrame->GetYaxis()->SetRangeUser(-9, 9);
  pullFrame->GetYaxis()->CenterTitle();

  pullFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  pullFrame->GetXaxis()->SetTitleOffset(1.55);
  pullFrame->GetXaxis()->SetLabelOffset(0.04);
  pullFrame->GetXaxis()->SetLabelSize(0.08);
  pullFrame->GetXaxis()->SetTitleSize(0.08);
  pullFrame->GetXaxis()->CenterTitle();

  pullFrame->GetYaxis()->SetTickSize(0.04);
  pullFrame->GetYaxis()->SetNdivisions(404);
  pullFrame->GetXaxis()->SetTickSize(0.03);
  pullFrame->Draw();

  // ===== horizontoal lines on pull hist ===== //
  // maybe we can make function for them...
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

  // negative
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

  c_mass->Draw();
  c_mass->SaveAs(("figs/" + outputPath_ + ".png").c_str());
}