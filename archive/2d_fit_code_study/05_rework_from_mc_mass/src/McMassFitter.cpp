#include "McMassFitter.h"

#include <iostream>
#include <algorithm> // std::min
#include <sstream>
#include "TStopwatch.h"
#include "TSystem.h"
#include "TFile.h"
// #include "TROOT.h" // gROOT
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TH1.h"
#include "TMath.h"

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooHist.h"
#include "RooMsgService.h"

#include "../../../headers/polarizationUtilities.h"

using namespace RooFit;

McMassFitter::McMassFitter(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh, float cosLow, float cosHigh)
    : ptLow_(ptLow), ptHigh_(ptHigh), yLow_(yLow), yHigh_(yHigh), cLow_(cLow), cHigh_(cHigh), cosLow_(cosLow), cosHigh_(cosHigh)
{
  // Reduce RooFit message
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Caching);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Integration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::InputArguments);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  std::cout << "McMassFitter object is created." << std::endl;
}

McMassFitter::~McMassFitter()
{
  delete f_mc_;
  
  delete mass_pdf_mc_;
  delete fit_result_mc_;

  delete mean_mc_;
  delete sigma1_mc_;
  delete sigma2_mc_;
  delete alpha1_mc_;
  delete power1_mc_;
  delete fracG1_mc_;
  delete nSig_mc_;

  // don' touch ds_tmp_, ds_mc_, mass_mc_
  std::cout << "McMassFitter object is destroyed." << std::endl;
}

void McMassFitter::run()
{
  TStopwatch t;
  t.Start();
  std::cout << "\n ===== Start MC mass fit =====\n";

  setupInputOutput();
  prepareDataset();
  buildModel();
  performFit();
  plotResults();
  saveResults();

  std::cout << "\n ===== Finish MC mass fit =====\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime()); 
}

void McMassFitter::setupInputOutput()
{
  std::ostringstream out_file_path;
  out_file_path << "mc_mass/mc_mass_pt" << ptLow_ << "-" << ptHigh_
                << "_cent" << cLow_ << "-" << cHigh_ << "_absY" << yLow_ << "-" << yHigh_
                << "_cos" << cosLow_ << "_" << cosHigh_;
  outpath_ = out_file_path.str();

  gSystem->mkdir("roots/mc_mass/", kTRUE);
  gSystem->mkdir("figs/mc_mass", kTRUE);
}

void McMassFitter::prepareDataset()
{
  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f && cBin>=%d && cBin<%d", ptLow_, ptHigh_, yLow_, yHigh_, massLow_, massHigh_, cLow_, cHigh_);
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";
  TString OS = "recoQQsign==0 &&";
  TString angle_cut = Form("&& cos_ep>%.2f && cos_ep<%.2f", cosLow_, cosHigh_);
  TString final_cut = OS + accCut + kineCut + angle_cut;

  f_mc_ = new TFile(inputPath_.c_str());
  ds_tmp_ = (RooDataSet *)f_mc_->Get("dataset");

  ds_mc_ = (RooDataSet *)ds_tmp_->reduce(final_cut.Data());
  ds_mc_->SetName("ds_mc_");

  mass_mc_ = (RooRealVar *)ds_mc_->get(0)->find("mass");
  mass_mc_->setRange(massLow_, massHigh_);
  mass_mc_->setRange("mcMassFit", massFitMin_, massFitMax_);
}

void McMassFitter::buildModel()
{
  mean_mc_ = new RooRealVar("mean_mc_", "", 3.096, 3.086, 3.106);
  sigma1_mc_ = new RooRealVar("sigma1_mc_", "", 0.01, 0.001, 0.1);
  sigma2_mc_ = new RooRealVar("sigma2_mc_", "", 0.03, 0.001, 0.1);
  alpha1_mc_ = new RooRealVar("alpha1_mc_", "", 2, 1, 5);
  power1_mc_ = new RooRealVar("power1_mc_", "", 2, 1, 5);
  fracG1_mc_ = new RooRealVar("fracG1_mc_", "", 0.6, 0.05, 0.95);

  auto mc_G1Sig = new RooGaussian("mc_G1Sig", "", *mass_mc_, *mean_mc_, *sigma1_mc_);
  auto mc_CB1Sig = new RooCBShape("mc_CB1Sig", "", *mass_mc_, *mean_mc_, *sigma2_mc_, *alpha1_mc_, *power1_mc_);

  auto mc_G1CB1Sig = new RooAddPdf("mc_G1CB1Sig", "", RooArgList(*mc_G1Sig, *mc_CB1Sig), RooArgList(*fracG1_mc_));

  nSig_mc_ = new RooRealVar("nSig_mc_", "", 20000, 1, 50000);
  mass_pdf_mc_ = new RooAddPdf("mass_pdf_mc_", "", RooArgList(*mc_G1CB1Sig), RooArgList(*nSig_mc_));
}

void McMassFitter::performFit()
{
  std::cout << "Is weighted: " << ds_mc_->isWeighted() << std::endl;
  fit_result_mc_ = mass_pdf_mc_->fitTo(*ds_mc_, Extended(kTRUE), Save(kTRUE), SumW2Error(kFALSE), NumCPU(24), Strategy(2), Range("mcMassFit"), Timer(kTRUE), PrintLevel(-1));
  fit_result_mc_->Print("V");
}

void McMassFitter::saveResults()
{
  TFile *outFile = new TFile(("roots/" + outpath_ + ".root").c_str(), "RECREATE");
  auto ws_mc = new RooWorkspace("ws_mc");
  mass_pdf_mc_->SetName("mass_pdf_mc");
  ws_mc->import(*mass_pdf_mc_);

  ws_mc->Write();
  fit_result_mc_->Write();

  outFile->Close();
  delete outFile;
}

void McMassFitter::plotResults()
{
  TCanvas *c_mc = new TCanvas("c_mc_mass", "", 800, 800);
  c_mc->cd();
  TPad *massPad = new TPad("massPad", "massPad", 0, 0.25, 0.98, 1.0);
  massPad->SetTicks(1, 1);
  massPad->Draw();
  massPad->cd();
  gPad->SetLogy();

  // ===== draw mass dist. ===== //
  
  auto massFrame = (RooPlot *)mass_mc_->frame(Title(""), Bins(nMassBin_), Range(massLow_, massHigh_));
  
  ds_mc_->plotOn(massFrame, Name("ds_mc_"), DataError(RooAbsData::SumW2), MarkerSize(.7));
  mass_pdf_mc_->plotOn(massFrame, Name("mass_pdf_mc_"), Range("mcMassFit"), NormRange("mcMassFit"), LineColor(kBlack));
  mass_pdf_mc_->plotOn(massFrame, LineStyle(kDashed), Components("mc_CB1Sig"), Name("mc_CB1Sig"), LineColor(44), LineWidth(2));
  mass_pdf_mc_->plotOn(massFrame, LineStyle(kDashed), Components("mc_G1Sig"), Name("mc_G1Sig"), LineColor(8), LineWidth(2));

  // set y plot range
  // RooPlot *massFrame = mass_mc_->frame(nMassBin); // bins
  TH1 *h_tmp = ds_mc_->createHistogram("hist", *mass_mc_, Binning(massFrame->GetNbinsX(), massFrame->GetXaxis()->GetXmin(), massFrame->GetXaxis()->GetXmax()));
  Double_t YMax = h_tmp->GetBinContent(h_tmp->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= h_tmp->GetNbinsX(); i++)
    if (h_tmp->GetBinContent(i) > 0)
      YMin = std::min(YMin, h_tmp->GetBinContent(i));
  
  double Ydown;
  double Yup;
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.1 / (1.0 - 0.1 - 0.4))));
  Yup = YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.1 - 0.4)));
  if (YMin <= 0)
    YMin = 0.01;
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
  // mass_pdf_mc_->plotOn(massFrame, Name("pdfMASS_Tot"), LineColor(kBlack));
  TLegend *leg_pdfs = new TLegend(text_x + 0.45, text_y - 0.2, text_x + 0.65, text_y);
  leg_pdfs->SetTextSize(text_size);
  leg_pdfs->SetTextFont(43);
  leg_pdfs->SetBorderSize(0);
  leg_pdfs->AddEntry(massFrame->findObject("ds_mc_"), "MC", "pe");
  leg_pdfs->AddEntry(massFrame->findObject("mass_pdf_mc_"), "Total", "l");
  leg_pdfs->AddEntry(massFrame->findObject("mc_CB1Sig"), "CB1", "l");
  leg_pdfs->AddEntry(massFrame->findObject("mc_G1Sig"), "Gauss1", "l");
  leg_pdfs->Draw("same");

  // ===== draw legends ===== //
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c ; Cent. %d - %d%s; %.2f < cos#theta_{EP} < %.2f", ptLow_, ptHigh_, cLow_ / 2, cHigh_ / 2, "%", cosLow_, cosHigh_), text_x, text_y, text_color, text_size);
  if (yLow_ == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh_), text_x, text_y - y_diff, text_color, text_size);
  else if (yLow_ != 0)
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow_, yHigh_), text_x, text_y - y_diff, text_color, text_size);
  // drawText(Form("Cent. %d - %d%s", cLow/2, cHigh/2, "%"),text_x,text_y-y_diff*2,text_color,text_size);
  drawText(Form("N_{J/#psi} = %.f #pm %.f", nSig_mc_->getVal(), nSig_mc_->getError()), text_x, text_y - y_diff * 2, text_color, text_size);
  drawText(Form("m_{J/#psi} = %.4f #pm %.4f", mean_mc_->getVal(), mean_mc_->getError()), text_x, text_y - y_diff * 3, text_color, text_size);
  drawText(Form("#alpha_{J/#psi} = %.4f #pm %.4f", alpha1_mc_->getVal(), alpha1_mc_->getError()), text_x, text_y - y_diff * 4, text_color, text_size);
  drawText(Form("f_{G1} = %.4f #pm %.4f", fracG1_mc_->getVal(), fracG1_mc_->getError()), text_x, text_y - y_diff * 5, text_color, text_size);
  drawText(Form("n_{J/#psi} = %.4f #pm %.4f", power1_mc_->getVal(), power1_mc_->getError()), text_x, text_y - y_diff * 6, text_color, text_size);
  drawText(Form("#sigma1_{J/#psi} = %.2f #pm %.2f MeV/c^{2}", (sigma1_mc_->getVal()) * 1000, (sigma1_mc_->getError()) * 1000), text_x, text_y - y_diff * 7, text_color, text_size);
  drawText(Form("#sigma2_{J/#psi} = %.2f #pm %.2f MeV/c^{2}", (sigma2_mc_->getVal()) * 1000, (sigma2_mc_->getError()) * 1000), text_x, text_y - y_diff * 8, text_color, text_size);

  
  // ===== draw pull dist. ===== //
  TPad *pullPad = new TPad("pullPad", "", 0, 0.001, 0.98, 0.32);
  c_mc->cd();
  pullPad->Draw();
  pullPad->cd();
  pullPad->SetTopMargin(0); // Upper and lower plot are joined
  pullPad->SetBottomMargin(0.67);
  pullPad->SetBottomMargin(0.4);
  pullPad->SetFillStyle(4000);
  pullPad->SetFrameFillStyle(4000);
  pullPad->SetTicks(1, 1);

  auto massPull = massFrame->pullHist("ds_mc_", "mass_pdf_mc_", true);
  // The error (TAxis::TAxis::SetRangeUser:0: RuntimeWarning: ulast > fXmax, fXmax is used) appears
  // because range of ds_mc_ is wider than mass_pdf_mc_(using massFitRange)
  auto pullFrame = mass_mc_->frame(Title(";;Pull Distribution"), Range(massLow_, massHigh_));
  massPull->SetMarkerSize(0.8);
  pullFrame->addPlotable(massPull, "P");

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

  c_mc->Draw();
  c_mc->SaveAs(("figs/" + outpath_ + ".png").c_str());

  delete c_mc;
}