#include "JpsiFitter.h"
#include "TFile.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include <iostream>
#include "RooCmdArg.h"
#include "TCanvas.h"
#include "TPad.h" // TPad, gPad
#include "TLine.h"
#include "TLegend.h"
#include "TMath.h" // min, Power

// sPlot
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "TH1D.h"


JpsiFitter::JpsiFitter(RooWorkspace &ws) : m_ws(ws) {}
JpsiFitter::~JpsiFitter() {}

void JpsiFitter::processTree(const std::string &filePath, const std::string &dsname, const std::string &obs, const std::string &cut){
  auto inputFile = std::unique_ptr<TFile>(TFile::Open(filePath.c_str()));
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Error: Can't open the input file: " << filePath << "\n";
    return;
  }


  RooDataSet *full_ds = (RooDataSet *)inputFile->Get("dataset");
  if (!full_ds) {
    std::cerr << "Error: Can't find 'dataset' in the inputfile" << "\n";
    return;
  }

  // apply cuts
  if (!cut.empty()) {
    RooDataSet *reduced_ds = (RooDataSet*)full_ds->reduce(cut.c_str());

    reduced_ds->SetName(dsname.c_str());
    m_ws.import(*reduced_ds);

    delete reduced_ds;
  } else {
    full_ds->SetName(dsname.c_str());
    m_ws.import(*full_ds);
  }
}

void JpsiFitter::print()
{
  // show the contents of workspace
  m_ws.Print("V");
}

void JpsiFitter::loadMassResult() {
  // get mass file
  std::string filePath = "roots/mass.root";

  auto inputFile = std::unique_ptr<TFile>(TFile::Open(filePath.c_str()));
  if (!inputFile || inputFile->IsZombie())
  {
    std::cerr << "Error: Can't open the input file: " << filePath << "\n";
    return;
  } else 
    std::cout << "SPlot: read mass fit results\n";

  // get workspace
  auto inputWs = (RooWorkspace*)inputFile->Get("wsMy");
  if (!inputWs) {
    std::cerr << "SPlot: Can't find 'wsMy' in the input file\n";
    return;
  }

  m_ws.import(*(inputWs->pdf("massModel")));
  m_ws.Print("V");
}

void JpsiFitter::doSplot() {
  // need orignal RooDataSet
  auto inputFile = std::unique_ptr<TFile>(TFile::Open("../../../files_roodata/RooDataSet_miniAOD_isMC0_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root"));
  inputFile->SetName("inputFile");
  RooDataSet *ds_tmp = (RooDataSet *)inputFile->Get("dataset");

  // nSig, nBkg from mass fit
  RooRealVar *sigYield = m_ws.var("nSig");
  RooRealVar *bkgYield = m_ws.var("nBkg");

  if (sigYield && bkgYield) {
    // set lower limit to 0 (SPlot requires it)
    sigYield->setRange(0, sigYield->getMax());
    bkgYield->setRange(0, bkgYield->getMax());
  } else {
    std::cerr << "SPlot: Can't find nSig or nBkg\n";
  }

  RooArgList yieldList;
  yieldList.add(*m_ws.var("nSig"));
  yieldList.add(*m_ws.var("nBkg"));

  std::cout << "Sig Yield: " << sigYield->getVal() << " +/- " << sigYield->getError() << "\n";
  std::cout << "Bkg Yield: " << bkgYield->getVal() << " +/- " << bkgYield->getError() << "\n";

  // Signal and Bkg PDFs of mass fit will be used
  // code below is messy, I keep it because it's legacy
  RooDataSet *ds_splot = (RooDataSet *)ds_tmp->Clone("ds_splot");
  // RooDataSet *ds_splot = (RooDataSet *)m_ws.data("ds_splot")->Clone("ds_splot");
  RooArgSet *cloneSet = (RooArgSet *)RooArgSet(*m_ws.pdf("massModel"), "massModel").snapshot(kTRUE);
  auto clone_mass_pdf = (RooAbsPdf *)cloneSet->find("massModel");
  clone_mass_pdf->setOperMode(RooAbsArg::ADirty, kTRUE);

  // ===== sPlot fit ===== //
  // sPlot: fds_splotit variable O1 to get the weightings of sig and bkg. We will aplly the weightings to variable O2 later.
  RooStats::SPlot sData = RooStats::SPlot("sData", "An SPlot", *ds_splot, clone_mass_pdf, yieldList);

  // the total yield and total weight should be same
  std::cout << "[INFO] Jpsi yield -> Mass Fit:" << m_ws.var("nSig")->getVal() << ", sWeights :" << sData.GetYieldFromSWeight("nSig") << "\n";
  std::cout << "[INFO] Bkg  yield -> Mass Fit:" << m_ws.var("nBkg")->getVal() << ", sWeights :" << sData.GetYieldFromSWeight("nBkg") << "\n";

  // individual events have different weights
  for (int i = 0; i < 10; i++)
  {
    // check weights of some events
    std::cout << "nSig weight: " << i << std::right << std::setw(12) << sData.GetSWeight(i, "nSig") << " /  nBkg weight: "
              << std::setw(12) << sData.GetSWeight(i, "nBkg") << " / Total Weight" << std::setw(12) << sData.GetSumOfEventSWeight(i)
              << "\n";
  }

  m_ws.import(*ds_splot);
}

void JpsiFitter::makeSplotPdfs() {
  // ===== build RooPlot and histograms ===== //
  // a lot of legacy codes

  double ctauErrLow = 0, ctauErrHigh = 10;
  int nBins = 100;

  std::cout << "SPlot: make SPlot PDFs\n";
  RooPlot *err_frame1 = m_ws.var("ctau3DErr")->frame(RooFit::Range(ctauErrLow, ctauErrHigh));

  TH1D *hTot = (TH1D *)m_ws.data("ds_splot")->createHistogram(("hTot"), *m_ws.var("ctau3DErr"), RooFit::Binning(err_frame1->GetNbinsX(), err_frame1->GetXaxis()->GetXmin(), err_frame1->GetXaxis()->GetXmax()));
  TH1D *hTot_M = (TH1D *)m_ws.data("ds_splot")->createHistogram(("hTot_M"), *m_ws.var("ctau3DErr"), RooFit::Binning(nBins, ctauErrLow, ctauErrHigh));

  // ===== build temporal sPlot weighted datasets ===== //
  std::cout << "SPlot: build temporal weighted datasets\n";
  RooDataSet *dataw_Sig_b = new RooDataSet("dataw_Sig_b", "TMP_SIG_DATA", (RooDataSet *)m_ws.data("ds_splot"),
                                           RooArgSet(*m_ws.var("ctau3DErr"), *m_ws.var("nSig_sw"), *m_ws.var("ctau3DRes"), *m_ws.var("ctau3D"), *m_ws.var("mass")), 0, "nSig_sw");
  TH1D *hSig = (TH1D *)dataw_Sig_b->createHistogram(("hSig"), *m_ws.var("ctau3DErr"), RooFit::Binning(nBins, ctauErrLow, ctauErrHigh));

  // Bkg
  RooDataSet *dataw_Bkg_b = new RooDataSet("dataw_Bkg_b", "TMP_BKG_DATA", (RooDataSet *)m_ws.data("ds_splot"),
                                           RooArgSet(*m_ws.var("ctau3DErr"), *m_ws.var("nBkg_sw"), *m_ws.var("ctau3DRes"), *m_ws.var("ctau3D"), *m_ws.var("mass")), 0, "nBkg_sw");
  TH1D *hBkg = (TH1D *)dataw_Bkg_b->createHistogram(("hBkg"), *m_ws.var("ctau3DErr"), RooFit::Binning(nBins, ctauErrLow, ctauErrHigh));

  // ===== set new ranges of ctau3DErr ===== //
  double ctauErrMin = 0;
  double ctauErrMax = 10;
  // new ErrMin
  for (int i = 0; i < hTot_M->GetNbinsX() / 2; i++)
  {
    if (hTot_M->GetBinContent(i) > 1)
    {
      ctauErrMin = hTot_M->GetBinLowEdge(i + 1);
      break;
    }
  }
  // new ErrMax
  for (int i = 0; i < hTot_M->GetNbinsX(); i++)
  {
    if (hTot_M->GetBinContent(i) >= 1 && hTot_M->GetBinContent(i + 1) < 1 && hTot_M->GetBinContent(i + 2) < 1)
    {
      // ctauErrMax = hSig->GetBinLowEdge(i)+hSig->GetBinWidth(i);
      ctauErrMax = hTot_M->GetBinLowEdge(i) + hTot_M->GetBinWidth(i);
      break;
    }
    // else { ctauErrMax = hSig->GetBinLowEdge(i); }
    else
    {
      ctauErrMax = hTot_M->GetBinLowEdge(i);
    }
  }

  std::cout << "\nSet new bin width\n";
  double BinWidth = (ctauErrHigh - ctauErrLow) / nBins;
  std::cout << "BinWidth : " << BinWidth << "\n";
  int newBins = (ctauErrMax - ctauErrMin) / BinWidth;
  std::cout << "newBins : " << newBins << "\n\n";


  // ===== build sPlot weighted datasets and RooHistPdfs ===== //
  // results of SPlot process
  TH1D *h_tot_test = (TH1D *)m_ws.data("ds_splot")->createHistogram(("h_tot_test"), *m_ws.var("ctau3DErr"), RooFit::Binning(newBins, ctauErrMin, ctauErrMax));
  // TH1D* h_tot_test = (TH1D*)m_ws.data("ds_splot")->createHistogram(("h_tot_test"), *m_ws.var("ctau3DErr"),RooFit::Binning(nBins,ctauErrMin,ctauErrMax));
  RooDataHist *totHist = new RooDataHist("ds_splot", "", *m_ws.var("ctau3DErr"), h_tot_test);
  RooHistPdf *err_tot_pdf = new RooHistPdf("err_tot_pdf", "Total RooHistPdf of ctauErr", *m_ws.var("ctau3DErr"), *totHist);

  // bkg ds and pdf
  TH1D *hBkg_w = (TH1D *)dataw_Bkg_b->createHistogram(("hBkg_w"), *m_ws.var("ctau3DErr"), RooFit::Binning(newBins, ctauErrMin, ctauErrMax));
  RooDataHist *bkgHist = new RooDataHist("ds_splot", "", *m_ws.var("ctau3DErr"), hBkg_w);
  RooHistPdf *err_bkg_pdf = new RooHistPdf("err_bkg_pdf", "hist pdf", *m_ws.var("ctau3DErr"), *bkgHist);

  // sig ds and pdf
  TH1D *hSig_w = (TH1D *)dataw_Sig_b->createHistogram(("hSig_w"), *m_ws.var("ctau3DErr"), RooFit::Binning(newBins, ctauErrMin, ctauErrMax));
  RooDataHist *sigHist = new RooDataHist("ds_splot", "", *m_ws.var("ctau3DErr"), hSig_w);
  RooHistPdf *err_sig_pdf = new RooHistPdf("err_sig_pdf", "hist pdf", *m_ws.var("ctau3DErr"), *sigHist);

  // Re-build datasets
  RooDataSet *ds_bkg_err = new RooDataSet("ds_bkg_err", "", (RooDataSet *)m_ws.data("ds_splot"),
                                          RooArgSet(*m_ws.var("ctau3DErr"), *m_ws.var("nBkg_sw"), *m_ws.var("ctau3DRes"), *m_ws.var("ctau3D"), *m_ws.var("mass")), 0, "nBkg_sw");
  RooDataSet *ds_sig_err = new RooDataSet("ds_sig_err", "", (RooDataSet *)m_ws.data("ds_splot"),
                                          RooArgSet(*m_ws.var("ctau3DErr"), *m_ws.var("nSig_sw"), *m_ws.var("ctau3DRes"), *m_ws.var("ctau3D"), *m_ws.var("mass")), 0, "nSig_sw");

  // ===== import sPlot results ===== //
  m_ws.import(*ds_sig_err);
  m_ws.import(*ds_bkg_err);
  m_ws.import(*err_tot_pdf);
  m_ws.import(*err_sig_pdf);
  m_ws.import(*err_bkg_pdf);


  // ===== this part should be moved to drawSplot() =====
  double minRange = (double)(floor(ctauErrMin * 100.) / 100.);
  double maxRange = (double)(ceil(ctauErrMax * 100.) / 100.);
  m_ws.var("ctau3DErr")->setRange("ctauErrWindow", ctauErrMin, ctauErrMax);
  err_frame1 = m_ws.var("ctau3DErr")->frame(RooFit::Bins(nBins), RooFit::Range(minRange - 0.01, maxRange + 0.01)); // modified

  std::cout << "Check number of events\n";
  std::cout << "ds_data: " << m_ws.data("ds_splot")->numEntries() << "\n";
  std::cout << "ds_splot: " << m_ws.data("ds_splot")->numEntries() << "\n\n";

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
  m_ws.data("ds_splot")->plotOn(err_frame2, RooFit::Name("dataCTAUERR_Tot"), RooFit::MarkerSize(.7), RooFit::Binning(newBins)); // Normalization(m_ws.->ds_splot("reducedDS_MC")->sumEntries()
  m_ws.pdf("err_tot_pdf")->plotOn(err_frame2, RooFit::Name("err_tot_pdf"), RooFit::LineColor(kGreen + 1), RooFit::Range(ctauErrMin, ctauErrMax), RooFit::LineWidth(2), RooFit::Normalization(m_ws.data("ds_splot")->sumEntries(), RooAbsReal::NumEvent));
  m_ws.data("ds_sig_err")->plotOn(err_frame2, RooFit::Name("dataHist_Sig"), RooFit::MarkerSize(.7), RooFit::LineColor(kRed + 2), RooFit::MarkerColor(kRed + 2), RooFit::Binning(newBins));
  m_ws.pdf("err_sig_pdf")->plotOn(err_frame2, RooFit::Name("err_sig_pdf"), RooFit::LineColor(kRed + 2), RooFit::LineWidth(2), RooFit::Range(ctauErrMin, ctauErrMax));
  m_ws.data("ds_bkg_err")->plotOn(err_frame2, RooFit::Name("dataHist_Bkg"), RooFit::MarkerSize(.7), RooFit::LineColor(kBlue + 2), RooFit::MarkerColor(kBlue + 2), RooFit::Binning(newBins));
  m_ws.pdf("err_bkg_pdf")->plotOn(err_frame2, RooFit::Name("err_bkg_pdf"), RooFit::LineColor(kBlue + 2), RooFit::LineWidth(2), RooFit::Range(ctauErrMin, ctauErrMax));

  // set range of y axis
  Double_t YMax = hTot->GetBinContent(hTot->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= hTot->GetNbinsX(); i++)
    if (hTot->GetBinContent(i) > 0)
      YMin = TMath::Min(YMin, hTot->GetBinContent(i));
  Double_t Yup(0.), Ydown(0.);
  Yup = YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.4 - 0.3)));
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.3 / (1.0 - 0.4 - 0.3))));
  err_frame2->GetYaxis()->SetRangeUser(Ydown, Yup);

  // ===== print lost events due to new ctau3DErr range ===== //
  std::cout << "final ctau3DErr range\n";
  std::cout << m_ws.var("ctau3DErr")->getMin() << ", " << m_ws.var("ctau3DErr")->getMax() << "\n\n";

  RooDataSet *ctauResCutDS = (RooDataSet *)ds_sig_err->reduce(RooArgSet(*(m_ws.var("ctau3DRes")), *(m_ws.var("ctau3D")), *(m_ws.var("ctau3DErr"))), Form("ctau3DErr>%f&&ctau3DErr<%f", ctauErrMin, ctauErrMax));
  ctauResCutDS->SetName("ctauResCutDS"); // used only for comparison
  m_ws.import(*ctauResCutDS);
  Double_t outTot = m_ws.data("ds_splot")->numEntries();
  Double_t outRes = m_ws.data("ctauResCutDS")->numEntries();
  std::cout << "Check how many events we lost due to new range cut\n";
  std::cout << "Total evt: " << outTot << "\n";
  std::cout << "Residual evt: " << outRes << "\n";
  std::cout << "lost evt: " << ((outTot - outRes) * 100) / outTot << " %\n\n";

  // ===== draw vertical lines to show ctau3DErr range ===== //
  TLine *minline = new TLine(ctauErrMin, 0.0, ctauErrMin, (Ydown * TMath::Power((Yup / Ydown), 0.5)));
  minline->SetLineStyle(2);
  minline->SetLineColor(1);
  minline->SetLineWidth(3);
  err_frame2->addObject(minline);
  TLine *maxline = new TLine(ctauErrMax, 0.0, ctauErrMax, (Ydown * TMath::Power((Yup / Ydown), 0.5)));
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
  // TLegend *leg_err = new TLegend(text_x + 0.5, text_y - 0.2, text_x + 0.7, text_y);
  // leg_err->SetTextSize(text_size);
  // leg_err->SetTextFont(43);
  // leg_err->SetBorderSize(0);
  // leg_err->AddEntry(err_frame2->findObject("dataCTAUERR_Tot"), "Data", "pe");
  // leg_err->AddEntry(err_frame2->findObject("err_tot_pdf"), "Total PDF", "l");
  // leg_err->AddEntry(err_frame2->findObject("err_sig_pdf"), "Signal", "l");
  // leg_err->AddEntry(err_frame2->findObject("err_bkg_pdf"), "Background", "l");
  // leg_err->Draw("same");
  // drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x, text_y, text_color, text_size);
  // if (yLow == 0)
  //   drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x, text_y - y_diff, text_color, text_size);
  // else if (yLow != 0)
  //   drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x, text_y - y_diff, text_color, text_size);
  // drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x, text_y - y_diff * 2, text_color, text_size);
  // drawText(Form("%.2f < cos#theta_{EP} < %.2f", cos_low, cos_high), text_x, text_y - y_diff * 3, text_color, text_size);
  // drawText(Form("Loss: (%.4f%s) %.f evts", (outTot - outRes) * 100 / outTot, "%", outTot - outRes), text_x, text_y - y_diff * 4, text_color, text_size);

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
  RooPlot *pull_frame = m_ws.var("ctau3DErr")->frame(RooFit::Bins(nBins), RooFit::Range(minRange - 0.01, maxRange + 0.01));
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
  c_err->SaveAs("figs/splot.png");
}

void JpsiFitter::drawSplot() {
  
}