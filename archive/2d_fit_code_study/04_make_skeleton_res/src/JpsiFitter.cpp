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

void JpsiFitter::processTree(const std::string &filePath, const std::string &dsname, const std::string &obs, const std::string &cut)
{
  auto inputFile = std::unique_ptr<TFile>(TFile::Open(filePath.c_str()));
  if (!inputFile || inputFile->IsZombie())
  {
    std::cerr << "Error: Can't open the input file: " << filePath << "\n";
    return;
  }

  RooDataSet *full_ds = (RooDataSet *)inputFile->Get("dataset");
  if (!full_ds)
  {
    std::cerr << "Error: Can't find 'dataset' in the inputfile" << "\n";
    return;
  }

  // apply cuts
  if (!cut.empty()) {
    RooDataSet *reduced_ds = (RooDataSet *)full_ds->reduce(cut.c_str());

    reduced_ds->SetName(dsname.c_str());
    m_ws.import(*reduced_ds);

    delete reduced_ds;
  }
  else {
    full_ds->SetName(dsname.c_str());
    m_ws.import(*full_ds);
  }
}

void JpsiFitter::print()
{
  // show the contents of workspace
  m_ws.Print("V");
}

void JpsiFitter::loadMassResult(const std::string &filePath)
{
  auto inputFile = std::unique_ptr<TFile>(TFile::Open(filePath.c_str()));
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Error: Can't open the input file: " << filePath << "\n";
    return;
  } else
    std::cout << "SPlot: read mass fit results\n";

  // get workspace
  auto inputWs = (RooWorkspace *)inputFile->Get("wsMy");
  if (!inputWs) {
    std::cerr << "SPlot: Can't find 'wsMy' in the input file\n";
    return;
  }

  m_ws.import(*(inputWs->pdf("massModel")));
  m_ws.Print("V");
}

void JpsiFitter::doSplot(const std::string &filePath, const std::string &cut)
{
  // need orignal RooDataSet
  auto inputFile = std::unique_ptr<TFile>(TFile::Open(filePath.c_str()));
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
  RooDataSet *ds_splot_tmp = (RooDataSet *)ds_tmp->Clone("ds_splot");
  RooDataSet *ds_splot = (RooDataSet *)ds_splot_tmp->reduce(cut.c_str());

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
  for (int i = 0; i < 10; i++) {
    // check weights of some events
    std::cout << "nSig weight: " << i << std::right << std::setw(12) << sData.GetSWeight(i, "nSig") << " /  nBkg weight: "
              << std::setw(12) << sData.GetSWeight(i, "nBkg") << " / Total Weight" << std::setw(12) << sData.GetSumOfEventSWeight(i)
              << "\n";
  }

  m_ws.import(*ds_splot);
}

void JpsiFitter::makeSplotPdfs(bool isForcedMax, double ctauErrMax_forced) {
  // ===== build RooPlot and histograms ===== //
  // a lot of legacy codes

  double ctauErrLow = 0, ctauErrHigh = 0.25;
  m_nBins = 100;

  std::cout << "SPlot: make SPlot PDFs\n";
  err_frame1 = m_ws.var("ctau3DErr")->frame(RooFit::Range(ctauErrLow, ctauErrHigh));

  m_hTot = (TH1D *)m_ws.data("ds_splot")->createHistogram(("hTot"), *m_ws.var("ctau3DErr"), RooFit::Binning(err_frame1->GetNbinsX(), err_frame1->GetXaxis()->GetXmin(), err_frame1->GetXaxis()->GetXmax()));
  TH1D *hTot_M = (TH1D *)m_ws.data("ds_splot")->createHistogram(("hTot_M"), *m_ws.var("ctau3DErr"), RooFit::Binning(m_nBins, ctauErrLow, ctauErrHigh));

  // ===== build temporal sPlot weighted datasets ===== //
  std::cout << "SPlot: build temporal weighted datasets\n";
  RooDataSet *dataw_Sig_b = new RooDataSet("dataw_Sig_b", "TMP_SIG_DATA", (RooDataSet *)m_ws.data("ds_splot"),
                                           RooArgSet(*m_ws.var("ctau3DErr"), *m_ws.var("nSig_sw"), *m_ws.var("ctau3DRes"), *m_ws.var("ctau3D"), *m_ws.var("mass")), 0, "nSig_sw");
  TH1D *hSig = (TH1D *)dataw_Sig_b->createHistogram(("hSig"), *m_ws.var("ctau3DErr"), RooFit::Binning(m_nBins, ctauErrLow, ctauErrHigh));

  // Bkg
  RooDataSet *dataw_Bkg_b = new RooDataSet("dataw_Bkg_b", "TMP_BKG_DATA", (RooDataSet *)m_ws.data("ds_splot"),
                                           RooArgSet(*m_ws.var("ctau3DErr"), *m_ws.var("nBkg_sw"), *m_ws.var("ctau3DRes"), *m_ws.var("ctau3D"), *m_ws.var("mass")), 0, "nBkg_sw");
  TH1D *hBkg = (TH1D *)dataw_Bkg_b->createHistogram(("hBkg"), *m_ws.var("ctau3DErr"), RooFit::Binning(m_nBins, ctauErrLow, ctauErrHigh));

  // ===== set new ranges of ctau3DErr ===== //
  double ctauErrMin_calc = 0;
  double ctauErrMax_calc = 10;
  // new ErrMin
  for (int i = 0; i < hTot_M->GetNbinsX() / 2; i++) {
    if (hTot_M->GetBinContent(i) > 1) {
      ctauErrMin_calc = hTot_M->GetBinLowEdge(i + 1);
      break;
    }
  }
  // new ErrMax
  for (int i = 0; i < hTot_M->GetNbinsX(); i++) {
    if (hTot_M->GetBinContent(i) >= 1 && hTot_M->GetBinContent(i + 1) < 1 && hTot_M->GetBinContent(i + 2) < 1)
    {
      // ctauErrMax_calc = hSig->GetBinLowEdge(i)+hSig->GetBinWidth(i);
      ctauErrMax_calc = hTot_M->GetBinLowEdge(i) + hTot_M->GetBinWidth(i);
      break;
    } else {
      ctauErrMax_calc = hTot_M->GetBinLowEdge(i);
    }
  }

  m_ctauErrMin = ctauErrMin_calc;
  m_ctauErrMax = ctauErrMax_calc;

  if (isForcedMax && m_ctauErrMin < ctauErrMax_forced)
    m_ctauErrMax = ctauErrMax_forced;
  else if(isForcedMax && m_ctauErrMin > ctauErrMax_forced) {
    std::cout << "Error: ctauErrMax_forced < m_ctauErrMin. Ignore forced_ctauErrMax\n";
  }

  std::cout << "\nSet new bin width\n";
  double BinWidth = (ctauErrHigh - ctauErrLow) / m_nBins;
  std::cout << "BinWidth : " << BinWidth << "\n";
  m_newBins = (m_ctauErrMax - m_ctauErrMin) / BinWidth;
  std::cout << "m_newBins : " << m_newBins << "\n\n";

  // ===== build sPlot weighted datasets and RooHistPdfs ===== //
  // results of SPlot process
  TH1D *h_tot_test = (TH1D *)m_ws.data("ds_splot")->createHistogram(("h_tot_test"), *m_ws.var("ctau3DErr"), RooFit::Binning(m_newBins, m_ctauErrMin, m_ctauErrMax));
  // TH1D* h_tot_test = (TH1D*)m_ws.data("ds_splot")->createHistogram(("h_tot_test"), *m_ws.var("ctau3DErr"),RooFit::Binning(m_nBins,m_ctauErrMin,m_ctauErrMax));
  RooDataHist *totHist = new RooDataHist("ds_splot", "", *m_ws.var("ctau3DErr"), h_tot_test);
  RooHistPdf *errTotModel = new RooHistPdf("errTotModel", "Total RooHistPdf of ctauErr", *m_ws.var("ctau3DErr"), *totHist);

  // bkg ds and pdf
  TH1D *hBkg_w = (TH1D *)dataw_Bkg_b->createHistogram(("hBkg_w"), *m_ws.var("ctau3DErr"), RooFit::Binning(m_newBins, m_ctauErrMin, m_ctauErrMax));
  RooDataHist *bkgHist = new RooDataHist("ds_splot", "", *m_ws.var("ctau3DErr"), hBkg_w);
  RooHistPdf *errBkgModel = new RooHistPdf("errBkgModel", "hist pdf", *m_ws.var("ctau3DErr"), *bkgHist);

  // sig ds and pdf
  TH1D *hSig_w = (TH1D *)dataw_Sig_b->createHistogram(("hSig_w"), *m_ws.var("ctau3DErr"), RooFit::Binning(m_newBins, m_ctauErrMin, m_ctauErrMax));
  RooDataHist *sigHist = new RooDataHist("ds_splot", "", *m_ws.var("ctau3DErr"), hSig_w);
  RooHistPdf *errSigModel = new RooHistPdf("errSigModel", "hist pdf", *m_ws.var("ctau3DErr"), *sigHist);

  // Re-build datasets
  RooDataSet *ds_errBkg = new RooDataSet("ds_errBkg", "", (RooDataSet *)m_ws.data("ds_splot"),
                                          RooArgSet(*m_ws.var("ctau3DErr"), *m_ws.var("nBkg_sw"), *m_ws.var("ctau3DRes"), *m_ws.var("ctau3D"), *m_ws.var("mass")), 0, "nBkg_sw");
  RooDataSet *ds_errSig = new RooDataSet("ds_errSig", "", (RooDataSet *)m_ws.data("ds_splot"),
                                          RooArgSet(*m_ws.var("ctau3DErr"), *m_ws.var("nSig_sw"), *m_ws.var("ctau3DRes"), *m_ws.var("ctau3D"), *m_ws.var("mass")), 0, "nSig_sw");

  // ===== import sPlot results ===== //
  m_ws.import(*ds_errSig);
  m_ws.import(*ds_errBkg);
  m_ws.import(*errTotModel);
  m_ws.import(*errSigModel);
  m_ws.import(*errBkgModel);


  m_minRange = (double)(floor(m_ctauErrMin * 100.) / 100.);
  m_maxRange = (double)(ceil(m_ctauErrMax * 100.) / 100.);
  m_ws.var("ctau3DErr")->setRange("ctauErrWindow", m_ctauErrMin, m_ctauErrMax);
  m_nBins = (m_maxRange - m_minRange) / BinWidth;
  err_frame1 = m_ws.var("ctau3DErr")->frame(RooFit::Bins(m_nBins), RooFit::Range(m_minRange - 0.01, m_maxRange + 0.01)); // modified
}

void JpsiFitter::drawSplot(const std::string &outPlotPath)
{
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

  // ===== print lost events due to new ctau3DErr range ===== //
  std::cout << "final ctau3DErr range\n";
  std::cout << m_ws.var("ctau3DErr")->getMin() << ", " << m_ws.var("ctau3DErr")->getMax() << "\n\n";

  RooDataSet *ctauResCutDS = (RooDataSet *)m_ws.data("ds_splot")->reduce(RooArgSet(*(m_ws.var("ctau3DRes")), *(m_ws.var("ctau3D")), *(m_ws.var("ctau3DErr"))), Form("ctau3DErr>%f&&ctau3DErr<%f", m_ctauErrMin, m_ctauErrMax));
  ctauResCutDS->SetName("ctauResCutDS"); // used only for comparison
  m_ws.import(*ctauResCutDS);

  Double_t outTot = m_ws.data("ds_splot")->numEntries();
  Double_t outRes = m_ws.data("ctauResCutDS")->numEntries(); // use it when you apply forced ctauErrMax
  std::cout << "Check how many events we lost due to new range cut\n";
  std::cout << "Total evt: " << outTot << "\n";
  std::cout << "Residual evt: " << outRes << "\n";
  std::cout << "lost evt: " << ((outTot - outRes) * 100) / outTot << " %\n\n";

  // ===== continue plotting =====
  // draw dist. on pad1
  // keep the legacy Name("legacy") -> Need to track in legend and pull.
  RooPlot *errFrame2 = (RooPlot *)err_frame1->Clone("errFrame2");

  // think about conditions wrt ctauResCutDS and ds_splot.
  m_ws.data("ctauResCutDS")->plotOn(errFrame2, RooFit::Name("dataCTAUERR_Tot"), RooFit::MarkerSize(.7), RooFit::Binning(m_newBins)); // Normalization(m_ws.->ds_splot("reducedDS_MC")->sumEntries()
  m_ws.pdf("errTotModel")->plotOn(errFrame2, RooFit::Name("errTotModel"), RooFit::LineColor(kGreen + 1), RooFit::Range(m_ctauErrMin, m_ctauErrMax), RooFit::LineWidth(2), RooFit::Normalization(m_ws.data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent));

  m_ws.data("ds_errSig")->plotOn(errFrame2, RooFit::Name("dataHist_Sig"), RooFit::MarkerSize(.7), RooFit::LineColor(kRed + 2), RooFit::MarkerColor(kRed + 2), RooFit::Binning(m_newBins));
  m_ws.pdf("errSigModel")->plotOn(errFrame2, RooFit::Name("errSigModel"), RooFit::LineColor(kRed + 2), RooFit::LineWidth(2), RooFit::Range(m_ctauErrMin, m_ctauErrMax));
  m_ws.data("ds_errBkg")->plotOn(errFrame2, RooFit::Name("dataHist_Bkg"), RooFit::MarkerSize(.7), RooFit::LineColor(kBlue + 2), RooFit::MarkerColor(kBlue + 2), RooFit::Binning(m_newBins));
  m_ws.pdf("errBkgModel")->plotOn(errFrame2, RooFit::Name("errBkgModel"), RooFit::LineColor(kBlue + 2), RooFit::LineWidth(2), RooFit::Range(m_ctauErrMin, m_ctauErrMax));

  // set range of y axis
  Double_t YMax = m_hTot->GetBinContent(m_hTot->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= m_hTot->GetNbinsX(); i++)
    if (m_hTot->GetBinContent(i) > 0)
      YMin = TMath::Min(YMin, m_hTot->GetBinContent(i));
  Double_t Yup(0.), Ydown(0.);
  Yup = YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.4 - 0.3)));
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.3 / (1.0 - 0.4 - 0.3))));
  errFrame2->GetYaxis()->SetRangeUser(Ydown, Yup);



  // ===== draw vertical lines to show ctau3DErr range ===== //
  TLine *minline = new TLine(m_ctauErrMin, 0.0, m_ctauErrMin, (Ydown * TMath::Power((Yup / Ydown), 0.5)));
  minline->SetLineStyle(2);
  minline->SetLineColor(1);
  minline->SetLineWidth(3);
  errFrame2->addObject(minline);
  TLine *maxline = new TLine(m_ctauErrMax, 0.0, m_ctauErrMax, (Ydown * TMath::Power((Yup / Ydown), 0.5)));
  maxline->SetLineStyle(2);
  maxline->SetLineColor(1);
  maxline->SetLineWidth(3);
  errFrame2->addObject(maxline);

  errFrame2->GetXaxis()->CenterTitle();
  errFrame2->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} Error (mm)");
  errFrame2->SetFillStyle(4000);
  errFrame2->GetYaxis()->SetTitleOffset(1.43);
  errFrame2->GetXaxis()->SetLabelSize(0);
  errFrame2->GetXaxis()->SetTitleSize(0);
  errFrame2->Draw();

  // ===== draw legends ===== //
  // TLegend *leg_err = new TLegend(text_x + 0.5, text_y - 0.2, text_x + 0.7, text_y);
  // leg_err->SetTextSize(text_size);
  // leg_err->SetTextFont(43);
  // leg_err->SetBorderSize(0);
  // leg_err->AddEntry(errFrame2->findObject("dataCTAUERR_Tot"), "Data", "pe");
  // leg_err->AddEntry(errFrame2->findObject("errTotModel"), "Total PDF", "l");
  // leg_err->AddEntry(errFrame2->findObject("errSigModel"), "Signal", "l");
  // leg_err->AddEntry(errFrame2->findObject("errBkgModel"), "Background", "l");
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
  TPad *pullPad = new TPad("pullPad", "", 0, 0.006, 0.98, 0.227);
  c_err->cd();
  pullPad->Draw();
  pullPad->cd();
  pullPad->SetTopMargin(0); // Upper and lower plot are joined
  pullPad->SetBottomMargin(0.67);
  pullPad->SetBottomMargin(0.4);
  pullPad->SetFillStyle(4000);
  pullPad->SetFrameFillStyle(4000);
  pullPad->SetTicks(1, 1);

  RooPlot *pullTmp = (RooPlot *)errFrame2->Clone("pullTmp");
  RooHist *pullHist = pullTmp->pullHist("dataCTAUERR_Tot", "errTotModel");
  pullHist->SetMarkerSize(0.8);
  RooPlot *pullFrame = m_ws.var("ctau3DErr")->frame(RooFit::Bins(m_nBins), RooFit::Range(m_minRange - 0.01, m_maxRange + 0.01));
  pullFrame->addPlotable(pullHist, "PX");
  pullFrame->SetTitle("");
  pullFrame->SetTitleSize(0);
  pullFrame->GetYaxis()->SetTitleOffset(0.3);
  pullFrame->GetYaxis()->SetTitle("Pull");
  pullFrame->GetYaxis()->SetTitleSize(0.15);
  pullFrame->GetYaxis()->SetLabelSize(0.15);
  pullFrame->GetYaxis()->SetRangeUser(-3.8, 3.8);
  pullFrame->GetYaxis()->CenterTitle();

  pullFrame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} Error (mm)");
  pullFrame->GetXaxis()->SetTitleOffset(1.05);
  pullFrame->GetXaxis()->SetLabelOffset(0.04);
  pullFrame->GetXaxis()->SetLabelSize(0.15);
  pullFrame->GetXaxis()->SetTitleSize(0.15);
  pullFrame->GetXaxis()->CenterTitle();

  pullFrame->GetYaxis()->SetTickSize(0.04);
  pullFrame->GetYaxis()->SetNdivisions(404);
  pullFrame->GetXaxis()->SetTickSize(0.03);
  pullFrame->Draw();

  // draw horizontal line at center -> It's not a result of fit so pull points should be 0
  TLine *l_pull = new TLine(m_minRange - 0.01, 0, m_maxRange + 0.01, 0);
  l_pull->SetLineStyle(1);
  l_pull->Draw("same");
  pullPad->Update();

  // ===== draw plots ===== //
  c_err->Update();
  c_err->SaveAs(outPlotPath.c_str());

  TFile f_out_test("roots/splot.root", "recreate");
  f_out_test.cd();
  m_ws.Write();
  f_out_test.Close();
}

bool JpsiFitter::loadSPlotResult(const std::string &filePath) {
  auto f_splot = std::unique_ptr<TFile>(TFile::Open(filePath.c_str()));
  if (!f_splot || f_splot->IsZombie())
  {
    std::cerr << "ERROR: Can't open the sPlot result file: " << filePath << "\n";
    return false;
  }
  std::cout << "INFO: sPlot result file '" << filePath << "' loaded.\n";

  RooWorkspace *wsMy = (RooWorkspace *)f_splot->Get("wsMy");
  RooDataSet *ds_errSig = (RooDataSet *)wsMy->data("ds_errSig");
  if (!ds_errSig)
  {
    std::cerr << "ERROR: Can't find 'ds_errSig' in the sPlot result file.\n";
    return false;
  }
  m_ws.import(*ds_errSig);
  return true;
}


void JpsiFitter::performResFit() {
  std::cout << "INFO: Performing Resolution Fit...\n";

  // setup variables
  RooAbsPdf *pdf = m_ws.pdf("ctauResModel");
  RooDataSet *data = (RooDataSet *)m_ws.data("ds_errSig");
  RooRealVar *obs = m_ws.var("ctau3DRes");
  if (!pdf || !data || !obs)
  {
    std::cerr << "ERROR: Missing components for resolution fit.\n";
    return;
  }

  // set fit range
  TH1 *h_tmp = data->createHistogram("h_tmp_res", *obs);
  double ctauResMin = -10; // default
  for (int i = 1; i <= h_tmp->GetNbinsX() / 2; i++)
  {
    if (h_tmp->GetBinContent(i) <= 0 && h_tmp->GetBinContent(i + 1) <= 0)
    {
      ctauResMin = h_tmp->GetBinLowEdge(i + 2);
    }
  }
  delete h_tmp; 

  obs->setRange("ctauResWindow", ctauResMin, 0); // nominal. fit only half
  std::cout << "INFO: Fit Range for ctau3DRes set to: [" << ctauResMin << ", 0]\n";

  // make a new dataset with new ctauRes range
  RooDataSet *dataToFit = (RooDataSet *)data->reduce(RooFit::CutRange("ctauResWindow"))->Clone("dataToFit");

  RooFitResult *fit_res = pdf->fitTo(*dataToFit, RooFit::Save(true), RooFit::SumW2Error(data->isWeighted()), RooFit::Extended(kTRUE), RooFit::NumCPU(4));

  std::cout << "INFO: Resolution fit result:\n";
  fit_res->Print("V");

  // fix parameters
  m_ws.var("s1_CtauRes")->setConstant(kTRUE);
  m_ws.var("rS21_CtauRes")->setConstant(kTRUE);
  m_ws.var("f_CtauRes")->setConstant(kTRUE);

  m_ws.import(*dataToFit);
}

void JpsiFitter::drawResPlot(const std::string &outputPath) {
  using namespace RooFit;
  std::cout << "INFO: Drawing Resolution Fit results...\n";

  RooRealVar *obs = m_ws.var("ctau3DRes");
  RooDataSet *data = (RooDataSet *)m_ws.data("dataToFit");
  RooAbsPdf *pdf = m_ws.pdf("ctauResModel");
  if (!obs || !data || !pdf)
  {
    std::cerr << "ERROR: Missing components for drawing resolution fit.\n";
    return;
  }

  TCanvas *c_res = new TCanvas("c_res", "Resolution Fit Plot", 600, 600);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
  pad1->SetBottomMargin(0.02);
  pad1->SetLogy();
  pad1->Draw();

  c_res->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.3);
  pad2->Draw();

  // res pad
  pad1->cd();
  RooPlot *frame = obs->frame(Title("Resolution Fit"));

  // draw dataset and pdfs
  data->plotOn(frame, Name("data"));
  pdf->plotOn(frame, Name("total_pdf"), LineColor(kBlack));
  pdf->plotOn(frame, Name("gauss1"), Components("GaussModel1_ctauRes"), LineColor(kGreen + 2), LineStyle(kDashed));
  pdf->plotOn(frame, Name("gauss2"), Components("GaussModel2_ctauRes"), LineColor(kRed + 2), LineStyle(kDashed));

  // set Y-axis
  TH1 *h_tmp = data->createHistogram("h_tmp_draw", *obs);
  double yMax = h_tmp->GetBinContent(h_tmp->GetMaximumBin());
  frame->GetYaxis()->SetRangeUser(0.1, yMax * 1.5);
  delete h_tmp;

  frame->GetXaxis()->SetLabelSize(0);
  frame->Draw();


  // draw pull
  pad2->cd();
  RooHist *pullHist = frame->pullHist("data", "total_pdf");
  RooPlot *pullFrame = obs->frame(Title("Pull Distribution"));
  pullFrame->addPlotable(pullHist, "P");
  pullFrame->GetYaxis()->SetRangeUser(-5, 5);
  pullFrame->GetYaxis()->SetTitle("Pull");
  pullFrame->GetYaxis()->CenterTitle();
  pullFrame->Draw();

  TLine *zeroLine = new TLine(obs->getMin(), 0, obs->getMax(), 0);
  zeroLine->SetLineColor(kRed);
  zeroLine->SetLineStyle(kDashed);
  zeroLine->Draw("same");

  c_res->SaveAs(outputPath.c_str());
}