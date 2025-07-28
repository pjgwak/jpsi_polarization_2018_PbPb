#include "CtauErrFit.h"
#include <iostream>

#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"

#include "RooWorkspace.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooStats/SPlot.h"
#include "../headers/polarizationUtilities.h"

using namespace RooFit;

CtauErrFit::CtauErrFit(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
                       float cosLow, float cosHigh, int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP)
    : ptLow(ptLow), ptHigh(ptHigh), yLow(yLow), yHigh(yHigh), cLow(cLow), cHigh(cHigh),
      cosLow(cosLow), cosHigh(cosHigh), PR(PR), PRw(PRw), fEffW(fEffW), fAccW(fAccW), isPtW(isPtW), isTnP(isTnP)
{
  gStyle->SetEndErrorSize(0);

  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
}

CtauErrFit::~CtauErrFit()
{
  RooMsgService::instance().getStream(0).addTopic(Caching);
  RooMsgService::instance().getStream(1).addTopic(Caching);
  RooMsgService::instance().getStream(0).addTopic(Plotting);
  RooMsgService::instance().getStream(1).addTopic(Plotting);
  RooMsgService::instance().getStream(0).addTopic(Integration);
  RooMsgService::instance().getStream(1).addTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
  
  delete fInputData;
  delete fMass;
  delete ws;
}

void CtauErrFit::run()
{
  setLabels();
  openInputFiles();
  setupWorkspaceAndData();
  performSPlot();
  setVariableRanges();
  buildPdfAndData();
  makePlot();
  saveResults();
}

void CtauErrFit::setLabels()
{
  std::cout << "===== setLabels() =====\n\n";
  DATE = "No_Weight_2";
  gStyle->SetEndErrorSize(0);
  gSystem->mkdir(Form("roots/2DFit_%s/", DATE.c_str()), kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/", DATE.c_str()), kTRUE);
  gSystem->mkdir(Form("../figs/2DFit_%s/CtauErr", DATE.c_str()), kTRUE);

  if (PRw == 1) fname = "PR";
  else if (PRw == 2) fname = "NP";

  kineLabel = std::string(getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh).Data());
}

void CtauErrFit::openInputFiles()
{
  std::cout << "===== openInputFiles() =====\n\n";
  fInputData = new TFile(inputFilePath.c_str());
  fMass = new TFile(Form("roots/2DFit_%s/MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));

  if (!fInputData || fInputData->IsZombie())
  {
    std::cerr << "CRITICAL: Input data file could not be opened. Aborting.\n";
    exit(1);
  }
  if (!fMass || fMass->IsZombie())
  {
    std::cerr << "CRITICAL: Mass fit result file could not be opened. Aborting.\n";
    exit(1);
  }
}

void CtauErrFit::setupWorkspaceAndData()
{
  std::cout << "===== setupWorkspaceAndData() =====\n\n";
  ws = new RooWorkspace("workspace");
  RooDataSet *dataset = (RooDataSet *)fInputData->Get("dataset");
  RooDataSet *datasetMass = (RooDataSet *)fMass->Get("datasetMass");
  RooAddPdf *pdfMASS_Tot = (RooAddPdf *)fMass->Get("pdfMASS_Tot");

  ws->import(*dataset);
  ws->import(*datasetMass);
  ws->import(*pdfMASS_Tot);
  // TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f && cBin>=%d && cBin<%d", ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh); // PbPb
  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f", ptLow, ptHigh, yLow, yHigh, massLow, massHigh); // OO
  TString accCut = "(((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2) && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)))"; // 2018 acceptance cut
  TString OSCut = "recoQQsign==0";
  // TString angleCut = Form("cos_ep>%.2f && cos_ep<%.2f", cosLow, cosHigh);
  TString angleCut = "true";
  TString finalCut = TString::Format("%s && %s && %s && %s", OSCut.Data(), accCut.Data(), kineCut.Data(), angleCut.Data());

  RooArgSet *argSet = new RooArgSet(*(ws->var("ctau3D")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y")), *(ws->var("weight")), *(ws->var("ctau3DRes")), *(ws->var("ctau3DErr")));
  argSet->add(*(ws->var("pt1")));
  argSet->add(*(ws->var("pt2")));
  argSet->add(*(ws->var("eta1")));
  argSet->add(*(ws->var("eta2")));
  argSet->add(*(ws->var("recoQQsign")));
  // argSet->add(*(ws->var("cBin"))); // Todo: turn on it for PbPb - OO

  // Don't use MC weight when use the sPlot technique (?) -> check legacy AN code too.
  RooDataSet *dsAB = (RooDataSet *)ws->data("dataset")->reduce(*argSet, finalCut.Data());
  ws->import(*dsAB, Rename("dsAB"));

  delete argSet;
}



void CtauErrFit::performSPlot()
{
  std::cout << "===== performSPlot() =====\n\n";
  RooRealVar *sigYield = ws->var("N_Jpsi");
  RooRealVar *bkgYield = ws->var("N_Bkg");

  sigYield->setMin(0); // sPlot requires a range: 0 ~ max
  bkgYield->setMin(0);

  RooArgList yieldList(*sigYield, *bkgYield);

  RooDataSet *data = (RooDataSet *)ws->data("dsAB")->Clone("TMP_DATA");
  
  // legacy codes.
  RooArgSet *cloneSet = (RooArgSet *)RooArgSet(*ws->pdf("pdfMASS_Tot"), "pdfMASS_Tot").snapshot(kTRUE);
  auto clone_mass_pdf = (RooAbsPdf *)cloneSet->find("pdfMASS_Tot");
  clone_mass_pdf->setOperMode(RooAbsArg::ADirty, kTRUE);

  RooStats::SPlot sData = RooStats::SPlot("sData", "An SPlot", *data, clone_mass_pdf, yieldList);

  ws->import(*data, Rename("dataset_SPLOT"));
  std::cout << "[INFO] Jpsi yield -> Mass Fit:" << ws->var("N_Jpsi")->getVal() << ", sWeights :" << sData.GetYieldFromSWeight("N_Jpsi") << "\n";
  std::cout << "[INFO] Bkg  yield -> Mass Fit:" << ws->var("N_Bkg")->getVal() << ", sWeights :" << sData.GetYieldFromSWeight("N_Bkg") << "\n";

  // make sWeighted dataset
  RooDataSet *ds_sigAllRange = new RooDataSet("ds_sigAllRange", "", (RooDataSet *)ws->data("dataset_SPLOT"),
                                           RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Jpsi_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("mass")), 0, "N_Jpsi_sw");
  

  RooDataSet *ds_bkgAllRange = new RooDataSet("ds_bkgAllRange", "TMP_BKG_DATA", (RooDataSet *)ws->data("dataset_SPLOT"),
                                           RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Bkg_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("mass")), 0, "N_Bkg_sw");

  ws->import(*ds_sigAllRange);
  ws->import(*ds_bkgAllRange);

  delete data;
}

void CtauErrFit::setVariableRanges()
{
  std::cout << "===== setVariableRanges() =====\n\n";
  ws->var("ctau3DErr")->setRange(ctauErrLow, ctauHigh);

  int nBinInit = std::min(int(round((ws->var("ctau3DErr")->getMax() - ws->var("ctau3DErr")->getMin()) / 0.0025)), 100);

  // TH1D *hSig = (TH1D *)ds_sigAllRange->createHistogram(("hSig"), *ws->var("ctau3DErr"), Binning(nBinInit, ctauErrLow, ctauErrHigh));
  TH1D *hTmp = (TH1D *)ws->data("dsAB")->createHistogram(("hTmp"), *ws->var("ctau3DErr"), Binning(nBinInit, ctauErrLow, ctauErrHigh));
  for (int i = 1; i < hTmp->GetNbinsX() / 2; i++)
  {
    if (hTmp->GetBinContent(i) > 1)
    { 
      ctauErrMin = hTmp->GetBinLowEdge(i + 1);
      break;
    }
  }

  if (ctauErrMax == -1) { // if user doesn't set the max, find a max automatically
    for (int i = 1; i < hTmp->GetNbinsX(); i++)
    {
      if (hTmp->GetBinContent(i) >= 1 && hTmp->GetBinContent(i + 1) < 1 && hTmp->GetBinContent(i + 2) < 1)
      {
        ctauErrMax = hTmp->GetBinLowEdge(i) + hTmp->GetBinWidth(i);
        break;
      }
      else
        ctauErrMax = hTmp->GetBinLowEdge(i);
    }
  }
  
  double binWidth = (ctauErrHigh - ctauErrLow) / nBinInit;
  nBins = (ctauErrMax - ctauErrMin) / binWidth;

  ws->var("ctau3DErr")->setRange(ctauErrMin, ctauErrMax);
  ws->var("ctau3DErr")->setRange("sPlotRange", ctauErrMin, ctauErrMax);
  std::cout << "ctau3DErr range: [" << ctauErrMin << ", " << ctauErrMax << "]\n";

  delete hTmp;
}

void CtauErrFit::buildPdfAndData()
{
  std::cout << "===== buildModels() =====\n\n";  
  // total
  TH1D *hTot = (TH1D *)ws->data("dsAB")->createHistogram(("hTot"), *ws->var("ctau3DErr"), Binning(nBins, ctauErrMin, ctauErrMax));
  RooDataHist *totHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), hTot);
  RooHistPdf *pdfCTAUERR_Tot = new RooHistPdf("pdfCTAUERR_Tot", "hist pdf", *ws->var("ctau3DErr"), *totHist);

  // bkg
  TH1D *hBkg_w = (TH1D *)ws->data("ds_bkgAllRange")->createHistogram(("hBkg_w"), *ws->var("ctau3DErr"), Binning(nBins, ctauErrMin, ctauErrMax));
  RooDataHist *bkgHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), hBkg_w);
  RooHistPdf *pdfCTAUERR_Bkg = new RooHistPdf("pdfCTAUERR_Bkg", "hist pdf", *ws->var("ctau3DErr"), *bkgHist);

  // sig
  TH1D *hSig_w = (TH1D *)ws->data("ds_sigAllRange")->createHistogram(("hSig_w"), *ws->var("ctau3DErr"), Binning(nBins, ctauErrMin, ctauErrMax));
  RooDataHist *sigHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), hSig_w);
  RooHistPdf *pdfCTAUERR_Jpsi = new RooHistPdf("pdfCTAUERR_Jpsi", "hist pdf", *ws->var("ctau3DErr"), *sigHist);

  // dataset for next steps
  RooDataSet *dataw_Bkg = new RooDataSet("dataw_Bkg", "TMP_BKG_DATA", (RooDataSet *)ws->data("dataset_SPLOT"),
                                         RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Bkg_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("mass")), 0, "N_Bkg_sw");
  RooDataSet *dataw_Sig = new RooDataSet("dataw_Sig", "TMP_SIG_DATA", (RooDataSet *)ws->data("dataset_SPLOT"),
                                         RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Jpsi_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("mass")), 0, "N_Jpsi_sw");

  ws->import(*dataw_Sig);
  ws->import(*dataw_Bkg);
  ws->import(*pdfCTAUERR_Tot);
  ws->import(*pdfCTAUERR_Jpsi);
  ws->import(*pdfCTAUERR_Bkg);

  delete hTot;
  delete hBkg_w;
  delete hSig_w;

  delete dataw_Bkg;
  delete dataw_Sig;
}


void CtauErrFit::makePlot()
{
  std::cout << "===== makePlot() =====\n\n";
  // double minRange = (double)(floor(ctauErrMin * 100.) / 100.);
  // double maxRange = (double)(ceil(ctauErrMax * 100.) / 100.);
  double minRange = ctauErrMin;
  double maxRange = ctauErrMax;

  ws->var("ctau3DErr")->setRange("sPlotRange", ctauErrMin, ctauErrMax);

  TCanvas *myCanvas = new TCanvas("myCanvas", "My plots", 554, 4, 550, 520);
  myCanvas->cd();

  TPad *pad1 = new TPad("pad1", "", 0, 0.16, 0.98, 1.0);
  pad1->SetTicks(1, 1);
  pad1->Draw();
  pad1->cd();
  gPad->SetLogy();

  // RooPlot *ctauFrame = ws->var("ctau3DErr")->frame(Range(ctauErrLow, ctauErrHigh));
  RooPlot *ctauFrame = ws->var("ctau3DErr")->frame(Bins(nBins), Range(minRange - 0.01, maxRange + 0.01)); // modified
  ctauFrame->SetTitle("");


  ws->data("dsAB")->plotOn(ctauFrame, Name("dataCTAUERR_Tot"), MarkerSize(.7), Binning(nBins)); // Normalization(wsmc->data("reducedDS_MC")->sumEntries()
  ws->pdf("pdfCTAUERR_Tot")->plotOn(ctauFrame, Name("pdfCTAUERR_Tot"), LineColor(kGreen + 1), Range(ctauErrMin, ctauErrMax), LineWidth(2), Normalization(ws->data("dsAB")->sumEntries(), RooAbsReal::NumEvent));
  ws->data("dataw_Sig")->plotOn(ctauFrame, Name("dataHist_Sig"), MarkerSize(.7), LineColor(kRed + 2), MarkerColor(kRed + 2), Binning(nBins));
  ws->pdf("pdfCTAUERR_Jpsi")->plotOn(ctauFrame, Name("pdfCTAUERR_Jpsi"), LineColor(kRed + 2), LineWidth(2), Range(ctauErrMin, ctauErrMax));
  ws->data("dataw_Bkg")->plotOn(ctauFrame, Name("dataHist_Bkg"), MarkerSize(.7), LineColor(kBlue + 2), MarkerColor(kBlue + 2), Binning(nBins));
  ws->pdf("pdfCTAUERR_Bkg")->plotOn(ctauFrame, Name("pdfCTAUERR_Bkg"), LineColor(kBlue + 2), LineWidth(2), Range(ctauErrMin, ctauErrMax));

  TH1D *hTmp = (TH1D *)ws->data("dsAB")->createHistogram(("hTmp"), *ws->var("ctau3DErr"), Binning(nBins, ctauErrMin, ctauErrMax));
  Double_t YMax = hTmp->GetBinContent(hTmp->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= hTmp->GetNbinsX(); i++)
    if (hTmp->GetBinContent(i) > 0)
      YMin = std::min(YMin, hTmp->GetBinContent(i));
  Double_t Yup(0.), Ydown(0.);
  Yup = YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.4 - 0.3)));
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.3 / (1.0 - 0.4 - 0.3))));
  delete hTmp;
  ctauFrame->GetYaxis()->SetRangeUser(Ydown, Yup);

  // cout<<ctauErrLow<<", "<<ctauErrHigh<<endl;
  std::cout << ws->var("ctau3DErr")->getMin() << ", " << ws->var("ctau3DErr")->getMax() << "\n";
  RooDataSet *ctauResCutDS = (RooDataSet *)ws->data("dataw_Sig")->reduce(RooArgSet(*(ws->var("ctau3DRes")), *(ws->var("ctau3D")), *(ws->var("ctau3DErr"))), Form("ctau3DErr>%f&&ctau3DErr<%f", ctauErrMin, ctauErrMax));
  ctauResCutDS->SetName("ctauResCutDS");
  ws->import(*ctauResCutDS);
  Double_t outTot = ws->data("dsAB")->numEntries();
  Double_t outRes = ws->data("ctauResCutDS")->numEntries();
  std::cout << "Tot evt: " << outTot << "" << "\n";
  std::cout << "Res evt: " << outRes << "" << "\n";
  std::cout << "lost evt: " << outTot - outRes << " evt (" << ((outTot - outRes) * 100) / outTot << " %)" << "\n";

  TLine *minline = new TLine(ctauErrMin, 0.0, ctauErrMin, (Ydown * TMath::Power((Yup / Ydown), 0.5)));
  minline->SetLineStyle(2);
  minline->SetLineColor(1);
  minline->SetLineWidth(3);
  ctauFrame->addObject(minline);
  TLine *maxline = new TLine(ctauErrMax, 0.0, ctauErrMax, (Ydown * TMath::Power((Yup / Ydown), 0.5)));
  maxline->SetLineStyle(2);
  maxline->SetLineColor(1);
  maxline->SetLineWidth(3);
  ctauFrame->addObject(maxline);

  ctauFrame->GetXaxis()->CenterTitle();
  ctauFrame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} Error (mm)");
  ctauFrame->SetFillStyle(4000);
  ctauFrame->GetYaxis()->SetTitleOffset(1.43);
  ctauFrame->GetXaxis()->SetLabelSize(0);
  ctauFrame->GetXaxis()->SetTitleSize(0);
  ctauFrame->Draw();
  // Double_t outTot = ws->data("dsAB")->numEntries();
  // Double_t outErr = ws->data("dsAB")->reduce(Form("(ctauErr>=%.6f || ctauErr<=%.6f)", ctauErrHigh, ctauErrLow))->numEntries();
  // cout<<(outErr*100)/outTot<<endl;
  TLegend *leg_B = new TLegend(text_x + 0.5, text_y - 0.2, text_x + 0.7, text_y);
  leg_B->SetTextSize(text_size);
  leg_B->SetTextFont(43);
  leg_B->SetBorderSize(0);
  leg_B->AddEntry(ctauFrame->findObject("dataCTAUERR_Tot"), "Data", "pe");
  leg_B->AddEntry(ctauFrame->findObject("pdfCTAUERR_Tot"), "Total PDF", "l");
  leg_B->AddEntry(ctauFrame->findObject("pdfCTAUERR_Jpsi"), "Signal", "l");
  leg_B->AddEntry(ctauFrame->findObject("pdfCTAUERR_Bkg"), "Background", "l");
  leg_B->Draw("same");
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x, text_y, text_color, text_size);
  if (yLow == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x, text_y - y_diff, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x, text_y - y_diff, text_color, text_size);
  drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x, text_y - y_diff * 2, text_color, text_size);
  // drawText(Form("n_{J/#psi} = %.f #pm %.f",sData.GetYieldFromSWeight("N_Jpsi"),ws->var("N_Jpsi")->getError()),text_x,text_y-y_diff*3,text_color,text_size);
  // drawText(Form("n_{Bkg} = %.f #pm %.f",sData.GetYieldFromSWeight("N_Bkg"), ws->var("N_Bkg")->getError()),text_x,text_y-y_diff*4,text_color,text_size);
  drawText(Form("Loss: (%.4f%s) %.f evts", (outTot - outRes) * 100 / outTot, "%", outTot - outRes), text_x, text_y - y_diff * 3, text_color, text_size);
  TPad *pad_B_2 = new TPad("pad_B_2", "pad_B_2", 0, 0.006, 0.98, 0.227);
  myCanvas->cd();
  pad_B_2->Draw();
  pad_B_2->cd();
  pad_B_2->SetTopMargin(0); // Upper and lower plot are joined
  pad_B_2->SetBottomMargin(0.67);
  pad_B_2->SetBottomMargin(0.4);
  pad_B_2->SetFillStyle(4000);
  pad_B_2->SetFrameFillStyle(4000);
  pad_B_2->SetTicks(1, 1);

  RooPlot *frameTMP = (RooPlot *)ctauFrame->Clone("frameTMP");
  RooHist *pullPad = frameTMP->pullHist("dataCTAUERR_Tot", "pdfCTAUERR_Tot");
  pullPad->SetMarkerSize(0.8);
  RooPlot *pullFrame = ws->var("ctau3DErr")->frame(Bins(nBins), Range(minRange - 0.01, maxRange + 0.01));
  pullFrame->addPlotable(pullPad, "PX");
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

  TLine *lB = new TLine(minRange - 0.01, 0, maxRange + 0.01, 0);
  lB->SetLineStyle(1);
  lB->Draw("same");
  pad_B_2->Update();


  myCanvas->Update();
  myCanvas->SaveAs(Form("../figs/2DFit_%s/CtauErr/CtauErr_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));
  myCanvas->SaveAs(Form("figs/2DFit_%s/CtauErr_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.png", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));


}

void CtauErrFit::saveResults()
{
  std::cout << "===== saveResults() =====\n\n";
  TFile *outputFile = new TFile(Form("roots/2DFit_%s/CtauErrResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP), "recreate");

  ws->data("dataw_Sig")->Write();
  ws->data("dataw_Bkg")->Write();
  ws->pdf("pdfCTAUERR_Jpsi")->Write();
  ws->pdf("pdfCTAUERR_Bkg")->Write();
  ws->pdf("pdfCTAUERR_Tot")->Write();

  outputFile->Close();
  delete outputFile;
}