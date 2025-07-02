#include "CtauErrFit.h"
#include <iostream>
#include <algorithm> // min()
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TMath.h"
#include "TLegend.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooMsgService.h"
#include "RooStats/SPlot.h"
#include "../../../headers/polarizationUtilities.h"

using namespace RooFit;

CtauErrFit::CtauErrFit(float ptLow, float ptHigh,
                             float yLow, float yHigh,
                             int cLow, int cHigh,
                             float cosLow, float cosHigh,
                             int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP)
    : ptLow(ptLow), ptHigh(ptHigh),
      yLow(yLow), yHigh(yHigh),
      cLow(cLow), cHigh(cHigh),
      cosLow(cosLow), cosHigh(cosHigh),
      PR(PR), PRw(PRw),
      fEffW(fEffW), fAccW(fAccW), isPtW(isPtW), isTnP(isTnP)
{
  gStyle->SetEndErrorSize(0);
  kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling);
  // RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
}

// Destructor
CtauErrFit::~CtauErrFit()
{
  RooMsgService::instance().getStream(1).addTopic(RooFit::ObjectHandling);
  delete ws;
  delete sData;
}

void CtauErrFit::run()
{
  std::cout << "===== Start run() =====\n";
  setLabels();
  prepareDataset();
  performSPlot();
  buildCtauErrModel();
  drawResult();
  saveResult();
}
void CtauErrFit::setLabels()
{
  std::cout << "===== Start setLabels() =====\n";
  DATE = "No_Weight_2";
  gSystem->mkdir(Form("roots/2DFit_%s/CtauErr", DATE.Data()), kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/CtauErr", DATE.Data()), kTRUE);

  fname;
  if (PRw == 1)
    fname = "PR";
  else if (PRw == 2)
    fname = "NP";
}
void CtauErrFit::prepareDataset()
{
  std::cout << "===== Start prepareDataset() =====\n";
  TFile *f1 = new TFile("../../../files_roodata/RooDataSet_miniAOD_isMC0_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root");
  TFile *fMass = new TFile(Form("roots/2DFit_%s/Mass/Mass_FixedFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));

  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f && cBin>=%d && cBin<%d", ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&"; // 2018 acceptance cut
  TString OS = "recoQQsign==0 &&";
  TString angle_cut = Form("&& cos_ep>%.2f && cos_ep<%.2f", cosLow, cosHigh);
  TString final_cut = OS + accCut + kineCut + angle_cut;

  RooDataSet *dataset = (RooDataSet *)f1->Get("dataset");
  RooWorkspace *ws_mass = (RooWorkspace *)fMass->Get("ws_mass");
  RooAddPdf *pdfMASS_Tot = (RooAddPdf *)ws_mass->pdf("pdfMASS_Tot");

  ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  ws->import(*pdfMASS_Tot);

  RooArgSet *argSetWo = new RooArgSet(*(ws->var("ctau3D")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y")), *(ws->var("ctau3DRes")), *(ws->var("ctau3DErr")));
  argSetWo->add(*(ws->var("pt1")));
  argSetWo->add(*(ws->var("pt2")));
  argSetWo->add(*(ws->var("eta1")));
  argSetWo->add(*(ws->var("eta2")));
  argSetWo->add(*(ws->var("recoQQsign")));
  argSetWo->add(*(ws->var("cBin")));

  RooArgSet *argSetW = new RooArgSet(*(ws->var("ctau3D")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y")), *(ws->var("weight")), *(ws->var("ctau3DRes")), *(ws->var("ctau3DErr")));
  argSetW->add(*(ws->var("pt1")));
  argSetW->add(*(ws->var("pt2")));
  argSetW->add(*(ws->var("eta1")));
  argSetW->add(*(ws->var("eta2")));
  argSetW->add(*(ws->var("recoQQsign")));
  argSetW->add(*(ws->var("cBin")));

  RooDataSet *datasetW = new RooDataSet("datasetW", "A sample",
                                        *argSetW,
                                        Import(*dataset), WeightVar(*ws->var("weight")));
  RooDataSet *datasetWo = new RooDataSet("datasetWo", "A sample",
                                         *argSetWo,
                                         Import(*dataset)); //, WeightVar(*ws->var("weight")));
  ws->import(*datasetW);
  ws->import(*datasetWo);

  RooDataSet *dsAB = (RooDataSet *)datasetW->reduce(kineCut.Data());
  RooDataSet *dsAB1 = (RooDataSet *)datasetWo->reduce(kineCut.Data());
  dsAB->SetName("dsAB");
  dsAB1->SetName("dsAB1");
  ws->import(*dsAB);
  ws->import(*dsAB1);

  std::cout << "Weight : " << ws->var("weight")->getVal() << "\n";
  std::cout << "pt: " << ptLow << "-" << ptHigh << ", y: " << yLow << "-" << yHigh << ", Cent: " << cLow << "-" << cHigh << "%" << "\n";
  std::cout << "####################################" << "\n";

  ws->var("ctau3DErr")->setRange(ctauErrLow, ctauHigh);

  delete dsAB;
  delete dsAB1;
  delete datasetW;
  delete datasetWo;
  delete argSetW;
  delete argSetWo;
  delete f1;
  delete fMass;
}

void CtauErrFit::performSPlot()
{
  std::cout << "===== Start performSPlot() =====\n";
  RooRealVar *sigYield = ws->var("N_Jpsi");
  RooRealVar *bkgYield = ws->var("N_Bkg");

  sigYield->setMin(0); // sPlot requires a range: 0 ~ max
  bkgYield->setMin(0);

  RooArgList yieldList(*sigYield, *bkgYield);

  std::cout << "Sig Yield: " << sigYield->getVal() << " +/- " << sigYield->getError() << "\n";
  std::cout << "Bkg Yield: " << bkgYield->getVal() << " +/- " << bkgYield->getError() << "\n";

  RooDataSet *dataForSPlot = (RooDataSet *)ws->data("dsAB");

  sData = new RooStats::SPlot("sData", "", *dataForSPlot, ws->pdf("pdfMASS_Tot"), yieldList);
  ws->import(*dataForSPlot, Rename("dataset_SPLOT"));

  std::cout << "[INFO] Jpsi yield -> Mass Fit:" << ws->var("N_Jpsi")->getVal() << ", sWeights :" << sData->GetYieldFromSWeight("N_Jpsi") << "\n";
  std::cout << "[INFO] Bkg  yield -> Mass Fit:" << ws->var("N_Bkg")->getVal() << ", sWeights :" << sData->GetYieldFromSWeight("N_Bkg") << "\n";

  delete sigYield;
  delete bkgYield;
  // delete dataForSPlot;
}

void CtauErrFit::buildCtauErrModel()
{
  std::cout << "===== Start buildCtauErrModel() =====\n";
  // int nBins = std::min(int(round((ws->var("ctau3DErr")->getMax() - ws->var("ctau3DErr")->getMin()) / 0.0025)), 100);
  int nBins = nCtauErrBins;

  RooDataSet *dataw_Bkg_b = new RooDataSet("dataw_Bkg_b", "TMP_BKG_DATA", (RooDataSet *)ws->data("dataset_SPLOT"),
                                           RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Bkg_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("mass")), 0, "N_Bkg_sw");
  RooDataSet *dataw_Sig_b = new RooDataSet("dataw_Sig_b", "TMP_SIG_DATA", (RooDataSet *)ws->data("dataset_SPLOT"),
                                           RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Jpsi_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("mass")), 0, "N_Jpsi_sw");
  // TH1D *hSig = (TH1D *)dataw_Sig_b->createHistogram(("hSig"), *ws->var("ctau3DErr"), Binning(nBins, ctauErrLow, ctauErrHigh));

  // automatic min, max range
  TH1D *h_tmp = (TH1D *)ws->data("dsAB")->createHistogram(("h_tmp"), *ws->var("ctau3DErr"), Binning(nBins, ctauErrLow, ctauErrHigh));
  for (int i = 0; i < h_tmp->GetNbinsX() / 2; i++)
  {
    // if(hSig->GetBinContent(i)<=0&&hSig->GetBinContent(i+1)<=0&&hSig->GetBinContent(i+2)>=1&&hSig->GetBinContent(i+3)>=1){
    if (h_tmp->GetBinContent(i) > 1)
    { // pt 7-7.5
      // if(ptLow==3&&ptHigh==4.5&&cLow==20&&cHigh==120) ctauErrMin = hSig->GetBinLowEdge(i);
      // else ctauErrMin = hSig->GetBinLowEdge(i+1);
      ctauErrMin = h_tmp->GetBinLowEdge(i + 1);
      break;
    }
  }

  for (int i = 0; i < h_tmp->GetNbinsX(); i++)
  {
    if (h_tmp->GetBinContent(i) >= 1 && h_tmp->GetBinContent(i + 1) < 1 && h_tmp->GetBinContent(i + 2) < 1)
    {
      // ctauErrMax = hSig->GetBinLowEdge(i)+hSig->GetBinWidth(i);
      ctauErrMax = h_tmp->GetBinLowEdge(i) + h_tmp->GetBinWidth(i);
      break;
    }
    // else { ctauErrMax = hSig->GetBinLowEdge(i); }
    else
    {
      ctauErrMax = h_tmp->GetBinLowEdge(i);
    }
  }

  // Todo: user custom max.
  ctauErrMax = 0.1203;
  std::cout << "ctauErrMax : " << ctauErrMax << " ctauErrMin : " << ctauErrMin << "\n";

  double BinWidth = (ctauErrHigh - ctauErrLow) / nBins;
  std::cout << " BinWidth : " << BinWidth << "\n";
  newBins = (ctauErrMax - ctauErrMin) / BinWidth;
  std::cout << " newBins : " << newBins << "\n";

  // ===== make a datasets and HistPdf with new ctau3DErr range =====
  // total
  TH1D *hTot_w = (TH1D *)ws->data("dsAB")->createHistogram(("hTot_w"), *ws->var("ctau3DErr"), Binning(newBins, ctauErrMin, ctauErrMax));
  TH1D *hBkg_w = (TH1D *)dataw_Bkg_b->createHistogram(("hBkg_w"), *ws->var("ctau3DErr"), Binning(newBins, ctauErrMin, ctauErrMax));
  TH1D *hSig_w = (TH1D *)dataw_Sig_b->createHistogram(("hSig_w"), *ws->var("ctau3DErr"), Binning(newBins, ctauErrMin, ctauErrMax));

  // temporal RooDataHist
  RooDataHist *totHist = new RooDataHist("totHist", "", *ws->var("ctau3DErr"), hTot_w);
  RooDataHist *bkgHist = new RooDataHist("bkgHist", "", *ws->var("ctau3DErr"), hBkg_w);
  RooDataHist *sigHist = new RooDataHist("sigHist", "", *ws->var("ctau3DErr"), hSig_w);

  RooHistPdf *pdfCTAUERR_Tot = new RooHistPdf("pdfCTAUERR_Tot", "", *ws->var("ctau3DErr"), *totHist);
  RooHistPdf *pdfCTAUERR_Jpsi = new RooHistPdf("pdfCTAUERR_Jpsi", "", *ws->var("ctau3DErr"), *sigHist);
  RooHistPdf *pdfCTAUERR_Bkg = new RooHistPdf("pdfCTAUERR_Bkg", "", *ws->var("ctau3DErr"), *bkgHist);

  RooDataSet *dataw_Bkg = new RooDataSet("dataw_Bkg", "TMP_BKG_DATA", (RooDataSet *)ws->data("dataset_SPLOT"),
                                         RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Bkg_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("mass")), 0, "N_Bkg_sw");
  RooDataSet *dataw_Sig = new RooDataSet("dataw_Sig", "TMP_SIG_DATA", (RooDataSet *)ws->data("dataset_SPLOT"),
                                         RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Jpsi_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("mass")), 0, "N_Jpsi_sw");

  ws->import(*dataw_Sig);
  ws->import(*dataw_Bkg);
  ws->import(*pdfCTAUERR_Tot);
  ws->import(*pdfCTAUERR_Jpsi);
  ws->import(*pdfCTAUERR_Bkg);

  delete dataw_Bkg;
  delete dataw_Sig;
  delete pdfCTAUERR_Tot;
  delete pdfCTAUERR_Jpsi;
  delete pdfCTAUERR_Bkg;
  delete sigHist;
  delete bkgHist;
  delete totHist;
  delete hSig_w;
  delete hBkg_w;
  delete hTot_w;
  delete h_tmp;
  delete dataw_Sig_b;
  delete dataw_Bkg_b;
}

void CtauErrFit::drawResult()
{
  std::cout << "===== Start drawResult() =====\n";
  int nBins = nCtauErrBins;

  RooPlot *myPlot_B = ws->var("ctau3DErr")->frame(Range(ctauErrLow, ctauErrHigh));
  TH1D *hTot = (TH1D *)ws->data("dsAB")->createHistogram(("hTot"), *ws->var("ctau3DErr"),
                                                         // Binning(myPlot_B->GetNbinsX(), myPlot_B->GetXaxis()->GetXmin(), myPlot_B->GetXaxis()->GetXmax()));
                                                         Binning(myPlot_B->GetNbinsX(), myPlot_B->GetXaxis()->GetXmin(), myPlot_B->GetXaxis()->GetXmax()));

  double minRange = (double)(floor(ctauErrMin * 100.) / 100.);
  double maxRange = (double)(ceil(ctauErrMax * 100.) / 100.);
  ws->var("ctau3DErr")->setRange("ctauErrWindow", ctauErrMin, ctauErrMax);
  myPlot_B = ws->var("ctau3DErr")->frame(Bins(nBins), Range(minRange - 0.01, maxRange + 0.01)); // modified

  std::cout << ws->data("dsAB")->numEntries() << "\n";
  std::cout << ws->data("dataset_SPLOT")->numEntries() << "\n";

  TCanvas *c_B = new TCanvas("canvas_B", "My plots", 554, 4, 550, 520);
  c_B->cd();
  TPad *pad_B_1 = new TPad("pad_B_1", "pad_B_1", 0, 0.16, 0.98, 1.0);
  pad_B_1->SetTicks(1, 1);
  pad_B_1->Draw();
  pad_B_1->cd();
  myPlot_B->SetTitle("");

  c_B->cd();
  c_B->SetLogy();

  pad_B_1->cd();
  gPad->SetLogy();
  RooPlot *myPlot2_B = (RooPlot *)myPlot_B->Clone();
  ws->data("dsAB")->plotOn(myPlot2_B, Name("dataCTAUERR_Tot"), MarkerSize(.7), Binning(newBins)); // Normalization(wsmc->data("reducedDS_MC")->sumEntries()
  ws->pdf("pdfCTAUERR_Tot")->plotOn(myPlot2_B, Name("pdfCTAUERR_Tot"), LineColor(kGreen + 1), Range(ctauErrMin, ctauErrMax), LineWidth(2), Normalization(ws->data("dsAB")->sumEntries(), RooAbsReal::NumEvent));
  // ws->data("dataw_Sig")->plotOn(myPlot2_B,Name("dataHist_Sig"), MarkerSize(.7), LineColor(kRed+2), MarkerColor(kRed+2), Binning(nBins));
  ws->data("dataw_Sig")->plotOn(myPlot2_B, Name("dataHist_Sig"), MarkerSize(.7), LineColor(kRed + 2), MarkerColor(kRed + 2), Binning(newBins));
  ws->pdf("pdfCTAUERR_Jpsi")->plotOn(myPlot2_B, Name("pdfCTAUERR_Jpsi"), LineColor(kRed + 2), LineWidth(2), Range(ctauErrMin, ctauErrMax));
  // ws->data("dataw_Bkg")->plotOn(myPlot2_B,Name("dataHist_Bkg"), MarkerSize(.7), LineColor(kBlue+2), MarkerColor(kBlue+2), Binning(nBins));
  ws->data("dataw_Bkg")->plotOn(myPlot2_B, Name("dataHist_Bkg"), MarkerSize(.7), LineColor(kBlue + 2), MarkerColor(kBlue + 2), Binning(newBins));
  ws->pdf("pdfCTAUERR_Bkg")->plotOn(myPlot2_B, Name("pdfCTAUERR_Bkg"), LineColor(kBlue + 2), LineWidth(2), Range(ctauErrMin, ctauErrMax));
  Double_t YMax = hTot->GetBinContent(hTot->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= hTot->GetNbinsX(); i++)
    if (hTot->GetBinContent(i) > 0)
      YMin = std::min(YMin, hTot->GetBinContent(i));
  Double_t Yup(0.), Ydown(0.);
  Yup = YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.4 - 0.3)));
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.3 / (1.0 - 0.4 - 0.3))));
  myPlot2_B->GetYaxis()->SetRangeUser(Ydown, Yup);
  // std::cout<<ctauErrLow<<", "<<ctauErrHigh<<"\n";
  std::cout << ws->var("ctau3DErr")->getMin() << ", " << ws->var("ctau3DErr")->getMax() << "\n";
  RooDataSet *ctauResCutDS = (RooDataSet *)ws->data("dataw_Sig")->reduce(RooArgSet(*(ws->var("ctau3DRes")), *(ws->var("ctau3D")), *(ws->var("ctau3DErr"))), Form("ctau3DErr>%f&&ctau3DErr<%f", ctauErrMin, ctauErrMax));
  ctauResCutDS->SetName("ctauResCutDS");
  ws->import(*ctauResCutDS);
  Double_t outTot = ws->data("dsAB")->numEntries();
  Double_t outRes = ws->data("ctauResCutDS")->numEntries();
  std::cout << "Tot evt: (" << outTot << ")" << "\n";
  std::cout << "Res evt: (" << outRes << ")" << "\n";
  std::cout << "lost evt: " << (outTot - outRes) << " evts (" << ((outTot - outRes) * 100) / outTot << " %)" << "\n";

  TLine *minline = new TLine(ctauErrMin, 0.0, ctauErrMin, (Ydown * TMath::Power((Yup / Ydown), 0.5)));
  minline->SetLineStyle(2);
  minline->SetLineColor(1);
  minline->SetLineWidth(3);
  myPlot2_B->addObject(minline);
  TLine *maxline = new TLine(ctauErrMax, 0.0, ctauErrMax, (Ydown * TMath::Power((Yup / Ydown), 0.5)));
  maxline->SetLineStyle(2);
  maxline->SetLineColor(1);
  maxline->SetLineWidth(3);
  myPlot2_B->addObject(maxline);

  myPlot2_B->GetXaxis()->CenterTitle();
  myPlot2_B->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} Error (mm)");
  myPlot2_B->SetFillStyle(4000);
  myPlot2_B->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_B->GetXaxis()->SetLabelSize(0);
  myPlot2_B->GetXaxis()->SetTitleSize(0);
  myPlot2_B->Draw();
  // Double_t outTot = ws->data("dsAB")->numEntries();
  // Double_t outErr = ws->data("dsAB")->reduce(Form("(ctauErr>=%.6f || ctauErr<=%.6f)", ctauErrHigh, ctauErrLow))->numEntries();
  // std::cout<<(outErr*100)/outTot<<"\n";
  TLegend *leg_B = new TLegend(text_x + 0.5, text_y - 0.2, text_x + 0.7, text_y);
  leg_B->SetTextSize(text_size);
  leg_B->SetTextFont(43);
  leg_B->SetBorderSize(0);
  leg_B->AddEntry(myPlot2_B->findObject("dataCTAUERR_Tot"), "Data", "pe");
  leg_B->AddEntry(myPlot2_B->findObject("pdfCTAUERR_Tot"), "Total PDF", "l");
  leg_B->AddEntry(myPlot2_B->findObject("pdfCTAUERR_Jpsi"), "Signal", "l");
  leg_B->AddEntry(myPlot2_B->findObject("pdfCTAUERR_Bkg"), "Background", "l");
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
  c_B->cd();
  pad_B_2->Draw();
  pad_B_2->cd();
  pad_B_2->SetTopMargin(0); // Upper and lower plot are joined
  pad_B_2->SetBottomMargin(0.67);
  pad_B_2->SetBottomMargin(0.4);
  pad_B_2->SetFillStyle(4000);
  pad_B_2->SetFrameFillStyle(4000);
  pad_B_2->SetTicks(1, 1);

  RooPlot *frameTMP_B = (RooPlot *)myPlot2_B->Clone("TMP");
  RooHist *hpull_B = frameTMP_B->pullHist("dataCTAUERR_Tot", "pdfCTAUERR_Tot");
  // RooHist* hpull_B = frameTMP_B->pullHist("dataCTAUERR_Tot","pdfCTAUERR_Tot");
  // RooHist* hpull_B = frameTMP_B->pullHist("dataHist_Sig","pdfCTAUERR_Jpsi");
  hpull_B->SetMarkerSize(0.8);
  // RooPlot* pullFrame_B = ws->var("ctau3DErr")->frame(Title("Pull Distribution")) ;
  RooPlot *pullFrame_B = ws->var("ctau3DErr")->frame(Bins(nBins), Range(minRange - 0.01, maxRange + 0.01));
  pullFrame_B->addPlotable(hpull_B, "PX");
  pullFrame_B->SetTitle("");
  pullFrame_B->SetTitleSize(0);
  pullFrame_B->GetYaxis()->SetTitleOffset(0.3);
  pullFrame_B->GetYaxis()->SetTitle("Pull");
  pullFrame_B->GetYaxis()->SetTitleSize(0.15);
  pullFrame_B->GetYaxis()->SetLabelSize(0.15);
  pullFrame_B->GetYaxis()->SetRangeUser(-3.8, 3.8);
  pullFrame_B->GetYaxis()->CenterTitle();

  pullFrame_B->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} Error (mm)");
  pullFrame_B->GetXaxis()->SetTitleOffset(1.05);
  pullFrame_B->GetXaxis()->SetLabelOffset(0.04);
  pullFrame_B->GetXaxis()->SetLabelSize(0.15);
  pullFrame_B->GetXaxis()->SetTitleSize(0.15);
  pullFrame_B->GetXaxis()->CenterTitle();

  pullFrame_B->GetYaxis()->SetTickSize(0.04);
  pullFrame_B->GetYaxis()->SetNdivisions(404);
  pullFrame_B->GetXaxis()->SetTickSize(0.03);
  pullFrame_B->Draw();

  TLine *lB = new TLine(minRange - 0.01, 0, maxRange + 0.01, 0);
  lB->SetLineStyle(1);
  lB->Draw("same");
  pad_B_2->Update();

  c_B->Update();
  c_B->SaveAs(Form("figs/2DFit_%s/CtauErr/ctauErr_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.png", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));

  delete ctauResCutDS;
  delete hTot;
  delete lB;
  // delete leg_B;
  // delete pullFrame_B;
  // delete hpull_B;
  delete frameTMP_B;
  delete myPlot2_B; // get maxline, minline
  delete myPlot_B;
  delete c_B;
}

void CtauErrFit::saveResult()
{
  std::cout << "===== Start saveResult() =====\n";
  TFile *outFile = new TFile(Form("roots/2DFit_%s/CtauErr/CtauErrResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP), "recreate");

  ws->SetName("old_ws");
  
  RooWorkspace newWs("ws_err", "");
  // newWs.import(*ws->pdf("pdfCTAUERR_Tot"));
  newWs.import(*ws->pdf("pdfCTAUERR_Bkg"));
  newWs.import(*ws->pdf("pdfCTAUERR_Jpsi"));
  newWs.import(*ws->data("dataw_Sig"));
  newWs.import(*ws->data("dataw_Bkg"));  
  newWs.Write();

  outFile->Close();
  delete outFile;
}