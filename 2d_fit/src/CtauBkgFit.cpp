#include "CtauBkgFit.h"
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH1D.h"

#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

#include "../headers/polarizationUtilities.h"

using namespace RooFit;

CtauBkgFit::CtauBkgFit(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
                       float cosLow, float cosHigh, int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP)
    : ptLow(ptLow), ptHigh(ptHigh), yLow(yLow), yHigh(yHigh), cLow(cLow), cHigh(cHigh),
      cosLow(cosLow), cosHigh(cosHigh), PR(PR), PRw(PRw),
      fEffW(fEffW), fAccW(fAccW), isPtW(isPtW), isTnP(isTnP)
{
  gStyle->SetEndErrorSize(0);

  RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling);

  // to avoid ternimal breaking due to too long Tracing Error message
  // It's okay as long as minimizer find good fitting point
  // Always check fit results.
  // RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Tracing);
  RooMsgService::instance().getStream(1).removeTopic(Tracing);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
}

CtauBkgFit::~CtauBkgFit()
{
  RooMsgService::instance().getStream(1).addTopic(RooFit::ObjectHandling);

  RooMsgService::instance().getStream(0).addTopic(Caching);
  RooMsgService::instance().getStream(1).addTopic(Caching);
  RooMsgService::instance().getStream(0).addTopic(Plotting);
  RooMsgService::instance().getStream(1).addTopic(Plotting);
  RooMsgService::instance().getStream(0).addTopic(Tracing);
  RooMsgService::instance().getStream(1).addTopic(Tracing);
  RooMsgService::instance().getStream(0).addTopic(Integration);
  RooMsgService::instance().getStream(1).addTopic(Integration);
  RooMsgService::instance().getStream(1).addTopic(NumIntegration);

  delete fMass;
  delete fCErr;
  delete fCRes;
  delete ws;
  delete fitResult;
}

void CtauBkgFit::run()
{
  setLabels();
  openInputFile();
  setupWorkspaceAndData();
  setVariableRanges();
  defineModel();
  performFit();
  drawPullPlot();
  drawRatioPlot();
  saveResults();
}

void CtauBkgFit::setLabels()
{
  std::cout << "===== setLabels() =====\n\n";
  DATE = "No_Weight_2";
  gSystem->mkdir(Form("roots/2DFit_%s/CtauBkg", DATE.c_str()), kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/CtauBkg", DATE.c_str()), kTRUE);

  kineLabel = std::string(getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh).Data());

  if (PRw == 1)
    fname = "PR";
  else if (PRw == 2)
    fname = "NP";
}

void CtauBkgFit::openInputFile()
{
  std::cout << "===== openInputFile() =====\n\n";

  fMass = new TFile(Form("roots/2DFit_%s/Mass/Mass_FixedFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));
  fCErr = new TFile(Form("roots/2DFit_%s/CtauErr/CtauErrResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));
  fCRes = new TFile(Form("roots/2DFit_%s/CtauRes/CtauResResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));

  if (!fMass || fMass->IsZombie())
  {
    std::cerr << "CRITICAL: Mass fit result file could not be opened. Aborting.\n";
    exit(1);
  }
  if (!fCErr || fCErr->IsZombie())
  {
    std::cerr << "CRITICAL: Ctau error result file could not be opened. Aborting.\n";
    exit(1);
  }
  if (!fCRes || fCRes->IsZombie())
  {
    std::cerr << "CRITICAL: Ctau resolution result file could not be opened. Aborting.\n";
    exit(1);
  }
}

void CtauBkgFit::setupWorkspaceAndData()
{
  std::cout << "===== setupWorkspaceAndData() =====\n\n";

  RooDataSet *datasetMass = (RooDataSet *)fMass->Get("datasetMass");
  RooAddPdf *pdfMASS_Tot = (RooAddPdf *)fMass->Get("pdfMASS_Tot");
  RooDataSet *dataw_Bkg = (RooDataSet *)fCErr->Get("dataw_Bkg");
  RooAddPdf *GaussModel_Tot = (RooAddPdf *)fCRes->Get("GaussModel_Tot");
  RooAddPdf *pdfCTAUERR_Bkg = (RooAddPdf *)fCErr->Get("pdfCTAUERR_Bkg");

  ws = new RooWorkspace("workspace");

  ws->import(*datasetMass);
  ws->import(*pdfMASS_Tot);
  ws->import(*dataw_Bkg);
  ws->import(*GaussModel_Tot);
  ws->import(*pdfCTAUERR_Bkg);
}

void CtauBkgFit::setVariableRanges()
{
  std::cout << "===== setVariableRanges() =====\n\n";
  double ctauErrMin, ctauErrMax;
  ws->data("dataw_Bkg")->getRange(*ws->var("ctau3DErr"), ctauErrMin, ctauErrMax);

  ws->var("ctau3D")->setRange(ctauLow, ctauHigh);
  ws->var("ctau3D")->setRange("bkgPlotRange", ctauLow, ctauHigh);
  
  ws->var("ctau3DErr")->setRange(ctauErrMin, ctauErrMax);
  ws->var("ctau3DErr")->setRange("bkgErrRange", ctauErrMin, ctauErrMax);
}

void CtauBkgFit::defineModel()
{
  std::cout << "===== defineModel() =====\n\n";
  // set parameters
  ws->factory("zeroMean[0.0]");
  ws->factory("b_Bkg[0.4, 0.1, 1]"); // NP fraction for bkg
  ws->factory("fDFSS[0.3, 0.05, 1.]");
  ws->factory("fDLIV[0.3, 0.05, 1.0]");
  ws->factory("lambdaDDS_Bkg[1, 0.001, 10]");
  ws->factory("lambdaDF_Bkg1[1, 0.001, 10]");
  ws->factory("lambdaDF_Bkg2[1, 0.001, 10]");
  ws->factory("lambdaDSS_Bkg1[1, 0.001, 10]");
  ws->factory("lambdaDSS_Bkg2[1, 0.001, 10]");
  ws->factory("fDSS12[0.3, 0.01, 1.]");
  ws->factory("fDF12[0.04, 0.01, 1.]");

  // Todo: remove this part. Res code will fix the parameters after fit
  ws->var("ctau1_CtauRes")->setConstant(kTRUE);
  ws->var("s1_CtauRes")->setConstant(kTRUE);
  ws->var("ctau2_CtauRes")->setConstant(kTRUE);
  ws->var("rS21_CtauRes")->setConstant(kTRUE);
  ws->var("f_CtauRes")->setConstant(kTRUE);

  // ===== resolution model =====
  ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes1", "ctau3D",
                   "ctau1_CtauRes", //"ctau1_CtauRes",
                   "s1_CtauRes",
                   "zeroMean",
                   "ctau3DErr"));
  ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes2", "ctau3D",
                   "ctau2_CtauRes", //"ctau2_CtauRes",
                   "s2_CtauRes",
                   "zeroMean",
                   "ctau3DErr"));
  if (nGauss == 3)
  {
    ws->var("ctau3_CtauRes")->setConstant(kTRUE);
    ws->var("rS32_CtauRes")->setConstant(kTRUE);
    ws->var("f2_CtauRes")->setConstant(kTRUE);
    ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes3", "ctau3D",
                     "ctau3_CtauRes", //"ctau3_CtauRes",
                     "s3_CtauRes",
                     "zeroMean",
                     "ctau3DErr"));
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "ctauRes32", "ctauRes3", "ctauRes2", "f2_CtauRes"));
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "pdfCTAURES", "ctauRes1", "ctauRes32", "f_CtauRes"));
  }
  else
  {
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "pdfCTAURES", "ctauRes1", "ctauRes2", "f_CtauRes"));
  }

  // ===== lifetime model =====
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUDSS1", "ctau3D", "lambdaDSS_Bkg1", "pdfCTAURES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUDSS2", "ctau3D", "lambdaDSS_Bkg2", "pdfCTAURES"));
  // ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", "pdfCTAUDF", "ctau3D", "lambdaDF_Bkg1", "pdfCTAURES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", "pdfCTAUDF1", "ctau3D", "lambdaDF_Bkg1", "pdfCTAURES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", "pdfCTAUDF2", "ctau3D", "lambdaDF_Bkg2", "pdfCTAURES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", "pdfCTAUDDS", "ctau3D", "lambdaDDS_Bkg", "pdfCTAURES"));

  // combine
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUDSS", "fDSS12", "pdfCTAUDSS1", "pdfCTAUDSS2"));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUDF", "fDF12", "pdfCTAUDF1", "pdfCTAUDF2"));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAU1", "fDFSS", "pdfCTAUDSS", "pdfCTAUDF"));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUCOND_BkgNoPR", "fDLIV", "pdfCTAU1", "pdfCTAUDDS")); // NP
  ws->factory(Form("SUM::%s(%s)", "pdfCTAUCOND_BkgPR", "pdfCTAURES"));                               // PR
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUCOND_Bkg", "b_Bkg", "pdfCTAUCOND_BkgNoPR", "pdfCTAUCOND_BkgPR"));

  // model for fit
  ws->factory(Form("RooExtendPdf::%s(%s,%s)", "pdfTot_Bkg", "pdfCTAUCOND_Bkg", "N_Bkg")); // N_Bkg is number of bkg from dataw_Bkg

  // ===== conditional model =====
  // To check goodness of fit -> NOT USED FOR FIT HERE.
  RooProdPdf pdfPR("pdfCTAU_BkgPR", "", *ws->pdf("pdfCTAUERR_Bkg"), Conditional(*ws->pdf("pdfCTAUCOND_BkgPR"), RooArgList(*ws->var("ctau3D"))));
  ws->import(pdfPR);
  
  RooProdPdf pdfNoPR("pdfCTAU_BkgNoPR", "", *ws->pdf("pdfCTAUERR_Bkg"), Conditional(*ws->pdf("pdfCTAUCOND_BkgNoPR"), RooArgList(*ws->var("ctau3D"))));
  ws->import(pdfNoPR);

  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAU_Bkg", "b_Bkg", "pdfCTAU_BkgNoPR", "pdfCTAU_BkgPR"));
  RooAbsPdf *pdfCTAU_Bkg_Tot = new RooAddPdf("pdfCTAU_Bkg_Tot", "pdfCTAU_Bkg_Tot", RooArgList(*ws->pdf("pdfCTAU_Bkg")), RooArgList(*ws->var("N_Bkg")));

  ws->import(*pdfCTAU_Bkg_Tot);

  delete pdfCTAU_Bkg_Tot;
}

void CtauBkgFit::performFit()
{
  std::cout << "===== performFit() =====\n\n";
  
  // set proper fit range
  TH1D *hTot = (TH1D *)ws->data("dataw_Bkg")->createHistogram(("hTot"), *ws->var("ctau3D"), Binning(nCtauBins, ctauLow, ctauHigh));
  if (ctauMin == -100)
    ctauMin = hTot->GetBinLowEdge(hTot->FindFirstBinAbove(1, 1));
  if (ctauMax == 100)
    ctauMax = hTot->GetBinLowEdge(hTot->FindLastBinAbove(2, 1)) + hTot->GetBinWidth(hTot->FindLastBinAbove(2, 1));
  delete hTot;
  ws->var("ctau3D")->setRange("bkgFitRange", ctauMin, ctauMax);

  // make a dataset for fitting
  RooDataSet *dataToFit = (RooDataSet *)ws->data("dataw_Bkg")->reduce(CutRange("bkgFitRange"));
  ws->import(*dataToFit, Rename("dataToFit"));

  // fit
  bool isWeighted = ws->data("dataToFit")->isWeighted();
  fitResult = ws->pdf("pdfTot_Bkg")->fitTo(*ws->data("dataToFit"), Save(), Range("bkgFitRange"), Extended(kTRUE), NumCPU(nCPU), PrintLevel(0), SumW2Error(isWeighted), RecoverFromUndefinedRegions(1.), PrintEvalErrors(-1), Strategy(2));
  fitResult->Print("V");

  std::cout << "fixing floating parameters to constant (after fitting)\n";
  const RooArgList &floatted_params = fitResult->floatParsFinal();
  for (auto arg : floatted_params)
  {
    RooRealVar *param = dynamic_cast<RooRealVar *>(arg);
    if (param)
    {
      std::cout << "Fixing parameter: " << param->GetName() << "\n";
      ws->var(param->GetName())->setConstant(kTRUE);
    }
  }

  delete dataToFit;
}

void CtauBkgFit::drawPullPlot()
{
  std::cout << "===== Start drawPullPlot() =====\n";
  TCanvas *myCanvas = new TCanvas("myCanvas", "", 800, 800);
  myCanvas->cd();
  TPad *ctauPad = new TPad("ctauPad", "Main plot pad", 0.0, 0.3, 1.0, 1.0);
  ctauPad->SetTicks(1, 1);
  ctauPad->SetLeftMargin(0.15);
  ctauPad->SetRightMargin(0.07);
  ctauPad->SetTopMargin(0.1);
  ctauPad->SetBottomMargin(0.02);
  ctauPad->Draw();
  ctauPad->cd();
  gPad->SetLogy();

  RooPlot *ctauFrame = ws->var("ctau3D")->frame(Bins(nCtauBins), Range(ctauLow, ctauHigh)); // bins
  ctauFrame->SetTitle("");

  // Todo: remove legacy?
  double normDSTot = 1.0;
  if (ws->data("dataToFit"))
  {
    normDSTot = ws->data("dataToFit")->sumEntries() / ws->data("dataToFit")->sumEntries();
  }
  double normBkg = 1.0;
  if (ws->data("dataw_Bkg"))
  {
    normBkg = ws->data("dataToFit")->sumEntries() * normDSTot / ws->data("dataw_Bkg")->sumEntries();
  }

  ws->pdf("pdfCTAU_Bkg_Tot")->setNormRange("ctaubkgFitRange");
  ctauFrame->updateNormVars(RooArgSet(*ws->var("ctau3D")));
  // ctauFrame->updateNormVars(RooArgSet(*ws->var("mass"), *ws->var("ctau3D"), *ws->var("ctau3DErr")));

  // plotting
  ws->data("dataToFit")->plotOn(ctauFrame, Name("data_points"), DataError(RooAbsData::SumW2)); // for normalization
  ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("total_pdf"), ProjWData(*ws->data("dataw_Bkg")));
  ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("prompt_pdf"), Components(RooArgSet(*ws->pdf("pdfCTAU_BkgPR"))), LineStyle(kDashed), LineColor(kGreen + 2), LineWidth(2), LineStyle(kDashed));
  ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("nonprompt_pdf"), Components(RooArgSet(*ws->pdf("pdfCTAU_BkgNoPR"))), LineStyle(kDashed), LineColor(kRed + 1), LineWidth(2), LineStyle(kDashed));
  
  ws->data("dataToFit")->plotOn(ctauFrame, Name("data_points"), DataError(RooAbsData::SumW2)); // to draw data points over PDFs

  // set Y range
  TH1 *h_tmp = ws->data("dataw_Bkg")->createHistogram("h_tmp", *ws->var("ctau3D"), Binning(ctauFrame->GetNbinsX(), ctauFrame->GetXaxis()->GetXmin(), ctauFrame->GetXaxis()->GetXmax()));
  Double_t YMax = h_tmp->GetBinContent(h_tmp->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= h_tmp->GetNbinsX(); i++)
    if (h_tmp->GetBinContent(i) > 0)
      YMin = std::min(YMin, h_tmp->GetBinContent(i));
  Double_t Yup(0.), Ydown(0.);
  // Yup = YMax*TMath::Power((YMax/0.1), 0.5);
  Yup = YMax * TMath::Power((YMax / 0.01), 0.5);
  Ydown = 0.01;

  // ctauFrame->GetYaxis()->SetRangeUser(10e-2, 10e7);
  ctauFrame->GetYaxis()->SetRangeUser(Ydown, Yup);
  ctauFrame->GetXaxis()->SetRangeUser(ctauLow, ctauHigh);
  ctauFrame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  ctauFrame->SetFillStyle(4000);
  ctauFrame->GetYaxis()->SetTitleOffset(1.43);
  ctauFrame->GetXaxis()->SetLabelSize(0);
  ctauFrame->GetXaxis()->SetTitleSize(0);
  ctauFrame->Draw();

  TLine *minline = new TLine(ctauMin, 0.0, ctauMin, Ydown * TMath::Power((Yup / Ydown), 0.4));
  minline->SetLineStyle(2);
  minline->SetLineColor(1);
  minline->SetLineWidth(3);
  ctauFrame->addObject(minline);
  TLine *maxline = new TLine(ctauMax, 0.0, ctauMax, Ydown * TMath::Power((Yup / Ydown), 0.4));
  maxline->SetLineStyle(2);
  maxline->SetLineColor(1);
  maxline->SetLineWidth(3);
  ctauFrame->addObject(maxline);

  // ctauFrame->Draw();
  TLegend *leg_E = new TLegend(text_x + 0.3, text_y + 0.04, text_x + 0.43, text_y - 0.15);
  leg_E->SetTextSize(text_size);
  leg_E->SetTextFont(43);
  leg_E->SetBorderSize(0);
  leg_E->AddEntry(ctauFrame->findObject("data_points"), "Data_Bkg", "pe");
  leg_E->AddEntry(ctauFrame->findObject("total_pdf"), "Total PDF", "fl");
  // leg_E->AddEntry(ctauFrame->findObject("test"),"?? PDF","l");
  leg_E->Draw("same");
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x + 0.05, text_y, text_color, text_size);
  if (yLow == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x + 0.05, text_y - y_diff, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x + 0.05, text_y - y_diff, text_color, text_size);
  drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x + 0.05, text_y - y_diff * 2, text_color, text_size);

  drawText(Form("N_{Bkg} = %.f #pm %.f", ws->var("N_Bkg")->getVal(), ws->var("N_Bkg")->getError()), text_x + 0.55, text_y, text_color, text_size);
  drawText(Form("b_{Bkg} = %.4f #pm %.4f", ws->var("b_Bkg")->getVal(), ws->var("b_Bkg")->getError()), text_x + 0.55, text_y - y_diff * 1, text_color, text_size);
  drawText(Form("fDFSS = %.4f #pm %.4f", ws->var("fDFSS")->getVal(), ws->var("fDFSS")->getError()), text_x + 0.55, text_y - y_diff * 2, text_color, text_size);
  drawText(Form("fDLIV = %.4f #pm %.4f", ws->var("fDLIV")->getVal(), ws->var("fDLIV")->getError()), text_x + 0.55, text_y - y_diff * 3, text_color, text_size);
  drawText(Form("#lambdaDDS_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDDS_Bkg")->getVal(), ws->var("lambdaDDS_Bkg")->getError()), text_x + 0.55, text_y - y_diff * 4, text_color, text_size);
  drawText(Form("#lambdaDF1_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDF_Bkg1")->getVal(), ws->var("lambdaDF_Bkg1")->getError()), text_x + 0.55, text_y - y_diff * 5, text_color, text_size);
  drawText(Form("#lambdaDF2_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDF_Bkg2")->getVal(), ws->var("lambdaDF_Bkg2")->getError()), text_x + 0.55, text_y - y_diff * 6, text_color, text_size);
  drawText(Form("#lambdaDSS1_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDSS_Bkg1")->getVal(), ws->var("lambdaDSS_Bkg1")->getError()), text_x + 0.55, text_y - y_diff * 7, text_color, text_size);
  drawText(Form("#lambdaDSS2_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDSS_Bkg2")->getVal(), ws->var("lambdaDSS_Bkg2")->getError()), text_x + 0.55, text_y - y_diff * 8, text_color, text_size);

  TPad *pullPad = new TPad("pullPad", "Pull plot pad", 0.0, 0.0, 1.0, 0.3);
  myCanvas->cd();
  pullPad->SetLeftMargin(0.15);
  pullPad->SetRightMargin(0.07);
  pullPad->SetTopMargin(0); // Upper and lower plot are joined
  // pullPad->SetTopMargin(0.05);
  pullPad->SetBottomMargin(0.35);
  pullPad->SetTicks(1, 1);
  pullPad->SetFillStyle(4000);
  pullPad->SetFrameFillStyle(4000);
  pullPad->Draw();
  pullPad->cd();

  RooPlot *frameTMP = (RooPlot *)ctauFrame->Clone("frameTMP");
  RooHist *pullHist = frameTMP->pullHist("data_points", "total_pdf", true);
  pullHist->SetMarkerSize(0.8);
  RooPlot *pullFrame = ws->var("ctau3D")->frame(Title("Pull Distribution"));
  pullFrame->addPlotable(pullHist, "P");
  pullFrame->SetTitle("");
  pullFrame->SetTitleSize(0);
  pullFrame->GetYaxis()->SetTitleOffset(0.3);
  pullFrame->GetYaxis()->SetTitle("Pull");
  pullFrame->GetYaxis()->SetTitleSize(0.08);
  pullFrame->GetYaxis()->SetLabelSize(0.08);
  pullFrame->GetYaxis()->SetRangeUser(-10, 10);
  pullFrame->GetYaxis()->CenterTitle();

  pullFrame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  pullFrame->GetXaxis()->SetTitleOffset(1.05);
  pullFrame->GetXaxis()->SetLabelOffset(0.04);
  pullFrame->GetXaxis()->SetLabelSize(0.08);
  pullFrame->GetXaxis()->SetTitleSize(0.08);
  pullFrame->GetXaxis()->CenterTitle();

  pullFrame->GetYaxis()->SetTickSize(0.04);
  pullFrame->GetYaxis()->SetNdivisions(404);
  pullFrame->GetXaxis()->SetTickSize(0.03);
  pullFrame->Draw();

  TLine *l1 = new TLine(ctauLow, 0, ctauHigh, 0);
  l1->SetLineStyle(1);
  l1->Draw("same");

  printChi2(ws, pullPad, frameTMP, fitResult, "ctau3D", "data_points", "total_pdf", nCtauBins, false);

  Double_t theNLL = fitResult->minNll();
  std::cout << " *** NLL : " << std::setprecision(15) << theNLL << "\n";

  myCanvas->Update();
  myCanvas->Draw();
  myCanvas->SaveAs(Form("figs/2DFit_%s/CtauBkg/Bkg_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.png", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));

  delete h_tmp;
  delete frameTMP;
  delete myCanvas;
}

void CtauBkgFit::drawRatioPlot()
{
  std::cout << "===== Start drawRatioPlot() =====\n";
  TCanvas *myCanvas = new TCanvas("myCanvas", "", 800, 800);
  myCanvas->cd();
  TPad *ctauPad = new TPad("ctauPad", "Main plot pad", 0.0, 0.3, 1.0, 1.0);
  ctauPad->SetTicks(1, 1);
  ctauPad->SetLeftMargin(0.15);
  ctauPad->SetRightMargin(0.07);
  ctauPad->SetTopMargin(0.1);
  ctauPad->SetBottomMargin(0.02);
  ctauPad->Draw();
  ctauPad->cd();
  gPad->SetLogy();

  RooPlot *ctauFrame = ws->var("ctau3D")->frame(Bins(nCtauBins), Range(ctauLow, ctauHigh)); // bins
  ctauFrame->SetTitle("");

  // Todo: remove legacy?
  double normDSTot = 1.0;
  if (ws->data("dataToFit"))
  {
    normDSTot = ws->data("dataToFit")->sumEntries() / ws->data("dataToFit")->sumEntries();
  }
  double normBkg = 1.0;
  if (ws->data("dataw_Bkg"))
  {
    normBkg = ws->data("dataToFit")->sumEntries() * normDSTot / ws->data("dataw_Bkg")->sumEntries();
  }

  ws->pdf("pdfCTAU_Bkg_Tot")->setNormRange("ctaubkgFitRange");
  ctauFrame->updateNormVars(RooArgSet(*ws->var("ctau3D")));
  // ctauFrame->updateNormVars(RooArgSet(*ws->var("mass"), *ws->var("ctau3D"), *ws->var("ctau3DErr")));

  // plotting
  ws->data("dataToFit")->plotOn(ctauFrame, Name("data_points"), DataError(RooAbsData::SumW2));
  ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("total_pdf"), ProjWData(*ws->data("dataw_Bkg")));
  ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("prompt_pdf"), Components(RooArgSet(*ws->pdf("pdfCTAU_BkgPR"))), LineStyle(kDashed), LineColor(kGreen + 2), LineWidth(2), LineStyle(kDashed));
  ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(ctauFrame, Name("nonprompt_pdf"), Components(RooArgSet(*ws->pdf("pdfCTAU_BkgNoPR"))), LineStyle(kDashed), LineColor(kRed + 1), LineWidth(2), LineStyle(kDashed));

  ws->data("dataToFit")->plotOn(ctauFrame, Name("data_points"), DataError(RooAbsData::SumW2)); // redraw data points on the pdfs

  // set Y range
  TH1 *h_tmp = ws->data("dataw_Bkg")->createHistogram("h_tmp", *ws->var("ctau3D"), Binning(ctauFrame->GetNbinsX(), ctauFrame->GetXaxis()->GetXmin(), ctauFrame->GetXaxis()->GetXmax()));
  Double_t YMax = h_tmp->GetBinContent(h_tmp->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= h_tmp->GetNbinsX(); i++)
    if (h_tmp->GetBinContent(i) > 0)
      YMin = std::min(YMin, h_tmp->GetBinContent(i));
  Double_t Yup(0.), Ydown(0.);
  // Yup = YMax*TMath::Power((YMax/0.1), 0.5);
  Yup = YMax * TMath::Power((YMax / 0.01), 0.5);
  Ydown = 0.01;

  // ctauFrame->GetYaxis()->SetRangeUser(10e-2, 10e7);
  ctauFrame->GetYaxis()->SetRangeUser(Ydown, Yup);
  ctauFrame->GetXaxis()->SetRangeUser(ctauLow, ctauHigh);
  ctauFrame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  ctauFrame->SetFillStyle(4000);
  ctauFrame->GetYaxis()->SetTitleOffset(1.43);
  ctauFrame->GetXaxis()->SetLabelSize(0);
  ctauFrame->GetXaxis()->SetTitleSize(0);
  ctauFrame->Draw();

  TLine *minline = new TLine(ctauMin, 0.0, ctauMin, Ydown * TMath::Power((Yup / Ydown), 0.4));
  minline->SetLineStyle(2);
  minline->SetLineColor(1);
  minline->SetLineWidth(3);
  ctauFrame->addObject(minline);
  TLine *maxline = new TLine(ctauMax, 0.0, ctauMax, Ydown * TMath::Power((Yup / Ydown), 0.4));
  maxline->SetLineStyle(2);
  maxline->SetLineColor(1);
  maxline->SetLineWidth(3);
  ctauFrame->addObject(maxline);

  // ctauFrame->Draw();
  TLegend *leg_E = new TLegend(text_x + 0.3, text_y + 0.04, text_x + 0.43, text_y - 0.15);
  leg_E->SetTextSize(text_size);
  leg_E->SetTextFont(43);
  leg_E->SetBorderSize(0);
  leg_E->AddEntry(ctauFrame->findObject("data_points"), "Data_Bkg", "pe");
  leg_E->AddEntry(ctauFrame->findObject("total_pdf"), "Total PDF", "fl");
  // leg_E->AddEntry(ctauFrame->findObject("test"),"?? PDF","l");
  leg_E->Draw("same");
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x + 0.05, text_y, text_color, text_size);
  if (yLow == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x + 0.05, text_y - y_diff, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x + 0.05, text_y - y_diff, text_color, text_size);
  drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x + 0.05, text_y - y_diff * 2, text_color, text_size);

  drawText(Form("N_{Bkg} = %.f #pm %.f", ws->var("N_Bkg")->getVal(), ws->var("N_Bkg")->getError()), text_x + 0.55, text_y, text_color, text_size);
  drawText(Form("b_{Bkg} = %.4f #pm %.4f", ws->var("b_Bkg")->getVal(), ws->var("b_Bkg")->getError()), text_x + 0.55, text_y - y_diff * 1, text_color, text_size);
  drawText(Form("fDFSS = %.4f #pm %.4f", ws->var("fDFSS")->getVal(), ws->var("fDFSS")->getError()), text_x + 0.55, text_y - y_diff * 2, text_color, text_size);
  drawText(Form("fDLIV = %.4f #pm %.4f", ws->var("fDLIV")->getVal(), ws->var("fDLIV")->getError()), text_x + 0.55, text_y - y_diff * 3, text_color, text_size);
  drawText(Form("#lambdaDDS_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDDS_Bkg")->getVal(), ws->var("lambdaDDS_Bkg")->getError()), text_x + 0.55, text_y - y_diff * 4, text_color, text_size);
  drawText(Form("#lambdaDF1_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDF_Bkg1")->getVal(), ws->var("lambdaDF_Bkg1")->getError()), text_x + 0.55, text_y - y_diff * 5, text_color, text_size);
  drawText(Form("#lambdaDF2_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDF_Bkg2")->getVal(), ws->var("lambdaDF_Bkg2")->getError()), text_x + 0.55, text_y - y_diff * 6, text_color, text_size);
  drawText(Form("#lambdaDSS1_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDSS_Bkg1")->getVal(), ws->var("lambdaDSS_Bkg1")->getError()), text_x + 0.55, text_y - y_diff * 7, text_color, text_size);
  drawText(Form("#lambdaDSS2_{Bkg} = %.4f #pm %.4f", ws->var("lambdaDSS_Bkg2")->getVal(), ws->var("lambdaDSS_Bkg2")->getError()), text_x + 0.55, text_y - y_diff * 8, text_color, text_size);

  TPad *ratioPad = new TPad("ratioPad", "", 0.0, 0.0, 1.0, 0.3);
  myCanvas->cd();
  ratioPad->SetLeftMargin(0.15);
  ratioPad->SetRightMargin(0.07);
  ratioPad->SetTopMargin(0); // Upper and lower plot are joined
  // ratioPad->SetTopMargin(0.05);
  ratioPad->SetBottomMargin(0.35);
  ratioPad->SetTicks(1, 1);
  ratioPad->SetFillStyle(4000);
  ratioPad->SetFrameFillStyle(4000);
  ratioPad->Draw();
  ratioPad->cd();

  RooHist *h_dataPoints = (RooHist *)ctauFrame->findObject("data_points");
  RooCurve *curveRatio = (RooCurve *)ctauFrame->findObject("total_pdf");

  RooPlot *ratioFrame = ws->var("ctau3D")->frame(Title("Ratio"), Range(ctauLow, ctauHigh));
  ratioFrame->SetTitle("");

  TGraphAsymmErrors *g_ratio = new TGraphAsymmErrors();
  int point_idx = 0;
  for (int i = 0; i < h_dataPoints->GetN(); ++i)
  {
    double x, y;
    h_dataPoints->GetPoint(i, x, y);
    double model_val = curveRatio->Eval(x);
    if (model_val > 1e-9 && y > 0)
    {
      double ratio = y / model_val;
      g_ratio->SetPoint(point_idx, x, ratio);
      double err_y_low = h_dataPoints->GetErrorYlow(i);
      double err_y_high = h_dataPoints->GetErrorYhigh(i);
      g_ratio->SetPointError(point_idx, 0, 0, err_y_low / model_val, err_y_high / model_val);
      point_idx++;
    }
  }

  ratioFrame->SetTitle("");
  ratioFrame->SetTitleSize(0);
  ratioFrame->GetYaxis()->SetTitleOffset(0.3);
  ratioFrame->GetYaxis()->SetTitle("Data / Fit");
  ratioFrame->GetYaxis()->SetTitleSize(0.08);
  ratioFrame->GetYaxis()->SetLabelSize(0.08);
  ratioFrame->GetYaxis()->SetRangeUser(0, 3);
  ratioFrame->GetYaxis()->CenterTitle();

  ratioFrame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  ratioFrame->GetXaxis()->SetTitleOffset(1.15);
  ratioFrame->GetXaxis()->SetLabelOffset(0.04);
  ratioFrame->GetXaxis()->SetLabelSize(0.08);
  ratioFrame->GetXaxis()->SetTitleSize(0.08);
  ratioFrame->GetXaxis()->CenterTitle();

  ratioFrame->GetYaxis()->SetTickSize(0.04);
  ratioFrame->GetYaxis()->SetNdivisions(404);
  ratioFrame->GetXaxis()->SetTickSize(0.03);
  ratioFrame->Draw();

  TLine *line_at_1 = new TLine(ctauLow, 1.0, ctauHigh, 1.0);
  line_at_1->SetLineColor(kRed);
  line_at_1->SetLineStyle(kDashed);
  line_at_1->Draw("same");

  g_ratio->SetMarkerStyle(kFullCircle);
  g_ratio->SetMarkerSize(0.7);
  g_ratio->Draw("P SAME");

  myCanvas->Update();
  myCanvas->Draw();
  myCanvas->SaveAs(Form("figs/2DFit_%s/CtauBkg/Bkg_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_ratio.png", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));

  delete h_tmp;
  delete myCanvas;
}

void CtauBkgFit::saveResults()
{
  std::cout << "===== saveResults() =====\n\n";
  TFile *outputFile = new TFile(Form("roots/2DFit_%s/CtauBkg/CtauBkgResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP), "recreate");

  fitResult->Write();
  ws->pdf("pdfCTAU_Bkg_Tot")->Write("pdfCTAU_Bkg_Tot");

  outputFile->Close();
  delete outputFile;
}