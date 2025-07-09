#include "CtauResFit.h"
#include <iostream>
#include <algorithm>

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
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooGaussModel.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooFormulaVar.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooAbsReal.h"
#include "RooAbsArg.h"
#include "RooGaussModel.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "../headers/polarizationUtilities.h"

using namespace RooFit;

CtauResFit::CtauResFit(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
                       float cosLow, float cosHigh, int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP)
    : ptLow(ptLow), ptHigh(ptHigh), yLow(yLow), yHigh(yHigh), cLow(cLow), cHigh(cHigh),
      cosLow(cosLow), cosHigh(cosHigh), PR(PR), PRw(PRw),
      fEffW(fEffW), fAccW(fAccW), isPtW(isPtW), isTnP(isTnP)
{
  RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling);
  gStyle->SetEndErrorSize(0);
}

CtauResFit::~CtauResFit()
{
  RooMsgService::instance().getStream(1).addTopic(RooFit::ObjectHandling);
  if (fMass) { fMass->Close(); delete fMass; }
  if (fCErr) { fCErr->Close(); delete fCErr; }
  delete ws;
  delete fitResult;
}

void CtauResFit::run()
{
  std::cout << "===== run() =====\n\n";
  setLabels();
  openInputFile();
  setupWorkspaceAndData();
  setVariableRanges();
  defineModel();
  performFit();
  makePlot();
  saveResults();
}

void CtauResFit::setLabels()
{
  std::cout << "===== setLabels() =====\n\n";
  DATE = "No_Weight_2";
  gSystem->mkdir(Form("roots/2DFit_%s/CtauRes", DATE.c_str()), kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/CtauRes", DATE.c_str()), kTRUE);

  kineLabel = std::string(getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh).Data());

  if (PRw == 1)
    fname = "PR";
  else if (PRw == 2)
    fname = "NP";
}

void CtauResFit::openInputFile()
{
  std::cout << "===== openInputFile() =====\n\n";
  fMass = new TFile(Form("roots/2DFit_%s/Mass/Mass_FixedFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));
  fCErr = new TFile(Form("roots/2DFit_%s/CtauErr/CtauErrResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));

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
}

void CtauResFit::setupWorkspaceAndData()
{
  std::cout << "===== setupWorkspaceAndData() =====\n\n";
  RooAddPdf *pdfMASS_Tot = (RooAddPdf *)fMass->Get("pdfMASS_Tot");
  RooDataSet *dataw_Sig = (RooDataSet *)fCErr->Get("dataw_Sig");

  ws = new RooWorkspace("workspace");
  ws->import(*pdfMASS_Tot);
  ws->import(*dataw_Sig);

  double ctauErrDataMin, ctauErrDataMax;
  dataw_Sig->getRange(*ws->var("ctau3DErr"), ctauErrDataMin, ctauErrDataMax);

  TString cutTmp = Form("ctau3DRes<0 && ctau3DErr>%f && ctau3DErr<%f", ctauErrDataMin, ctauErrDataMax);
  RooDataSet *ctauResCutDS = (RooDataSet *)dataw_Sig->reduce(cutTmp.Data());
  ctauResCutDS->SetName("ctauResCutDS");
  ws->import(*ctauResCutDS);
  delete ctauResCutDS;
}

void CtauResFit::setVariableRanges()
{
  std::cout << "===== setVariableRanges() =====\n\n";

  TH1D *hTot = (TH1D *)ws->data("dataw_Sig")->createHistogram(("hTot"), *ws->var("ctau3DRes"), Binning(nCtauResBins, ctauResLow, ctauResHigh));

  if (ctauResMin == -100) {
    ctauResMin = hTot->GetBinLowEdge(hTot->FindFirstBinAbove(1, 1));
    if (ctauResMin < -10)
      ctauResMin = -10;
  }

  if (ctauResMax == 100) {
    ctauResMax = hTot->GetBinLowEdge(hTot->FindLastBinAbove(0, 1)) + hTot->GetBinWidth(hTot->FindLastBinAbove(0, 1));
    if (ctauResMax > -10)
      ctauResMax = 10;
  }
  delete hTot;

  ws->var("ctau3DRes")->setRange("resFitRange", ctauResMin, ctauResMax);
}

void CtauResFit::defineModel()
{
  std::cout << "===== defineModel() =====\n\n";
  // variables
  ws->factory("One[1.0]");
  ws->factory("ctauRes_mean[0.0]");

  ws->factory("ctau1_CtauRes[0.]");
  ws->factory("ctau2_CtauRes[0.]");
  ws->factory("ctau3_CtauRes[0.]");
  ws->factory("ctau4_CtauRes[0.]");

  ws->factory("s1_CtauRes[0.7, 0.1, 1.1]");
  ws->factory("rS21_CtauRes[1.5, 1.0, 3.0]");
  ws->factory("rS32_CtauRes[1.5, 1.0, 3.0]");
  ws->factory("rS43_CtauRes[1.5, 1.0, 3.0]");

  ws->factory("f_CtauRes[0.344, 0.01, 1.]");
  ws->factory("f2_CtauRes[0.322, 0.01, 1.]");
  ws->factory("f3_CtauRes[0.5, 1e-6, 1.]");

  ws->factory("RooFormulaVar::s2_CtauRes('@0*@1',{rS21_CtauRes,s1_CtauRes})");
  ws->factory("RooFormulaVar::s3_CtauRes('@0*@1',{rS32_CtauRes,s2_CtauRes})");
  ws->factory("RooFormulaVar::s4_CtauRes('@0*@1',{rS43_CtauRes,s3_CtauRes})");

  // build models
  TString varName = "ctau3DRes";
  ws->factory(Form("GaussModel::%s(%s, %s, %s,One,One)", "GaussModel1_ctauRes", varName.Data(),
                   "ctau1_CtauRes", //"ctau1_CtauRes",
                   "s1_CtauRes"));
  ws->factory(Form("GaussModel::%s(%s, %s, %s,One,One)", "GaussModel2_ctauRes", varName.Data(),
                   "ctau2_CtauRes", //"ctau2_CtauRes",
                   "s2_CtauRes"));
  ws->factory(Form("GaussModel::%s(%s, %s, %s,One,One)", "GaussModel3_ctauRes", varName.Data(),
                   "ctau3_CtauRes", //"ctau3_CtauRes",
                   "s3_CtauRes"));
  ws->factory(Form("GaussModel::%s(%s, %s, %s,One,One)", "GaussModel4_ctauRes", varName.Data(),
                   "ctau4_CtauRes", //"ctau3_CtauRes",
                   "s4_CtauRes"));

  // combine models
  if (nGauss == 4)
  {
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModel43_ctauRes", "GaussModel4_ctauRes", "GaussModel3_ctauRes", "f3_CtauRes"));
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModel23_ctauRes", "GaussModel43_ctauRes", "GaussModel2_ctauRes", "f2_CtauRes"));
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModelCOND_ctauRes", "GaussModel1_ctauRes", "GaussModel23_ctauRes", "f_CtauRes"));
  }
  else if (nGauss == 3)
  {
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModel23_ctauRes", "GaussModel3_ctauRes", "GaussModel2_ctauRes", "f2_CtauRes"));
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModelCOND_ctauRes", "GaussModel1_ctauRes", "GaussModel23_ctauRes", "f_CtauRes"));
  }
  else if (nGauss == 2)
  {
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModelCOND_ctauRes", "GaussModel1_ctauRes", "GaussModel2_ctauRes", "f_CtauRes"));
  }

  RooExtendPdf GaussModel_Tot("GaussModel_Tot", "Total Resolution Model",
                          *ws->pdf("GaussModelCOND_ctauRes"), *ws->var("N_Jpsi"));
  ws->import(GaussModel_Tot);

  // RooAddPdf *GaussModel_Tot = new RooAddPdf("GaussModel_Tot");
  // RooAbsPdf *ctauResModel = ctauResModel = new RooAddPdf("GaussModel_Tot", "GaussModel_Tot", *ws->pdf("GaussModelCOND_ctauRes"), *ws->var("N_Jpsi"));
  // ws->import(*ctauResModel);
}

void CtauResFit::performFit()
{
  std::cout << "===== performFit() =====\n\n";
  RooAbsData *initial_data = ws->data("ctauResCutDS");
  RooDataSet *dataToFit = (RooDataSet *)initial_data->reduce(CutRange("resFitRange"));
  ws->import(*dataToFit, Rename("dataToFit"));
  delete dataToFit;

  RooAbsPdf *model = ws->pdf("GaussModel_Tot");
  RooAbsData *finalData = ws->data("dataToFit"); // don't use the dataToFit local variable.

  bool isWeighted = finalData->isWeighted();
  fitResult = model->fitTo(*finalData,
                           Save(),
                           Extended(true),
                           NumCPU(nCPU),
                           PrintLevel(-1),
                           SumW2Error(isWeighted),
                           Strategy(2));

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
}

void CtauResFit::makePlot()
{
  std::cout << "===== makePlot() =====\n\n";
  TCanvas *myCanvas = new TCanvas("canvas_C", "", 800, 800);

  TPad *pad1 = new TPad("pad1", "", 0, 0.16, 0.98, 1.0);
  pad1->SetTicks(1, 1);
  pad1->Draw();
  pad1->cd();
  gPad->SetLogy();

  RooPlot *resFrame = ws->var("ctau3DRes")->frame(Bins(nCtauResBins), Range(ctauResLow, ctauResHigh)); // bins
  resFrame->SetTitle("");

  ws->data("ctauResCutDS")->plotOn(resFrame, Name("dataHist_ctauRes"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7), LineColor(kRed + 2), MarkerColor(kRed + 2));
  ws->pdf("GaussModel_Tot")->plotOn(resFrame, Name("modelHist_ctauRes"), NormRange("resFitRange"), LineColor(kBlack));
  ws->pdf("GaussModel_Tot")->plotOn(resFrame, Name("modelHist_gm1"), NormRange("resFitRange"), Components(*ws->pdf("GaussModel1_ctauRes")), LineColor(kGreen + 2));
  ws->pdf("GaussModel_Tot")->plotOn(resFrame, Name("modelHist_gm2"), NormRange("resFitRange"), Components(*ws->pdf("GaussModel2_ctauRes")), LineColor(kRed + 2));
  if (nGauss >=3)
    ws->pdf("GaussModel_Tot")->plotOn(resFrame, Name("modelHist_gm3"), NormRange("resFitRange"), Components(*ws->pdf("GaussModel3_ctauRes")), LineColor(kBlue + 2));
  if (nGauss >= 4)
    ws->pdf("GaussModel_Tot")->plotOn(resFrame, Name("modelHist_gm4"), NormRange("resFitRange"), Components(*ws->pdf("GaussModel4_ctauRes")), LineColor(kMagenta + 2));

  // y-axis range
  TH1 *hTmp = ws->data("ctauResCutDS")->createHistogram("hTmp", *ws->var("ctau3DRes"), Binning(resFrame->GetNbinsX(), resFrame->GetXaxis()->GetXmin(), resFrame->GetXaxis()->GetXmax()));
  double YMax = hTmp->GetBinContent(hTmp->GetMaximumBin());
  double YMin = 1e99;
  for (int i = 1; i <= hTmp->GetNbinsX(); i++)
    if (hTmp->GetBinContent(i) > 0)
      YMin = std::min(YMin, hTmp->GetBinContent(i));
  delete hTmp;
  double Yup(0.), Ydown(0.);
  Yup = YMax * TMath::Power((YMax / YMin), (0.5 / (1.0 - 0.5 - 0.2)));
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.2 / (1.0 - 0.5 - 0.2))));
  resFrame->GetYaxis()->SetRangeUser(Ydown, Yup);

  double outTot = ws->data("ctauResCutDS")->numEntries();
  double outRes = ws->data("dataToFit")->numEntries();
  std::cout << "Tot evt: " << outTot << "\n";
  std::cout << "Res evt: " << outRes << "\n";
  std::cout << "lost evt: " << (outTot - outRes) << " events (" << (outTot - outRes) * 100 / outTot << " %)\n";

  if (outRes > 0.0)
  {
    // TLine   *minline = new TLine(rangeErr[0], 0.0, rangeErr[0], Ydown*TMath::Power((Yup/Ydown),0.4));
    TLine *minline = new TLine(ctauResMin, 0.0, ctauResMin, Ydown * TMath::Power((Yup / Ydown), 0.4));
    minline->SetLineStyle(2);
    minline->SetLineColor(1);
    minline->SetLineWidth(3);
    resFrame->addObject(minline);
    TLine *maxline = new TLine(ctauResMax, 0.0, ctauResMax, Ydown * TMath::Power((Yup / Ydown), 0.4));
    maxline->SetLineStyle(2);
    maxline->SetLineColor(1);
    maxline->SetLineWidth(3);
    resFrame->addObject(maxline);
  }
  resFrame->GetXaxis()->CenterTitle();
  resFrame->GetXaxis()->SetTitle("#frac{#font[12]{l}_{J/#psi}}{#sigma_{#font[12]{l}_{J/#psi}}}");
  resFrame->SetFillStyle(4000);
  resFrame->GetYaxis()->SetTitleOffset(1.43);
  resFrame->GetXaxis()->SetLabelSize(0);
  resFrame->GetXaxis()->SetTitleSize(0);
  resFrame->Draw();

  TLegend *leg_C = new TLegend(text_x + 0.29, text_y + 0.03, text_x + 0.39, text_y - 0.17);
  leg_C->SetTextSize(text_size);
  leg_C->SetTextFont(43);
  leg_C->SetBorderSize(0);
  leg_C->AddEntry(resFrame->findObject("dataHist_ctauRes"), "Data", "pe");
  leg_C->AddEntry(resFrame->findObject("modelHist_ctauRes"), "Total PDF", "l");
  leg_C->AddEntry(resFrame->findObject("modelHist_gm1"), "Gauss 1", "l");
  leg_C->AddEntry(resFrame->findObject("modelHist_gm2"), "Gauss 2", "l");
  if (nGauss >= 3)
    leg_C->AddEntry(resFrame->findObject("modelHist_gm3"), "Gauss 3", "l");
  if (nGauss >= 4)
    leg_C->AddEntry(resFrame->findObject("modelHist_gm4"), "Gauss 4", "l");
  leg_C->Draw("same");

  // std::cout<<"s2/s1: "<<ws->var("s2_CtauRes")->getVal()/ws->var("s1_CtauRes")->getVal()<<"\n";
  // std::cout<<"s3/s2: "<<ws->var("s3_CtauRes")->getVal()/ws->var("s2_CtauRes")->getVal()<<"\n";
  // std::cout << "s1: " << ws->var("s1_CtauRes")->getVal() << "\n";

  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x, text_y, text_color, text_size);
  if (yLow == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x, text_y - y_diff, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x, text_y - y_diff, text_color, text_size);
  drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x, text_y - y_diff * 2, text_color, text_size);
  drawText(Form("Loss: %.f evts (%.4f %s)", outTot - outRes, (outTot - outRes) * 100 / outTot, "%"), text_x, text_y - y_diff * 3, text_color, text_size);
  // std::cout<<"lost evt: ("<<(outRes*100)/outTot<<")%, "<<outRes<<"evts"<<"\n";

  drawText(Form("N_{J/#psi} = %.f #pm %.f", ws->var("N_Jpsi")->getVal(), ws->var("N_Jpsi")->getError()), text_x + 0.5, text_y, text_color, text_size);
  drawText(Form("s1_{Res} = %.4f #pm %.4f", ws->var("s1_CtauRes")->getVal(), ws->var("s1_CtauRes")->getError()), text_x + 0.5, text_y - y_diff * 1, text_color, text_size);
  drawText(Form("(s2/s1)_{Res} = %.4f #pm %.4f", ws->var("rS21_CtauRes")->getVal(), ws->var("rS21_CtauRes")->getError()), text_x + 0.5, text_y - y_diff * 2, text_color, text_size);
  if (nGauss == 3)
  {
    drawText(Form("(s3/s2)_{Res} = %.4f #pm %.4f", ws->var("rS32_CtauRes")->getVal(), ws->var("rS32_CtauRes")->getError()), text_x + 0.5, text_y - y_diff * 3, text_color, text_size);
    drawText(Form("f_{Res} = %.4f #pm %.4f", ws->var("f_CtauRes")->getVal(), ws->var("f_CtauRes")->getError()), text_x + 0.5, text_y - y_diff * 4, text_color, text_size);
    drawText(Form("f2_{Res} = %.4f #pm %.4f", ws->var("f2_CtauRes")->getVal(), ws->var("f2_CtauRes")->getError()), text_x + 0.5, text_y - y_diff * 5, text_color, text_size);
  }
  else if (nGauss == 2)
  {
    drawText(Form("f_{Res} = %.4f #pm %.4f", ws->var("f_CtauRes")->getVal(), ws->var("f_CtauRes")->getError()), text_x + 0.5, text_y - y_diff * 3, text_color, text_size);
  }

  // pull pad
  TPad *pad2 = new TPad("pad2", "", 0, 0.006, 0.98, 0.227);
  myCanvas->cd();
  pad2->Draw();
  pad2->cd();
  pad2->SetTopMargin(0); // Upper and lower plot are joined
  pad2->SetBottomMargin(0.67);
  pad2->SetBottomMargin(0.4);
  pad2->SetFillStyle(4000);
  pad2->SetFrameFillStyle(4000);
  pad2->SetTicks(1, 1);

  // ignore empty bin warning of residual pull. - comment them out when you want to test the codes.
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);

  RooPlot *frameTMP = (RooPlot *)resFrame->Clone("frameTMP");
  RooHist *pullHist = frameTMP->pullHist("dataHist_ctauRes", "modelHist_ctauRes", true);
  pullHist->SetMarkerSize(0.8);
  // RooPlot* pullFrame = ws->var("ctau3DRes")->frame(Title("Pull Distribution"), Bins(nCtauResBins), Range(ctauResLow, ctauResHigh)) ;
  RooPlot *pullFrame = ws->var("ctau3DRes")->frame(Title("Pull Distribution"), Bins(resFrame->GetNbinsX()), Range(resFrame->GetXaxis()->GetXmin(), resFrame->GetXaxis()->GetXmax()));
  pullFrame->addPlotable(pullHist, "PX");
  pullFrame->SetTitle("");
  pullFrame->SetTitleSize(0);
  pullFrame->GetYaxis()->SetTitleOffset(0.3);
  pullFrame->GetYaxis()->SetTitle("Pull");
  pullFrame->GetYaxis()->SetTitleSize(0.15);
  pullFrame->GetYaxis()->SetLabelSize(0.15);
  pullFrame->GetYaxis()->SetRangeUser(-3.8, 3.8);
  pullFrame->GetYaxis()->CenterTitle();

  pullFrame->GetXaxis()->SetTitle("#frac{#font[12]{l}_{J/#psi}}{#sigma_{#font[12]{l}_{J/#psi}}}");
  pullFrame->GetXaxis()->SetTitleOffset(1.05);
  pullFrame->GetXaxis()->SetLabelOffset(0.04);
  pullFrame->GetXaxis()->SetLabelSize(0.15);
  pullFrame->GetXaxis()->SetTitleSize(0.15);
  pullFrame->GetXaxis()->CenterTitle();

  pullFrame->GetYaxis()->SetTickSize(0.04);
  pullFrame->GetYaxis()->SetNdivisions(404);
  pullFrame->GetXaxis()->SetTickSize(0.03);
  pullFrame->Draw();

  TLine *lC = new TLine(ctauResLow, 0, ctauResHigh, 0);
  lC->SetLineStyle(1);
  lC->Draw("same");

  printChi2(ws, pad2, frameTMP, fitResult, "ctau3DRes", "dataHist_ctauRes", "modelHist_ctauRes", nCtauResBins, false);
  pad2->Update();

  myCanvas->Update();
  myCanvas->SaveAs(Form("figs/2DFit_%s/CtauRes/CtauRes_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.png", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP));

  // turn on alarm
  RooMsgService::instance().getStream(0).addTopic(Plotting);
  RooMsgService::instance().getStream(1).addTopic(Plotting);

  delete frameTMP;
  delete resFrame;
  delete pullFrame;
  delete myCanvas;
}

void CtauResFit::saveResults()
{
  std::cout << "===== saveResults() =====\n\n";
  TFile *outputFile = new TFile(Form("roots/2DFit_%s/CtauRes/CtauResResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.c_str(), kineLabel.c_str(), fname.c_str(), fEffW, fAccW, isPtW, isTnP), "recreate");

  ws->pdf("GaussModel_Tot")->Write("GaussModel_Tot");
  fitResult->Write();

  outputFile->Close();
  delete outputFile;
}