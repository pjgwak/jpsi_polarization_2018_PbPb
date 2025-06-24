#include <iostream>
#include <string>
#include <RooGaussian.h>

using namespace std;
using namespace RooFit;

// global variables - will be member variables
// got from JpsiUtility.h
const int nCtauResBins = 72;
const float ctauResLow = -10, ctauResHigh = 10;


void setupWorkspace(RooWorkspace *ws, const std::string &mass_path, const std::string &err_path);
double processDataset(RooWorkspace *ws);
void buildResModel(RooWorkspace *ws, int nGauss);
RooFitResult *perfromResFit(RooWorkspace *ws);
void drawResPlot(RooWorkspace *ws, double ctauResMin, double ctauResMax, int nGauss, const std::string &outpath); // Todo: use outPath
void saveOutput(RooWorkspace *ws, const std::string &outpath);

void ref_res()
{
  gStyle->SetEndErrorSize(0);

  // setupWorkspace
  RooWorkspace *ws_res = new RooWorkspace("ws_res");
  setupWorkspace(ws_res, "roots/mass.root", "roots/splot.root");
  
  // processDataset
  double fitRangeMin = processDataset(ws_res);
  double fitRangeMax = 0;

  // buildResModel
  int nGauss = 3;
  buildResModel(ws_res, nGauss);

  // perfromResFit
  RooFitResult *fit_result = perfromResFit(ws_res);
  

  // drawResPlot
  drawResPlot(ws_res, fitRangeMin, fitRangeMax, nGauss, "figs/res.png");

  // saveOutput
  saveOutput(ws_res, "roots/res.root");
  // ws_res->Print();

  delete fit_result;
  delete ws_res;
  std::cout << "===== Finish ref_res.C =====\n";
}

// RooWorkspace* setupWorkspace(TFile *f_mass, TFile *f_err)
void setupWorkspace(RooWorkspace *ws, const std::string &mass_path, const std::string &err_path)
{
  std::cout << "===== start setupWorkspace() =====\n";
  // get 2 input files
  TFile *f_mass = new TFile(mass_path.c_str());
  TFile *f_err = new TFile(err_path.c_str());

  if (!f_mass || !f_err) {
    std::cout << "Error: Could not open one of the files.\n";
    return;
  }

  // import nSig from mass
  RooWorkspace *ws_mass = (RooWorkspace*)f_mass->Get("wsMy");
  RooRealVar *nSig = (RooRealVar *)ws_mass->var("nSig");

  // import splot workspace
  RooWorkspace *ws_err = (RooWorkspace *)f_err->Get("wsMy");
  RooDataSet *ds_errSig = (RooDataSet *)ws_err->data("ds_errSig");

  ws->import(*nSig);
  ws->import(*ds_errSig);

  // cout << "\nnSig: " << ws->var("nSig")->getVal() << "+/-" << ws->var("nSig")->getError() << "\n\n";

  delete f_mass;
  delete f_err;
  std::cout << "===== finish setupWorkspace() =====\n";
}

double processDataset(RooWorkspace *ws)
{
  std::cout << "===== start processDataset() =====\n";
  // ===== make an elementary dataset =====
  // pick 3 variables up and apply a ctauErr range cut. Will be used till the end of Res fit
  // in splot "ctauResCutDS" was used for the dataset applied ctau3DErr range cut, same?
  double ctauErrMin = ws->var("ctau3DErr")->getMin();
  double ctauErrMax = ws->var("ctau3DErr")->getMax();

  // also apply Res < 0 cut. Nominal only use negative values
  RooDataSet *ctauResCutDS = (RooDataSet *)ws->data("ds_errSig")->reduce(RooArgSet(*(ws->var("ctau3DRes")), *(ws->var("ctau3D")), *(ws->var("ctau3DErr"))), Form("ctau3DRes<0&&ctau3DErr>%f&&ctau3DErr<%f", ctauErrMin, ctauErrMax));
  ctauResCutDS->SetName("ctauResCutDS");
  ws->import(*ctauResCutDS);

  // ===== decide ctauResMin range =====
  // and nominal only use negative values
  double ctauResMin = -10;
  double ctauResMax = 0;
  
  // find a proper ctauResMin
  TH1D *h_tmp = (TH1D *)ws->data("ds_errSig")->createHistogram(("h_tmp"), *ws->var("ctau3DRes"), Binning(nCtauResBins, ctauResLow, ctauResHigh));
  // double ctauResMax = h_tmp->GetBinCenter(h_tmp->FindLastBinAbove(1,1));
  for (int i = 0; i < h_tmp->GetNbinsX() / 2; i++)
  {
    // cout<<"Content: "<<h_tmp->GetBinContent(i)<<endl;
    if (h_tmp->GetBinContent(i) <= 0 && h_tmp->GetBinContent(i + 1) <= 0)
    {
      // cout<<"#####"<<i<<": "<<h_tmp->GetBinLowEdge(i+2)<<endl;
      ctauResMin = h_tmp->GetBinLowEdge(i + 2);
    }
    // if(h_tmp->GetBinContent(i)>1)ctauResMax = h_tmp->GetBinCenter(i)+h_tmp->GetBinWidth(i);
  }
  delete h_tmp;
  
  ws->var("ctau3DRes")->setRange("ctauResWindow", ctauResMin, 0);
  std::cout << "ctauRes Range: [" << ctauResMin << ", 0]\n";

  // === make a dataset with ctauRes range
  RooDataSet *ds_res = (RooDataSet *)ctauResCutDS->reduce(Form("ctau3DRes>=%.f&&ctau3DRes<=0", ctauResMin))->Clone("ds_res"); // original name: dataToFit
  ws->import(*ds_res);
  std::cout << "===== finish processDataset() =====\n";

  return ctauResMin;
}

void buildResModel(RooWorkspace *ws, int nGauss)
{
  std::cout << "===== start buildResModel() =====\n";

  // Todo: use config to set up the variable
  ws->factory("ctauRes_mean[0.]"); // fixed moean of gauss

  // sigmas
  ws->factory("s1_CtauRes[0.5, 1e-6, 1.0]");
  ws->factory("rS21_CtauRes[1.5, 1.0, 5.0]");  // s2 = s1 * rS21. To make sigma2 get bigger than sigma1
  ws->factory("rS32_CtauRes[2.5, 1.0, 5.0]");  // s3 = s2 * rS32
  ws->factory("rS43_CtauRes[1.5, 1.0, 10.0]"); // s4 = s3 * rS43

  // fractions
  ws->factory("f_CtauRes[0.2, 0., 1.]");
  ws->factory("f2_CtauRes[0.2, 0., 1.]");
  ws->factory("f3_CtauRes[0.5, 0., 1.]");

  // sigmas formulars
  ws->factory("RooFormulaVar::s2_CtauRes('@0*@1',{s1_CtauRes, rS21_CtauRes})");
  ws->factory("RooFormulaVar::s3_CtauRes('@0*@1',{s2_CtauRes, rS32_CtauRes})");
  ws->factory("RooFormulaVar::s4_CtauRes('@0*@1',{s3_CtauRes, rS43_CtauRes})");

  // create gauss
  ws->factory("RooGaussian::GaussModel1_ctauRes(ctau3DRes, ctauRes_mean, s1_CtauRes)");
  ws->factory("RooGaussian::GaussModel2_ctauRes(ctau3DRes, ctauRes_mean, s2_CtauRes)");
  ws->factory("RooGaussian::GaussModel3_ctauRes(ctau3DRes, ctauRes_mean, s3_CtauRes)");
  ws->factory("RooGaussian::GaussModel4_ctauRes(ctau3DRes, ctauRes_mean, s4_CtauRes)");

  // combine gauss - use the same names in old codes to avoid problems
  std::string finalModelName = "GaussModelCOND_ctauRes";
  if (nGauss == 2)
  {
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", finalModelName.c_str(), "GaussModel1_ctauRes", "GaussModel2_ctauRes", "f_CtauRes"));
  }
  else if (nGauss == 3)
  {
    ws->factory(Form("AddModel::GaussModel23_ctauRes({%s, %s}, {%s})", "GaussModel3_ctauRes", "GaussModel2_ctauRes", "f2_CtauRes"));
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", finalModelName.c_str(), "GaussModel1_ctauRes", "GaussModel23_ctauRes", "f_CtauRes"));
  }
  else if (nGauss == 4)
  {
    ws->factory(Form("AddModel::GaussModel43_ctauRes({%s, %s}, {%s})", "GaussModel4_ctauRes", "GaussModel3_ctauRes", "f3_CtauRes"));
    ws->factory(Form("AddModel::GaussModel23_ctauRes({%s, %s}, {%s})", "GaussModel43_ctauRes", "GaussModel2_ctauRes", "f2_CtauRes"));
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", finalModelName.c_str(), "GaussModel1_ctauRes", "GaussModel23_ctauRes", "f_CtauRes"));
  }
  else
  {
    std::cerr << "ERROR: nGauss should be 2, 3, or 4.\n";
    return;
  }

  ws->factory(Form("SUM::gaus_res_model(nSig*%s)", finalModelName.c_str()));

  // old method making fit model.
  // RooAddPdf *gaus_res_model = new RooAddPdf("gaus_res_model");
  // RooAbsPdf *ctauResModel = ctauResModel = new RooAddPdf("gaus_res_model", "gaus_res_model", *ws->pdf("GaussModelCOND_ctauRes"), *ws->var("nSig"));
  // ws->import(*ctauResModel);

  std::cout << "===== finish buildResModel() =====\n";
}

RooFitResult *perfromResFit(RooWorkspace *ws)
{
  std::cout << "===== start perfromResFit() =====\n";
  
  // fitTo
  bool isWeighted = ws->data("ds_res")->isWeighted();
  RooFitResult *fit_res = ws->pdf("gaus_res_model")->fitTo(
    *ws->data("ds_res"),
    Save(),
    SumW2Error(isWeighted),
    Extended(kTRUE),
    NumCPU(8),
    PrintLevel(-1));
  fit_res->Print("V");
  ws->import(*fit_res, "fit_res");

  // fix parameters for next steps
  std::cout << "fixing floating parameters to constant (after fitting)\n";
  const RooArgList &floatted_params = fit_res->floatParsFinal();
  for (auto arg : floatted_params) {
    RooRealVar *param = dynamic_cast<RooRealVar *>(arg);
    if (param)
    {
      std::cout << "Fixing parameter: " << param->GetName() << "\n";
      ws->var(param->GetName())->setConstant(kTRUE);
    }
  }

  return fit_res;
  std::cout << "===== finish perfromResFit() =====\n";
}

void drawResPlot(RooWorkspace *ws, double ctauResMin, double ctauResMax, int nGauss, const std::string &outpath)
{
  std::cout << "===== start drawResPlot() =====\n";
  // use default plotting setting.
  // if some kinematic bin needs to be redrawn -> use python plotting code

  // make a main canvas
  TCanvas *c_res = new TCanvas("c_res", "", 800, 600);
  c_res->cd();

  // make a pad1 - distribution
  TPad *pad1 = new TPad("pad1", "", 0, 0.25, 0.98, 1.0);
  pad1->SetTicks(1, 1);
  pad1->Draw();
  pad1->cd();

  // make a frame1 - distribution
  RooPlot *frame1 = ws->var("ctau3DRes")->frame(Bins(nCtauResBins), Range(ctauResLow, ctauResHigh)); // bins
  frame1->SetTitle("");
  pad1->cd();
  gPad->SetLogy();

  // ===== draw - dataset and models =====
  // dataset used for fitting
  ws->data("ds_res")->plotOn(frame1, Name("dataHist_ctauRes"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7), LineColor(kBlack), MarkerColor(kBlack));

  ws->pdf("gaus_res_model")->plotOn(frame1, Name("modelHist_ctauRes"), Precision(1e-4), NormRange("ctauResWindow"), Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), LineColor(kBlack));
  ws->pdf("gaus_res_model")->plotOn(frame1, Name("modelHist_gm1"), Precision(1e-4), NormRange("ctauResWindow"), Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel1_ctauRes")), LineColor(kGreen + 2));
  ws->pdf("gaus_res_model")->plotOn(frame1, Name("modelHist_gm2"), Precision(1e-4), NormRange("ctauResWindow"), Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel2_ctauRes")), LineColor(kRed + 2));

  if (nGauss >= 3)
  {
  ws->pdf("gaus_res_model")->plotOn(frame1, Name("modelHist_gm3"), Precision(1e-4), NormRange("ctauResWindow"), Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel3_ctauRes")), LineColor(kBlue + 2));
  }
  if (nGauss >= 4)
  {
    ws->pdf("gaus_res_model")->plotOn(frame1, Name("modelHist_gm4"), Precision(1e-4), NormRange("ctauResWindow"), Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel4_ctauRes")), LineColor(kMagenta + 2));
  }

  // To show the orignal data points - see the plot. red points exist beyond vertical line.
  ws->data("ctauResCutDS")->plotOn(frame1, Name("dataHist_ctauRes"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7), LineColor(kRed + 2), MarkerColor(kRed + 2));

  // ===== set good Y range =====
  TH1 *h = ws->data("ctauResCutDS")->createHistogram("hist", *ws->var("ctau3DRes"), Binning(frame1->GetNbinsX(), frame1->GetXaxis()->GetXmin(), frame1->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= h->GetNbinsX(); i++)
    if (h->GetBinContent(i) > 0)
      YMin = min(YMin, h->GetBinContent(i));
  Double_t Yup(0.), Ydown(0.);
  Yup = YMax * TMath::Power((YMax / YMin), (0.5 / (1.0 - 0.5 - 0.2)));
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.2 / (1.0 - 0.5 - 0.2))));
  frame1->GetYaxis()->SetRangeUser(Ydown, Yup);

  // ===== check lost event due to new Res range cut =====
  Double_t outTot = ws->data("ctauResCutDS")->numEntries();
  Double_t outRes = ws->data("ds_res")->numEntries();
  std::cout << "Tot evt: " << outTot << "\n";
  std::cout << "Fitted evt: " << outRes << "\n";
  std::cout << "lost evt: " << (outTot - outRes) << "(" << (outTot - outRes) * 100 / outTot << ")%" << "\n\n";

  // draw vertical lines to show a fit range
  if (outRes > 0.0)
  {
    // TLine   *minline = new TLine(rangeErr[0], 0.0, rangeErr[0], Ydown*TMath::Power((Yup/Ydown),0.4));
    TLine *minline = new TLine(ctauResMin, 0.0, ctauResMin, Ydown * TMath::Power((Yup / Ydown), 0.4));
    minline->SetLineStyle(2);
    minline->SetLineColor(1);
    minline->SetLineWidth(3);
    frame1->addObject(minline);
    TLine *maxline = new TLine(ctauResMax, 0.0, ctauResMax, Ydown * TMath::Power((Yup / Ydown), 0.4));
    maxline->SetLineStyle(2);
    maxline->SetLineColor(1);
    maxline->SetLineWidth(3);
    frame1->addObject(maxline);
  }

  // frame1 cosmetics
  frame1->GetXaxis()->CenterTitle();
  frame1->GetXaxis()->SetTitle("#frac{#font[12]{l}_{J/#psi}}{#sigma_{#font[12]{l}_{J/#psi}}}");
  frame1->SetFillStyle(4000);
  frame1->GetYaxis()->SetTitleOffset(2);
  frame1->GetXaxis()->SetLabelSize(0);
  frame1->GetXaxis()->SetTitleSize(0);
  frame1->Draw();

  // // ===== legend ===== //
  // skip now -> need to text_x, text_y etc
  // TLegend *leg_C = new TLegend(text_x + 0.29, text_y + 0.03, text_x + 0.39, text_y - 0.17);
  // leg_C->SetTextSize(text_size);
  // leg_C->SetTextFont(43);
  // leg_C->SetBorderSize(0);
  // leg_C->AddEntry(frame1->findObject("dataHist_ctauRes"), "Data", "pe");
  // leg_C->AddEntry(frame1->findObject("modelHist_ctauRes"), "Total PDF", "l");
  // leg_C->AddEntry(frame1->findObject("modelHist_gm1"), "Gauss 1", "l");
  // leg_C->AddEntry(frame1->findObject("modelHist_gm2"), "Gauss 2", "l");
  // if (nGauss >= 3)
  //   leg_C->AddEntry(frame1->findObject("modelHist_gm3"), "Gauss 3", "l");
  // if (nGauss == 4)
  //   leg_C->AddEntry(frame1->findObject("modelHist_gm4"), "Gauss 4", "l");
  // leg_C->Draw("same");
  // // cout<<"s2/s1: "<<ws->var("s2_CtauRes")->getVal()/ws->var("s1_CtauRes")->getVal()<<endl;
  // // cout<<"s3/s2: "<<ws->var("s3_CtauRes")->getVal()/ws->var("s2_CtauRes")->getVal()<<endl;
  // // cout << "s1: " << ws->var("s1_CtauRes")->getVal() << endl;

  // // ===== draw latex ===== //
  // drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x, text_y, text_color, text_size);
  // if (yLow == 0)
  //   drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x, text_y - y_diff, text_color, text_size);
  // else if (yLow != 0)
  //   drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x, text_y - y_diff, text_color, text_size);
  // drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x, text_y - y_diff * 2, text_color, text_size);
  // drawText(Form("Loss: (%.4f%s) %.f evts", (outTot - outRes) * 100 / outTot, "%", outTot - outRes), text_x, text_y - y_diff * 3, text_color, text_size);
  // // cout<<"lost evt: ("<<(outRes*100)/outTot<<")%, "<<outRes<<"evts"<<endl;
  // drawText(Form("N_{J/#psi} = %.f #pm %.f", ws->var("nSig")->getVal(), ws->var("nSig")->getError()), text_x + 0.5, text_y, text_color, text_size);
  // drawText(Form("s1_{Res} = %.4f #pm %.3f", ws->var("s1_CtauRes")->getVal(), ws->var("s1_CtauRes")->getError()), text_x + 0.5, text_y - y_diff * 1, text_color, text_size);
  // drawText(Form("(s2/s1)_{Res} = %.3f #pm %.3f", ws->var("rS21_CtauRes")->getVal(), ws->var("rS21_CtauRes")->getError()), text_x + 0.5, text_y - y_diff * 2, text_color, text_size);
  // drawText(Form("f_{Res} = %.4f #pm %.4f", ws->var("f_CtauRes")->getVal(), ws->var("f_CtauRes")->getError()), text_x + 0.5, text_y - y_diff * 3, text_color, text_size);
  // if (nGauss >= 3)
  //   drawText(Form("(s3/s2)_{Res} = %.4f #pm %.3f", ws->var("rS32_CtauRes")->getVal(), ws->var("rS32_CtauRes")->getError()), text_x + 0.5, text_y - y_diff * 4, text_color, text_size);

  // if (nGauss >= 3)
  //   drawText(Form("f2_{Res} = %.4f #pm %.3f", ws->var("f2_CtauRes")->getVal(), ws->var("f2_CtauRes")->getError()), text_x + 0.5, text_y - y_diff * 5, text_color, text_size);


  // ===== draw pull =====
  TPad *pullPad = new TPad("pullPad", "", 0, 0.001, 0.98, 0.32);
  c_res->cd();
  pullPad->Draw();
  pullPad->cd();
  pullPad->SetTopMargin(0); // Upper and lower plot are joined
  pullPad->SetBottomMargin(0.67);
  pullPad->SetBottomMargin(0.4);
  pullPad->SetFillStyle(4000);
  pullPad->SetFrameFillStyle(4000);
  pullPad->SetTicks(1, 1);

  RooPlot *pull_tmp = (RooPlot *)frame1->Clone("TMP");
  RooHist *pull_hist = pull_tmp->pullHist("dataHist_ctauRes", "modelHist_ctauRes", true);
  pull_hist->SetMarkerSize(0.8);
  RooPlot *pullFrame = ws->var("ctau3DRes")->frame(Title("Pull Distribution"), Bins(nCtauResBins), Range(ctauResLow, ctauResHigh));
  pullFrame->addPlotable(pull_hist, "PX");

  // cosmetics
  pullFrame->SetTitle("");
  pullFrame->SetTitleSize(0);
  pullFrame->GetYaxis()->SetTitleOffset(0.3);
  pullFrame->GetYaxis()->SetTitle("Pull");
  pullFrame->GetYaxis()->SetTitleSize(0.08);
  pullFrame->GetYaxis()->SetLabelSize(0.08);
  pullFrame->GetYaxis()->SetRangeUser(-9, 9);
  pullFrame->GetYaxis()->CenterTitle();

  pullFrame->GetXaxis()->SetTitle("#frac{#font[12]{l}_{J/#psi}}{#sigma_{#font[12]{l}}}");
  pullFrame->GetXaxis()->SetTitleOffset(1.55);
  pullFrame->GetXaxis()->SetLabelOffset(0.04);
  pullFrame->GetXaxis()->SetLabelSize(0.08);
  pullFrame->GetXaxis()->SetTitleSize(0.08);
  pullFrame->GetXaxis()->CenterTitle();

  pullFrame->GetYaxis()->SetTickSize(0.04);
  pullFrame->GetYaxis()->SetNdivisions(404);
  pullFrame->GetXaxis()->SetTickSize(0.03);
  pullFrame->Draw();

  // draw center line
  TLine *lC = new TLine(ctauResLow, 0, ctauResHigh, 0);
  lC->SetLineStyle(1);
  lC->Draw("same");

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

  // ===== draw canvas ===== //
  c_res->Update();
  c_res->Draw();
  
  // gSystem->mkdir(Form("figs/res/"), kTRUE);
  c_res->SaveAs(outpath.c_str());

  ws->import(*c_res);
  
  delete c_res;
  std::cout << "===== finish drawResPlot() =====\n";
}

void saveOutput(RooWorkspace *ws, const std::string &outpath)
{
  std::cout << "===== start saveOutput() =====\n";
  // gSystem->mkdir(Form("roots/res"), kTRUE);
  auto outfile = new TFile(outpath.c_str(), "recreate");
  // ctauResModel->Write(); // will be saved as gaus_res_model in root file
  ws->Write();
  outfile->Close();
  
  delete outfile;

  std::cout << "output saved to: " << outpath << "\n";
  std::cout << "===== finish saveOutput() =====\n";
}