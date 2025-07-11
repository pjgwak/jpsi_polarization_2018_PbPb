#include "MassFit.h"
#include "../headers/polarizationUtilities.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <filesystem>

#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMath.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1D.h"

#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooHist.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

using namespace RooFit;

MassFit::MassFit(const std::string &global_config_path, const std::string &local_config_path)
{
  RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling);
  loadConfiguration(global_config_path, local_config_path);
}

MassFit::~MassFit()
{
  RooMsgService::instance().getStream(1).addTopic(RooFit::ObjectHandling);
  delete ws;
  delete fitResult;
  if (fInputData)
  {
    if (fInputData->IsOpen())
      fInputData->Close();
    delete fInputData;
  }
  if (fMcParams)
  {
    if (fMcParams->IsOpen())
      fMcParams->Close();
    delete fMcParams;
  }
}

void MassFit::run()
{
  std::cout << "===== MassFit::run() =====" << "\n";
  setupWorkspaceAndData();
  setVariableRanges();
  initializeParameters();
  defineModel();
  performFit();
  makePlot();
  saveResults();
  std::cout << "===== MassFit finished =====" << "\n";
}

void MassFit::loadConfiguration(const std::string &global_path, const std::string &local_path)
{
  // read global and local config files
  // fill the member variables

  std::ifstream f_global(global_path);
  if (!f_global.is_open())
  {
    std::cerr << "FATAL: Global config file could not be opened: " << global_path << "\n";
    exit(1);
  }

  std::ifstream f_local(local_path);
  if (!f_local.is_open())
  {
    std::cerr << "FATAL: Local config file could not be opened: " << local_path << "\n";
    exit(1);
  }

  global_config = json::parse(f_global);
  local_config = json::parse(f_local);
  massfit_config = local_config.at("MassFit");

  // --- labeling ---
  std::string roodata_dir = global_config.at("roodata_dir");
  full_input_path = roodata_dir + "/" + massfit_config.at("input_file").get<std::string>();

  std::string date_tag = global_config.value("date_tag", "test");
  std::string base_output_roots = global_config.value("base_output_dir_roots", "roots/");
  std::string base_output_figs = global_config.value("base_output_dir_figs", "fig/");

  // --- set variables ---
  json kine = local_config.at("kinematics");
  ptLow = kine["ptLow"];
  ptHigh = kine["ptHigh"];
  yLow = kine["yLow"];
  yHigh = kine["yHigh"];
  cLow = kine["cLow"];
  cHigh = kine["cHigh"];
  cosLow = kine["cosLow"];
  cosHigh = kine["cosHigh"];

  json weights = local_config.at("weights");
  PR = weights["PR"];
  PRw = weights["PRw"];
  fEffW = weights["fEffW"];
  fAccW = weights["fAccW"];
  isPtW = weights["isPtW"];
  isTnP = weights["isTnP"];
  fname = PRw == 1 ? "PR" : "NP";

  // make output file path
  TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0, cLow, cHigh);
  std::string out_dir = Form("2DFit_%s/Mass", date_tag.c_str());

  output_root_path = base_output_roots + out_dir;
  output_fig_path = base_output_figs + out_dir;
  output_base_name = Form("Mass_%s", kineLabel.Data());

  gSystem->mkdir(output_root_path.c_str(), kTRUE);
  gSystem->mkdir(output_fig_path.c_str(), kTRUE);

  // --- get MC mass fit result ---
  std::string mc_mass_result_dir = base_output_roots + Form("2DFit_%s/mc_Mass", date_tag.c_str());
  std::string mc_mass_filename = Form("mc_Mass_%s_FitResult.root", kineLabel.Data());
  std::string full_mc_param_path = mc_mass_result_dir + "/" + mc_mass_filename;
  fMcParams = new TFile(full_mc_param_path.c_str());
  if (!fMcParams || fMcParams->IsZombie())
  {
    std::cerr << "FATAL: MassFit could not open required input from McMassFit: " << full_mc_param_path << "\n";
    exit(1);
  }

  // --- other setting ---
  plotting_style = massfit_config.value("plotting", json::object());
  nMassBin = massfit_config.at("plot_bins");
  massPlotMin = massfit_config["plot_range"]["min"];
  massPlotMax = massfit_config["plot_range"]["max"];
  massFitMin = massfit_config["fit_range"]["min"];
  massFitMax = massfit_config["fit_range"]["max"];

  ws = new RooWorkspace("workspace", "Mass Fit Workspace");
}

void MassFit::setupWorkspaceAndData()
{
  std::cout << "===== MassFit::setupWorkspaceAndData() =====" << "\n";

  // load the roodataset
  fInputData = new TFile(full_input_path.c_str());
  if (!fInputData || fInputData->IsZombie())
  {
    std::cerr << "FATAL: Input data file could not be opened: " << full_input_path << "\n";
    exit(1);
  }
  RooDataSet *dataset = (RooDataSet *)fInputData->Get("dataset");
  if (!dataset)
  {
    std::cerr << "FATAL: RooDataSet 'dataset' not found in input file: " << full_input_path << "\n";
    exit(1);
  }
  ws->import(*dataset, Rename("original_dataset"));

  // --- cuts ---
  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && cBin>=%d && cBin<%d", ptLow, ptHigh, yLow, yHigh, cLow, cHigh);
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )";
  TString osCut = "recoQQsign==0";
  TString angleCut = Form("cos_ep>%.2f && cos_ep<%.2f", cosLow, cosHigh);
  TString finalCut = TString::Format("%s && %s && %s && %s", osCut.Data(), accCut.Data(), kineCut.Data(), angleCut.Data());

  // --- set weight if need ---
  bool apply_weights = massfit_config.value("apply_event_weights", false);
  RooDataSet *ds_to_reduce = (RooDataSet *)ws->data("original_dataset");
  if (!ds_to_reduce)
  {
    std::cerr << "FATAL: 'original_dataset' not found in RooWorkspace.\n";
    exit(1);
  }

  RooDataSet *ds_weighted_tmp = nullptr;
  if (apply_weights)
  {
    std::cout << "INFO: Applying event weights for MassFit." << "\n";

    ds_weighted_tmp = new RooDataSet("ds_weighted_tmp", "", ds_to_reduce, *ds_to_reduce->get(), finalCut.Data(), "weight");
    ds_to_reduce = ds_weighted_tmp;
  }

  RooDataSet *dsAB = (RooDataSet *)ds_to_reduce->reduce(finalCut.Data());
  if (!dsAB)
  {
    std::cerr << "FATAL: Final reduced dataset 'dsAB' could not be created.\n";
    exit(1);
  }
  ws->import(*dsAB, Rename("dsAB"));

  if (ds_weighted_tmp)
    delete ds_weighted_tmp;
  delete dsAB;
}

void MassFit::setVariableRanges()
{
  std::cout << "===== MassFit::setVariableRanges() =====" << "\n";

  RooRealVar *mass = ws->var("mass");
  if (!mass)
  {
    std::cerr << "FATAL: RooRealVar 'mass' not found in workspace.\n";
    exit(1);
  }

  mass->setRange(massPlotMin, massPlotMax);
  mass->setRange("massPlotRange", massPlotMin, massPlotMax);
  // mass->setRange("massFitRange", massFitMin, massFitMax); // both are same
}

void MassFit::initializeParameters()
{
  std::cout << "===== MassFit::initializeParameters() =====\n";
  // --- get RooArgSet from datasetMass ---
  RooDataSet *dataset_fit = (RooDataSet *)fMcParams->Get("datasetMass");
  if (!dataset_fit)
  {
    std::cerr << "FATAL: datasetMass not found in MC fit result.\n";
    fMcParams->ls();
    exit(1);
  }
  RooWorkspace *ws_fit = new RooWorkspace("workspace_fit");
  ws_fit->import(*dataset_fit); // import mc parameters

  // --- get parameters according to pdf ---
  std::string model_type = massfit_config.value("signal_model_type", "DoubleCB");

  if (model_type == "DoubleCB")
  {
    mc_alpha = ws_fit->var("alpha_1_A")->getVal();
    mc_n = ws_fit->var("n_1_A")->getVal();
    mc_sigma = ws_fit->var("sigma_1_A")->getVal();
    mc_x = ws_fit->var("x_A")->getVal();
    mc_f = ws_fit->var("f")->getVal();

    if (mc_f > 1.0)
      mc_f = 1.0;
    if (mc_f < 0.0)
      mc_f = 0.0;
  }
  else if (model_type == "SingleGauss")
  {
    mc_sigma = ws_fit->var("sigma")->getVal();
  }
  else if (model_type == "CBG")
  {
    mc_sigma = ws_fit->var("sigma_cb")->getVal();
    mc_alpha = ws_fit->var("alpha_cb")->getVal();
    mc_n = ws_fit->var("n_cb")->getVal();
    mc_f = ws_fit->var("frac_cb")->getVal();

    if (mc_f > 1.0)
      mc_f = 1.0;
    if (mc_f < 0.0)
      mc_f = 0.0;
  }
  else
  {
    std::cerr << "FATAL: Unknown signal_model_type: " << model_type << "\n";
    exit(1);
  }
}

void MassFit::defineModel()
{
  std::cout << "===== MassFit::defineModel() =====" << "\n";

  // build signal
  RooRealVar mean_jpsi("mean_jpsi", "m_{J/#psi}", 3.096, 3.0, 3.2);
  RooRealVar sigma_1_A("sigma_1_A", "#sigma_{1}", mc_sigma, 0.01, 0.1);
  RooRealVar x_A("x_A", "#sigma_{2}/#sigma_{1}", mc_x);
  RooRealVar alpha_1_A("alpha_1_A", "#alpha", mc_alpha);
  RooRealVar n_1_A("n_1_A", "n", mc_n);
  RooRealVar f("f", "f_{CB1}", mc_f);

  ws->import({mean_jpsi, sigma_1_A, x_A, alpha_1_A, n_1_A, f});

  // fix parameters in json
  std::vector<std::string> fixed_params = massfit_config.value("fixed_params_from_mc", std::vector<std::string>{});
  for (const auto &pname : fixed_params)
  {
    if (ws->var(pname.c_str()))
    {
      std::cout << "INFO: Fixing parameter: " << pname << "\n";
      ws->var(pname.c_str())->setConstant(kTRUE);
    }
  }

  // DoubleCB
  ws->factory("RooFormulaVar::sigma_2_A('@0*@1',{sigma_1_A, x_A})");
  ws->factory("RooFormulaVar::alpha_2_A('1.0*@0',{alpha_1_A})");
  ws->factory("RooFormulaVar::n_2_A('1.0*@0',{n_1_A})");

  RooCBShape cb1("cball_1_A", "CB1", *ws->var("mass"), *ws->var("mean_jpsi"),
                 *ws->var("sigma_1_A"), *ws->var("alpha_1_A"), *ws->var("n_1_A"));
  RooCBShape cb2("cball_2_A", "CB2", *ws->var("mass"), *ws->var("mean_jpsi"),
                 *ws->function("sigma_2_A"), *ws->function("alpha_2_A"), *ws->function("n_2_A"));

  RooAddPdf pdf_sig("pdfMASS_Jpsi", "Signal Shape", {cb1, cb2}, *ws->var("f"));
  ws->import(pdf_sig, RecycleConflictNodes());

  // --- build bkg model ---
  json bkg_cfg = massfit_config.at("background_model");
  std::string bkg_type = bkg_cfg.at("type");

  if (bkg_type == "Chebychev")
  {
    RooArgList cheby_params;
    for (auto &[name, p] : bkg_cfg.at("params").items())
    {
      RooRealVar *param = new RooRealVar(name.c_str(), p.at("title").get<std::string>().c_str(),
                                         p.at("val").get<double>(), p.at("min").get<double>(), p.at("max").get<double>());
      cheby_params.add(*param);
      ws->import(*param);
    }
    RooChebychev bkg_pdf("pdfMASS_bkg", "Background Shape", *ws->var("mass"), cheby_params);
    ws->import(bkg_pdf, RecycleConflictNodes());
  }
  else
  {
    std::cerr << "FATAL: Unknown background model: " << bkg_type << "\n";
    exit(1);
  }

  // ---  yields and final model ---
  for (auto &[name, p] : massfit_config.at("yields").items())
  {
    ws->factory(Form("%s[%.1f, %.1f, %.1f]", name.c_str(),
                     p.at("val").get<double>(), p.at("min").get<double>(), p.at("max").get<double>()));
    ws->var(name.c_str())->SetTitle(p.at("title").get<std::string>().c_str());
  }

  RooAddPdf total_pdf("pdfMASS_Tot", "Signal + Background",
                      {*ws->pdf("pdfMASS_Jpsi"), *ws->pdf("pdfMASS_bkg")},
                      {*ws->var("N_Jpsi"), *ws->var("N_Bkg")});
  ws->import(total_pdf, RecycleConflictNodes());
}

void MassFit::performFit()
{
  std::cout << "===== MassFit::performFit() =====" << "\n";

  RooAbsPdf *pdf = ws->pdf("pdfMASS_Tot");
  RooAbsData *data = ws->data("dsAB");

  bool isWeighted = data->isWeighted();

  fitResult = pdf->fitTo(*data,
                         Save(),
                         Extended(kTRUE),
                        //  Range("massPlotRange"),
                         NumCPU(nCPU),
                         SumW2Error(isWeighted),
                         Strategy(2),
                         PrintLevel(0));

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

void MassFit::makePlot()
{
  std::cout << "===== MassFit::makePlot() =====" << "\n";

  auto plot_cfg = massfit_config.at("plotting");
  std::vector<std::string> params_to_plot = massfit_config.value("plot_params", std::vector<std::string>{});

  // set canvas and pad
  int canvas_w = plot_cfg["canvas_size"].value("width", 800);
  int canvas_h = plot_cfg["canvas_size"].value("height", 800);
  double pad_split = plot_cfg.value("pad_split_ratio", 0.3);
  auto canvas = new TCanvas("mass_canvas", "Mass Fit", canvas_w, canvas_h);

  auto pad_main = new TPad("pad_main", "", 0, pad_split, 1, 1);
  auto pad_pull = new TPad("pad_pull", "", 0, 0, 1, pad_split);
  pad_main->SetBottomMargin(0.02);
  pad_pull->SetTopMargin(0);
  pad_pull->SetBottomMargin(0.35);
  pad_main->Draw();
  pad_pull->Draw();

  // --- main plot ---
  pad_main->cd();
  pad_main->SetLogy();

  RooRealVar *mass = ws->var("mass");
  RooAbsData *data = ws->data("dsAB");
  RooAbsPdf *pdf = ws->pdf("pdfMASS_Tot");

  RooPlot *frame = mass->frame(Title(""), Bins(nMassBin), Range("massPlotRange"));
  data->plotOn(frame, Name("dataOS"));
  pdf->plotOn(frame, Name("pdfMASS_Tot"), LineColor(kRed), Range("massPlotRange"), NormRange("massPlotRange"));

  // signals
  auto *sig_pdf = dynamic_cast<RooAddPdf *>(ws->pdf("pdfMASS_Jpsi"));
  std::vector<std::string> sig_comp_names;

  if (sig_pdf)
  {
    const RooArgList &comps = sig_pdf->pdfList();

    std::vector<int> sig_colors = plot_cfg["signal_components"].value("colors", std::vector<int>{kBlue, kGreen + 2});
    std::vector<int> sig_styles = plot_cfg["signal_components"].value("styles", std::vector<int>{2, 4});

    for (int i = 0; i < comps.getSize(); ++i)
    {
      RooAbsPdf *signal_comp = dynamic_cast<RooAbsPdf *>(&comps[i]);
      if (!signal_comp)
        continue;

      TString name = signal_comp->GetName();
      TString plot_name = Form("bkg + %s", name.Data());
      sig_comp_names.push_back(plot_name.Data());

      RooArgSet comp_set;
      comp_set.add(*ws->pdf("pdfMASS_bkg"));
      comp_set.add(*signal_comp);

      ws->pdf("pdfMASS_Tot")->plotOn(frame, Components(comp_set), Name(plot_name), LineColor(sig_colors[i % sig_colors.size()]), LineStyle(sig_styles[i % sig_styles.size()]), Range("massPlotRange"), NormRange("massPlotRange"));
    }
  }

  // bkg model
  auto bkg_cfg = plot_cfg.value("background_component", json::object());
  pdf->plotOn(frame,
              Components(*ws->pdf("pdfMASS_bkg")),
              Name("pdfMASS_bkg"),
              LineColor(bkg_cfg.value("color", 920)),
              LineStyle(bkg_cfg.value("style", 2)),
              Range("massPlotRange"), NormRange("massPlotRange"));

  // Y-range
  auto hist_tmp = data->createHistogram("htmp", *mass, Binning(nMassBin, massPlotMin, massPlotMax));
  double ymax = hist_tmp->GetMaximum();
  double ymin = 1e9;
  for (int i = 1; i <= hist_tmp->GetNbinsX(); ++i)
    if (hist_tmp->GetBinContent(i) > 0)
      ymin = std::min(ymin, hist_tmp->GetBinContent(i));
  double ydown = ymin / std::pow(ymax / ymin, 0.1 / 0.5);
  double yup = 2*ymax * std::pow(ymax / ymin, 0.4 / 0.5);
  frame->GetYaxis()->SetRangeUser(std::max(ydown, 0.1), yup);
  delete hist_tmp;

  frame->GetYaxis()->SetTitle(Form("Events / ( %.3f GeV/c^{2} )", (massPlotMax - massPlotMin) / nMassBin));
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetXaxis()->SetTitleSize(0);
  frame->GetYaxis()->SetTitleOffset(1.4);
  frame->Draw();

  // --- print parameters ---
  json kine_cfg = plot_cfg.value("kinematics_text", json::object());

  double x_kine = kine_cfg["x"];
  double y_kine = kine_cfg["y"];
  double dy_kine = kine_cfg["dy"];
  int size_kine = kine_cfg["size"];

  drawText(Form("%.1f < p_{T} < %.1f GeV/c", ptLow, ptHigh), x_kine, y_kine, kBlack, size_kine);
  drawText(Form("%.1f < |y| < %.1f", yLow, yHigh), x_kine, y_kine - dy_kine, kBlack, size_kine);
  drawText(Form("Cent. %d-%d %%", cLow / 2, cHigh / 2), x_kine, y_kine - 2 * dy_kine, kBlack, size_kine);

  json param_cfg = plot_cfg.value("parameters_text", json::object());
  float px = param_cfg.value("x", 0.6f);
  float py = param_cfg.value("y", 0.85f);
  float dy = param_cfg.value("dy", 0.055f);
  int size = param_cfg.value("size", 18);

  for (const auto &pname_json : params_to_plot)
  {
    std::string pname = pname_json;
    RooRealVar *param = ws->var(pname.c_str());
    if (!param)
      continue;

    TString title = param->GetTitle();
    double val = param->getVal();
    double err = param->getError();

    if (TString(pname).BeginsWith("N_"))
    {
      drawText(Form("%s = %.0f #pm %.0f", title.Data(), val, err), px, py, kBlack, size);
    }
    else
    {
      int precision = 3;
      double min_error = TMath::Power(10, -precision);
      if (param->isConstant())
      {
        drawText(Form("%s = %.*f (fixed)", title.Data(), precision, val), px, py, kBlack, size);
      }
      else if (err < min_error)
      {
        drawText(Form("%s = %.*f #pm < %.*f", title.Data(), precision, val, precision, min_error), px, py, kBlack, size);
      }
      else
      {
        drawText(Form("%s = %.*f #pm %.*f", title.Data(), precision, val, precision, err), px, py, kBlack, size);
      }
    }

    py -= dy;
  }

  // --- legend ---
  auto leg_cfg = plot_cfg["legend"];
  auto legend = new TLegend(leg_cfg["x1"], leg_cfg["y1"], leg_cfg["x2"], leg_cfg["y2"]);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(leg_cfg.value("text_size", 0.04));
  legend->AddEntry(frame->findObject("dataOS"), "Data", "pe");
  // legend->AddEntry((TObject *)nullptr, "", "");
  legend->AddEntry(frame->findObject("pdfMASS_Tot"), "Total Fit", "l");
  legend->AddEntry(frame->findObject("pdfMASS_bkg"), "Background", "l");

  // signal + bkg component
  for (const auto &plot_name : sig_comp_names)
  {
    TObject *obj = frame->findObject(plot_name.c_str());
    if (obj) {
      legend->AddEntry(obj, plot_name.c_str(), "l");
    }
      
  }

  legend->Draw("same");

  // --- Pull plot ---
  pad_pull->cd();
  RooPlot *frameTmp = (RooPlot *)frame->Clone("frameTMP");
  RooPlot *pullFrame = mass->frame(Title(" "), Bins(nMassBin), Range("massPlotRange"));
  RooHist *pullHist = frameTmp->pullHist("dataOS", "pdfMASS_Tot", true);
  pullFrame->addPlotable(pullHist, "P");

  double pullMin = plot_cfg["pull_plot"].value("y_range_min", -5.0);
  double pullMax = plot_cfg["pull_plot"].value("y_range_max", 5.0);
  pullFrame->GetYaxis()->SetRangeUser(pullMin, pullMax);
  pullFrame->GetYaxis()->SetTitle("Pull");
  pullFrame->GetYaxis()->SetTitleSize(0.12);
  pullFrame->GetYaxis()->SetLabelSize(0.12);
  pullFrame->GetYaxis()->SetTitleOffset(0.35);
  pullFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  pullFrame->GetXaxis()->SetTitleSize(0.12);
  pullFrame->GetXaxis()->SetLabelSize(0.12);
  pullFrame->GetXaxis()->SetTitleOffset(1.15);
  pullFrame->GetXaxis()->CenterTitle();
  pullFrame->Draw();

  TLine *l0 = new TLine(massLow, 0, massHigh, 0);
  l0->SetLineStyle(kDashed);
  l0->Draw("same");

  printChi2(ws, pad_pull, frameTmp, fitResult, "mass", "dataOS", "pdfMASS_Tot", nMassBin, false);

  std::string out_path = Form("%s/%s.png", output_fig_path.c_str(), output_base_name.c_str());
  canvas->SaveAs(out_path.c_str());
  std::cout << "Plot saved to: " << out_path << "\n";

  delete frameTmp;
  delete pullFrame;
  delete canvas;
}

void MassFit::saveResults()
{
  std::cout << "===== MassFit::saveResults() =====" << "\n";

  std::string root_save_path = Form("%s/%s_FitResult.root", output_root_path.c_str(), output_base_name.c_str());
  TFile *outputFile = new TFile(root_save_path.c_str(), "RECREATE");

  if (!outputFile || outputFile->IsZombie())
  {
    std::cerr << "ERROR: Could not create output file: " << root_save_path << '\n';
    return;
  }

  outputFile->cd();

  if (ws->pdf("pdfMASS_Tot"))
    ws->pdf("pdfMASS_Tot")->Write("pdfMASS_Tot");
  if (ws->pdf("pdfMASS_bkg"))
    ws->pdf("pdfMASS_bkg")->Write("pdfMASS_bkg");

  if (fitResult)
    fitResult->Write("fitResult");

  RooArgSet *fitargs = new RooArgSet();
  fitargs->add(fitResult->floatParsFinal());
  RooDataSet *datasetMass = new RooDataSet("datasetMass", "dataset with Mass Fit result", *fitargs);
  datasetMass->add(*fitargs);
  datasetMass->Write("datasetMass");

  outputFile->Close();
  std::cout << "Results saved to " << root_save_path << '\n';

  delete fitargs;
  delete datasetMass;
  delete outputFile;
}