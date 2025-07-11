#include "McMassFit.h"
#include <filesystem>
#include <fstream>
#include <iostream>

#include "TLatex.h"
#include "TFile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"

#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLine.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "TCanvas.h"
#include "TMath.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooMsgService.h"

#include "../headers/polarizationUtilities.h"

using namespace RooFit;

McMassFit::McMassFit(const std::string &global_config_path, const std::string &local_config_path)
{
  RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling);
  loadConfiguration(global_config_path, local_config_path);
}

McMassFit::~McMassFit()
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
}

void McMassFit::run()
{
  std::cout << "===== McMassFit::run() =====\n\n";
  setupWorkspaceAndData();
  setVariableRanges();
  defineModel();
  performFit();
  makePlot();
  saveResults();
  std::cout << "===== McMassFit finished =====\n";
}

void McMassFit::loadConfiguration(const std::string &global_path_str, const std::string &local_path_str)
{
  try
  {
    // load global config file
    std::ifstream f_global(global_path_str);
    if (!f_global.is_open())
    {
      std::cerr << "FATAL: Global config file could not be opened: " << global_path_str << "\n";
      exit(1);
    }
    global_config = json::parse(f_global);

    // load local config file
    std::ifstream f_local(local_path_str);
    if (!f_local.is_open())
    {
      std::cerr << "FATAL: Local config file could not be opened: " << local_path_str << "\n";
      exit(1);
    }
    local_config = json::parse(f_local);

    // --- set variables from json file ---
    std::string roodata_dir = global_config.at("roodata_dir");
    json mcMassConfig = local_config.at("McMassFit");
    full_input_path = roodata_dir + "/" + mcMassConfig.at("input_file").get<std::string>();

    // nCPU = global_config.value("num_cpu", 8);

    json kine_settings = local_config.at("kinematics");
    ptLow = kine_settings.at("ptLow");
    ptHigh = kine_settings.at("ptHigh");
    yLow = kine_settings.at("yLow");
    yHigh = kine_settings.at("yHigh");
    cLow = kine_settings.at("cLow");
    cHigh = kine_settings.at("cHigh");
    cosLow = kine_settings.at("cosLow");
    cosHigh = kine_settings.at("cosHigh");

    json flag_settings = local_config.at("weights");
    PRw = flag_settings.at("PRw");
    PR = flag_settings.at("PR");
    fEffW = flag_settings.at("fEffW");
    fAccW = flag_settings.at("fAccW");
    isPtW = flag_settings.at("isPtW");
    isTnP = flag_settings.at("isTnP");
    fname = (PRw == 1) ? "PR" : "NP";

    // model info
    yields = mcMassConfig.at("yields");
    signal_model_type = mcMassConfig.at("signal_model").at("type");
    signal_model_params = mcMassConfig.at("signal_model").at("params");

    // plotting info
    params_to_plot = mcMassConfig.value("plot_params", std::vector<std::string>{});
    plotting_style = mcMassConfig.value("plotting", json::object());
    nMassBin = mcMassConfig.at("plot_bins");
    massPlotMin = mcMassConfig.at("plot_range").at("min");
    massPlotMax = mcMassConfig.at("plot_range").at("max");
    massFitMin = mcMassConfig.at("fit_range").at("min");
    massFitMax = mcMassConfig.at("fit_range").at("max");

    std::string date_tag = global_config.value("date_tag", "test");
    std::string base_output_roots = global_config.value("base_output_dir_roots", "roots/");
    std::string base_output_figs = global_config.value("base_output_dir_figs", "fig/");

    std::string out_dir_name = Form("2DFit_%s/mc_Mass", date_tag.c_str());
    output_root_path = base_output_roots + out_dir_name;
    output_fig_path = base_output_figs + out_dir_name;

    TString kineLabel_tmp = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);
    output_base_name = Form("mc_Mass_%s", kineLabel_tmp.Data());

    gSystem->mkdir(output_root_path.c_str(), kTRUE);
    gSystem->mkdir(output_fig_path.c_str(), kTRUE);

    ws = new RooWorkspace("workspace", "MC Mass Fit Workspace");
  }
  catch (const json::parse_error &e)
  {
    std::cerr << "\nFATAL: JSON parsing error in " << e.what() << '\n'
              << "Message: " << e.what() << '\n'
              << "Location: Byte " << e.byte << " in one of the config files." << '\n'
              << "HINT: Please check for missing commas, extra commas, or mismatched brackets in your JSON files.\n";
    exit(1);
  }
  catch (const json::out_of_range &e)
  {
    std::cerr << "\nFATAL: JSON key not found: " << e.what() << '\n'
              << "HINT: Please check if all required keys (like 'McMassFit', 'kinematics', etc.) exist in your config files.\n";
    exit(1);
  }
}

void McMassFit::setupWorkspaceAndData()
{
  std::cout << "===== McMassFit::setupWorkspaceAndData() =====" << "\n";
  fInputData = new TFile(full_input_path.c_str());
  if (!fInputData || fInputData->IsZombie())
  {
    std::cerr << "FATAL: Input data file could not be opened: " << full_input_path << "\n";
    exit(1);
  }

  RooDataSet *ds_reduce = (RooDataSet *)fInputData->Get("dataset");
  ws->import(*ds_reduce, Rename("ds_reduce"));

  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && cBin>=%d && cBin<%d", ptLow, ptHigh, yLow, yHigh, cLow, cHigh);
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )";
  TString osCut = "recoQQsign==0";
  TString angleCut = Form("cos_ep>%.2f && cos_ep<%.2f", cosLow, cosHigh);
  TString finalCut = TString::Format("%s && %s && %s && %s", osCut.Data(), accCut.Data(), kineCut.Data(), angleCut.Data());

  RooArgSet varsToKeep;
  varsToKeep.add(*ws->var("mass"));
  varsToKeep.add(*ws->var("pt"));
  varsToKeep.add(*ws->var("y"));
  varsToKeep.add(*ws->var("cBin"));
  varsToKeep.add(*ws->var("recoQQsign"));
  varsToKeep.add(*ws->var("pt1"));
  varsToKeep.add(*ws->var("pt2"));
  varsToKeep.add(*ws->var("eta1"));
  varsToKeep.add(*ws->var("eta2"));
  varsToKeep.add(*ws->var("cos_ep"));

  bool apply_weights = local_config.at("McMassFit").value("apply_event_weights", false);
  RooDataSet *ds_weighted = nullptr;
  if (apply_weights)
  {
    std::cout << "INFO: Applying event weights." << "\n";
    ds_weighted = new RooDataSet("ds_weighted", "", ds_reduce, *ds_reduce->get(), finalCut.Data(), "weight");
    ds_reduce = ds_weighted;
  }

  RooDataSet *dsAB = (RooDataSet *)ds_reduce->reduce(varsToKeep, finalCut.Data());
  dsAB->SetName("dsAB");
  ws->import(*dsAB);

  if (ds_weighted)
    delete ds_weighted;
  delete dsAB;
}

void McMassFit::setVariableRanges()
{
  std::cout << "===== McMassFit::setVariableRanges() =====\n";
  ws->var("mass")->setRange("mcPlotRange", massPlotMin, massPlotMax);
  ws->var("mass")->setRange("mcFitRange", massFitMin, massFitMax);
}

void McMassFit::defineModel()
{
  std::cout << "===== McMassFit::defineModel() for type: " << signal_model_type << " =====" << "\n";

  // parameters
  for (auto const &[name, p] : yields.items())
  {
    ws->factory(Form("%s[%.4f, %.4f, %.4f]", name.c_str(),
                     p.at("val").get<double>(),
                     p.at("min").get<double>(),
                     p.at("max").get<double>()));
    ws->var(name.c_str())->SetTitle(p.at("title").get<std::string>().c_str());
  }

  for (auto const &[name, p] : signal_model_params.items())
  {
    ws->factory(Form("%s[%.4f, %.4f, %.4f]", name.c_str(),
                     p.at("val").get<double>(),
                     p.at("min").get<double>(),
                     p.at("max").get<double>()));
    ws->var(name.c_str())->SetTitle(p.at("title").get<std::string>().c_str());
  }

  std::vector<std::string> required_vars;
  if (signal_model_type == "DoubleCB")
  {
    required_vars = {"mass", "mean_jpsi", "sigma_1_A", "x_A", "alpha_1_A", "n_1_A", "f"};
  }
  else if (signal_model_type == "SingleGauss")
  {
    required_vars = {"mass", "mean", "sigma"};
  }

  for (const auto &var_name : required_vars)
  {
    if (!ws->var(var_name.c_str()))
    {
      std::cerr << "FATAL: Required variable '" << var_name << "' not found in workspace before building model!" << '\n';
      std::cerr << "HINT: Check if the variable is defined in the input RooDataSet or in the 'params' section of your JSON config." << '\n';
      ws->Print();
      exit(1);
    }
  }

  // build models
  signal_component_names.clear();
  bkg_component_names.clear();

  if (signal_model_type == "DoubleCB")
  {
    ws->factory("RooFormulaVar::sigma_2_A('@0*@1',{sigma_1_A, x_A})");
    ws->factory("RooFormulaVar::alpha_2_A('1.0*@0',{alpha_1_A})");
    ws->factory("RooFormulaVar::n_2_A('1.0*@0',{n_1_A})");

    RooCBShape cb_1_A("cball_1_A", "CB1", *ws->var("mass"), *ws->var("mean_jpsi"), *ws->var("sigma_1_A"), *ws->var("alpha_1_A"), *ws->var("n_1_A"));
    RooCBShape cb_2_A("cball_2_A", "CB2", *ws->var("mass"), *ws->var("mean_jpsi"), *ws->function("sigma_2_A"), *ws->function("alpha_2_A"), *ws->function("n_2_A"));

    ws->import(cb_1_A);
    ws->import(cb_2_A);
    signal_component_names.push_back("cball_1_A");
    signal_component_names.push_back("cball_2_A");

    RooAddPdf pdfMASS_Jpsi("pdfMASS_Jpsi", "Signal Shape", {cb_1_A, cb_2_A}, *ws->var("f"));
    ws->import(pdfMASS_Jpsi, RecycleConflictNodes());
  }
  else if (signal_model_type == "SingleGauss")
  {
    RooGaussian gauss("pdfMASS_Jpsi", "Signal Gauss", *ws->var("mass"), *ws->var("mean"), *ws->var("sigma"));

    ws->import(gauss);
    signal_component_names.push_back("pdfMASS_Jpsi");
  }
  else if (signal_model_type == "CBG")
  {
    RooCBShape cb_shape("cb_shape", "Crystal Ball Component", *ws->var("mass"), *ws->var("mean_jpsi"), *ws->var("sigma_cb"), *ws->var("alpha_cb"), *ws->var("n_cb"));
    RooGaussian gauss_shape("gauss_shape", "Gaussian Component", *ws->var("mass"), *ws->var("mean_jpsi"), *ws->var("sigma_gauss"));

    ws->import(cb_shape);
    ws->import(gauss_shape);
    signal_component_names.push_back("cb_shape");
    signal_component_names.push_back("gauss_shape");

    RooAddPdf pdfMASS_Jpsi("pdfMASS_Jpsi", "Signal Shape (CB+Gauss)", {cb_shape, gauss_shape}, *ws->var("frac_cb"));
    ws->import(pdfMASS_Jpsi, RooFit::RecycleConflictNodes());
  }
  else
  {
    std::cerr << "FATAL: Unknown signal model type '" << signal_model_type << "'" << "\n";
    exit(1);
  }

  RooAddPdf pdfMASS_Tot("pdfMASS_Tot", "Total MC PDF",
                        RooArgList(*ws->pdf("pdfMASS_Jpsi")),
                        RooArgList(*ws->var("N_Jpsi")));

  ws->import(pdfMASS_Tot, RooFit::RecycleConflictNodes());
}

void McMassFit::performFit()
{
  std::cout << "===== McMassFit::performFit() =====" << "\n";
  RooAbsPdf *pdf = ws->pdf("pdfMASS_Tot");
  RooAbsData *data = ws->data("dsAB");

  bool isWeighted = data->isWeighted();
  fitResult = pdf->fitTo(*data,
                         Save(),
                         Extended(kTRUE),
                         Range("mcFitRange"),
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

void McMassFit::makePlot()
{
  std::cout << "===== McMassFit::makePlot() =====" << '\n';

  // read plot setting
  json plot_cfg = local_config.at("McMassFit").at("plotting");
  int canvas_w = plot_cfg.value("canvas_size", json::object()).value("width", 800);
  int canvas_h = plot_cfg.value("canvas_size", json::object()).value("height", 800);
  double pad_split_ratio = plot_cfg.value("pad_split_ratio", 0.3);

  TCanvas *myCanvas = new TCanvas("myCanvas", "MC Mass Fit", canvas_w, canvas_h);

  // pad 1, 2
  TPad *pad1 = new TPad("pad1", "Main pad", 0.0, pad_split_ratio, 1.0, 1.0);
  pad1->SetBottomMargin(0.02);
  pad1->SetTicks(1, 1);
  pad1->Draw();

  TPad *pad2 = new TPad("pad2", "Pull pad", 0.0, 0.0, 1.0, pad_split_ratio);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.35);
  pad2->SetTicks(1, 1);
  pad2->Draw();

  // --- pad1 ---
  pad1->cd();
  pad1->SetLogy();

  // add frame, draw data and plots
  RooPlot *massFrame = ws->var("mass")->frame(Title(""), Bins(nMassBin), Range(massPlotMin, massPlotMax));

  ws->data("dsAB")->plotOn(massFrame, Name("dataOS"));
  ws->pdf("pdfMASS_Tot")->plotOn(massFrame, Name("pdfMASS_Tot"), LineColor(kBlack), Range("mcFitRange"), NormRange("mcFitRange"));

  json sig_style = plot_cfg.value("signal_components", json::object());
  std::vector<int> sig_colors = sig_style.value("colors", std::vector<int>{4, 2});
  std::vector<int> sig_styles = sig_style.value("styles", std::vector<int>{2, 7});
  for (size_t i = 0; i < signal_component_names.size(); ++i)
  {
    const std::string &name = signal_component_names[i];
    ws->pdf("pdfMASS_Tot")->plotOn(massFrame, Components(name.c_str()), Name(name.c_str()), LineColor(sig_colors[i % sig_colors.size()]), LineStyle(sig_styles[i % sig_styles.size()]), Range("mcFitRange"), NormRange("mcFitRange"));
  }

  // set y range
  TH1 *hTmp = ws->data("dsAB")->createHistogram("hTmp", *ws->var("mass"), Binning(nMassBin, massPlotMin, massPlotMax));
  Double_t YMax = hTmp->GetBinContent(hTmp->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= hTmp->GetNbinsX(); i++)
  {
    if (hTmp->GetBinContent(i) > 0)
    {
      YMin = std::min(YMin, hTmp->GetBinContent(i));
    }
  }
  double Ydown = YMin / (TMath::Power((YMax / YMin), (0.1 / (1.0 - 0.1 - 0.4))));
  double Yup = 2 * YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.1 - 0.4)));

  if (Ydown <= 0)
    Ydown = 0.1;

  massFrame->GetYaxis()->SetRangeUser(Ydown, Yup);
  delete hTmp;

  // draw frame
  massFrame->GetYaxis()->SetTitle(Form("Events / ( %.3f GeV/c^{2} )", (massPlotMax - massPlotMin) / nMassBin));
  massFrame->GetYaxis()->SetTitleOffset(1.4);
  massFrame->GetXaxis()->SetLabelSize(0);
  massFrame->GetXaxis()->SetTitleSize(0);
  massFrame->Draw();

  // draw parameters
  json kine_text_style = plot_cfg.value("kinematics_text", json::object());
  drawText(Form("%.1f < p_{T} < %.1f GeV/c", ptLow, ptHigh), kine_text_style.value("x", 0.18f), kine_text_style.value("y", 0.85f), kBlack, kine_text_style.value("size", 20));
  drawText(Form("%.1f < |y| < %.1f", yLow, yHigh), kine_text_style.value("x", 0.18f), kine_text_style.value("y", 0.85f) - kine_text_style.value("dy", 0.06f), kBlack, kine_text_style.value("size", 20));
  drawText(Form("Cent. %d-%d %%", cLow / 2, cHigh / 2), kine_text_style.value("x", 0.18f), kine_text_style.value("y", 0.85f) - 2 * kine_text_style.value("dy", 0.06f), kBlack, kine_text_style.value("size", 20));

  json param_text_style = plot_cfg.value("parameters_text", json::object());
  float param_y = param_text_style.value("y", 0.85f);
  float param_dy = param_text_style.value("dy", 0.06f);

  for (const auto &param_name : params_to_plot)
  {
    RooRealVar *param = ws->var(param_name.c_str());
    if (param)
    {
      TString title = param->GetTitle();
      double val = param->getVal();
      double err = param->getError();

      if (TString(param->GetName()).BeginsWith("N_"))
      {
        drawText(Form("%s = %.0f #pm %.0f", title.Data(), val, err),
                 param_text_style.value("x", 0.65f), param_y, kBlack, param_text_style.value("size", 20));
      }
      else
      {
        int precision = 3;
        double min_displayable_error = TMath::Power(10, -precision); // 0.001

        if (err > 0 && err < min_displayable_error)
        {

          drawText(Form("%s = %.*f #pm < %.*f", title.Data(), precision, val, precision, min_displayable_error),
                   param_text_style.value("x", 0.65f), param_y, kBlack, param_text_style.value("size", 20));
        }
        else
        {
          drawText(Form("%s = %.*f #pm %.*f", title.Data(), precision, val, precision, err),
                   param_text_style.value("x", 0.65f), param_y, kBlack, param_text_style.value("size", 20));
        }
      }
      param_y -= param_dy;
    }
  }

  json legend_style = plot_cfg.value("legend", json::object());
  TLegend *leg = new TLegend(legend_style.value("x1", 0.2f), legend_style.value("y1", 0.65f), legend_style.value("x2", 0.5f), legend_style.value("y2", 0.82f));
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->AddEntry(massFrame->findObject("dataOS"), "Data", "pe");
  leg->AddEntry(massFrame->findObject("pdfMASS_Tot"), "Total Fit", "l");
  for (const auto &name : signal_component_names)
  {
    leg->AddEntry(massFrame->findObject(name.c_str()), ws->pdf(name.c_str())->GetTitle(), "l");
  }
  leg->Draw("same");

  // --- pull pad ---
  pad2->cd();
  RooPlot *frameTMP = (RooPlot *)massFrame->Clone("frameTMP");
  RooPlot *pullFrame = ws->var("mass")->frame(Title(" "), Bins(nMassBin), Range(massPlotMin, massPlotMax));
  RooHist *pullHist = frameTMP->pullHist("dataOS", "pdfMASS_Tot", true);
  pullFrame->addPlotable(pullHist, "P");

  double pull_min = plot_cfg.value("pull_plot", json::object()).value("y_range_min", -5.0);
  double pull_max = plot_cfg.value("pull_plot", json::object()).value("y_range_max", 5.0);
  pullFrame->GetYaxis()->SetRangeUser(pull_min, pull_max);
  pullFrame->GetYaxis()->SetTitle("Pull");
  pullFrame->GetYaxis()->SetTitleSize(0.12);
  pullFrame->GetYaxis()->SetLabelSize(0.12);
  pullFrame->GetYaxis()->SetTitleOffset(0.35);
  pullFrame->GetYaxis()->SetNdivisions(505);
  pullFrame->GetYaxis()->CenterTitle();

  pullFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  pullFrame->GetXaxis()->SetTitleSize(0.12);
  pullFrame->GetXaxis()->SetLabelSize(0.12);
  pullFrame->GetXaxis()->SetTitleOffset(1.15);
  pullFrame->GetXaxis()->CenterTitle();
  pullFrame->Draw();

  TLine *l1 = new TLine(massPlotMin, 0, massPlotMax, 0);
  l1->SetLineStyle(kDashed);
  l1->Draw("same");
  printChi2(ws, pad2, frameTMP, fitResult, "mass", "dataOS", "pdfMASS_Tot", nMassBin, false);

  std::string fig_save_path = Form("%s/%s.png", output_fig_path.c_str(), output_base_name.c_str());
  myCanvas->SaveAs(fig_save_path.c_str());
  std::cout << "Plot saved to " << fig_save_path << '\n';

  delete frameTMP;
  delete myCanvas;
  delete massFrame;
  delete pullFrame;
}

void McMassFit::saveResults()
{
  std::cout << "===== McMassFit::saveResults() =====" << '\n';

  std::string root_save_path = Form("%s/%s_FitResult.root", output_root_path.c_str(), output_base_name.c_str());
  TFile *outputFile = new TFile(root_save_path.c_str(), "RECREATE");
  if (!outputFile || outputFile->IsZombie())
  {
    std::cerr << "ERROR: Could not create output file: " << root_save_path << '\n';
    return;
  }
  outputFile->cd();

  ws->pdf("pdfMASS_Tot")->Write();
  fitResult->Write("fitResult");

  // dataset noticing the list of parameters -> Do we need it?
  RooArgSet *fitargs = new RooArgSet();
  fitargs->add(fitResult->floatParsFinal());
  RooDataSet *datasetMass = new RooDataSet("datasetMass", "dataset with Mass Fit result", *fitargs);
  datasetMass->Write();

  outputFile->Close();

  std::cout << "Results saved to " << root_save_path << '\n';

  delete fitargs;
  delete datasetMass;
  delete outputFile;
}