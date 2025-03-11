#include <iostream>
#include <RooGaussian.h>
#include <RooFormulaVar.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"
#include "../headers/rootFitHeaders.h"
#include "../headers/commonUtility.h"
#include "../headers/JpsiUtility.h"
#include "../headers/cutsAndBin.h"
#include "../headers/CMS_lumi_v2mass.C"
#include "../headers/tdrstyle.C"



void bkg()
{
    TStopwatch *t = new TStopwatch;
    t->Start();
    cout << "\n=================================\n";
    cout << "\n Start ctau bkg fit\n";
    cout << "\n=================================\n";

    double nCtauBins = 200;
    // set macro configure
    using namespace std;
    using namespace RooFit;
    gStyle->SetEndErrorSize(0);

    // input parameters
    float cos_low = 0.0, cos_high = 0.1;
    float ptLow = 3;
    float ptHigh = 6.5;
    float yLow = 1.6;
    float yHigh = 2.4;
    int cLow = 60;
    int cHigh = 180;

    // Usually not used
    int PR = 0; // 0=PR, 1=NP, 2=Inc.
    int PRw = 1;
    bool fEffW = false;
    bool fAccW = false;
    bool isPtW = false;
    bool isTnP = false;
    double massLow = 2.6;
    double massHigh = 3.5; // Jpsi mass range

    // be quiet please
    // RooMsgService::instance().getStream(0).removeTopic(Caching); // in RooFit::
    // RooMsgService::instance().getStream(1).removeTopic(Caching);
    RooMsgService::instance().getStream(0).removeTopic(Plotting);
    RooMsgService::instance().getStream(1).removeTopic(Plotting);
    // RooMsgService::instance().getStream(0).removeTopic(Integration);
    // RooMsgService::instance().getStream(1).removeTopic(Integration);
    // RooMsgService::instance().getStream(0).removeTopic(InputArguments);
    // RooMsgService::instance().getStream(1).removeTopic(InputArguments);
    // RooMsgService::instance().getStream(0).removeTopic(NumIntegration);
    // RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
    // RooMsgService::instance().getStream(0).removeTopic(Minimization);
    // RooMsgService::instance().getStream(1).removeTopic(Minimization);
    RooMsgService::instance().setGlobalKillBelow(WARNING);


    // ===== labeling ===== //
    std::ostringstream out_file_path;
    out_file_path << "bkg/bkg_pt" << ptLow << "-" << ptHigh
                  << "_cent" << cLow << "-" << cHigh << "_absY" << yLow << "-" << yHigh
                  << "_cos" << cos_low << "_" << cos_high;
    string out_ss = out_file_path.str();

    std::ostringstream input_file;
    input_file << "pt" << ptLow << "-" << ptHigh
               << "_cent" << cLow << "-" << cHigh << "_absY" << yLow << "-" << yHigh
               << "_cos" << cos_low << "_" << cos_high
               << ".root";
    string in_ss = input_file.str();

    // PR, NP name label
    // TString fname = "";
    // if (PRw == 1) fname = "PR";
    // else if (PRw == 2) fname = "NP";

    // TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);


    // ===== make output folders ===== //
    gSystem->mkdir(Form("roots/bkg/"), kTRUE);
    gSystem->mkdir(Form("figs/bkg"), kTRUE);


    // ===== kinematic cut ===== //
    // cut is applied in f_err dataset


    // ===== import inputs ===== //
    auto f_mass = new TFile(("roots/mass/mass_" + in_ss).c_str());
    auto f_err = new TFile(("roots/err/err_" + in_ss).c_str());
    auto f_res = new TFile(("roots/res/res_" + in_ss).c_str());

    auto ds_red_mass = (RooDataSet *)f_mass->Get("ds_red_mass");
    auto mass_pdf = (RooAddPdf *)f_mass->Get("mass_pdf");
    auto ds_bkg_err = (RooDataSet *)f_err->Get("ds_bkg_err");
    auto err_bkg_pdf = (RooAddPdf *)f_err->Get("err_bkg_pdf");
    auto gaus_res_model = (RooAddPdf *)f_res->Get("gaus_res_model");

    RooWorkspace *ws_bkg = new RooWorkspace("ws_bkg");
    ws_bkg->import(*ds_red_mass);
    ws_bkg->import(*mass_pdf);
    ws_bkg->import(*ds_bkg_err);
    ws_bkg->import(*gaus_res_model);
    ws_bkg->import(*err_bkg_pdf);
    // ws_bkg->Print("V")


    // ===== set ranges ===== //
    double ctauErrMin, ctauErrMax;
    ws_bkg->data("ds_bkg_err")->getRange(*ws_bkg->var("ctau3DErr"), ctauErrMin, ctauErrMax);
    ctauLow=-3, ctauHigh=4;
    ws_bkg->var("ctau3D")->setRange(ctauLow, ctauHigh);
    ws_bkg->var("ctau3D")->setRange("ctauRange", ctauLow, ctauHigh);
    ws_bkg->var("ctau3DErr")->setRange(ctauErrMin, ctauErrMax);
    ws_bkg->var("ctau3DErr")->setRange("ctauErrRange", ctauErrMin, ctauErrMax);
    // ws_bkg->var("ctau3D")->Print();
    // ws_bkg->var("ctau3DErr")->Print();


    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // *-* Build model *-*//
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    // ===== make res model ===== //
    // fix Res parameters
    ws_bkg->var("ctau1_CtauRes")->setConstant(kTRUE);
    ws_bkg->var("s1_CtauRes")->setConstant(kTRUE);
    ws_bkg->var("ctau2_CtauRes")->setConstant(kTRUE);
    ws_bkg->var("rS21_CtauRes")->setConstant(kTRUE);
    // ws_bkg->var("ctau3_CtauRes")->setConstant(kTRUE);	ws_bkg->var("rS32_CtauRes")->setConstant(kTRUE);
    // ws_bkg->var("f2_CtauRes")->setConstant(kTRUE);
    ws_bkg->var("f_CtauRes")->setConstant(kTRUE);
    
    cout << "\nn_bkg: " << ws_bkg->var("n_bkg")->getVal() << "+/-" << ws_bkg->var("n_bkg")->getError() << endl;
    cout << "Bias of gaus1: " << ws_bkg->var("ctau1_CtauRes")->getVal() << endl;
    cout << "Fraction of gaus1: " << ws_bkg->var("f_CtauRes")->getVal() << "+/-" << ws_bkg->var("f_CtauRes")->getError() << endl;
    cout << "Sigma of gaus1: " << ws_bkg->var("s1_CtauRes")->getVal() << endl;

    ws_bkg->factory("zeroMean[0.0]");
    ws_bkg->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes1", "ctau3D",
                     "ctau1_CtauRes", //"ctau1_CtauRes",
                     "s1_CtauRes",
                     "One",
                     "ctau3DErr"));
    ws_bkg->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes2", "ctau3D",
                     "ctau2_CtauRes", //"ctau2_CtauRes",
                     "s2_CtauRes",
                     "One",
                     "ctau3DErr"));
    // ws_bkg->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes3", "ctau3D",
    //       "ctau3_CtauRes", //"ctau3_CtauRes",
    //       "s3_CtauRes",
    //       "zeroMean",
    //       "ctau3DErr"
    //       ));
    // ws_bkg->factory(Form("AddModel::%s({%s, %s}, {%s})", "ctauRes32", "ctauRes3", "ctauRes2", "f2_CtauRes"));
    ws_bkg->factory(Form("AddModel::%s({%s, %s}, {%s})", "pdfCTAURES", "ctauRes1", "ctauRes2", "f_CtauRes"));
    

    // ===== buld bkg decay model ===== //
    ws_bkg->factory("b_Bkg[0.1, 1e-3, 1]"); // NP fraction for bkg - Not final b fraciton
    ws_bkg->factory("fDFSS[0.1, 1e-3, 1.]");
    ws_bkg->factory("fDLIV[0.1, 1e-3, 1.]");
    ws_bkg->factory("lambdaDDS_Bkg[0.5, 0, 10]");
    ws_bkg->factory("lambdaDF_Bkg[ 0.5, 0, 10]");
    ws_bkg->factory("lambdaDSS_Bkg[0.5, 0, 10]");
    
    // make 3 exp
    ws_bkg->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUDSS1", "ctau3D", "lambdaDSS_Bkg", "pdfCTAURES"));
    ws_bkg->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUDSS2", "ctau3D", "lambdaDSS_Bkg2[0.1,0.01,1]", "pdfCTAURES"));
    ws_bkg->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUDSS", "fSS[0.5,0.01,1]", "pdfCTAUDSS1", "pdfCTAUDSS2"));

    ws_bkg->factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", "pdfCTAUDF1", "ctau3D", "lambdaDF_Bkg", "pdfCTAURES"));
    ws_bkg->factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", "pdfCTAUDF2", "ctau3D", "lambdaDF_Bkg2[0.1,0.01,1]", "pdfCTAURES"));
    ws_bkg->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUDF", "fFF[0.5,0.01,1]", "pdfCTAUDF1", "pdfCTAUDF2"));
    // ws_bkg->factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", "pdfCTAUDF1", "ctau3D", "lambdaDF_Bkg2[0.1,0.01,1]", "pdfCTAURES"));
    // ws_bkg->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUDF", "fFS[0.5,0.01,1]", "pdfCTAUDF1", "pdfCTAUDF1"));

    ws_bkg->factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", "pdfCTAUDDS", "ctau3D", "lambdaDDS_Bkg", "pdfCTAURES"));
    // ws_bkg->factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", "pdfCTAUDDS2", "ctau3D", "lambdaDDS_Bkg2[0.1,0.01,1]", "pdfCTAURES"));
    // ws_bkg->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUDDS", "fDDS[0.5,0.01,1]", "pdfCTAUDDS1", "pdfCTAUDDS2"));

    ws_bkg->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAU1", "fDFSS", "pdfCTAUDSS", "pdfCTAUDF"));
    ws_bkg->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUCOND_BkgNoPR", "fDLIV", "pdfCTAU1", "pdfCTAUDDS")); // NP

    // ctauBkgPR
    ws_bkg->factory(Form("SUM::%s(%s)", "pdfCTAUCOND_BkgPR", "pdfCTAURES"));
    ws_bkg->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUCOND_Bkg", "b_Bkg", "pdfCTAUCOND_BkgNoPR", "pdfCTAUCOND_BkgPR"));

    // total
    ws_bkg->factory(Form("RooExtendPdf::%s(%s,%s)", "pdfTot_Bkg", "pdfCTAUCOND_Bkg", "n_bkg")); // n_bkg is number of bkg from ds_bkg_err

    RooProdPdf pdfPR("pdfCTAU_BkgPR", "", *ws_bkg->pdf("err_bkg_pdf"), Conditional(*ws_bkg->pdf("pdfCTAUCOND_BkgPR"), RooArgList(*ws_bkg->var("ctau3D"))));
    ws_bkg->import(pdfPR);
    RooProdPdf pdfNoPR("pdfCTAU_BkgNoPR", "", *ws_bkg->pdf("err_bkg_pdf"), Conditional(*ws_bkg->pdf("pdfCTAUCOND_BkgNoPR"), RooArgList(*ws_bkg->var("ctau3D"))));
    ws_bkg->import(pdfNoPR);
    ws_bkg->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAU_Bkg", "b_Bkg", "pdfCTAU_BkgNoPR", "pdfCTAU_BkgPR"));
    RooAbsPdf *bkg_fit_pdf = new RooAddPdf("bkg_fit_pdf", "bkg_fit_pdf", RooArgList(*ws_bkg->pdf("pdfCTAU_Bkg")), RooArgList(*ws_bkg->var("n_bkg")));
    ws_bkg->import(*bkg_fit_pdf);
    // RooAbsPdf* ctauBkgModel = ctauBkgModel = new RooAddPdf("pdfTot_Bkg","pdfTot_Bkg",*ws_bkg->pdf("pdfCTAUCOND_Bkg"), *ws_bkg->var("n_bkg"));
    // ws_bkg->import(*pdfCTAUCOND_Bkg);


    // ===== preparation for drawing ===== //
    TH1D *h_tmp = (TH1D *)ws_bkg->data("ds_bkg_err")->createHistogram(("h_tmp"), *ws_bkg->var("ctau3D"), Binning(nCtauBins, ctauLow, ctauHigh));
    double ctauMin = h_tmp->GetBinLowEdge(h_tmp->FindFirstBinAbove(1, 1));
    if (cLow == 80 && cHigh == 100)
        ctauMin = -1.0;
    else if (cLow == 100 && cHigh == 180)
        ctauMin = -0.6;
    double ctauMax = h_tmp->GetBinLowEdge(h_tmp->FindLastBinAbove(2, 1)) + h_tmp->GetBinWidth(h_tmp->FindLastBinAbove(2, 1));
    // ws_bkg->var("ctau3D")->setRange("ctauFitRange", ctauLow, ctauHigh);

    TCanvas *c_bkg = new TCanvas("c_bkg", "My plots", 800, 800);
    c_bkg->cd();
    TPad *pad_bkg1 = new TPad("pad_bkg1", "pad_bkg1", 0, 0.16, 1, 1.0);
    pad_bkg1->SetTicks(1, 1);
    pad_bkg1->Draw();
    pad_bkg1->cd();
    RooPlot *ctau_frame1 = ws_bkg->var("ctau3D")->frame(Bins(nCtauBins), Range(ctauLow, ctauHigh)); // bins
    ctau_frame1->SetTitle("");

    ws_bkg->pdf("bkg_fit_pdf")->setNormRange("ctauWindow");


    // ===== fit here ===== //
    // // buld dataset for fit
    // RooDataSet* dataToFit = (RooDataSet*)ds_bkg_err->reduce(Form("ctau3D>=%.f&&ctau3D<=%.f",ctauMin, ctauMax))->Clone("ds_bkg_err");
    RooDataSet *dataToFit = (RooDataSet *)ds_bkg_err->reduce(Form("ctau3D>=%.f&&ctau3D<=%.f", ctauLow, ctauHigh))->Clone("dataToFit");
    ws_bkg->import(*dataToFit);

    pad_bkg1->cd();
    gPad->SetLogy();
    RooPlot *ctau_frame2 = (RooPlot *)ctau_frame1->Clone("ctau_frame2");

    double normDSTot = 1.0;
    if (ws_bkg->data("dataToFit"))
    {
        normDSTot = ws_bkg->data("dataToFit")->sumEntries() / ws_bkg->data("dataToFit")->sumEntries();
    }
    double normBkg = 1.0;
    if (ws_bkg->data("ds_bkg_err"))
    {
        normBkg = ws_bkg->data("dataToFit")->sumEntries() * normDSTot / ws_bkg->data("ds_bkg_err")->sumEntries();
    }


    // to avoid ternimal breaking due to too long Tracing Error message
    // It's okay as long as minimizer find good fitting point
    // Always check fit results.
    RooMsgService::instance().getStream(0).removeTopic(Tracing);

    bool isWeighted = ws_bkg->data("ds_bkg_err")->isWeighted();
    // RooFitResult* fit_bkg = ws_bkg->pdf("bkg_fit_pdf")->fitTo(*ds_bkg_err, Save(), Range("ctauRange"), Extended(kTRUE), NumCPU(12));


    RooFitResult *fit_bkg = ws_bkg->pdf("bkg_fit_pdf")->fitTo(*dataToFit, Save(), Extended(1), NumCPU(12), SumW2Error(isWeighted), Strategy(2), RecoverFromUndefinedRegions(1.), PrintEvalErrors(-1), PrintLevel(0));
    
    ws_bkg->import(*fit_bkg, "fit_bkg");

    // fis parameters
    ws_bkg->var("b_Bkg")->setConstant(kTRUE);
    ws_bkg->var("fDFSS")->setConstant(kTRUE);
    ws_bkg->var("fDLIV")->setConstant(kTRUE);
    ws_bkg->var("lambdaDDS_Bkg")->setConstant(kTRUE);
    ws_bkg->var("lambdaDF_Bkg")->setConstant(kTRUE);
    ws_bkg->var("lambdaDSS_Bkg")->setConstant(kTRUE);
    ws_bkg->var("fSS")->setConstant(kTRUE);

    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // *-* Plotting - very long but simple *-*//
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    ctau_frame2->updateNormVars(RooArgSet(*ws_bkg->var("ctau3D")));
    // ctau_frame2->updateNormVars(RooArgSet(*ws_bkg->var("ctau3D")));

    // ws_bkg->data("ds_bkg_err")->plotOn(ctau_frame2, Name("data_ctauBkg"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.7)); // to show full ranges

    RooDataSet *expDataY = (RooDataSet *)dataToFit->reduce(*ws_bkg->var("ctau3DErr"));
    RooAbsData *binnedDataY = expDataY->binnedClone();

    // draw data point to provide normalization info to PDFs
    ws_bkg->data("dataToFit")->plotOn(ctau_frame2, Name("data_ctauBkg"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.7)); // data points inside range

    // ws_bkg->pdf("bkg_fit_pdf")->plotOn(ctau_frame2, ProjWData(*binnedDataY), LineColor(kCyan), LineStyle(kDotted));
    // ws_bkg->pdf("bkg_fit_pdf")->plotOn(ctau_frame2, LineColor(kCyan), LineStyle(kDotted));
    // ws_bkg->pdf("bkg_fit_pdf")->forceNumInt(kTRUE);
    ws_bkg->pdf("bkg_fit_pdf")->plotOn(ctau_frame2, Name("ctauBkg_Tot"), NumCPU(12), LineColor(kBlack));

    // draw data point again - to bring point forward
    // ws_bkg->data("dataToFit")->plotOn(ctau_frame2, Name("data_ctauBkg"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.7)); // data points inside range

    if (ws_bkg->pdf("pdfCTAU_BkgPR"))
    {
        ws_bkg->pdf("bkg_fit_pdf")->plotOn(ctau_frame2, Name("BKGPR"), Components(RooArgSet(*ws_bkg->pdf("pdfCTAU_BkgPR"))), LineColor(44), NumCPU(12), Range("ctauRange"), LineStyle(kDashed), LineWidth(2));
    }
    if (ws_bkg->pdf("pdfCTAU_BkgNoPR"))
    {
        ws_bkg->pdf("bkg_fit_pdf")->plotOn(ctau_frame2, Name("BKGNoPR"), Components(RooArgSet(*ws_bkg->pdf("pdfCTAU_BkgNoPR"))), LineColor(8), NumCPU(12), Range("ctauRange"), LineStyle(kDashed), LineWidth(2));
    }
    // ws_bkg->pdf("bkg_fit_pdf")->plotOn(ctau_frame2, Name("PDF"), , LineColor(kBlack), Precision(1e-4), Range("ctauRange"));
    // Normalization(normBkg, RooAbsReal::NumEvent), NumCPU(12), ProjWData(RooArgSet(*ws_bkg->var("ctau3DErr")), *ws_bkg->data("ds_bkg_err"), kTRUE)

    // set proper y range
    ctau_frame2->GetYaxis()->SetRangeUser(10e-2, 10e7);
    TH1 *h = ws_bkg->data("ds_bkg_err")->createHistogram("hist", *ws_bkg->var("ctau3D"), Binning(ctau_frame1->GetNbinsX(), ctau_frame1->GetXaxis()->GetXmin(), ctau_frame1->GetXaxis()->GetXmax()));
    Double_t YMax = h->GetBinContent(h->GetMaximumBin());
    Double_t YMin = 1e99;
    for (int i = 1; i <= h->GetNbinsX(); i++)
        if (h->GetBinContent(i) > 0)
            YMin = min(YMin, h->GetBinContent(i));
    Double_t Yup(0.), Ydown(0.);
    // Yup = YMax*TMath::Power((YMax/0.1), 0.5);
    Yup = YMax * TMath::Power((YMax / 0.01), 0.5);
    Ydown = 0.01;
    ctau_frame2->GetYaxis()->SetRangeUser(Ydown, Yup);
    ctau_frame2->GetXaxis()->SetRangeUser(-4, 7);
    ctau_frame2->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
    ctau_frame2->SetFillStyle(4000);
    ctau_frame2->GetYaxis()->SetTitleOffset(1.43);
    ctau_frame2->GetXaxis()->SetLabelSize(0);
    ctau_frame2->GetXaxis()->SetTitleSize(0);


    // ===== draw lines to present Res rantge ===== //
    TLine *minline = new TLine(ctauMin, 0.0, ctauMin, Ydown * TMath::Power((Yup / Ydown), 0.4));
    minline->SetLineStyle(2);
    minline->SetLineColor(1);
    minline->SetLineWidth(3);
    ctau_frame2->addObject(minline);
    TLine *maxline = new TLine(ctauMax, 0.0, ctauMax, Ydown * TMath::Power((Yup / Ydown), 0.4));
    maxline->SetLineStyle(2);
    maxline->SetLineColor(1);
    maxline->SetLineWidth(3);
    ctau_frame2->addObject(maxline);


    // ===== legend ===== //
    ctau_frame2->Draw();
    TLegend *leg_bkg = new TLegend(text_x + 0.25, text_y + 0.04, text_x + 0.38, text_y - 0.15);
    leg_bkg->SetTextSize(text_size);
    leg_bkg->SetTextFont(43);
    leg_bkg->SetBorderSize(0);
    leg_bkg->AddEntry(ctau_frame2->findObject("data_ctauBkg"), "Data_Bkg", "pe");
    // leg_bkg->AddEntry(ctau_frame2->findObject("BKGPR"), "Bkg PR", "fl");
    // leg_bkg->AddEntry(ctau_frame2->findObject("BKGNoPR"), "Bkg NonPR", "fl");
    leg_bkg->Draw("same");


    // ===== draw latex ===== //
    drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x, text_y, text_color, text_size);
    if (yLow == 0)
        drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x, text_y - y_diff, text_color, text_size);
    else if (yLow != 0)
        drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x, text_y - y_diff, text_color, text_size);
    drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x, text_y - y_diff * 2, text_color, text_size);
    drawText(Form("N_{Bkg} = %.f #pm %.f", ws_bkg->var("n_bkg")->getVal(), ws_bkg->var("n_bkg")->getError()), text_x + 0.5, text_y, text_color, text_size);
    drawText(Form("b_{Bkg} = %.4f #pm %.4f", ws_bkg->var("b_Bkg")->getVal(), ws_bkg->var("b_Bkg")->getError()), text_x + 0.5, text_y - y_diff * 1, text_color, text_size);
    drawText(Form("fDFSS = %.4f #pm %.4f", ws_bkg->var("fDFSS")->getVal(), ws_bkg->var("fDFSS")->getError()), text_x + 0.5, text_y - y_diff * 2, text_color, text_size);
    drawText(Form("fDLIV = %.4f #pm %.4f", ws_bkg->var("fDLIV")->getVal(), ws_bkg->var("fDLIV")->getError()), text_x + 0.5, text_y - y_diff * 3, text_color, text_size);
    drawText(Form("#lambdaDDS_{Bkg} = %.4f #pm %.4f", ws_bkg->var("lambdaDDS_Bkg")->getVal(), ws_bkg->var("lambdaDDS_Bkg")->getError()), text_x + 0.5, text_y - y_diff * 4, text_color, text_size);
    drawText(Form("#lambdaDF_{Bkg} = %.4f #pm %.4f", ws_bkg->var("lambdaDF_Bkg")->getVal(), ws_bkg->var("lambdaDF_Bkg")->getError()), text_x + 0.5, text_y - y_diff * 5, text_color, text_size);
    drawText(Form("#lambdaDSS_{Bkg} = %.4f #pm %.4f", ws_bkg->var("lambdaDSS_Bkg")->getVal(), ws_bkg->var("lambdaDSS_Bkg")->getError()), text_x + 0.5, text_y - y_diff * 6, text_color, text_size);
    
 
    // ===== draw pull ===== //
    TPad *pull_pad = new TPad("pull_pad", "pull_pad", 0, 0.006, 1, 0.227);
    c_bkg->cd();
    pull_pad->Draw();
    pull_pad->cd();
    pull_pad->SetTopMargin(0); // Upper and lower plot are joined
    pull_pad->SetBottomMargin(0.67);
    pull_pad->SetBottomMargin(0.4);
    pull_pad->SetFillStyle(4000);
    pull_pad->SetFrameFillStyle(4000);
    pull_pad->SetTicks(1, 1);

    RooPlot *tmp_pull_frame = (RooPlot *)ctau_frame2->Clone("tmp_pull_frame");
    RooHist *pull_hist = tmp_pull_frame->pullHist("data_ctauBkg", "ctauBkg_Tot", true);
    // RooHist *pull_hist = tmp_pull_frame->pullHist("data_ctauBkg", "ctauBkg_Tot", true); // flow code
    pull_hist->SetMarkerSize(0.8);
    RooPlot *pull_frame = ws_bkg->var("ctau3D")->frame(Title("Pull Distribution"), Bins(nCtauBins), Range(ctauLow, ctauHigh));
    pull_frame->addPlotable(pull_hist, "PX");
    pull_frame->SetTitle("");
    pull_frame->SetTitleSize(0);
    pull_frame->GetYaxis()->SetTitleOffset(0.3);
    pull_frame->GetYaxis()->SetTitle("Pull");
    pull_frame->GetYaxis()->SetTitleSize(0.08);
    pull_frame->GetYaxis()->SetLabelSize(0.08);
    pull_frame->GetYaxis()->SetRangeUser(-50, 50);
    pull_frame->GetYaxis()->CenterTitle();
    pull_frame->GetYaxis()->SetTickSize(0.04);
    pull_frame->GetYaxis()->SetNdivisions(404);
    pull_frame->GetXaxis()->SetTickSize(0.03);

    pull_frame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
    pull_frame->GetXaxis()->SetTitleOffset(1.55);
    pull_frame->GetXaxis()->SetLabelOffset(0.04);
    pull_frame->GetXaxis()->SetLabelSize(0.08);
    pull_frame->GetXaxis()->SetTitleSize(0.08);
    pull_frame->GetXaxis()->CenterTitle();
    pull_frame->Draw();

    TLine *lD = new TLine(ctauLow, 0, ctauHigh, 0);
    lD->SetLineStyle(1);
    lD->Draw("same");


    // ===== horizontoal lines on pull hist ===== //
    // maybe we can make function for them...
    TLine *line_at_2 = new TLine(pull_frame->GetXaxis()->GetXmin(), 2, pull_frame->GetXaxis()->GetXmax(), 2);
    line_at_2->SetLineColor(kBlack);
    line_at_2->SetLineStyle(3);
    line_at_2->SetLineWidth(1);
    line_at_2->Draw();

    TLine *line_at_4 = new TLine(pull_frame->GetXaxis()->GetXmin(), 4, pull_frame->GetXaxis()->GetXmax(), 4);
    line_at_4->SetLineColor(kBlack);
    line_at_4->SetLineStyle(3);
    line_at_4->SetLineWidth(2);
    line_at_4->Draw();

    TLine *line_at_8 = new TLine(pull_frame->GetXaxis()->GetXmin(), 8, pull_frame->GetXaxis()->GetXmax(), 8);
    line_at_8->SetLineColor(kBlack);
    line_at_8->SetLineStyle(3);
    line_at_8->SetLineWidth(1);
    line_at_8->Draw();

    // negative
    TLine *line_at_n2 = new TLine(pull_frame->GetXaxis()->GetXmin(), -2, pull_frame->GetXaxis()->GetXmax(), -2);
    line_at_n2->SetLineColor(kBlack);
    line_at_n2->SetLineStyle(3);
    line_at_n2->SetLineWidth(1);
    line_at_n2->Draw();

    TLine *line_at_n4 = new TLine(pull_frame->GetXaxis()->GetXmin(), -4, pull_frame->GetXaxis()->GetXmax(), -4);
    line_at_n4->SetLineColor(kBlack);
    line_at_n4->SetLineStyle(3);
    line_at_n4->SetLineWidth(2);
    line_at_n4->Draw();

    TLine *line_at_n8 = new TLine(pull_frame->GetXaxis()->GetXmin(), -8, pull_frame->GetXaxis()->GetXmax(), -8);
    line_at_n8->SetLineColor(kBlack);
    line_at_n8->SetLineStyle(3);
    line_at_n8->SetLineWidth(1);
    line_at_n8->Draw();


    // ===== chi2 ===== //
    printChi2(ws_bkg, pull_pad, tmp_pull_frame, fit_bkg, "ctau3D", "data_ctauBkg", "ctauBkg_Tot", nCtauBins, false);

    pull_pad->Update();


    // ===== draw canvas ===== //
    c_bkg->Update();
    c_bkg->SaveAs(("figs/" + out_ss + ".png").c_str());
    c_bkg->SaveAs(("figs/" + out_ss + ".pdf").c_str());

    // ===== test2 ===== //
    // // Get ROOT histograms
    // TH1 *h_func = ws_bkg->pdf("bkg_fit_pdf")->createHistogram("h_func", *ws_bkg->var("ctau3D"), Binning(ctau_frame1->GetNbinsX(), ctau_frame1->GetXaxis()->GetXmin(), ctau_frame1->GetXaxis()->GetXmax()));
    // TH1 *h_data = dataToFit->createHistogram("h_data", *ws_bkg->var("ctau3D"), Binning(ctau_frame1->GetNbinsX(), ctau_frame1->GetXaxis()->GetXmin(), ctau_frame1->GetXaxis()->GetXmax()));
    // double KStest = h_data->KolmogorovTest(h_func);
    // TLatex *latex = new TLatex(); // prepare text in LaTeX format
    // latex->SetTextSize(0.035);
    // latex->SetNDC();
    // latex->DrawLatex(0.25, 0.75, Form("K-S test = %.2f", KStest)); // K-S test = 1 means very high probability of data coming from the distribution described by the model
    // h_data->SetMarkerColor(kRed);
    // TCanvas *cnv_KS = new TCanvas();
    // cnv_KS->cd();
    // gPad->SetLogy();
    // h_data->Draw();
    // h_func->Draw("same");
    // cout << "KStest: " << KStest << endl;

    // ===== ratio test cnavas ===== //
    auto c = new TCanvas();
    c->cd();
    TH1 *h2 = dataToFit->createHistogram("h2", *ws_bkg->var("ctau3D"), Binning(ctau_frame1->GetNbinsX(), ctau_frame1->GetXaxis()->GetXmin(), ctau_frame1->GetXaxis()->GetXmax()));
    // auto h1 = data->createHistogram("x");

    auto f1 = ws_bkg->pdf("bkg_fit_pdf")->asTF(*ws_bkg->var("ctau3D"), RooArgSet(), *ws_bkg->var("ctau3D"));
    auto fcor = new TF1("fcor", [&](double *x, double *p)
                        { return f1->EvalPar(x, p) * h2->GetSumOfWeights() * h2->GetBinWidth(1); }, ctau_frame1->GetXaxis()->GetXmin(), ctau_frame1->GetXaxis()->GetXmax(), 0);
    auto f2 = (TF1 *)fcor->Clone("f2"); // for drawing
    double f2_max = f2->Eval(0);
    cout << "Eval f2 at ctau3D = 0: " << f2->Eval(0) << endl;
    h2->SetMaximum(f2_max * 5); // usually, f2's maximum is bigger than a hist's
    gPad->SetLogy();
    h2->Draw("");
    f2->Draw("same");

    int n_tmp = h->GetNbinsX();
    for (int bin = 1; bin <= n_tmp; ++bin)
    {
        Double_t binCenter = h2->GetBinCenter(bin);
        Double_t result = f2->Eval(binCenter);
        Double_t binContent = h2->GetBinContent(bin);

        cout << "Bin center: " << binCenter
             << ", Bin content: " << binContent
             << ", f(" << binCenter << ") = " << result << std::endl;
    }

    Double_t chi2 = h2->Chisquare(f2, "L");
    cout << "Chi2: " << chi2 << endl;
    c->Draw();
    c->SaveAs(("figs/" + out_ss + "_TF1.png").c_str());

    // ===== Export results ===== //
    auto out_file = new TFile(("roots/" + out_ss + ".root").c_str(), "recreate");

    fit_bkg->Write();
    ws_bkg->Write();
    out_file->Close();

    fit_bkg->Print("V");

    cout << "\n=================================\n";
    cout << "\n Finish ctau bkg fit\n";
    cout << "\n=================================\n";
    t->Stop();
    printf("RealTime=%f seconds, CpuTime=%f seconds\n", t->RealTime(), t->CpuTime());
}