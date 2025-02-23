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
#include "../headers/cutsAndBin.h"
#include "../headers/CMS_lumi_v2mass.C"
#include "../headers/tdrstyle.C"
#include "../headers/rootFitHeaders.h"
#include "../headers/commonUtility.h"
#include "../headers/JpsiUtility.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"


void err()
{
    TStopwatch *t = new TStopwatch;
    t->Start();
    cout << "\n=================================\n";
    cout << "\n Start ctau3DErr sPlot work\n";
    cout << "\n=================================\n";

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
    RooMsgService::instance().getStream(0).removeTopic(Caching); // in RooFit::
    RooMsgService::instance().getStream(1).removeTopic(Caching);
    RooMsgService::instance().getStream(0).removeTopic(Plotting);
    RooMsgService::instance().getStream(1).removeTopic(Plotting);
    RooMsgService::instance().getStream(0).removeTopic(Integration);
    RooMsgService::instance().getStream(1).removeTopic(Integration);
    RooMsgService::instance().getStream(0).removeTopic(InputArguments);
    RooMsgService::instance().getStream(1).removeTopic(InputArguments);
    RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
    RooMsgService::instance().getStream(1).removeTopic(Minimization);
    RooMsgService::instance().setGlobalKillBelow(WARNING);


    // ===== labeling ===== //
    std::ostringstream out_file_path;
    out_file_path << "err/err_pt" << ptLow << "-" << ptHigh
                  << "_cent" << cLow << "-" << cHigh << "_absY" << yLow << "-" << yHigh
                  << "_cos" << cos_low << "_" << cos_high;
    string out_ss = out_file_path.str();

    // PR, NP name label
    // TString fname = "";
    // if (PRw == 1) fname = "PR";
    // else if (PRw == 2) fname = "NP";

    // TString kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);


    // ===== make output folders ===== //
    gSystem->mkdir(Form("roots/err/"), kTRUE);
    gSystem->mkdir(Form("figs/err"), kTRUE);


    // ===== kinematic cut ===== //
    TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f && cBin>=%d && cBin<%d", ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);
    TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&"; // 2018 acceptance cut
    TString OS = "recoQQsign==0 &&";

    TString angle_cut = Form("&& cos_ep>%.2f && cos_ep<%.2f", cos_low, cos_high);

    TString final_cut = OS + accCut + kineCut + angle_cut;


    // ===== import inputs ===== //
    // mass ds_splot
    auto f_data = new TFile("../files_roodata/RooDataSet_miniAOD_isMC0_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root");
    f_data->SetName("f_data");
    RooDataSet *ds_tmp = (RooDataSet *)f_data->Get("dataset");

    // import mass fit result
    std::ostringstream in_mass_fit;
    in_mass_fit << "roots/mass/mass_pt" << ptLow << "-" << ptHigh
              << "_cent" << cLow << "-" << cHigh << "_absY" << yLow << "-" << yHigh
              << "_cos" << cos_low << "_" << cos_high
              << ".root";
    string ss = in_mass_fit.str();
    auto f_mass = new TFile(ss.c_str());

    auto ds_mass = (RooDataSet *)f_mass->Get("ds_red_mass");
    auto mass_pdf = (RooAddPdf *)f_mass->Get("mass_pdf");


    // ===== make datasets and workspaces ===== //
    RooWorkspace *ws_tmp = new RooWorkspace("ws_tmp");
    ws_tmp->import(*ds_tmp);

    RooArgSet *argSet = new RooArgSet(*(ws_tmp->var("ctau3D")), *(ws_tmp->var("mass")), *(ws_tmp->var("pt")), *(ws_tmp->var("y")), *(ws_tmp->var("weight")), *(ws_tmp->var("ctau3DRes")), *(ws_tmp->var("ctau3DErr")));
    argSet->add(*(ws_tmp->var("pt1")));
    argSet->add(*(ws_tmp->var("pt2")));
    argSet->add(*(ws_tmp->var("eta1")));
    argSet->add(*(ws_tmp->var("eta2")));
    argSet->add(*(ws_tmp->var("recoQQsign")));
    argSet->add(*(ws_tmp->var("cBin")));
    argSet->add(*(ws_tmp->var("cos_ep")));

    RooDataSet *ds_tmpW = new RooDataSet("ds_tmpW", "weighted",
                                          *argSet,
                                          Import(*ds_tmp), WeightVar(*ws_tmp->var("weight")));
    // not used?
    // RooDataSet *datasetWo = new RooDataSet("datasetWo", "A sample",
    //                                        *argSet,
    //                                        Import(*ds_tmp));
    // RooDataSet *ds_data1 = (RooDataSet *)datasetWo->reduce(*argSet, kineCut.Data());

    // it will be cloned to ds_splot
    RooDataSet *ds_data = (RooDataSet *)ds_tmpW->reduce(*argSet, kineCut.Data());
    ds_data->SetName("ds_data");

    // workspace to use from now on
    RooWorkspace *ws_err = new RooWorkspace("ws_err");
    ws_err->import(*ds_mass);
    ws_err->import(*mass_pdf);
    ws_err->import(*ds_data);

    // set ctau3DErr range and bins - from legacy code
    // ws_err->var("ctau3DErr")->setRange(0,0.25);
    ws_err->var("ctau3DErr")->setRange(ctauErrLow, ctauHigh);
    int nBins = min(int(round((ws_err->var("ctau3DErr")->getMax() - ws_err->var("ctau3DErr")->getMin()) / 0.0025)), 100);
    cout << "nBin : " << nBins << endl;


    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // *-* Apply sPlot technique *-*//
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    
    // n_sig, n_bkg will be used for sPlot
    RooArgList yieldList;
    RooRealVar *sigYield = ws_err->var("n_sig");
    RooRealVar *bkgYield = ws_err->var("n_bkg");
    yieldList.add(*ws_err->var("n_sig"));
    yieldList.add(*ws_err->var("n_bkg"));

    cout << "Sig Yield: " << sigYield->getVal() << " +/- " << sigYield->getError() << endl;
    cout << "Bkg Yield: " << bkgYield->getVal() << " +/- " << bkgYield->getError() << endl;

    // Signal and Bkg PDFs of mass fit will be used
    // code below is messy, I keep it because it's legacy
    RooDataSet *ds_splot = (RooDataSet *)ws_err->data("ds_data")->Clone("ds_splot");
    RooArgSet *cloneSet = (RooArgSet *)RooArgSet(*ws_err->pdf("mass_pdf"), "mass_pdf").snapshot(kTRUE);
    auto clone_mass_pdf = (RooAbsPdf *)cloneSet->find("mass_pdf");
    clone_mass_pdf->setOperMode(RooAbsArg::ADirty, kTRUE);


    // ===== sPlot fit ===== //
    //sPlot: fds_splotit variable O1 to get the weightings of sig and bkg. We will aplly the weightings to variable O2 later.
    RooStats::SPlot sData = RooStats::SPlot("sData", "An SPlot", *ds_splot, clone_mass_pdf, yieldList);

    // the total yield and total weight should be same
    cout << "[INFO] Jpsi yield -> Mass Fit:" << ws_err->var("n_sig")->getVal() << ", sWeights :" << sData.GetYieldFromSWeight("n_sig") << endl;
    cout << "[INFO] Bkg  yield -> Mass Fit:" << ws_err->var("n_bkg")->getVal() << ", sWeights :" << sData.GetYieldFromSWeight("n_bkg") << endl;

    // individual events have different weights
    for (int i = 0; i < 10; i++)
    {
        // check weights of some events
        std::cout << "n_sig weight: " << i << std::right << std::setw(12) << sData.GetSWeight(i, "n_sig") << " /  n_bkg weight: "
                  << std::setw(12) << sData.GetSWeight(i, "n_bkg") << " / Total Weight" << std::setw(12) << sData.GetSumOfEventSWeight(i)
                  << std::endl;
    }

    ws_err->import(*ds_splot);
    

    
    // ===== build RooPlot and histograms ===== //
    // a lot of legacy codes
    RooPlot *err_frame1 = ws_err->var("ctau3DErr")->frame(Range(ctauErrLow, ctauErrHigh));

    TH1D *hTot = (TH1D *)ws_err->data("ds_data")->createHistogram(("hTot"), *ws_err->var("ctau3DErr"), Binning(err_frame1->GetNbinsX(), err_frame1->GetXaxis()->GetXmin(), err_frame1->GetXaxis()->GetXmax()));
    TH1D *hTot_M = (TH1D *)ws_err->data("ds_data")->createHistogram(("hTot_M"), *ws_err->var("ctau3DErr"), Binning(nBins, ctauErrLow, ctauErrHigh));


    // ===== build temporal sPlot weighted datasets ===== //
    RooDataSet *dataw_Sig_b = new RooDataSet("dataw_Sig_b", "TMP_SIG_DATA", (RooDataSet *)ws_err->data("ds_splot"),
                                             RooArgSet(*ws_err->var("ctau3DErr"), *ws_err->var("n_sig_sw"), *ws_err->var("ctau3DRes"), *ws_err->var("ctau3D"), *ws_err->var("mass")), 0, "n_sig_sw");
    TH1D *hSig = (TH1D *)dataw_Sig_b->createHistogram(("hSig"), *ws_err->var("ctau3DErr"), Binning(nBins, ctauErrLow, ctauErrHigh));

    // Bkg
    RooDataSet *dataw_Bkg_b = new RooDataSet("dataw_Bkg_b", "TMP_BKG_DATA", (RooDataSet *)ws_err->data("ds_splot"),
                                             RooArgSet(*ws_err->var("ctau3DErr"), *ws_err->var("n_bkg_sw"), *ws_err->var("ctau3DRes"), *ws_err->var("ctau3D"), *ws_err->var("mass")), 0, "n_bkg_sw");
    TH1D *hBkg = (TH1D *)dataw_Bkg_b->createHistogram(("hBkg"), *ws_err->var("ctau3DErr"), Binning(nBins, ctauErrLow, ctauErrHigh));


    // ===== set new ranges of ctau3DErr ===== //
    double ctauErrMin;
    double ctauErrMax;
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

    cout << "\nSet new bin width\n";
    double BinWidth = (ctauErrHigh - ctauErrLow) / nBins;
    cout << "BinWidth : " << BinWidth << endl;
    int newBins = (ctauErrMax - ctauErrMin) / BinWidth;
    cout << "newBins : " << newBins << "\n\n";


    // ===== build sPlot weighted datasets and RooHistPdfs ===== //
    // total ds and pdf
    TH1D *h_tot_test = (TH1D *)ws_err->data("ds_data")->createHistogram(("h_tot_test"), *ws_err->var("ctau3DErr"), Binning(newBins, ctauErrMin, ctauErrMax));
    // TH1D* h_tot_test = (TH1D*)ws_err->data("ds_data")->createHistogram(("h_tot_test"), *ws_err->var("ctau3DErr"),Binning(nBins,ctauErrMin,ctauErrMax));
    RooDataHist *totHist = new RooDataHist("ds_data", "", *ws_err->var("ctau3DErr"), h_tot_test);
    RooHistPdf *err_tot_pdf = new RooHistPdf("err_tot_pdf", "Total RooHistPdf of ctauErr", *ws_err->var("ctau3DErr"), *totHist);

    // bkg ds and pdf
    TH1D *hBkg_w = (TH1D *)dataw_Bkg_b->createHistogram(("hBkg_w"), *ws_err->var("ctau3DErr"), Binning(newBins, ctauErrMin, ctauErrMax));
    RooDataHist *bkgHist = new RooDataHist("ds_data", "", *ws_err->var("ctau3DErr"), hBkg_w);
    RooHistPdf *err_bkg_pdf = new RooHistPdf("err_bkg_pdf", "hist pdf", *ws_err->var("ctau3DErr"), *bkgHist);
    
    // sig ds and pdf
    TH1D *hSig_w = (TH1D *)dataw_Sig_b->createHistogram(("hSig_w"), *ws_err->var("ctau3DErr"), Binning(newBins, ctauErrMin, ctauErrMax));
    RooDataHist *sigHist = new RooDataHist("ds_data", "", *ws_err->var("ctau3DErr"), hSig_w);
    RooHistPdf *err_sig_pdf = new RooHistPdf("err_sig_pdf", "hist pdf", *ws_err->var("ctau3DErr"), *sigHist);
    
    // Re-build datasets
    RooDataSet *ds_bkg_err = new RooDataSet("ds_bkg_err", "", (RooDataSet *)ws_err->data("ds_splot"),
                                           RooArgSet(*ws_err->var("ctau3DErr"), *ws_err->var("n_bkg_sw"), *ws_err->var("ctau3DRes"), *ws_err->var("ctau3D"), *ws_err->var("mass")), 0, "n_bkg_sw");
    RooDataSet *ds_sig_err = new RooDataSet("ds_sig_err", "", (RooDataSet *)ws_err->data("ds_splot"),
                                           RooArgSet(*ws_err->var("ctau3DErr"), *ws_err->var("n_sig_sw"), *ws_err->var("ctau3DRes"), *ws_err->var("ctau3D"), *ws_err->var("mass")), 0, "n_sig_sw");

    // ===== import sPlot results ===== //
    ws_err->import(*ds_sig_err);
    ws_err->import(*ds_bkg_err);
    ws_err->import(*err_tot_pdf);
    ws_err->import(*err_sig_pdf);
    ws_err->import(*err_bkg_pdf);


    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    // *-* Plotting - very long but simple *-*//
    // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    double minRange = (double)(floor(ctauErrMin * 100.) / 100.);
    double maxRange = (double)(ceil(ctauErrMax * 100.) / 100.);
    ws_err->var("ctau3DErr")->setRange("ctauErrWindow", ctauErrMin, ctauErrMax);
    err_frame1 = ws_err->var("ctau3DErr")->frame(Bins(nBins), Range(minRange - 0.01, maxRange + 0.01)); // modified

    cout << "Check number of events\n";
    cout << "ds_data: " << ws_err->data("ds_data")->numEntries() << "\n";
    cout << "ds_splot: " <<ws_err->data("ds_splot")->numEntries() << "\n\n";

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
    ws_err->data("ds_data")->plotOn(err_frame2, Name("dataCTAUERR_Tot"), MarkerSize(.7), Binning(newBins)); // Normalization(ws_errmc->ds_splot("reducedDS_MC")->sumEntries()
    ws_err->pdf("err_tot_pdf")->plotOn(err_frame2, Name("err_tot_pdf"), LineColor(kGreen + 1), Range(ctauErrMin, ctauErrMax), LineWidth(2), Normalization(ws_err->data("ds_data")->sumEntries(), RooAbsReal::NumEvent));
    ws_err->data("ds_sig_err")->plotOn(err_frame2, Name("dataHist_Sig"), MarkerSize(.7), LineColor(kRed + 2), MarkerColor(kRed + 2), Binning(newBins));
    ws_err->pdf("err_sig_pdf")->plotOn(err_frame2, Name("err_sig_pdf"), LineColor(kRed + 2), LineWidth(2), Range(ctauErrMin, ctauErrMax));
    ws_err->data("ds_bkg_err")->plotOn(err_frame2, Name("dataHist_Bkg"), MarkerSize(.7), LineColor(kBlue + 2), MarkerColor(kBlue + 2), Binning(newBins));
    ws_err->pdf("err_bkg_pdf")->plotOn(err_frame2, Name("err_bkg_pdf"), LineColor(kBlue + 2), LineWidth(2), Range(ctauErrMin, ctauErrMax));

    // set range of y axis 
    Double_t YMax = hTot->GetBinContent(hTot->GetMaximumBin());
    Double_t YMin = 1e99;
    for (int i = 1; i <= hTot->GetNbinsX(); i++)
        if (hTot->GetBinContent(i) > 0)
            YMin = min(YMin, hTot->GetBinContent(i));
    Double_t Yup(0.), Ydown(0.);
    Yup = YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.4 - 0.3)));
    Ydown = YMin / (TMath::Power((YMax / YMin), (0.3 / (1.0 - 0.4 - 0.3))));
    err_frame2->GetYaxis()->SetRangeUser(Ydown, Yup);
    

    // ===== print lost events due to new ctau3DErr range ===== //
    cout << "final ctau3DErr range\n";
    cout << ws_err->var("ctau3DErr")->getMin() << ", " << ws_err->var("ctau3DErr")->getMax() << "\n\n";

    RooDataSet *ctauResCutDS = (RooDataSet *)ds_sig_err->reduce(RooArgSet(*(ws_err->var("ctau3DRes")), *(ws_err->var("ctau3D")), *(ws_err->var("ctau3DErr"))), Form("ctau3DErr>%f&&ctau3DErr<%f", ctauErrMin, ctauErrMax));
    ctauResCutDS->SetName("ctauResCutDS"); // used only for comparison
    ws_err->import(*ctauResCutDS);
    Double_t outTot = ws_err->data("ds_data")->numEntries();
    Double_t outRes = ws_err->data("ctauResCutDS")->numEntries();
    cout << "Check how many events we lost due to new range cut\n";
    cout << "Total evt: " << outTot << "" << endl;
    cout << "Residual evt: " << outRes << "" << endl;
    cout << "lost evt: " << ((outTot - outRes) * 100) / outTot << " %\n\n";

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
    TLegend *leg_err = new TLegend(text_x + 0.5, text_y - 0.2, text_x + 0.7, text_y);
    leg_err->SetTextSize(text_size);
    leg_err->SetTextFont(43);
    leg_err->SetBorderSize(0);
    leg_err->AddEntry(err_frame2->findObject("dataCTAUERR_Tot"), "Data", "pe");
    leg_err->AddEntry(err_frame2->findObject("err_tot_pdf"), "Total PDF", "l");
    leg_err->AddEntry(err_frame2->findObject("err_sig_pdf"), "Signal", "l");
    leg_err->AddEntry(err_frame2->findObject("err_bkg_pdf"), "Background", "l");
    leg_err->Draw("same");
    drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x, text_y, text_color, text_size);
    if (yLow == 0)
        drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x, text_y - y_diff, text_color, text_size);
    else if (yLow != 0)
        drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x, text_y - y_diff, text_color, text_size);
    drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x, text_y - y_diff * 2, text_color, text_size);
    drawText(Form("%.2f < cos#theta_{EP} < %.2f", cos_low, cos_high), text_x, text_y - y_diff * 3, text_color, text_size);
    drawText(Form("Loss: (%.4f%s) %.f evts", (outTot - outRes) * 100 / outTot, "%", outTot - outRes), text_x, text_y - y_diff * 4, text_color, text_size);


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
    RooPlot *pull_frame = ws_err->var("ctau3DErr")->frame(Bins(nBins), Range(minRange - 0.01, maxRange + 0.01));
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
    c_err->SaveAs(("figs/" + out_ss + ".png").c_str());


    // ===== Export results ===== //
    // TFile *outFile = new TFile(Form("roots/2DFit_%s/CtauErr_v2Bins/CtauErrResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_v2_%.2f.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP, v2), "recreate");
    // ds_bkg_err->Write();
    // ds_sig_err->Write();
    // err_tot_pdf->Write();
    // err_bkg_pdf->Write();
    // err_sig_pdf->Write();
    // outFile->Close();

    cout << "\n=================================\n";
    cout << "\n Finish ctau3DErr sPlot work\n";
    cout << "\n=================================\n";
    t->Stop();
    printf("RealTime=%f seconds, CpuTime=%f seconds\n", t->RealTime(), t->CpuTime());
}