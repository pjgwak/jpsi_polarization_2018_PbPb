# include "polarizationUtilities.h"

#include <cmath>

#include "RooWorkspace.h"
#include "TPad.h"
#include "RooPlot.h"
#include "TString.h"
#include "RooFitResult.h"
#include "RooHist.h"
#include "RooRealVar.h"
#include "TLatex.h"
#include "TColor.h"

// ===== from JpsiUtility.h =====
// external variables
int nCPU = 24;

float text_x = 0.15, text_y = 0.816;
float y_diff = 0.05; // 0.05
float text_size = 13; // 13
int text_color = 1;

float massLow = 2.6, massHigh = 3.5;
int nMassBin = 36;

float ctauLow = -4.0, ctauHigh = 7.0;
double nCtauBins = 200; // ctau3D

double ctauErrLow = 1e-6, ctauErrHigh = 0.22;
int nCtauErrBins = 144;
double binWidth = 0.0025;

float ctauResLow = -10, ctauResHigh = 10;
int nCtauResBins = 72;

int nCtauTrueBins = 50;

ParticleMass pdgMass = {3.096, 3.686, 9.460, 10.023, 10.355, 91.188, 0.139570, 0.49367};

void printChi2(RooWorkspace *myws, TPad *Pad, RooPlot *frame, RooFitResult *fitRes, std::string varLabel, std::string dataLabel, std::string pdfLabel, int nBins, bool useDefaultName)
{
  // updated: filter only non-zero bins for calcualtion to avoid inf numerator - pjgwak
  // Use: hpull->GetN() instead of nBins
  // Todo: So the parameter nBins doesn't do anything! but I leave it to avoid chaing every code...

  double chi2 = 0;
  unsigned int ndof = 0;
  Pad->cd();
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextSize(0.1);
  // unsigned int nFitPar = myws->pdf(pdfLabel.c_str())->getParameters(*myws->data(dataLabel.c_str()))->selectByAttrib("Constant",kFALSE)->getSize();
  unsigned int nFitPar = fitRes->floatParsFinal().getSize();
  RooHist *hpull = frame->pullHist(dataLabel.c_str(), pdfLabel.c_str(), true);
  double *ypulls = hpull->GetY();
  unsigned int nFullBins = 0;
  for (int i = 0; i < hpull->GetN(); i++)
  {
    if (ypulls[i] == 0)
      continue;
    chi2 += ypulls[i] * ypulls[i];
    nFullBins++;
  }
  ndof = nFullBins - nFitPar;
  t->DrawLatex(0.7, 0.85, Form("#chi^{2}/ndof = %.0f / %d ", chi2, ndof));

  if (useDefaultName)
  {
    RooRealVar chi2Var("chi2", "chi2", chi2);
    RooRealVar ndofVar("ndof", "ndof", ndof);
    myws->import(chi2Var);
    myws->import(ndofVar);
  }
  else
  {
    RooRealVar chi2Var((std::string("chi2_") + pdfLabel).c_str(), (std::string("chi2_") + pdfLabel).c_str(), chi2);
    RooRealVar ndofVar((std::string("ndof_") + pdfLabel).c_str(), (std::string("ndof_") + pdfLabel).c_str(), ndof);
    myws->import(chi2Var);
    myws->import(ndofVar);
  }
  delete hpull;
};

// ===== from commonUtility.h =====
void drawText(const char *text, float xp, float yp, int textColor = kBlack, int textSize)
{

  TLatex *tex = new TLatex(xp, yp, text);
  tex->SetTextFont(43);
  //   if(bold)tex->SetTextFont(43);
  tex->SetTextSize(textSize);
  tex->SetTextColor(textColor);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}


// ===== from cutsAndBin.h =====
TString getKineLabel(float ptLow, float ptHigh, float yLow, float yHigh, float muPtCut_, int cLow, int cHigh)
{
  TString kineLabel = Form("pt%.1f-%.1f_y%.1f-%.1f_muPt%.1f", ptLow, ptHigh, yLow, yHigh, (float)muPtCut_);
  kineLabel = kineLabel + Form("_centrality%d-%d", (int)cLow, (int)cHigh);
  return kineLabel;
}

TString getKineLabelpp(float ptLow, float ptHigh, float yLow, float yHigh, float muPtCut_)
{
  TString kineLabel = Form("pt%.1f-%.1f_y%.1f-%.1f_muPt%.1f", ptLow, ptHigh, yLow, yHigh, (float)muPtCut_);
  return kineLabel;
}