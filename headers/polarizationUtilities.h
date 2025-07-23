#pragma once

#include <string>

class RooWorkspace;
class TPad;
class RooPlot;
class RooFitResult;
class TString;

// ===== from JpsiUtility.h =====
extern int nCPU;

extern float text_x;
extern float text_y;
extern float y_diff;
extern float text_size;
extern int text_color;

extern float massLow, massHigh;
extern int nMassBin;

extern float ctauLow, ctauHigh;
extern double nCtauBins;

extern double ctauErrLow;
extern double ctauErrHigh;
extern int nCtauErrBins;
extern double binWidth;

extern float ctauResLow, ctauResHigh;
extern int nCtauResBins;

extern int nCtauTrueBins;

struct ParticleMass
{
  double JPsi, Psi2S, Y1S, Y2S, Y3S, Z, PiPlus, KaPlus;
};
extern ParticleMass pdgMass;

// functions
void printChi2(RooWorkspace *myws, TPad *Pad, RooPlot *frame, RooFitResult *fitRes, std::string varLabel, std::string dataLabel, std::string pdfLabel, int nBins, bool useDefaultName = true);

//
//
//

// ===== from commonUtility.h =====
void drawText(const char *text, float xp, float yp, int textColor, int textSize = 18);


// ===== from cutsAndBin.h =====
TString getKineLabel(float ptLow, float ptHigh, float yLow, float yHigh, float muPtCut_, int cLow, int cHigh);

TString getKineLabelpp(float ptLow, float ptHigh, float yLow, float yHigh, float muPtCut_);