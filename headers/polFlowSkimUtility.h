#pragma once
// It collects the necessary values and functiosn from legacy CMS HI headers
// for the process of "Onia to FlowSkim" 


// =======================================
// --- from Raa cutsAndBin.h ---
// =======================================
const int nBins = 200; // table of bin edges

extern const double binTable[nBins+1];
extern const double binTable_up[nBins+1];
extern const double binTable_down[nBins+1];

int getHiBinFromhiHF(const double hiHF);
int getHiBinFromhiHF_Up(const double hiHF);
int getHiBinFromhiHF_Down(const double hiHF);
double findNcoll(int hiBin);


// =======================================
// --- from Raa tnp_weight_lowptPbPb.h ---
// =======================================
#include "TMath.h"

// IN THIS FILE YOU WILL FIND:
// +++++++++++++++++++++++++++++++++++++++
// - Trigger: (tnp_weight_trg_pbpb)
//   * filterId = 0: Jpsi L2 filter
//   * filterId = 1: Jpsi L3 filter
//   * filterId = 2: Upsi L2 filter
//   * filterId = 3: Upsi L3 filter
//   * idx = 0:  nominal
//   * idx = -1: syst variation, +1 sigma
//   * idx = -2: syst variation, -1 sigma
//   * idx = +1: stat variation, +1 sigma
//   * idx = +2: stat variation, -1 sigma
// - MuID: (tnp_weight_muid_pbpb)
//   * idx = 0:  nominal
//   * idx = -1: syst variation, +1 sigma
//   * idx = -2: syst variation, -1 sigma
//   * idx = +1: stat variation, +1 sigma
//   * idx = +2: stat variation, -1 sigma
// - Inner tracking: (tnp_weight_trk_pbpb)
//   * idx = 0:  nominal
//   * idx = -1: syst variation, +1 sigma
//   * idx = -2: syst variation, -1 sigma
//   * idx = +1: stat variation, +1 sigma
//   * idx = +2: stat variation, -1 sigma
// +++++++++++++++++++++++++++++++++++++++

double tnp_weight_muid_pbpb(double pt, double eta, int idx=0);
double tnp_weight_trk_pbpb(double eta, int idx=0);
double tnp_weight_trg_pbpb(double pt, double eta, int filterId=0,int idx=0, bool getRawDen= false);