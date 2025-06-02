#pragma once

#include "TString.h"
#include "TH1D.h"
#include "TTree.h"
#include "Data.h"
#include "TVector3.h"
#include "TLorentzVector.h"

class Algorithm {
  public:
    /**
     * @brief construct a Alogrithm instance
     */
    Algorithm();

    /**
     * @brief initialize the algorithm
     * @param data - pointer to the instance of Data class
     */
    void initialize(Data* data);

    /**
     * @brief execute logics
     */
    void execute(Long64_t irqq);

    /**
     * @brief pointer to the TH1D class instance
     */
    TH1D* h_runNb_m = nullptr;

    /**
     * @brief pointer to the output TTree class
     */
    TTree* m_myTree = nullptr;

    /**
     * @brief turn on(off) selection cuts
     */
    bool cut_centrality0_180 = false;
    bool cut_jpsiMass = false;
    bool cut_softMuons = false;
    bool cut_vtxProbability = false;
    bool cut_oppositeSign = false;
    bool cut_singleMuonAcc = false;

    /**
     * @brief methods to make the outputs
     */
    void fillTree();

  protected:
    /**
     * @brief apply the selection
     */
    bool passedSelection(Long64_t irqq);
    bool passedSoftMuonCut(Long64_t irqq);
    bool passedMuonAcc2018();
    bool checkMuonAcc2018(double muPt, double muEta);

    /**
     * @brief fill the histograms and TTree
     */
    void standByFilling(Long64_t irqq);
    void fillHists();

    /**
     * @brief variables for polarization
     */
    TVector3 muplHX;
    TVector3 muplCS;

    /**
     * @brief polarization transform functions
     */
    void MuPlusVector_Helicity(const TLorentzVector &QQLV_Lab, const TLorentzVector &MuPlusLV_Lab);
    void MuPlusVector_CollinsSoper(const TLorentzVector &QQLV_Lab);

    /**
     * @brief instance of the Data class
     */
    Data *m_data = nullptr;

    /**
     * @brief set output branches
     */
    void setBranches();

    /**
     * @brief
     */

    /**
     * @brief variables of m_myTree
     */
    int evt;
    int runN;
    int lumi;
    int cBin;
    int nDimu;
    float vz;
    float mass[1000];
    float pt[1000];
    float pt1[1000];
    float pt2[1000];
    float y[1000];
    float phi[1000];
    float phi1[1000];
    float phi2[1000];
    float eta[1000];
    float eta1[1000];
    float eta2[1000];
    float weight0[1000];
    float weight1[1000];
    int recoQQsign[1000];
    float ctau3D[1000];
    float ctau3DErr[1000];
    float ctau3DTrue[1000];
    float ctau3DRes[1000];
    double TnPweight[1000] = {1.};
    float cosHX[1000];
    float phiHX[1000];
    float cosCS[1000];
    float phiCS[1000];
    double weight = 1;

    // // for polarization
    // float cos_theta[1000];
    // float cos_theta1[1000];
    // float cos_cs[1000];
    // float phi_cs[1000];
    // float cos_hx[1000];
    // float phi_hx[1000];
    // float cos_px[1000];
    // float phi_px[1000];
    // float cos_ep[1000];
    // float phi_ep[1000];
    // float cos_lab[1000];
    // float phi_lab[1000];
};