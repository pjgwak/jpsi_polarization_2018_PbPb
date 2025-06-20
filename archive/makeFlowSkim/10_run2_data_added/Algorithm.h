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
    void initialize(Data* data, bool isRun2);

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
    bool cut_HLTriggerPbPbJpsi2018 = false;
    bool cut_recoQQTrigger = false;
    bool cut_L2L3FilterPbPbJpsi2018 = false;
    bool cut_runNb327123 = false;

    /**
     * @brief methods to make the outputs
     */
    void fillTree();

  protected:
    /**
     * Physics trigger numbers
     */
    int kTrigJpsi2018 = 12;
    int kL2filterJpsi2018 = 16;
    int kL3filterJpsi2018 = 17;

    int kTrigJpsipp2015 = 3;
    

    /**
     * @brief apply the selection
     */
    bool passedSelection(Long64_t irqq);
    bool passedSoftMuonCut(Long64_t irqq);
    bool passedMuonAcc2018();
    bool checkMuonAcc2018(double muPt, double muEta);
    bool passedRunNb327123();
    bool passedHLTriggerJpsiPbPb2018();
    bool passedRecoQQTrigger(Long64_t irqq);
    bool passedHLFilterJpsiPbPb2018(Long64_t irqq);

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
     * and run2 sample flag
     */
    Data *m_data = nullptr;
    bool m_isRun2 = false;

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
};


// // ====== Childern classes ===== //
// class AlgorithmRun2 {
// 	public:
// 		AlgorithmRun2();
// 		~AlgorithmRun2() override = default;

// 		/**
// 		 * @brief initialize the event loop
// 		 */
// 		void initialize()override;
// 		void execute() override;

// 		/**
// 		 * @brief instance of the Data class
// 		 */
// 		DataRun2 *m_data = nullptr;
// };

// class AlgorithmRun3 {
// 	public:
// 		AlgorithmRun3();
// 		~AlgorithmRun3() override = default;

// 		void initialize() override;
// 		void execute() override;

// 		DataRun3 *m_data = nullptr;
// };