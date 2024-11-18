#include <ctime>
#include <TLorentzVector.h>
#include "../cms_headers/commonUtility.h"
#include "../cms_headers/HiEvtPlaneList.h"
#include "../cms_headers/cutsAndBin.h"
// #include "../cms_headers/tnp_weight_lowptPbPb.h" -> 차이점??
#include "../cms_headers/tnp_weight_lowptPbPb_num_den_new.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

static const long MAXTREESIZE = 1000000000000;
double getAccWeight(TH1D *h = 0, double pt = 0);
double getEffWeight(TH1D *h = 0, double pt = 0);

// nevt: -1 = Loop all events, isMC: Data or MC, MCtype: 1(Signal) or 2(NP)
// kTrigSel: Apply run2 particle trigger,
// hiHFBinEdge: Heavy ion group Forward Calorimeter systematics info. None(0), Up(1), Down(2),
// PDtype: Central(1) or Peripheral(2)
void onia_to_skim_jpsi_mc(int nevt = -1, bool isMC = true, int MCtype = 1, int kTrigSel = kTrigJpsi, int hiHFBinEdge = 0, int PDtype = 1, bool isPtWeight = true)
{
	using namespace std;
	using namespace hi;

	TString date_label = "241010";

	// Example of using event plane namespace
	cout
		<< " Index of " << EPNames[HFm2] << " = " << HFm2 << endl;
	cout << " Index of " << EPNames[HFp2] << " = " << HFp2 << endl;
	cout << " Index of " << EPNames[trackmid2] << " = " << trackmid2 << endl;

	// Onia input path
	TString fnameDataReReco = "/disk1/Oniatree/miniAOD/OniaTree_miniAOD_HIDoubleMuonPD_addQVect_merged.root";
	TString fnameDataReRecoPeri = "/disk1/Oniatree/miniAOD/OniaTree_miniAOD_HIDoubleMuonPsiPeriPD_addQVect_merged.root";
	TString fnameMC = "/disk1/Oniatree/miniAOD/OniaTreeMC_miniAOD_Jpsi_HydjetMB_5p02TeV_merged.root";
	TString fnameMC_BtoJ = "/disk1/Oniatree/miniAOD/Oniatree_MC_BToJpsi_pThat-2_TuneCP5-EvtGen_HydjetDrumMB_5p02TeV_pythia8_miniAOD.root";

	// pT weighting function inputs
	TFile *fPtW1;
	TFile *fPtW2;
	if (PDtype == 1 || MCtype == 1)
	{
		fPtW1 = new TFile("local_cos_eff_acc_study/pt_weight_inputs/ratioDataMC_AA_Jpsi_DATA_y0_1p6_211201.root", "read");
		fPtW2 = new TFile("local_cos_eff_acc_study/pt_weight_inputs/ratioDataMC_AA_Jpsi_DATA_Forward_y_211218.root", "read");
	}
	else if (PDtype == 2 || MCtype == 2)
	{
		fPtW1 = new TFile("./compareDataToMC/WeightedFcN_fit/ratioDataMC_AA_BtoJpsi_DATA_All_y.root", "read");
		fPtW2 = new TFile("./compareDataToMC/WeightedFcN_fit/ratioDataMC_AA_BtoJpsi_DATA_Forward_y.root", "read");
	}
	TF1 *fptw1 = (TF1 *)fPtW1->Get("dataMC_Ratio1");
	TF1 *fptw2 = (TF1 *)fPtW2->Get("dataMC_Ratio1");

	// Label of output
	TString fPD;
	if (PDtype == 1)
		fPD = "DB";
	else if (PDtype == 2)
		fPD = "DBPeri";

	TString fMCtype;
	if (MCtype == 1)
		fMCtype = "Signal";
	else if (MCtype == 2)
		fMCtype = "NPOnly";

	// Read tree inside root file by using TChain
	TChain *mytree = new TChain("hionia/myTree");
	if (!isMC)
	{
		if (PDtype == 1)
			mytree->Add(fnameDataReReco.Data());
		else if (PDtype == 2)
			mytree->Add(fnameDataReRecoPeri.Data());
	}
	else if (isMC)
	{
		if (MCtype == 1)
			mytree->Add(fnameMC.Data());
		else
			mytree->Add(fnameMC_BtoJ.Data());
	}

	// ROOT file has limit for # of branch
	const int maxBranchSize = 500;

	// Prepare variables
	UInt_t runNb;
	UInt_t eventNb, LS;
	float zVtx;
	Int_t Centrality;
	ULong64_t HLTriggers;
	Float_t SumET_HF;
	Short_t Reco_QQ_size;
	Short_t Reco_mu_size;
	TClonesArray *Reco_QQ_4mom = nullptr;
	TClonesArray *Reco_mu_4mom = nullptr;
	ULong64_t Reco_QQ_trig[maxBranchSize];	//[Reco_QQ_size]
	ULong64_t Reco_mu_trig[maxBranchSize];	//[Reco_QQ_size]
	Float_t Reco_QQ_VtxProb[maxBranchSize]; //[Reco_QQ_size]
	TBranch *b_runNb;
	TBranch *b_eventNb;
	TBranch *b_LS;
	TBranch *b_zVtx;
	TBranch *b_Centrality;
	TBranch *b_HLTriggers;
	TBranch *b_SumET_HF;
	TBranch *b_Reco_QQ_size;
	TBranch *b_Reco_mu_size;
	TBranch *b_Reco_QQ_4mom;
	TBranch *b_Reco_mu_4mom;
	TBranch *b_Reco_QQ_trig;
	TBranch *b_Reco_mu_trig;
	TBranch *b_Reco_QQ_VtxProb;
	Bool_t Reco_mu_highPurity[maxBranchSize]; //[Reco_QQ_size]
	TBranch *b_Reco_mu_highPurity;

	// Connect to branch
	if (!isMC)
	{
		mytree->SetBranchAddress("runNb", &runNb, &b_runNb);
		mytree->SetBranchAddress("LS", &LS, &b_LS);
	}

	mytree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
	mytree->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
	mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
	mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
	mytree->SetBranchAddress("SumET_HF", &SumET_HF, &b_SumET_HF);
	mytree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
	mytree->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
	mytree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
	mytree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
	mytree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
	mytree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
	mytree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
	mytree->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);

	//  muon id
	Short_t Reco_QQ_mupl_idx[maxBranchSize];
	Short_t Reco_QQ_mumi_idx[maxBranchSize];
	TBranch *b_Reco_QQ_mupl_idx;
	TBranch *b_Reco_QQ_mumi_idx;
	mytree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
	mytree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);

	Int_t Reco_mu_nTrkHits[maxBranchSize]; //[Reco_mu_size]
	TBranch *b_Reco_mu_nTrkHits;
	mytree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
	Float_t Reco_mu_normChi2_global[maxBranchSize]; //[Reco_mu_size]
	TBranch *b_Reco_mu_normChi2_global;
	mytree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
	Int_t Reco_mu_nMuValHits[maxBranchSize]; //[Reco_mu_size]
	TBranch *b_Reco_mu_nMuValHits;
	mytree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
	Int_t Reco_mu_StationsMatched[maxBranchSize]; //[Reco_mu_size]
	TBranch *b_Reco_mu_StationsMatched;
	mytree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
	Float_t Reco_mu_dxy[maxBranchSize];	   //[Reco_mu_size]
	Float_t Reco_mu_dxyErr[maxBranchSize]; //[Reco_mu_size]
	TBranch *b_Reco_mu_dxy;
	TBranch *b_Reco_mu_dxyErr;
	mytree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
	mytree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
	Float_t Reco_mu_dz[maxBranchSize];	  //[Reco_mu_size]
	Float_t Reco_mu_dzErr[maxBranchSize]; //[Reco_mu_size]
	TBranch *b_Reco_mu_dz;
	TBranch *b_Reco_mu_dzErr;
	mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
	mytree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
	Int_t Reco_mu_nTrkWMea[maxBranchSize]; //[Reco_mu_size]
	TBranch *b_Reco_mu_nTrkWMea;
	mytree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
	Bool_t Reco_mu_TMOneStaTight[maxBranchSize]; //[Reco_mu_size]
	TBranch *b_Reco_mu_TMOneStaTight;
	if (!isMC)
	{
		mytree->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
	}

	Int_t Reco_mu_nPixWMea[maxBranchSize]; //[Reco_mu_size]
	TBranch *b_Reco_mu_nPixWMea;
	mytree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
	Short_t Reco_QQ_sign[maxBranchSize]; //[Reco_QQ_size]
	TBranch *b_Reco_QQ_sign;
	mytree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
	Float_t rpAng[29]; //[nEP]
	TBranch *b_rpAng;
	// mytree->SetBranchAddress("rpAng", rpAng, &b_rpAng);

	Int_t Reco_mu_nPixValHits[maxBranchSize]; //[Reco_QQ_size]
	TBranch *b_Reco_mu_nPixValHits;
	mytree->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
	// Float_t         Reco_mu_ptErr_global[maxBranchSize];   //[Reco_QQ_size]
	// TBranch        *b_Reco_mu_ptErr_global;
	// mytree->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);

	Int_t Reco_mu_SelectionType[maxBranchSize];
	TBranch *b_Reco_mu_SelectionType;
	mytree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);

	Float_t Reco_QQ_ctau3D[maxBranchSize];
	Float_t Reco_QQ_ctauErr3D[maxBranchSize];
	TBranch *b_Reco_QQ_ctau3D;
	TBranch *b_Reco_QQ_ctauErr3D;
	mytree->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
	mytree->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);

	Short_t Reco_mu_whichGen[maxBranchSize];
	TBranch *b_Reco_mu_whichGen;
	Float_t Gen_weight;
	TBranch *b_Gen_weight;

	if (isMC)
	{
		mytree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen, &b_Reco_mu_whichGen);
		mytree->SetBranchAddress("Gen_weight", &Gen_weight, &b_Gen_weight);
	}

	// Particle tirgger selection
	// Refer to "Eff_Acc/cutsAndBin.h:int kTrigJpsi = 12;"
	int trigIndx = 0;
	if (kTrigSel == kTrigJpsi)
		trigIndx = 0;
	else if (kTrigSel == kTrigUps)
		trigIndx = 1;
	else if (kTrigSel == kTrigL1DBOS40100)
		trigIndx = 2;
	else if (kTrigSel == kTrigL1DB50100)
		trigIndx = 3;

	// Initialize variables for tnp weighting
	double tnp_weight = 1;
	double tnp_trig_weight = 1;
	double tnp_trig_weight_muplL2_num = -1;
	double tnp_trig_weight_muplL3_num = -1;
	double tnp_trig_weight_mumiL2_num = -1;
	double tnp_trig_weight_mumiL3_num = -1;
	double tnp_trig_weight_muplL2_den = -1;
	double tnp_trig_weight_muplL3_den = -1;
	double tnp_trig_weight_mumiL2_den = -1;
	double tnp_trig_weight_mumiL3_den = -1;
	double tnp_trig_weight_num = 1;
	double tnp_trig_weight_den = 1;
	double tnp_trig_weight_mupl = -1;
	double tnp_trig_weight_mumi = -1;
	double tnp_trig_dimu = -1;

	// L2, L3 trigger - Don't touch for Jpsi
	int kL2filter = 16;
	int kL3filter = 17;

	// To check number of muons
	int count = 0;
	int counttnp = 0;

	// HF up, down
	TString fCentSelHF = "HFNom";
	if (hiHFBinEdge == 1)
		fCentSelHF = "HFUp";
	else if (hiHFBinEdge == -1)
		fCentSelHF = "HFDown";
	TFile *out_file;
	if (isMC)
	{
		out_file = new TFile(Form("skimmed_files/OniaFlowSkim_%s_miniAOD_isMC%d_%s_%s_%s.root", fTrigName[trigIndx].Data(), isMC, fMCtype.Data(), fCentSelHF.Data(), date_label.Data()), "recreate");
	}
	else
	{
		out_file = new TFile(Form("skimmed_files/OniaFlowSkim_%sTrig_%sPD_miniAOD_isMC%d_%s_%s.root", fTrigName[trigIndx].Data(), fPD.Data(), isMC, fCentSelHF.Data(), date_label.Data()), "recreate");
	}

	// Variables for output
	const static int nMaxDimu = 1000;
	int evt;
	int runN;
	int lumi;
	int cBin;
	int nDimu;
	float vz;
	float mass[nMaxDimu];
	float pt[nMaxDimu];
	float y[nMaxDimu];
	float phi[nMaxDimu];
	float eta[nMaxDimu];
	float eta1[nMaxDimu];
	float eta2[nMaxDimu];
	float phi1[nMaxDimu];
	float phi2[nMaxDimu];
	float pt1[nMaxDimu];
	float pt2[nMaxDimu];
	float weight0[nMaxDimu];
	float weight1[nMaxDimu];
	int recoQQsign[nMaxDimu];
	float ctau3D[nMaxDimu];
	float ctau3DErr[nMaxDimu];
	float ctau3DRes[nMaxDimu];
	float ctau3D2S[nMaxDimu];
	float ctau3DErr2S[nMaxDimu];
	float ctau3DRes2S[nMaxDimu];
	double TnPweight[nMaxDimu] = {1.};
	double weight[nMaxDimu] = {1.};
	double pt_weight = 1;
	double cos_theta_hx[nMaxDimu];

	// Connec to output tree branch
	// muon+muon pair event tree
	TTree *mmevttree = new TTree("mmepevt", "dimuonAndEventPlanes in event based");
	mmevttree->SetMaxTreeSize(MAXTREESIZE);
	mmevttree->Branch("event", &evt, "event/I");
	mmevttree->Branch("runN", &runN, "runN/I");
	mmevttree->Branch("lumi", &lumi, "lumi/I");
	mmevttree->Branch("cBin", &cBin, "cBin/I");
	mmevttree->Branch("vz", &vz, "vz/F");
	mmevttree->Branch("nDimu", &nDimu, "nDimu/I");
	mmevttree->Branch("mass", mass, "mass[nDimu]/F");
	mmevttree->Branch("y", y, "y[nDimu]/F");
	mmevttree->Branch("pt", pt, "pt[nDimu]/F");
	mmevttree->Branch("pt1", pt1, "pt1[nDimu]/F");
	mmevttree->Branch("pt2", pt2, "pt2[nDimu]/F");
	mmevttree->Branch("eta", eta, "eta[nDimu]/F");
	mmevttree->Branch("eta1", eta1, "eta1[nDimu]/F");
	mmevttree->Branch("eta2", eta2, "eta2[nDimu]/F");
	mmevttree->Branch("recoQQsign", recoQQsign, "recoQQsign[nDimu]/I");
	mmevttree->Branch("ctau3D", ctau3D, "ctau3D[nDimu]/F");
	mmevttree->Branch("ctau3DErr", ctau3DErr, "ctau3DErr[nDimu]/F");
	mmevttree->Branch("ctau3DRes", ctau3DRes, "ctau3DRes[nDimu]/F");
	mmevttree->Branch("ctau3D2S", ctau3D2S, "ctau3D2S[nDimu]/F");
	mmevttree->Branch("ctau3DErr2S", ctau3DErr2S, "ctau3DErr2S[nDimu]/F");
	mmevttree->Branch("ctau3DRes2S", ctau3DRes2S, "ctau3DRes2S[nDimu]/F");
	mmevttree->Branch("weight", &weight, "weight[nDimu]/D");
	mmevttree->Branch("TnPweight", TnPweight, "TnPweight[nDimu]/D");
	mmevttree->Branch("cos_theta_hx", cos_theta_hx, "cos_theta_hx[nDimu]/D");
	////////////////////////////////////////////////////////////////////////
	////////////////// TLorentzVector dummies
	////////////////////////////////////////////////////////////////////////
	TLorentzVector *JP_Reco = new TLorentzVector;
	TLorentzVector *mupl_Reco = new TLorentzVector;
	TLorentzVector *mumi_Reco = new TLorentzVector;

	int counts[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	// Event loop start
	if (nevt == -1)
		nevt = mytree->GetEntries();

	cout << "Total events = " << nevt << endl;
	double weight_ = 1;
	for (int iev = 0; iev < nevt; ++iev)
	{
		// if (iev == 1627597 || iev == 5734759 || iev == 6916044 || iev == 11512879 || iev == 14288366)
		//{
		//	continue;
		// }
		//  cout << iev << endl;
		//   Progress announcement
		if (iev % 100000 == 0)
		{
			cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() << " (" << (int)(100. * iev / mytree->GetEntries()) << "%)" << endl;
			cout << "# of dimuon pair: " << count << endl;
		}

		// Get values from input tree
		mytree->GetEntry(iev);

		nDimu = 0;

		// Need for PD file - PD file has both PD and Peri runs. Only run < 327123: central collision
		if (!isMC && PDtype == 1 && runNb >= 327123)
			continue;

		// weight 1 means no weight - usually data events have no weight.
		weight_ = 1.;
		if (isMC)
		{
			weight_ = findNcoll(Centrality) * Gen_weight;
		}
		counts[0]++;

		// HLT trigger - Must match to Jpsi
		if (!((HLTriggers & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel))))
			continue;

		counts[1]++;
		counts[2] += Reco_QQ_size;

		// cout << "Reco_QQ_size : " << Reco_QQ_size << endl;
		for (Int_t irqq = 0; irqq < Reco_QQ_size; ++irqq)
		{

			// Loop inside the event. One event can have reconsdtruced particles.
			// Is this particle Jpsi satisfying cuts?
			// irqq: Just index for interation. Maybe named from "iterator of qq vector"?
			// cout << "irqq : " << irqq << endl;
			runN = runNb;
			evt = eventNb;
			lumi = LS;
			cBin = -999;
			if (hiHFBinEdge == 0)
				cBin = getHiBinFromhiHF(SumET_HF);
			else if (hiHFBinEdge == 1)
				cBin = getHiBinFromhiHF_Up(SumET_HF);
			else if (hiHFBinEdge == -1)
				cBin = getHiBinFromhiHF_Down(SumET_HF);
			if (cBin == -999)
			{
				cout << "ERROR!!! No HF Centrality Matching!!" << endl;
				return;
			}
			vz = zVtx;

			JP_Reco = (TLorentzVector *)Reco_QQ_4mom->At(irqq);
			mupl_Reco = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
			mumi_Reco = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);

			// Acceptance cut
			bool muplAcc = (((TMath::Abs(mupl_Reco->Eta()) <= 1.2) && (mupl_Reco->Pt() >= 3.5)) ||
							((TMath::Abs(mupl_Reco->Eta()) > 1.2) && (TMath::Abs(mupl_Reco->Eta()) <= 2.1) && (mupl_Reco->Pt() >= 5.47 - 1.89 * (TMath::Abs(mupl_Reco->Eta())))) ||
							((TMath::Abs(mupl_Reco->Eta()) > 2.1) && (TMath::Abs(mupl_Reco->Eta()) <= 2.4) && (mupl_Reco->Pt() >= 1.5)));
			bool mumiAcc = (((TMath::Abs(mumi_Reco->Eta()) <= 1.2) && (mumi_Reco->Pt() >= 3.5)) ||
							((TMath::Abs(mumi_Reco->Eta()) > 1.2) && (TMath::Abs(mumi_Reco->Eta()) <= 2.1) && (mumi_Reco->Pt() >= 5.47 - 1.89 * (TMath::Abs(mumi_Reco->Eta())))) ||
							((TMath::Abs(mumi_Reco->Eta()) > 2.1) && (TMath::Abs(mumi_Reco->Eta()) <= 2.4) && (mumi_Reco->Pt() >= 1.5)));
			counts[3]++;
			if (!(muplAcc && mumiAcc))
				continue;

			// Particle musth match to Jpsi trigger
			if (!((Reco_QQ_trig[irqq] & ((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel))))
				continue;

			if (isMC)
			{
				if (Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1)
					continue;
				if (Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1)
					continue;
			}

			bool passMuonTypePl = true;
			passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]] & ((int)pow(2, 1)));
			passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]] & ((int)pow(2, 3)));

			bool passMuonTypeMi = true;
			passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]] & ((int)pow(2, 1)));
			passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]] & ((int)pow(2, 3)));

			// Soft muon cut for positive muon
			bool muplSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]]==true) &&
				(Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]] > 5) &&
				(Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]] > 0) &&
				(fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]]) < 0.3) &&
				(fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]]) < 20.) &&
				passMuonTypePl //                       &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true)
			);

			// Soft muon cut for negative muon
			bool mumiSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]]==true) &&
				(Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5) &&
				(Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0) &&
				(fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]]) < 0.3) &&
				(fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]]) < 20.) &&
				passMuonTypeMi //                        &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true)
			);

			counts[4]++;
			// muons must pass soft muon cut
			if (!(muplSoft && mumiSoft))
				continue;

			counts[5]++;
			// Vertext chi^2 probability must > 1 %
			// Roughly speaking, it's he probability muons came from vertex
			if (Reco_QQ_VtxProb[irqq] < 0.01)
				continue;

			// Muon pair must have opposite sign -> This cut will be applied in skim to RooDataset step.
			recoQQsign[irqq] = Reco_QQ_sign[irqq];
			count++;

			// Calculate TnP weighting -> Used to get acc and eff weightings in other codes
			if (isMC)
			{
				tnp_weight = 1;
				tnp_trig_weight = 1;
				tnp_trig_weight_mupl = -1;
				tnp_trig_weight_mumi = -1;
				tnp_trig_weight_muplL2_num = -1;
				tnp_trig_weight_muplL3_num = -1;
				tnp_trig_weight_mumiL2_num = -1;
				tnp_trig_weight_mumiL3_num = -1;
				tnp_trig_weight_muplL2_den = -1;
				tnp_trig_weight_muplL3_den = -1;
				tnp_trig_weight_mumiL2_den = -1;
				tnp_trig_weight_mumiL3_den = -1;
				tnp_trig_weight_num = 1;
				tnp_trig_weight_den = 1;

				tnp_weight = tnp_weight * tnp_weight_muid_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 0) * tnp_weight_muid_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 0); // mu id
				tnp_weight = tnp_weight * tnp_weight_trk_pbpb(mupl_Reco->Eta(), 0) * tnp_weight_trk_pbpb(mumi_Reco->Eta(), 0);									   // inner tracker

				counts[6]++;
				// Trigger part
				if (!((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) && (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))))
				{
					//         cout << "irqq : " << irqq << " - iev : " << iev << endl;
					//        cout << "TnP ERROR !!!! ::: No matched L2 filter1 " << endl;
					continue;
				}
				bool mupl_L2Filter = ((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))) ? true : false;
				bool mupl_L3Filter = ((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter))) ? true : false;
				bool mumi_L2Filter = ((Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))) ? true : false;
				bool mumi_L3Filter = ((Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter))) ? true : false;
				if (mupl_L2Filter == false || mumi_L2Filter == false)
				{
					cout << "TnP ERROR !!!! ::: No matched L2 filter2 " << endl;
					cout << endl;
				}

				bool mupl_isL2 = (mupl_L2Filter && !mupl_L3Filter) ? true : false;
				bool mupl_isL3 = (mupl_L2Filter && mupl_L3Filter) ? true : false;
				bool mumi_isL2 = (mumi_L2Filter && !mumi_L3Filter) ? true : false;
				bool mumi_isL3 = (mumi_L2Filter && mumi_L3Filter) ? true : false;
				bool SelDone = false;

				if (mupl_isL2 && mumi_isL3)
				{
					tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 0, 0);
					tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 1, 0);
					SelDone = true;
					tnp_trig_weight = tnp_trig_weight_mupl * tnp_trig_weight_mumi;
				}
				else if (mupl_isL3 && mumi_isL2)
				{
					tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 1, 0);
					tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 0, 0);
					SelDone = true;
					tnp_trig_weight = tnp_trig_weight_mupl * tnp_trig_weight_mumi;
				}
				else if (mupl_isL3 && mumi_isL3)
				{

					tnp_trig_weight_muplL2_num = tnp_weight_trg_pbpb_num(mupl_Reco->Pt(), mupl_Reco->Eta(), 0, 0);
					tnp_trig_weight_muplL3_num = tnp_weight_trg_pbpb_num(mupl_Reco->Pt(), mupl_Reco->Eta(), 1, 0);
					tnp_trig_weight_mumiL2_num = tnp_weight_trg_pbpb_num(mumi_Reco->Pt(), mumi_Reco->Eta(), 0, 0);
					tnp_trig_weight_mumiL3_num = tnp_weight_trg_pbpb_num(mumi_Reco->Pt(), mumi_Reco->Eta(), 1, 0);

					tnp_trig_weight_muplL2_den = tnp_weight_trg_pbpb_den(mupl_Reco->Pt(), mupl_Reco->Eta(), 0, 0);
					tnp_trig_weight_muplL3_den = tnp_weight_trg_pbpb_den(mupl_Reco->Pt(), mupl_Reco->Eta(), 1, 0);
					tnp_trig_weight_mumiL2_den = tnp_weight_trg_pbpb_den(mumi_Reco->Pt(), mumi_Reco->Eta(), 0, 0);
					tnp_trig_weight_mumiL3_den = tnp_weight_trg_pbpb_den(mumi_Reco->Pt(), mumi_Reco->Eta(), 1, 0);

					tnp_trig_weight_num = tnp_trig_weight_muplL2_num * tnp_trig_weight_mumiL3_num + tnp_trig_weight_mumiL2_num * tnp_trig_weight_muplL3_num - tnp_trig_weight_muplL3_num * tnp_trig_weight_mumiL3_num;
					tnp_trig_weight_den = tnp_trig_weight_muplL2_den * tnp_trig_weight_mumiL3_den + tnp_trig_weight_mumiL2_den * tnp_trig_weight_muplL3_den - tnp_trig_weight_muplL3_den * tnp_trig_weight_mumiL3_den;
					tnp_trig_weight = tnp_trig_weight_num / tnp_trig_weight_den;
				}
				tnp_weight = tnp_weight * tnp_trig_weight;
				counttnp++;
			}
			// End of selection cuts

			// Lorentz transfomration
			// Define 4-momentum in Lab frame
			TLorentzVector jpsi_lab(JP_Reco->Px(), JP_Reco->Py(),
									JP_Reco->Pz(), JP_Reco->Energy());
			TLorentzVector mupl_lab(mupl_Reco->Px(), mupl_Reco->Py(),
									mupl_Reco->Pz(), mupl_Reco->Energy());
			TLorentzVector mumi_lab(mumi_Reco->Px(), mumi_Reco->Py(),
									mumi_Reco->Pz(), mumi_Reco->Energy());

			// Copy 4-momentum
			TLorentzVector *mupl_QRF = new TLorentzVector(mupl_lab);
			TLorentzVector *mumi_QRF = new TLorentzVector(mumi_lab);

			// Boost to Quarkonia rest frame
			mupl_QRF->Boost(-jpsi_lab.BoostVector());
			mumi_QRF->Boost(-jpsi_lab.BoostVector());

			pt_weight = 1;
			if (isPtWeight && fabs(JP_Reco->Rapidity()) < 1.6)
				pt_weight = fptw1->Eval(JP_Reco->Pt());
			if (isPtWeight && fabs(JP_Reco->Rapidity()) > 1.6)
				pt_weight = fptw2->Eval(JP_Reco->Pt());

			// Push the values into the branches
			if (isMC)
			{
				TnPweight[nDimu] = tnp_weight;
				weight[nDimu] = weight_ * tnp_weight * pt_weight;
			}

			mass[nDimu] = JP_Reco->M();
			phi[nDimu] = JP_Reco->Phi();
			phi1[nDimu] = mupl_Reco->Phi();
			phi2[nDimu] = mumi_Reco->Phi();
			eta[nDimu] = JP_Reco->Eta();
			y[nDimu] = JP_Reco->Rapidity();
			pt[nDimu] = JP_Reco->Pt();
			pt1[nDimu] = mupl_Reco->Pt();
			pt2[nDimu] = mumi_Reco->Pt();
			eta1[nDimu] = mupl_Reco->Eta();
			eta2[nDimu] = mumi_Reco->Eta();
			ctau3D[nDimu] = Reco_QQ_ctau3D[irqq];
			ctau3DErr[nDimu] = Reco_QQ_ctauErr3D[irqq];
			ctau3DRes[nDimu] = (Reco_QQ_ctau3D[irqq]) / (Reco_QQ_ctauErr3D[irqq]);
			ctau3D2S[nDimu] = ctau3D[nDimu] * (pdgMass.Psi2S / pdgMass.JPsi);
			ctau3DErr2S[nDimu] = ctau3DErr[nDimu] * (pdgMass.Psi2S / pdgMass.JPsi);
			ctau3DRes2S[nDimu] = ctau3DRes[nDimu] * (pdgMass.Psi2S / pdgMass.JPsi);
			cos_theta_hx[irqq] = TMath::Cos(jpsi_lab.Angle(mupl_QRF->Vect()));
			nDimu++;

			// Delocate memories of object
			mupl_QRF->Delete();
			mumi_QRF->Delete();
		} // end of dimuon loop

		if (nDimu > 0)
			mmevttree->Fill();

	} // end of event loop
	cout << "count " << count << endl;
	cout << "counttnp " << counttnp << endl;
	cout << "0 : " << counts[0] << endl;
	cout << "1 : " << counts[1] << endl;
	cout << "2 : " << counts[2] << endl;
	cout << "3 : " << counts[3] << endl;
	cout << "4 : " << counts[4] << endl;
	cout << "5 : " << counts[5] << endl;
	cout << "6 : " << counts[6] << endl;

	out_file->cd();
	mmevttree->Write();
	out_file->Close();
}
