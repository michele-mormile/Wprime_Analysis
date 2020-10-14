#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include <vector>
#include <assert.h>
#include <TMVA/Reader.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <sstream>
#include <string>
#include "TFileCollection.h"
#include "THashList.h"
#include "TBenchmark.h"
#include "TF1.h"
//#include "/afs/cern.ch/user/f/fcarneva/CMSSW_10_2_0/src/Tprime/TprimeAnalysis/src/DMTopVariables.h"
#include "/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/interface/DMTopVariables.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/interface/Weights.h"
//#include "Weights.h"
#include "/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/interface/MT2Utility.h"
#include "/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/interface/mt2w_bisect.h"
#include "/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/interface/mt2bl_bisect.h"
#include "/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/interface/Mt2Com_bisect.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "/home/mmormile/CMSSW_9_4_9/src/Wprime/WprimeAnalysis/interface/topTagging.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

using namespace std;

bool sortByVal(const pair<int, float> &a,
               const pair<int, float> &b)
{
    return (a.second > b.second);
}

int main(int argc, char **argv)
{
    string path(argv[1]); // list of root file

    string outFile(argv[2]);

    TFileCollection fc("FileCollection", "FileCollection", path.c_str());

    TString treePath = "Events"; //QUI PATH

    TChain chain(treePath);
    chain.AddFileInfoList(fc.GetList());

    //----------------Devo vedere quali di queste variabili servono davvero-------------//

    int sizeMax = 30;
    int PFCandsSizeMax = 500;
    Int_t jetSize(0.), fatjetSize(0.), muonsTight_size(0.), muonsLoose_size = (0.), muons_size(0.), genP_size(0.);

    float nLepTop(0.);

    float genP_Mass[sizeMax];
    float genP_Pt[sizeMax];
    float genP_Phi[sizeMax];
    float genP_Eta[sizeMax];
    float genP_flavor[sizeMax];
    Float_t genP_mom1I[sizeMax], genP_mom2I[sizeMax];

    float jetmass[sizeMax];
    float jetpt[sizeMax];
    float jetphi[sizeMax];
    float jeteta[sizeMax];
    float jetcsv[sizeMax];
    float jetDeepCSV[sizeMax];
    float jetisDeepCSVL[sizeMax];
    float jetisDeepCSVM[sizeMax];
    float jetisDeepCSVT[sizeMax];
    float jetId[sizeMax];

    float jethadronAK8flavour[sizeMax];
    float jetpartonflavour[sizeMax];

    int nJetPFCands;
    int Jet_nConstituents[sizeMax];
    int JetPFCands_candIdx[PFCandsSizeMax];
    int JetPFCands_jetIdx[PFCandsSizeMax];
    float JetPFCands_mass[PFCandsSizeMax];
    float JetPFCands_pt[PFCandsSizeMax];
    float JetPFCands_eta[PFCandsSizeMax];
    float JetPFCands_phi[PFCandsSizeMax];
    float JetPFCands_d0[PFCandsSizeMax];
    float JetPFCands_dz[PFCandsSizeMax];
    float JetPFCands_pdgId[PFCandsSizeMax];

    float fatjetmass[sizeMax];
    float fatjetpt[sizeMax];
    float fatjeteta[sizeMax];
    float fatjetphi[sizeMax];
    float SDmass[sizeMax];
    float Prunedmass[sizeMax];
    float fatjettopmass[sizeMax];

    float fatjettau3[sizeMax];
    float fatjettau2[sizeMax];
    float fatjettau1[sizeMax];
    float fatjetsubjets_Ncsv[sizeMax];
    float fatjetsubjets_Ncsv_tm[sizeMax];
    float fatjet_nDeepCSVsubj[sizeMax];
    float fatjetsubjetsN[sizeMax];
    float fatjettau3OVER2[sizeMax], fatjettau2OVER1[sizeMax];

    float metptCorr[1];
    float metphiCorr[1];
    float muons_IsGlobalMuon[sizeMax], muons_IsTrackerMuon[sizeMax], muons_NumberMatchedStations[sizeMax], muons_NumberOfValidTrackerHits[sizeMax];

    float nLooseMuons[1], nLooseElectrons[1], nTightMuons[1], nTightElectrons[1], nMuons[1];
    float muonsLoose_Pt[sizeMax], muonsTight_Pt[sizeMax], muonsLoose_Phi[sizeMax], muonsTight_Phi[sizeMax], muonsTight_Eta[sizeMax], muons_DB[sizeMax], muons_Dxyerr[sizeMax], muons_Dz[sizeMax], muons_Dxy[sizeMax], muons_MiniIso[sizeMax];
    float muonsCharge[sizeMax], muon_IsHighPtMuon[sizeMax], muon_IsTightPtMuon[sizeMax], muon_IsLoosePtMuon[sizeMax];
    float muons_Pt[sizeMax], muons_Phi[sizeMax], muons_Eta[sizeMax], muons_Mass[sizeMax], muons_Iso04[sizeMax];

    vector<float> muonsHigh_Pt, muonsHigh_E, muonsHigh_Phi, muonsHigh_Eta, elHEEP_Pt, elHEEP_Phi, elHEEP_Eta, elHEEP_E, Lep_E, Lep_Pt, Lep_Phi, Lep_Eta;

    //-------------------------Devo vedere quali branch servono davvero --------------------//

    // chain.SetBranchAddress("Event_nLepTop", &nLepTop);

    chain.SetBranchAddress("MET_pt", metptCorr);
    chain.SetBranchAddress("MET_phi", metphiCorr);

    chain.SetBranchAddress("nGenPart", &genP_size);
    chain.SetBranchAddress("GenPart_pt", genP_Pt);
    chain.SetBranchAddress("GenPart_phi", genP_Phi);
    chain.SetBranchAddress("GenPart_eta", genP_Eta);
    chain.SetBranchAddress("GenPart_mass", genP_Mass);
    chain.SetBranchAddress("GenPart_pdgId", genP_flavor);
    chain.SetBranchAddress("GenPart_genPartIdxMother", genP_mom1I);

    //Jet AK4
    chain.SetBranchAddress("Jet_mass", jetmass);
    chain.SetBranchAddress("Jet_pt", jetpt);
    chain.SetBranchAddress("Jet_phi", jetphi);
    chain.SetBranchAddress("Jet_eta", jeteta);
    chain.SetBranchAddress("Jet_btagCSVV2", jetcsv);
    // chain.SetBranchAddress("Jet_IsCSVM", jetiscsvm_pre);
    chain.SetBranchAddress("nJet", &jetSize);
    chain.SetBranchAddress("Jet_jetId", jetId);
    chain.SetBranchAddress("Jet_partonFlavour", jetpartonflavour);
    //chain.SetBranchAddress("Jet_IsCSVL", jetisCSVL);
    chain.SetBranchAddress("Jet_btagDeepB", jetDeepCSV);
    // chain.SetBranchAddress("Jet_IsDeepCSVL", jetisDeepCSVL);
    // chain.SetBranchAddress("Jet_IsDeepCSVM", jetisDeepCSVM);
    // chain.SetBranchAddress("Jet_IsDeepCSVT", jetisDeepCSVT);
    chain.SetBranchAddress("Jet_nConstituents", Jet_nConstituents);
    chain.SetBranchAddress("JetPFCandsAK4_candIdx", JetPFCands_candIdx);
    chain.SetBranchAddress("JetPFCandsAK4_jetIdx", JetPFCands_jetIdx);
    chain.SetBranchAddress("nJetPFCandsAK4", &nJetPFCands);
    chain.SetBranchAddress("JetPFCands_mass", JetPFCands_mass);
    chain.SetBranchAddress("JetPFCands_pt", JetPFCands_pt);
    chain.SetBranchAddress("JetPFCands_eta", JetPFCands_eta);
    chain.SetBranchAddress("JetPFCands_phi", JetPFCands_phi);
    chain.SetBranchAddress("JetPFCands_d0", JetPFCands_d0);
    chain.SetBranchAddress("JetPFCands_dz", JetPFCands_dz);
    chain.SetBranchAddress("JetPFCands_pdgId", JetPFCands_pdgId);

    // //JET AK8
    // chain.SetBranchAddress("jetsAK8Puppi_E", fatjete);
    // chain.SetBranchAddress("jetsAK8Puppi_Pt", fatjetpt);
    // chain.SetBranchAddress("jetsAK8Puppi_Eta", fatjeteta);
    // chain.SetBranchAddress("jetsAK8Puppi_Phi", fatjetphi);
    // chain.SetBranchAddress("jetsAK8Puppi_TopMass", fatjettopmass);
    // chain.SetBranchAddress("jetsAK8Puppi_PartonFlavour", jethadronAK8flavour); //Lorenzo
    // chain.SetBranchAddress("jetsAK8Puppi_nDeepCSVsubj_l", fatjet_nDeepCSVsubj);
    // //chain.SetBranchAddress("jetsAK8Puppi_nCSVsubj_ l", fatjet_nCSVsubj);
    // //chain.SetBranchAddress("subjetsAK8Puppi_DeepCSV", subj_DeepCSV_fatjet);
    // chain.SetBranchAddress("jetsAK8Puppi_tau3OVERtau2", fatjettau3OVER2);
    // chain.SetBranchAddress("jetsAK8Puppi_tau2OVERtau1", fatjettau2OVER1);
    // chain.SetBranchAddress("jetsAK8Puppi_tau3", fatjettau3);
    // chain.SetBranchAddress("jetsAK8Puppi_tau2", fatjettau2);
    // chain.SetBranchAddress("jetsAK8Puppi_tau1", fatjettau1);
    // //chain.SetBranchAddress("jetsAK8Puppi_prunedMassPuppi", fatjetprunedmass);
    // chain.SetBranchAddress("jetsAK8Puppi_prunedMass", Prunedmass);
    // chain.SetBranchAddress("jetsAK8Puppi_softDropMass", SDmass);
    // //chain.SetBranchAddress("jetsAK8Puppi_prunedMass", SDmass);
    // chain.SetBranchAddress("jetsAK8Puppi_size", &fatjetSize);
    // chain.SetBranchAddress("jetsAK8Puppi_nCSVsubj", fatjetsubjets_Ncsv);
    // chain.SetBranchAddress("jetsAK8Puppi_nCSVsubj_tm", fatjetsubjets_Ncsv_tm);
    // chain.SetBranchAddress("jetsAK8Puppi_nJ", fatjetsubjetsN);

    //MUONS
    // chain.SetBranchAddress("Event_nLooseMuons", nLooseMuons);
    // chain.SetBranchAddress("Event_nLooseElectrons", nLooseElectrons);
    chain.SetBranchAddress("Muon_looseId", muon_IsLoosePtMuon);
    chain.SetBranchAddress("Muon_tightId", muon_IsTightPtMuon);
    chain.SetBranchAddress("Muon_highPtId", muon_IsHighPtMuon); //#1<<<<
    chain.SetBranchAddress("Muon_pfRelIso04_all", muons_Iso04);
    chain.SetBranchAddress("Muon_pt", muons_Pt);     //FIXME
    chain.SetBranchAddress("Muon_mass", muons_Mass); //FIXME
    chain.SetBranchAddress("nMuon", &muons_size);
    chain.SetBranchAddress("Muon_phi", muons_Phi); //Lorenzo
    chain.SetBranchAddress("Muon_eta", muons_Eta); // Lorenzo
    chain.SetBranchAddress("Muon_charge", muonsCharge);
    // chain.SetBranchAddress("Event_nMuonsSF", &nMuons);
    // chain.SetBranchAddress("muonsTight_Pt", muonsTight_Pt); //FIXME
    // chain.SetBranchAddress("muonsTight_size", &muonsTight_size);
    // chain.SetBranchAddress("muonsTight_Phi", muonsTight_Phi); //Lorenzo
    // chain.SetBranchAddress("muonsTight_Eta", muonsTight_Eta); // Lorenzo
    // chain.SetBranchAddress("Event_nTightMuons", nTightMuons);
    // chain.SetBranchAddress("muonsLoose_Pt", muonsLoose_Pt);
    // chain.SetBranchAddress("muonsLoose_size", &muonsLoose_size); //
    // chain.SetBranchAddress("muonsLoose_Phi", muonsLoose_Phi);
    // chain.SetBranchAddress("muons_DB", muons_DB);
    // chain.SetBranchAddress("muons_DBerr", muons_DBerr);
    chain.SetBranchAddress("Muon_dz", muons_Dz);
    chain.SetBranchAddress("Muon_dxy", muons_Dxy);
    chain.SetBranchAddress("Muon_dxyErr", muons_Dxyerr);
    chain.SetBranchAddress("Muon_isGlobal", muons_IsGlobalMuon);
    chain.SetBranchAddress("Muon_isTracker", muons_IsTrackerMuon);
    chain.SetBranchAddress("Muon_nStations", muons_NumberMatchedStations);
    chain.SetBranchAddress("Muon_miniPFRelIso_all", muons_MiniIso);
    chain.SetBranchAddress("Muon_nTrackerLayers", muons_NumberOfValidTrackerHits);

    //-------Qui consideriamo i trees che vanno nell'output file----------//

    TLorentzVector muon_vect, jet_vect, muon_vect_gen, b_vect, met_vect, el_vect_gen, muon_vect_boosted, jet_vect_boosted;
    float deltaRTemp;
    double costhetap;

    TopUtilities Top1, top2;

    int PFCandsSizeMax_output = 5;

    float mu_pt_merged = 0.;
    float mu_e_merged = 0.;
    float mu_phi_merged = 0.;
    float mu_eta_merged = 0.;
    float mu_ch_merged = 0.;
    // float mu_DB_merged = 0.;
    // float mu_DBerr_merged = 0.;
    float mu_Dz_merged = 0.;
    float mu_Dxy_merged = 0.;
    float mu_Dxyerr_merged = 0.;
    float mu_MiniIso_merged = 0.;
    float mu_Iso_merged = 0.;
    float mu_IsGlobal_merged = 0.;
    float mu_IsTracker_merged = 0.;
    float mu_NumberMatchedStations_merged = 0.;
    float mu_NumberOfValidTrackerHits_merged = 0.;
    float mu_Dxy_fract_merged = 0.;
    float pt_rel_merged = 0.;
    float mu_isHigh_merged = 0.;
    float mu_isTight_merged = 0.;
    float mu_isLoose_merged = 0.;
    int mu_high_truth_merged = 0;
    int tau_high_truth_merged = 0;
    int musize_merged = 0;
    float nHadZ_merged = 0;
    float nLepTop_merged = 0;

    float nHadZ_merged_reco = 0;
    float nLepTop_merged_reco = 0;

    float jet_isDeepCSVL_merged = 0;
    float jet_isDeepCSVM_merged = 0;
    float jet_isDeepCSVT_merged = 0;
    float jet_DeepCSV_merged = 0;
    float jet_partonFlavour_merged = 0;

    float jet_pt_merged = 0.;
    float jet_e_merged = 0.;
    float jet_phi_merged = 0.;
    float jet_eta_merged = 0.;
    int jet_high_truth_merged = 0;
    int jetsize_merged = 0;

    float JetPFCands_mass_merged[PFCandsSizeMax_output] = {0};
    float JetPFCands_pt_merged[PFCandsSizeMax_output] = {0};
    float JetPFCands_eta_merged[PFCandsSizeMax_output] = {0};
    float JetPFCands_phi_merged[PFCandsSizeMax_output] = {0};
    float JetPFCands_d0_merged[PFCandsSizeMax_output] = {0};
    float JetPFCands_dz_merged[PFCandsSizeMax_output] = {0};
    float JetPFCands_pdgId_merged[PFCandsSizeMax_output] = {0};

    int jet_hasPromptLep_merged = 0;

    float top_pt_merged = 0.;
    float top_e_merged = 0.;
    float top_phi_merged = 0.;
    float top_eta_merged = 0.;
    float top_M_merged = 0.;
    float top_Mt_merged = 0.;
    float rate_top_M_pT_merged = 0.;
    float top_M_corr_merged = 0.;

    int top_high_truth_merged = 0;
    int topsize_merged = 0;
    int top_category_merged = 0;

    float mu_boosted_pt_merged = 0.;
    float mu_boosted_e_merged = 0.;
    float mu_boosted_phi_merged = 0.;
    float mu_boosted_eta_merged = 0.;

    float jet_boosted_pt_merged = 0.;
    float jet_boosted_e_merged = 0.;
    float jet_boosted_phi_merged = 0.;
    float jet_boosted_eta_merged = 0.;

    float top_nu_pt_merged = 0.;
    float top_nu_e_merged = 0.;
    float top_nu_phi_merged = 0.;
    float top_nu_eta_merged = 0.;
    float top_nu_M_merged = 0.;

    float met_pt_merged = 0;
    float met_phi_merged = 0;
    int neutrino_number_merged = 0;

    double costheta_merged = 0.;

    float mu_pt_resolved = 0.;
    float mu_e_resolved = 0.;
    float mu_phi_resolved = 0.;
    float mu_eta_resolved = 0.;
    float mu_ch_resolved = 0.;
    // float mu_DB_resolved = 0.;
    // float mu_DBerr_resolved = 0.;
    float mu_Dz_resolved = 0.;
    float mu_Dxy_resolved = 0.;
    float mu_Dxyerr_resolved = 0.;
    float mu_MiniIso_resolved = 0.;
    float mu_Iso_resolved = 0.;
    float mu_IsGlobal_resolved = 0.;
    float mu_IsTracker_resolved = 0.;
    float mu_Dxy_fract_resolved = 0.;
    float pt_rel_resolved = 0.;
    float mu_isHigh_resolved = 0.;
    float mu_isTight_resolved = 0.;
    float mu_isLoose_resolved = 0.;
    float mu_NumberMatchedStations_resolved = 0.;
    float mu_NumberOfValidTrackerHits_resolved = 0.;

    int mu_high_truth_resolved = 0;
    int tau_high_truth_resolved = 0;
    int musize_resolved = 0;
    float nHadZ_resolved = 0;
    float nLepTop_resolved = 0;

    float nHadZ_resolved_reco = 0;
    float nLepTop_resolved_reco = 0;

    float jet_DeepCSV_resolved = 0;
    float jet_isDeepCSVL_resolved = 0;
    float jet_isDeepCSVM_resolved = 0;
    float jet_isDeepCSVT_resolved = 0;
    float jet_partonFlavour_resolved = 0;

    float jet_pt_resolved = 0.;
    float jet_e_resolved = 0.;
    float jet_phi_resolved = 0.;
    float jet_eta_resolved = 0.;
    int jet_high_truth_resolved = 0;
    int jetsize_resolved = 0;
    int jet_hasPromptLep_resolved = 0;

    float JetPFCands_mass_resolved[PFCandsSizeMax_output] = {0};
    float JetPFCands_pt_resolved[PFCandsSizeMax_output] = {0};
    float JetPFCands_eta_resolved[PFCandsSizeMax_output] = {0};
    float JetPFCands_phi_resolved[PFCandsSizeMax_output] = {0};
    float JetPFCands_d0_resolved[PFCandsSizeMax_output] = {0};
    float JetPFCands_dz_resolved[PFCandsSizeMax_output] = {0};
    float JetPFCands_pdgId_resolved[PFCandsSizeMax_output] = {0};

    float top_pt_resolved = 0.;
    float top_e_resolved = 0.;
    float top_phi_resolved = 0.;
    float top_eta_resolved = 0.;
    float top_M_resolved = 0.;
    float top_Mt_resolved = 0.;
    float rate_top_M_pT_resolved = 0.;
    float top_M_corr_resolved = 0.;

    int top_high_truth_resolved = 0;
    int topsize_resolved = 0;
    int top_category_resolved = 0;

    float mu_boosted_pt_resolved = 0.;
    float mu_boosted_e_resolved = 0.;
    float mu_boosted_phi_resolved = 0.;
    float mu_boosted_eta_resolved = 0.;

    float jet_boosted_pt_resolved = 0.;
    float jet_boosted_e_resolved = 0.;
    float jet_boosted_phi_resolved = 0.;
    float jet_boosted_eta_resolved = 0.;

    float top_nu_pt_resolved = 0.;
    float top_nu_e_resolved = 0.;
    float top_nu_phi_resolved = 0.;
    float top_nu_eta_resolved = 0.;
    float top_nu_M_resolved = 0.;

    float met_pt_resolved = 0;
    float met_phi_resolved = 0;
    int neutrino_number_resolved = 0;

    double costheta_resolved = 0.;

    TFile *f2 = new TFile((outFile).c_str(), "RECREATE");

    TTree *is_merged = new TTree("is_merged", "is_merged");
    //branch 1 top merged
    is_merged->Branch("Muon_Pt", &mu_pt_merged);
    is_merged->Branch("Muon_Phi", &mu_phi_merged);
    is_merged->Branch("Muon_Eta", &mu_eta_merged);
    is_merged->Branch("Muon_Mass", &mu_e_merged);
    is_merged->Branch("Muon_Size", &musize_merged);
    is_merged->Branch("Muon_Charge", &mu_ch_merged);
    // is_merged->Branch("Muon_DB", &mu_DB_merged);
    // is_merged->Branch("Muon_DBerr", &mu_DBerr_merged);
    is_merged->Branch("Muon_Dz", &mu_Dz_merged);
    is_merged->Branch("Muon_MiniIso", &mu_MiniIso_merged);
    is_merged->Branch("Muon_Iso", &mu_Iso_merged);
    is_merged->Branch("Muon_IsTrackerMuon", &mu_IsTracker_merged);
    is_merged->Branch("Muon_NumberMatchedStations", &mu_NumberMatchedStations_merged);
    is_merged->Branch("Muon_NumberOfValidTrackerHits", &mu_NumberOfValidTrackerHits_merged);
    is_merged->Branch("Muon_Dxy", &mu_Dxy_merged);
    is_merged->Branch("Muon_Dxyerr", &mu_Dxyerr_merged);
    is_merged->Branch("Muon_IsGlobalMuon", &mu_IsGlobal_merged);
    is_merged->Branch("Muon_High_Truth", &mu_high_truth_merged);
    is_merged->Branch("Muon_Dxy_fract", &mu_Dxy_fract_merged);
    is_merged->Branch("Muon_Pt_Rel", &pt_rel_merged);
    is_merged->Branch("Muon_isHigh", &mu_isHigh_merged);
    is_merged->Branch("Muon_isTight", &mu_isTight_merged);
    is_merged->Branch("Muon_isLoose", &mu_isLoose_merged);

    is_merged->Branch("Jet_High_Truth", &jet_high_truth_merged);
    is_merged->Branch("Top_High_Truth", &top_high_truth_merged);
    is_merged->Branch("Tau_High_Truth", &tau_high_truth_merged);
    is_merged->Branch("Top_Category", &top_category_merged);

    is_merged->Branch("Jet_Pt", &jet_pt_merged);
    is_merged->Branch("Jet_Phi", &jet_phi_merged);
    is_merged->Branch("Jet_Eta", &jet_eta_merged);
    is_merged->Branch("Jet_E", &jet_e_merged);
    is_merged->Branch("Jet_Size", &jetsize_merged);

    is_merged->Branch("JetPFCands_mass", JetPFCands_mass_merged);
    is_merged->Branch("JetPFCands_pt", JetPFCands_pt_merged);
    is_merged->Branch("JetPFCands_eta", JetPFCands_eta_merged);
    is_merged->Branch("JetPFCands_phi", JetPFCands_phi_merged);
    is_merged->Branch("JetPFCands_d0", JetPFCands_d0_merged);
    is_merged->Branch("JetPFCands_dz", JetPFCands_dz_merged);
    is_merged->Branch("JetPFCands_pdgId", JetPFCands_pdgId_merged);

    is_merged->Branch("top_Pt", &top_pt_merged);
    is_merged->Branch("top_Phi", &top_phi_merged);
    is_merged->Branch("top_Eta", &top_eta_merged);
    is_merged->Branch("top_E", &top_e_merged);
    is_merged->Branch("top_M", &top_M_merged);
    is_merged->Branch("top_Mt", &top_Mt_merged);
    is_merged->Branch("rate_top_M_pT", &rate_top_M_pT_merged);
    is_merged->Branch("top_M_corr", &top_M_corr_merged);
    is_merged->Branch("top_Size", &topsize_merged);
    is_merged->Branch("Event_nHadZ", &nHadZ_merged);
    is_merged->Branch("Event_nLepTop", &nLepTop_merged);

    is_merged->Branch("Event_nHadZ_reco", &nHadZ_merged_reco);
    is_merged->Branch("Event_nLepTop_reco", &nLepTop_merged_reco);

    is_merged->Branch("Jet_isDeepCSVL", &jet_isDeepCSVL_merged);
    is_merged->Branch("Jet_isDeepCSVM", &jet_isDeepCSVM_merged);
    is_merged->Branch("Jet_isDeepCSVT", &jet_isDeepCSVT_merged);
    is_merged->Branch("Jet_DCSV", &jet_DeepCSV_merged);
    is_merged->Branch("Jet_PartonFlavour", &jet_partonFlavour_merged);
    is_merged->Branch("Jet_hasPromptLepton", &jet_hasPromptLep_merged);

    is_merged->Branch("Muon_Boosted_Pt", &mu_boosted_pt_merged);
    is_merged->Branch("Muon_Boosted_Phi", &mu_boosted_phi_merged);
    is_merged->Branch("Muon_Boosted_Eta", &mu_boosted_eta_merged);
    is_merged->Branch("Muon_Boosted_E", &mu_boosted_e_merged);

    is_merged->Branch("Jet_Boosted_Pt", &jet_boosted_pt_merged);
    is_merged->Branch("Jet_Boosted_Phi", &jet_boosted_phi_merged);
    is_merged->Branch("Jet_Boosted_Eta", &jet_boosted_eta_merged);
    is_merged->Branch("Jet_Boosted_E", &jet_boosted_e_merged);

    is_merged->Branch("top_nu_Pt", &top_nu_pt_merged);
    is_merged->Branch("top_nu_Phi", &top_nu_phi_merged);
    is_merged->Branch("top_nu_Eta", &top_nu_eta_merged);
    is_merged->Branch("top_nu_E", &top_nu_e_merged);
    is_merged->Branch("top_nu_M", &top_nu_M_merged);

    is_merged->Branch("MET_Pt", &met_pt_merged);
    is_merged->Branch("MET_Phi", &met_phi_merged);
    is_merged->Branch("Neutrino_number", &neutrino_number_merged);

    is_merged->Branch("Costheta", &costheta_merged);

    TTree *is_resolved = new TTree("is_resolved", "is_resolved");
    //branch 1 top resolved
    is_resolved->Branch("Muon_Pt", &mu_pt_resolved);
    is_resolved->Branch("Muon_Phi", &mu_phi_resolved);
    is_resolved->Branch("Muon_Eta", &mu_eta_resolved);
    is_resolved->Branch("Muon_E", &mu_e_resolved);
    is_resolved->Branch("Muon_Size", &musize_resolved);
    is_resolved->Branch("Muon_Charge", &mu_ch_resolved);
    // is_resolved->Branch("Muon_DB", &mu_DB_resolved);
    // is_resolved->Branch("Muon_DBerr", &mu_DBerr_resolved);
    is_resolved->Branch("Muon_Dz", &mu_Dz_resolved);
    is_resolved->Branch("Muon_MiniIso", &mu_MiniIso_resolved);
    is_resolved->Branch("Muon_Iso", &mu_Iso_resolved);
    is_resolved->Branch("Muon_IsTrackerMuon", &mu_IsTracker_resolved);
    is_resolved->Branch("Muon_NumberMatchedStations", &mu_NumberMatchedStations_resolved);
    is_resolved->Branch("Muon_NumberOfValidTrackerHits", &mu_NumberOfValidTrackerHits_resolved);
    is_resolved->Branch("Muon_Dxy", &mu_Dxy_resolved);
    is_resolved->Branch("Muon_Dxyerr", &mu_Dxyerr_resolved);
    is_resolved->Branch("Muon_IsGlobalMuon", &mu_IsGlobal_resolved);
    is_resolved->Branch("Muon_Dxy_fract", &mu_Dxy_fract_resolved);
    is_resolved->Branch("Muon_Pt_Rel", &pt_rel_resolved);
    is_resolved->Branch("Muon_isHigh", &mu_isHigh_resolved);
    is_resolved->Branch("Muon_isTight", &mu_isTight_resolved);
    is_resolved->Branch("Muon_isLoose", &mu_isLoose_resolved);
    is_resolved->Branch("Muon_High_Truth", &mu_high_truth_resolved);
    is_resolved->Branch("Jet_High_Truth", &jet_high_truth_resolved);
    is_resolved->Branch("Top_High_Truth", &top_high_truth_resolved);
    is_resolved->Branch("Tau_High_Truth", &tau_high_truth_resolved);
    is_resolved->Branch("Top_Category", &top_category_resolved);

    is_resolved->Branch("Jet_Pt", &jet_pt_resolved);
    is_resolved->Branch("Jet_Phi", &jet_phi_resolved);
    is_resolved->Branch("Jet_Eta", &jet_eta_resolved);
    is_resolved->Branch("Jet_E", &jet_e_resolved);
    is_resolved->Branch("Jet_Size", &jetsize_resolved);

    is_resolved->Branch("JetPFCands_mass", JetPFCands_mass_resolved);
    is_resolved->Branch("JetPFCands_pt", JetPFCands_pt_resolved);
    is_resolved->Branch("JetPFCands_eta", JetPFCands_eta_resolved);
    is_resolved->Branch("JetPFCands_phi", JetPFCands_phi_resolved);
    is_resolved->Branch("JetPFCands_d0", JetPFCands_d0_resolved);
    is_resolved->Branch("JetPFCands_dz", JetPFCands_dz_resolved);
    is_resolved->Branch("JetPFCands_pdgId", JetPFCands_pdgId_resolved);

    is_resolved->Branch("top_Pt", &top_pt_resolved);
    is_resolved->Branch("top_Phi", &top_phi_resolved);
    is_resolved->Branch("top_Eta", &top_eta_resolved);
    is_resolved->Branch("top_E", &top_e_resolved);
    is_resolved->Branch("top_M", &top_M_resolved);
    is_resolved->Branch("top_Mt", &top_Mt_resolved);
    is_resolved->Branch("rate_top_M_pT", &rate_top_M_pT_resolved);
    is_resolved->Branch("top_M_corr", &top_M_corr_resolved);
    is_resolved->Branch("top_Size", &topsize_resolved);
    is_resolved->Branch("Event_nHadZ", &nHadZ_resolved);
    is_resolved->Branch("Event_nLepTop", &nLepTop_resolved);
    is_resolved->Branch("Event_nHadZ_reco", &nHadZ_resolved_reco);
    is_resolved->Branch("Event_nLepTop_reco", &nLepTop_resolved_reco);

    is_resolved->Branch("Jet_DCSV", &jet_DeepCSV_resolved);
    is_resolved->Branch("Jet_isDeepCSVL", &jet_isDeepCSVL_resolved);
    is_resolved->Branch("Jet_isDeepCSVM", &jet_isDeepCSVM_resolved);
    is_resolved->Branch("Jet_isDeepCSVT", &jet_isDeepCSVT_resolved);
    is_resolved->Branch("Jet_PartonFlavour", &jet_partonFlavour_resolved);
    is_resolved->Branch("Jet_hasPromptLepton", &jet_hasPromptLep_resolved);

    is_resolved->Branch("Muon_Boosted_Pt", &mu_boosted_pt_resolved);
    is_resolved->Branch("Muon_Boosted_Phi", &mu_boosted_phi_resolved);
    is_resolved->Branch("Muon_Boosted_Eta", &mu_boosted_eta_resolved);
    is_resolved->Branch("Muon_Boosted_E", &mu_boosted_e_resolved);

    is_resolved->Branch("Jet_Boosted_Pt", &jet_boosted_pt_resolved);
    is_resolved->Branch("Jet_Boosted_Phi", &jet_boosted_phi_resolved);
    is_resolved->Branch("Jet_Boosted_Eta", &jet_boosted_eta_resolved);
    is_resolved->Branch("Jet_Boosted_E", &jet_boosted_e_resolved);

    is_resolved->Branch("top_nu_Pt", &top_nu_pt_resolved);
    is_resolved->Branch("top_nu_Phi", &top_nu_phi_resolved);
    is_resolved->Branch("top_nu_Eta", &top_nu_eta_resolved);
    is_resolved->Branch("top_nu_E", &top_nu_e_resolved);
    is_resolved->Branch("top_nu_M", &top_nu_M_resolved);

    is_resolved->Branch("MET_Pt", &met_pt_resolved);
    is_resolved->Branch("MET_Phi", &met_phi_resolved);
    is_resolved->Branch("Neutrino_number", &neutrino_number_resolved);

    is_resolved->Branch("Costheta", &costheta_resolved);
    is_resolved->Branch("Muon_Pt_Rel", &pt_rel_merged);
    float Event_Number = 0.;
    float Event_run = 0.;

    vector<int> *Z_index = new vector<int>;
    vector<int> *Z_b_index = new vector<int>;

    int counter_merged = 0, counter_resolved = 0;

    TH2F *histo_DP_DR = new TH2F("histo_DR_DP", "histo_DR_DP", 100, 0, 0.8, 100, -0.4, 0.4);

    for (int i = 0; i < chain.GetEntries(); i++)
    {
        chain.GetEntry(i);

        cout << "Event: " << i << endl;

        mu_pt_resolved = 0;
        mu_e_resolved = 0;
        mu_phi_resolved = 0;
        mu_eta_resolved = 0;
        mu_ch_resolved = 0;
        // mu_DB_resolved = 0;
        // mu_DBerr_resolved = 0;
        mu_Dz_resolved = 0;
        mu_Dxy_resolved = 0;
        mu_Dxyerr_resolved = 0;
        mu_MiniIso_resolved = 0;
        mu_Iso_resolved = 0;
        mu_IsGlobal_resolved = 0;
        mu_IsTracker_resolved = 0;
        mu_NumberMatchedStations_resolved = 0;
        mu_NumberOfValidTrackerHits_resolved = 0;
        mu_Dxy_fract_resolved = 0.;
        pt_rel_resolved = 0.;
        mu_isHigh_resolved = 0.;
        mu_isTight_resolved = 0.;
        mu_isLoose_resolved = 0.;
        mu_high_truth_resolved = 0;
        jet_high_truth_resolved = 0;
        top_category_resolved = 0;

        top_high_truth_resolved = 0;
        nHadZ_resolved = 0;
        nLepTop_resolved = 0;
        nHadZ_resolved_reco = 0;
        nLepTop_resolved_reco = 0;

        jet_isDeepCSVL_resolved = 0;
        jet_isDeepCSVM_resolved = 0;
        jet_isDeepCSVT_resolved = 0;
        jet_DeepCSV_resolved = 0;
        jet_partonFlavour_resolved = 0;
        jet_hasPromptLep_resolved = 0;

        jet_pt_resolved = 0;
        jet_e_resolved = 0;
        jet_phi_resolved = 0;
        jet_eta_resolved = 0;

        for (int y = 0; y < PFCandsSizeMax_output; y++)
        {
            JetPFCands_mass_resolved[y] = 0;
            JetPFCands_pt_resolved[y] = 0;
            JetPFCands_eta_resolved[y] = 0;
            JetPFCands_phi_resolved[y] = 0;
            JetPFCands_d0_resolved[y] = 0;
            JetPFCands_dz_resolved[y] = 0;
            JetPFCands_pdgId_resolved[y] = 0;
        }
        top_pt_resolved = 0;
        top_e_resolved = 0;
        top_phi_resolved = 0;
        top_eta_resolved = 0;
        top_M_resolved = 0;
        top_Mt_resolved = 0;
        top_M_corr_resolved = 0;
        rate_top_M_pT_resolved = 0;

        mu_boosted_pt_resolved = 0.;
        mu_boosted_e_resolved = 0.;
        mu_boosted_phi_resolved = 0.;
        mu_boosted_eta_resolved = 0.;

        tau_high_truth_resolved = 0.;

        jet_boosted_pt_resolved = 0.;
        jet_boosted_e_resolved = 0.;
        jet_boosted_phi_resolved = 0.;
        jet_boosted_eta_resolved = 0.;

        top_nu_pt_resolved = 0.;
        top_nu_e_resolved = 0.;
        top_nu_phi_resolved = 0.;
        top_nu_eta_resolved = 0.;
        top_nu_M_resolved = 0.;

        met_pt_resolved = 0;
        met_phi_resolved = 0;

        costheta_resolved = 0.;

        mu_pt_merged = 0;
        mu_e_merged = 0;
        mu_phi_merged = 0;
        mu_eta_merged = 0;
        mu_ch_merged = 0;
        // mu_DB_merged = 0;
        mu_Dxy_fract_merged = 0;
        mu_isHigh_merged = 0.;
        mu_isTight_merged = 0.;
        mu_isLoose_merged = 0.;
        pt_rel_merged = 0.;
        mu_Dz_merged = 0;
        mu_Dxy_merged = 0;
        mu_Dxyerr_merged = 0;
        mu_MiniIso_merged = 0;
        mu_Iso_merged = 0;
        mu_IsGlobal_merged = 0;
        mu_IsTracker_merged = 0;
        mu_NumberMatchedStations_merged = 0;
        mu_NumberOfValidTrackerHits_merged = 0;
        mu_high_truth_merged = 0;
        jet_high_truth_merged = 0;
        top_high_truth_merged = 0;
        tau_high_truth_merged = 0.;
        top_category_merged = 0;

        nHadZ_merged = 0;
        nLepTop_merged = 0;
        nHadZ_merged_reco = 0;
        nLepTop_merged_reco = 0;

        jet_isDeepCSVL_merged = 0;
        jet_isDeepCSVM_merged = 0;
        jet_isDeepCSVT_merged = 0;
        jet_DeepCSV_merged = 0;
        jet_partonFlavour_merged = 0;
        jet_hasPromptLep_merged = 0;

        jet_pt_merged = 0;
        jet_e_merged = 0;
        jet_phi_merged = 0;
        jet_eta_merged = 0;

        for (int y = 0; y < PFCandsSizeMax_output; y++)
        {
            JetPFCands_mass_merged[y] = 0;
            JetPFCands_pt_merged[y] = 0;
            JetPFCands_eta_merged[y] = 0;
            JetPFCands_phi_merged[y] = 0;
            JetPFCands_d0_merged[y] = 0;
            JetPFCands_dz_merged[y] = 0;
            JetPFCands_pdgId_merged[y] = 0;
        }

        top_pt_merged = 0;
        top_e_merged = 0;
        top_phi_merged = 0;
        top_eta_merged = 0;
        top_M_merged = 0;
        top_Mt_merged = 0;
        top_M_corr_merged = 0;
        rate_top_M_pT_merged = 0;

        mu_boosted_pt_merged = 0.;
        mu_boosted_e_merged = 0.;
        mu_boosted_phi_merged = 0.;
        mu_boosted_eta_merged = 0.;

        jet_boosted_pt_merged = 0.;
        jet_boosted_e_merged = 0.;
        jet_boosted_phi_merged = 0.;
        jet_boosted_eta_merged = 0.;

        top_nu_pt_merged = 0.;
        top_nu_e_merged = 0.;
        top_nu_phi_merged = 0.;
        top_nu_eta_merged = 0.;
        top_nu_M_merged = 0.;

        met_pt_merged = 0;
        met_phi_merged = 0;

        costheta_merged = 0.;
        Z_index->clear();
        Z_b_index->clear();

        float top_rec[100];
        int muon_high_truth[100];
        int top_high_truth[100];
        int jet_high_truth[100];
        int tau_high_truth[100];

        int jet_hasPromptLep[100];

        int Z_tag = 0;
        int nZ_had = 0;

        met_vect.SetPtEtaPhiE(metptCorr[0], 0, metphiCorr[0], metptCorr[0]);

        // for (int y = 0; y < fatjetSize; y++)
        // {
        // TLorentzVector fatjet_vect;
        // fatjet_vect.SetPtEtaPhiE(fatjetpt[y], fatjeteta[y], fatjetphi[y], fatjete[y]);

        //     if (SDmass[y] < 105 && SDmass[y] > 60 && fatjettau2OVER1[y] < 0.45)
        //     {
        //         if (fatjetpt[y] > 200)
        //         {

        //             for (int k = 0; k < jetSize; k++)
        //             {
        //                 if (jetisDeepCSVL[k] == 1)
        //                 {
        //                     jet_vect.SetPtEtaPhiE(jetpt[k], jeteta[k], jetphi[k], jete[k]);
        //                     float delta_R_true_corrected = sqrt(pow(fatjeteta[y] - jeteta[k], 2) + pow(fatjet_vect.DeltaPhi(jet_vect), 2));
        //                     if (delta_R_true_corrected > 1.2)
        //                     {
        //                         Z_tag++;
        //                         Z_index->push_back(y);
        //                         Z_b_index->push_back(k);
        //                     }
        //                 }
        //             }
        //         }
        //     }
        //     if (Z_tag != 0)
        //     {
        //         nZ_had++;
        //     }
        // }

        // RICHIEDERE UN MUONE DI ALTO PT

        /// RICOSTRUZIONE MUONI

        // SELEZIONE MUONI TRUE E FALSE HIGH PT

        topsize_merged = 0;
        topsize_resolved = 0;

        neutrino_number_merged = 0;
        neutrino_number_resolved = 0;

        for (int u = 0; u < genP_size; u++)
        {
            int index1 = 0;
            index1 = (int)(genP_mom1I[u] - 1);

            if ((index1 > -1) && ((abs(genP_flavor[u]) == 12) || (abs(genP_flavor[u]) == 14) || (abs(genP_flavor[u]) == 16)) && ((abs(genP_flavor[index1]) == 24) || (abs(genP_flavor[index1]) == 23)))
            {
                neutrino_number_merged++;
                neutrino_number_resolved++;
            }
        }

        for (int k = 0; k < jetSize; k++) // Unpackiamo la size di muoni e jet
        {
            if (jetDeepCSV[k] != -2 && metptCorr[0] > 70)
            {
                for (int j = 0; j < muons_size; j++)
                {
                    jet_vect.SetPtEtaPhiM(jetpt[k], jeteta[k], jetphi[k], jetmass[k]);
                    if (muons_Pt[j] > 10)
                    {
                        muon_vect.SetPtEtaPhiM(muons_Pt[j], muons_Eta[j], muons_Phi[j], muons_Mass[j]);
                        float deltaR = sqrt(pow(jeteta[k] - muons_Eta[j], 2) + pow(muon_vect.DeltaPhi(jet_vect), 2));
                        if (deltaR < 2 && deltaR > .4)
                        {
                            topsize_resolved++;
                        }
                        if (deltaR < .4)
                        {
                            topsize_merged++;
                        }
                    }
                }
            }
        }

        for (int k = 0; k < jetSize; k++)
        {
            jet_high_truth[k] = 0;
            jet_hasPromptLep[k] = 0;
            int Z_b_match = 0;
            for (unsigned int f = 0; f < (Z_index->size()); f++)
            {
                if (k == Z_b_index->at(f))
                {
                    Z_b_match = 1;
                }
            }
            Z_b_match = 1;

            jet_vect.SetPtEtaPhiM(jetpt[k], jeteta[k], jetphi[k], jetmass[k]);

            if (jetDeepCSV[k] != -2 && Z_b_match != 0 && metptCorr[0] > 70)
            { //&& nTop_Lep>0

                for (int u = 0; u < genP_size; u++)
                {
                    int index1 = 0;
                    index1 = (int)(genP_mom1I[u] - 1);

                    if (index1 > -1)
                    {
                        if (abs(genP_flavor[u]) == 5 && (genP_flavor[u] * jetpartonflavour[k] > 0.) && abs(genP_flavor[index1]) == 6)
                        { //
                            TLorentzVector jet_vect_gen;
                            jet_vect_gen.SetPtEtaPhiM(genP_Pt[u], genP_Eta[u], genP_Phi[u], genP_Mass[u]);
                            float delta_R_true_corrected = sqrt(pow(genP_Eta[u] - jeteta[k], 2) + pow(jet_vect_gen.DeltaPhi(jet_vect), 2));
                            if (delta_R_true_corrected < 0.4)
                            {
                                jet_high_truth[k] = 1;
                            } // MC match
                        }     // part giusta
                    }
                } //loop generatore

                //loop per selezionare i PFCands a pt più alto
                vector<pair<int, float>> sorted_PFCands_indexes_pt;

                for (int j = 0; j < nJetPFCands; j++)
                {
                    if (JetPFCands_jetIdx[j] == k)
                    {
                        sorted_PFCands_indexes_pt.push_back(make_pair(JetPFCands_candIdx[j], JetPFCands_pt[JetPFCands_candIdx[j]]));
                    }
                }

                sort(sorted_PFCands_indexes_pt.begin(), sorted_PFCands_indexes_pt.end(), sortByVal);

                //loop generatore per controllare presenza prompt lepton in jet
                for (int u = 0; u < genP_size; u++)
                {
                    int indexLepton = 0;
                    indexLepton = (int)(genP_mom1I[u] - 1);

                    if (indexLepton > -1)
                    {
                        if (((abs(genP_flavor[u]) == 11) || (abs(genP_flavor[u]) == 13) || (abs(genP_flavor[u]) == 15)) && (abs(genP_flavor[indexLepton]) == 24))
                        {
                            TLorentzVector lepton_vect_gen;
                            lepton_vect_gen.SetPtEtaPhiM(genP_Pt[u], genP_Eta[u], genP_Phi[u], genP_Mass[u]);
                            float delta_R_lep_jet = sqrt(pow(genP_Eta[u] - jeteta[k], 2) + pow(lepton_vect_gen.DeltaPhi(jet_vect), 2));
                            if (delta_R_lep_jet < .4)
                                jet_hasPromptLep[k] = 1;
                        }
                    }
                }

                for (int j = 0; j < muons_size; j++)
                {
                    muon_high_truth[j] = 0;
                    tau_high_truth[j] = 0;
                    top_rec[j] = 0;
                    if (muons_Pt[j] > 10)
                    {
                        muon_vect.SetPtEtaPhiM(muons_Pt[j], muons_Eta[j], muons_Phi[j], muons_Mass[j]);

                        deltaRTemp = sqrt(pow(jeteta[k] - muons_Eta[j], 2) + pow(muon_vect.DeltaPhi(jet_vect), 2));

                        if (deltaRTemp < 2. && deltaRTemp > 0.4)
                        { //&& mu_isHigh->at(j)==1  && mu_Iso->at(j)<0.2

                            top_high_truth[j] = 0;
                            top_rec[j] = 1;
                            TLorentzVector tot_vect_3 = jet_vect + muon_vect;
                            math::PtEtaPhiELorentzVector tot_vect = Top1.top4Momentum(muon_vect, jet_vect, met_vect.Px(), met_vect.Py());
                            TLorentzVector tot_vect_1, tot_vect_2;
                            tot_vect_1.SetPxPyPzE(-tot_vect.Px(), -tot_vect.Py(), -tot_vect.Pz(), tot_vect.E());

                            top_nu_e_resolved = (tot_vect.E());
                            top_nu_phi_resolved = (tot_vect.Phi());
                            top_nu_pt_resolved = (tot_vect.Pt());
                            top_nu_eta_resolved = (tot_vect.Eta());
                            top_nu_M_resolved = (sqrt(tot_vect.M2()));

                            top_e_resolved = (tot_vect_3.E());
                            top_phi_resolved = (tot_vect_3.Phi());
                            top_pt_resolved = (tot_vect_3.Pt());
                            top_eta_resolved = (tot_vect_3.Eta());

                            top_M_resolved = (sqrt(tot_vect_3.M2()));
                            top_Mt_resolved = (tot_vect_3 + met_vect).Mt();
                            rate_top_M_pT_resolved = top_M_resolved / top_pt_resolved;
                            top_M_corr_resolved = sqrt(abs(top_M_resolved * top_M_resolved + 2 * muon_vect.P() * tot_vect_3.E() - 2 * top_pt_resolved * metptCorr[0] * cos(tot_vect_3.DeltaPhi(met_vect)) - 2 * tot_vect_3.Pz() * muon_vect.Pz()));
                            jet_high_truth_resolved = (jet_high_truth[k]);
                            tot_vect_2.SetPxPyPzE(tot_vect.Px(), tot_vect.Py(), tot_vect.Pz(), tot_vect.E());
                            costhetap = top2.costhetapol(muon_vect, jet_vect, tot_vect_2);

                            mu_pt_resolved = muons_Pt[j];
                            mu_e_resolved = muon_vect.E();
                            mu_phi_resolved = muons_Phi[j];
                            mu_eta_resolved = muons_Eta[j];
                            mu_ch_resolved = muonsCharge[j];

                            jet_pt_resolved = jetpt[k];
                            jet_e_resolved = jet_vect.E();
                            jet_phi_resolved = jetphi[k];
                            jet_eta_resolved = jeteta[k];

                            if (jetDeepCSV[k] > 0.1522)
                            {
                                jet_isDeepCSVL_resolved = 1;
                            }
                            if (jetDeepCSV[k] > 0.4941)
                            {
                                jet_isDeepCSVM_resolved = 1;
                            }
                            if (jetDeepCSV[k] > 0.8001)
                            {
                                jet_isDeepCSVT_resolved = 1;
                            }
                            jet_DeepCSV_resolved = jetDeepCSV[k];
                            jet_partonFlavour_resolved = jetpartonflavour[k];
                            jet_hasPromptLep_resolved = jet_hasPromptLep[k];

                            cout << "Event: " << i << "   Jet: " << k << "   Components: " << Jet_nConstituents[k] << endl; 

                            for (int y = 0; y < std::min(PFCandsSizeMax_output, Jet_nConstituents[k]); y++)
                            {
                                cout << "Index PFCand: " << sorted_PFCands_indexes_pt[y].first << "    PT PFCand: " << sorted_PFCands_indexes_pt[y].second << endl;
                                JetPFCands_mass_resolved[y] = JetPFCands_mass[sorted_PFCands_indexes_pt[y].first];
                                JetPFCands_pt_resolved[y] = JetPFCands_pt[sorted_PFCands_indexes_pt[y].first];
                                JetPFCands_eta_resolved[y] = JetPFCands_eta[sorted_PFCands_indexes_pt[y].first];
                                JetPFCands_phi_resolved[y] = JetPFCands_phi[sorted_PFCands_indexes_pt[y].first];
                                JetPFCands_d0_resolved[y] = JetPFCands_d0[sorted_PFCands_indexes_pt[y].first];
                                JetPFCands_dz_resolved[y] = JetPFCands_dz[sorted_PFCands_indexes_pt[y].first];
                                JetPFCands_pdgId_resolved[y] = JetPFCands_pdgId[sorted_PFCands_indexes_pt[y].first];
                            }

                            // mu_DB_resolved = muons_DB[j];
                            // mu_DBerr_resolved = muons_DBerr[j];
                            mu_Dz_resolved = muons_Dz[j];
                            mu_Dxy_resolved = muons_Dxy[j];
                            mu_Dxyerr_resolved = muons_Dxyerr[j];
                            mu_MiniIso_resolved = muons_MiniIso[j];
                            mu_Iso_resolved = muons_Iso04[j];
                            mu_IsGlobal_resolved = muons_IsGlobalMuon[j];
                            mu_IsTracker_resolved = muons_IsTrackerMuon[j];
                            mu_NumberMatchedStations_resolved = muons_NumberMatchedStations[j];
                            mu_NumberOfValidTrackerHits_resolved = muons_NumberOfValidTrackerHits[j];
                            mu_isHigh_resolved = muon_IsHighPtMuon[j];
                            mu_isTight_resolved = muon_IsTightPtMuon[j];
                            mu_isLoose_resolved = muon_IsLoosePtMuon[j];
                            mu_Dxy_fract_resolved = muons_Dxy[j] / (muons_Dxyerr[j]);
                            pt_rel_resolved = ((muon_vect.Vect()).Cross(jet_vect.Vect())).Mag() / ((jet_vect.Vect()).Mag());

                            met_pt_resolved = metptCorr[0];
                            met_phi_resolved = metptCorr[0];

                            for (int u = 0; u < genP_size; u++)
                            {
                                int index1 = genP_mom1I[u] - 1;
                                if (abs(genP_flavor[u]) == 13 && (genP_flavor[u] * muonsCharge[j] < 0.) && abs(genP_flavor[index1]) == 24)
                                {
                                    //                 if(abs(gP_flavor->at(u))==13 && (gP_flavor->at(u)*mu_ch->at(j)<0.) && abs(gP_flavor->at(index1))==24){
                                    muon_vect_gen.SetPtEtaPhiM(genP_Pt[u], genP_Eta[u], genP_Phi[u], genP_Mass[u]);
                                    float delta_R_true_corrected = sqrt(pow(genP_Eta[u] - muons_Eta[j], 2) + pow(muon_vect_gen.DeltaPhi(muon_vect), 2));
                                    if (delta_R_true_corrected < 0.1)
                                    {
                                        muon_high_truth[j] = 1;
                                    } // MC match
                                }     // part giusta
                            }         //loop generatore

                            for (int u = 0; u < genP_size; u++)
                            {
                                int index1 = genP_mom1I[u] - 1;
                                if (abs(genP_flavor[u]) == 15 && (genP_flavor[u] * muonsCharge[j] < 0.) && abs(genP_flavor[index1]) == 24)
                                {
                                    //                 if(abs(gP_flavor->at(u))==13 && (gP_flavor->at(u)*mu_ch->at(j)<0.) && abs(gP_flavor->at(index1))==24){
                                    muon_vect_gen.SetPtEtaPhiM(genP_Pt[u], genP_Eta[u], genP_Phi[u], genP_Mass[u]);
                                    float delta_R_true_corrected = sqrt(pow(genP_Eta[u] - muons_Eta[j], 2) + pow(muon_vect_gen.DeltaPhi(muon_vect), 2));
                                    if (delta_R_true_corrected < 0.05)
                                    {
                                        tau_high_truth[j] = 1;
                                    } // MC match
                                }     // part giusta
                            }         //loop generatore
                            muon_vect.Boost(tot_vect_1.BoostVector());
                            jet_vect.Boost(tot_vect_1.BoostVector());

                            jet_boosted_pt_resolved = (jet_vect.Pt());
                            jet_boosted_e_resolved = (jet_vect.E());
                            jet_boosted_phi_resolved = (jet_vect.Phi());
                            jet_boosted_eta_resolved = (jet_vect.Eta());

                            mu_boosted_pt_resolved = (muon_vect.Pt());
                            mu_boosted_e_resolved = (muon_vect.E());
                            mu_boosted_phi_resolved = (muon_vect.Phi());
                            mu_boosted_eta_resolved = (muon_vect.Eta());
                            costheta_resolved = (costhetap);

                            jet_vect.SetPtEtaPhiM(jetpt[k], jeteta[k], jetphi[k], jetmass[k]);
                            muon_vect.SetPtEtaPhiM(muons_Pt[j], muons_Eta[j], muons_Phi[j], muons_Mass[j]);

                            if ((muon_high_truth[j] == 1 || tau_high_truth[j] == 1) && jet_high_truth[k] == 1 && muonsCharge[j] * jetpartonflavour[k] > 0)
                            {
                                top_high_truth[j] = 1;
                                top_category_resolved = 3;
                            }

                            if ((muon_high_truth[j] == 0 && tau_high_truth[j] == 0) && jet_high_truth[k] == 0)
                            {
                                top_category_resolved = 0;
                            }

                            if ((muon_high_truth[j] == 0 && tau_high_truth[j] == 0) && jet_high_truth[k] == 1)
                            {
                                top_category_resolved = 1;
                            }

                            if ((muon_high_truth[j] == 1 || tau_high_truth[j] == 1) && jet_high_truth[k] == 0)
                            {
                                top_category_resolved = 2;
                            }

                            if ((muon_high_truth[j] == 1 || tau_high_truth[j] == 1) && jet_high_truth[k] == 1 && muonsCharge[j] * jetpartonflavour[k] < 0)
                            {
                                top_category_resolved = 4;
                            }

                            mu_high_truth_resolved = (muon_high_truth[j]);
                            tau_high_truth_resolved = tau_high_truth[j];
                            top_high_truth_resolved = (top_high_truth[j]);
                            counter_resolved = counter_resolved + 1;
                            nLepTop_resolved = nLepTop;
                            is_resolved->Fill();
                        } // selezione delta resolved

                        else if (deltaRTemp < 0.4)
                        { // && mu_isHigh->at(j)==1 && abs( mu_DB->at(j)/(mu_DBerr->at(j)))<2.5 && (mu_pt->at(j)/(jet_pt->at(k)))>0.1
                            top_high_truth[j] = 0;
                            top_rec[j] = -1;
                            if (false)
                                cout << " top_rec " << top_rec[j] << endl;
                            TLorentzVector tot_vect_3, jet_vect_tmp = jet_vect;
                            tot_vect_3 = jet_vect;
                            jet_vect = jet_vect - muon_vect;

                            math::PtEtaPhiELorentzVector tot_vect = Top1.top4Momentum(muon_vect, jet_vect, met_vect.Px(), met_vect.Py());

                            TLorentzVector tot_vect_1, tot_vect_2;
                            tot_vect_1.SetPxPyPzE(-tot_vect.Px(), -tot_vect.Py(), -tot_vect.Pz(), tot_vect.E());
                            top_nu_e_merged = (tot_vect.E());
                            top_nu_phi_merged = (tot_vect.Phi());
                            top_nu_pt_merged = (tot_vect.Pt());
                            top_nu_eta_merged = (tot_vect.Eta());

                            //top_reco_merged=(top_rec[j]);
                            top_nu_M_merged = sqrt(tot_vect.M2());
                            //muon_index_merged=(j);
                            //jet_index_merged=(k);
                            jet_high_truth_merged = (jet_high_truth[k]);
                            top_e_merged = (tot_vect_3.E());
                            top_phi_merged = (tot_vect_3.Phi());
                            top_pt_merged = (tot_vect_3.Pt());
                            top_eta_merged = (tot_vect_3.Eta());

                            top_M_merged = sqrt(tot_vect_3.M2());
                            top_Mt_merged = (tot_vect_3 + met_vect).Mt();
                            rate_top_M_pT_merged = top_M_merged / top_pt_merged;
                            top_M_corr_merged = sqrt(abs(top_M_merged * top_M_merged + 2 * muon_vect.P() * tot_vect_3.E() - 2 * top_pt_merged * metptCorr[0] * cos(tot_vect_3.DeltaPhi(met_vect)) - 2 * tot_vect_3.Pz() * muon_vect.Pz()));
                            tot_vect_2.SetPxPyPzE(tot_vect.Px(), tot_vect.Py(), tot_vect.Pz(), tot_vect.E());
                            costhetap = top2.costhetapol(muon_vect, jet_vect, tot_vect_2);

                            if (jetDeepCSV[k] > 0.1522)
                            {
                                jet_isDeepCSVL_merged = 1;
                            }
                            if (jetDeepCSV[k] > 0.4941)
                            {
                                jet_isDeepCSVM_merged = 1;
                            }
                            if (jetDeepCSV[k] > 0.8001)
                            {
                                jet_isDeepCSVT_merged = 1;
                            }
                            jet_DeepCSV_merged = jetDeepCSV[k];
                            jet_partonFlavour_merged = jetpartonflavour[k];
                            jet_hasPromptLep_merged = jet_hasPromptLep[k];

                            for (int y = 0; y < std::min(PFCandsSizeMax_output, Jet_nConstituents[k]); y++)
                            {
                                JetPFCands_mass_merged[y] = JetPFCands_mass[sorted_PFCands_indexes_pt[y].first];
                                JetPFCands_pt_merged[y] = JetPFCands_pt[sorted_PFCands_indexes_pt[y].first];
                                JetPFCands_eta_merged[y] = JetPFCands_eta[sorted_PFCands_indexes_pt[y].first];
                                JetPFCands_phi_merged[y] = JetPFCands_phi[sorted_PFCands_indexes_pt[y].first];
                                JetPFCands_d0_merged[y] = JetPFCands_d0[sorted_PFCands_indexes_pt[y].first];
                                JetPFCands_dz_merged[y] = JetPFCands_dz[sorted_PFCands_indexes_pt[y].first];
                                JetPFCands_pdgId_merged[y] = JetPFCands_pdgId[sorted_PFCands_indexes_pt[y].first];
                            }

                            mu_pt_merged = muons_Pt[j];
                            mu_e_merged = muon_vect.E();
                            mu_phi_merged = muons_Phi[j];
                            mu_eta_merged = muons_Eta[j];
                            mu_ch_merged = muonsCharge[j];
                            pt_rel_merged = ((muon_vect.Vect()).Cross(tot_vect_3.Vect())).Mag() / ((tot_vect_3.Vect()).Mag());

                            jet_pt_merged = jetpt[k];
                            jet_e_merged = tot_vect_3.E();
                            jet_phi_merged = jetphi[k];
                            jet_eta_merged = jeteta[k];

                            // mu_DB_merged = muons_DB[j];
                            // mu_DBerr_merged = muons_DBerr[j];
                            mu_Dz_merged = muons_Dz[j];
                            mu_Dxy_merged = muons_Dxy[j];
                            mu_Dxyerr_merged = muons_Dxyerr[j];
                            mu_MiniIso_merged = muons_MiniIso[j];
                            mu_Iso_merged = muons_Iso04[j];
                            mu_IsGlobal_merged = muons_IsGlobalMuon[j];
                            mu_IsTracker_merged = muons_IsTrackerMuon[j];
                            mu_NumberMatchedStations_merged = muons_NumberMatchedStations[j];
                            mu_NumberOfValidTrackerHits_merged = muons_NumberOfValidTrackerHits[j];
                            mu_isHigh_merged = muon_IsHighPtMuon[j];
                            mu_isTight_merged = muon_IsTightPtMuon[j];
                            mu_isLoose_merged = muon_IsLoosePtMuon[j];
                            mu_Dxy_fract_merged = muons_Dxy[j] / (muons_Dxyerr[j]);

                            met_pt_merged = metptCorr[0];
                            met_phi_merged = metptCorr[0];

                            for (int u = 0; u < genP_size; u++)
                            {
                                int index1 = genP_mom1I[u] - 1;
                                if (abs(genP_flavor[u]) == 13 && (genP_flavor[u] * muonsCharge[j] < 0.) && abs(genP_flavor[index1]) == 24)
                                {
                                    //                 if(abs(gP_flavor->at(u))==13 && (gP_flavor->at(u)*mu_ch->at(j)<0.) && abs(gP_flavor->at(index1))==24){
                                    muon_vect_gen.SetPtEtaPhiM(genP_Pt[u], genP_Eta[u], genP_Phi[u], genP_Mass[u]);
                                    float delta_R_true_corrected = sqrt(pow(genP_Eta[u] - muons_Eta[j], 2) + pow(muon_vect_gen.DeltaPhi(muon_vect), 2));
                                    if (delta_R_true_corrected < 0.1)
                                    {
                                        muon_high_truth[j] = 1;
                                    } // MC match
                                }     // part giusta
                            }         //loop generatore

                            for (int u = 0; u < genP_size; u++)
                            {
                                int index1 = genP_mom1I[u] - 1;
                                if (abs(genP_flavor[u]) == 15 && (genP_flavor[u] * muonsCharge[j] < 0.) && abs(genP_flavor[index1]) == 24)
                                {
                                    //                 if(abs(gP_flavor->at(u))==13 && (gP_flavor->at(u)*mu_ch->at(j)<0.) && abs(gP_flavor->at(index1))==24){
                                    muon_vect_gen.SetPtEtaPhiM(genP_Pt[u], genP_Eta[u], genP_Phi[u], genP_Mass[u]);
                                    float delta_R_true_corrected = sqrt(pow(genP_Eta[u] - muons_Eta[j], 2) + pow(muon_vect_gen.DeltaPhi(muon_vect), 2));
                                    if (delta_R_true_corrected < 0.05)
                                    {
                                        tau_high_truth[j] = 1;

                                    } // MC match
                                }     // part giusta
                            }         //loop generatore

                            jet_vect = jet_vect + muon_vect;
                            muon_vect.Boost(tot_vect_1.BoostVector());
                            jet_vect.Boost(tot_vect_1.BoostVector());

                            jet_boosted_pt_merged = (jet_vect.Pt());
                            jet_boosted_e_merged = (jet_vect.E());
                            jet_boosted_phi_merged = (jet_vect.Phi());
                            jet_boosted_eta_merged = (jet_vect.Eta());

                            mu_boosted_pt_merged = (muon_vect.Pt());
                            mu_boosted_e_merged = (muon_vect.E());
                            mu_boosted_phi_merged = (muon_vect.Phi());
                            mu_boosted_eta_merged = (muon_vect.Eta());
                            costheta_merged = (costhetap);

                            jet_vect.SetPtEtaPhiM(jetpt[k], jeteta[k], jetphi[k], jetmass[k]);
                            muon_vect.SetPtEtaPhiM(muons_Pt[j], muons_Eta[j], muons_Phi[j], muons_Mass[j]);

                            if ((muon_high_truth[j] == 1 || tau_high_truth[j] == 1) && jet_high_truth[k] == 1 && muonsCharge[j] * jetpartonflavour[k] > 0)
                            {
                                top_high_truth[j] = 1;
                                top_category_merged = 3;
                            }

                            if ((muon_high_truth[j] == 0 && tau_high_truth[j] == 0) && jet_high_truth[k] == 0)
                            {
                                top_category_merged = 0;
                            }

                            if ((muon_high_truth[j] == 0 && tau_high_truth[j] == 0) && jet_high_truth[k] == 1)
                            {
                                top_category_merged = 1;
                            }

                            if ((muon_high_truth[j] == 1 || tau_high_truth[j] == 1) && jet_high_truth[k] == 0)
                            {
                                top_category_merged = 2;
                            }

                            if ((muon_high_truth[j] == 1 || tau_high_truth[j] == 1) && jet_high_truth[k] == 1 && muonsCharge[j] * jetpartonflavour[k] < 0)
                            {
                                top_category_merged = 4;
                            }

                            mu_high_truth_merged = (muon_high_truth[j]);
                            top_high_truth_merged = (top_high_truth[j]);
                            tau_high_truth_merged = tau_high_truth[j];
                            counter_merged = counter_merged + 1;
                            nLepTop_merged = nLepTop;
                            is_merged->Fill();
                            jet_vect = jet_vect_tmp;
                        } // selezione delta merged
                    }     // muone giusto
                }         // muon loop
            }             // b tagged
        }                 //loop jet

        for (int u = 0; u < genP_size; u++)
        {
            for (int j = 0; j < muons_size; j++)
            {
                int index1 = genP_mom1I[u] - 1;
                if (muons_Pt[j] > 10)
                {
                    if (abs(genP_flavor[u]) == 13 && (genP_flavor[u] * muonsCharge[j] < 0.) && abs(genP_flavor[index1]) == 24)
                    {
                        muon_vect.SetPtEtaPhiM(muons_Pt[j], muons_Eta[j], muons_Phi[j], muons_Mass[j]);
                        muon_vect_gen.SetPtEtaPhiM(genP_Pt[u], genP_Eta[u], genP_Phi[u], genP_Mass[u]);

                        float delta_R_MC = sqrt(pow(genP_Eta[u] - muons_Eta[j], 2) + pow(muon_vect_gen.DeltaPhi(muon_vect), 2));
                        float delta_pT_frac = (muons_Pt[j] - genP_Pt[u]) / genP_Pt[u];

                        histo_DP_DR->Fill(delta_R_MC, delta_pT_frac);
                        histo_DP_DR->GetXaxis()->SetTitle("DR");
                        histo_DP_DR->GetYaxis()->SetTitle("DPt/Pt");
                    }
                }
            }
        }
    }

    //end of loop entries

    is_merged->Write();

    is_resolved->Write();

    histo_DP_DR->Write();

    f2->Close();

    return 1;
}