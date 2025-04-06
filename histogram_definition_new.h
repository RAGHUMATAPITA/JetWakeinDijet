#include "call_libraries.h"

TH1D* hEvents = new TH1D("hEvents", "", 10, 0, 10);
TH1D* hpthat = new TH1D("hpthat", "", 200, 0, 1200.);
TH1D* hpthatW = new TH1D("hpthatW", "", 100, 0, 1.);
TH1D* hCent = new TH1D("hCent", "", 200, 0, 200);
TH1D* hCent_MB = new TH1D("hCent_MB", "", 200, 0, 200);
TH1D* hZvtx = new TH1D("hZvtx", "", 200, -20, 20);

TH1D* hZvtx_ldsld = new TH1D("hZvtx_ldsld", "", 200, -20, 20);
TH1D* hZvtx_ldsld_Gen = new TH1D("hZvtx_ldsld_Gen", "", 200, -20, 20);

// for Centrality bin
//const int NCentbin = 5;
//double CentbinEdge[NCentbin+1] = {0., 20., 60., 100., 160., 200};
const int NCentbin = 6;
double CentbinEdge[NCentbin+1] = {0., 20., 60., 100., 140., 180., 200};
TH1D* hCentBin = new TH1D("hCentBin", "", NCentbin, CentbinEdge);

const int jtpT_nbins = 7;
double jtpT_bins_edge[jtpT_nbins+1] = {40., 120., 150., 190., 230., 300., 500., 5020.};
TH1D* hJetpTBin = new TH1D("hJetpTBin", "", jtpT_nbins, jtpT_bins_edge);

//const int trkpT_nbins = 5;
//double trkpT_bins_edge[trkpT_nbins+1] = {0.9, 2.0, 4.0, 10.0, 20.0, 5020.0};
const int trkpT_nbins = 3;
double trkpT_bins_edge[trkpT_nbins+1] = {0.9, 2.0, 4.0, 10.1};
TH1D* hTrkpTBin = new TH1D("hTrkpTBin", "", trkpT_nbins, trkpT_bins_edge);

// for JES and JER hist
int    bins4D_jer_jes[4]   =   { 400 ,  250  , NCentbin        , 8};
double xmin4D_jer_jes[4]   =   { 0. ,   0.   , 0               , 0};
double xmax4D_jer_jes[4]   =   { 4.  ,  5000., (double)NCentbin, 8};

THnSparseD * hJer_Jes_CorrpT_refpT_ctbin_flavour_W = new THnSparseD("hJer_Jes_CorrpT_refpT_ctbin_flavour_W", "", 4, bins4D_jer_jes, xmin4D_jer_jes, xmax4D_jer_jes);
THnSparseD * hJer_Jes_ldCorrpT_ldrefpT_ctbin_flavour_W = new THnSparseD("hJer_Jes_ldCorrpT_ldrefpT_ctbin_flavour_W", "", 4, bins4D_jer_jes, xmin4D_jer_jes, xmax4D_jer_jes);
THnSparseD * hJer_Jes_sldCorrpT_sldrefpT_ctbin_flavour_W = new THnSparseD("hJer_Jes_sldCorrpT_sldrefpT_ctbin_flavour_W", "", 4, bins4D_jer_jes, xmin4D_jer_jes, xmax4D_jer_jes);

// for jet histograms
int    bins4D_jet[5]   =   { 200  , 50   , 64           , NCentbin        , 8};
double xmin4D_jet[5]   =   { 0.   , -2.5 , -TMath::Pi() , 0.              , 0};
double xmax4D_jet[5]   =   { 1000.,  2.5 , TMath::Pi()  , (double)NCentbin, 8};

// reco/data/reco matched
// inclusive jet
THnSparseD * hJet_CorrpT_Eta_Phi_ctbin_flavour_pTCut_W = new THnSparseD("hJet_CorrpT_Eta_Phi_ctbin_flavour_pTCut_W", "", 5, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_Matched_CorrpT_Eta_Phi_ctbin_flavour_pTCut_W = new THnSparseD("hJet_Matched_CorrpT_Eta_Phi_ctbin_flavour_pTCut_W", "", 5, bins4D_jet, xmin4D_jet, xmax4D_jet);

// leading jet
THnSparseD * hJet_ldCorrpT_Eta_Phi_ctbin_flavour_pTCut_W = new THnSparseD("hJet_ldCorrpT_Eta_Phi_ctbin_flavour_pTCut_W", "", 5, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_Matched_ldCorrpT_Eta_Phi_ctbin_flavour_pTCut_W = new THnSparseD("hJet_Matched_ldCorrpT_Eta_Phi_ctbin_flavour_pTCut_W", "", 5, bins4D_jet, xmin4D_jet, xmax4D_jet);

// subleading jet
THnSparseD * hJet_sldCorrpT_Eta_Phi_ctbin_flavour_pTCut_W = new THnSparseD("hJet_sldCorrpT_Eta_Phi_ctbin_flavour_pTCut_W", "", 5, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD * hJet_Matched_sldCorrpT_Eta_Phi_ctbin_flavour_pTCut_W = new THnSparseD("hJet_Matched_sldCorrpT_Eta_Phi_ctbin_flavour_pTCut_W", "", 5, bins4D_jet, xmin4D_jet, xmax4D_jet);

// gen
// inclusive jet 
THnSparseD * hJet_GenpT_Eta_Phi_ctbin_flavour_pTCut_W = new THnSparseD("hJet_GenpT_Eta_Phi_ctbin_flavour_pTCut_W", "", 5, bins4D_jet, xmin4D_jet, xmax4D_jet);

// leading jet
THnSparseD * hJet_ldGenpT_Eta_Phi_ctbin_flavour_pTCut_W = new THnSparseD("hJet_ldGenpT_Eta_Phi_ctbin_flavour_pTCut_W", "", 5, bins4D_jet, xmin4D_jet, xmax4D_jet);

// subleading jet
THnSparseD * hJet_sldGenpT_Eta_Phi_ctbin_flavour_pTCut_W = new THnSparseD("hJet_sldGenpT_Eta_Phi_ctbin_flavour_pTCut_W", "", 5, bins4D_jet, xmin4D_jet, xmax4D_jet);

// for track histograms
int    bins4D_trk[4]   =   { 50  , 25   , 32           , NCentbin         };
double xmin4D_trk[4]   =   { 0.  , -2.5 , -TMath::Pi() , 0.               };
double xmax4D_trk[4]   =   { 20. ,  2.5 , TMath::Pi()  , (double)NCentbin };

// reco/data
THnSparseD * hTrk_pT_Eta_Phi_ctbin_pTCut_W = new THnSparseD("hTrk_pT_Eta_Phi_ctbin_pTCut_W", "", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
// corr pt
THnSparseD * hTrk_CorrpT_Eta_Phi_ctbin_pTCut_W = new THnSparseD("hTrk_CorrpT_Eta_Phi_ctbin_pTCut_W", "", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
// gen
THnSparseD * hTrk_GenpT_Eta_Phi_ctbin_pTCut_W = new THnSparseD("hTrk_GenpT_Eta_Phi_ctbin_pTCut_W", "", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);

// all sube
THnSparseD * hTrk_Signal_GenpT_Eta_Phi_ctbin_pTCut_W = new THnSparseD("hTrk_Signal_GenpT_Eta_Phi_ctbin_pTCut_W", "", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD * hTrk_Mixing_GenpT_Eta_Phi_ctbin_pTCut_W = new THnSparseD("hTrk_Mixing_GenpT_Eta_Phi_ctbin_pTCut_W", "", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);

// for sube = 0
THnSparseD * hTrk_sube0_Signal_GenpT_Eta_Phi_ctbin_pTCut_W = new THnSparseD("hTrk_sube0_Signal_GenpT_Eta_Phi_ctbin_pTCut_W", "", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD * hTrk_sube0_Mixing_GenpT_Eta_Phi_ctbin_pTCut_W = new THnSparseD("hTrk_sube0_Mixing_GenpT_Eta_Phi_ctbin_pTCut_W", "", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
// for sube = 1
THnSparseD * hTrk_sube1_Signal_GenpT_Eta_Phi_ctbin_pTCut_W = new THnSparseD("hTrk_sube1_Signal_GenpT_Eta_Phi_ctbin_pTCut_W", "", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD * hTrk_sube1_Mixing_GenpT_Eta_Phi_ctbin_pTCut_W = new THnSparseD("hTrk_sube1_Mixing_GenpT_Eta_Phi_ctbin_pTCut_W", "", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);

// for ntracks histograms
TH1D* hntrk_Signal_0 = new TH1D("hntrk_Signal_0", "", 2400, 0, 2400);
TH1D* hntrk_Signal_1 = new TH1D("hntrk_Signal_1", "", 2400, 0, 2400);
TH1D* hntrk_Signal_2 = new TH1D("hntrk_Signal_2", "", 2400, 0, 2400);
TH1D* hntrk_Signal_3 = new TH1D("hntrk_Signal_3", "", 2400, 0, 2400);
TH1D* hntrk_Signal_4 = new TH1D("hntrk_Signal_4", "", 2400, 0, 2400);

TH1D* hntrk_Signal_Check_0 = new TH1D("hntrk_Signal_Check_0", "", 2400, 0, 2400);
TH1D* hntrk_Signal_Check_1 = new TH1D("hntrk_Signal_Check_1", "", 2400, 0, 2400);
TH1D* hntrk_Signal_Check_2 = new TH1D("hntrk_Signal_Check_2", "", 2400, 0, 2400);
TH1D* hntrk_Signal_Check_3 = new TH1D("hntrk_Signal_Check_3", "", 2400, 0, 2400);
TH1D* hntrk_Signal_Check_4 = new TH1D("hntrk_Signal_Check_4", "", 2400, 0, 2400);

TH1D* hntrkoffline_Signal_0 = new TH1D("hntrkoffline_Signal_0", "", 2400, 0, 2400);
TH1D* hntrkoffline_Signal_1 = new TH1D("hntrkoffline_Signal_1", "", 2400, 0, 2400);
TH1D* hntrkoffline_Signal_2 = new TH1D("hntrkoffline_Signal_2", "", 2400, 0, 2400);
TH1D* hntrkoffline_Signal_3 = new TH1D("hntrkoffline_Signal_3", "", 2400, 0, 2400);
TH1D* hntrkoffline_Signal_4 = new TH1D("hntrkoffline_Signal_4", "", 2400, 0, 2400);

TH1D* hntrk_Mixing_0 = new TH1D("hntrk_Mixing_0", "", 2400, 0, 2400);
TH1D* hntrk_Mixing_1 = new TH1D("hntrk_Mixing_1", "", 2400, 0, 2400);
TH1D* hntrk_Mixing_2 = new TH1D("hntrk_Mixing_2", "", 2400, 0, 2400);
TH1D* hntrk_Mixing_3 = new TH1D("hntrk_Mixing_3", "", 2400, 0, 2400);
TH1D* hntrk_Mixing_4 = new TH1D("hntrk_Mixing_4", "", 2400, 0, 2400);

// for gen sube0
TH1D* hntrk_gen_sube0_Signal_0 = new TH1D("hntrk_gen_sube0_Signal_0", "", 2400, 0, 2400);
TH1D* hntrk_gen_sube0_Signal_1 = new TH1D("hntrk_gen_sube0_Signal_1", "", 2400, 0, 2400);
TH1D* hntrk_gen_sube0_Signal_2 = new TH1D("hntrk_gen_sube0_Signal_2", "", 2400, 0, 2400);
TH1D* hntrk_gen_sube0_Signal_3 = new TH1D("hntrk_gen_sube0_Signal_3", "", 2400, 0, 2400);
TH1D* hntrk_gen_sube0_Signal_4 = new TH1D("hntrk_gen_sube0_Signal_4", "", 2400, 0, 2400);

// for gen sube1
TH1D* hntrk_gen_sube1_Signal_0 = new TH1D("hntrk_gen_sube1_Signal_0", "", 2400, 0, 2400);
TH1D* hntrk_gen_sube1_Signal_1 = new TH1D("hntrk_gen_sube1_Signal_1", "", 2400, 0, 2400);
TH1D* hntrk_gen_sube1_Signal_2 = new TH1D("hntrk_gen_sube1_Signal_2", "", 2400, 0, 2400);
TH1D* hntrk_gen_sube1_Signal_3 = new TH1D("hntrk_gen_sube1_Signal_3", "", 2400, 0, 2400);
TH1D* hntrk_gen_sube1_Signal_4 = new TH1D("hntrk_gen_sube1_Signal_4", "", 2400, 0, 2400);

// for ntracks histograms
TH1D* hvtxz_Signal_0 = new TH1D("hvtxz_Signal_0", "", 200, -20, 20);
TH1D* hvtxz_Signal_1 = new TH1D("hvtxz_Signal_1", "", 200, -20, 20);
TH1D* hvtxz_Signal_2 = new TH1D("hvtxz_Signal_2", "", 200, -20, 20);
TH1D* hvtxz_Signal_3 = new TH1D("hvtxz_Signal_3", "", 200, -20, 20);
TH1D* hvtxz_Signal_4 = new TH1D("hvtxz_Signal_4", "", 200, -20, 20);

// For leading/subleading jet - tracks correlation
double etaW = (4.0*trk_eta_cut)/nBins_deta;
double phiW = (2.0*TMath::Pi())/nBins_dphi;

double eta_upedge = (2.0*trk_eta_cut)+(etaW/2.0);
double eta_lowedge = -(2.0*trk_eta_cut)-(etaW/2.0);

double phi_upedge = ((3.0*TMath::Pi()) - phiW)/2.0;
double phi_lowedge = -(TMath::Pi() - phiW)/2.0;

int    bins_corr[5]   =   { nBins_deta+1, nBins_dphi-1,  trkpT_nbins,         8,  NCentbin         };
double xmin_corr[5]   =   { eta_lowedge,  phi_lowedge,   0,                   0,  0.               };
double xmax_corr[5]   =   { eta_upedge,   phi_upedge,    (double)trkpT_nbins, 8,  (double)NCentbin };

int    bins_pair[3]   =   { 1, 8,  NCentbin         };
double xmin_pair[3]   =   { 0, 0,  0,               };
double xmax_pair[3]   =   { 1, 8,  (double)NCentbin };

// signal
// reco/data
// for jet pairs
THnSparseD* hldsld_Jet_pair_Signal_RapAsym = new THnSparseD("hldsld_Jet_pair_Signal_RapAsym", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_Signal_RapAsym_1 = new THnSparseD("hldsld_Jet_pair_Signal_RapAsym_1", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_Signal_RapAsym_2 = new THnSparseD("hldsld_Jet_pair_Signal_RapAsym_2", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_Signal_RapAsym_3 = new THnSparseD("hldsld_Jet_pair_Signal_RapAsym_3", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_Signal_RapAsym_sh = new THnSparseD("hldsld_Jet_pair_Signal_RapAsym_sh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_Signal_RapAsym_oh = new THnSparseD("hldsld_Jet_pair_Signal_RapAsym_oh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_Signal = new THnSparseD("hldsld_Jet_pair_Signal", "", 3, bins_pair, xmin_pair, xmax_pair);
// reco jet - gen tracks
// for sube = 0
THnSparseD* hldsld_Jet_pair_sube0_Signal_RapAsym = new THnSparseD("hldsld_Jet_pair_sube0_Signal_RapAsym", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube0_Signal_RapAsym_1 = new THnSparseD("hldsld_Jet_pair_sube0_Signal_RapAsym_1", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube0_Signal_RapAsym_2 = new THnSparseD("hldsld_Jet_pair_sube0_Signal_RapAsym_2", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube0_Signal_RapAsym_3 = new THnSparseD("hldsld_Jet_pair_sube0_Signal_RapAsym_3", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube0_Signal_RapAsym_sh = new THnSparseD("hldsld_Jet_pair_sube0_Signal_RapAsym_sh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube0_Signal_RapAsym_oh = new THnSparseD("hldsld_Jet_pair_sube0_Signal_RapAsym_oh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube0_Signal = new THnSparseD("hldsld_Jet_pair_sube0_Signal", "", 3, bins_pair, xmin_pair, xmax_pair);
// for sube = 1
THnSparseD* hldsld_Jet_pair_sube1_Signal_RapAsym = new THnSparseD("hldsld_Jet_pair_sube1_Signal_RapAsym", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube1_Signal_RapAsym_1 = new THnSparseD("hldsld_Jet_pair_sube1_Signal_RapAsym_1", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube1_Signal_RapAsym_2 = new THnSparseD("hldsld_Jet_pair_sube1_Signal_RapAsym_2", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube1_Signal_RapAsym_3 = new THnSparseD("hldsld_Jet_pair_sube1_Signal_RapAsym_3", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube1_Signal_RapAsym_sh = new THnSparseD("hldsld_Jet_pair_sube1_Signal_RapAsym_sh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube1_Signal_RapAsym_oh = new THnSparseD("hldsld_Jet_pair_sube1_Signal_RapAsym_oh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube1_Signal = new THnSparseD("hldsld_Jet_pair_sube1_Signal", "", 3, bins_pair, xmin_pair, xmax_pair);

//leading jet - trk
THnSparseD* hldJet_Trk_Signal_RapAsym_1 = new THnSparseD("hldJet_Trk_Signal_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_Signal_RapAsym_2 = new THnSparseD("hldJet_Trk_Signal_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_Signal_RapAsym_3 = new THnSparseD("hldJet_Trk_Signal_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_Signal_RapAsym_sh = new THnSparseD("hldJet_Trk_Signal_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_Signal_RapAsym_oh = new THnSparseD("hldJet_Trk_Signal_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_Signal = new THnSparseD("hldJet_Trk_Signal", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_Signal_RapAsym = new THnSparseD("hldJet_Trk_Signal_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
// reco jet - gen tracks
// for sube = 0
THnSparseD* hldJet_Trk_sube0_Signal_RapAsym = new THnSparseD("hldJet_Trk_sube0_Signal_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube0_Signal_RapAsym_1 = new THnSparseD("hldJet_Trk_sube0_Signal_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube0_Signal_RapAsym_2 = new THnSparseD("hldJet_Trk_sube0_Signal_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube0_Signal_RapAsym_3 = new THnSparseD("hldJet_Trk_sube0_Signal_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube0_Signal_RapAsym_sh = new THnSparseD("hldJet_Trk_sube0_Signal_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube0_Signal_RapAsym_oh = new THnSparseD("hldJet_Trk_sube0_Signal_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube0_Signal = new THnSparseD("hldJet_Trk_sube0_Signal", "", 5, bins_corr, xmin_corr, xmax_corr);
// for sube = 1
THnSparseD* hldJet_Trk_sube1_Signal_RapAsym = new THnSparseD("hldJet_Trk_sube1_Signal_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube1_Signal_RapAsym_1 = new THnSparseD("hldJet_Trk_sube1_Signal_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube1_Signal_RapAsym_2 = new THnSparseD("hldJet_Trk_sube1_Signal_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube1_Signal_RapAsym_3 = new THnSparseD("hldJet_Trk_sube1_Signal_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube1_Signal_RapAsym_sh = new THnSparseD("hldJet_Trk_sube1_Signal_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube1_Signal_RapAsym_oh = new THnSparseD("hldJet_Trk_sube1_Signal_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube1_Signal = new THnSparseD("hldJet_Trk_sube1_Signal", "", 5, bins_corr, xmin_corr, xmax_corr);

//subleading jet - trk
THnSparseD* hsldJet_Trk_Signal_RapAsym = new THnSparseD("hsldJet_Trk_Signal_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldJet_Trk_Signal_RapAsym_1 = new THnSparseD("hsldJet_Trk_Signal_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldJet_Trk_Signal_RapAsym_2 = new THnSparseD("hsldJet_Trk_Signal_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldJet_Trk_Signal_RapAsym_3 = new THnSparseD("hsldJet_Trk_Signal_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldJet_Trk_Signal_RapAsym_sh = new THnSparseD("hsldJet_Trk_Signal_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldJet_Trk_Signal_RapAsym_oh = new THnSparseD("hsldJet_Trk_Signal_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldJet_Trk_Signal = new THnSparseD("hsldJet_Trk_Signal", "", 5, bins_corr, xmin_corr, xmax_corr);

// gen
// for jet pairs
THnSparseD* hldsld_GenJet_pair_Signal_RapAsym = new THnSparseD("hldsld_GenJet_pair_Signal_RapAsym", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_Signal_RapAsym_1 = new THnSparseD("hldsld_GenJet_pair_Signal_RapAsym_1", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_Signal_RapAsym_2 = new THnSparseD("hldsld_GenJet_pair_Signal_RapAsym_2", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_Signal_RapAsym_3 = new THnSparseD("hldsld_GenJet_pair_Signal_RapAsym_3", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_Signal_RapAsym_sh = new THnSparseD("hldsld_GenJet_pair_Signal_RapAsym_sh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_Signal_RapAsym_oh = new THnSparseD("hldsld_GenJet_pair_Signal_RapAsym_oh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_Signal = new THnSparseD("hldsld_GenJet_pair_Signal", "", 3, bins_pair, xmin_pair, xmax_pair);
// for sube = 0
THnSparseD* hldsld_GenJet_pair_sube0_Signal_RapAsym = new THnSparseD("hldsld_GenJet_pair_sube0_Signal_RapAsym", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube0_Signal_RapAsym_1 = new THnSparseD("hldsld_GenJet_pair_sube0_Signal_RapAsym_1", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube0_Signal_RapAsym_2 = new THnSparseD("hldsld_GenJet_pair_sube0_Signal_RapAsym_2", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube0_Signal_RapAsym_3 = new THnSparseD("hldsld_GenJet_pair_sube0_Signal_RapAsym_3", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube0_Signal_RapAsym_sh = new THnSparseD("hldsld_GenJet_pair_sube0_Signal_RapAsym_sh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube0_Signal_RapAsym_oh = new THnSparseD("hldsld_GenJet_pair_sube0_Signal_RapAsym_oh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube0_Signal = new THnSparseD("hldsld_GenJet_pair_sube0_Signal", "", 3, bins_pair, xmin_pair, xmax_pair);
// for sube = 1
THnSparseD* hldsld_GenJet_pair_sube1_Signal_RapAsym = new THnSparseD("hldsld_GenJet_pair_sube1_Signal_RapAsym", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube1_Signal_RapAsym_1 = new THnSparseD("hldsld_GenJet_pair_sube1_Signal_RapAsym_1", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube1_Signal_RapAsym_2 = new THnSparseD("hldsld_GenJet_pair_sube1_Signal_RapAsym_2", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube1_Signal_RapAsym_3 = new THnSparseD("hldsld_GenJet_pair_sube1_Signal_RapAsym_3", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube1_Signal_RapAsym_sh = new THnSparseD("hldsld_GenJet_pair_sube1_Signal_RapAsym_sh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube1_Signal_RapAsym_oh = new THnSparseD("hldsld_GenJet_pair_sube1_Signal_RapAsym_oh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube1_Signal = new THnSparseD("hldsld_GenJet_pair_sube1_Signal", "", 3, bins_pair, xmin_pair, xmax_pair);

//leading jet - trk
THnSparseD* hldGenJet_Trk_Signal_RapAsym = new THnSparseD("hldGenJet_Trk_Signal_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_Signal_RapAsym_1 = new THnSparseD("hldGenJet_Trk_Signal_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_Signal_RapAsym_2 = new THnSparseD("hldGenJet_Trk_Signal_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_Signal_RapAsym_3 = new THnSparseD("hldGenJet_Trk_Signal_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_Signal_RapAsym_sh = new THnSparseD("hldGenJet_Trk_Signal_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_Signal_RapAsym_oh = new THnSparseD("hldGenJet_Trk_Signal_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_Signal = new THnSparseD("hldGenJet_Trk_Signal", "", 5, bins_corr, xmin_corr, xmax_corr);
// for sube =0
THnSparseD* hldGenJet_Trk_sube0_Signal_RapAsym = new THnSparseD("hldGenJet_Trk_sube0_Signal_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube0_Signal_RapAsym_1 = new THnSparseD("hldGenJet_Trk_sube0_Signal_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube0_Signal_RapAsym_2 = new THnSparseD("hldGenJet_Trk_sube0_Signal_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube0_Signal_RapAsym_3 = new THnSparseD("hldGenJet_Trk_sube0_Signal_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube0_Signal_RapAsym_sh = new THnSparseD("hldGenJet_Trk_sube0_Signal_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube0_Signal_RapAsym_oh = new THnSparseD("hldGenJet_Trk_sube0_Signal_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube0_Signal = new THnSparseD("hldGenJet_Trk_sube0_Signal", "", 5, bins_corr, xmin_corr, xmax_corr);
// for sube =1
THnSparseD* hldGenJet_Trk_sube1_Signal_RapAsym = new THnSparseD("hldGenJet_Trk_sube1_Signal_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube1_Signal_RapAsym_1 = new THnSparseD("hldGenJet_Trk_sube1_Signal_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube1_Signal_RapAsym_2 = new THnSparseD("hldGenJet_Trk_sube1_Signal_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube1_Signal_RapAsym_3 = new THnSparseD("hldGenJet_Trk_sube1_Signal_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube1_Signal_RapAsym_sh = new THnSparseD("hldGenJet_Trk_sube1_Signal_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube1_Signal_RapAsym_oh = new THnSparseD("hldGenJet_Trk_sube1_Signal_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube1_Signal = new THnSparseD("hldGenJet_Trk_sube1_Signal", "", 5, bins_corr, xmin_corr, xmax_corr);

//subleading jet - trk
THnSparseD* hsldGenJet_Trk_Signal_RapAsym = new THnSparseD("hsldGenJet_Trk_Signal_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_Signal_RapAsym_1 = new THnSparseD("hsldGenJet_Trk_Signal_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_Signal_RapAsym_2 = new THnSparseD("hsldGenJet_Trk_Signal_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_Signal_RapAsym_3 = new THnSparseD("hsldGenJet_Trk_Signal_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_Signal_RapAsym_sh = new THnSparseD("hsldGenJet_Trk_Signal_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_Signal_RapAsym_oh = new THnSparseD("hsldGenJet_Trk_Signal_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_Signal = new THnSparseD("hsldGenJet_Trk_Signal", "", 5, bins_corr, xmin_corr, xmax_corr);
// for sube =0
THnSparseD* hsldGenJet_Trk_sube0_Signal_RapAsym = new THnSparseD("hsldGenJet_Trk_sube0_Signal_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube0_Signal_RapAsym_1 = new THnSparseD("hsldGenJet_Trk_sube0_Signal_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube0_Signal_RapAsym_2 = new THnSparseD("hsldGenJet_Trk_sube0_Signal_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube0_Signal_RapAsym_3 = new THnSparseD("hsldGenJet_Trk_sube0_Signal_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube0_Signal_RapAsym_sh = new THnSparseD("hsldGenJet_Trk_sube0_Signal_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube0_Signal_RapAsym_oh = new THnSparseD("hsldGenJet_Trk_sube0_Signal_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube0_Signal = new THnSparseD("hsldGenJet_Trk_sube0_Signal", "", 5, bins_corr, xmin_corr, xmax_corr);
// for sube =1
THnSparseD* hsldGenJet_Trk_sube1_Signal_RapAsym = new THnSparseD("hsldGenJet_Trk_sube1_Signal_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube1_Signal_RapAsym_1 = new THnSparseD("hsldGenJet_Trk_sube1_Signal_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube1_Signal_RapAsym_2 = new THnSparseD("hsldGenJet_Trk_sube1_Signal_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube1_Signal_RapAsym_3 = new THnSparseD("hsldGenJet_Trk_sube1_Signal_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube1_Signal_RapAsym_sh = new THnSparseD("hsldGenJet_Trk_sube1_Signal_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube1_Signal_RapAsym_oh = new THnSparseD("hsldGenJet_Trk_sube1_Signal_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube1_Signal = new THnSparseD("hsldGenJet_Trk_sube1_Signal", "", 5, bins_corr, xmin_corr, xmax_corr);

// mixing
// for jet pairs
THnSparseD* hldsld_Jet_pair_Mixing_RapAsym = new THnSparseD("hldsld_Jet_pair_Mixing_RapAsym", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_Mixing_RapAsym_1 = new THnSparseD("hldsld_Jet_pair_Mixing_RapAsym_1", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_Mixing_RapAsym_2 = new THnSparseD("hldsld_Jet_pair_Mixing_RapAsym_2", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_Mixing_RapAsym_3 = new THnSparseD("hldsld_Jet_pair_Mixing_RapAsym_3", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_Mixing_RapAsym_sh = new THnSparseD("hldsld_Jet_pair_Mixing_RapAsym_sh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_Mixing_RapAsym_oh = new THnSparseD("hldsld_Jet_pair_Mixing_RapAsym_oh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_Mixing = new THnSparseD("hldsld_Jet_pair_Mixing", "", 3, bins_pair, xmin_pair, xmax_pair);
// reco jet - gen tracks
// for sube = 0
THnSparseD* hldsld_Jet_pair_sube0_Mixing_RapAsym = new THnSparseD("hldsld_Jet_pair_sube0_Mixing_RapAsym", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube0_Mixing_RapAsym_1 = new THnSparseD("hldsld_Jet_pair_sube0_Mixing_RapAsym_1", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube0_Mixing_RapAsym_2 = new THnSparseD("hldsld_Jet_pair_sube0_Mixing_RapAsym_2", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube0_Mixing_RapAsym_3 = new THnSparseD("hldsld_Jet_pair_sube0_Mixing_RapAsym_3", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube0_Mixing_RapAsym_sh = new THnSparseD("hldsld_Jet_pair_sube0_Mixing_RapAsym_sh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube0_Mixing_RapAsym_oh = new THnSparseD("hldsld_Jet_pair_sube0_Mixing_RapAsym_oh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube0_Mixing = new THnSparseD("hldsld_Jet_pair_sube0_Mixing", "", 3, bins_pair, xmin_pair, xmax_pair);
// for sube = 1
THnSparseD* hldsld_Jet_pair_sube1_Mixing_RapAsym = new THnSparseD("hldsld_Jet_pair_sube1_Mixing_RapAsym", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube1_Mixing_RapAsym_1 = new THnSparseD("hldsld_Jet_pair_sube1_Mixing_RapAsym_1", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube1_Mixing_RapAsym_2 = new THnSparseD("hldsld_Jet_pair_sube1_Mixing_RapAsym_2", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube1_Mixing_RapAsym_3 = new THnSparseD("hldsld_Jet_pair_sube1_Mixing_RapAsym_3", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube1_Mixing_RapAsym_sh = new THnSparseD("hldsld_Jet_pair_sube1_Mixing_RapAsym_sh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube1_Mixing_RapAsym_oh = new THnSparseD("hldsld_Jet_pair_sube1_Mixing_RapAsym_oh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_Jet_pair_sube1_Mixing = new THnSparseD("hldsld_Jet_pair_sube1_Mixing", "", 3, bins_pair, xmin_pair, xmax_pair);

// reco/data
//leading jet - trk
THnSparseD* hldJet_Trk_Mixing_RapAsym = new THnSparseD("hldJet_Trk_Mixing_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_Mixing_RapAsym_1 = new THnSparseD("hldJet_Trk_Mixing_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_Mixing_RapAsym_2 = new THnSparseD("hldJet_Trk_Mixing_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_Mixing_RapAsym_3 = new THnSparseD("hldJet_Trk_Mixing_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_Mixing_RapAsym_sh = new THnSparseD("hldJet_Trk_Mixing_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_Mixing_RapAsym_oh = new THnSparseD("hldJet_Trk_Mixing_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_Mixing = new THnSparseD("hldJet_Trk_Mixing", "", 5, bins_corr, xmin_corr, xmax_corr);
// reco jet - gen tracks
// for sube = 0
THnSparseD* hldJet_Trk_sube0_Mixing_RapAsym = new THnSparseD("hldJet_Trk_sube0_Mixing_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube0_Mixing_RapAsym_1 = new THnSparseD("hldJet_Trk_sube0_Mixing_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube0_Mixing_RapAsym_2 = new THnSparseD("hldJet_Trk_sube0_Mixing_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube0_Mixing_RapAsym_3 = new THnSparseD("hldJet_Trk_sube0_Mixing_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube0_Mixing_RapAsym_sh = new THnSparseD("hldJet_Trk_sube0_Mixing_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube0_Mixing_RapAsym_oh = new THnSparseD("hldJet_Trk_sube0_Mixing_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube0_Mixing = new THnSparseD("hldJet_Trk_sube0_Mixing", "", 5, bins_corr, xmin_corr, xmax_corr);
// for sube = 1
THnSparseD* hldJet_Trk_sube1_Mixing_RapAsym = new THnSparseD("hldJet_Trk_sube1_Mixing_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube1_Mixing_RapAsym_1 = new THnSparseD("hldJet_Trk_sube1_Mixing_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube1_Mixing_RapAsym_2 = new THnSparseD("hldJet_Trk_sube1_Mixing_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube1_Mixing_RapAsym_3 = new THnSparseD("hldJet_Trk_sube1_Mixing_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube1_Mixing_RapAsym_sh = new THnSparseD("hldJet_Trk_sube1_Mixing_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube1_Mixing_RapAsym_oh = new THnSparseD("hldJet_Trk_sube1_Mixing_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldJet_Trk_sube1_Mixing = new THnSparseD("hldJet_Trk_sube1_Mixing", "", 5, bins_corr, xmin_corr, xmax_corr);

//subleading jet - trk
THnSparseD* hsldJet_Trk_Mixing_RapAsym = new THnSparseD("hsldJet_Trk_Mixing_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldJet_Trk_Mixing_RapAsym_1 = new THnSparseD("hsldJet_Trk_Mixing_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldJet_Trk_Mixing_RapAsym_2 = new THnSparseD("hsldJet_Trk_Mixing_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldJet_Trk_Mixing_RapAsym_3 = new THnSparseD("hsldJet_Trk_Mixing_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldJet_Trk_Mixing_RapAsym_sh = new THnSparseD("hsldJet_Trk_Mixing_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldJet_Trk_Mixing_RapAsym_oh = new THnSparseD("hsldJet_Trk_Mixing_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldJet_Trk_Mixing = new THnSparseD("hsldJet_Trk_Mixing", "", 5, bins_corr, xmin_corr, xmax_corr);

// gen
// for jet pairs
THnSparseD* hldsld_GenJet_pair_Mixing_RapAsym = new THnSparseD("hldsld_GenJet_pair_Mixing_RapAsym", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_Mixing_RapAsym_1 = new THnSparseD("hldsld_GenJet_pair_Mixing_RapAsym_1", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_Mixing_RapAsym_2 = new THnSparseD("hldsld_GenJet_pair_Mixing_RapAsym_2", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_Mixing_RapAsym_3 = new THnSparseD("hldsld_GenJet_pair_Mixing_RapAsym_3", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_Mixing_RapAsym_sh = new THnSparseD("hldsld_GenJet_pair_Mixing_RapAsym_sh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_Mixing_RapAsym_oh = new THnSparseD("hldsld_GenJet_pair_Mixing_RapAsym_oh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_Mixing = new THnSparseD("hldsld_GenJet_pair_Mixing", "", 3, bins_pair, xmin_pair, xmax_pair);
// for sube = 0
THnSparseD* hldsld_GenJet_pair_sube0_Mixing_RapAsym = new THnSparseD("hldsld_GenJet_pair_sube0_Mixing_RapAsym", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube0_Mixing_RapAsym_1 = new THnSparseD("hldsld_GenJet_pair_sube0_Mixing_RapAsym_1", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube0_Mixing_RapAsym_2 = new THnSparseD("hldsld_GenJet_pair_sube0_Mixing_RapAsym_2", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube0_Mixing_RapAsym_3 = new THnSparseD("hldsld_GenJet_pair_sube0_Mixing_RapAsym_3", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube0_Mixing_RapAsym_sh = new THnSparseD("hldsld_GenJet_pair_sube0_Mixing_RapAsym_sh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube0_Mixing_RapAsym_oh = new THnSparseD("hldsld_GenJet_pair_sube0_Mixing_RapAsym_oh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube0_Mixing = new THnSparseD("hldsld_GenJet_pair_sube0_Mixing", "", 3, bins_pair, xmin_pair, xmax_pair);
// for sube = 1
THnSparseD* hldsld_GenJet_pair_sube1_Mixing_RapAsym = new THnSparseD("hldsld_GenJet_pair_sube1_Mixing_RapAsym", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube1_Mixing_RapAsym_1 = new THnSparseD("hldsld_GenJet_pair_sube1_Mixing_RapAsym_1", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube1_Mixing_RapAsym_2 = new THnSparseD("hldsld_GenJet_pair_sube1_Mixing_RapAsym_2", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube1_Mixing_RapAsym_3 = new THnSparseD("hldsld_GenJet_pair_sube1_Mixing_RapAsym_3", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube1_Mixing_RapAsym_sh = new THnSparseD("hldsld_GenJet_pair_sube1_Mixing_RapAsym_sh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube1_Mixing_RapAsym_oh = new THnSparseD("hldsld_GenJet_pair_sube1_Mixing_RapAsym_oh", "", 3, bins_pair, xmin_pair, xmax_pair);
THnSparseD* hldsld_GenJet_pair_sube1_Mixing = new THnSparseD("hldsld_GenJet_pair_sube1_Mixing", "", 3, bins_pair, xmin_pair, xmax_pair);

//leading jet - trk
THnSparseD* hldGenJet_Trk_Mixing_RapAsym = new THnSparseD("hldGenJet_Trk_Mixing_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_Mixing_RapAsym_1 = new THnSparseD("hldGenJet_Trk_Mixing_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_Mixing_RapAsym_2 = new THnSparseD("hldGenJet_Trk_Mixing_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_Mixing_RapAsym_3 = new THnSparseD("hldGenJet_Trk_Mixing_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_Mixing_RapAsym_sh = new THnSparseD("hldGenJet_Trk_Mixing_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_Mixing_RapAsym_oh = new THnSparseD("hldGenJet_Trk_Mixing_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_Mixing = new THnSparseD("hldGenJet_Trk_Mixing", "", 5, bins_corr, xmin_corr, xmax_corr);
// for sube =0
THnSparseD* hldGenJet_Trk_sube0_Mixing_RapAsym = new THnSparseD("hldGenJet_Trk_sube0_Mixing_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube0_Mixing_RapAsym_1 = new THnSparseD("hldGenJet_Trk_sube0_Mixing_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube0_Mixing_RapAsym_2 = new THnSparseD("hldGenJet_Trk_sube0_Mixing_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube0_Mixing_RapAsym_3 = new THnSparseD("hldGenJet_Trk_sube0_Mixing_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube0_Mixing_RapAsym_sh = new THnSparseD("hldGenJet_Trk_sube0_Mixing_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube0_Mixing_RapAsym_oh = new THnSparseD("hldGenJet_Trk_sube0_Mixing_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube0_Mixing = new THnSparseD("hldGenJet_Trk_sube0_Mixing", "", 5, bins_corr, xmin_corr, xmax_corr);
// for sube =1
THnSparseD* hldGenJet_Trk_sube1_Mixing_RapAsym = new THnSparseD("hldGenJet_Trk_sube1_Mixing_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube1_Mixing_RapAsym_1 = new THnSparseD("hldGenJet_Trk_sube1_Mixing_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube1_Mixing_RapAsym_2 = new THnSparseD("hldGenJet_Trk_sube1_Mixing_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube1_Mixing_RapAsym_3 = new THnSparseD("hldGenJet_Trk_sube1_Mixing_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube1_Mixing_RapAsym_sh = new THnSparseD("hldGenJet_Trk_sube1_Mixing_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube1_Mixing_RapAsym_oh = new THnSparseD("hldGenJet_Trk_sube1_Mixing_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hldGenJet_Trk_sube1_Mixing = new THnSparseD("hldGenJet_Trk_sube1_Mixing", "", 5, bins_corr, xmin_corr, xmax_corr);

//subleading jet - trk
THnSparseD* hsldGenJet_Trk_Mixing_RapAsym = new THnSparseD("hsldGenJet_Trk_Mixing_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_Mixing_RapAsym_1 = new THnSparseD("hsldGenJet_Trk_Mixing_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_Mixing_RapAsym_2 = new THnSparseD("hsldGenJet_Trk_Mixing_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_Mixing_RapAsym_3 = new THnSparseD("hsldGenJet_Trk_Mixing_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_Mixing_RapAsym_sh = new THnSparseD("hsldGenJet_Trk_Mixing_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_Mixing_RapAsym_oh = new THnSparseD("hsldGenJet_Trk_Mixing_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_Mixing = new THnSparseD("hsldGenJet_Trk_Mixing", "", 5, bins_corr, xmin_corr, xmax_corr);
// for sube =0
THnSparseD* hsldGenJet_Trk_sube0_Mixing_RapAsym = new THnSparseD("hsldGenJet_Trk_sube0_Mixing_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube0_Mixing_RapAsym_1 = new THnSparseD("hsldGenJet_Trk_sube0_Mixing_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube0_Mixing_RapAsym_2 = new THnSparseD("hsldGenJet_Trk_sube0_Mixing_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube0_Mixing_RapAsym_3 = new THnSparseD("hsldGenJet_Trk_sube0_Mixing_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube0_Mixing_RapAsym_sh = new THnSparseD("hsldGenJet_Trk_sube0_Mixing_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube0_Mixing_RapAsym_oh = new THnSparseD("hsldGenJet_Trk_sube0_Mixing_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube0_Mixing = new THnSparseD("hsldGenJet_Trk_sube0_Mixing", "", 5, bins_corr, xmin_corr, xmax_corr);
// for sube =1
THnSparseD* hsldGenJet_Trk_sube1_Mixing_RapAsym = new THnSparseD("hsldGenJet_Trk_sube1_Mixing_RapAsym", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube1_Mixing_RapAsym_1 = new THnSparseD("hsldGenJet_Trk_sube1_Mixing_RapAsym_1", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube1_Mixing_RapAsym_2 = new THnSparseD("hsldGenJet_Trk_sube1_Mixing_RapAsym_2", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube1_Mixing_RapAsym_3 = new THnSparseD("hsldGenJet_Trk_sube1_Mixing_RapAsym_3", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube1_Mixing_RapAsym_sh = new THnSparseD("hsldGenJet_Trk_sube1_Mixing_RapAsym_sh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube1_Mixing_RapAsym_oh = new THnSparseD("hsldGenJet_Trk_sube1_Mixing_RapAsym_oh", "", 5, bins_corr, xmin_corr, xmax_corr);
THnSparseD* hsldGenJet_Trk_sube1_Mixing = new THnSparseD("hsldGenJet_Trk_sube1_Mixing", "", 5, bins_corr, xmin_corr, xmax_corr);

void sumw2()
{
  //Event histograms
  hpthat->Sumw2();
  hpthatW->Sumw2();
  hEvents->Sumw2();
  hCent->Sumw2();
  hCent_MB->Sumw2();
  hZvtx->Sumw2();
  hZvtx_ldsld->Sumw2();
  hZvtx_ldsld_Gen->Sumw2();

  //~~~~~~~~~~~~~JES and JER histogram
  hJer_Jes_CorrpT_refpT_ctbin_flavour_W->Sumw2();
  hJer_Jes_ldCorrpT_ldrefpT_ctbin_flavour_W->Sumw2();
  hJer_Jes_sldCorrpT_sldrefpT_ctbin_flavour_W->Sumw2();
  
  //~~~~~~~~~~~~~reco/data/reco matched Jet histograms
  // inclusive jets
  hJet_CorrpT_Eta_Phi_ctbin_flavour_pTCut_W->Sumw2();
  hJet_Matched_CorrpT_Eta_Phi_ctbin_flavour_pTCut_W->Sumw2();
  
  // leading jets
  hJet_ldCorrpT_Eta_Phi_ctbin_flavour_pTCut_W->Sumw2();
  hJet_Matched_ldCorrpT_Eta_Phi_ctbin_flavour_pTCut_W->Sumw2();

  // subleading jets
  hJet_sldCorrpT_Eta_Phi_ctbin_flavour_pTCut_W->Sumw2();
  hJet_Matched_sldCorrpT_Eta_Phi_ctbin_flavour_pTCut_W->Sumw2();

  //~~~~~~~~~~~~~~gen Jet histograms
  // inclusive jets
  hJet_GenpT_Eta_Phi_ctbin_flavour_pTCut_W->Sumw2();

  // leading jets
  hJet_ldGenpT_Eta_Phi_ctbin_flavour_pTCut_W->Sumw2();

  // subleading jets
  hJet_sldGenpT_Eta_Phi_ctbin_flavour_pTCut_W->Sumw2();

  //~~~~~~~~~~~~~~~reco/data Trk histograms
  hTrk_pT_Eta_Phi_ctbin_pTCut_W->Sumw2();
  // corr pT
  hTrk_CorrpT_Eta_Phi_ctbin_pTCut_W->Sumw2();
  
  // gen Trk histograms
  hTrk_GenpT_Eta_Phi_ctbin_pTCut_W->Sumw2();

  // all sube
  hTrk_Signal_GenpT_Eta_Phi_ctbin_pTCut_W->Sumw2();
  hTrk_Mixing_GenpT_Eta_Phi_ctbin_pTCut_W->Sumw2();
  // for sube = 0
  hTrk_sube0_Signal_GenpT_Eta_Phi_ctbin_pTCut_W->Sumw2();
  hTrk_sube0_Mixing_GenpT_Eta_Phi_ctbin_pTCut_W->Sumw2();
  // for sube = 1
  hTrk_sube1_Signal_GenpT_Eta_Phi_ctbin_pTCut_W->Sumw2();
  hTrk_sube1_Mixing_GenpT_Eta_Phi_ctbin_pTCut_W->Sumw2();

  // ntrks histograms
  hntrk_Signal_0->Sumw2();
  hntrk_Signal_1->Sumw2();
  hntrk_Signal_2->Sumw2();
  hntrk_Signal_3->Sumw2();
  hntrk_Signal_4->Sumw2();


  hntrk_Signal_Check_0->Sumw2();
  hntrk_Signal_Check_1->Sumw2();
  hntrk_Signal_Check_2->Sumw2();
  hntrk_Signal_Check_3->Sumw2();
  hntrk_Signal_Check_4->Sumw2();

  hntrkoffline_Signal_0->Sumw2();
  hntrkoffline_Signal_1->Sumw2();
  hntrkoffline_Signal_2->Sumw2();
  hntrkoffline_Signal_3->Sumw2();
  hntrkoffline_Signal_4->Sumw2();

  hntrk_Mixing_0->Sumw2();
  hntrk_Mixing_1->Sumw2();
  hntrk_Mixing_2->Sumw2();
  hntrk_Mixing_3->Sumw2();
  hntrk_Mixing_4->Sumw2();

  // ntrks gen sube0 histograms
  hntrk_gen_sube0_Signal_0->Sumw2();
  hntrk_gen_sube0_Signal_1->Sumw2();
  hntrk_gen_sube0_Signal_2->Sumw2();
  hntrk_gen_sube0_Signal_3->Sumw2();
  hntrk_gen_sube0_Signal_4->Sumw2();

  // ntrks gen sube1 histograms
  hntrk_gen_sube1_Signal_0->Sumw2();
  hntrk_gen_sube1_Signal_1->Sumw2();
  hntrk_gen_sube1_Signal_2->Sumw2();
  hntrk_gen_sube1_Signal_3->Sumw2();
  hntrk_gen_sube1_Signal_4->Sumw2();
  
  // vz histograms
  hvtxz_Signal_0->Sumw2();
  hvtxz_Signal_1->Sumw2();
  hvtxz_Signal_2->Sumw2();
  hvtxz_Signal_3->Sumw2();
  hvtxz_Signal_4->Sumw2();
  //~~~~~~~~~~~~~~signal
  // for jet pairs
  hldsld_Jet_pair_Signal_RapAsym->Sumw2();
  hldsld_Jet_pair_Signal_RapAsym_1->Sumw2();
  hldsld_Jet_pair_Signal_RapAsym_2->Sumw2();
  hldsld_Jet_pair_Signal_RapAsym_3->Sumw2();
  hldsld_Jet_pair_Signal_RapAsym_sh->Sumw2();
  hldsld_Jet_pair_Signal_RapAsym_oh->Sumw2();
  hldsld_Jet_pair_Signal->Sumw2();
  // reco jet - gen tracks
  //for sube = 0
  hldsld_Jet_pair_sube0_Signal_RapAsym->Sumw2();
  hldsld_Jet_pair_sube0_Signal_RapAsym_1->Sumw2();
  hldsld_Jet_pair_sube0_Signal_RapAsym_2->Sumw2();
  hldsld_Jet_pair_sube0_Signal_RapAsym_3->Sumw2();
  hldsld_Jet_pair_sube0_Signal_RapAsym_sh->Sumw2();
  hldsld_Jet_pair_sube0_Signal_RapAsym_oh->Sumw2();
  hldsld_Jet_pair_sube0_Signal->Sumw2();
  //for sube = 1
  hldsld_Jet_pair_sube1_Signal_RapAsym->Sumw2();
  hldsld_Jet_pair_sube1_Signal_RapAsym_1->Sumw2();
  hldsld_Jet_pair_sube1_Signal_RapAsym_2->Sumw2();
  hldsld_Jet_pair_sube1_Signal_RapAsym_3->Sumw2();
  hldsld_Jet_pair_sube1_Signal_RapAsym_sh->Sumw2();
  hldsld_Jet_pair_sube1_Signal_RapAsym_oh->Sumw2();
  hldsld_Jet_pair_sube1_Signal->Sumw2();

  // reco/data correlation histograms
  // leading jet - trk
  hldJet_Trk_Signal_RapAsym->Sumw2();
  hldJet_Trk_Signal_RapAsym_1->Sumw2();
  hldJet_Trk_Signal_RapAsym_2->Sumw2();
  hldJet_Trk_Signal_RapAsym_3->Sumw2();
  hldJet_Trk_Signal_RapAsym_sh->Sumw2();
  hldJet_Trk_Signal_RapAsym_oh->Sumw2();
  hldJet_Trk_Signal->Sumw2();
  // reco jet - gen tracks
  //for sube = 0
  hldJet_Trk_sube0_Signal_RapAsym->Sumw2();
  hldJet_Trk_sube0_Signal_RapAsym_1->Sumw2();
  hldJet_Trk_sube0_Signal_RapAsym_2->Sumw2();
  hldJet_Trk_sube0_Signal_RapAsym_3->Sumw2();
  hldJet_Trk_sube0_Signal_RapAsym_sh->Sumw2();
  hldJet_Trk_sube0_Signal_RapAsym_oh->Sumw2();
  hldJet_Trk_sube0_Signal->Sumw2();
  //for sube = 1
  hldJet_Trk_sube1_Signal_RapAsym->Sumw2();
  hldJet_Trk_sube1_Signal_RapAsym_1->Sumw2();
  hldJet_Trk_sube1_Signal_RapAsym_2->Sumw2();
  hldJet_Trk_sube1_Signal_RapAsym_3->Sumw2();
  hldJet_Trk_sube1_Signal_RapAsym_sh->Sumw2();
  hldJet_Trk_sube1_Signal_RapAsym_oh->Sumw2();
  hldJet_Trk_sube1_Signal->Sumw2();

  // subleading jet - trk
  hsldJet_Trk_Signal_RapAsym->Sumw2();
  hsldJet_Trk_Signal_RapAsym_1->Sumw2();
  hsldJet_Trk_Signal_RapAsym_2->Sumw2();
  hsldJet_Trk_Signal_RapAsym_3->Sumw2();
  hsldJet_Trk_Signal_RapAsym_sh->Sumw2();
  hsldJet_Trk_Signal_RapAsym_oh->Sumw2();
  hsldJet_Trk_Signal->Sumw2();
  
  // gen correlation histograms
  // for jet pairs
  hldsld_GenJet_pair_Signal_RapAsym->Sumw2();
  hldsld_GenJet_pair_Signal_RapAsym_1->Sumw2();
  hldsld_GenJet_pair_Signal_RapAsym_2->Sumw2();
  hldsld_GenJet_pair_Signal_RapAsym_3->Sumw2();
  hldsld_GenJet_pair_Signal_RapAsym_sh->Sumw2();
  hldsld_GenJet_pair_Signal_RapAsym_oh->Sumw2();
  hldsld_GenJet_pair_Signal->Sumw2();
  // for sube = 0
  hldsld_GenJet_pair_sube0_Signal_RapAsym->Sumw2();
  hldsld_GenJet_pair_sube0_Signal_RapAsym_1->Sumw2();
  hldsld_GenJet_pair_sube0_Signal_RapAsym_2->Sumw2();
  hldsld_GenJet_pair_sube0_Signal_RapAsym_3->Sumw2();
  hldsld_GenJet_pair_sube0_Signal_RapAsym_sh->Sumw2();
  hldsld_GenJet_pair_sube0_Signal_RapAsym_oh->Sumw2();
  hldsld_GenJet_pair_sube0_Signal->Sumw2();
  // for sube = 1
  hldsld_GenJet_pair_sube1_Signal_RapAsym->Sumw2();
  hldsld_GenJet_pair_sube1_Signal_RapAsym_1->Sumw2();
  hldsld_GenJet_pair_sube1_Signal_RapAsym_2->Sumw2();
  hldsld_GenJet_pair_sube1_Signal_RapAsym_3->Sumw2();
  hldsld_GenJet_pair_sube1_Signal_RapAsym_sh->Sumw2();
  hldsld_GenJet_pair_sube1_Signal_RapAsym_oh->Sumw2();
  hldsld_GenJet_pair_sube1_Signal->Sumw2();

  // leading jet - trk
  hldGenJet_Trk_Signal_RapAsym->Sumw2();
  hldGenJet_Trk_Signal_RapAsym_1->Sumw2();
  hldGenJet_Trk_Signal_RapAsym_2->Sumw2();
  hldGenJet_Trk_Signal_RapAsym_3->Sumw2();
  hldGenJet_Trk_Signal_RapAsym_sh->Sumw2();
  hldGenJet_Trk_Signal_RapAsym_oh->Sumw2();
  hldGenJet_Trk_Signal->Sumw2();
  // for sube = 0
  hldGenJet_Trk_sube0_Signal_RapAsym->Sumw2();
  hldGenJet_Trk_sube0_Signal_RapAsym_1->Sumw2();
  hldGenJet_Trk_sube0_Signal_RapAsym_2->Sumw2();
  hldGenJet_Trk_sube0_Signal_RapAsym_3->Sumw2();
  hldGenJet_Trk_sube0_Signal_RapAsym_sh->Sumw2();
  hldGenJet_Trk_sube0_Signal_RapAsym_oh->Sumw2();
  hldGenJet_Trk_sube0_Signal->Sumw2();
  // for sube = 1
  hldGenJet_Trk_sube1_Signal_RapAsym->Sumw2();
  hldGenJet_Trk_sube1_Signal_RapAsym_1->Sumw2();
  hldGenJet_Trk_sube1_Signal_RapAsym_2->Sumw2();
  hldGenJet_Trk_sube1_Signal_RapAsym_3->Sumw2();
  hldGenJet_Trk_sube1_Signal_RapAsym_sh->Sumw2();
  hldGenJet_Trk_sube1_Signal_RapAsym_oh->Sumw2();
  hldGenJet_Trk_sube1_Signal->Sumw2();
  
  // subleading jet - trk
  hsldGenJet_Trk_Signal_RapAsym->Sumw2();
  hsldGenJet_Trk_Signal_RapAsym_1->Sumw2();
  hsldGenJet_Trk_Signal_RapAsym_2->Sumw2();
  hsldGenJet_Trk_Signal_RapAsym_3->Sumw2();
  hsldGenJet_Trk_Signal_RapAsym_sh->Sumw2();
  hsldGenJet_Trk_Signal_RapAsym_oh->Sumw2();
  hsldGenJet_Trk_Signal->Sumw2();
  // for sube = 0
  hsldGenJet_Trk_sube0_Signal_RapAsym->Sumw2();
  hsldGenJet_Trk_sube0_Signal_RapAsym_1->Sumw2();
  hsldGenJet_Trk_sube0_Signal_RapAsym_2->Sumw2();
  hsldGenJet_Trk_sube0_Signal_RapAsym_3->Sumw2();
  hsldGenJet_Trk_sube0_Signal_RapAsym_sh->Sumw2();
  hsldGenJet_Trk_sube0_Signal_RapAsym_oh->Sumw2();
  hsldGenJet_Trk_sube0_Signal->Sumw2();
  // for sube = 1
  hsldGenJet_Trk_sube1_Signal_RapAsym->Sumw2();
  hsldGenJet_Trk_sube1_Signal_RapAsym_1->Sumw2();
  hsldGenJet_Trk_sube1_Signal_RapAsym_2->Sumw2();
  hsldGenJet_Trk_sube1_Signal_RapAsym_3->Sumw2();
  hsldGenJet_Trk_sube1_Signal_RapAsym_sh->Sumw2();
  hsldGenJet_Trk_sube1_Signal_RapAsym_oh->Sumw2();
  hsldGenJet_Trk_sube1_Signal->Sumw2();
  
  //~~~~~~~~~~~~~mixing
    // for jet pairs
  hldsld_Jet_pair_Mixing_RapAsym->Sumw2();
  hldsld_Jet_pair_Mixing_RapAsym_1->Sumw2();
  hldsld_Jet_pair_Mixing_RapAsym_2->Sumw2();
  hldsld_Jet_pair_Mixing_RapAsym_3->Sumw2();
  hldsld_Jet_pair_Mixing_RapAsym_sh->Sumw2();
  hldsld_Jet_pair_Mixing_RapAsym_oh->Sumw2();
  hldsld_Jet_pair_Mixing->Sumw2();
  // reco jet - gen tracks
  //for sube = 0
  hldsld_Jet_pair_sube0_Mixing_RapAsym->Sumw2();
  hldsld_Jet_pair_sube0_Mixing_RapAsym_1->Sumw2();
  hldsld_Jet_pair_sube0_Mixing_RapAsym_2->Sumw2();
  hldsld_Jet_pair_sube0_Mixing_RapAsym_3->Sumw2();
  hldsld_Jet_pair_sube0_Mixing_RapAsym_sh->Sumw2();
  hldsld_Jet_pair_sube0_Mixing_RapAsym_oh->Sumw2();
  hldsld_Jet_pair_sube0_Mixing->Sumw2();
  //for sube = 1
  hldsld_Jet_pair_sube1_Mixing_RapAsym->Sumw2();
  hldsld_Jet_pair_sube1_Mixing_RapAsym_1->Sumw2();
  hldsld_Jet_pair_sube1_Mixing_RapAsym_2->Sumw2();
  hldsld_Jet_pair_sube1_Mixing_RapAsym_3->Sumw2();
  hldsld_Jet_pair_sube1_Mixing_RapAsym_sh->Sumw2();
  hldsld_Jet_pair_sube1_Mixing_RapAsym_oh->Sumw2();
  hldsld_Jet_pair_sube1_Mixing->Sumw2();

  // reco/data correlation histograms
  // leading jet - trk
  hldJet_Trk_Mixing_RapAsym->Sumw2();
  hldJet_Trk_Mixing_RapAsym_1->Sumw2();
  hldJet_Trk_Mixing_RapAsym_2->Sumw2();
  hldJet_Trk_Mixing_RapAsym_3->Sumw2();
  hldJet_Trk_Mixing_RapAsym_sh->Sumw2();
  hldJet_Trk_Mixing_RapAsym_oh->Sumw2();
  hldJet_Trk_Mixing->Sumw2();
  // reco jet - gen tracks
  //for sube = 0
  hldJet_Trk_sube0_Mixing_RapAsym->Sumw2();
  hldJet_Trk_sube0_Mixing_RapAsym_1->Sumw2();
  hldJet_Trk_sube0_Mixing_RapAsym_2->Sumw2();
  hldJet_Trk_sube0_Mixing_RapAsym_3->Sumw2();
  hldJet_Trk_sube0_Mixing_RapAsym_sh->Sumw2();
  hldJet_Trk_sube0_Mixing_RapAsym_oh->Sumw2();
  hldJet_Trk_sube0_Mixing->Sumw2();
  //for sube = 1
  hldJet_Trk_sube1_Mixing_RapAsym->Sumw2();
  hldJet_Trk_sube1_Mixing_RapAsym_1->Sumw2();
  hldJet_Trk_sube1_Mixing_RapAsym_2->Sumw2();
  hldJet_Trk_sube1_Mixing_RapAsym_3->Sumw2();
  hldJet_Trk_sube1_Mixing_RapAsym_sh->Sumw2();
  hldJet_Trk_sube1_Mixing_RapAsym_oh->Sumw2();
  hldJet_Trk_sube1_Mixing->Sumw2();
  
  // subleading jet - trk
  hsldJet_Trk_Mixing_RapAsym->Sumw2();
  hsldJet_Trk_Mixing_RapAsym_1->Sumw2();
  hsldJet_Trk_Mixing_RapAsym_2->Sumw2();
  hsldJet_Trk_Mixing_RapAsym_3->Sumw2();
  hsldJet_Trk_Mixing_RapAsym_sh->Sumw2();
  hsldJet_Trk_Mixing_RapAsym_oh->Sumw2();
  hsldJet_Trk_Mixing->Sumw2();
  
  // gen correlation histograms
  // for jet pairs
  hldsld_GenJet_pair_Mixing_RapAsym->Sumw2();
  hldsld_GenJet_pair_Mixing_RapAsym_1->Sumw2();
  hldsld_GenJet_pair_Mixing_RapAsym_2->Sumw2();
  hldsld_GenJet_pair_Mixing_RapAsym_3->Sumw2();
  hldsld_GenJet_pair_Mixing_RapAsym_sh->Sumw2();
  hldsld_GenJet_pair_Mixing_RapAsym_oh->Sumw2();
  hldsld_GenJet_pair_Mixing->Sumw2();
  // for sube = 0
  hldsld_GenJet_pair_sube0_Mixing_RapAsym->Sumw2();
  hldsld_GenJet_pair_sube0_Mixing_RapAsym_1->Sumw2();
  hldsld_GenJet_pair_sube0_Mixing_RapAsym_2->Sumw2();
  hldsld_GenJet_pair_sube0_Mixing_RapAsym_3->Sumw2();
  hldsld_GenJet_pair_sube0_Mixing_RapAsym_sh->Sumw2();
  hldsld_GenJet_pair_sube0_Mixing_RapAsym_oh->Sumw2();
  hldsld_GenJet_pair_sube0_Mixing->Sumw2();
  // for sube = 1
  hldsld_GenJet_pair_sube1_Mixing_RapAsym->Sumw2();
  hldsld_GenJet_pair_sube1_Mixing_RapAsym_1->Sumw2();
  hldsld_GenJet_pair_sube1_Mixing_RapAsym_2->Sumw2();
  hldsld_GenJet_pair_sube1_Mixing_RapAsym_3->Sumw2();
  hldsld_GenJet_pair_sube1_Mixing_RapAsym_sh->Sumw2();
  hldsld_GenJet_pair_sube1_Mixing_RapAsym_oh->Sumw2();
  hldsld_GenJet_pair_sube1_Mixing->Sumw2();
  
  // leading jet - trk
  hldGenJet_Trk_Mixing_RapAsym->Sumw2();
  hldGenJet_Trk_Mixing_RapAsym_1->Sumw2();
  hldGenJet_Trk_Mixing_RapAsym_2->Sumw2();
  hldGenJet_Trk_Mixing_RapAsym_3->Sumw2();
  hldGenJet_Trk_Mixing_RapAsym_sh->Sumw2();
  hldGenJet_Trk_Mixing_RapAsym_oh->Sumw2();
  hldGenJet_Trk_Mixing->Sumw2();
  // for sube = 0
  hldGenJet_Trk_sube0_Mixing_RapAsym->Sumw2();
  hldGenJet_Trk_sube0_Mixing_RapAsym_1->Sumw2();
  hldGenJet_Trk_sube0_Mixing_RapAsym_2->Sumw2();
  hldGenJet_Trk_sube0_Mixing_RapAsym_3->Sumw2();
  hldGenJet_Trk_sube0_Mixing_RapAsym_sh->Sumw2();
  hldGenJet_Trk_sube0_Mixing_RapAsym_oh->Sumw2();
  hldGenJet_Trk_sube0_Mixing->Sumw2();
  // for sube = 1
  hldGenJet_Trk_sube1_Mixing_RapAsym->Sumw2();
  hldGenJet_Trk_sube1_Mixing_RapAsym_1->Sumw2();
  hldGenJet_Trk_sube1_Mixing_RapAsym_2->Sumw2();
  hldGenJet_Trk_sube1_Mixing_RapAsym_3->Sumw2();
  hldGenJet_Trk_sube1_Mixing_RapAsym_sh->Sumw2();
  hldGenJet_Trk_sube1_Mixing_RapAsym_oh->Sumw2();
  hldGenJet_Trk_sube1_Mixing->Sumw2();
  
  // subleading jet - trk
  hsldGenJet_Trk_Mixing_RapAsym->Sumw2();
  hsldGenJet_Trk_Mixing_RapAsym_1->Sumw2();
  hsldGenJet_Trk_Mixing_RapAsym_2->Sumw2();
  hsldGenJet_Trk_Mixing_RapAsym_3->Sumw2();
  hsldGenJet_Trk_Mixing_RapAsym_sh->Sumw2();
  hsldGenJet_Trk_Mixing_RapAsym_oh->Sumw2();
  hsldGenJet_Trk_Mixing->Sumw2();
  // for sube = 0
  hsldGenJet_Trk_sube0_Mixing_RapAsym->Sumw2();
  hsldGenJet_Trk_sube0_Mixing_RapAsym_1->Sumw2();
  hsldGenJet_Trk_sube0_Mixing_RapAsym_2->Sumw2();
  hsldGenJet_Trk_sube0_Mixing_RapAsym_3->Sumw2();
  hsldGenJet_Trk_sube0_Mixing_RapAsym_sh->Sumw2();
  hsldGenJet_Trk_sube0_Mixing_RapAsym_oh->Sumw2();
  hsldGenJet_Trk_sube0_Mixing->Sumw2();
  // for sube = 1
  hsldGenJet_Trk_sube1_Mixing_RapAsym->Sumw2();
  hsldGenJet_Trk_sube1_Mixing_RapAsym_1->Sumw2();
  hsldGenJet_Trk_sube1_Mixing_RapAsym_2->Sumw2();
  hsldGenJet_Trk_sube1_Mixing_RapAsym_3->Sumw2();
  hsldGenJet_Trk_sube1_Mixing_RapAsym_sh->Sumw2();
  hsldGenJet_Trk_sube1_Mixing_RapAsym_oh->Sumw2();
  hsldGenJet_Trk_sube1_Mixing->Sumw2();

}

void Write_Event_hist(const bool& is_MC)
{
  hEvents->Write();
  if(is_MC)
    {
      hpthat->Write();
      hpthatW->Write();
    }
  hCent->Write();
  hCent_MB->Write();
  hZvtx->Write();
  hZvtx_ldsld->Write();
  if(is_MC)
    {
      hZvtx_ldsld_Gen->Write();
    }
}

void Write_Jet_QA_hist(const bool& is_MC, const bool& is_JES_JER)
{
  // reco/data/reco matched Jet histograms
  // inclusive jets
  hJet_CorrpT_Eta_Phi_ctbin_flavour_pTCut_W->Write();
  // leading jets
  hJet_ldCorrpT_Eta_Phi_ctbin_flavour_pTCut_W->Write();
  // subleading jets
  hJet_sldCorrpT_Eta_Phi_ctbin_flavour_pTCut_W->Write();
  
  if(is_MC)
    {
      // inclusive jets
      hJet_Matched_CorrpT_Eta_Phi_ctbin_flavour_pTCut_W->Write();
      // leading jets
      hJet_Matched_ldCorrpT_Eta_Phi_ctbin_flavour_pTCut_W->Write();
      // subleading jets
      hJet_Matched_sldCorrpT_Eta_Phi_ctbin_flavour_pTCut_W->Write();
      
      // gen histograms
      // inclusive jets
      hJet_GenpT_Eta_Phi_ctbin_flavour_pTCut_W->Write();
      // leading jets
      hJet_ldGenpT_Eta_Phi_ctbin_flavour_pTCut_W->Write();
      // subleading jets
      hJet_sldGenpT_Eta_Phi_ctbin_flavour_pTCut_W->Write();

      if(is_JES_JER)
	{
	  hJer_Jes_CorrpT_refpT_ctbin_flavour_W->Write();
	  hJer_Jes_ldCorrpT_ldrefpT_ctbin_flavour_W->Write();
	  hJer_Jes_sldCorrpT_sldrefpT_ctbin_flavour_W->Write();
	}
    }
}

void Write_Trk_QA_hist(const bool& is_MC)
{
  /*
  // reco/data Trk histograms
  hTrk_pT_Eta_Phi_ctbin_pTCut_W->Write();
  //corr pT
  hTrk_CorrpT_Eta_Phi_ctbin_pTCut_W->Write();
  
  hTrk_Signal_GenpT_Eta_Phi_ctbin_pTCut_W->Write();
  hTrk_Mixing_GenpT_Eta_Phi_ctbin_pTCut_W->Write();
  */
  
  if(is_MC)
    {
      /*
      // gen Trk histograms
      hTrk_GenpT_Eta_Phi_ctbin_pTCut_W->Write();
      */

      /*
      // for sube = 0
      hTrk_sube0_Signal_GenpT_Eta_Phi_ctbin_pTCut_W->Write();
      hTrk_sube0_Mixing_GenpT_Eta_Phi_ctbin_pTCut_W->Write();
      // for sube = 1
      hTrk_sube1_Signal_GenpT_Eta_Phi_ctbin_pTCut_W->Write();
      hTrk_sube1_Mixing_GenpT_Eta_Phi_ctbin_pTCut_W->Write();
      */
      
      // ntrks gen sube0 histograms
      hntrk_gen_sube0_Signal_0->Write();
      hntrk_gen_sube0_Signal_1->Write();
      hntrk_gen_sube0_Signal_2->Write();
      hntrk_gen_sube0_Signal_3->Write();
      hntrk_gen_sube0_Signal_4->Write();
      
      // ntrks gen sube1 histograms
      hntrk_gen_sube1_Signal_0->Write();
      hntrk_gen_sube1_Signal_1->Write();
      hntrk_gen_sube1_Signal_2->Write();
      hntrk_gen_sube1_Signal_3->Write();
      hntrk_gen_sube1_Signal_4->Write();
    }
  
  // ntrks histograms
  hntrk_Signal_0->Write();
  hntrk_Signal_1->Write();
  hntrk_Signal_2->Write();
  hntrk_Signal_3->Write();
  hntrk_Signal_4->Write();

  hntrk_Signal_Check_0->Write();
  hntrk_Signal_Check_1->Write();
  hntrk_Signal_Check_2->Write();
  hntrk_Signal_Check_3->Write();
  hntrk_Signal_Check_4->Write();

  hntrkoffline_Signal_0->Write();
  hntrkoffline_Signal_1->Write();
  hntrkoffline_Signal_2->Write();
  hntrkoffline_Signal_3->Write();
  hntrkoffline_Signal_4->Write();
  
  hntrk_Mixing_0->Write();
  hntrk_Mixing_1->Write();
  hntrk_Mixing_2->Write();
  hntrk_Mixing_3->Write();
  hntrk_Mixing_4->Write();
  
  hvtxz_Signal_0->Write();
  hvtxz_Signal_1->Write();
  hvtxz_Signal_2->Write();
  hvtxz_Signal_3->Write();
  hvtxz_Signal_4->Write();
}

void Write_Jet_Trk_Corr_hist(const bool& is_MC, const bool& isRcJetGnTrk)
{
  // signal
  // for jet pairs
  hldsld_Jet_pair_Signal_RapAsym->Write();
  hldsld_Jet_pair_Signal_RapAsym_1->Write();
  hldsld_Jet_pair_Signal_RapAsym_2->Write();
  hldsld_Jet_pair_Signal_RapAsym_3->Write();
  hldsld_Jet_pair_Signal_RapAsym_sh->Write();
  hldsld_Jet_pair_Signal_RapAsym_oh->Write();
  hldsld_Jet_pair_Signal->Write();

  if(isRcJetGnTrk)
    {
      // reco jet - gen tracks
      //for sube = 0
      hldsld_Jet_pair_sube0_Signal_RapAsym->Write();
      hldsld_Jet_pair_sube0_Signal_RapAsym_1->Write();
      hldsld_Jet_pair_sube0_Signal_RapAsym_2->Write();
      hldsld_Jet_pair_sube0_Signal_RapAsym_3->Write();
      hldsld_Jet_pair_sube0_Signal_RapAsym_sh->Write();
      hldsld_Jet_pair_sube0_Signal_RapAsym_oh->Write();
      hldsld_Jet_pair_sube0_Signal->Write();
      //for sube = 1
      hldsld_Jet_pair_sube1_Signal_RapAsym->Write();
      hldsld_Jet_pair_sube1_Signal_RapAsym_1->Write();
      hldsld_Jet_pair_sube1_Signal_RapAsym_2->Write();
      hldsld_Jet_pair_sube1_Signal_RapAsym_3->Write();
      hldsld_Jet_pair_sube1_Signal_RapAsym_sh->Write();
      hldsld_Jet_pair_sube1_Signal_RapAsym_oh->Write();
      hldsld_Jet_pair_sube1_Signal->Write();
    }
  
  // reco/data correlation histograms
  // leading jet - trk
  hldJet_Trk_Signal_RapAsym->Write();
  hldJet_Trk_Signal_RapAsym_1->Write();
  hldJet_Trk_Signal_RapAsym_2->Write();
  hldJet_Trk_Signal_RapAsym_3->Write();
  hldJet_Trk_Signal_RapAsym_sh->Write();
  hldJet_Trk_Signal_RapAsym_oh->Write();
  hldJet_Trk_Signal->Write();
  
  if(isRcJetGnTrk)
    {
      // reco jet - gen tracks
      //for sube = 0
      hldJet_Trk_sube0_Signal_RapAsym->Write();
      hldJet_Trk_sube0_Signal_RapAsym_1->Write();
      hldJet_Trk_sube0_Signal_RapAsym_2->Write();
      hldJet_Trk_sube0_Signal_RapAsym_3->Write();
      hldJet_Trk_sube0_Signal_RapAsym_sh->Write();
      hldJet_Trk_sube0_Signal_RapAsym_oh->Write();
      hldJet_Trk_sube0_Signal->Write();
      //for sube = 1
      hldJet_Trk_sube1_Signal_RapAsym->Write();
      hldJet_Trk_sube1_Signal_RapAsym_1->Write();
      hldJet_Trk_sube1_Signal_RapAsym_2->Write();
      hldJet_Trk_sube1_Signal_RapAsym_3->Write();
      hldJet_Trk_sube1_Signal_RapAsym_sh->Write();
      hldJet_Trk_sube1_Signal_RapAsym_oh->Write();
      hldJet_Trk_sube1_Signal->Write();
    }
  
  /*
  // subleading jet - trk
  hsldJet_Trk_Signal_RapAsym->Write();
  hsldJet_Trk_Signal_RapAsym_1->Write();
  hsldJet_Trk_Signal_RapAsym_2->Write();
  hsldJet_Trk_Signal_RapAsym_3->Write();
  hsldJet_Trk_Signal_RapAsym_sh->Write();
  hsldJet_Trk_Signal_RapAsym_oh->Write();
  hsldJet_Trk_Signal->Write();
  */
  
  // mixing
  // for jet pairs
  hldsld_Jet_pair_Mixing_RapAsym->Write();
  hldsld_Jet_pair_Mixing_RapAsym_1->Write();
  hldsld_Jet_pair_Mixing_RapAsym_2->Write();
  hldsld_Jet_pair_Mixing_RapAsym_3->Write();
  hldsld_Jet_pair_Mixing_RapAsym_sh->Write();
  hldsld_Jet_pair_Mixing_RapAsym_oh->Write();
  hldsld_Jet_pair_Mixing->Write();

  if(isRcJetGnTrk)
    {
      // reco jet - gen tracks
      //for sube = 0
      hldsld_Jet_pair_sube0_Mixing_RapAsym->Write();
      hldsld_Jet_pair_sube0_Mixing_RapAsym_1->Write();
      hldsld_Jet_pair_sube0_Mixing_RapAsym_2->Write();
      hldsld_Jet_pair_sube0_Mixing_RapAsym_3->Write();
      hldsld_Jet_pair_sube0_Mixing_RapAsym_sh->Write();
      hldsld_Jet_pair_sube0_Mixing_RapAsym_oh->Write();
      hldsld_Jet_pair_sube0_Mixing->Write();
      //for sube = 1
      hldsld_Jet_pair_sube1_Mixing_RapAsym->Write();
      hldsld_Jet_pair_sube1_Mixing_RapAsym_1->Write();
      hldsld_Jet_pair_sube1_Mixing_RapAsym_2->Write();
      hldsld_Jet_pair_sube1_Mixing_RapAsym_3->Write();
      hldsld_Jet_pair_sube1_Mixing_RapAsym_sh->Write();
      hldsld_Jet_pair_sube1_Mixing_RapAsym_oh->Write();
      hldsld_Jet_pair_sube1_Mixing->Write();
    }
  
  // reco/data correlation histograms
  // leading jet - trk
  hldJet_Trk_Mixing_RapAsym->Write();
  hldJet_Trk_Mixing_RapAsym_1->Write();
  hldJet_Trk_Mixing_RapAsym_2->Write();
  hldJet_Trk_Mixing_RapAsym_3->Write();
  hldJet_Trk_Mixing_RapAsym_sh->Write();
  hldJet_Trk_Mixing_RapAsym_oh->Write();
  hldJet_Trk_Mixing->Write();

  if(isRcJetGnTrk)
    {
      // reco jet - gen tracks
      //for sube = 0
      hldJet_Trk_sube0_Mixing_RapAsym->Write();
      hldJet_Trk_sube0_Mixing_RapAsym_1->Write();
      hldJet_Trk_sube0_Mixing_RapAsym_2->Write();
      hldJet_Trk_sube0_Mixing_RapAsym_3->Write();
      hldJet_Trk_sube0_Mixing_RapAsym_sh->Write();
      hldJet_Trk_sube0_Mixing_RapAsym_oh->Write();
      hldJet_Trk_sube0_Mixing->Write();
      //for sube = 1
      hldJet_Trk_sube1_Mixing_RapAsym->Write();
      hldJet_Trk_sube1_Mixing_RapAsym_1->Write();
      hldJet_Trk_sube1_Mixing_RapAsym_2->Write();
      hldJet_Trk_sube1_Mixing_RapAsym_3->Write();
      hldJet_Trk_sube1_Mixing_RapAsym_sh->Write();
      hldJet_Trk_sube1_Mixing_RapAsym_oh->Write();
      hldJet_Trk_sube1_Mixing->Write();
    }
  
  /*
  // subleading jet - trk
  hsldJet_Trk_Mixing_RapAsym->Write();
  hsldJet_Trk_Mixing_RapAsym_1->Write();
  hsldJet_Trk_Mixing_RapAsym_2->Write();
  hsldJet_Trk_Mixing_RapAsym_3->Write();
  hsldJet_Trk_Mixing_RapAsym_sh->Write();
  hsldJet_Trk_Mixing_RapAsym_oh->Write();
  hsldJet_Trk_Mixing->Write();
  */
  
  if(is_MC)
    {
      // signal
      // for jet pairs
      hldsld_GenJet_pair_Signal_RapAsym->Write();
      hldsld_GenJet_pair_Signal_RapAsym_1->Write();
      hldsld_GenJet_pair_Signal_RapAsym_2->Write();
      hldsld_GenJet_pair_Signal_RapAsym_3->Write();
      hldsld_GenJet_pair_Signal_RapAsym_sh->Write();
      hldsld_GenJet_pair_Signal_RapAsym_oh->Write();
      hldsld_GenJet_pair_Signal->Write();
      
      if(!isRcJetGnTrk)
	{
	  // for sube = 0
	  hldsld_GenJet_pair_sube0_Signal_RapAsym->Write();
	  hldsld_GenJet_pair_sube0_Signal_RapAsym_1->Write();
	  hldsld_GenJet_pair_sube0_Signal_RapAsym_2->Write();
	  hldsld_GenJet_pair_sube0_Signal_RapAsym_3->Write();
	  hldsld_GenJet_pair_sube0_Signal_RapAsym_sh->Write();
	  hldsld_GenJet_pair_sube0_Signal_RapAsym_oh->Write();
	  hldsld_GenJet_pair_sube0_Signal->Write();
	  // for sube = 1
	  hldsld_GenJet_pair_sube1_Signal_RapAsym->Write();
	  hldsld_GenJet_pair_sube1_Signal_RapAsym_1->Write();
	  hldsld_GenJet_pair_sube1_Signal_RapAsym_2->Write();
	  hldsld_GenJet_pair_sube1_Signal_RapAsym_3->Write();
	  hldsld_GenJet_pair_sube1_Signal_RapAsym_sh->Write();
	  hldsld_GenJet_pair_sube1_Signal_RapAsym_oh->Write();
	  hldsld_GenJet_pair_sube1_Signal->Write();
	}
      
      // gen correlation histograms
      // leading jet - trk
      hldGenJet_Trk_Signal_RapAsym->Write();
      hldGenJet_Trk_Signal_RapAsym_1->Write();
      hldGenJet_Trk_Signal_RapAsym_2->Write();
      hldGenJet_Trk_Signal_RapAsym_3->Write();
      hldGenJet_Trk_Signal_RapAsym_sh->Write();
      hldGenJet_Trk_Signal_RapAsym_oh->Write();
      hldGenJet_Trk_Signal->Write();
      
      if(!isRcJetGnTrk)
	{
	  // for sube = 0
	  hldGenJet_Trk_sube0_Signal_RapAsym->Write();
	  hldGenJet_Trk_sube0_Signal_RapAsym_1->Write();
	  hldGenJet_Trk_sube0_Signal_RapAsym_2->Write();
	  hldGenJet_Trk_sube0_Signal_RapAsym_3->Write();
	  hldGenJet_Trk_sube0_Signal_RapAsym_sh->Write();
	  hldGenJet_Trk_sube0_Signal_RapAsym_oh->Write();
	  hldGenJet_Trk_sube0_Signal->Write();
	  // for sube = 1
	  hldGenJet_Trk_sube1_Signal_RapAsym->Write();
	  hldGenJet_Trk_sube1_Signal_RapAsym_1->Write();
	  hldGenJet_Trk_sube1_Signal_RapAsym_2->Write();
	  hldGenJet_Trk_sube1_Signal_RapAsym_3->Write();
	  hldGenJet_Trk_sube1_Signal_RapAsym_sh->Write();
	  hldGenJet_Trk_sube1_Signal_RapAsym_oh->Write();
	  hldGenJet_Trk_sube1_Signal->Write();
	}
      /*
      // subleading jet - trk
      hsldGenJet_Trk_Signal_RapAsym->Write();
      hsldGenJet_Trk_Signal_RapAsym_1->Write();
      hsldGenJet_Trk_Signal_RapAsym_2->Write();
      hsldGenJet_Trk_Signal_RapAsym_3->Write();
      hsldGenJet_Trk_Signal_RapAsym_sh->Write();
      hsldGenJet_Trk_Signal_RapAsym_oh->Write();
      hsldGenJet_Trk_Signal->Write();
      // for sube = 0
      hsldGenJet_Trk_sube0_Signal_RapAsym->Write();
      hsldGenJet_Trk_sube0_Signal_RapAsym_1->Write();
      hsldGenJet_Trk_sube0_Signal_RapAsym_2->Write();
      hsldGenJet_Trk_sube0_Signal_RapAsym_3->Write();
      hsldGenJet_Trk_sube0_Signal_RapAsym_sh->Write();
      hsldGenJet_Trk_sube0_Signal_RapAsym_oh->Write();
      hsldGenJet_Trk_sube0_Signal->Write();
      // for sube = 1
      hsldGenJet_Trk_sube1_Signal_RapAsym->Write();
      hsldGenJet_Trk_sube1_Signal_RapAsym_1->Write();
      hsldGenJet_Trk_sube1_Signal_RapAsym_2->Write();
      hsldGenJet_Trk_sube1_Signal_RapAsym_3->Write();
      hsldGenJet_Trk_sube1_Signal_RapAsym_sh->Write();
      hsldGenJet_Trk_sube1_Signal_RapAsym_oh->Write();
      hsldGenJet_Trk_sube1_Signal->Write();
      */
      
      // mixing
      // for jet pairs
      hldsld_GenJet_pair_Mixing_RapAsym->Write();
      hldsld_GenJet_pair_Mixing_RapAsym_1->Write();
      hldsld_GenJet_pair_Mixing_RapAsym_2->Write();
      hldsld_GenJet_pair_Mixing_RapAsym_3->Write();
      hldsld_GenJet_pair_Mixing_RapAsym_sh->Write();
      hldsld_GenJet_pair_Mixing_RapAsym_oh->Write();
      hldsld_GenJet_pair_Mixing->Write();

      if(!isRcJetGnTrk)
	{
	  // for sube = 0
	  hldsld_GenJet_pair_sube0_Mixing_RapAsym->Write();
	  hldsld_GenJet_pair_sube0_Mixing_RapAsym_1->Write();
	  hldsld_GenJet_pair_sube0_Mixing_RapAsym_2->Write();
	  hldsld_GenJet_pair_sube0_Mixing_RapAsym_3->Write();
	  hldsld_GenJet_pair_sube0_Mixing_RapAsym_sh->Write();
	  hldsld_GenJet_pair_sube0_Mixing_RapAsym_oh->Write();
	  hldsld_GenJet_pair_sube0_Mixing->Write();
	  // for sube = 1
	  hldsld_GenJet_pair_sube1_Mixing_RapAsym->Write();
	  hldsld_GenJet_pair_sube1_Mixing_RapAsym_1->Write();
	  hldsld_GenJet_pair_sube1_Mixing_RapAsym_2->Write();
	  hldsld_GenJet_pair_sube1_Mixing_RapAsym_3->Write();
	  hldsld_GenJet_pair_sube1_Mixing_RapAsym_sh->Write();
	  hldsld_GenJet_pair_sube1_Mixing_RapAsym_oh->Write();
	  hldsld_GenJet_pair_sube1_Mixing->Write();
	}
      
      // gen correlation histograms
      // leading jet - trk
      hldGenJet_Trk_Mixing_RapAsym->Write();
      hldGenJet_Trk_Mixing_RapAsym_1->Write();
      hldGenJet_Trk_Mixing_RapAsym_2->Write();
      hldGenJet_Trk_Mixing_RapAsym_3->Write();
      hldGenJet_Trk_Mixing_RapAsym_sh->Write();
      hldGenJet_Trk_Mixing_RapAsym_oh->Write();
      hldGenJet_Trk_Mixing->Write();

      if(!isRcJetGnTrk)
	{
	  // for sube = 0
	  hldGenJet_Trk_sube0_Mixing_RapAsym->Write();
	  hldGenJet_Trk_sube0_Mixing_RapAsym_1->Write();
	  hldGenJet_Trk_sube0_Mixing_RapAsym_2->Write();
	  hldGenJet_Trk_sube0_Mixing_RapAsym_3->Write();
	  hldGenJet_Trk_sube0_Mixing_RapAsym_sh->Write();
	  hldGenJet_Trk_sube0_Mixing_RapAsym_oh->Write();
	  hldGenJet_Trk_sube0_Mixing->Write();
	  // for sube = 1
	  hldGenJet_Trk_sube1_Mixing_RapAsym->Write();
	  hldGenJet_Trk_sube1_Mixing_RapAsym_1->Write();
	  hldGenJet_Trk_sube1_Mixing_RapAsym_2->Write();
	  hldGenJet_Trk_sube1_Mixing_RapAsym_3->Write();
	  hldGenJet_Trk_sube1_Mixing_RapAsym_sh->Write();
	  hldGenJet_Trk_sube1_Mixing_RapAsym_oh->Write();
	  hldGenJet_Trk_sube1_Mixing->Write();
	}
      
      /*
      // subleading jet - trk
      hsldGenJet_Trk_Mixing_RapAsym->Write();
      hsldGenJet_Trk_Mixing_RapAsym_1->Write();
      hsldGenJet_Trk_Mixing_RapAsym_2->Write();
      hsldGenJet_Trk_Mixing_RapAsym_3->Write();
      hsldGenJet_Trk_Mixing_RapAsym_sh->Write();
      hsldGenJet_Trk_Mixing_RapAsym_oh->Write();
      hsldGenJet_Trk_Mixing->Write();
      // for sube = 0
      hsldGenJet_Trk_sube0_Mixing_RapAsym->Write();
      hsldGenJet_Trk_sube0_Mixing_RapAsym_1->Write();
      hsldGenJet_Trk_sube0_Mixing_RapAsym_2->Write();
      hsldGenJet_Trk_sube0_Mixing_RapAsym_3->Write();
      hsldGenJet_Trk_sube0_Mixing_RapAsym_sh->Write();
      hsldGenJet_Trk_sube0_Mixing_RapAsym_oh->Write();
      hsldGenJet_Trk_sube0_Mixing->Write();
      // for sube = 1
      hsldGenJet_Trk_sube1_Mixing_RapAsym->Write();
      hsldGenJet_Trk_sube1_Mixing_RapAsym_1->Write();
      hsldGenJet_Trk_sube1_Mixing_RapAsym_2->Write();
      hsldGenJet_Trk_sube1_Mixing_RapAsym_3->Write();
      hsldGenJet_Trk_sube1_Mixing_RapAsym_sh->Write();
      hsldGenJet_Trk_sube1_Mixing_RapAsym_oh->Write();
      hsldGenJet_Trk_sube1_Mixing->Write();
      */
    }
}
