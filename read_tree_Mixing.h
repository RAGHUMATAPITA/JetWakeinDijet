#include "call_libraries.h"  // call libraries from ROOT and C++

// declare variables

//~~~~~~~~~~~~~~~~~~~~~ event quantities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
float vertexz_MB; // event z vertex
int hiBin_MB; // event centrality (used if use_centrality = true in input_variables.h)
ULong64_t event_MB;
UInt_t run_MB;
//~~~~~~~~~~~~~~~~~~~~~track quantities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// reco jets
const int nmax_MB = 99999;

// reco tracks
int ntrk_MB;                 // number of track
float trkpt_MB[nmax_MB];       // track pT
float trketa_MB[nmax_MB];      // track eta
float trkphi_MB[nmax_MB];      // track phi

// gen tracks
std::vector<float> *gen_trkpt_MB = 0;  // gen particle pT
std::vector<float> *gen_trketa_MB = 0; // gen particle eta
std::vector<float> *gen_trkphi_MB = 0; // gen particle phi
std::vector<int> *gen_trksube_MB = 0;   // gen particle pid

void read_tree_Mixing(TChain *tree, bool is_MC, TString colliding_system)
{
  tree->SetBranchStatus("*", 0); // disable all branches - this is important while reading big files
  
  // enable branches of interest -> see definition of each variables above
  
  // event quantities
  
  tree->SetBranchStatus("vz", 1);
  tree->SetBranchAddress("vz", &vertexz_MB);
  
  tree->SetBranchStatus("evt", 1);
  tree->SetBranchAddress("evt", &event_MB);
  
  tree->SetBranchStatus("run", 1);
  tree->SetBranchAddress("run", &run_MB);

  if(colliding_system=="PbPb")
    {
      tree->SetBranchStatus("hiBin", 1);
      tree->SetBranchAddress("hiBin", &hiBin_MB);
    }
  
  // reco or data track quantities
  
  tree->SetBranchStatus("nTrk", 1);
  tree->SetBranchAddress("nTrk", &ntrk_MB);
  
  tree->SetBranchStatus("trkPt", 1);
  tree->SetBranchAddress("trkPt", &trkpt_MB);
  
  tree->SetBranchStatus("trkEta", 1);
  tree->SetBranchAddress("trkEta", &trketa_MB);
  
  tree->SetBranchStatus("trkPhi", 1);
  tree->SetBranchAddress("trkPhi", &trkphi_MB);
  
  if (is_MC)
    {
      // gen particle quantities
      tree->SetBranchStatus("pt", 1);
      tree->SetBranchAddress("pt", &gen_trkpt_MB);
      
      tree->SetBranchStatus("eta", 1);
      tree->SetBranchAddress("eta", &gen_trketa_MB);
      
      tree->SetBranchStatus("phi", 1);
      tree->SetBranchAddress("phi", &gen_trkphi_MB);

      tree->SetBranchStatus("sube", 1);
      tree->SetBranchAddress("sube", &gen_trksube_MB);
    }
}
