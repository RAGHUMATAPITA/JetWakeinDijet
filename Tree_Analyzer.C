#include "call_libraries.h"  // call libraries from ROOT and C++
#include "input_variables.h" // input variables 
#include "histogram_definition_new.h" // define histograms
#include "read_tree.h" // read the TChains
#include "read_tree_Mixing.h" // read Tree for MB sample for Mixing
#include "JetCorrector.h" // reader for JEC
#include "JetUncertainty.h" // reader for JEU
#include "function_defination.h" // function defined here
 
void Tree_Analyzer(TString input_file, int itxtoutFile, TString out_file, TString colliding_system, int isMC)
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Strat Initializing PbPb and pp used input variables~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(isMC == 1)
    {
      is_MC = true;
    }
  else 
    {
      is_MC = false;
    }

  // for Centrality and Vz weight
  TFile* hCentVzFile;
  TF1* fVzweight;
  TF1* fCentweight;
  
  // for pt weight 
  TFile* fptFile_pp;
  TFile* fptFile_PbPb[NCentbin];
  TF1* fptWeight_PbPb[NCentbin];
  TF1* fptWeight;

  // for JER correction from JER fit, eta dependent SF
  TFile* fJERFile_pp;
  TFile* fJERFile_PbPb[NCentbin];
  TF1* fJERWeight_PbPb[NCentbin];
  TF1* fJERWeight;

  // for trk efficiecny
  TFile* ftrk_eff;
  TFile* ftrk_fak; // only for PbPb
  
  TH2F* htrk_eff_pp;
  TH2F* htrk_fak_pp;
  TH2F* htrk_sec_pp;

  TH3F* htrk_eff_PbPb;
  TH3F* htrk_fak_PbPb;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~:Setting for pp ref data and MC:~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TString colliding_system_filename;
  
  if(colliding_system == "pp")
    {
      isMBEventsMixing = false; // always false for pp
      colliding_system_filename = "pp";
      use_cent = false;
      doEventFilter = false;
      
      if(is_MC)
	{
	  colliding_system_filename = "Pythia";
	  is_CentBin_and_VZ_Weight = true;
	  is_ptWeight = false;
	  is_JES_JER = false;
	  is_Gen_Reco_Correlation = false;
	  is_JER_Correction = false;
	  isRcJetGnTrk = false; 
	}

      //if(!is_MC)
      if(!is_MC || is_MC)
	{
          is_JEU = false;
          is_JEU_up = false;
          is_JEU_down = false;
	}

      if(is_MC)
        {
          is_JER = true;
          is_JER_nom = true;
          is_JER_up = false;
          is_JER_down = false;
        }

      /*
      if(!is_MC){is_JetTrigger = true;} // for data
      else {is_JetTrigger = false;}
      */
      
      is_JetTrigger = true;
      jet_trigger = pp_jet_trigger;

      if(is_MC && is_CentBin_and_VZ_Weight)
	{
	  CentBin_and_VZ_Fun = pp_CentBin_and_VZ_Fun;
	  hCentVzFile = TFile::Open(Form("ReweightFile/%s", pp_CentBin_and_VZ_Fun.Data()));
	  fVzweight = (TF1*)hCentVzFile->Get("hvzFun");
	}

      if(is_MC && is_ptWeight)
	{
	  ptWeightFun = pp_ptWeight;
	  fptFile_pp = TFile::Open(Form("ReweightFile/%s.root", ptWeightFun.Data()));
	  fptWeight = (TF1*)fptFile_pp->Get("hptFun");
	}

      if(is_MC && is_JER)
	{
	  JER_File = JER_File_pp;
          fJERFile_pp = TFile::Open(Form("ReweightFile/%s.root", JER_File.Data()));
          fJERWeight = (TF1*)fJERFile_pp->Get("JER_Fit");
	}

      jet_collection = pp_jet_collection;
      JEC_file = pp_JEC_file;

      //if(!is_MC)
      if(!is_MC || is_MC)
	{
	  JEC_file_data = pp_JEC_file_data;
	  if(is_JEU)	  
	    {
	      JEU_file_data = pp_JEU_file_data;
	    }
	}
      
      if(is_MC && is_JER) 
	{
	  JER_file_data = pp_JER_file_data;
	}

      // for trk eff
      trk_eff_file = trk_eff_file_pp;
      ftrk_eff = TFile::Open(Form("Eff_File/%s", trk_eff_file.Data()));
      htrk_eff_pp = (TH2F*)ftrk_eff->Get("rEff");
      htrk_fak_pp = (TH2F*)ftrk_eff->Get("rFak");
      htrk_sec_pp = (TH2F*)ftrk_eff->Get("rSec");
      
      event_filter_str.resize(0);
      event_filter_str.push_back("pBeamScrapingFilter");
      event_filter_str.push_back("pPAprimaryVertexFilter");
      event_filter_str.push_back("HBHENoiseFilterResultRun2Loose");
      
    }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~:Setting for PbPb data and MC:~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(colliding_system == "PbPb")
    {
      isMBEventsMixing = false; // true/ false
      colliding_system_filename = "PbPb";
      
      use_cent = true;
      
      if(is_MC) 
	{
	  colliding_system_filename = "PH";
	  
	  doEventFilter = true; // default true
	  is_CentBin_and_VZ_Weight = true;
	  is_ptWeight = false;
	  is_JES_JER = false;
	  is_Gen_Reco_Correlation = false;
	  is_JER_Correction = false;
	  isRcJetGnTrk = false; 
	}
      
      //if(!is_MC)
      if(!is_MC || is_MC)
	{
	  is_JEU = false;
	  is_JEU_up = false;
	  is_JEU_down = false;
	}

      if(is_MC)
	{
	  is_JER = true;
	  is_JER_nom = true;
	  is_JER_up = false;
	  is_JER_down = false;
	}

      /*
      if(!is_MC){is_JetTrigger = true;} // for data
      else {is_JetTrigger = false;}
      */
      
      is_JetTrigger = true;
      jet_trigger = PbPb_jet_trigger;

      if(is_MC && is_CentBin_and_VZ_Weight)
	{
	  CentBin_and_VZ_Fun = PbPb_CentBin_and_VZ_Fun;
	  hCentVzFile = TFile::Open(Form("ReweightFile/%s", PbPb_CentBin_and_VZ_Fun.Data()));
	  fVzweight = (TF1*)hCentVzFile->Get("hvzFun");
	  fCentweight = (TF1*)hCentVzFile->Get("hCentFun");
	}

      if(is_MC && is_ptWeight)
        {
          ptWeightFun = PbPb_ptWeight;
	  for(int ictt = 0; ictt < NCentbin-1; ictt++)
            {
              fptFile_PbPb[ictt] = TFile::Open(Form("ReweightFile/%s_%d.root", ptWeightFun.Data(), ictt));
	      fptWeight_PbPb[ictt] = (TF1*)fptFile_PbPb[ictt]->Get("hptFun");
            }
        }

      if(is_MC && is_JER)
        {
          JER_File = JER_File_PbPb;
	  for(int ictt = 0; ictt < NCentbin-1; ictt++)
            {
	      fJERFile_PbPb[ictt] = TFile::Open(Form("ReweightFile/%s%d.root", JER_File.Data(), ictt));
	      fJERWeight_PbPb[ictt] = (TF1*)fJERFile_PbPb[ictt]->Get("JER_Fit");
	    }
        }

      jet_collection = PbPb_jet_collection;
      JEC_file = PbPb_JEC_file;
      
      //if(!is_MC)
      if(!is_MC || is_MC)
	{
	  JEC_file_data = PbPb_JEC_file_data;
	  if(is_JEU)
            {
              JEU_file_data = PbPb_JEU_file_data;
            }
	}

      if(is_MC && is_JER) 
	{
	  JER_file_data = PbPb_JER_file_data;
	}

       // for trk eff
      trk_eff_file = trk_eff_file_PbPb;
      trk_fak_file = trk_fak_file_PbPb;
      ftrk_eff = TFile::Open(Form("Eff_File/%s", trk_eff_file.Data()));
      ftrk_fak = TFile::Open(Form("Eff_File/%s", trk_fak_file.Data()));
      htrk_eff_PbPb = (TH3F*)ftrk_eff->Get("Eff3D");
      htrk_fak_PbPb = (TH3F*)ftrk_fak->Get("Fak3D");

      event_filter_str.resize(0);
      event_filter_str.push_back("pprimaryVertexFilter");
      event_filter_str.push_back("HBHENoiseFilterResultRun2Loose");
      event_filter_str.push_back("collisionEventSelectionAOD");
      event_filter_str.push_back("phfCoincFilter2Th4");
      event_filter_str.push_back("pclusterCompatibilityFilter");

      /*
      if(is_MC && colliding_system == "PbPb")
	{	  
	  // AOD MC
	  event_filter_str.push_back("pprimaryVertexFilter");
	  event_filter_str.push_back("HBHENoiseFilterResultRun2Loose");
	  event_filter_str.push_back("collisionEventSelectionAOD");
	  event_filter_str.push_back("phfCoincFilter2Th4");
	  event_filter_str.push_back("pclusterCompatibilityFilter");
	}
      else if(!is_MC && colliding_system == "PbPb")
	{
	  // MAOD data
	  event_filter_str.push_back("pclusterCompatibilityFilter");
	  event_filter_str.push_back("pprimaryVertexFilter");
	  event_filter_str.push_back("pphfCoincFilter2Th4");
	}
      */
    }

  TFile *xjetFrac = new TFile();
  TH1D* hXjetFrac_vs_genpT = new TH1D();
  TF1* fXjetFrac_vs_genpT = new TF1();

  if(is_MC && doEventFilter && colliding_system == "PbPb")
    {
      xjetFrac = TFile::Open(Form("ReweightFile/%s", xJetsw.Data()));
      hXjetFrac_vs_genpT = (TH1D*)xjetFrac->Get("hUM_ld_CorrpT_gt_noW"); // from histogram
      //fXjetFrac_vs_genpT = (TF1*)xjetFrac->Get("fgaus"); // from fitting
    }

  
  //~~~~~~~~~~~~~~~~~~~~~End Initializing PbPb and pp used input variables~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Start printing the input informations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  std::cout<<endl;
  std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Start printing the input informations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  std::cout<<endl;
  std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~Event quantities~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  std::cout<<endl;
  
  std::cout<<"colliding_system is: "<<colliding_system.Data()<<std::endl;
  std::cout<<"is it Monte Carlo: "<<std::boolalpha<<is_MC<<std::endl;
  std::cout<<"use Centrality: "<<std::boolalpha<<use_cent<<std::endl;
  std::cout<<"doEventFilter: "<<std::boolalpha<<doEventFilter<<std::endl;
  std::cout<<"skimed event_filter_str size is: "<<event_filter_str.size()<<std::endl;
  std::cout<<"skimed event filters are: "<<std::endl;
  
  for(unsigned int ifl = 0; ifl < event_filter_str.size(); ifl++)
    {
      std::cout<<event_filter_str[ifl].Data()<<", ";
    }
  std::cout<<std::endl;
  
  std::cout<<"vertex z cut: "<<vz_cut_min<<" to "<<vz_cut_max<<" cm"<<std::endl;
  
  if(colliding_system == "PbPb")
    {
      std::cout<<"Maximum hiBin cut: "<<centCut<<std::endl;
    }
  
  if(is_MC)
    {
      std::cout<<"Minimum pThat cut: "<<pthat_cut<<" GeV"<<std::endl;   
    }

  std::cout<<"is_CentBin_and_VZ_Weight: "<<std::boolalpha<<is_CentBin_and_VZ_Weight<<std::endl;
  if(is_MC && is_CentBin_and_VZ_Weight)
    {
      std::cout<<"Cent bin and Vertex z weight file is: "<<hCentVzFile->GetName()<<std::endl;
    }
  
  std::cout<<endl;

  std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~Jet quantities~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  std::cout<<endl;

  std::cout<<"Jet Collection: "<<jet_collection.Data()<<std::endl;
  
  std::cout<<"is_JetTrigger: "<<std::boolalpha<<is_JetTrigger<<std::endl;
  if(is_JetTrigger)
    {
      std::cout<<"Applied jet trigger: "<<jet_trigger.Data()<<std::endl;
    }

  std::cout<<"JEC File used: "<<JEC_file.Data()<<std::endl;
  if(!is_MC)
  {
    std::cout<<"JEC File for data used: "<<JEC_file_data.Data()<<std::endl;
  }

  std::cout<<"is_JEU: "<<std::boolalpha<<is_JEU<<std::endl;
  if((!is_MC || is_MC) && is_JEU) 
    {
      std::cout<<"JEU File for data used: "<<JEU_file_data.c_str()<<std::endl; 
    }

  std::cout<<"is_JES_JER: "<<std::boolalpha<<is_JES_JER<<std::endl;

  if(is_MC && doEventFilter && colliding_system == "PbPb")
    {
      std::cout<<"xJetsw file is: "<<xjetFrac->GetName()<<std::endl;
    }
  
  std::cout<<"jet Eta cut: "<<jet_eta_min_cut<<" to "<<jet_eta_max_cut<<std::endl;
  std::cout<<"jet pT cut: "<<jet_pt_min_cut<<" to "<<jet_pt_max_cut<<" GeV"<<std::endl;
  std::cout<<"jet leading pT cut: "<<leading_pT_min_cut<<" to "<<leading_pT_max_cut<<" GeV"<<std::endl;
  std::cout<<"jet subleading pT cut: "<<subleading_pT_min_cut<<" to "<<subleading_pT_max_cut<<" GeV"<<std::endl;
  std::cout<<"is_ptWeight : "<<std::boolalpha<<is_ptWeight<<std::endl;
  
  if(is_MC && is_ptWeight)
    {
      if(colliding_system == "pp")
	{
	  std::cout<<"pt weight file is : "<<fptFile_pp->GetName()<<std::endl;
	}
      else if(colliding_system == "PbPb")
	{
	  for(int ictt = 0; ictt < NCentbin-1; ictt++)
	    {
	      std::cout<<"pt weight file is : "<<fptFile_PbPb[ictt]->GetName()<<std::endl;
	    }
	}
    }
  
  std::cout<<"is_JER_Correction : "<<std::boolalpha<<is_JER_Correction<<std::endl;
  std::cout<<"is_JER, is_JER_nom, is_JER_up, and is_JER_down: "<<std::boolalpha<<is_JER<<"  "<<is_JER_nom<<"  "<<is_JER_up<<"  "<<is_JER_down<<std::endl;

  if(is_MC && is_JER)
    {
      if(colliding_system == "pp")
	{
	  std::cout<<"JER weight file is : "<<fJERFile_pp->GetName()<<std::endl;
	}
      
      if(colliding_system == "PbPb")
	{
	  for(int ictt = 0; ictt < NCentbin-1; ictt++)
	    {
	      std::cout<<"JER weight file is : "<<fJERFile_PbPb[ictt]->GetName()<<std::endl;
	    }
	} 
    }

  std::cout<<"is_Gen_Reco_Correlation : "<<std::boolalpha<<is_Gen_Reco_Correlation<<std::endl;
  //std::cout<<"is_EWTA_Correction : "<<std::boolalpha<<is_EWTA_Correction<<std::endl;

  std::cout<<endl;
  std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~Track quantities~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  std::cout<<endl;
  std::cout<<"For reco/data, Only HP tracks are considered"<<std::endl;
  std::cout<<"Only charged tracks are considered"<<std::endl;
  std::cout<<"trk min. & max. pT cut: "<<trk_pt_min_cut<<" GeV,  "<<trk_pt_max_cut<<" GeV"<<std::endl;
  std::cout<<"trk eta cut: "<<trk_eta_cut<<std::endl;
  std::cout<<"trk DCA XY and Z sig. cut: "<<trk_dca_xy_cut<<",  "<<trk_dca_z_cut<<std::endl;
  std::cout<<"trk pT Reso. cut: "<<trk_pt_resolution_cut<<std::endl;
  if(colliding_system == "PbPb")
    {
      std::cout<<"trk nhits cut: "<<nhits<<std::endl;
      std::cout<<"trk chi2/ndf/nlayes cut: "<<chi2_ndf_nlayer_cut<<std::endl;
      std::cout<<"if(trk_algo == 6 && trk_mva < 0.98) continue; cut applied"<<std::endl;
    }
  std::cout<<"Efficiecny correction file: "<<ftrk_eff->GetName()<<std::endl;
  if(colliding_system == "PbPb") std::cout<<"Fake correction file: "<<ftrk_fak->GetName()<<std::endl;
  
  std::cout<<endl;
  std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~Mixing quantities~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  std::cout<<endl;
  std::cout<<"Number of Mixing events: "<<bkgFactor<<std::endl;
  std::cout<<"Z vertex between two events: "<<(double)Dvz_cut<<" cm"<<std::endl;
  if(colliding_system == "PbPb") std::cout<<"Centrality difference between two events: "<<(Dhibin_cut)/2.<<"%"<<std::endl;
  else if(colliding_system == "pp") std::cout<<"Ntrack difference between two events: "<<(int)Dntrk_cut<<std::endl;
  std::cout<<endl;

  std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~Please note that~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  std::cout<<endl;
  if(isMBEventsMixing){std::cout<<"isMBEventsMixing: "<<std::boolalpha<<isMBEventsMixing<<std::endl; std::cout<<"Mixing will be done using MB events"<<std::endl;}
  else {std::cout<<"isMBEventsMixing: "<<std::boolalpha<<isMBEventsMixing<<std::endl; std::cout<<"Mixing will be done using dijet events"<<std::endl;}
  
  if(isRcJetGnTrk){std::cout<<"isRcJetGnTrk: "<<std::boolalpha<<isRcJetGnTrk<<std::endl; std::cout<<"Correlation will be done between reco jets and gen tracks and vice varsa"<<std::endl;}
  else {std::cout<<"isRcJetGnTrk: "<<std::boolalpha<<isRcJetGnTrk<<std::endl; std::cout<<"Correlation will be done between reco jets and reco tracks and vice varsa"<<std::endl;}

  std::cout<<endl;
  std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~End printing the input informations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~End printing the input informations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  std::cout<<endl;

  //~~~~~~~~~~~~~~calling sumw2 for all the histograms~~~~~~~~~~~~~~~~~~~
  sumw2(); 
  TH1::SetDefaultSumw2();
  
  //~~~~~~~~~~~~~~input JEC file~~~~~~~~~~~~~~~~~~~~~~~~~
  vector<string> JECFiles;
  JECFiles.push_back(Form("JEC_files/%s", JEC_file.Data()));
  if(!is_MC)
    {
      JECFiles.push_back(Form("JEC_files/%s", JEC_file_data.Data())); // for data only
      std::cout<<"It is runing for data (L2L3 residual file is attached)"<<std::endl;
      std::cout<<endl;
    }

  JetCorrector JEC(JECFiles);
  
  // for uncertainity in JEC in data (Used in systematics)
  /*
  JetUncertainty* JEU_ld = nullptr;
  JetUncertainty* JEU_sld = nullptr;
  */
  JetUncertainty* JEU = nullptr;
  
  //if(!is_MC && is_JEU)
  if((!is_MC || is_MC)  && is_JEU)
    {
      /*
      JEU_ld = new JetUncertainty(JEU_file_data);
      JEU_sld = new JetUncertainty(JEU_file_data);
      */
      JEU = new JetUncertainty(JEU_file_data); 
    }

  // for uncertainity in JER
  if(is_MC && is_JER)
    {
      if(is_JER_nom){SFCondition = 0;}
      else if(is_JER_up){SFCondition = 1;}
      else if(is_JER_down){SFCondition = 2;}
    }

  //  open input forest/skim file
  fstream openInputFile;
  openInputFile.open(Form("%s",input_file.Data()), ios::in);
  if(!openInputFile.is_open())
    {
      cout << "List of input files not founded!" << endl;
      return;
    }

  // Make a chain and a vector of file names
  std::vector<TString> file_name_vector;
  string file_chain;
  while(getline(openInputFile, file_chain))
    {
      //file_name_vector.push_back(Form("root://cmsxrootd.fnal.gov/%s", file_chain.c_str()));
      if(colliding_system == "PbPb")
	{
	  //file_name_vector.push_back(Form("root://cmsxrootd.fnal.gov:1094/%s", file_chain.c_str()));
	  //file_name_vector.push_back(Form("davs://xrootd-vanderbilt.sites.opensciencegrid.org:1094/%s", file_chain.c_str()));
	  //file_name_vector.push_back(Form("root://xrootd-vanderbilt.sites.opensciencegrid.org/%s", file_chain.c_str()));
	  file_name_vector.push_back(file_chain.c_str());
	}
      else if(colliding_system == "pp")
	{
	  if(is_MC)
	    {
	      file_name_vector.push_back(file_chain.c_str());
	    }
	  else file_name_vector.push_back(Form("root://cmsxrootd.fnal.gov:1094/%s", file_chain.c_str()));
	}
      //file_name_vector.push_back(file_chain.c_str());
    }
  openInputFile.close();

  // Read the trees to be added in the Chain
  TChain *hlt_tree = new TChain("hltanalysis/HltTree");
  TChain *hea_tree = new TChain("hiEvtAnalyzer/HiTree");
  TChain *ski_tree = new TChain("skimanalysis/HltTree");
  TChain *jet_tree = new TChain(Form("%s/t",jet_collection.Data()));
  TChain *trk_tree = new TChain("ppTrack/trackTree");
  TChain *gen_tree = nullptr;

  if(is_MC)
    {
      gen_tree = new TChain("HiGenParticleAna/hi");
    }
  
  for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++)
    {
      /*
      // Convert TString to std::string
      std::string filePath = std::string(*listIterator);
      // Check if file exists on disk using `gfal-ls`
      std::string command = "gfal-ls " + filePath + " > /dev/null 2>&1";
      bool isFileOnDisk = (gSystem->Exec(command.c_str()) == 0);
      if (!isFileOnDisk) {
        std::cout << "File not found on disk: " << filePath << std::endl;
        continue;
      }
      */

      std::string filePath = std::string(*listIterator);

      TFile *testfile = TFile::Open(*listIterator);

      if(!testfile || testfile->IsZombie() || testfile->TestBit(TFile::kRecovered))
	{
	  cout << "File: " << *listIterator << " failed to open" << endl;
	  delete testfile;  // Cleanup memory to avoid leaks
	  testfile = nullptr;
	  continue;
	}
      else
	{
	  if(gSystem->AccessPathName(filePath.c_str()) == 0)
	    {
	      cout << "Adding file:--- " << *listIterator << "--- to the chains" << endl;

	      hlt_tree->Add(*listIterator);
	      hea_tree->Add(*listIterator);
	      ski_tree->Add(*listIterator);
	      jet_tree->Add(*listIterator);
	      trk_tree->Add(*listIterator);
	      if(is_MC)
		{
		  gen_tree->Add(*listIterator);
		}
	    }
	  delete testfile; // Close file after use
	  testfile = nullptr;
	}
    }

  hlt_tree->AddFriend(hea_tree);
  hlt_tree->AddFriend(ski_tree);	
  hlt_tree->AddFriend(jet_tree);
  hlt_tree->AddFriend(trk_tree);
  if(is_MC)
    {
      hlt_tree->AddFriend(gen_tree);
    }

  // calling function to read forest/skim tree
  read_tree(hlt_tree, is_MC, jet_trigger, colliding_system, event_filter_str); // access the tree informations
  
  int nevents = hlt_tree->GetEntries(); // number of events

  std::cout<<std::endl;
  cout << "Total number of events in those files: "<< nevents << endl;
  std::cout<<std::endl;
  
  //~~~~~~~~~~~~Define vectors use for jet tracks Signal and Mixed event correlation~~~~~~~~~~~~~~~~~~~~~~
  // 1D event vectors
  //for gen
  std::vector<double> Gen_Evtw_vec_1D;                  // event weight
  std::vector<double> Gen_Vertexz_vec_1D;                   // vertex z
  std::vector<Long64_t> Gen_Evtno_vec_1D;                   // vertex z
  std::vector<int> Gen_EvtCount_vec_1D;                     // event count no 
  std::vector<int> Gen_HiBin_vec_1D;                        // centrality
  std::vector<int> Gen_HiBinValue_vec_1D;                   // centrality
    
  //for reco or data
  std::vector<double> Reco_Evtw_vec_1D;                 // event weight
  std::vector<double> Vertexz_vec_1D;                   // vertex z
  std::vector<Long64_t> Reco_Evtno_vec_1D;                   // vertex z
  std::vector<int> Reco_EvtCount_vec_1D;                     // event count no 
  std::vector<int> HiBin_vec_1D;                        // centrality
  std::vector<int> HiBinValue_vec_1D;                   // centrality
  
  // 1D leading jet vectors (Event wise)
  //for gen
  std::vector<TVector3> Gen_Filtered_ldJet_CorrpT_vec_1D;      // leading corrected jet
  std::vector<int> Gen_Filtered_ldrefpartonB_vec_1D;           // leading jet flavor
  std::vector<double> Gen_Filtered_ldJetW_vec_1D;              // leading jet weight
  //for reco or data
  std::vector<TVector3> Reco_Filtered_ldJet_CorrpT_vec_1D;      // leading corrected jet
  std::vector<int> Reco_Filtered_ldrefpartonB_vec_1D;           // leading jet flavor
  std::vector<double> Reco_Filtered_ldJetW_vec_1D;              // leading jet weight

  //1D subleading jet vectors (Event wise)
  //for gen
  std::vector<TVector3> Gen_Filtered_sldJet_CorrpT_vec_1D;     // sub leading corrected jet
  std::vector<int> Gen_Filtered_sldrefpartonB_vec_1D;          // sub leading jet flavor
  std::vector<double> Gen_Filtered_sldJetW_vec_1D;             // sub leading jet weight
  //for reco or data
  std::vector<TVector3> Reco_Filtered_sldJet_CorrpT_vec_1D;     // sub leading corrected jet
  std::vector<int> Reco_Filtered_sldrefpartonB_vec_1D;          // sub leading jet flavor
  std::vector<double> Reco_Filtered_sldJetW_vec_1D;              // sub leading jet weight

  //2D tracks vectors (Event wise and track wise)
  //gen vector
  std::vector<std::vector<TVector3>> Gen_Filtered_InclTrk_pT_vec_2D;            // trk vec
  std::vector<std::vector<double>> Gen_Filtered_InclTrkW_vec_2D;                // trk weight 
  std::vector<std::vector<int>> Gen_Filtered_InclTrkCharge_vec_2D;              // trk charge
  std::vector<std::vector<int>> Gen_Filtered_InclTrkSube_vec_2D;                // trk sube
  //gen jets and reco tracks
  std::vector<std::vector<TVector3>> GenReco_Filtered_InclTrk_pT_vec_2D;            // trk vec
  std::vector<std::vector<double>> GenReco_Filtered_InclTrkW_vec_2D;                // trk weight 
  std::vector<std::vector<int>> GenReco_Filtered_InclTrkCharge_vec_2D;              // trk charge
  std::vector<std::vector<int>> GenReco_Filtered_InclTrkSube_vec_2D;                // trk sube
  
  //reco/data vector
  std::vector<std::vector<TVector3>> Reco_Filtered_InclTrk_pT_vec_2D;           // trk vec
  std::vector<std::vector<double>> Reco_Filtered_InclTrkW_vec_2D;               // trk weight
  std::vector<std::vector<int>> Reco_Filtered_InclTrkCharge_vec_2D;             // trk charge
  std::vector<std::vector<int>> Reco_Filtered_InclTrkSube_vec_2D;                // trk sube
  //reco jets and gen tracks
  std::vector<std::vector<TVector3>> RecoGen_Filtered_InclTrk_pT_vec_2D;           // trk vec
  std::vector<std::vector<double>> RecoGen_Filtered_InclTrkW_vec_2D;               // trk weight
  std::vector<std::vector<int>> RecoGen_Filtered_InclTrkCharge_vec_2D;             // trk charge
  std::vector<std::vector<int>> RecoGen_Filtered_InclTrkSube_vec_2D;                // trk sube
  //~~~~~~~~~~~~End Defining vectors use for jet tracks Signal and Mixed event correlation~~~~~~~~~~~~~~~~~~~~~~
  
  int count = 0, count_wcut = 0; // to counts the jets

  int evtcount_Reco = 0; // required for mixing
  int evtcount_Gen = 0; // required for mixing

  for(int i = 0; i < nevents; i++) //event loop //start
  //for(int i = 0; i < 1000; i++) //event loop start
    {
      hlt_tree->GetEntry(i);
      
      if(i%10000 == 0)
	{
	  std::cout<<i<<"  events running"<<std::endl;
	}
      
      count = count+nref; // count total number of jets in the events without any cuts
      
      if(vertexz <= vz_cut_min || vertexz >= vz_cut_max) continue; // apply z vertex cut

      hEvents->AddBinContent(1,1); // after VZ cut

      // determine centrality here
      int centbin;
      int ctbin = -1;
      
      //int hiCentBin = getHiBinFromhiHF(hiHF, 0); // 0 = nominal, 1 = up, 2 = down // for systematics

      if(use_cent) // if you use centrality
	{
	  if(is_MC)
            {
	      if(hiBin <= 9) continue; 
	      centbin = hiBin - 10; // match MC multiplicity with data multiplicity
	    }
          else
            {
              centbin = hiBin;
	      //centbin = hiCentBin; //for systematics 
            }

	  if(centbin >= centCut || centbin < 0) continue;
	  //if(centbin >= centCut || centbin < 100) continue;
	  ctbin = hCentBin->FindBin(centbin) - 1;
	}
      else // if you use multiplicity 
	{
	  centbin = -1;
	  if(centbin >= mult_Cut) continue;
	  ctbin = 1;
	}

      hEvents->AddBinContent(2,1); // after cent cut

      // for pT weight to match MC pT with the data pT
      if(is_MC && is_ptWeight && colliding_system == "PbPb")
        {
	  if(ctbin == 0){fptWeight = fptWeight_PbPb[0];}
	  else if(ctbin == 1){fptWeight = fptWeight_PbPb[1];}
	  else if(ctbin == 2){fptWeight = fptWeight_PbPb[2];}
	  else if(ctbin == 3){fptWeight = fptWeight_PbPb[3];}
	  else if(ctbin == 4){fptWeight = fptWeight_PbPb[3];}
	  else{std::cout<<"Centrality bins are more than 4"<<std::endl; break;}
	}

      // for JER correction (To smear dist. according to JER)
      if(is_MC && is_JER && colliding_system == "PbPb")
        {
	  if(ctbin == 0){fJERWeight = fJERWeight_PbPb[0];}
	  else if(ctbin == 1){fJERWeight = fJERWeight_PbPb[1];}
	  else if(ctbin == 2){fJERWeight = fJERWeight_PbPb[2];}
	  else if(ctbin == 3){fJERWeight = fJERWeight_PbPb[3];}
	  else if(ctbin == 4){fJERWeight = fJERWeight_PbPb[4];}
	  else{std::cout<<"Centrality bins are more than 4"<<std::endl; break;}
	}

      double ptHatw = 1.;
      double pTHat = 0.;
      if(is_MC)
	{
	  ptHatw = weight;
	  pTHat = pthat;
	  if(pTHat == 0 || pTHat <= pthat_cut) continue; // apply pTHat cut
	}
	 
      hEvents->AddBinContent(3,1); // after pT hat cut
      
      bool skimmed_evtfilter = false;

      for(int ii = 0; ii < (int) event_filter_str.size(); ii++)
	{
	  if (event_filter_bool[ii] != 1) // condition for the skimmed event filters
	    {
	      skimmed_evtfilter = true;
	      break;
	    }
	}

      if(skimmed_evtfilter) continue; // apply the skimmed event filters

      hEvents->AddBinContent(4,1); // after event cuts

      if(is_JetTrigger)
	{
	  if(jet_trigger_bit != 1) continue; // apply jet trigger
	}

      hEvents->AddBinContent(5,1); // after jet trigger cut

      if(nref <= 0) continue; // if there is no jets in an event
      if(is_MC && ngen <= 0) continue; // if there is no jets in an event
			       
      hEvents->AddBinContent(6,1); // after nref and ngen <= 1 cut

      // determine reco/data event weight (pTHatw, centrality, and VZ weight) here
      double Evtw = ptHatw; // Event weight
            
      if(is_MC && is_CentBin_and_VZ_Weight)
	{
	  if(colliding_system == "pp")
	    {
	      Evtw = ptHatw*(fVzweight->Eval(vertexz));
	    }
	  else if(colliding_system == "PbPb")
	    {
	      Evtw = ptHatw*(fVzweight->Eval(vertexz))*(fCentweight->Eval(centbin));
	    }
	}

      // prepare for reco event filter based on unmatched leading pt > ptHat
      double mum_ld_CorrpT = -999., mum_ld_CorrEta = -999., mum_ld_CorrPhi = -999.;// leading matched unmatched jet quantities   
      int mum_ld_index = -999, mum_ld_falvorB = -9999; double mum_ld_refpt = -9999.;

      std::vector<double> highest_GenpT;
      double high_GenpT = -999.;

      bool excludeEvents = false;

      double genEvtw = Evtw;

      if(is_MC && doEventFilter && colliding_system == "PbPb")
        {
          for(int iref = 0; iref < (int)nref; iref++) // reco jet loop to calculate leading jet for the filter
            {
              if(trackMax[iref]/rawpt[iref] < 0.01)continue; // Cut for jets for very low maxium pT track                  
              if(trackMax[iref]/rawpt[iref] > 0.98)continue; // Cut for jets where all the pT is taken by one track        
	      
              JEC.SetJetPT(rawpt[iref]);
              JEC.SetJetEta(jteta[iref]);
              JEC.SetJetPhi(jtphi[iref]);
              float jet_ptcorr = JEC.GetCorrectedPT();

	      int refpartonnn = -99;
              if(fabs(refparton_flavor[iref]) >= 1 && fabs(refparton_flavor[iref]) <= 6)
                {
                  refpartonnn = fabs(refparton_flavor[iref]);
                }
              else if(fabs(refparton_flavor[iref]) == 21)
                {
                  refpartonnn = 7;
                }
	      else 
		{
		  refpartonnn = 0;
		}
	      
	      int refpartonnnB = -99;
              if(fabs(refparton_flavorForB[iref]) >= 1 && fabs(refparton_flavorForB[iref]) <= 6)
                {
                  refpartonnnB = fabs(refparton_flavorForB[iref]);
                }
              else if(fabs(refparton_flavorForB[iref]) == 21)
                {
                  refpartonnnB = 7;
                }
	      else
		{
		  refpartonnnB = 0;
		}

	      if(refpartonnn == -99 || refpartonnnB == -99) continue;

	      find_leading_jets(jet_ptcorr, jteta[iref], jtphi[iref], iref, refpt[iref], refpartonnnB, mum_ld_CorrpT, mum_ld_CorrEta,mum_ld_CorrPhi, mum_ld_index, mum_ld_refpt, mum_ld_falvorB);

	    } // reco jet loop end
			
	  for (int igen = 0; igen < (int)ngen; igen++) // gen jet loop to calculate leading jet for the filter
	    {
	      highest_GenpT.push_back(gen_jtpt[igen]);
	    }

	  std::sort(highest_GenpT.begin(), highest_GenpT.end()); // sort the vector
	  
	  if(highest_GenpT.size() >= 1)
	    {
	      high_GenpT = highest_GenpT[highest_GenpT.size()-1];
	    }

	  highest_GenpT.clear(); // clear the gen pT vector
	  	  
          if(mum_ld_index >= 0)
            {
              if((TMath::Abs(mum_ld_CorrEta) < jet_eta_max_cut ) && (mum_ld_CorrpT > leading_pT_min_cut && mum_ld_CorrpT < leading_pT_max_cut))
		{
                  if((int)refpt[mum_ld_index] == -999 && mum_ld_CorrpT > pTHat)
                    {
                      excludeEvents = true;
                    }
                  else excludeEvents = false;
                }
            }
          if(!excludeEvents && (int)high_GenpT != -999)
            {
              double genw = hXjetFrac_vs_genpT->GetBinContent(hXjetFrac_vs_genpT->FindBin(high_GenpT));    
              //double genw = fXjetFrac_vs_genpT->Eval(high_GenpT);
              genEvtw = (Evtw)*(1./(1. - genw));
            }
        } // is_MC, doEventFilter , and PbPb condition

      //std::cout<<"after excludeEvents: "<<(int)excludeEvents<<"  "<<genEvtw<<"  "<<mum_ld_CorrpT<<"  "<<pTHat<<std::endl;
            
      if(is_MC && doEventFilter && colliding_system == "PbPb")
	{
	  if(excludeEvents) continue; // filter the events based on leading unmacthed Jet pT > pThat
	}
      
      hEvents->AddBinContent(7,1); // after excludeEvents

      // Fill the bascis event quantities histograms
      if(is_MC)
	{
	  hpthat->Fill(pTHat, Evtw);
	  hpthatW->Fill(ptHatw);
	}

      hCent->Fill(centbin, Evtw);
      hZvtx->Fill(vertexz, Evtw);

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:Reco/Data jet loop start:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      //Jet loop start
      
      double ld_RawpT=-999.;   int ld_RawIndex=-999;      // leading Raw jet quantities                          
      double sld_RawpT=-999.;  int sld_RawIndex=-999;     // subleading Raw jet quantities                       
      double ld_CorrpT=-999.;  int ld_CorrIndex=-999;     // leading corrected jet quantities                       
      double sld_CorrpT=-999.; int sld_CorrIndex=-999;    // subleading corrected jet quantities
      double ld_Matched_CorrpT=-999.;  int ld_Matched_CorrIndex=-999;     // leading corrected jet quantities                       
      double sld_Matched_CorrpT=-999.; int sld_Matched_CorrIndex=-999;    // subleading corrected jet quantities    

      double ld_refpT=-999., sld_refpT=-999., ld_Matched_refpT=-999., sld_Matched_refpT=-999.;

      for (int j = 0; j < nref; j++) //Jet loop start
	{
	  if(trackMax[j]/rawpt[j] < 0.01)continue; // Cut for jets for very low maxium pT track
	  if(trackMax[j]/rawpt[j] > 0.98)continue; // Cut for jets where all the pT is taken by one track
	  
	  double jet_pt_raw = rawpt[j];
	  double jet_eta = jteta[j];
	  double jet_phi = jtphi[j];
	  double ref_jet_pt = refpt[j];

	  // JEC correction
	  JEC.SetJetPT(jet_pt_raw); 
	  JEC.SetJetEta(jet_eta); 
	  JEC.SetJetPhi(jet_phi);

	  float jet_pt_corr = JEC.GetCorrectedPT();
	  float jet_pt_corr_before = JEC.GetCorrectedPT();

	  // JEU correction
	  //if(!is_MC && is_JEU)
	  if((!is_MC || is_MC)  && is_JEU)
	    {
	      JEU->SetJetPT(jet_pt_corr);
	      JEU->SetJetEta(jet_eta);
	      JEU->SetJetPhi(jet_phi);
	      
	      if(is_JEU_down && !is_JEU_up)
		{
		  jet_pt_corr = jet_pt_corr*(1. - (JEU->GetUncertainty().first));
		}
	      else if(is_JEU_up && !is_JEU_down)
		{
		  jet_pt_corr = jet_pt_corr*(1. + (JEU->GetUncertainty().second));
		}
	    }
	  
	  // JER eta dependent scale factor
	  if(is_MC && is_JER)
	    {
	      double resolution_factor = 1.;
	      if(jet_eta <= jet_eta_min_cut || jet_eta >= jet_eta_max_cut)
		{
		  resolution_factor = 1.;
		}
	      else
		{
		  resolution_factor = getJERSFFromEta(jet_eta, SFCondition, colliding_system);
		}

	      double extraResolution = TMath::Sqrt(TMath::Max(resolution_factor*resolution_factor - 1.0, 0.0)); // from JetMET
	      double sigma_smear  = extraResolution*fJERWeight->Eval(jet_pt_corr); // some % worst --> from JetMET times JER
              if(jet_pt_corr < 0.) sigma_smear = extraResolution*fJERWeight->Eval(0.1);
              if(jet_pt_corr > 4999.) sigma_smear = extraResolution*fJERWeight->Eval(4998.9);

	      /*
	      gRandom->SetSeed(0);
	      double JER_smear = gRandom->Gaus(1,sigma_smear);
	      while( JER_smear < 0 ){ JER_smear = gRandom->Gaus(1,sigma_smear); }
	      */
	      
	      TRandom3 random(12345);
	      double JER_smear = random.Gaus(1,sigma_smear);
	      while( JER_smear < 0 ){ JER_smear = random.Gaus(1,sigma_smear); }
	      	      
	      jet_pt_corr = jet_pt_corr*JER_smear;
	      //std::cout<<"Jet pt after JERSF: "<<jet_pt_corr<<"  "<<ref_jet_pt<<"  "<<fJERWeight->Eval(ref_jet_pt)<<"  "<<jet_eta<<"  "<<resolution_factor<<"  "<<std::endl;
	    }
	  //std::cout<<"Jet pt before and after JERSF: "<<jet_pt_corr_before<<"  "<<jet_pt_corr<<std::endl;
	  
	  // jet flavour for in MC
	  int refpartonB = 0, refpartonB_PM = 0;
	  if(is_MC)
            {
	      if(fabs(refparton_flavorForB[j]) >= 1 && fabs(refparton_flavorForB[j]) <= 6)
		{
		  refpartonB = fabs(refparton_flavorForB[j]);
		  
		  if(refparton_flavorForB[j] > 0)
		    {
		      refpartonB_PM = fabs(refparton_flavorForB[j]);
		    }
		  else if(refparton_flavorForB[j] < 0)
		    {
		      refpartonB_PM = fabs(refparton_flavorForB[j])+8;
		    }
		}
	      else if(fabs(refparton_flavorForB[j]) == 21)
		{
		  refpartonB = 7;
		  refpartonB_PM = 7;
		}
	      else
		{
		  refpartonB = 0; 
		  refpartonB_PM = 0;
		}
	    }
	  else // for data
	    {
	      refpartonB = 1;
	      refpartonB_PM = 1;
	    }

	  // find leading subleading reco/data jets
	  //find_leading_subleading_jets(jet_pt_raw, j, ld_RawpT, ld_RawIndex , sld_RawpT, sld_RawIndex);       // for Raw pt 
	  //find_leading_subleading_jets(jet_pt_corr, j, ld_CorrpT, ld_CorrIndex, sld_CorrpT, sld_CorrIndex);   // for corrected pt
	  find_leading_subleading_jets(jet_pt_corr, ref_jet_pt, j, ld_CorrpT, ld_refpT, ld_CorrIndex, sld_CorrpT, sld_refpT, sld_CorrIndex);   // for corrected pt

	  // find leading subleading reco matched jets
	  if(is_MC)
	    {
	      if(ref_jet_pt > 0. && jet_pt_raw > 0. && jet_pt_corr > 0.)  // discrad the unmacthed jet (refpt == -999.)                       
		{
		  //find_leading_subleading_jets(jet_pt_corr, j, ld_Matched_CorrpT, ld_Matched_CorrIndex, sld_Matched_CorrpT, sld_Matched_CorrIndex);
		  find_leading_subleading_jets(jet_pt_corr, ref_jet_pt, j, ld_Matched_CorrpT, ld_Matched_refpT, ld_Matched_CorrIndex, sld_Matched_CorrpT, sld_Matched_refpT, sld_Matched_CorrIndex);
		}
	    }
	  
	  if(is_MC && is_JES_JER)
	    {
	      if(jet_eta > jet_eta_min_cut && jet_eta < jet_eta_max_cut)
		{
		  double JES_ratio_corrpT_vs_refpT = jet_pt_corr/ref_jet_pt;
		  double Jet_CorrpT_refpT_ctbin_flavour[4] = {JES_ratio_corrpT_vs_refpT, ref_jet_pt, (double)ctbin, (double)refpartonB};
		  //hJer_Jes_CorrpT_refpT_ctbin_flavour_W->Fill(Jet_CorrpT_refpT_ctbin_flavour, (Evtw));
		  hJer_Jes_CorrpT_refpT_ctbin_flavour_W->Fill(Jet_CorrpT_refpT_ctbin_flavour, ptHatw);
		}
	    }
	  
	  if(jet_eta > jet_eta_min_cut && jet_eta < jet_eta_max_cut)
	    {
	      if(jet_pt_corr > jet_pt_min_cut && jet_pt_corr < jet_pt_max_cut)
		{
		  count_wcut ++; // to count jet after cut
		  
		  // jet pT weight in MC
		  double Reco_jetpT_W = 1.; // for data 1.
		  
		  if(is_MC && is_ptWeight)
		    {
		      Reco_jetpT_W = 1./(fptWeight->Eval(jet_pt_corr));
		    }

		  // Filled corrected reco/data jet THnSparse
		  double Jet_CorrpT_Eta_Phi_ctbin_flavour[5] = {jet_pt_corr, jet_eta, jet_phi, (double)ctbin, (double)refpartonB};
		  hJet_CorrpT_Eta_Phi_ctbin_flavour_pTCut_W->Fill(Jet_CorrpT_Eta_Phi_ctbin_flavour, (Evtw)*(Reco_jetpT_W));
		  
		  if(is_MC)
		    {
		      if(ref_jet_pt > 0. && jet_pt_raw > 0. && jet_pt_corr > 0.)  // discrad the unmacthed jet (refpt == -999.)
			{	      
			  // Filled reco matched corrected jet THnSparse
			  double Jet_Matched_CorrpT_Eta_Phi_ctbin_flavour[5] = {jet_pt_corr, jet_eta, jet_phi, (double)ctbin, (double)refpartonB};
			  hJet_Matched_CorrpT_Eta_Phi_ctbin_flavour_pTCut_W->Fill(Jet_Matched_CorrpT_Eta_Phi_ctbin_flavour, (Evtw)*(Reco_jetpT_W));
			}
		    }  // is_MC condition (matched)
		} // jet eta cut
	    } // jet pt cut
	} // reco/ref Jet loop end~~~~~~~~~~~~~~~~~~
    
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:Reco/Data jets loop end:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      // leading subleading jet conditions for reco/data
      //reco or data
      if(nref > 1) // condition for at least 2 jets in an event to find leading and subleading jets
	{
	  hEvents->AddBinContent(8,1); // after nref > 1 cut
	  
	  if(ld_CorrpT != -999. && ld_CorrIndex != -999 && sld_CorrpT != -999. && sld_CorrIndex != -999)       
	    {
	      double ld_jet_pt  = ld_CorrpT;
	      double ld_jet_eta = jteta[ld_CorrIndex];
	      double ld_jet_phi = jtphi[ld_CorrIndex];
	      double ld_refjet_pt = refpt[ld_CorrIndex];
	      int ld_jet_refparonB = (int)refparton_flavorForB[ld_CorrIndex];
	      
	      double sld_jet_pt  = sld_CorrpT;
	      double sld_jet_eta = jteta[sld_CorrIndex];
	      double sld_jet_phi = jtphi[sld_CorrIndex];
	      double sld_refjet_pt = refpt[sld_CorrIndex];
	      int sld_jet_refparonB = (int)refparton_flavorForB[sld_CorrIndex];

	      //std::cout<<ld_refjet_pt<<"  "<<ld_refpT<<"  "<<ld_jet_pt<<"  "<<sld_refjet_pt<<"  "<<sld_refpT<<"  "<<sld_jet_pt<<std::endl;
	      /*
	      // JEU correction
	      //if(!is_MC && is_JEU)
	      if((!is_MC || is_MC)  && is_JEU)
		{
		  // for ld jet
		  JEU_ld->SetJetPT(ld_jet_pt);
		  JEU_ld->SetJetEta(ld_jet_eta);
		  JEU_ld->SetJetPhi(ld_jet_phi);

		  // for sld jet
		  JEU_sld->SetJetPT(sld_jet_pt);
		  JEU_sld->SetJetEta(sld_jet_eta);
		  JEU_sld->SetJetPhi(sld_jet_phi);
		  
		  if(is_JEU_down && !is_JEU_up)
		    {
		      ld_jet_pt = ld_jet_pt*(1. - (JEU_ld->GetUncertainty().first)); // for ld jet
		      sld_jet_pt = sld_jet_pt*(1. - (JEU_sld->GetUncertainty().first)); // for sld jet
		    }
		  else if(is_JEU_up && !is_JEU_down)
		    {
		      ld_jet_pt = ld_jet_pt*(1. + (JEU_ld->GetUncertainty().second)); // for ld jet
		      sld_jet_pt = sld_jet_pt*(1. + (JEU_sld->GetUncertainty().second)); // for sld jet
		    }
		}
	      
	      if(is_MC && is_JER)
		{
		  // for leading jet
		  double resolution_factor_ld = 1.;
		  if(ld_jet_eta <= jet_eta_min_cut || ld_jet_eta >= jet_eta_max_cut)
		    {
		      resolution_factor_ld = 1.;
		    }
		  else
		    {
		      resolution_factor_ld = getJERSFFromEta(ld_jet_eta, SFCondition, colliding_system);
		    }
		  
		  double extraResolution_ld = TMath::Sqrt(TMath::Max(resolution_factor_ld*resolution_factor_ld - 1.0, 0.0)); // found jet resolution
		  double sigma_smear_ld  = extraResolution_ld*fJERWeight->Eval(ld_refjet_pt); // some % worst --> from JetMET
		  if(ld_refjet_pt < 0.) sigma_smear_ld = extraResolution_ld*fJERWeight->Eval(0.1);
		  if(ld_refjet_pt > 5020.) sigma_smear_ld = extraResolution_ld*fJERWeight->Eval(5019.9);
		  
		  gRandom->SetSeed(0);
		  
		  double JER_smear_ld = gRandom->Gaus(1.,sigma_smear_ld);
		  while( JER_smear_ld < 0 ){ JER_smear_ld = gRandom->Gaus(1.,sigma_smear_ld); }
		  ld_jet_pt = ld_jet_pt*JER_smear_ld;

		  // for subleading jet
		  double resolution_factor_sld = 1.;
		  if(sld_jet_eta <= jet_eta_min_cut || sld_jet_eta >= jet_eta_max_cut)
		    {
		      resolution_factor_sld = 1.;
		    }
		  else
		    {
		      resolution_factor_sld = getJERSFFromEta(sld_jet_eta, SFCondition, colliding_system);
		    }
		  
		  double extraResolution_sld = TMath::Sqrt(TMath::Max(resolution_factor_sld*resolution_factor_sld - 1.0, 0.0)); // found jet resolution
		  double sigma_smear_sld  = extraResolution_sld*fJERWeight->Eval(sld_refjet_pt); // some % worst --> from JetMET
		  if(sld_refjet_pt < 0.) sigma_smear_sld = extraResolution_sld*fJERWeight->Eval(0.1);
		  if(sld_refjet_pt > 5020.) sigma_smear_sld = extraResolution_sld*fJERWeight->Eval(5019.9);
		  
		  gRandom->SetSeed(0);
		  
		  double JER_smear_sld = gRandom->Gaus(1.,sigma_smear_sld);
		  while( JER_smear_sld < 0 ){ JER_smear_sld = gRandom->Gaus(1.,sigma_smear_sld); }
		  sld_jet_pt = sld_jet_pt*JER_smear_sld;
		}
	      */
	      
	      // to estimate JES and JER for ld and subld jet
	      if(is_MC && is_JES_JER)
		{
		  if(ld_jet_eta > jet_eta_min_cut && ld_jet_eta < jet_eta_max_cut)
		    {
		      int ld_jet_refparonB_jes = 0;
		      if(fabs(ld_jet_refparonB) >= 1 && fabs(ld_jet_refparonB) <= 6) {ld_jet_refparonB_jes = fabs(ld_jet_refparonB);}
		      else if(fabs(ld_jet_refparonB) == 21){ld_jet_refparonB_jes = 7;}
		      else ld_jet_refparonB_jes = 0;
		      
		      double JES_ratio_ldcorrpT_vs_ldrefpT = ld_jet_pt/ld_refjet_pt;
		      double Jet_ldCorrpT_ldrefpT_ctbin_flavour[4] = {JES_ratio_ldcorrpT_vs_ldrefpT, ld_refjet_pt, (double)ctbin, (double)ld_jet_refparonB_jes};
		      //hJer_Jes_ldCorrpT_ldrefpT_ctbin_flavour_W->Fill(Jet_ldCorrpT_ldrefpT_ctbin_flavour, (Evtw)*(Reco_ldjetpT_W));
		      hJer_Jes_ldCorrpT_ldrefpT_ctbin_flavour_W->Fill(Jet_ldCorrpT_ldrefpT_ctbin_flavour, ptHatw);
		    }
		  if(sld_jet_eta > jet_eta_min_cut && sld_jet_eta < jet_eta_max_cut)
		    {
		      int sld_jet_refparonB_jes = 0;
		      if(fabs(sld_jet_refparonB) >= 1 && fabs(sld_jet_refparonB) <= 6) {sld_jet_refparonB_jes = fabs(sld_jet_refparonB);}
		      else if(fabs(sld_jet_refparonB) == 21){sld_jet_refparonB_jes = 7;}
		      else sld_jet_refparonB_jes = 0;
		      
		      double JES_ratio_sldcorrpT_vs_sldrefpT = sld_jet_pt/sld_refjet_pt;
		      double Jet_sldCorrpT_sldrefpT_ctbin_flavour[4] = {JES_ratio_sldcorrpT_vs_sldrefpT, sld_refjet_pt, (double)ctbin, (double)sld_jet_refparonB_jes};
		      //hJer_Jes_sldCorrpT_sldrefpT_ctbin_flavour_W->Fill(Jet_sldCorrpT_sldrefpT_ctbin_flavour, (Evtw)*(Reco_sldjetpT_W));
		      hJer_Jes_sldCorrpT_sldrefpT_ctbin_flavour_W->Fill(Jet_sldCorrpT_sldrefpT_ctbin_flavour, ptHatw);
		    }
		}
	      
	      double Dphi_ld_sld = TVector2::Phi_mpi_pi(ld_jet_phi - sld_jet_phi);

	      if(fabs(Dphi_ld_sld) > leading_subleading_deltaphi_min)
		{
		  if((ld_jet_eta > jet_eta_min_cut && ld_jet_eta < jet_eta_max_cut) && (sld_jet_eta > jet_eta_min_cut && sld_jet_eta < jet_eta_max_cut))
		    {
		      if((ld_jet_pt > leading_pT_min_cut && ld_jet_pt < leading_pT_max_cut) && (sld_jet_pt > subleading_pT_min_cut && sld_jet_pt < subleading_pT_max_cut))
		      //if((ld_jet_pt > leading_pT_min_cut && ld_jet_pt < leading_pT_max_cut) && (sld_jet_pt > subleading_pT_min_cut && sld_jet_pt < subleading_pT_max_cut) && (ld_refjet_pt > leading_pT_min_cut && ld_refjet_pt < leading_pT_max_cut) && (sld_refjet_pt > leading_pT_min_cut && sld_refjet_pt < leading_pT_max_cut)) // for spil over correction
			{
			  // flavour determination of leading and subleading jets
			  int ld_jet_refparonBB = 0, ld_jet_refparonBB_PM = 0;
			  if(is_MC)
			    {
			      if(fabs(ld_jet_refparonB) >= 1 && fabs(ld_jet_refparonB) <= 6)
				{
				  ld_jet_refparonBB = fabs(ld_jet_refparonB);
				  
				  if(ld_jet_refparonB > 0)
				    {
				      ld_jet_refparonBB_PM = fabs(ld_jet_refparonB);
				    }
				  else if(ld_jet_refparonB < 0)
				    {
				      ld_jet_refparonBB_PM = fabs(ld_jet_refparonB)+8;
				    }
				}
			      else if(fabs(ld_jet_refparonB) == 21)
				{
				  ld_jet_refparonBB = 7;
				  ld_jet_refparonBB_PM = 7;
				}
			      else
				{
				  ld_jet_refparonBB = 0;
				  ld_jet_refparonBB_PM = 0;
				}
			    }
			  else // for data
			    {
			      ld_jet_refparonBB = 1;
			      ld_jet_refparonBB_PM = 1;
			    }
						  
			  int sld_jet_refparonBB = 0, sld_jet_refparonBB_PM = 0;
			  if(is_MC)
			    {
			      if(fabs(sld_jet_refparonB) >= 1 && fabs(sld_jet_refparonB) <= 6)
				{
				  sld_jet_refparonBB = fabs(sld_jet_refparonB);
				  
				  if(sld_jet_refparonB > 0)
				    {
				      sld_jet_refparonBB_PM = fabs(sld_jet_refparonB);
				    }
				  else if(sld_jet_refparonB < 0)
				    {
				      sld_jet_refparonBB_PM = fabs(sld_jet_refparonB)+8;
				    }
				}
			      else if(fabs(sld_jet_refparonB) == 21)
				{
				  sld_jet_refparonBB = 7;
				  sld_jet_refparonBB_PM = 7;
				}
			      else
				{
				  sld_jet_refparonBB = 0;
				  sld_jet_refparonBB_PM = 0;
				}
			    }
			  else // for data
			    {
			      sld_jet_refparonBB = 1;
			      sld_jet_refparonBB_PM = 1;
			    }

			  // Jet pT weight for MC
			  double Reco_ldjetpT_W = 1., Reco_sldjetpT_W = 1.; // initialize with 1 , for data 1.
			  if(is_MC && is_ptWeight)
			    {
			      Reco_ldjetpT_W = 1./(fptWeight->Eval(ld_jet_pt));
			      Reco_sldjetpT_W = 1./(fptWeight->Eval(sld_jet_pt));
			    }
			  
			  // Fill TVector3 for leading and subleading jets
			  TVector3 Reco_Filtered_ldJet_CorrpT, Reco_Filtered_sldJet_CorrpT;
			  Reco_Filtered_ldJet_CorrpT.SetPtEtaPhi(ld_jet_pt, ld_jet_eta, ld_jet_phi);
			  Reco_Filtered_sldJet_CorrpT.SetPtEtaPhi(sld_jet_pt, sld_jet_eta, sld_jet_phi);
			  
			  // Filled leading subleading reco/data corrected jet THnSparse
			  double Jet_ldCorrpT_Eta_Phi_ctbin_flavour[5] = {ld_jet_pt, ld_jet_eta, ld_jet_phi, (double)ctbin, (double)ld_jet_refparonBB};
			  hJet_ldCorrpT_Eta_Phi_ctbin_flavour_pTCut_W->Fill(Jet_ldCorrpT_Eta_Phi_ctbin_flavour, (Evtw)*(Reco_ldjetpT_W));
			  
			  double Jet_sldCorrpT_Eta_Phi_ctbin_flavour[5] = {sld_jet_pt, sld_jet_eta, sld_jet_phi, (double)ctbin, (double)sld_jet_refparonBB};
			  hJet_sldCorrpT_Eta_Phi_ctbin_flavour_pTCut_W->Fill(Jet_sldCorrpT_Eta_Phi_ctbin_flavour, (Evtw)*(Reco_sldjetpT_W));
			  
			  // Fill 1D leading and subleading jet vector for each event. Those event does not have dijet, it will fill null
			  
			  // Reco or data
			  //1D leading jet
			  Reco_Filtered_ldJet_CorrpT_vec_1D.push_back(Reco_Filtered_ldJet_CorrpT);
			  Reco_Filtered_ldrefpartonB_vec_1D.push_back(ld_jet_refparonBB);
			  Reco_Filtered_ldJetW_vec_1D.push_back(Reco_ldjetpT_W);
			  //1D subleading jet
			  Reco_Filtered_sldJet_CorrpT_vec_1D.push_back(Reco_Filtered_sldJet_CorrpT);
			  Reco_Filtered_sldrefpartonB_vec_1D.push_back(sld_jet_refparonBB);
			  Reco_Filtered_sldJetW_vec_1D.push_back(Reco_sldjetpT_W);
			  
			  //1D event vector
			  hZvtx_ldsld->Fill(vertexz, Evtw);
			  Reco_Evtw_vec_1D.push_back(Evtw);
			  Reco_Evtno_vec_1D.push_back(event);
			  Reco_EvtCount_vec_1D.push_back(evtcount_Reco);
			  Vertexz_vec_1D.push_back(vertexz);
			  HiBin_vec_1D.push_back(ctbin);
			  HiBinValue_vec_1D.push_back(centbin);

			  evtcount_Reco++;
			  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~select reco/data tracks for selected leading jets~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			  //Define 1D inclusive jet vector used to fill 2D inclusive jet vector
			  std::vector<TVector3> Reco_Filtered_InclTrk_pT_vec_1D;           // trk vec
			  std::vector<double> Reco_Filtered_InclTrkW_vec_1D;               // trk weight
			  std::vector<int> Reco_Filtered_InclTrkCharge_vec_1D;             // trk charge
			  std::vector<int> Reco_Filtered_InclTrkSube_vec_1D;               // trk sube

			  double NtrkCount = 0.;
			  
			  //start reco/data tracks loop
			  for(int irctrk = 0; irctrk < ntrk; irctrk++)
			    {
			      float trk_pt = (float)trkpt[irctrk];
			      float trk_eta = (float)trketa[irctrk];
			      float trk_phi = (float)trkphi[irctrk];
			      bool trk_hp = (bool)highpur[irctrk];
			      int trk_chg = (int)trkcharge[irctrk];
			      int trk_nhits = (int)trknhits[irctrk];
			      float trk_chi2 = (float)trkchi2[irctrk];
			      int trk_ndf = (int)trkndof[irctrk];
			      int trk_nlayers = (int)trknlayer[irctrk];
			      float trk_pterr = (float)trkpterr[irctrk];
			      float trk_dxy = (float)trkdcaxy[irctrk];
			      float trk_dxyerr = (float)trkdcaxyerr[irctrk];
			      float trk_dz = (float)trkdcaz[irctrk];
			      float trk_dzerr = (float)trkdcazerr[irctrk];
			      float ETCalo = (((float)pfEcal[irctrk] + (float)pfHcal[irctrk])/(TMath::CosH(trk_eta)));
			      int trk_algo; float trk_mva;
			      if(colliding_system == "PbPb")
				{
				  trk_algo = (int)trkalgo[irctrk];
				  trk_mva = (float)trkmva[irctrk];
				}
			      
			      if(!trk_hp) continue;
			      if(trk_chg == 0) continue;
			      if(fabs(trk_pterr/trk_pt) > trk_pt_resolution_cut) continue;
			      if(fabs(trk_dxy/trk_dxyerr) > trk_dca_xy_cut) continue;
			      if(fabs(trk_dz/trk_dzerr) > trk_dca_z_cut) continue;
			      if(colliding_system == "PbPb")
				{
				  if(trk_nhits < nhits) continue;
				  if(trk_algo == 6 && trk_mva < 0.98) continue;
				  if(trk_chi2/trk_ndf/trk_nlayers >= chi2_ndf_nlayer_cut) continue;
				  if(trk_pt > 20 && fabs(ETCalo/trk_pt) < calo_matching) continue;
				}
			      if(trk_pt <= trk_pt_min_cut || trk_pt >= trk_pt_max_cut) continue;
			      //if(trk_pt < trk_pt_min_cut) continue;
			      if(fabs(trk_eta) >= trk_eta_cut) continue;
			      
			      double trkeff = 1.;
			      double trkfak = 1.;
			      double trksec = 1.;
			      double trk_wt = 1.;
			      
			      if(colliding_system == "PbPb")
				{
				  if(is_MC)
				    {
				      //trkeff = htrk_eff_PbPb->GetBinContent(htrk_eff_PbPb->FindBin(trk_eta, trk_pt, centbin+10));
				      //trkfak = htrk_fak_PbPb->GetBinContent(htrk_fak_PbPb->FindBin(trk_eta, trk_pt, centbin+10));
				      trkeff = htrk_eff_PbPb->GetBinContent(htrk_eff_PbPb->FindBin(trk_eta, trk_pt, centbin));
                                      trkfak = htrk_fak_PbPb->GetBinContent(htrk_fak_PbPb->FindBin(trk_eta, trk_pt, centbin));
				    }
				  else
				    {
				      trkeff = htrk_eff_PbPb->GetBinContent(htrk_eff_PbPb->FindBin(trk_eta, trk_pt, centbin));
				      trkfak = htrk_fak_PbPb->GetBinContent(htrk_fak_PbPb->FindBin(trk_eta, trk_pt, centbin));
				    }
				}
			      else if(colliding_system == "pp")
				{
				  trkeff = (htrk_eff_pp->GetBinContent(htrk_eff_pp->FindBin(trk_eta, trk_pt)))*0.979; //0.979 is scale factor from D mesons
				  trkfak = htrk_fak_pp->GetBinContent(htrk_fak_pp->FindBin(trk_eta, trk_pt));
				  trksec = htrk_sec_pp->GetBinContent(htrk_sec_pp->FindBin(trk_eta, trk_pt));
				}
			      if(trkeff > 0.001)
				{
				  if(colliding_system == "PbPb")
				    {
				      trk_wt = (1. - trkfak)/trkeff;
				    }
				  else if(colliding_system == "pp")
				    {
				      //trk_wt = (1. - trkfak)/trkeff;
				      trk_wt = ((1. - trkfak)*(1. - trksec))/trkeff;
				    }
				}
			      else trk_wt = 1.;
			      if((trk_pt < 0 || trk_pt > 500) || (fabs(trk_eta) > 2.4)) trk_wt = 0.; // check bound
			      //if(trk_pt < 0 || fabs(trk_eta) > 2.4) trk_wt = 0.; // check bound
			      //if(trk_pt > 500) trk_wt = 1.;
			      
			      //std::cout<<"Trk pt, eta, and phi, eff : "<<trk_pt<<"  "<<trk_eta<<"  "<<trk_phi<<"  "<<trk_wt<<std::endl;

			      NtrkCount += trk_wt; // corrected ntrk count
			      
			      TVector3 Reco_Filtered_InclTrk_pT;
			      Reco_Filtered_InclTrk_pT.SetPtEtaPhi(trk_pt, trk_eta, trk_phi);
			      
			      Reco_Filtered_InclTrk_pT_vec_1D.push_back(Reco_Filtered_InclTrk_pT);
			      Reco_Filtered_InclTrkW_vec_1D.push_back(trk_wt);
			      Reco_Filtered_InclTrkCharge_vec_1D.push_back(trk_chg);
			      Reco_Filtered_InclTrkSube_vec_1D.push_back(1);
			      
			      double Trk_pT_Eta_Phi_ctbin[4] = {trk_pt, trk_eta, trk_phi, (double)ctbin};
			      hTrk_pT_Eta_Phi_ctbin_pTCut_W->Fill(Trk_pT_Eta_Phi_ctbin, 1.);
			      hTrk_CorrpT_Eta_Phi_ctbin_pTCut_W->Fill(Trk_pT_Eta_Phi_ctbin, trk_wt);
			    } // reco track loop end

			  // Fill ntrk count for check
			  if(ctbin == 0) hntrk_Signal_Check_0->Fill(NtrkCount, Evtw);
			  else if(ctbin == 1) hntrk_Signal_Check_1->Fill(NtrkCount, Evtw);
			  else if(ctbin == 2) hntrk_Signal_Check_2->Fill(NtrkCount, Evtw);
			  else if(ctbin == 3) hntrk_Signal_Check_3->Fill(NtrkCount, Evtw);
			  else if(ctbin == 4) hntrk_Signal_Check_4->Fill(NtrkCount, Evtw);
			  
			  // filled reco/data 2D trk vector
			  Reco_Filtered_InclTrk_pT_vec_2D.push_back(Reco_Filtered_InclTrk_pT_vec_1D);
			  Reco_Filtered_InclTrkW_vec_2D.push_back(Reco_Filtered_InclTrkW_vec_1D);
			  Reco_Filtered_InclTrkCharge_vec_2D.push_back(Reco_Filtered_InclTrkCharge_vec_1D);
			  Reco_Filtered_InclTrkSube_vec_2D.push_back(Reco_Filtered_InclTrkSube_vec_1D);
			  
			  // clear reco/data 1D trk vector used to fill 2D vector
			  Reco_Filtered_InclTrk_pT_vec_1D.clear();
			  Reco_Filtered_InclTrkW_vec_1D.clear();
			  Reco_Filtered_InclTrkCharge_vec_1D.clear();
			  Reco_Filtered_InclTrkSube_vec_1D.clear();

			  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~select gen tracks for selected leading reco/data jets~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			  // for reco jets and gen tracks
			  std::vector<TVector3> RecoGen_Filtered_InclTrk_pT_vec_1D;           // trk vec
			  std::vector<double> RecoGen_Filtered_InclTrkW_vec_1D;               // trk weight
			  std::vector<int> RecoGen_Filtered_InclTrkCharge_vec_1D;             // trk charge
			  std::vector<int> RecoGen_Filtered_InclTrkSube_vec_1D;               // trk sube
			  
			  if(isRcJetGnTrk)
			    {
			      //start gen track loop in reco jet loop                                                                                                  
			      for(int ircgntrk = 0; ircgntrk < gen_trkpt->size(); ircgntrk++)
				{
				  float recogentrk_pt = float((*gen_trkpt)[ircgntrk]);
				  float recogentrk_eta = float((*gen_trketa)[ircgntrk]);
				  float recogentrk_phi = float((*gen_trkphi)[ircgntrk]);
				  int recogentrk_chg = int((*gen_trkchg)[ircgntrk]);
				  int recogentrk_sube = int((*gen_trksube)[ircgntrk]);
				  
				  if(recogentrk_chg == 0) continue;
				  if(recogentrk_pt <= trk_pt_min_cut || recogentrk_pt >= trk_pt_max_cut) continue;
				  //if(recogentrk_pt < trk_pt_min_cut) continue;
				  if(fabs(recogentrk_eta) >= trk_eta_cut) continue;
				  
				  //std::cout<<"recogentrk_sube is: "<<recogentrk_sube<<std::endl;
				  double recogentrk_wt = 1.;
				  
				  TVector3 RecoGen_Filtered_InclTrk_pT;
				  RecoGen_Filtered_InclTrk_pT.SetPtEtaPhi(recogentrk_pt, recogentrk_eta, recogentrk_phi);
				  
				  RecoGen_Filtered_InclTrk_pT_vec_1D.push_back(RecoGen_Filtered_InclTrk_pT);
				  RecoGen_Filtered_InclTrkW_vec_1D.push_back(recogentrk_wt);
				  RecoGen_Filtered_InclTrkCharge_vec_1D.push_back(recogentrk_chg);
				  RecoGen_Filtered_InclTrkSube_vec_1D.push_back(recogentrk_sube);
				}
			      // filled gen 2D trk vector
			      RecoGen_Filtered_InclTrk_pT_vec_2D.push_back(RecoGen_Filtered_InclTrk_pT_vec_1D);
			      RecoGen_Filtered_InclTrkW_vec_2D.push_back(RecoGen_Filtered_InclTrkW_vec_1D);
			      RecoGen_Filtered_InclTrkCharge_vec_2D.push_back(RecoGen_Filtered_InclTrkCharge_vec_1D);
			      RecoGen_Filtered_InclTrkSube_vec_2D.push_back(RecoGen_Filtered_InclTrkSube_vec_1D);
			      // clear gen 1D trk vector used to fill 2D vector
			      RecoGen_Filtered_InclTrk_pT_vec_1D.clear();
			      RecoGen_Filtered_InclTrkW_vec_1D.clear();
			      RecoGen_Filtered_InclTrkCharge_vec_1D.clear();
			      RecoGen_Filtered_InclTrkSube_vec_1D.clear();
			    } // if(isRcJetGnTrk)
			} // leading subleadin phi condition
		    } // leading subleading jet pt cut
		} // leading subleading jet eta cut
	    } // -999 condition
      	} // nref > 1 condition
    
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:Reco/Data tracks loop end:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:Event loop continue:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:Gen jet loop start:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // gen jet loop
      if(is_MC)
        {
	  double ld_GenpT=-999.; int ld_GenIndex=-999;     // leading gen jet quantities
	  double sld_GenpT=-999.; int sld_GenIndex=-999;   // subleading gen jet quantities   
	  int ld_RefIndex = -999, sld_RefIndex = -999;     // leading subleading ref index
	  
	  for (int ign = 0; ign < (int)ngen; ign++) //gen Jet loop start                                        
            {
	      int index_ref = -1;
	      int index_gen = -1;
	      
	      double min_Dpt = 999999999.;
	      
	      for (int irf = 0; irf < (int)nref; irf++) //ref Jet loop start to calculate the parton flavour
		{
		  //if(fabs(gen_jtpt[ign] - refpt[irf]) < std::numeric_limits<float>::epsilon()) // match gen jet with ref jet to find ref eta and phi
		  if(fabs(gen_jtpt[ign] - refpt[irf]) < min_Dpt)
		    {
		      min_Dpt = fabs(gen_jtpt[ign] - refpt[irf]);
		      index_ref = irf;
		      index_gen = ign;
		    }
		}
	      
	      if(index_ref >= 0 && index_gen >= 0)
		{
		  double gen_jet_pt = gen_jtpt[index_gen];
		  double gen_jet_eta = gen_jteta[index_gen];
		  double gen_jet_phi = gen_jtphi[index_gen];
		  
		  int refpartonB_gn = 0, refpartonB_gn_PM = 0;
		  if(fabs(refparton_flavorForB[index_ref]) >= 1 && fabs(refparton_flavorForB[index_ref]) <= 6)
		    {
		      refpartonB_gn = fabs(refparton_flavorForB[index_ref]);
		      
		      if(refparton_flavorForB[index_ref] > 0)
			{
			  refpartonB_gn_PM = fabs(refparton_flavorForB[index_ref]);
			}
		      else if(refparton_flavorForB[index_ref] < 0)
			{
			  refpartonB_gn_PM = fabs(refparton_flavorForB[index_ref])+8;
			}
		    }
		  else if(fabs(refparton_flavorForB[index_ref]) == 21)
		    {
		      refpartonB_gn = 7;
		      refpartonB_gn_PM = 7;
		    }
		  else
		    {
		      refpartonB_gn = 0;
		      refpartonB_gn_PM = 0;
		    }
		  
		  // find leading subleading gen jets
		  //find_leading_subleading_jets(gen_jet_pt, index_gen, ld_GenpT, ld_GenIndex, sld_GenpT, sld_GenIndex);   // for gen jet pt
		  find_leading_subleading_gen_jets(gen_jet_pt, index_gen, index_ref, ld_GenpT, ld_GenIndex, ld_RefIndex, sld_GenpT, sld_GenIndex, sld_RefIndex);
		  
		  if(gen_jet_eta > jet_eta_min_cut && gen_jet_eta < jet_eta_max_cut)
		    {
		      if(gen_jet_pt > jet_pt_min_cut && gen_jet_pt < jet_pt_max_cut)
			{
			  // jet pT weight in MC
			  double Gen_jetpT_W = 1.;
			  if(is_MC && is_ptWeight)
			    {
			      Gen_jetpT_W = 1./(fptWeight->Eval(gen_jet_pt));
			    }
			  // Filled corrected gen jet THnSparse
			  double Jet_GenpT_Eta_Phi_ctbin_flavour[5] = {gen_jet_pt, gen_jet_eta, gen_jet_phi, (double)ctbin, (double)refpartonB_gn};
			  hJet_GenpT_Eta_Phi_ctbin_flavour_pTCut_W->Fill(Jet_GenpT_Eta_Phi_ctbin_flavour, (genEvtw)*(Gen_jetpT_W));
			} // gen jet pt cut
		    } // gen jet eta cut
		} // if(index_ref >= 0 && index_gen >= 0)
	    } // gen jet loop end~~~~~~~~~~~~~~~~~~~~~~~~~
	
	  //filling leading and subleading gen jet quantities
	  	  
	  if(ngen > 1) // condition for at least 2 jets in an event to find leading and subleading jets                    
	    {
	      //if(ld_GenpT != -999. && ld_GenIndex != -999 && sld_GenpT != -999. && sld_GenIndex != -999)
	      if(ld_GenpT != -999. && ld_GenIndex != -999 && sld_GenpT != -999. && sld_GenIndex != -999 && ld_RefIndex != -999 && sld_RefIndex != -999)       
		{
		  double ld_Genjet_pt  = ld_GenpT;
		  double ld_Genjet_eta = gen_jteta[ld_GenIndex];
		  double ld_Genjet_phi = gen_jtphi[ld_GenIndex];
		  int ld_Genjet_refparonB = (int)refparton_flavorForB[ld_RefIndex];
		  
		  double sld_Genjet_pt  = sld_GenpT;
		  double sld_Genjet_eta = gen_jteta[sld_GenIndex];
		  double sld_Genjet_phi = gen_jtphi[sld_GenIndex];
		  int sld_Genjet_refparonB = (int)refparton_flavorForB[sld_RefIndex];
		  
		  double Dphi_ld_sld_Gen = TVector2::Phi_mpi_pi(ld_Genjet_phi - sld_Genjet_phi);
		  
		  if(fabs(Dphi_ld_sld_Gen) > leading_subleading_deltaphi_min)
		    {
		      if((ld_Genjet_eta > jet_eta_min_cut && ld_Genjet_eta < jet_eta_max_cut) && (sld_Genjet_eta > jet_eta_min_cut && sld_Genjet_eta < jet_eta_max_cut))
			{
			  if((ld_Genjet_pt > leading_pT_min_cut && ld_Genjet_pt < leading_pT_max_cut) && (sld_Genjet_pt > subleading_pT_min_cut && sld_Genjet_pt < subleading_pT_max_cut))
			    {
			      // determine the flavour of leading and subleadit jets
			      int ld_Genjet_refparonBB = 0, ld_Genjet_refparonBB_PM = 0;
			      if(fabs(ld_Genjet_refparonB) >= 1 && fabs(ld_Genjet_refparonB) <= 6)
				{
				  ld_Genjet_refparonBB = fabs(ld_Genjet_refparonB);
				  
				  if(ld_Genjet_refparonB > 0)
				    {
				      ld_Genjet_refparonBB_PM = fabs(ld_Genjet_refparonB);
				    }
				  else if(ld_Genjet_refparonB < 0)
				    {
				      ld_Genjet_refparonBB_PM = fabs(ld_Genjet_refparonB)+8;
				    }
				}
			      else if(fabs(ld_Genjet_refparonB) == 21)
				{
				  ld_Genjet_refparonBB = 7;
				  ld_Genjet_refparonBB_PM = 7;
				}
			      else
				{
				  ld_Genjet_refparonBB = 0;
				  ld_Genjet_refparonBB_PM = 0;
				}
			      
			      int sld_Genjet_refparonBB = 0, sld_Genjet_refparonBB_PM = 0;
			      if(fabs(sld_Genjet_refparonB) >= 1 && fabs(sld_Genjet_refparonB) <= 6)
				{
				  sld_Genjet_refparonBB = fabs(sld_Genjet_refparonB);
				  
				  if(sld_Genjet_refparonB > 0)
				    {
				      sld_Genjet_refparonBB_PM = fabs(sld_Genjet_refparonB);
				    }
				  else if(sld_Genjet_refparonB < 0)
				    {
				      sld_Genjet_refparonBB_PM = fabs(sld_Genjet_refparonB)+8;
				    }
				}
			      else if(fabs(sld_Genjet_refparonB) == 21)
				{
				  sld_Genjet_refparonBB = 7;
				  sld_Genjet_refparonBB_PM = 7;
				}
			      else
				{
				  sld_Genjet_refparonBB = 0;
				  sld_Genjet_refparonBB_PM = 0;
				}
			      
			      // ldjet and sldjet gen jet pT weight in MC
			      double Gen_ldjetpT_W = 1., Gen_sldjetpT_W = 1.; // Initialize with 1
			      if(is_MC && is_ptWeight)
				{
				  Gen_ldjetpT_W = 1./(fptWeight->Eval(ld_Genjet_pt));
				  Gen_sldjetpT_W = 1./(fptWeight->Eval(sld_Genjet_pt));
				}
			      
			      // fill TVector3
			      TVector3 Gen_Filtered_ldJet_CorrpT, Gen_Filtered_sldJet_CorrpT;
			      Gen_Filtered_ldJet_CorrpT.SetPtEtaPhi(ld_Genjet_pt, ld_Genjet_eta, ld_Genjet_phi);
			      Gen_Filtered_sldJet_CorrpT.SetPtEtaPhi(sld_Genjet_pt, sld_Genjet_eta, sld_Genjet_phi);
			      
			      // Filled leading subleading gen jet THnSparse
			      double Jet_ldGenpT_Eta_Phi_ctbin_flavour[5] = {ld_Genjet_pt, ld_Genjet_eta, ld_Genjet_phi, (double)ctbin, (double)ld_Genjet_refparonBB};
			      hJet_ldGenpT_Eta_Phi_ctbin_flavour_pTCut_W->Fill(Jet_ldGenpT_Eta_Phi_ctbin_flavour, (genEvtw)*(Gen_ldjetpT_W));
			      
			      double Jet_sldGenpT_Eta_Phi_ctbin_flavour[5] = {sld_Genjet_pt, sld_Genjet_eta, sld_Genjet_phi, (double)ctbin, (double)sld_Genjet_refparonBB};
			      hJet_sldGenpT_Eta_Phi_ctbin_flavour_pTCut_W->Fill(Jet_sldGenpT_Eta_Phi_ctbin_flavour, (genEvtw)*(Gen_sldjetpT_W));
			      
			      
			      // Fill 1D leading and subleading jet vector for each event. Those event does not have dijet, it will fill null
			      //1D leading jets
			      Gen_Filtered_ldJet_CorrpT_vec_1D.push_back(Gen_Filtered_ldJet_CorrpT);
			      Gen_Filtered_ldrefpartonB_vec_1D.push_back(ld_Genjet_refparonBB);
			      Gen_Filtered_ldJetW_vec_1D.push_back(Gen_ldjetpT_W);
			      
			      //1D subleading jets
			      Gen_Filtered_sldJet_CorrpT_vec_1D.push_back(Gen_Filtered_sldJet_CorrpT);
			      Gen_Filtered_sldrefpartonB_vec_1D.push_back(sld_Genjet_refparonBB);
			      Gen_Filtered_sldJetW_vec_1D.push_back(Gen_ldjetpT_W);
			      
			      // Filled event quantities vector
			      hZvtx_ldsld_Gen->Fill(vertexz, Evtw);
			      Gen_Evtw_vec_1D.push_back(genEvtw);
			      Gen_Evtno_vec_1D.push_back(event);
			      Gen_EvtCount_vec_1D.push_back(evtcount_Gen);
			      Gen_Vertexz_vec_1D.push_back(vertexz);
			      Gen_HiBin_vec_1D.push_back(ctbin);
			      Gen_HiBinValue_vec_1D.push_back(centbin);

			      evtcount_Gen++;
			      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:Gen jet loop end:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			      
			      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~select gen tracks for selected leading jets~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			      //Define 1D tracks vectors used to fill 2D trk vectors
			      std::vector<TVector3> Gen_Filtered_InclTrk_pT_vec_1D;           // trk vec
			      std::vector<double> Gen_Filtered_InclTrkW_vec_1D;               // trk weight
			      std::vector<int> Gen_Filtered_InclTrkCharge_vec_1D;             // trk charge
			      std::vector<int> Gen_Filtered_InclTrkSube_vec_1D;             // trk sube
			      
			      //start gen track loop                                                                                                  
			      for(int igntrk = 0; igntrk < gen_trkpt->size(); igntrk++)
				{
				  float gentrk_pt = float((*gen_trkpt)[igntrk]);
				  float gentrk_eta = float((*gen_trketa)[igntrk]);
				  float gentrk_phi = float((*gen_trkphi)[igntrk]);
				  int gentrk_chg = int((*gen_trkchg)[igntrk]);
				  int gentrk_sube = int((*gen_trksube)[igntrk]);

				  if(gentrk_chg == 0) continue;
				  if(gentrk_pt <= trk_pt_min_cut || gentrk_pt >= trk_pt_max_cut) continue;
				  //if(gentrk_pt < trk_pt_min_cut) continue;
				  if(fabs(gentrk_eta) >= trk_eta_cut) continue;
				  
				  //std::cout<<"gentrk_sube is: "<<gentrk_sube<<std::endl;
				  double gentrk_wt = 1.;
				  
				  TVector3 Gen_Filtered_InclTrk_pT;
				  Gen_Filtered_InclTrk_pT.SetPtEtaPhi(gentrk_pt, gentrk_eta, gentrk_phi);
				  
				  Gen_Filtered_InclTrk_pT_vec_1D.push_back(Gen_Filtered_InclTrk_pT);
				  Gen_Filtered_InclTrkW_vec_1D.push_back(gentrk_wt);
				  Gen_Filtered_InclTrkCharge_vec_1D.push_back(gentrk_chg);
				  Gen_Filtered_InclTrkSube_vec_1D.push_back(gentrk_sube);
				    
				  double Trk_GenpT_Eta_Phi_ctbin[4] = {gentrk_pt, gentrk_eta, gentrk_phi, (double)ctbin};
				  hTrk_GenpT_Eta_Phi_ctbin_pTCut_W->Fill(Trk_GenpT_Eta_Phi_ctbin, gentrk_wt);
				}
			      
			      // filled gen 2D trk vector
			      Gen_Filtered_InclTrk_pT_vec_2D.push_back(Gen_Filtered_InclTrk_pT_vec_1D);
			      Gen_Filtered_InclTrkW_vec_2D.push_back(Gen_Filtered_InclTrkW_vec_1D);
			      Gen_Filtered_InclTrkCharge_vec_2D.push_back(Gen_Filtered_InclTrkCharge_vec_1D);
			      Gen_Filtered_InclTrkSube_vec_2D.push_back(Gen_Filtered_InclTrkSube_vec_1D);
			      
			      // clear gen 1D trk vector used to fill 2D vector
			      Gen_Filtered_InclTrk_pT_vec_1D.clear();
			      Gen_Filtered_InclTrkW_vec_1D.clear();
			      Gen_Filtered_InclTrkCharge_vec_1D.clear();
			      Gen_Filtered_InclTrkSube_vec_1D.clear();
			      
			      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~select reco tracks for selected leading gen jets~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			      // for gen jets and reco tracks
			      std::vector<TVector3> GenReco_Filtered_InclTrk_pT_vec_1D;           // trk vec
			      std::vector<double> GenReco_Filtered_InclTrkW_vec_1D;               // trk weight
			      std::vector<int> GenReco_Filtered_InclTrkCharge_vec_1D;             // trk charge
			      std::vector<int> GenReco_Filtered_InclTrkSube_vec_1D;             // trk sube

			      if(isRcJetGnTrk)
				{
				  //start reco/data tracks loop in gen jet loop
				  for(int ignrctrk = 0; ignrctrk < ntrk; ignrctrk++)
				    {
				      float genrecotrk_pt = (float)trkpt[ignrctrk];
				      float genrecotrk_eta = (float)trketa[ignrctrk];
				      float genrecotrk_phi = (float)trkphi[ignrctrk];
				      bool genrecotrk_hp = (bool)highpur[ignrctrk];
				      int genrecotrk_chg = (int)trkcharge[ignrctrk];
				      int genrecotrk_nhits = (int)trknhits[ignrctrk];
				      float genrecotrk_chi2 = (float)trkchi2[ignrctrk];
				      int genrecotrk_ndf = (int)trkndof[ignrctrk];
				      int genrecotrk_nlayers = (int)trknlayer[ignrctrk];
				      float genrecotrk_pterr = (float)trkpterr[ignrctrk];
				      float genrecotrk_dxy = (float)trkdcaxy[ignrctrk];
				      float genrecotrk_dxyerr = (float)trkdcaxyerr[ignrctrk];
				      float genrecotrk_dz = (float)trkdcaz[ignrctrk];
				      float genrecotrk_dzerr = (float)trkdcazerr[ignrctrk];
				      float ETCalo = (((float)pfEcal[ignrctrk] + (float)pfHcal[ignrctrk])/(TMath::CosH(genrecotrk_eta)));
				      int genrecotrk_algo; float genrecotrk_mva;
				      if(colliding_system == "PbPb")
					{
					  genrecotrk_algo = (int)trkalgo[ignrctrk];
					  genrecotrk_mva = (float)trkmva[ignrctrk];
					}
				      
				      if(!genrecotrk_hp) continue;
				      if(genrecotrk_chg == 0) continue;
				      if(fabs(genrecotrk_pterr/genrecotrk_pt) > trk_pt_resolution_cut) continue;
				      if(fabs(genrecotrk_dxy/genrecotrk_dxyerr) > trk_dca_xy_cut) continue;
				      if(fabs(genrecotrk_dz/genrecotrk_dzerr) > trk_dca_z_cut) continue;
				      if(colliding_system == "PbPb")
					{
					  if(genrecotrk_nhits < nhits) continue;
					  if(genrecotrk_algo == 6 && genrecotrk_mva < 0.98) continue;
					  if(genrecotrk_chi2/genrecotrk_ndf/genrecotrk_nlayers >= chi2_ndf_nlayer_cut) continue;
					  if(genrecotrk_pt > 20 && fabs(ETCalo/genrecotrk_pt) < calo_matching) continue;
					}
				      if(genrecotrk_pt <= trk_pt_min_cut || genrecotrk_pt >= trk_pt_max_cut) continue;
				      //if(genrecotrk_pt < trk_pt_min_cut) continue;
				      if(fabs(genrecotrk_eta) >= trk_eta_cut) continue;
				      
				      double genrecotrkeff = 1.;
				      double genrecotrkfak = 1.;
				      double genrecotrksec = 1.;
				      double genrecotrk_wt = 1.;
				      
				      if(colliding_system == "PbPb")
					{
					  if(is_MC)
					    {
					      //genrecotrkeff = htrk_eff_PbPb->GetBinContent(htrk_eff_PbPb->FindBin(genrecotrk_eta, genrecotrk_pt, centbin+10));
					      //genrecotrkfak = htrk_fak_PbPb->GetBinContent(htrk_fak_PbPb->FindBin(genrecotrk_eta, genrecotrk_pt, centbin+10));
					      genrecotrkeff = htrk_eff_PbPb->GetBinContent(htrk_eff_PbPb->FindBin(genrecotrk_eta, genrecotrk_pt, centbin));
                                              genrecotrkfak = htrk_fak_PbPb->GetBinContent(htrk_fak_PbPb->FindBin(genrecotrk_eta, genrecotrk_pt, centbin));
					    }
					  else
					    {
					      genrecotrkeff = htrk_eff_PbPb->GetBinContent(htrk_eff_PbPb->FindBin(genrecotrk_eta, genrecotrk_pt, centbin));
					      genrecotrkfak = htrk_fak_PbPb->GetBinContent(htrk_fak_PbPb->FindBin(genrecotrk_eta, genrecotrk_pt, centbin));
					    }
					}
				      else if(colliding_system == "pp")
					{
					  genrecotrkeff = (htrk_eff_pp->GetBinContent(htrk_eff_pp->FindBin(genrecotrk_eta, genrecotrk_pt)))*0.979; //0.979 is scale factor from D mesons
					  genrecotrkfak = htrk_fak_pp->GetBinContent(htrk_fak_pp->FindBin(genrecotrk_eta, genrecotrk_pt));
					  genrecotrksec = htrk_sec_pp->GetBinContent(htrk_sec_pp->FindBin(genrecotrk_eta, genrecotrk_pt));
					}
				      if(genrecotrkeff > 0.001)
					{
					  if(colliding_system == "PbPb")
					    {
					      genrecotrk_wt = (1. - genrecotrkfak)/genrecotrkeff;
					    }
					  else if(colliding_system == "pp")
					    {
					      //genrecotrk_wt = (1. - trkfak)/trkeff;
					      genrecotrk_wt = ((1. - genrecotrkfak)*(1. - genrecotrksec))/genrecotrkeff;
					    }
					}
				      else genrecotrk_wt = 1.;
				      if((genrecotrk_pt < 0 || genrecotrk_pt > 500) || (fabs(genrecotrk_eta) > 2.4)) genrecotrk_wt = 0.; // check bound
				      //if(genrecotrk_pt < 0 || fabs(genrecotrk_eta) > 2.4) genrecotrk_wt = 0.; // check bound
				      //if(genrecotrk_pt > 500) genrecotrk_wt = 1.;
				      
				      TVector3 GenReco_Filtered_InclTrk_pT;
				      GenReco_Filtered_InclTrk_pT.SetPtEtaPhi(genrecotrk_pt, genrecotrk_eta, genrecotrk_phi);
				      
				      GenReco_Filtered_InclTrk_pT_vec_1D.push_back(GenReco_Filtered_InclTrk_pT);
				      GenReco_Filtered_InclTrkW_vec_1D.push_back(genrecotrk_wt);
				      GenReco_Filtered_InclTrkCharge_vec_1D.push_back(genrecotrk_chg);
				      GenReco_Filtered_InclTrkSube_vec_1D.push_back(1);
				    } // reco track loop end
				  // filled reco/data 2D trk vector
				  GenReco_Filtered_InclTrk_pT_vec_2D.push_back(GenReco_Filtered_InclTrk_pT_vec_1D);
				  GenReco_Filtered_InclTrkW_vec_2D.push_back(GenReco_Filtered_InclTrkW_vec_1D);
				  GenReco_Filtered_InclTrkCharge_vec_2D.push_back(GenReco_Filtered_InclTrkCharge_vec_1D);
				  GenReco_Filtered_InclTrkSube_vec_2D.push_back(GenReco_Filtered_InclTrkSube_vec_1D);
				  
				  // clear reco/data 1D trk vector used to fill 2D vector
				  GenReco_Filtered_InclTrk_pT_vec_1D.clear();
				  GenReco_Filtered_InclTrkW_vec_1D.clear();
				  GenReco_Filtered_InclTrkCharge_vec_1D.clear();
				  GenReco_Filtered_InclTrkSube_vec_1D.clear();
				} //if(isRcJetGnTrk)
			    } // Gen leading subleadin phi condition
			} // Gen leading subleading jet pt cut
		    } // Gen matchd leading subleading jet eta cut
		} // Gen -999 condition
	    } //if(ngen > 1)
	} // is_MC condition
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:Gen tracks loop end:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    } // events loop end

  std::cout<<"Total nref and nref with cuts are: "<<count<<"  "<<count_wcut<<std::endl;
  std::cout<<std::endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:Events loop end:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(colliding_system == "PbPb" && isMBEventsMixing)
    {
      std::cout<<"~~~~~~~~~~~Now started for MinimunBias(MB) events~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
      std::cout<<std::endl;
    }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~start event loop for MB events used for mixing~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  //~~~~~~~~~~~~Define vectors use for jet tracks Mixed event correlation from MB sample~~~~~~~~~~~~~~~~~~~~~~
  // 1D event vectors
  std::vector<double> Vertexz_vec_MB_1D;                   // vertex z
  std::vector<Long64_t> Evtno_vec_MB_1D;                   // event number
  std::vector<int> HiBin_vec_MB_1D;                        // centrality
  std::vector<int> HiBinValue_vec_MB_1D;                   // centrality
  std::vector<int> EvtCount_vec_MB_1D;                        // centrality
  
  //2D tracks vectors (Event wise and track wise) 
  //gen vector
  std::vector<std::vector<TVector3>> Gen_Filtered_InclTrk_pT_vec_MB_2D;            // trk vec
  std::vector<std::vector<double>> Gen_Filtered_InclTrkW_vec_MB_2D;                // trk weight
  std::vector<std::vector<int>> Gen_Filtered_InclTrkSube_vec_MB_2D;                // trk sube
  
  //reco/data vector
  std::vector<std::vector<TVector3>> Reco_Filtered_InclTrk_pT_vec_MB_2D;           // trk vec
  std::vector<std::vector<double>> Reco_Filtered_InclTrkW_vec_MB_2D;               // trk weight
  std::vector<std::vector<int>> Reco_Filtered_InclTrkSube_vec_MB_2D;                // trk sube
  //~~~~~~~~~~~~End Defining vectors use for jet track Mixed event correlation from MB sample~~~~~~~~~~~~~~~~~~~~~~

  TChain *hea_tree_MB = nullptr;
  TChain *trk_tree_MB = nullptr;
  TChain *gen_tree_MB = nullptr;
  
  if(colliding_system == "PbPb" && isMBEventsMixing)
    {
      //  open input forest/skim file
      TString input_file_MB;
      if(is_MC) input_file_MB = inputFileMC_MB;
      else if(!is_MC) input_file_MB  = inputFileData_MB;
      
      fstream openInputFile_MB;
      openInputFile_MB.open(Form("%s",input_file_MB.Data()), ios::in);
      if(!openInputFile_MB.is_open())
	{
	  cout << "List of input MB files not founded!" << endl;
	  return;
	}
      
      // Make a chain and a vector of file names
      std::vector<TString> file_name_vector_MB;
      string file_chain_MB;
      while(getline(openInputFile_MB, file_chain_MB))
	{
	  //file_name_vector.push_back(Form("root://cmsxrootd.fnal.gov/%s", file_chain.c_str()));
	  file_name_vector_MB.push_back(file_chain_MB.c_str());
	}
      openInputFile_MB.close();
      
      hea_tree_MB = new TChain("hiEvtAnalyzer/HiTree");
      trk_tree_MB = new TChain("ppTrack/trackTree");
      if(is_MC) gen_tree_MB = new TChain("HiGenParticleAna/hi");
      
      for (std::vector<TString>::iterator listIterator_MB = file_name_vector_MB.begin(); listIterator_MB != file_name_vector_MB.end(); listIterator_MB++)
	{
	  TFile *testfile_MB = TFile::Open(*listIterator_MB);
	  
	  if(!testfile_MB || testfile_MB->IsZombie() || testfile_MB->TestBit(TFile::kRecovered))
	    {
	      cout << "MB File: " << *listIterator_MB << " failed to open" << endl;
	      delete testfile_MB;
	      testfile_MB = nullptr;
	      continue;
	    }
	  else
	    {
	      cout << "Adding MB file:--- " << *listIterator_MB << "--- to the chains" << endl;
	      
	      hea_tree_MB->Add(*listIterator_MB);
	      trk_tree_MB->Add(*listIterator_MB);
	      if(is_MC) gen_tree_MB->Add(*listIterator_MB);
	      
	      delete testfile_MB; // Close file after use
	      testfile_MB = nullptr;
	    }
	}
      
      hea_tree_MB->AddFriend(trk_tree_MB);
      if(is_MC) hea_tree_MB->AddFriend(gen_tree_MB);

      // calling function to read forest/skim tree
      read_tree_Mixing(hea_tree_MB, is_MC, colliding_system); // access the tree informations
      
      int nevents_MB = hea_tree_MB->GetEntries();
      
      std::cout<<std::endl;
      cout << "Total number of MB events in those files: "<< nevents_MB << endl;
      std::cout<<std::endl;

      int evtcount = 0; // required for mixing
      
      //for(int imb = 0; imb < nevents_MB; imb++) //event loop //start
      //for(int imb = 0; imb < 1000; imb++) //event loop start
      for(int imb = 0; imb < (int)((double)(nevents_MB)/1.7); imb++) //event loop //start
	{
	  hea_tree_MB->GetEntry(imb);
	  
	  if(imb%50000 == 0)
	    {
	      std::cout<<imb<<"  MB events running"<<std::endl;
	    }
	  
	  if(vertexz_MB <= vz_cut_min || vertexz_MB >= vz_cut_max) continue; // apply z vertex cut
	  
	  // determine centrality here
	  int centbin_MB;
	  int ctbin_MB = -1;
	  
	  if(use_cent) // if you use centrality
	    {
	      if(is_MC)
		{
		  if(hiBin_MB <= 9) continue;
		  centbin_MB = hiBin_MB - 10; // match MC multiplicity with data multiplicity
		  //centbin_MB = hiBin_MB;
		}
	      else
		{
		  centbin_MB = hiBin_MB;
		}
	      
	      if(centbin_MB >= centCut || centbin_MB < 0) continue;
	      ctbin_MB = hCentBin->FindBin(centbin_MB) - 1;
	    }
	  else // if you use multiplicity
	    {
	      centbin_MB = -1;
	      if(centbin_MB >= mult_Cut) continue;
	      ctbin_MB = 1;
	    }

	  
	  // Filled event quantities vector                                                                                                    
	  hCent_MB->Fill(centbin_MB, 1);
	  Evtno_vec_MB_1D.push_back(event_MB);
	  Vertexz_vec_MB_1D.push_back(vertexz_MB);
	  HiBin_vec_MB_1D.push_back(ctbin_MB);
	  HiBinValue_vec_MB_1D.push_back(centbin_MB);
	  EvtCount_vec_MB_1D.push_back(evtcount);

	  evtcount++;

	  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~select reco/data tracks from MB events~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~         
	  //Define 1D inclusive trk vector used to fill 2D inclusive trk vector                                                               
	  std::vector<TVector3> Reco_Filtered_InclTrk_pT_vec_MB_1D;           // trk vec                                                         
	  std::vector<double> Reco_Filtered_InclTrkW_vec_MB_1D;               // trk weight                                                      
	  std::vector<int> Reco_Filtered_InclTrkSube_vec_MB_1D;               // trk sube

	  //start reco/data tracks loop                                                                                                            
	  for(int irctrk_MB = 0; irctrk_MB < ntrk_MB; irctrk_MB++)
	    {
	      float trk_pt_MB = (float)trkpt_MB[irctrk_MB];
	      float trk_eta_MB = (float)trketa_MB[irctrk_MB];
	      float trk_phi_MB = (float)trkphi_MB[irctrk_MB];

	      if(trk_pt_MB <= trk_pt_min_cut || trk_pt_MB >= trk_pt_max_cut) continue;
	      if(fabs(trk_eta_MB) >= trk_eta_cut) continue;

	      double trkeff_MB = 1., trkfak_MB = 1.;
	      if(is_MC)
		{
		  //trkeff_MB = htrk_eff_PbPb->GetBinContent(htrk_eff_PbPb->FindBin(trk_eta_MB, trk_pt_MB, centbin_MB+10));
		  //trkfak_MB = htrk_fak_PbPb->GetBinContent(htrk_fak_PbPb->FindBin(trk_eta_MB, trk_pt_MB, centbin_MB+10));
		  trkeff_MB = htrk_eff_PbPb->GetBinContent(htrk_eff_PbPb->FindBin(trk_eta_MB, trk_pt_MB, centbin_MB));
                  trkfak_MB = htrk_fak_PbPb->GetBinContent(htrk_fak_PbPb->FindBin(trk_eta_MB, trk_pt_MB, centbin_MB));
		}
	      else
		{
		  trkeff_MB = htrk_eff_PbPb->GetBinContent(htrk_eff_PbPb->FindBin(trk_eta_MB, trk_pt_MB, centbin_MB));
                  trkfak_MB = htrk_fak_PbPb->GetBinContent(htrk_fak_PbPb->FindBin(trk_eta_MB, trk_pt_MB, centbin_MB));
		}
	      double trk_wt_MB = 1.;
	      if(trkeff_MB > 0.001){trk_wt_MB = (1. - trkfak_MB)/trkeff_MB;}
	      if((trk_pt_MB < 0 || trk_pt_MB > 500) || (fabs(trk_eta_MB) > 2.4)) trk_wt_MB = 0.; // check bound

	      TVector3 Reco_Filtered_InclTrk_pT_MB;
	      Reco_Filtered_InclTrk_pT_MB.SetPtEtaPhi(trk_pt_MB, trk_eta_MB, trk_phi_MB);
	      
	      Reco_Filtered_InclTrk_pT_vec_MB_1D.push_back(Reco_Filtered_InclTrk_pT_MB);
	      Reco_Filtered_InclTrkW_vec_MB_1D.push_back(trk_wt_MB);
	      Reco_Filtered_InclTrkSube_vec_MB_1D.push_back(1);
	    } // track loop end
	  
	  // filled MB reco/data 2D trk vector                                                                                                        
	  Reco_Filtered_InclTrk_pT_vec_MB_2D.push_back(Reco_Filtered_InclTrk_pT_vec_MB_1D);
	  Reco_Filtered_InclTrkW_vec_MB_2D.push_back(Reco_Filtered_InclTrkW_vec_MB_1D);
	  Reco_Filtered_InclTrkSube_vec_MB_2D.push_back(Reco_Filtered_InclTrkSube_vec_MB_1D);
	  
	  // clear MB reco/data 1D trk vector used to fill 2D vector                                                                                  
	  Reco_Filtered_InclTrk_pT_vec_MB_1D.clear();
	  Reco_Filtered_InclTrkW_vec_MB_1D.clear();
	  Reco_Filtered_InclTrkSube_vec_MB_1D.clear();
	  
	  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~select MB gen tracks~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~          
	  // for reco jets and gen tracks                                                                                                          
	  std::vector<TVector3> Gen_Filtered_InclTrk_pT_vec_MB_1D;           // trk vec                                                           
	  std::vector<double> Gen_Filtered_InclTrkW_vec_MB_1D;               // trk weight                                                        
	  std::vector<int> Gen_Filtered_InclTrkSube_vec_MB_1D;               // trk sube                                                          
	  
	  if(is_MC)
	    {
	      //start MB gen track loop
	      for(int igntrk_MB = 0; igntrk_MB < gen_trkpt_MB->size(); igntrk_MB++)
		{
		  float gentrk_pt_MB = float((*gen_trkpt_MB)[igntrk_MB]);
		  float gentrk_eta_MB = float((*gen_trketa_MB)[igntrk_MB]);
		  float gentrk_phi_MB = float((*gen_trkphi_MB)[igntrk_MB]);
		  int gentrk_sube_MB = int((*gen_trksube_MB)[igntrk_MB]);

		  if(gentrk_pt_MB <= trk_pt_min_cut || gentrk_pt_MB >= trk_pt_max_cut) continue;
		  if(fabs(gentrk_eta_MB) >= trk_eta_cut) continue;
		  
		  double gentrk_wt_MB = 1.;

		  TVector3 Gen_Filtered_InclTrk_pT_MB;
		  Gen_Filtered_InclTrk_pT_MB.SetPtEtaPhi(gentrk_pt_MB, gentrk_eta_MB, gentrk_phi_MB);
		  
		  Gen_Filtered_InclTrk_pT_vec_MB_1D.push_back(Gen_Filtered_InclTrk_pT_MB);
		  Gen_Filtered_InclTrkW_vec_MB_1D.push_back(gentrk_wt_MB);
		  Gen_Filtered_InclTrkSube_vec_MB_1D.push_back(gentrk_sube_MB);
		} // gen trk loop end
	      // filled gen 2D trk vector                                                                                                          
	      Gen_Filtered_InclTrk_pT_vec_MB_2D.push_back(Gen_Filtered_InclTrk_pT_vec_MB_1D);
	      Gen_Filtered_InclTrkW_vec_MB_2D.push_back(Gen_Filtered_InclTrkW_vec_MB_1D);
	      Gen_Filtered_InclTrkSube_vec_MB_2D.push_back(Gen_Filtered_InclTrkSube_vec_MB_1D);
	      // clear gen 1D trk vector used to fill 2D vector                                                                                    
	      Gen_Filtered_InclTrk_pT_vec_MB_1D.clear();
	      Gen_Filtered_InclTrkW_vec_MB_1D.clear();
	      Gen_Filtered_InclTrkSube_vec_MB_1D.clear();
	    }// if(is_MC)
	}// MB event loop end
    }// if (colliding_system && isMBEventsMixing)
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~calling function for jet tracks signal and mixing correlation~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(!isRcJetGnTrk) // for gen jets gen tracks and reco jets reco tracks
    {
      std::cout<<"Correlation will be done between reco jets and reco tracks and vice varsa"<<std::endl;

      // signal
      // for reco/data
      Jet_Track_signal_corr_ldsld(colliding_system, Reco_Evtw_vec_1D, HiBin_vec_1D, HiBinValue_vec_1D, Vertexz_vec_1D, Reco_Filtered_ldJet_CorrpT_vec_1D, Reco_Filtered_ldrefpartonB_vec_1D, Reco_Filtered_ldJetW_vec_1D, Reco_Filtered_sldJet_CorrpT_vec_1D, Reco_Filtered_sldrefpartonB_vec_1D, Reco_Filtered_sldJetW_vec_1D, Reco_Filtered_InclTrk_pT_vec_2D, Reco_Filtered_InclTrkW_vec_2D, Reco_Filtered_InclTrkCharge_vec_2D, Reco_Filtered_InclTrkSube_vec_2D, true, false); // true for reco, false for sube

      // mixing
      // for reco/data
      if(colliding_system == "PbPb")
	{
	  if(!isMBEventsMixing)
	    {
	      std::cout<<"Mixed event is done with dijet data"<<std::endl;
	      std::cout<<endl;
	      Jet_Track_mixing_corr_ldsld(colliding_system, Reco_Evtw_vec_1D, HiBin_vec_1D, HiBinValue_vec_1D, Vertexz_vec_1D, Reco_Evtno_vec_1D, Reco_EvtCount_vec_1D, Reco_Filtered_ldJet_CorrpT_vec_1D, Reco_Filtered_ldrefpartonB_vec_1D, Reco_Filtered_ldJetW_vec_1D, Reco_Filtered_sldJet_CorrpT_vec_1D, Reco_Filtered_sldrefpartonB_vec_1D, Reco_Filtered_sldJetW_vec_1D, Reco_Filtered_InclTrk_pT_vec_2D, Reco_Filtered_InclTrkW_vec_2D, Reco_Filtered_InclTrkCharge_vec_2D, Reco_Filtered_InclTrkSube_vec_2D, true, false); // true for reco, false for sube
	    }
	  else if(isMBEventsMixing)
	    {
	      std::cout<<"Mixed event is done with MB sample"<<std::endl;
	      std::cout<<endl;
	      Jet_Track_mixing_withMB_corr_ldsld(colliding_system, Reco_Evtw_vec_1D, HiBin_vec_1D, HiBinValue_vec_1D, Vertexz_vec_1D, Reco_Evtno_vec_1D, Reco_Filtered_ldJet_CorrpT_vec_1D, Reco_Filtered_ldrefpartonB_vec_1D, Reco_Filtered_ldJetW_vec_1D, Reco_Filtered_sldJet_CorrpT_vec_1D, Reco_Filtered_sldrefpartonB_vec_1D, Reco_Filtered_sldJetW_vec_1D, Reco_Filtered_InclTrk_pT_vec_2D, Reco_Filtered_InclTrkW_vec_2D, Reco_Filtered_InclTrkCharge_vec_2D, Reco_Filtered_InclTrkSube_vec_2D, HiBin_vec_MB_1D, HiBinValue_vec_MB_1D, Vertexz_vec_MB_1D, Evtno_vec_MB_1D, EvtCount_vec_MB_1D, Reco_Filtered_InclTrk_pT_vec_MB_2D, Reco_Filtered_InclTrkW_vec_MB_2D, Reco_Filtered_InclTrkSube_vec_MB_2D, true, false); // true for reco, false for sube, 2nd false always
	    }
	}
      else if(colliding_system == "pp")
	{
	  std::cout<<"Mixed event is done with diejt data"<<std::endl;
	  std::cout<<endl;
	  Jet_Track_mixing_corr_ldsld(colliding_system, Reco_Evtw_vec_1D, HiBin_vec_1D, HiBinValue_vec_1D, Vertexz_vec_1D, Reco_Evtno_vec_1D, Reco_EvtCount_vec_1D, Reco_Filtered_ldJet_CorrpT_vec_1D, Reco_Filtered_ldrefpartonB_vec_1D, Reco_Filtered_ldJetW_vec_1D, Reco_Filtered_sldJet_CorrpT_vec_1D, Reco_Filtered_sldrefpartonB_vec_1D, Reco_Filtered_sldJetW_vec_1D, Reco_Filtered_InclTrk_pT_vec_2D, Reco_Filtered_InclTrkW_vec_2D, Reco_Filtered_InclTrkCharge_vec_2D, Reco_Filtered_InclTrkSube_vec_2D, true, false); // true for reco, false for sube
	}
      
      // for gen
      if(is_MC)
	{
	  // signal
	  if(colliding_system == "PbPb" && !isMBEventsMixing) // for sube 0
	    //if(colliding_system == "PbPb") // for sube 0
	    {
	      Jet_Track_signal_corr_ldsld(colliding_system, Gen_Evtw_vec_1D, Gen_HiBin_vec_1D, Gen_HiBinValue_vec_1D, Gen_Vertexz_vec_1D, Gen_Filtered_ldJet_CorrpT_vec_1D, Gen_Filtered_ldrefpartonB_vec_1D, Gen_Filtered_ldJetW_vec_1D, Gen_Filtered_sldJet_CorrpT_vec_1D, Gen_Filtered_sldrefpartonB_vec_1D, Gen_Filtered_sldJetW_vec_1D, Gen_Filtered_InclTrk_pT_vec_2D, Gen_Filtered_InclTrkW_vec_2D, Gen_Filtered_InclTrkCharge_vec_2D, Gen_Filtered_InclTrkSube_vec_2D, false, true); // false for gen, true for sube
	    }
	  else 
	    {
	      Jet_Track_signal_corr_ldsld(colliding_system, Gen_Evtw_vec_1D, Gen_HiBin_vec_1D, Gen_HiBinValue_vec_1D, Gen_Vertexz_vec_1D, Gen_Filtered_ldJet_CorrpT_vec_1D, Gen_Filtered_ldrefpartonB_vec_1D, Gen_Filtered_ldJetW_vec_1D, Gen_Filtered_sldJet_CorrpT_vec_1D, Gen_Filtered_sldrefpartonB_vec_1D, Gen_Filtered_sldJetW_vec_1D, Gen_Filtered_InclTrk_pT_vec_2D, Gen_Filtered_InclTrkW_vec_2D, Gen_Filtered_InclTrkCharge_vec_2D, Gen_Filtered_InclTrkSube_vec_2D, false, false); // false for gen, true for sube
	    }

	  // mixing
	  if(colliding_system == "PbPb")
	    {
	      if(!isMBEventsMixing)
		{
		  std::cout<<"Mixed event is done with diejt data"<<std::endl;
		  std::cout<<endl;
		  Jet_Track_mixing_corr_ldsld(colliding_system, Gen_Evtw_vec_1D, Gen_HiBin_vec_1D, Gen_HiBinValue_vec_1D, Gen_Vertexz_vec_1D, Gen_Evtno_vec_1D, Gen_EvtCount_vec_1D, Gen_Filtered_ldJet_CorrpT_vec_1D, Gen_Filtered_ldrefpartonB_vec_1D, Gen_Filtered_ldJetW_vec_1D, Gen_Filtered_sldJet_CorrpT_vec_1D, Gen_Filtered_sldrefpartonB_vec_1D, Gen_Filtered_sldJetW_vec_1D, Gen_Filtered_InclTrk_pT_vec_2D, Gen_Filtered_InclTrkW_vec_2D, Gen_Filtered_InclTrkCharge_vec_2D, Gen_Filtered_InclTrkSube_vec_2D, false, true); // false for gen, true for sube
		}
	      else if(isMBEventsMixing)
		{
		  std::cout<<"Mixed event is done with MB sample"<<std::endl;
		  std::cout<<endl;
		  Jet_Track_mixing_withMB_corr_ldsld(colliding_system, Gen_Evtw_vec_1D, Gen_HiBin_vec_1D, Gen_HiBinValue_vec_1D, Gen_Vertexz_vec_1D, Gen_Evtno_vec_1D, Gen_Filtered_ldJet_CorrpT_vec_1D, Gen_Filtered_ldrefpartonB_vec_1D, Gen_Filtered_ldJetW_vec_1D, Gen_Filtered_sldJet_CorrpT_vec_1D, Gen_Filtered_sldrefpartonB_vec_1D, Gen_Filtered_sldJetW_vec_1D, Gen_Filtered_InclTrk_pT_vec_2D, Gen_Filtered_InclTrkW_vec_2D, Gen_Filtered_InclTrkCharge_vec_2D, Gen_Filtered_InclTrkSube_vec_2D, HiBin_vec_MB_1D, HiBinValue_vec_MB_1D, Vertexz_vec_MB_1D, Evtno_vec_MB_1D, EvtCount_vec_MB_1D, Gen_Filtered_InclTrk_pT_vec_MB_2D, Gen_Filtered_InclTrkW_vec_MB_2D, Gen_Filtered_InclTrkSube_vec_MB_2D, false, false); // false for gen, false for sube, 2nd false always
		}
	    }
	  else if(colliding_system == "pp")
	    {
	      std::cout<<"Mixed event is done with diejt data"<<std::endl;
	      std::cout<<endl;
	      Jet_Track_mixing_corr_ldsld(colliding_system, Gen_Evtw_vec_1D, Gen_HiBin_vec_1D, Gen_HiBinValue_vec_1D, Gen_Vertexz_vec_1D, Gen_Evtno_vec_1D, Gen_EvtCount_vec_1D, Gen_Filtered_ldJet_CorrpT_vec_1D, Gen_Filtered_ldrefpartonB_vec_1D, Gen_Filtered_ldJetW_vec_1D, Gen_Filtered_sldJet_CorrpT_vec_1D, Gen_Filtered_sldrefpartonB_vec_1D, Gen_Filtered_sldJetW_vec_1D, Gen_Filtered_InclTrk_pT_vec_2D, Gen_Filtered_InclTrkW_vec_2D, Gen_Filtered_InclTrkCharge_vec_2D, Gen_Filtered_InclTrkSube_vec_2D, false, false); // false for gen, true for sube
	    }
	} // is MC
    }
  else if(isRcJetGnTrk) // for gen jets reco tracks and reco jets gen tracks
    {
      std::cout<<"Correlation will be done between reco jets and gen tracks and vice varsa"<<std::endl;
      
      // signal
      // for reco/data
      Jet_Track_signal_corr_ldsld(colliding_system, Reco_Evtw_vec_1D, HiBin_vec_1D, HiBinValue_vec_1D, Vertexz_vec_1D, Reco_Filtered_ldJet_CorrpT_vec_1D, Reco_Filtered_ldrefpartonB_vec_1D, Reco_Filtered_ldJetW_vec_1D, Reco_Filtered_sldJet_CorrpT_vec_1D, Reco_Filtered_sldrefpartonB_vec_1D, Reco_Filtered_sldJetW_vec_1D, RecoGen_Filtered_InclTrk_pT_vec_2D, RecoGen_Filtered_InclTrkW_vec_2D, RecoGen_Filtered_InclTrkCharge_vec_2D, RecoGen_Filtered_InclTrkSube_vec_2D, true, false); // true for reco, false for sube
      
      // mixing
      // for reco/data
      if(colliding_system == "PbPb")
        {
          if(!isMBEventsMixing)
            {
	      std::cout<<"Mixed event is done dijet data"<<std::endl;
	      std::cout<<endl;
	      Jet_Track_mixing_corr_ldsld(colliding_system, Reco_Evtw_vec_1D, HiBin_vec_1D, HiBinValue_vec_1D, Vertexz_vec_1D, Reco_Evtno_vec_1D, Reco_EvtCount_vec_1D, Reco_Filtered_ldJet_CorrpT_vec_1D, Reco_Filtered_ldrefpartonB_vec_1D, Reco_Filtered_ldJetW_vec_1D, Reco_Filtered_sldJet_CorrpT_vec_1D, Reco_Filtered_sldrefpartonB_vec_1D, Reco_Filtered_sldJetW_vec_1D, RecoGen_Filtered_InclTrk_pT_vec_2D, RecoGen_Filtered_InclTrkW_vec_2D, RecoGen_Filtered_InclTrkCharge_vec_2D, RecoGen_Filtered_InclTrkSube_vec_2D, true, false); // true for reco, false for sube
	    }
	  else if(isMBEventsMixing)
	    {
	      std::cout<<"Mixed event is done with MB sample"<<std::endl;
	      std::cout<<endl;
	      Jet_Track_mixing_withMB_corr_ldsld(colliding_system, Reco_Evtw_vec_1D, HiBin_vec_1D, HiBinValue_vec_1D, Vertexz_vec_1D, Reco_Evtno_vec_1D, Reco_Filtered_ldJet_CorrpT_vec_1D, Reco_Filtered_ldrefpartonB_vec_1D, Reco_Filtered_ldJetW_vec_1D, Reco_Filtered_sldJet_CorrpT_vec_1D, Reco_Filtered_sldrefpartonB_vec_1D, Reco_Filtered_sldJetW_vec_1D, RecoGen_Filtered_InclTrk_pT_vec_2D, RecoGen_Filtered_InclTrkW_vec_2D, RecoGen_Filtered_InclTrkCharge_vec_2D, RecoGen_Filtered_InclTrkSube_vec_2D, HiBin_vec_MB_1D, HiBinValue_vec_MB_1D, Vertexz_vec_MB_1D, Evtno_vec_MB_1D, EvtCount_vec_MB_1D, Gen_Filtered_InclTrk_pT_vec_MB_2D, Gen_Filtered_InclTrkW_vec_MB_2D, Gen_Filtered_InclTrkSube_vec_MB_2D, true, false); // true for reco, false for sube, 2nd false always
	    }
	}
      else if(colliding_system == "pp")
	{
	  std::cout<<"Mixed event is done dijet data"<<std::endl;
	  std::cout<<endl;
	  Jet_Track_mixing_corr_ldsld(colliding_system, Reco_Evtw_vec_1D, HiBin_vec_1D, HiBinValue_vec_1D, Vertexz_vec_1D, Reco_Evtno_vec_1D, Reco_EvtCount_vec_1D, Reco_Filtered_ldJet_CorrpT_vec_1D, Reco_Filtered_ldrefpartonB_vec_1D, Reco_Filtered_ldJetW_vec_1D, Reco_Filtered_sldJet_CorrpT_vec_1D, Reco_Filtered_sldrefpartonB_vec_1D, Reco_Filtered_sldJetW_vec_1D, RecoGen_Filtered_InclTrk_pT_vec_2D, RecoGen_Filtered_InclTrkW_vec_2D, RecoGen_Filtered_InclTrkCharge_vec_2D, RecoGen_Filtered_InclTrkSube_vec_2D, true, false); // true for reco, false for sube
	}
      
      // for gen
      if(is_MC)
	{
	  // signal
	  Jet_Track_signal_corr_ldsld(colliding_system, Gen_Evtw_vec_1D, Gen_HiBin_vec_1D, Gen_HiBinValue_vec_1D, Gen_Vertexz_vec_1D, Gen_Filtered_ldJet_CorrpT_vec_1D, Gen_Filtered_ldrefpartonB_vec_1D, Gen_Filtered_ldJetW_vec_1D, Gen_Filtered_sldJet_CorrpT_vec_1D, Gen_Filtered_sldrefpartonB_vec_1D, Gen_Filtered_sldJetW_vec_1D, GenReco_Filtered_InclTrk_pT_vec_2D, GenReco_Filtered_InclTrkW_vec_2D, GenReco_Filtered_InclTrkCharge_vec_2D, GenReco_Filtered_InclTrkSube_vec_2D, false, false); // false for reco, false for sube
	  
	  // mixing
	  if(colliding_system == "PbPb")
            {
              if(!isMBEventsMixing)
                {
                  std::cout<<"Mixed event is done with diejt data"<<std::endl;
                  std::cout<<endl;
		  Jet_Track_mixing_corr_ldsld(colliding_system, Gen_Evtw_vec_1D, Gen_HiBin_vec_1D, Gen_HiBinValue_vec_1D, Gen_Vertexz_vec_1D, Gen_Evtno_vec_1D, Gen_EvtCount_vec_1D, Gen_Filtered_ldJet_CorrpT_vec_1D, Gen_Filtered_ldrefpartonB_vec_1D, Gen_Filtered_ldJetW_vec_1D, Gen_Filtered_sldJet_CorrpT_vec_1D, Gen_Filtered_sldrefpartonB_vec_1D, Gen_Filtered_sldJetW_vec_1D, GenReco_Filtered_InclTrk_pT_vec_2D, GenReco_Filtered_InclTrkW_vec_2D, GenReco_Filtered_InclTrkCharge_vec_2D, GenReco_Filtered_InclTrkSube_vec_2D, false, false); // false for gen, false for sube
		}
	      else if(isMBEventsMixing)
		{
		  std::cout<<"Mixed event is done with MB sample"<<std::endl;
		  std::cout<<endl;
		  Jet_Track_mixing_withMB_corr_ldsld(colliding_system, Gen_Evtw_vec_1D, Gen_HiBin_vec_1D, Gen_HiBinValue_vec_1D, Gen_Vertexz_vec_1D, Gen_Evtno_vec_1D, Gen_Filtered_ldJet_CorrpT_vec_1D, Gen_Filtered_ldrefpartonB_vec_1D, Gen_Filtered_ldJetW_vec_1D, Gen_Filtered_sldJet_CorrpT_vec_1D, Gen_Filtered_sldrefpartonB_vec_1D, Gen_Filtered_sldJetW_vec_1D, GenReco_Filtered_InclTrk_pT_vec_2D, GenReco_Filtered_InclTrkW_vec_2D, GenReco_Filtered_InclTrkCharge_vec_2D, GenReco_Filtered_InclTrkSube_vec_2D, HiBin_vec_MB_1D, HiBinValue_vec_MB_1D, Vertexz_vec_MB_1D, Evtno_vec_MB_1D, EvtCount_vec_MB_1D, Reco_Filtered_InclTrk_pT_vec_MB_2D, Reco_Filtered_InclTrkW_vec_MB_2D, Reco_Filtered_InclTrkSube_vec_MB_2D, false, false); // false for gen, false for sube, 2nd false always
		}
	    }
	  else if(colliding_system == "pp")
	    {
	      std::cout<<"Mixed event is done with diejt data"<<std::endl;
	      std::cout<<endl;
	      Jet_Track_mixing_corr_ldsld(colliding_system, Gen_Evtw_vec_1D, Gen_HiBin_vec_1D, Gen_HiBinValue_vec_1D, Gen_Vertexz_vec_1D, Gen_Evtno_vec_1D, Gen_EvtCount_vec_1D, Gen_Filtered_ldJet_CorrpT_vec_1D, Gen_Filtered_ldrefpartonB_vec_1D, Gen_Filtered_ldJetW_vec_1D, Gen_Filtered_sldJet_CorrpT_vec_1D, Gen_Filtered_sldrefpartonB_vec_1D, Gen_Filtered_sldJetW_vec_1D, GenReco_Filtered_InclTrk_pT_vec_2D, GenReco_Filtered_InclTrkW_vec_2D, GenReco_Filtered_InclTrkCharge_vec_2D, GenReco_Filtered_InclTrkSube_vec_2D, false, false); // false for gen, false for sube
	    }
	} // is_MC
    } // isRcJetGnTrk
  
  // clear 1D event vector
  //for gen
  Gen_Evtw_vec_1D.clear();
  Gen_Evtno_vec_1D.clear();
  Gen_EvtCount_vec_1D.clear();
  Gen_Vertexz_vec_1D.clear();
  Gen_HiBin_vec_1D.clear();
  Gen_HiBinValue_vec_1D.clear();
  
  //for reco
  Reco_Evtw_vec_1D.clear();
  Reco_Evtno_vec_1D.clear();
  Reco_EvtCount_vec_1D.clear();
  Vertexz_vec_1D.clear();
  HiBin_vec_1D.clear();
  HiBinValue_vec_1D.clear();

  // clear 1D leading jet vector
  //for gen
  Gen_Filtered_ldJet_CorrpT_vec_1D.clear();
  Gen_Filtered_ldrefpartonB_vec_1D.clear();
  Gen_Filtered_ldJetW_vec_1D.clear();
    
  //for reco or data
  Reco_Filtered_ldJet_CorrpT_vec_1D.clear();
  Reco_Filtered_ldrefpartonB_vec_1D.clear();
  Reco_Filtered_ldJetW_vec_1D.clear();
    
  //clear 1D subleading jet vectors
  //for gen
  Gen_Filtered_sldJet_CorrpT_vec_1D.clear();
  Gen_Filtered_sldrefpartonB_vec_1D.clear();
  Gen_Filtered_sldJetW_vec_1D.clear();

  //for reco or data
  Reco_Filtered_sldJet_CorrpT_vec_1D.clear();
  Reco_Filtered_sldrefpartonB_vec_1D.clear();
  Reco_Filtered_ldJetW_vec_1D.clear();

  // clear 2D trk vector
  // for gen
  Gen_Filtered_InclTrk_pT_vec_2D.clear();
  Gen_Filtered_InclTrkW_vec_2D.clear();
  Gen_Filtered_InclTrkCharge_vec_2D.clear();
  Gen_Filtered_InclTrkSube_vec_2D.clear();
  // for gen jets and reco tracks
  GenReco_Filtered_InclTrk_pT_vec_2D.clear();
  GenReco_Filtered_InclTrkW_vec_2D.clear();
  GenReco_Filtered_InclTrkCharge_vec_2D.clear();
  GenReco_Filtered_InclTrkSube_vec_2D.clear();
  
  // for reco
  Reco_Filtered_InclTrk_pT_vec_2D.clear();
  Reco_Filtered_InclTrkW_vec_2D.clear();
  Reco_Filtered_InclTrkCharge_vec_2D.clear();
  Reco_Filtered_InclTrkSube_vec_2D.clear();
  // for reco jets and gen tracks
  RecoGen_Filtered_InclTrk_pT_vec_2D.clear();
  RecoGen_Filtered_InclTrkW_vec_2D.clear();
  RecoGen_Filtered_InclTrkCharge_vec_2D.clear();
  RecoGen_Filtered_InclTrkSube_vec_2D.clear();

  // close / delete all the tree, files, histograms, TF1 that are open
  delete hea_tree;
  delete hlt_tree;
  delete ski_tree;
  delete jet_tree;
  delete trk_tree;
  delete gen_tree;

  // clear 1D MB event vector
  Evtno_vec_MB_1D.clear();
  Vertexz_vec_MB_1D.clear();
  HiBin_vec_MB_1D.clear();
  HiBinValue_vec_MB_1D.clear();
  EvtCount_vec_MB_1D.clear();
  
  // clear 2D MB gen trk vector
  Gen_Filtered_InclTrk_pT_vec_MB_2D.clear();
  Gen_Filtered_InclTrkW_vec_MB_2D.clear();
  Gen_Filtered_InclTrkSube_vec_MB_2D.clear();

  // clear 2D MB gen trk vector
  Reco_Filtered_InclTrk_pT_vec_MB_2D.clear();
  Reco_Filtered_InclTrkW_vec_MB_2D.clear();
  Reco_Filtered_InclTrkSube_vec_MB_2D.clear();

  delete hea_tree_MB;
  delete trk_tree_MB;
  delete gen_tree_MB;
  
  if(is_MC && is_CentBin_and_VZ_Weight)
    {
      hCentVzFile->Close();
      delete hCentVzFile;
    }

  if(is_MC && colliding_system == "pp" && is_ptWeight)
    {
      fptFile_pp->Close();
      delete fptFile_pp;
    }

  if(is_MC && colliding_system == "PbPb" && is_ptWeight)
    {
      for(int ictt = 0; ictt < NCentbin-1; ictt++)
	{
	  fptFile_PbPb[ictt]->Close();
	  delete fptFile_PbPb[ictt];
	}
    }
  
  if(is_MC && colliding_system == "pp" && is_JER)
    {
      fJERFile_pp->Close();                                                                                                                 
      delete fJERFile_pp;
    }

  if(is_MC && colliding_system == "PbPb" && is_JER)
    {
      for(int ictt = 0; ictt < NCentbin-1; ictt++)
	{
	  fJERFile_PbPb[ictt]->Close();
	  delete fJERFile_PbPb[ictt];
	}
    }

  if(is_MC && is_JEU)
    {
      /*
      delete JEU_ld;
      JEU_ld = nullptr;
      
      delete JEU_sld;
      JEU_sld = nullptr;
      */
      delete JEU;
      JEU = nullptr;
    }
  if(is_MC && doEventFilter && colliding_system == "PbPb")
    {
      xjetFrac->Close();
      delete xjetFrac;
    }

  if(colliding_system == "pp")
    {
      delete htrk_eff_pp;
      delete htrk_fak_pp;
      delete htrk_sec_pp;
    }
  else if(colliding_system == "PbPb")
    {
      delete htrk_eff_PbPb;
      delete htrk_fak_PbPb;
    }
  
  ftrk_eff->Close();
  delete ftrk_eff;
  if(colliding_system == "PbPb")
    {
      ftrk_fak->Close();
      delete ftrk_fak;
    }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~define output file
  std::string outfilename = Form("%s/%s_Outfile_pTHat%1.1f_JetpT%1.1f_LdpT%1.1f_SldpT%1.1f_JetEta%1.1f_%d",out_file.Data(), colliding_system_filename.Data(), pthat_cut, jet_pt_min_cut, leading_pT_min_cut, subleading_pT_min_cut, jet_eta_max_cut, itxtoutFile);
  
  std::replace(outfilename.begin(), outfilename.end(), '.', 'p'); // replace . to p
  std::replace(outfilename.begin(), outfilename.end(), '-', 'N'); // replace - to N for negative

  TFile* fout = new TFile(Form("%s.root", outfilename.c_str()), "recreate");	 

  //TFile* fout = new TFile(Form("%s/%s_Outfile_NoRefpTCut_pTHat%1.1f_JetpT%1.1f_LdpT%1.1f_SldpT%1.1f_NoTrigger_%d.root",out_file.Data(), colliding_system.Data(), pthat_cut, jet_pt_min_cut, leading_pT_min, subleading_pT_min, itxtoutFile), "recreate");

  fout->mkdir("Event_Hist");
  fout->cd("Event_Hist");
  Write_Event_hist(is_MC);

  /*
  fout->mkdir("Jet_QA_Hist");
  fout->cd("Jet_QA_Hist");
  Write_Jet_QA_hist(is_MC, is_JES_JER);
  */

  fout->mkdir("Trk_QA_Hist");
  fout->cd("Trk_QA_Hist");
  Write_Trk_QA_hist(is_MC);
  
  fout->mkdir("Jet_Trk_Corr_Hist");
  fout->cd("Jet_Trk_Corr_Hist");
  Write_Jet_Trk_Corr_hist(is_MC);
    
  fout->Write();
  fout->Close();
  delete fout;
} // void Tree_Analyzer() end

// main program
int main(int argc, char **argv)
{
  clock_t sec_start, sec_end;
  sec_start = clock();

  //TDatime* date = new TDatime();

  print_start(); // start timing print

  using namespace std;

  TString inputfile;
  int itxtout;
  TString outfile;
  TString coll_sys;
  int ismc;

  if(argc == 1)
    {
      std::cout<<"You did not pass any argument to the code other than the program name"<<std::endl;
      std::cout<<"You need to pass 6 arguments including the program name"<<std::endl;
    }

  if(argc >=1 && argc <=5)
    {
      std::cout<<"Only "<<argc<<" arguments you have given including the program name"<<std::endl;
      std::cout<<"You need to pass 6 arguments including the program name"<<std::endl;
    }

  if(argc == 6)
    {
      std::cout<<std::endl;
      std::cout<<"You have given "<< argc <<" arguments including the program name;  Your program will run"<<std::endl;
      std::cout<<std::endl;

      inputfile = argv[1];
      itxtout = atoi(argv[2]);
      outfile = argv[3];
      coll_sys = argv[4];
      ismc = atoi(argv[5]);

      Tree_Analyzer(inputfile, itxtout, outfile, coll_sys, ismc);
    }

  sec_end = clock(); // stop time counting                                                                                                            
  cout << "========================================" << endl;
  cout << "Total running time: " << (double)(sec_end - sec_start) / CLOCKS_PER_SEC << " [s]" << endl;
  cout << "========================================" << endl;

  print_stop(); // Print time, date and hour when it stops
  
  return 0;
}
