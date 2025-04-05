//#include "call_libraries.h"

bool pTHatFilter(double jetPt, double pthat)
{
  bool result = false;
  
  if(jetPt < 1.702*pthat + 9.701) result = true;
  
  return result;
}

void find_leading_subleading_jets(double pt, double refpt, int index, double &leadpt, double &leadrefpt, int &leadindex, double &sublpt, double &sublrefpt, int &sublindex)
{
  if( pt > leadpt )
    {
      sublpt = leadpt;
      sublindex = leadindex;
      sublrefpt = leadrefpt;
      
      leadpt = pt;
      leadindex = index;
      leadrefpt = refpt;
    }
  else if( pt > sublpt)
    {
      sublpt = pt;
      sublindex = index;
      sublrefpt = refpt;
    }
}

void find_leading_subleading_gen_jets(double pt, int index, int ref_index, double &leadpt, int &leadindex, int &ref_leadindex, double &sublpt, int &sublindex, int &ref_subleadindex)
{
  if( pt > leadpt )
    {
      sublpt = leadpt;
      sublindex = leadindex;
      ref_subleadindex = ref_leadindex;
      
      leadpt = pt;
      leadindex = index;
      ref_leadindex = ref_index;
    }
  else if( pt > sublpt)
    {
      sublpt = pt;
      sublindex = index;
      ref_subleadindex = ref_index;
    }
}

void find_leading_subleading_tracks(double pt, int index, double &leadpt, int &leadindex, double &sublpt, int &sublindex)
{
  if( pt > leadpt )
    {
      sublpt = leadpt;
      sublindex = leadindex;
      
      leadpt = pt;
      leadindex = index;
    }
  else if( pt > sublpt)
    {
      sublpt = pt;
      sublindex = index;
    }
}

void find_leading_jets(double pt, double eta, double phi, int index, double refpt, int falvorB, double &leadpt, double &leadeta, double &leadphi, int &leadindex, double &leadrefpt, int &leadfalvorB)
{
  if( pt > leadpt )
    {
      leadpt = pt;
      leadeta = eta;
      leadphi = phi;
      leadindex = index;
      leadrefpt = refpt;
      leadfalvorB = falvorB;
    }
}

void Jet_Track_signal_corr_ldsld(const TString& colliding_system, const std::vector<double>& Evtw_vec_1D, const std::vector<int>& HiBin_vec_1D, const std::vector<int>& HiBinValue_vec_1D, const std::vector<double>& Vertexz_vec_1D, const std::vector<TVector3>& Filtered_ldjet_vec_1D, const std::vector<int>& Filtered_ldrefpartonB_vec_1D, const std::vector<double>& Filtered_ldJetW_vec_1D, const std::vector<TVector3>& Filtered_sldjet_vec_1D, const std::vector<int>& Filtered_sldrefpartonB_vec_1D, const std::vector<double>& Filtered_sldJetW_vec_1D, const std::vector<std::vector<TVector3>>& Filtered_Trk_pT_vec_2D, const std::vector<std::vector<double>>& Filtered_TrkW_vec_2D, const std::vector<std::vector<int>>& Filtered_TrkCharge_vec_2D, const std::vector<std::vector<int>>& Filtered_TrkSube_vec_2D, const bool& isrc, const bool& do_sube, const bool& do_sube_rcjetgntrk)
{
  std::cout<<endl;
  if(Evtw_vec_1D.size() != Filtered_ldjet_vec_1D.size()){std::cout<<"event numbers are not same, pleaes check"<<std::endl;}
  else {std::cout<<"Total events is: "<<Evtw_vec_1D.size()<<"  "<<Filtered_ldjet_vec_1D.size()<<"  "<<Filtered_Trk_pT_vec_2D.size()<<std::endl;}

  std::cout<<endl;

  if(isrc) {std::cout<<"~~~~~start leading subleading jet tracks signal correlation~~~~~~~~"<<std::endl;}
  else {std::cout<<"~~~~~start gen leading subleading jet tracks signal correlation~~~~~~~~"<<std::endl;}

  int TotaldijetEvent = 0;
  
  for(int ievt = 0; ievt < Vertexz_vec_1D.size(); ievt++)
    {
      TotaldijetEvent++;
	
      double evtw = Evtw_vec_1D[ievt];
      int ctbin = HiBin_vec_1D[ievt];
      double vz = Vertexz_vec_1D[ievt];

      TVector3 ldjet_vec = Filtered_ldjet_vec_1D[ievt];
      double ldJetW = Filtered_ldJetW_vec_1D[ievt];
      int ldJetflvor = Filtered_ldrefpartonB_vec_1D[ievt];
      double ldJet_pt = ldjet_vec.Pt();
      double ldJet_eta = ldjet_vec.Eta();
      double ldJet_phi = ldjet_vec.Phi();
      
      TVector3 sldjet_vec = Filtered_sldjet_vec_1D[ievt];
      double sldJetW = Filtered_sldJetW_vec_1D[ievt];
      int sldJetflvor = Filtered_sldrefpartonB_vec_1D[ievt];
      double sldJet_pt = sldjet_vec.Pt();
      double sldJet_eta = sldjet_vec.Eta();
      double sldJet_phi = sldjet_vec.Phi();

      double Deta_ldsldJet = fabs(ldJet_eta - sldJet_eta);

      	    
      if(ldJet_pt < leading_pT_min_cut || sldJet_pt < subleading_pT_min_cut || fabs(ldJet_eta) > jet_eta_max_cut || fabs(sldJet_eta) > jet_eta_max_cut) {std::cout<<"Leading and subleading jets are out of pT and eta cuts, please check"<<std::endl;}

      if(fabs(ldJet_eta) > 0.5) continue; // for test
      
      int ldjtTRkCorr = 0, ldjtTRkCorr_RapAsym = 0, ldjtTRkCorr_RapAsym_1 = 0, ldjtTRkCorr_RapAsym_2 = 0, ldjtTRkCorr_RapAsym_3 = 0, ldjtTRkCorr_RapAsym_sh = 0, ldjtTRkCorr_RapAsym_oh = 0; // for reco jet reco trks
      int rcJet_gnTrk_sube0_ldjtTRkCorr = 0, rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym = 0, rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_1 = 0, rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_2 = 0, rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_3 = 0, rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_sh = 0, rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_oh = 0; // reco jet gen trks sub == 0
      int rcJet_gnTrk_sube1_ldjtTRkCorr = 0, rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym = 0, rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_1 = 0, rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_2 = 0, rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_3 = 0, rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_sh = 0, rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_oh = 0; // reco jet gen trks sub > 0
      int gn_ldjtTRkCorr = 0, gn_ldjtTRkCorr_RapAsym = 0, gn_ldjtTRkCorr_RapAsym_1 = 0, gn_ldjtTRkCorr_RapAsym_2 = 0, gn_ldjtTRkCorr_RapAsym_3 = 0, gn_ldjtTRkCorr_RapAsym_sh = 0, gn_ldjtTRkCorr_RapAsym_oh = 0; // for gen jets gen trks
      int gn_sube0_ldjtTRkCorr = 0, gn_sube0_ldjtTRkCorr_RapAsym = 0, gn_sube0_ldjtTRkCorr_RapAsym_1 = 0, gn_sube0_ldjtTRkCorr_RapAsym_2 = 0, gn_sube0_ldjtTRkCorr_RapAsym_3 = 0, gn_sube0_ldjtTRkCorr_RapAsym_sh = 0, gn_sube0_ldjtTRkCorr_RapAsym_oh = 0; // for gen jets gen trks sube == 0
      int gn_sube1_ldjtTRkCorr = 0, gn_sube1_ldjtTRkCorr_RapAsym = 0, gn_sube1_ldjtTRkCorr_RapAsym_1 = 0, gn_sube1_ldjtTRkCorr_RapAsym_2 = 0, gn_sube1_ldjtTRkCorr_RapAsym_3 = 0, gn_sube1_ldjtTRkCorr_RapAsym_sh = 0, gn_sube1_ldjtTRkCorr_RapAsym_oh = 0; // for gen jets gen trks sube > 0

      double TotalNtrk = 0., TotalNtrk_gen_sube0 = 0., TotalNtrk_gen_sube1 = 0.;
      double TotalNtrkOffline = 0.;
      
      for(int itrk = 0; itrk < Filtered_Trk_pT_vec_2D[ievt].size(); itrk++)
	{
	  TVector3 trk_vec = Filtered_Trk_pT_vec_2D[ievt][itrk];
	  double trk_w = Filtered_TrkW_vec_2D[ievt][itrk];
	  double trk_pt = trk_vec.Pt();
	  double trk_eta = trk_vec.Eta();
	  double trk_phi = trk_vec.Phi();
	  int trk_chg = Filtered_TrkCharge_vec_2D[ievt][itrk];
	  int trk_sube = Filtered_TrkSube_vec_2D[ievt][itrk];
	  int trk_ptbin = hTrkpTBin->FindBin(trk_pt) - 1;

	  if(trk_pt < trk_pt_min_cut || trk_pt > trk_pt_max_cut || fabs(trk_eta) > trk_eta_cut || trk_chg == 0) {std::cout<<"Tracks are out of  pT and eta cuts, please check"<<std::endl;}

	  if(!isrc)
	    {
	      double Trk_pT_Eta_Phi_ctbin[4] = {trk_pt, trk_eta, trk_phi, (double)ctbin};
	      hTrk_Signal_GenpT_Eta_Phi_ctbin_pTCut_W->Fill(Trk_pT_Eta_Phi_ctbin, trk_w);
	      
	      if(do_sube)
		{
		  if(trk_sube == 0)
		    {
		      hTrk_sube0_Signal_GenpT_Eta_Phi_ctbin_pTCut_W->Fill(Trk_pT_Eta_Phi_ctbin, trk_w);
		      TotalNtrk_gen_sube0 += trk_w;
		    }
		  else
		    {
		      hTrk_sube1_Signal_GenpT_Eta_Phi_ctbin_pTCut_W->Fill(Trk_pT_Eta_Phi_ctbin, trk_w);
		      TotalNtrk_gen_sube1 += trk_w;
		    }
		}
	    }

	  if(isrc)
	    {
	      //if(trk_ptbin == 0 || trk_ptbin == 1) TotalNtrk += trk_w; // sum MB trks
	      TotalNtrk += trk_w;
	      TotalNtrkOffline++;
	    }
	  
	  // for leading jet - trk
	  double Deta_ldjet_trk = trk_eta - ldJet_eta;
	  double Dphi_ldjet_trk = trk_phi - ldJet_phi;

	  if(Dphi_ldjet_trk > 1.5*TMath::Pi())
	    {
	      Dphi_ldjet_trk = Dphi_ldjet_trk - 2.0*TMath::Pi();
	    }
	  else if(Dphi_ldjet_trk < -0.5*TMath::Pi())
	    {
	      Dphi_ldjet_trk = Dphi_ldjet_trk + 2.0*TMath::Pi();
	    }
	  
	  // for subleading jet - trk
	  double Deta_sldjet_trk = trk_eta - sldJet_eta;
	  double Dphi_sldjet_trk = trk_phi - sldJet_phi;

	  if(Dphi_sldjet_trk > 1.5*TMath::Pi())
	    {
	      Dphi_sldjet_trk = Dphi_sldjet_trk - 2.0*TMath::Pi();
	    }
	  else if(Dphi_sldjet_trk < -0.5*TMath::Pi())
	    {
	      Dphi_sldjet_trk = Dphi_sldjet_trk + 2.0*TMath::Pi();
	    }
	  
	  // leading jet - trk
	  double ldjet_trk_signal[5] = {Deta_ldjet_trk, Dphi_ldjet_trk, (double)trk_ptbin, (double)ldJetflvor, (double)ctbin};
	  
	  // subleading jet - trk
	  double sldjet_trk_signal[5] = {Deta_sldjet_trk, Dphi_sldjet_trk, (double)trk_ptbin, (double)sldJetflvor, (double)ctbin};
	  
	  if(isrc) // for reco/data
	    {
	      ldjtTRkCorr++;
	      
	      if(do_sube_rcjetgntrk) // for reco jet gen trk sube
		{
		  if(trk_sube == 0){rcJet_gnTrk_sube0_ldjtTRkCorr++;}
		  else {rcJet_gnTrk_sube1_ldjtTRkCorr++;}
		}
	      
	      // leading jet -trk
	      hldJet_Trk_Signal->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
	      
	      // subleading jet -trk
	      hsldJet_Trk_Signal->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));

	      if(do_sube_rcjetgntrk) // for sube
		{
		  if(trk_sube == 0)
		    {
		      // leading reco jet - gentrk
		      hldJet_sube0_Signal->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
		    }
		  else
		    {
		      // leading reco jet - gentrk
		      hldJet_sube1_Signal->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
		    }
		}
	      
	      if(ldJet_eta > sldJet_eta) // crucial for rapidity asymmetry
		{
		  ldjtTRkCorr_RapAsym++;

		  if(do_sube_rcjetgntrk) // for reco jet gen trk sube
		    {
		      if(trk_sube == 0){rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym++;}
		      else {rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym++;}
		    }
			      
		  // leading jet -trk
		  hldJet_Trk_Signal_RapAsym->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
		  
		  // subleading jet -trk
		  hsldJet_Trk_Signal_RapAsym->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));

		  if(do_sube_rcjetgntrk) // for sube
		    {
		      if(trk_sube == 0)
			{
			  // leading reco jet - gentrk
			  hldJet_sube0_Signal_RapAsym->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			}
		      else
			{
			  // leading reco jet - gentrk
			  hldJet_sube1_Signal_RapAsym->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			}
		    }
		  
		  if(Deta_ldsldJet < 0.5) // mid rapidity
		    {
		      ldjtTRkCorr_RapAsym_1++;

		      if(do_sube_rcjetgntrk) // for reco jet gen trk sube
			{
			  if(trk_sube == 0){rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_1++;}
			  else {rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_1++;}
			}

		      // leading jet -trk
		      hldJet_Trk_Signal_RapAsym_1->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
		      
		      // subleading jet -trk
		      hsldJet_Trk_Signal_RapAsym_1->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));

		      if(do_sube_rcjetgntrk) // for sube
			{
			  if(trk_sube == 0)
			    {
			      // leading reco jet - gentrk
			      hldJet_sube0_Signal_RapAsym_1->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			    }
			  else
			    {
			      // leading reco jet - gentrk
			      hldJet_sube1_Signal_RapAsym_1->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			    }
			}
		    }
		  else if(Deta_ldsldJet > 0.5 && Deta_ldsldJet < 1.0) // intermediate rapidity
		    {
		      ldjtTRkCorr_RapAsym_2++;

		      if(do_sube_rcjetgntrk) // for reco jet gen trk sube
			{
			  if(trk_sube == 0){rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_2++;}
			  else {rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_2++;}
			}
				      
		      // leading jet -trk
		      hldJet_Trk_Signal_RapAsym_2->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
		      
		      // subleading jet -trk
		      hsldJet_Trk_Signal_RapAsym_2->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));

		      if(do_sube_rcjetgntrk) // for sube
			{
			  if(trk_sube == 0)
			    {
			      // leading reco jet - gentrk
			      hldJet_sube0_Signal_RapAsym_2->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			    }
			  else
			    {
			      // leading reco jet - gentrk
			      hldJet_sube1_Signal_RapAsym_2->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			    }
			}
		    }
		  else if(Deta_ldsldJet > 1.0) // higher rapidity
		    {
		      ldjtTRkCorr_RapAsym_3++;
		      
		      if(do_sube_rcjetgntrk) // for reco jet gen trk sube
			{
			  if(trk_sube == 0){rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_3++;}
			  else {rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_3++;}
			}
		      
		      // leading jet -trk
		      hldJet_Trk_Signal_RapAsym_3->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
		      
		      // subleading jet -trk
		      hsldJet_Trk_Signal_RapAsym_3->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));
		      
		      if(do_sube_rcjetgntrk) // for sube
			{
			  if(trk_sube == 0)
			    {
			      // leading reco jet - gentrk
			      hldJet_sube0_Signal_RapAsym_3->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			    }
			  else
			    {
			      // leading reco jet - gentrk
			      hldJet_sube1_Signal_RapAsym_3->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			    }
			}
		    }
		  if((ldJet_eta)*(sldJet_eta) > 0) // same hemisphere 
		    {
		      ldjtTRkCorr_RapAsym_sh++;
		      
		      if(do_sube_rcjetgntrk) // for reco jet gen trk sube
			{
			  if(trk_sube == 0){rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_sh++;}
			  else {rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_sh++;}
			}

		      // leading jet -trk
		      hldJet_Trk_Signal_RapAsym_sh->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
		      
		      // subleading jet -trk
		      hsldJet_Trk_Signal_RapAsym_sh->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));

		      if(do_sube_rcjetgntrk) // for sube
			{
			  if(trk_sube == 0)
			    {
			      // leading reco jet - gentrk
			      hldJet_sube0_Signal_RapAsym_sh->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			    }
			  else
			    {
			      // leading reco jet - gentrk
			      hldJet_sube1_Signal_RapAsym_sh->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			    }
			}
		    }
		  else if((ldJet_eta)*(sldJet_eta) < 0) // opposite hemisphere
		    {
		      ldjtTRkCorr_RapAsym_oh++;
		      
		      if(do_sube_rcjetgntrk) // for reco jet gen trk sube
			{
			  if(trk_sube == 0){rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_oh++;}
			  else {rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_oh++;}
			}
				      
		      // leading jet -trk
		      hldJet_Trk_Signal_RapAsym_oh->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
		      
		      // subleading jet -trk
		      hsldJet_Trk_Signal_RapAsym_oh->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));
		      
		      if(do_sube_rcjetgntrk) // for sube
			{
			  if(trk_sube == 0)
			    {
			      // leading reco jet - gentrk
			      hldJet_sube0_Signal_RapAsym_oh->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			    }
			  else
			    {
			      // leading reco jet - gentrk
			      hldJet_sube1_Signal_RapAsym_oh->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			    }
			}
		    }
		}
	    } // for reco/data
	  else
	    {
	      gn_ldjtTRkCorr++;
	      if(do_sube) // for sube
		{
		  if(trk_sube == 0) {gn_sube0_ldjtTRkCorr++;}
		  else {gn_sube1_ldjtTRkCorr++;}
		}
	      // leading jet -trk
	      hldGenJet_Trk_Signal->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
	      
	      // subleading jet -trk
	      hsldGenJet_Trk_Signal->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));

	      if(do_sube) // for sube
		{
		  if(trk_sube == 0)
		    {
		      // subleading jet -trk
		      hldGenJet_Trk_sube0_Signal->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
		      // subleading jet -trk 
		      hsldGenJet_Trk_sube0_Signal->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));
		    }
		  else
		    {
		      // leading jet -trk
		      hldGenJet_Trk_sube1_Signal->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
		      // subleading jet -trk
		      hsldGenJet_Trk_sube1_Signal->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));
		    }
		}
	      
	      if(ldJet_eta > sldJet_eta) // crucial for rapidity asymmetry
		{
		  gn_ldjtTRkCorr_RapAsym++;
		  if(do_sube) // for sube
		    {
		      if(trk_sube == 0) {gn_sube0_ldjtTRkCorr_RapAsym++;}
		      else {gn_sube1_ldjtTRkCorr_RapAsym++;}
		    }
		  // leading jet -trk
		  hldGenJet_Trk_Signal_RapAsym->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			  
		  // subleading jet -trk
		  hsldGenJet_Trk_Signal_RapAsym->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));

		  if(do_sube) // for sube
		    {
		      if(trk_sube == 0)
			{
			  // subleading jet -trk
			  hldGenJet_Trk_sube0_Signal_RapAsym->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			  // subleading jet -trk 
			  hsldGenJet_Trk_sube0_Signal_RapAsym->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));
			}
		      else
			{
			  // leading jet -trk
			  hldGenJet_Trk_sube1_Signal_RapAsym->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			  // subleading jet -trk
			  hsldGenJet_Trk_sube1_Signal_RapAsym->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));
			}
		    }
		  
		  if(Deta_ldsldJet < 0.5) // mid rapidity
		    {
		      gn_ldjtTRkCorr_RapAsym_1++;
		      if(do_sube) // for sube
			{
			  if(trk_sube == 0) {gn_sube0_ldjtTRkCorr_RapAsym_1++;}
			  else {gn_sube1_ldjtTRkCorr_RapAsym_1++;}
			}
		      
		      // leading jet -trk
		      hldGenJet_Trk_Signal_RapAsym_1->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
		      
		      // subleading jet -trk
		      hsldGenJet_Trk_Signal_RapAsym_1->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));
		      
		      if(do_sube) // for sube
			{
			  if(trk_sube == 0)
			    {
			      // subleading jet -trk
			      hldGenJet_Trk_sube0_Signal_RapAsym_1->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			      // subleading jet -trk 
			      hsldGenJet_Trk_sube0_Signal_RapAsym_1->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));
			    }
			  else
			    {
			      // leading jet -trk
			      hldGenJet_Trk_sube1_Signal_RapAsym_1->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			      // subleading jet -trk
			      hsldGenJet_Trk_sube1_Signal_RapAsym_1->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));
			    }
			}
		    }
		  else if(Deta_ldsldJet > 0.5 && Deta_ldsldJet < 1.0) // intermediate rapidity
		    {
		      gn_ldjtTRkCorr_RapAsym_2++;
		      if(do_sube) // for sube
			{
			  if(trk_sube == 0) {gn_sube0_ldjtTRkCorr_RapAsym_2++;}
			  else {gn_sube1_ldjtTRkCorr_RapAsym_2++;}
			}
		      
		      // leading jet -trk
		      hldGenJet_Trk_Signal_RapAsym_2->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
		      
		      // subleading jet -trk
		      hsldGenJet_Trk_Signal_RapAsym_2->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));

		      if(do_sube) // for sube
			{
			  if(trk_sube == 0)
			    {
			      // subleading jet -trk
			      hldGenJet_Trk_sube0_Signal_RapAsym_2->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			      // subleading jet -trk 
			      hsldGenJet_Trk_sube0_Signal_RapAsym_2->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));
			    }
			  else
			    {
			      // leading jet -trk
			      hldGenJet_Trk_sube1_Signal_RapAsym_2->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			      // subleading jet -trk
			      hsldGenJet_Trk_sube1_Signal_RapAsym_2->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));
			    }
			}
		    }
		  else if(Deta_ldsldJet > 1.0) // higher rapidity
		    {
		      gn_ldjtTRkCorr_RapAsym_3++;
		      if(do_sube) // for sube
			{
			  if(trk_sube == 0) {gn_sube0_ldjtTRkCorr_RapAsym_3++;}
			  else {gn_sube1_ldjtTRkCorr_RapAsym_3++;}
			}
 
		      // leading jet -trk
		      hldGenJet_Trk_Signal_RapAsym_3->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
		      
		      // subleading jet -trk
		      hsldGenJet_Trk_Signal_RapAsym_3->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));

		      if(do_sube) // for sube
			{
			  if(trk_sube == 0)
			    {
			      // subleading jet -trk
			      hldGenJet_Trk_sube0_Signal_RapAsym_3->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			      // subleading jet -trk 
			      hsldGenJet_Trk_sube0_Signal_RapAsym_3->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));
			    }
			  else
			    {
			      // leading jet -trk
			      hldGenJet_Trk_sube1_Signal_RapAsym_3->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			      // subleading jet -trk
			      hsldGenJet_Trk_sube1_Signal_RapAsym_3->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));
			    }
			}
		    }
		  if((ldJet_eta)*(sldJet_eta) > 0) // same hemisphere
		    {
		      gn_ldjtTRkCorr_RapAsym_sh++;
		      if(do_sube) // for sube
			{
			  if(trk_sube == 0) {gn_sube0_ldjtTRkCorr_RapAsym_sh++;}
			  else {gn_sube1_ldjtTRkCorr_RapAsym_sh++;}
			}
		      // leading jet -trk
		      hldGenJet_Trk_Signal_RapAsym_sh->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
		      
		      // subleading jet -trk
		      hsldGenJet_Trk_Signal_RapAsym_sh->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));

		      if(do_sube) // for sube
			{
			  if(trk_sube == 0)
			    {
			      // subleading jet -trk
			      hldGenJet_Trk_sube0_Signal_RapAsym_sh->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			      // subleading jet -trk 
			      hsldGenJet_Trk_sube0_Signal_RapAsym_sh->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));
			    }
			  else
			    {
			      // leading jet -trk
			      hldGenJet_Trk_sube1_Signal_RapAsym_sh->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			      // subleading jet -trk
			      hsldGenJet_Trk_sube1_Signal_RapAsym_sh->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));
			    }
			}
		    }
		  else if((ldJet_eta)*(sldJet_eta) < 0) // opposite hemisphere
		    {
		      gn_ldjtTRkCorr_RapAsym_oh++;
		      if(do_sube) // for sube
			{
			  if(trk_sube == 0) {gn_sube0_ldjtTRkCorr_RapAsym_oh++;}
			  else {gn_sube1_ldjtTRkCorr_RapAsym_oh++;}
			}
				      
		      // leading jet -trk
		      hldGenJet_Trk_Signal_RapAsym_oh->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
		      
		      // subleading jet -trk
		      hsldGenJet_Trk_Signal_RapAsym_oh->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));

		      if(do_sube) // for sube
			{
			  if(trk_sube == 0)
			    {
			      // subleading jet -trk
			      hldGenJet_Trk_sube0_Signal_RapAsym_oh->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			      // subleading jet -trk 
			      hsldGenJet_Trk_sube0_Signal_RapAsym_oh->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));
			    }
			  else
			    {
			      // leading jet -trk
			      hldGenJet_Trk_sube1_Signal_RapAsym_oh->Fill(ldjet_trk_signal, (evtw*ldJetW*trk_w));
			      // subleading jet -trk
			      hsldGenJet_Trk_sube1_Signal_RapAsym_oh->Fill(sldjet_trk_signal, (evtw*sldJetW*trk_w));
			    }
			}
		    }
		} // if(ldjetEta > sldjetEta)
	    } // for gen
	}// Track loop end
    
      /*
      // ntrk sum
      if(colliding_system == "PbPb" && !isrc && do_sube)
	{
	  if(ctbin == 0) hntrk_Signal_0->Fill(TotalNtrk);
	  else if(ctbin == 1) hntrk_Signal_1->Fill(TotalNtrk);
	  else if(ctbin == 2) hntrk_Signal_2->Fill(TotalNtrk);
	  else if(ctbin == 3) hntrk_Signal_3->Fill(TotalNtrk);
	  else if(ctbin == 4) hntrk_Signal_4->Fill(TotalNtrk);
	}
      */

      if(isrc)
	{
	  if(ctbin == 0) {hntrk_Signal_0->Fill(TotalNtrk, evtw); hntrkoffline_Signal_0->Fill(TotalNtrkOffline, evtw);}
	  else if(ctbin == 1) {hntrk_Signal_1->Fill(TotalNtrk, evtw); hntrkoffline_Signal_1->Fill(TotalNtrkOffline, evtw);}
	  else if(ctbin == 2) {hntrk_Signal_2->Fill(TotalNtrk, evtw); hntrkoffline_Signal_2->Fill(TotalNtrkOffline, evtw);}
	  else if(ctbin == 3) {hntrk_Signal_3->Fill(TotalNtrk, evtw); hntrkoffline_Signal_3->Fill(TotalNtrkOffline, evtw);}
	  else if(ctbin == 4) {hntrk_Signal_4->Fill(TotalNtrk, evtw); hntrkoffline_Signal_4->Fill(TotalNtrkOffline, evtw);}

	  if(ctbin == 0) hvtxz_Signal_0->Fill(vz, evtw);
          else if(ctbin == 1) hvtxz_Signal_1->Fill(vz, evtw);
          else if(ctbin == 2) hvtxz_Signal_2->Fill(vz, evtw);
          else if(ctbin == 3) hvtxz_Signal_3->Fill(vz, evtw);
          else if(ctbin == 4) hvtxz_Signal_4->Fill(vz, evtw);
	}
      else if(!isrc)
	{
	  if(do_sube)
	    {
	      if(ctbin == 0)
		{
		  hntrk_gen_sube0_Signal_0->Fill(TotalNtrk_gen_sube0, evtw);
		  hntrk_gen_sube1_Signal_0->Fill(TotalNtrk_gen_sube1, evtw);
		}
	      else if(ctbin == 1)
		{
		  hntrk_gen_sube0_Signal_1->Fill(TotalNtrk_gen_sube0, evtw);
                  hntrk_gen_sube1_Signal_1->Fill(TotalNtrk_gen_sube1, evtw);
		}
	      else if(ctbin == 2)
		{
		 hntrk_gen_sube0_Signal_2->Fill(TotalNtrk_gen_sube0, evtw);
		 hntrk_gen_sube1_Signal_2->Fill(TotalNtrk_gen_sube1, evtw);
		}
	      else if(ctbin == 3)
		{
		  hntrk_gen_sube0_Signal_3->Fill(TotalNtrk_gen_sube0, evtw);
                  hntrk_gen_sube1_Signal_3->Fill(TotalNtrk_gen_sube1, evtw);
		}
	      else if(ctbin == 4)
		{
		  hntrk_gen_sube0_Signal_4->Fill(TotalNtrk_gen_sube0, evtw);
                  hntrk_gen_sube1_Signal_4->Fill(TotalNtrk_gen_sube1, evtw);
		}
	    }
	}
    
      // for jet pair
      double jet_pair[3] = {0, (double)ldJetflvor, (double)ctbin};
      if(isrc)
	{
	  if(ldjtTRkCorr > 0) {hldsld_Jet_pair_Signal->Fill(jet_pair, (evtw*ldJetW));}
	  if(do_sube_rcjetgntrk) // for sube
	    {
	      if(rcJet_gnTrk_sube0_ldjtTRkCorr > 0){hldsld_Jet_pair_sube0_Signal->Fill(jet_pair, (evtw*ldJetW));}
	      if(rcJet_gnTrk_sube1_ldjtTRkCorr > 0){hldsld_Jet_pair_sube1_Signal->Fill(jet_pair, (evtw*ldJetW));}
	    }
	  if(ldJet_eta > sldJet_eta) // crucial for rapidity asymmetry                                                                             
	    {
	      if(ldjtTRkCorr_RapAsym > 0){hldsld_Jet_pair_Signal_RapAsym->Fill(jet_pair, (evtw*ldJetW));}
	      if(do_sube_rcjetgntrk) // for sube
		{
		  if(rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym > 0){hldsld_Jet_pair_sube0_Signal_RapAsym->Fill(jet_pair, (evtw*ldJetW));}
		  if(rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym > 0){hldsld_Jet_pair_sube1_Signal_RapAsym->Fill(jet_pair, (evtw*ldJetW));}
		}
	      if(Deta_ldsldJet < 0.5) // mid rapidity
		{
		  if(ldjtTRkCorr_RapAsym_1 > 0){hldsld_Jet_pair_Signal_RapAsym_1->Fill(jet_pair, (evtw*ldJetW));}
		  if(do_sube_rcjetgntrk) // for sube
		    {
		      if(rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_1 > 0){hldsld_Jet_pair_sube0_Signal_RapAsym_1->Fill(jet_pair, (evtw*ldJetW));}
		      if(rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_1 > 0){hldsld_Jet_pair_sube1_Signal_RapAsym_1->Fill(jet_pair, (evtw*ldJetW));}
		    }
		}
	      else if(Deta_ldsldJet > 0.5 && Deta_ldsldJet < 1.0) // intermediate rapidity
		{
		  if(ldjtTRkCorr_RapAsym_2 > 0){hldsld_Jet_pair_Signal_RapAsym_2->Fill(jet_pair, (evtw*ldJetW));}
		  if(do_sube_rcjetgntrk) // for sube
		    {
		      if(rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_2 > 0){hldsld_Jet_pair_sube0_Signal_RapAsym_2->Fill(jet_pair, (evtw*ldJetW));}
		      if(rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_2 > 0){hldsld_Jet_pair_sube1_Signal_RapAsym_2->Fill(jet_pair, (evtw*ldJetW));}
		    }
		}
	      else if(Deta_ldsldJet > 1.0) // higher rapidity
		{
		  if(ldjtTRkCorr_RapAsym_3 > 0){hldsld_Jet_pair_Signal_RapAsym_3->Fill(jet_pair, (evtw*ldJetW));}
		  if(do_sube_rcjetgntrk) // for sube
		    {
		      if(rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_3 > 0){hldsld_Jet_pair_sube0_Signal_RapAsym_3->Fill(jet_pair, (evtw*ldJetW));}
		      if(rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_3 > 0){hldsld_Jet_pair_sube1_Signal_RapAsym_3->Fill(jet_pair, (evtw*ldJetW));}
		    }
		}
	      if((ldJet_eta)*(sldJet_eta) > 0) // same hemisphere
		{
		  if(ldjtTRkCorr_RapAsym_sh > 0){hldsld_Jet_pair_Signal_RapAsym_sh->Fill(jet_pair, (evtw*ldJetW));}
		  if(do_sube_rcjetgntrk) // for sube
		    {
		      if(rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_sh > 0){hldsld_Jet_pair_sube0_Signal_RapAsym_sh->Fill(jet_pair, (evtw*ldJetW));}
		      if(rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_sh > 0){hldsld_Jet_pair_sube1_Signal_RapAsym_sh->Fill(jet_pair, (evtw*ldJetW));}
		    }
		}
	      else if((ldJet_eta)*(sldJet_eta) < 0) // opposite hemisphere
		{
		  if(ldjtTRkCorr_RapAsym_oh > 0){hldsld_Jet_pair_Signal_RapAsym_oh->Fill(jet_pair, (evtw*ldJetW));}
		  if(do_sube_rcjetgntrk) // for sube
		    {
		      if(rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_oh > 0){hldsld_Jet_pair_sube0_Signal_RapAsym_oh->Fill(jet_pair, (evtw*ldJetW));}
		      if(rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_oh > 0){hldsld_Jet_pair_sube1_Signal_RapAsym_oh->Fill(jet_pair, (evtw*ldJetW));}
		    }		  
		}
	    }
	}
      else
	{
	  if(gn_ldjtTRkCorr > 0){hldsld_GenJet_pair_Signal->Fill(jet_pair, (evtw*ldJetW));}
	  if(do_sube) // for sube
	    {
	      if(gn_sube0_ldjtTRkCorr > 0){hldsld_GenJet_pair_sube0_Signal->Fill(jet_pair, (evtw*ldJetW));}
	      if(gn_sube1_ldjtTRkCorr > 0){hldsld_GenJet_pair_sube1_Signal->Fill(jet_pair, (evtw*ldJetW));} 
	    }
	  if(ldJet_eta > sldJet_eta) // crucial for rapidity asymmetry                                                                             
	    {
	      if(gn_ldjtTRkCorr_RapAsym > 0){hldsld_GenJet_pair_Signal_RapAsym->Fill(jet_pair, (evtw*ldJetW));}
	      if(do_sube) // for sube
		{
		  if(gn_sube0_ldjtTRkCorr_RapAsym > 0){hldsld_GenJet_pair_sube0_Signal_RapAsym->Fill(jet_pair, (evtw*ldJetW));}
		  if(gn_sube1_ldjtTRkCorr_RapAsym > 0){hldsld_GenJet_pair_sube1_Signal_RapAsym->Fill(jet_pair, (evtw*ldJetW));} 
		} 
	      if(Deta_ldsldJet < 0.5) // mid rapidity
		{
		  if(gn_ldjtTRkCorr_RapAsym_1 > 0){hldsld_GenJet_pair_Signal_RapAsym_1->Fill(jet_pair, (evtw*ldJetW));}
		  if(do_sube) // for sube
		    {
		      if(gn_sube0_ldjtTRkCorr_RapAsym_1 > 0){hldsld_GenJet_pair_sube0_Signal_RapAsym_1->Fill(jet_pair, (evtw*ldJetW));}
		      if(gn_sube1_ldjtTRkCorr_RapAsym_1 > 0){hldsld_GenJet_pair_sube1_Signal_RapAsym_1->Fill(jet_pair, (evtw*ldJetW));} 
		    } 
		}
	      else if(Deta_ldsldJet > 0.5 && Deta_ldsldJet < 1.0) // intermediate rapidity
		{
		  if(gn_ldjtTRkCorr_RapAsym_2 > 0){hldsld_GenJet_pair_Signal_RapAsym_2->Fill(jet_pair, (evtw*ldJetW));}
		  if(do_sube) // for sube
		    {
		      if(gn_sube0_ldjtTRkCorr_RapAsym_2 > 0){hldsld_GenJet_pair_sube0_Signal_RapAsym_2->Fill(jet_pair, (evtw*ldJetW));}
		      if(gn_sube1_ldjtTRkCorr_RapAsym_2 > 0){hldsld_GenJet_pair_sube1_Signal_RapAsym_2->Fill(jet_pair, (evtw*ldJetW));} 
		    } 
		}
	      else if(Deta_ldsldJet > 1.0) // higher rapidity
		{
		  if(gn_ldjtTRkCorr_RapAsym_3 > 0){hldsld_GenJet_pair_Signal_RapAsym_3->Fill(jet_pair, (evtw*ldJetW));}
		  if(do_sube) // for sube
		    {
		      if(gn_sube0_ldjtTRkCorr_RapAsym_3 > 0){hldsld_GenJet_pair_sube0_Signal_RapAsym_3->Fill(jet_pair, (evtw*ldJetW));}
		      if(gn_sube1_ldjtTRkCorr_RapAsym_3 > 0){hldsld_GenJet_pair_sube1_Signal_RapAsym_3->Fill(jet_pair, (evtw*ldJetW));} 
		    } 
		}
	      if((ldJet_eta)*(sldJet_eta) > 0) // same hemisphere
		{
		  if(gn_ldjtTRkCorr_RapAsym_sh > 0){hldsld_GenJet_pair_Signal_RapAsym_sh->Fill(jet_pair, (evtw*ldJetW));}
		  if(do_sube) // for sube
		    {
		      if(gn_sube0_ldjtTRkCorr_RapAsym_sh > 0){hldsld_GenJet_pair_sube0_Signal_RapAsym_sh->Fill(jet_pair, (evtw*ldJetW));}
		      if(gn_sube1_ldjtTRkCorr_RapAsym_sh > 0){hldsld_GenJet_pair_sube1_Signal_RapAsym_sh->Fill(jet_pair, (evtw*ldJetW));} 
		    } 
		}
	      else if((ldJet_eta)*(sldJet_eta) < 0) // opposite hemisphere
		{
		  if(gn_ldjtTRkCorr_RapAsym_oh > 0){hldsld_GenJet_pair_Signal_RapAsym_oh->Fill(jet_pair, (evtw*ldJetW));}
		  if(do_sube) // for sube
		    {
		      if(gn_sube0_ldjtTRkCorr_RapAsym_oh > 0){hldsld_GenJet_pair_sube0_Signal_RapAsym_oh->Fill(jet_pair, (evtw*ldJetW));}
		      if(gn_sube1_ldjtTRkCorr_RapAsym_oh > 0){hldsld_GenJet_pair_sube1_Signal_RapAsym_oh->Fill(jet_pair, (evtw*ldJetW));} 
		    } 
		}
	    } // if (ldeta > subld eta)
	} // for gen
    } // event/leading subleadin jet loop end

  std::cout<<endl;
  std::cout<<"TotaldijetEvent is: "<<TotaldijetEvent<<std::endl;
  std::cout<<endl;
  if(isrc) {std::cout<<"~~~~~end leading subleading jet tracks signal correlation~~~~~~~~"<<std::endl;}
  else {std::cout<<"~~~~~end gen leading subleading jet tracks signal correlation~~~~~~~~"<<std::endl;}
  std::cout<<endl;
}// function loop

void Jet_Track_mixing_corr_ldsld(const TString& colliding_system, const std::vector<double>& Evtw_vec_1D, const std::vector<int>& HiBin_vec_1D, const std::vector<int>& HiBinValue_vec_1D, const std::vector<double>& Vertexz_vec_1D, const std::vector<Long64_t>& Evtno_vec_1D, const std::vector<int>& EvtCount_vec_1D, const std::vector<TVector3>& Filtered_ldjet_vec_1D, const std::vector<int>& Filtered_ldrefpartonB_vec_1D, const std::vector<double>& Filtered_ldJetW_vec_1D, const std::vector<TVector3>& Filtered_sldjet_vec_1D, const std::vector<int>& Filtered_sldrefpartonB_vec_1D, const std::vector<double>& Filtered_sldJetW_vec_1D, const std::vector<std::vector<TVector3>>& Filtered_Trk_pT_vec_2D, const std::vector<std::vector<double>>& Filtered_TrkW_vec_2D, const std::vector<std::vector<int>>& Filtered_TrkCharge_vec_2D, const std::vector<std::vector<int>>& Filtered_TrkSube_vec_2D, const bool& isrc, const bool& do_sube, const bool& do_sube_rcjetgntrk)
{
  std::cout<<endl;

  if(isrc) {std::cout<<"~~~~~start leading subleading jet tracks mixing correlation~~~~~~~~"<<std::endl;}
  else {std::cout<<"~~~~~start gen leading subleading jet tracks mixing correlation~~~~~~~~"<<std::endl;}

  int TotaldijetEvent = 0;
  
  for(int ievt = 0; ievt < Vertexz_vec_1D.size(); ievt++)
    {
      TotaldijetEvent++;
      
      if(ievt%10000 == 0) std::cout<<ievt<<"  events running for mixing correlation of total events "<< Vertexz_vec_1D.size() <<std::endl;
      
      double evtw = Evtw_vec_1D[ievt];
      //double evtw_i = Evtw_vec_1D[ievt];
      int ctbin = HiBin_vec_1D[ievt];
      int centbin = HiBinValue_vec_1D[ievt];
      int ctbinn = hCentBin->FindBin(centbin) -1;
      double vz = Vertexz_vec_1D[ievt];
      Long64_t evtno = Evtno_vec_1D[ievt];
      
      if(colliding_system == "PbPb")
	{
	  if(ctbin != ctbinn) {std::cout<<"ctbin and ctbinn are not same in PbPb, please check"<<std::endl;}
	}
      
      TVector3 ldjet_vec = Filtered_ldjet_vec_1D[ievt];
      double ldJetW = Filtered_ldJetW_vec_1D[ievt];
      int ldJetflvor = Filtered_ldrefpartonB_vec_1D[ievt];
      double ldJet_pt = ldjet_vec.Pt();
      double ldJet_eta = ldjet_vec.Eta();
      double ldJet_phi = ldjet_vec.Phi();
      
      TVector3 sldjet_vec = Filtered_sldjet_vec_1D[ievt];
      double sldJetW = Filtered_sldJetW_vec_1D[ievt];
      int sldJetflvor = Filtered_sldrefpartonB_vec_1D[ievt];
      double sldJet_pt = sldjet_vec.Pt();
      double sldJet_eta = sldjet_vec.Eta();
      double sldJet_phi = sldjet_vec.Phi();

      int ntrk = Filtered_Trk_pT_vec_2D[ievt].size();
      
      double Deta_ldsldJet = fabs(ldJet_eta - sldJet_eta);

      if(fabs(ldJet_eta) > 0.5) continue; // for test
      
      /*
      // mixing algorithm
      int mixstart = ievt+1;
      int mixend   = (int)Vertexz_vec_1D.size();
      
      if(mixstart > (0.5*(Vertexz_vec_1D.size())))
	{
	  mixstart = 0;
	  //mixend   = (int)Vertexz_vec_1D.size();
	  mixend = ievt -1;
	}
      
      int nmix = 0;
      */

      std::srand(std::time(0));
      std::random_device rd;
      std::mt19937 g(rd());
      
      std::vector<int> SelectMBEventCountIDForMixing;
      
      // Shuffle File 2 events once before matching
      std::vector<int> shuffled_EvtCount_vec_1D = EvtCount_vec_1D;  
      std::shuffle(shuffled_EvtCount_vec_1D.begin(), shuffled_EvtCount_vec_1D.end(), g);

      for(int ievt_j = 0; ievt_j < shuffled_EvtCount_vec_1D.size(); ievt_j++) // loop over events to select desire mix events
	{
	  if(SelectMBEventCountIDForMixing.size() >= bkgFactor) break;  // Stop early if enough events are found
	  
	  int EventCountID = shuffled_EvtCount_vec_1D[ievt_j];
	  
	  int ctbin_j = HiBin_vec_1D[EventCountID];
          int centbin_j = HiBinValue_vec_1D[EventCountID];
          double vz_j = Vertexz_vec_1D[EventCountID];
	  int ntrk_j = Filtered_Trk_pT_vec_2D[EventCountID].size();
	  Long64_t evtno_j = Evtno_vec_1D[EventCountID];

	  if(colliding_system == "PbPb")
	    {
	      //if((fabs(vz - vz_j) < Dvz_cut) && (ctbin == ctbin_j) && (fabs(centbin - centbin_j) <= Dhibin_cut) && (evtno != evtno_j) && (fabs(ntrk - ntrk_j) <= Dntrk_cut))
	      if((fabs(vz - vz_j) <= Dvz_cut) && (ctbin == ctbin_j) && (fabs(centbin - centbin_j) <= Dhibin_cut) && (evtno != evtno_j))
		{
		  SelectMBEventCountIDForMixing.push_back(EventCountID); // push_back the selected event counts for mixing
		}
	    }
	  else if(colliding_system == "pp")
	    {
	      if((fabs(vz - vz_j) <= Dvz_cut) && (fabs(ntrk - ntrk_j) <= Dntrk_cut) && (evtno != evtno_j))
		{
		  SelectMBEventCountIDForMixing.push_back(EventCountID); // push_back the selected event counts for mixing
		}
	    }
	}
      
      if(ievt%10000 == 0)
	{
	  std::cout<<"number of mixing events found for "<<ievt<<"th dijet event is: "<<SelectMBEventCountIDForMixing.size()<<std::endl;
	  //std::cout<<"Selected events ID are: ";
	  //for(int jevt = 0; jevt < SelectMBEventCountIDForMixing.size(); jevt++){std::cout<<SelectMBEventCountIDForMixing[jevt]<<", "; }
	}

      //for(int jevt = mixstart; jevt < mixend; jevt++) // 2nd event loop
      //if(ievt == jevt) continue;

      for(int jevt = 0; jevt < SelectMBEventCountIDForMixing.size(); jevt++) // 2nd event loop
	{
	  int MixEvtID = SelectMBEventCountIDForMixing[jevt];
	  
	  double evtw_j = Evtw_vec_1D[MixEvtID];
	  
	  //if(ievt%200 == 0){std::cout<<ievt<<"  "<<MixEvtID<<"  "<<ctbin<<"  "<<HiBin_vec_1D[MixEvtID]<<"  "<<centbin<<"  "<<HiBinValue_vec_1D[MixEvtID]<<"  "<<vz<<"  "<<Vertexz_vec_1D[MixEvtID]<<"  "<<evtno<<"  "<<Evtno_vec_1D[MixEvtID]<<"  "<<ntrk<<"  "<<Filtered_Trk_pT_vec_2D[MixEvtID].size()<<std::endl;}

	  /*
	  double evtw = evtw_i;
	  if(colliding_system == "PbPb") {evtw = (evtw_i)*(evtw_j);}
	  */

	  int ldjtTRkCorr = 0, ldjtTRkCorr_RapAsym = 0, ldjtTRkCorr_RapAsym_1 = 0, ldjtTRkCorr_RapAsym_2 = 0, ldjtTRkCorr_RapAsym_3 = 0, ldjtTRkCorr_RapAsym_sh = 0, ldjtTRkCorr_RapAsym_oh = 0; // for reco jet reco trks
	  int rcJet_gnTrk_sube0_ldjtTRkCorr = 0, rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym = 0, rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_1 = 0, rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_2 = 0, rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_3 = 0, rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_sh = 0, rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_oh = 0; // reco jet gen trks sub == 0
	  int rcJet_gnTrk_sube1_ldjtTRkCorr = 0, rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym = 0, rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_1 = 0, rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_2 = 0, rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_3 = 0, rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_sh = 0, rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_oh = 0; // reco jet gen trks sub > 0
	  int gn_ldjtTRkCorr = 0, gn_ldjtTRkCorr_RapAsym = 0, gn_ldjtTRkCorr_RapAsym_1 = 0, gn_ldjtTRkCorr_RapAsym_2 = 0, gn_ldjtTRkCorr_RapAsym_3 = 0, gn_ldjtTRkCorr_RapAsym_sh = 0, gn_ldjtTRkCorr_RapAsym_oh = 0; // for gen jets gen trks
	  int gn_sube0_ldjtTRkCorr = 0, gn_sube0_ldjtTRkCorr_RapAsym = 0, gn_sube0_ldjtTRkCorr_RapAsym_1 = 0, gn_sube0_ldjtTRkCorr_RapAsym_2 = 0, gn_sube0_ldjtTRkCorr_RapAsym_3 = 0, gn_sube0_ldjtTRkCorr_RapAsym_sh = 0, gn_sube0_ldjtTRkCorr_RapAsym_oh = 0; // for gen jets gen trks sube == 0
	  int gn_sube1_ldjtTRkCorr = 0, gn_sube1_ldjtTRkCorr_RapAsym = 0, gn_sube1_ldjtTRkCorr_RapAsym_1 = 0, gn_sube1_ldjtTRkCorr_RapAsym_2 = 0, gn_sube1_ldjtTRkCorr_RapAsym_3 = 0, gn_sube1_ldjtTRkCorr_RapAsym_sh = 0, gn_sube1_ldjtTRkCorr_RapAsym_oh = 0; // for gen jets gen trks sube > 0
	  
	  double TotalNtrk = 0.;

	  for(int itrk = 0; itrk < Filtered_Trk_pT_vec_2D[MixEvtID].size(); itrk++)
	    {
	      TVector3 trk_vec = Filtered_Trk_pT_vec_2D[MixEvtID][itrk];
	      double trk_w = Filtered_TrkW_vec_2D[MixEvtID][itrk];
	      double trk_pt = trk_vec.Pt();
	      double trk_eta = trk_vec.Eta();
	      double trk_phi = trk_vec.Phi();
	      int trk_sube = Filtered_TrkSube_vec_2D[MixEvtID][itrk];
	      int trk_ptbin = hTrkpTBin->FindBin(trk_pt) - 1;

	      if(!isrc)
		{
		  double Trk_pT_Eta_Phi_ctbin[4] = {trk_pt, trk_eta, trk_phi, (double)ctbin};
		  hTrk_Mixing_GenpT_Eta_Phi_ctbin_pTCut_W->Fill(Trk_pT_Eta_Phi_ctbin, trk_w);
		  
		  if(do_sube)
		    {
		      if(trk_sube == 0){hTrk_sube0_Mixing_GenpT_Eta_Phi_ctbin_pTCut_W->Fill(Trk_pT_Eta_Phi_ctbin, trk_w);}
		      else {hTrk_sube1_Mixing_GenpT_Eta_Phi_ctbin_pTCut_W->Fill(Trk_pT_Eta_Phi_ctbin, trk_w);}
		      
		      //if(colliding_system == "PbPb" && trk_sube != 0 && trk_pt >= 1 && trk_pt < 2) TotalNtrk += trk_w; // sum MB trks
		    }
		}
	      if(isrc)
		{
		  //if(trk_ptbin == 0) TotalNtrk += trk_w;
		  TotalNtrk += trk_w;
		}
	      // for leading jet - trk
	      double Deta_ldjet_trk = trk_eta - ldJet_eta;
	      double Dphi_ldjet_trk = trk_phi - ldJet_phi;
	      
	      if(Dphi_ldjet_trk > 1.5*TMath::Pi())
		{
		  Dphi_ldjet_trk = Dphi_ldjet_trk - 2.0*TMath::Pi();
		}
	      else if(Dphi_ldjet_trk < -0.5*TMath::Pi())
		{
		  Dphi_ldjet_trk = Dphi_ldjet_trk + 2.0*TMath::Pi();
		}
	      
	      // for subleading jet - trk
	      double Deta_sldjet_trk = trk_eta - sldJet_eta;
	      double Dphi_sldjet_trk = trk_phi - sldJet_phi;
	      
	      if(Dphi_sldjet_trk > 1.5*TMath::Pi())
		{
		  Dphi_sldjet_trk = Dphi_sldjet_trk - 2.0*TMath::Pi();
		}
	      else if(Dphi_sldjet_trk < -0.5*TMath::Pi())
		{
		  Dphi_sldjet_trk = Dphi_sldjet_trk + 2.0*TMath::Pi();
		}
	      
	      // leading jet -trk
	      double ldjet_trk_mixing[5] = {Deta_ldjet_trk, Dphi_ldjet_trk, (double)trk_ptbin, (double)ldJetflvor, (double)ctbin};
	      
	      // subleading jet -trk
	      double sldjet_trk_mixing[5] = {Deta_sldjet_trk, Dphi_sldjet_trk, (double)trk_ptbin, (double)sldJetflvor, (double)ctbin};
	      
	      if(isrc) // for reco/data
		{
		  ldjtTRkCorr++;

		  if(do_sube_rcjetgntrk) // for reco jet gen trk sube
		    {
		      if(trk_sube == 0){rcJet_gnTrk_sube0_ldjtTRkCorr++;}
		      else {rcJet_gnTrk_sube1_ldjtTRkCorr++;}
		    }
		  
		  // leading jet -trk
		  hldJet_Trk_Mixing->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
		  
		  // subleading jet -trk
		  hsldJet_Trk_Mixing->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));

		  if(do_sube_rcjetgntrk) // for sube
		    {
		      if(trk_sube == 0)
			{
			  // leading reco jet - gentrk
			  hldJet_sube0_Mixing->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			}
		      else
			{
			  // leading reco jet - gentrk
			  hldJet_sube1_Mixing->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			}
		    }
		  if(ldJet_eta > sldJet_eta) // crucial for rapidity asymmetry
		    {
		      ldjtTRkCorr_RapAsym++;

		      if(do_sube_rcjetgntrk) // for reco jet gen trk sube
			{
			  if(trk_sube == 0){rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym++;}
			  else {rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym++;}
			}
		      
		      // leading jet -trk
		      hldJet_Trk_Mixing_RapAsym->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));

		      // subleading jet -trk
		      hsldJet_Trk_Mixing_RapAsym->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));

		      if(do_sube_rcjetgntrk) // for sube
			{
			  if(trk_sube == 0)
			    {
			      // leading reco jet - gentrk
			      hldJet_sube0_Mixing_RapAsym->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			    }
			  else
			    {
			      // leading reco jet - gentrk
			      hldJet_sube1_Mixing_RapAsym->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			    }
			}
		      if(Deta_ldsldJet < 0.5) // mid rapidity
			{
			  ldjtTRkCorr_RapAsym_1++;

			  if(do_sube_rcjetgntrk) // for reco jet gen trk sube
			    {
			      if(trk_sube == 0){rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_1++;}
			      else {rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_1++;}
			    }
			  
			  // leading jet -trk
			  hldJet_Trk_Mixing_RapAsym_1->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldJet_Trk_Mixing_RapAsym_1->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));

			  if(do_sube_rcjetgntrk) // for sube
			    {
			      if(trk_sube == 0)
				{
				  // leading reco jet - gentrk
				  hldJet_sube0_Mixing_RapAsym_1->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				}
			      else
				{
				  // leading reco jet - gentrk
				  hldJet_sube1_Mixing_RapAsym_1->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				}
			    }
			}
		      else if(Deta_ldsldJet > 0.5 && Deta_ldsldJet < 1.0) // intermediate rapidity
			{
			  ldjtTRkCorr_RapAsym_2++;
			  
			  if(do_sube_rcjetgntrk) // for reco jet gen trk sube
			    {
			      if(trk_sube == 0){rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_2++;}
			      else {rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_2++;}
			    }
			  
			  // leading jet -trk
			  hldJet_Trk_Mixing_RapAsym_2->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldJet_Trk_Mixing_RapAsym_2->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));

			  if(do_sube_rcjetgntrk) // for sube
			    {
			      if(trk_sube == 0)
				{
				  // leading reco jet - gentrk
				  hldJet_sube0_Mixing_RapAsym_2->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				}
			      else
				{
				  // leading reco jet - gentrk
				  hldJet_sube1_Mixing_RapAsym_2->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				}
			    }
			}
		      else if(Deta_ldsldJet > 1.0) // higher rapidity
			{
			  ldjtTRkCorr_RapAsym_3++;

			  if(do_sube_rcjetgntrk) // for reco jet gen trk sube
			    {
			      if(trk_sube == 0){rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_3++;}
			      else {rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_3++;}
			    }
			  
			  // leading jet -trk
			  hldJet_Trk_Mixing_RapAsym_3->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldJet_Trk_Mixing_RapAsym_3->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));

			  if(do_sube_rcjetgntrk) // for sube
			    {
			      if(trk_sube == 0)
				{
				  // leading reco jet - gentrk
				  hldJet_sube0_Mixing_RapAsym_3->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				}
			      else
				{
				  // leading reco jet - gentrk
				  hldJet_sube1_Mixing_RapAsym_3->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				}
			    }
			}
		      if((ldJet_eta)*(sldJet_eta) > 0) // same hemisphere 
			{
			  ldjtTRkCorr_RapAsym_sh++;
			  
			  if(do_sube_rcjetgntrk) // for reco jet gen trk sube
			    {
			      if(trk_sube == 0){rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_sh++;}
			      else {rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_sh++;}
			    }
			  
			  // leading jet -trk
			  hldJet_Trk_Mixing_RapAsym_sh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldJet_Trk_Mixing_RapAsym_sh->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));

			  if(do_sube_rcjetgntrk) // for sube
			    {
			      if(trk_sube == 0)
				{
				  // leading reco jet - gentrk
				  hldJet_sube0_Mixing_RapAsym_sh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				}
			      else
				{
				  // leading reco jet - gentrk
				  hldJet_sube1_Mixing_RapAsym_sh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				}
			    }
			}
		      else if((ldJet_eta)*(sldJet_eta) < 0) // opposite hemisphere
			{
			  ldjtTRkCorr_RapAsym_oh++;

			  if(do_sube_rcjetgntrk) // for reco jet gen trk sube
			    {
			      if(trk_sube == 0){rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_oh++;}
			      else {rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_oh++;}
			    }
						  
			  // leading jet -trk
			  hldJet_Trk_Mixing_RapAsym_oh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldJet_Trk_Mixing_RapAsym_oh->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
			  
			  if(do_sube_rcjetgntrk) // for sube
			    {
			      if(trk_sube == 0)
				{
				  // leading reco jet - gentrk
				  hldJet_sube0_Mixing_RapAsym_oh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				}
			      else
				{
				  // leading reco jet - gentrk
				  hldJet_sube1_Mixing_RapAsym_oh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				}
			    }
			}
		    }
		} // for reco/data
	      else
		{
		  gn_ldjtTRkCorr++;
		  if(do_sube) // for sube
		    {
		      if(trk_sube == 0) {gn_sube0_ldjtTRkCorr++;}
		      else {gn_sube1_ldjtTRkCorr++;}
		    }
		  
		  // leading jet -trk
		  hldGenJet_Trk_Mixing->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
		  
		  // subleading jet -trk
		  hsldGenJet_Trk_Mixing->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));

		  if(do_sube) // for sube
		    {
		      if(trk_sube == 0)
			{
			  // subleading jet -trk
			  hldGenJet_Trk_sube0_Mixing->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  // subleading jet -trk 
			  hsldGenJet_Trk_sube0_Mixing->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
			}
		      else
			{
			  // leading jet -trk
			  hldGenJet_Trk_sube1_Mixing->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  // subleading jet -trk
			  hsldGenJet_Trk_sube1_Mixing->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
			}
		    }
		  
		  if(ldJet_eta > sldJet_eta) // crucial for rapidity asymmetry
		    {
		      gn_ldjtTRkCorr_RapAsym++;
		      if(do_sube) // for sube
			{
			  if(trk_sube == 0) {gn_sube0_ldjtTRkCorr_RapAsym++;}
			  else {gn_sube1_ldjtTRkCorr_RapAsym++;}
			}

		      // leading jet -trk
		      hldGenJet_Trk_Mixing_RapAsym->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
		      
		      // subleading jet -trk
		      hsldGenJet_Trk_Mixing_RapAsym->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));

		      if(do_sube) // for sube
			{
			  if(trk_sube == 0)
			    {
			      // subleading jet -trk
			      hldGenJet_Trk_sube0_Mixing_RapAsym->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			      // subleading jet -trk 
			      hsldGenJet_Trk_sube0_Mixing_RapAsym->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
			    }
			  else
			    {
			      // leading jet -trk
			      hldGenJet_Trk_sube1_Mixing_RapAsym->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			      // subleading jet -trk
			      hsldGenJet_Trk_sube1_Mixing_RapAsym->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
			    }
			}
		      
		      if(Deta_ldsldJet < 0.5) // mid rapidity
			{
			  gn_ldjtTRkCorr_RapAsym_1++;
			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0) {gn_sube0_ldjtTRkCorr_RapAsym_1++;}
			      else {gn_sube1_ldjtTRkCorr_RapAsym_1++;}
			    }
			  
			  // leading jet -trk
			  hldGenJet_Trk_Mixing_RapAsym_1->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldGenJet_Trk_Mixing_RapAsym_1->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));

			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0)
				{
				  // subleading jet -trk
				  hldGenJet_Trk_sube0_Mixing_RapAsym_1->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk 
				  hsldGenJet_Trk_sube0_Mixing_RapAsym_1->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			      else
				{
				  // leading jet -trk
				  hldGenJet_Trk_sube1_Mixing_RapAsym_1->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk
				  hsldGenJet_Trk_sube1_Mixing_RapAsym_1->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			    }
			}
		      else if(Deta_ldsldJet > 0.5 && Deta_ldsldJet < 1.0) // intermediate rapidity
			{
			  gn_ldjtTRkCorr_RapAsym_2++;
			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0) {gn_sube0_ldjtTRkCorr_RapAsym_2++;}
			      else {gn_sube1_ldjtTRkCorr_RapAsym_2++;}
			    }
  
			  // leading jet -trk
			  hldGenJet_Trk_Mixing_RapAsym_2->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldGenJet_Trk_Mixing_RapAsym_2->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));

			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0)
				{
				  // subleading jet -trk
				  hldGenJet_Trk_sube0_Mixing_RapAsym_2->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk 
				  hsldGenJet_Trk_sube0_Mixing_RapAsym_2->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			      else
				{
				  // leading jet -trk
				  hldGenJet_Trk_sube1_Mixing_RapAsym_2->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk
				  hsldGenJet_Trk_sube1_Mixing_RapAsym_2->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			    }
			}
		      else if(Deta_ldsldJet > 1.0) // higher rapidity
			{
			  gn_ldjtTRkCorr_RapAsym_3++;
			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0) {gn_sube0_ldjtTRkCorr_RapAsym_3++;}
			      else {gn_sube1_ldjtTRkCorr_RapAsym_3++;}
			    }

			  // leading jet -trk
			  hldGenJet_Trk_Mixing_RapAsym_3->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldGenJet_Trk_Mixing_RapAsym_3->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
			  
			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0)
				{
				  // subleading jet -trk
				  hldGenJet_Trk_sube0_Mixing_RapAsym_3->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk 
				  hsldGenJet_Trk_sube0_Mixing_RapAsym_3->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			      else
				{
				  // leading jet -trk
				  hldGenJet_Trk_sube1_Mixing_RapAsym_3->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk
				  hsldGenJet_Trk_sube1_Mixing_RapAsym_3->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			    }
			}
		      if((ldJet_eta)*(sldJet_eta) > 0) // same hemisphere
			{
			  gn_ldjtTRkCorr_RapAsym_sh++;
			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0) {gn_sube0_ldjtTRkCorr_RapAsym_sh++;}
			      else {gn_sube1_ldjtTRkCorr_RapAsym_sh++;}
			    }

			  // leading jet -trk
			  hldGenJet_Trk_Mixing_RapAsym_sh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldGenJet_Trk_Mixing_RapAsym_sh->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));

			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0)
				{
				  // subleading jet -trk
				  hldGenJet_Trk_sube0_Mixing_RapAsym_sh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk 
				  hsldGenJet_Trk_sube0_Mixing_RapAsym_sh->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			      else
				{
				  // leading jet -trk
				  hldGenJet_Trk_sube1_Mixing_RapAsym_sh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk
				  hsldGenJet_Trk_sube1_Mixing_RapAsym_sh->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			    }
			}
		      else if((ldJet_eta)*(sldJet_eta) < 0) // opposite hemisphere
			{
			  gn_ldjtTRkCorr_RapAsym_oh++;
			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0) {gn_sube0_ldjtTRkCorr_RapAsym_oh++;}
			      else {gn_sube1_ldjtTRkCorr_RapAsym_oh++;}
			    }
			  
			  // leading jet -trk
			  hldGenJet_Trk_Mixing_RapAsym_oh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldGenJet_Trk_Mixing_RapAsym_oh->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
			  
			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0)
				{
				  // subleading jet -trk
				  hldGenJet_Trk_sube0_Mixing_RapAsym_oh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk 
				  hsldGenJet_Trk_sube0_Mixing_RapAsym_oh->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			      else
				{
				  // leading jet -trk
				  hldGenJet_Trk_sube1_Mixing_RapAsym_oh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk
				  hsldGenJet_Trk_sube1_Mixing_RapAsym_oh->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			    }
			}
		    } // if(ldjetEta > sldjetEta)
		} // for gen
	    }// Track loop end
	  /*
	  // ntrk sum
	  if(colliding_system == "PbPb" && !isrc && do_sube)
	    {
	      if(ctbin == 0) hntrk_Mixing_0->Fill(TotalNtrk);
	      else if(ctbin == 1) hntrk_Mixing_1->Fill(TotalNtrk);
	      else if(ctbin == 2) hntrk_Mixing_2->Fill(TotalNtrk);
	      else if(ctbin == 3) hntrk_Mixing_3->Fill(TotalNtrk);
	      else if(ctbin == 4) hntrk_Mixing_4->Fill(TotalNtrk);
	    }
	  */

	  if(isrc)
	    {
	      if(ctbin == 0) hntrk_Mixing_0->Fill(TotalNtrk, evtw);
	      else if(ctbin == 1) hntrk_Mixing_1->Fill(TotalNtrk, evtw);
	      else if(ctbin == 2) hntrk_Mixing_2->Fill(TotalNtrk, evtw);
	      else if(ctbin == 3) hntrk_Mixing_3->Fill(TotalNtrk, evtw);
	      else if(ctbin == 4) hntrk_Mixing_4->Fill(TotalNtrk, evtw);
	    }
	      
	  // for jet pair
	  double jet_pair[3] = {0, (double)ldJetflvor, (double)ctbin};
	  if(isrc)
	    {
	      if(ldjtTRkCorr > 0){hldsld_Jet_pair_Mixing->Fill(jet_pair, (evtw*ldJetW));}
	      if(do_sube_rcjetgntrk) // for sube
		{
		  if(rcJet_gnTrk_sube0_ldjtTRkCorr > 0){hldsld_Jet_pair_sube0_Mixing->Fill(jet_pair, (evtw*ldJetW));}
		  if(rcJet_gnTrk_sube1_ldjtTRkCorr > 0){hldsld_Jet_pair_sube1_Mixing->Fill(jet_pair, (evtw*ldJetW));}
		}
	      if(ldJet_eta > sldJet_eta) // crucial for rapidity asymmetry
		{
		  if(ldjtTRkCorr_RapAsym > 0){hldsld_Jet_pair_Mixing_RapAsym->Fill(jet_pair, (evtw*ldJetW));}
		  if(do_sube_rcjetgntrk) // for sube
		    {
		      if(rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym > 0){hldsld_Jet_pair_sube0_Mixing_RapAsym->Fill(jet_pair, (evtw*ldJetW));}
		      if(rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym > 0){hldsld_Jet_pair_sube1_Mixing_RapAsym->Fill(jet_pair, (evtw*ldJetW));}
		    }
		  if(Deta_ldsldJet < 0.5) // mid rapidity
		    {
		      if(ldjtTRkCorr_RapAsym_1 > 0){hldsld_Jet_pair_Mixing_RapAsym_1->Fill(jet_pair, (evtw*ldJetW));}
		      if(do_sube_rcjetgntrk) // for sube
			{
			  if(rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym > 0){hldsld_Jet_pair_sube0_Mixing_RapAsym_1->Fill(jet_pair, (evtw*ldJetW));}
			  if(rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym > 0){hldsld_Jet_pair_sube1_Mixing_RapAsym_1->Fill(jet_pair, (evtw*ldJetW));}
			}
		    }
		  else if(Deta_ldsldJet > 0.5 && Deta_ldsldJet < 1.0) // intermediate rapidity
		    {
		      if(ldjtTRkCorr_RapAsym_2 > 0){hldsld_Jet_pair_Mixing_RapAsym_2->Fill(jet_pair, (evtw*ldJetW));}
		      if(do_sube_rcjetgntrk) // for sube
			{
			  if(rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym > 0){hldsld_Jet_pair_sube0_Mixing_RapAsym_2->Fill(jet_pair, (evtw*ldJetW));}
			  if(rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym > 0){hldsld_Jet_pair_sube1_Mixing_RapAsym_2->Fill(jet_pair, (evtw*ldJetW));}
			}
		    }
		  else if(Deta_ldsldJet > 1.0) // higher rapidity
		    {
		      if(ldjtTRkCorr_RapAsym_3 > 0){hldsld_Jet_pair_Mixing_RapAsym_3->Fill(jet_pair, (evtw*ldJetW));}
		      if(do_sube_rcjetgntrk) // for sube
			{
			  if(rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym_3 > 0){hldsld_Jet_pair_sube0_Mixing_RapAsym_3->Fill(jet_pair, (evtw*ldJetW));}
			  if(rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym_3 > 0){hldsld_Jet_pair_sube1_Mixing_RapAsym_3->Fill(jet_pair, (evtw*ldJetW));}
			}
		    }
		  if((ldJet_eta)*(sldJet_eta) > 0) // same hemisphere
		    {
		      if(ldjtTRkCorr_RapAsym_sh > 0){hldsld_Jet_pair_Mixing_RapAsym_sh->Fill(jet_pair, (evtw*ldJetW));}
		      if(do_sube_rcjetgntrk) // for sube
			{
			  if(rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym > 0){hldsld_Jet_pair_sube0_Mixing_RapAsym_sh->Fill(jet_pair, (evtw*ldJetW));}
			  if(rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym > 0){hldsld_Jet_pair_sube1_Mixing_RapAsym_sh->Fill(jet_pair, (evtw*ldJetW));}
			}
		    }
		  else if((ldJet_eta)*(sldJet_eta) < 0) // opposite hemisphere
		    {
		      if(ldjtTRkCorr_RapAsym_oh > 0){hldsld_Jet_pair_Mixing_RapAsym_oh->Fill(jet_pair, (evtw*ldJetW));}
		      if(do_sube_rcjetgntrk) // for sube
			{
			  if(rcJet_gnTrk_sube0_ldjtTRkCorr_RapAsym > 0){hldsld_Jet_pair_sube0_Mixing_RapAsym_oh->Fill(jet_pair, (evtw*ldJetW));}
			  if(rcJet_gnTrk_sube1_ldjtTRkCorr_RapAsym > 0){hldsld_Jet_pair_sube1_Mixing_RapAsym_oh->Fill(jet_pair, (evtw*ldJetW));}
			}
		    }
		}
	    }
	  else
	    {
	      if(gn_ldjtTRkCorr > 0){hldsld_GenJet_pair_Mixing->Fill(jet_pair, (evtw*ldJetW));}
	      if(do_sube) // for sube
		{
		  if(gn_sube0_ldjtTRkCorr > 0){hldsld_GenJet_pair_sube0_Mixing->Fill(jet_pair, (evtw*ldJetW));}
		  if(gn_sube1_ldjtTRkCorr > 0){hldsld_GenJet_pair_sube1_Mixing->Fill(jet_pair, (evtw*ldJetW));} 
		} 
	      if(ldJet_eta > sldJet_eta) // crucial for rapidity asymmetry
		{
		  if(gn_ldjtTRkCorr_RapAsym > 0){hldsld_GenJet_pair_Mixing_RapAsym->Fill(jet_pair, (evtw*ldJetW));}
		  if(do_sube) // for sube
		    {
		      if(gn_sube0_ldjtTRkCorr_RapAsym > 0){hldsld_GenJet_pair_sube0_Mixing_RapAsym->Fill(jet_pair, (evtw*ldJetW));}
		      if(gn_sube1_ldjtTRkCorr_RapAsym > 0){hldsld_GenJet_pair_sube1_Mixing_RapAsym->Fill(jet_pair, (evtw*ldJetW));} 
		    }
		  
		  if(Deta_ldsldJet < 0.5) // mid rapidity
		    {
		      if(gn_ldjtTRkCorr_RapAsym_1 > 0){hldsld_GenJet_pair_Mixing_RapAsym_1->Fill(jet_pair, (evtw*ldJetW));}
		      if(do_sube) // for sube
			{
			  if(gn_sube0_ldjtTRkCorr_RapAsym_1 > 0){hldsld_GenJet_pair_sube0_Mixing_RapAsym_1->Fill(jet_pair, (evtw*ldJetW));}
			  if(gn_sube1_ldjtTRkCorr_RapAsym_1 > 0){hldsld_GenJet_pair_sube1_Mixing_RapAsym_1->Fill(jet_pair, (evtw*ldJetW));} 
			}
		    }
		  else if(Deta_ldsldJet > 0.5 && Deta_ldsldJet < 1.0) // intermediate rapidity
		    {
		      if(gn_ldjtTRkCorr_RapAsym_2 > 0){hldsld_GenJet_pair_Mixing_RapAsym_2->Fill(jet_pair, (evtw*ldJetW));}
		      if(do_sube) // for sube
			{
			  if(gn_sube0_ldjtTRkCorr_RapAsym_2 > 0){hldsld_GenJet_pair_sube0_Mixing_RapAsym_2->Fill(jet_pair, (evtw*ldJetW));}
			  if(gn_sube1_ldjtTRkCorr_RapAsym_2 > 0){hldsld_GenJet_pair_sube1_Mixing_RapAsym_2->Fill(jet_pair, (evtw*ldJetW));} 
			}
		    }
		  else if(Deta_ldsldJet > 1.0) // higher rapidity
		    {
		      if(gn_ldjtTRkCorr_RapAsym_3 > 0){hldsld_GenJet_pair_Mixing_RapAsym_3->Fill(jet_pair, (evtw*ldJetW));}
		      if(do_sube) // for sube
			{
			  if(gn_sube0_ldjtTRkCorr_RapAsym_3 > 0){hldsld_GenJet_pair_sube0_Mixing_RapAsym_3->Fill(jet_pair, (evtw*ldJetW));}
			  if(gn_sube1_ldjtTRkCorr_RapAsym_3 > 0){hldsld_GenJet_pair_sube1_Mixing_RapAsym_3->Fill(jet_pair, (evtw*ldJetW));} 
			}
		    }
		  if((ldJet_eta)*(sldJet_eta) > 0) // same hemisphere
		    {
		      if(gn_ldjtTRkCorr_RapAsym_sh > 0){hldsld_GenJet_pair_Mixing_RapAsym_sh->Fill(jet_pair, (evtw*ldJetW));}
		      if(do_sube) // for sube
			{
			  if(gn_sube0_ldjtTRkCorr_RapAsym_sh > 0){hldsld_GenJet_pair_sube0_Mixing_RapAsym_sh->Fill(jet_pair, (evtw*ldJetW));}
			  if(gn_sube1_ldjtTRkCorr_RapAsym_sh > 0){hldsld_GenJet_pair_sube1_Mixing_RapAsym_sh->Fill(jet_pair, (evtw*ldJetW));} 
			}
		    }
		  else if((ldJet_eta)*(sldJet_eta) < 0) // opposite hemisphere
		    {
		      if(gn_ldjtTRkCorr_RapAsym_oh > 0){hldsld_GenJet_pair_Mixing_RapAsym_oh->Fill(jet_pair, (evtw*ldJetW));}
		      if(do_sube) // for sube
			{
			  if(gn_sube0_ldjtTRkCorr_RapAsym_oh > 0){hldsld_GenJet_pair_sube0_Mixing_RapAsym_oh->Fill(jet_pair, (evtw*ldJetW));}
			  if(gn_sube1_ldjtTRkCorr_RapAsym_oh > 0){hldsld_GenJet_pair_sube1_Mixing_RapAsym_oh->Fill(jet_pair, (evtw*ldJetW));} 
			}
		    }
		}
	    } // for gen
	} // 2nd event loop end
    } // event/leading subleadin jet loop end
  std::cout<<endl;
  std::cout<<"TotaldijetEvent in Mxing loop is: "<<TotaldijetEvent<<std::endl;
  std::cout<<endl;
  if(isrc) {std::cout<<"~~~~~end leading subleading jet tracks mixing correlation~~~~~~~~"<<std::endl;}
  else {std::cout<<"~~~~~end gen leading subleading jet tracks mixing correlation~~~~~~~~"<<std::endl;}
  std::cout<<endl;
}// function loop


void Jet_Track_mixing_withMB_corr_ldsld(const TString& colliding_system, const std::vector<double>& Evtw_vec_1D, const std::vector<int>& HiBin_vec_1D, const std::vector<int>& HiBinValue_vec_1D, const std::vector<double>& Vertexz_vec_1D, const std::vector<Long64_t>& Evtno_vec_1D, const std::vector<TVector3>& Filtered_ldjet_vec_1D, const std::vector<int>& Filtered_ldrefpartonB_vec_1D, const std::vector<double>& Filtered_ldJetW_vec_1D, const std::vector<TVector3>& Filtered_sldjet_vec_1D, const std::vector<int>& Filtered_sldrefpartonB_vec_1D, const std::vector<double>& Filtered_sldJetW_vec_1D, const std::vector<std::vector<TVector3>>& Filtered_Trk_pT_vec_2D, const std::vector<std::vector<double>>& Filtered_TrkW_vec_2D, const std::vector<std::vector<int>>& Filtered_TrkCharge_vec_2D, const std::vector<std::vector<int>>& Filtered_TrkSube_vec_2D, const std::vector<int>& HiBin_vec_MB_1D, const std::vector<int>& HiBinValue_vec_MB_1D, const std::vector<double>& Vertexz_vec_MB_1D, const std::vector<Long64_t>& Evtno_vec_MB_1D, const std::vector<int>& EvtCount_vec_MB_1D, const std::vector<std::vector<TVector3>>& Filtered_Trk_pT_vec_MB_2D, const std::vector<std::vector<double>>& Filtered_TrkW_vec_MB_2D, const std::vector<std::vector<int>>& Filtered_TrkSube_vec_MB_2D, const bool& isrc, const bool& do_sube, const bool& do_sube_rcjetgntrk)
{
  std::cout<<endl;

  if(isrc) {std::cout<<"~~~~~start leading subleading jet tracks mixing correlation with MB events~~~~~~~~"<<std::endl;}
  else {std::cout<<"~~~~~start gen leading subleading jet tracks mixing correlation with MB events~~~~~~~~"<<std::endl;}
  std::cout<<endl;

  if(Evtno_vec_MB_1D.size() != Filtered_Trk_pT_vec_MB_2D.size()){std::cout<<"MB events selcted for event and track vectors are not matching; please check"<<std::endl;}
  std::cout<<"Total selected MB events are: "<<Evtno_vec_MB_1D.size()<<"  "<<EvtCount_vec_MB_1D.size()<<"  "<<Filtered_Trk_pT_vec_MB_2D.size()<<std::endl;
  std::cout<<endl;
  
  int TotaldijetEvent = 0;
  
  for(int ievt = 0; ievt < Vertexz_vec_1D.size(); ievt++)
    {
      TotaldijetEvent++;
      
      if(ievt%10000 == 0) std::cout<<ievt<<"  events running for mixing correlation of total events "<< Vertexz_vec_1D.size() <<std::endl;
      
      double evtw = Evtw_vec_1D[ievt];
      //double evtw_i = Evtw_vec_1D[ievt];
      int ctbin = HiBin_vec_1D[ievt];
      int centbin = HiBinValue_vec_1D[ievt];
      int ctbinn = hCentBin->FindBin(centbin) -1;
      double vz = Vertexz_vec_1D[ievt];
      Long64_t evtno = Evtno_vec_1D[ievt];
      
      if(colliding_system == "PbPb")
	{
	  if(ctbin != ctbinn) {std::cout<<"ctbin and ctbinn are not same in PbPb, please check"<<std::endl;}
	}
      
      TVector3 ldjet_vec = Filtered_ldjet_vec_1D[ievt];
      double ldJetW = Filtered_ldJetW_vec_1D[ievt];
      int ldJetflvor = Filtered_ldrefpartonB_vec_1D[ievt];
      double ldJet_pt = ldjet_vec.Pt();
      double ldJet_eta = ldjet_vec.Eta();
      double ldJet_phi = ldjet_vec.Phi();
      
      TVector3 sldjet_vec = Filtered_sldjet_vec_1D[ievt];
      double sldJetW = Filtered_sldJetW_vec_1D[ievt];
      int sldJetflvor = Filtered_sldrefpartonB_vec_1D[ievt];
      double sldJet_pt = sldjet_vec.Pt();
      double sldJet_eta = sldjet_vec.Eta();
      double sldJet_phi = sldjet_vec.Phi();

      int ntrk = Filtered_Trk_pT_vec_2D[ievt].size();
      
      double Deta_ldsldJet = fabs(ldJet_eta - sldJet_eta);

      /*
      // mixing algorithm
      int mixstart = ievt+1;
      int mixend   = (int)Vertexz_vec_1D.size();
      
      if(mixstart > (0.5*(Vertexz_vec_1D.size())))
	{
	  mixstart = 0;
	  //mixend   = (int)Vertexz_vec_1D.size();
	  mixend = ievt -1;
	}
      
      int nmix = 0;
      */
      
      std::srand(std::time(0));
      std::random_device rd;
      std::mt19937 g(rd());

      std::vector<int> SelectMBEventCountIDForMixing;

      // Shuffle File 2 events once before matching
      std::vector<int> shuffled_EvtCount_vec_MB_1D = EvtCount_vec_MB_1D;  
      std::shuffle(shuffled_EvtCount_vec_MB_1D.begin(), shuffled_EvtCount_vec_MB_1D.end(), g);
      
      for(int ievt_j = 0; ievt_j < shuffled_EvtCount_vec_MB_1D.size(); ievt_j++) // loop over MV events to select desire mix events
	{
	  if(SelectMBEventCountIDForMixing.size() >= bkgFactor) break;  // Stop early if enough events are found
	  
	  int EventCountID = shuffled_EvtCount_vec_MB_1D[ievt_j];
	  
	  int ctbin_j = HiBin_vec_MB_1D[EventCountID];
          int centbin_j = HiBinValue_vec_MB_1D[EventCountID];
          double vz_j = Vertexz_vec_MB_1D[EventCountID];
	  int ntrk_j = Filtered_Trk_pT_vec_MB_2D[EventCountID].size();
	  Long64_t evtno_j = Evtno_vec_MB_1D[EventCountID];
	  
	  //if((fabs(vz - vz_j) < Dvz_cut) && (ctbin == ctbin_j) && (fabs(centbin - centbin_j) <= Dhibin_cut) && (evtno != evtno_j) && (fabs(ntrk - ntrk_j) <= Dntrk_cut))
	  if((fabs(vz - vz_j) <= Dvz_cut) && (ctbin == ctbin_j) && (fabs(centbin - centbin_j) <= Dhibin_cut) && (evtno != evtno_j))
	    {
	      SelectMBEventCountIDForMixing.push_back(EventCountID); // push_back the selected event counts for mixing
	    }
	}

      if(ievt%10000 == 0)
	{
	  std::cout<<"number of MB mixing events found for "<<ievt<<"th dijet event is: "<<SelectMBEventCountIDForMixing.size()<<std::endl;
	  //std::cout<<"Selected events ID are: ";
	  //for(int jevt = 0; jevt < SelectMBEventCountIDForMixing.size(); jevt++){std::cout<<SelectMBEventCountIDForMixing[jevt]<<", "; }
	}

      //if(ievt%500 == 0){std::cout<<SelectMBEventCountIDForMixing.size()<<" events found for mixing for the ith event:"<<ievt<<std::endl;}
      
      for(int jevt = 0; jevt < SelectMBEventCountIDForMixing.size(); jevt++) // 2nd event loop
	{
	  int MixEvtID = SelectMBEventCountIDForMixing[jevt];
	  
	  //if(ievt%200 == 0 && ctbin == 4){std::cout<<ievt<<"  "<<MixEvtID<<"  "<<ctbin<<"  "<<HiBin_vec_MB_1D[MixEvtID]<<"  "<<centbin<<"  "<<HiBinValue_vec_MB_1D[MixEvtID]<<"  "<<vz<<"  "<<Vertexz_vec_MB_1D[MixEvtID]<<"  "<<evtno<<"  "<<Evtno_vec_MB_1D[MixEvtID]<<"  "<<ntrk<<"  "<<Filtered_Trk_pT_vec_MB_2D[MixEvtID].size()<<std::endl;}

	  
	  int ldjtTRkCorr = 0, ldjtTRkCorr_RapAsym = 0, ldjtTRkCorr_RapAsym_1 = 0, ldjtTRkCorr_RapAsym_2 = 0, ldjtTRkCorr_RapAsym_3 = 0, ldjtTRkCorr_RapAsym_sh = 0, ldjtTRkCorr_RapAsym_oh = 0;
	  int gn_ldjtTRkCorr = 0, gn_ldjtTRkCorr_RapAsym = 0, gn_ldjtTRkCorr_RapAsym_1 = 0, gn_ldjtTRkCorr_RapAsym_2 = 0, gn_ldjtTRkCorr_RapAsym_3 = 0, gn_ldjtTRkCorr_RapAsym_sh = 0, gn_ldjtTRkCorr_RapAsym_oh = 0;
	  int gn_sube0_ldjtTRkCorr = 0, gn_sube0_ldjtTRkCorr_RapAsym = 0, gn_sube0_ldjtTRkCorr_RapAsym_1 = 0, gn_sube0_ldjtTRkCorr_RapAsym_2 = 0, gn_sube0_ldjtTRkCorr_RapAsym_3 = 0, gn_sube0_ldjtTRkCorr_RapAsym_sh = 0, gn_sube0_ldjtTRkCorr_RapAsym_oh = 0;
	  int gn_sube1_ldjtTRkCorr = 0, gn_sube1_ldjtTRkCorr_RapAsym = 0, gn_sube1_ldjtTRkCorr_RapAsym_1 = 0, gn_sube1_ldjtTRkCorr_RapAsym_2 = 0, gn_sube1_ldjtTRkCorr_RapAsym_3 = 0, gn_sube1_ldjtTRkCorr_RapAsym_sh = 0, gn_sube1_ldjtTRkCorr_RapAsym_oh = 0;

	  //double TotalNtrk = 0.;

	  for(int itrk = 0; itrk < Filtered_Trk_pT_vec_MB_2D[MixEvtID].size(); itrk++)
	    {
	      TVector3 trk_vec = Filtered_Trk_pT_vec_MB_2D[MixEvtID][itrk];
	      double trk_w = Filtered_TrkW_vec_MB_2D[MixEvtID][itrk];
	      double trk_pt = trk_vec.Pt();
	      double trk_eta = trk_vec.Eta();
	      double trk_phi = trk_vec.Phi();
	      int trk_sube = Filtered_TrkSube_vec_MB_2D[MixEvtID][itrk];
	      int trk_ptbin = hTrkpTBin->FindBin(trk_pt) - 1;

	      double Trk_pT_Eta_Phi_ctbin[4] = {trk_pt, trk_eta, trk_phi, (double)ctbin};
	      hTrk_Mixing_GenpT_Eta_Phi_ctbin_pTCut_W->Fill(Trk_pT_Eta_Phi_ctbin, trk_w);

	      //if(colliding_system == "PbPb" && !isrc && trk_pt >= 1 && trk_pt < 2) TotalNtrk += trk_w;
	      if(!isrc && do_sube)
		{
		  double Trk_pT_Eta_Phi_ctbin[4] = {trk_pt, trk_eta, trk_phi, (double)ctbin};
		  if(trk_sube == 0){hTrk_sube0_Mixing_GenpT_Eta_Phi_ctbin_pTCut_W->Fill(Trk_pT_Eta_Phi_ctbin, trk_w);}
		  else {hTrk_sube1_Mixing_GenpT_Eta_Phi_ctbin_pTCut_W->Fill(Trk_pT_Eta_Phi_ctbin, trk_w);}
		}
	      
	      // for leading jet - trk
	      double Deta_ldjet_trk = trk_eta - ldJet_eta;
	      double Dphi_ldjet_trk = trk_phi - ldJet_phi;
	      
	      if(Dphi_ldjet_trk > 1.5*TMath::Pi())
		{
		  Dphi_ldjet_trk = Dphi_ldjet_trk - 2.0*TMath::Pi();
		}
	      else if(Dphi_ldjet_trk < -0.5*TMath::Pi())
		{
		  Dphi_ldjet_trk = Dphi_ldjet_trk + 2.0*TMath::Pi();
		}
	      
	      // for subleading jet - trk
	      double Deta_sldjet_trk = trk_eta - sldJet_eta;
	      double Dphi_sldjet_trk = trk_phi - sldJet_phi;
	      
	      if(Dphi_sldjet_trk > 1.5*TMath::Pi())
		{
		  Dphi_sldjet_trk = Dphi_sldjet_trk - 2.0*TMath::Pi();
		}
	      else if(Dphi_sldjet_trk < -0.5*TMath::Pi())
		{
		  Dphi_sldjet_trk = Dphi_sldjet_trk + 2.0*TMath::Pi();
		}
	      
	      // leading jet -trk
	      double ldjet_trk_mixing[5] = {Deta_ldjet_trk, Dphi_ldjet_trk, (double)trk_ptbin, (double)ldJetflvor, (double)ctbin};
	      
	      // subleading jet -trk
	      double sldjet_trk_mixing[5] = {Deta_sldjet_trk, Dphi_sldjet_trk, (double)trk_ptbin, (double)sldJetflvor, (double)ctbin};
	      
	      if(isrc) // for reco/data
		{
		  ldjtTRkCorr++;

		  // leading jet -trk
		  hldJet_Trk_Mixing->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
		  
		  // subleading jet -trk
		  hsldJet_Trk_Mixing->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));

		  if(ldJet_eta > sldJet_eta) // crucial for rapidity asymmetry
		    {
		      ldjtTRkCorr_RapAsym++;
		      		      
		      // leading jet -trk
		      hldJet_Trk_Mixing_RapAsym->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));

		      // subleading jet -trk
		      hsldJet_Trk_Mixing_RapAsym->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));

		      if(Deta_ldsldJet < 0.5) // mid rapidity
			{
			  ldjtTRkCorr_RapAsym_1++;
			  
			  // leading jet -trk
			  hldJet_Trk_Mixing_RapAsym_1->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldJet_Trk_Mixing_RapAsym_1->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
			}
		      else if(Deta_ldsldJet > 0.5 && Deta_ldsldJet < 1.0) // intermediate rapidity
			{
			  ldjtTRkCorr_RapAsym_2++;
									  
			  // leading jet -trk
			  hldJet_Trk_Mixing_RapAsym_2->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldJet_Trk_Mixing_RapAsym_2->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
			}
		      else if(Deta_ldsldJet > 1.0) // higher rapidity
			{
			  ldjtTRkCorr_RapAsym_3++;
						  
			  // leading jet -trk
			  hldJet_Trk_Mixing_RapAsym_3->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldJet_Trk_Mixing_RapAsym_3->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
			}
		      if((ldJet_eta)*(sldJet_eta) > 0) // same hemisphere 
			{
			  ldjtTRkCorr_RapAsym_sh++;
						  
			  // leading jet -trk
			  hldJet_Trk_Mixing_RapAsym_sh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldJet_Trk_Mixing_RapAsym_sh->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
			  
			}
		      else if((ldJet_eta)*(sldJet_eta) < 0) // opposite hemisphere
			{
			  ldjtTRkCorr_RapAsym_oh++;
			  
			  // leading jet -trk
			  hldJet_Trk_Mixing_RapAsym_oh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldJet_Trk_Mixing_RapAsym_oh->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
			}
		    }
		} // for reco/data
	      else
		{
		  gn_ldjtTRkCorr++;
		  if(do_sube) // for sube
		    {
		      if(trk_sube == 0) {gn_sube0_ldjtTRkCorr++;}
		      else {gn_sube1_ldjtTRkCorr++;}
		    }
		  
		  // leading jet -trk
		  hldGenJet_Trk_Mixing->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
		  
		  // subleading jet -trk
		  hsldGenJet_Trk_Mixing->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));

		  if(do_sube) // for sube
		    {
		      if(trk_sube == 0)
			{
			  // subleading jet -trk
			  hldGenJet_Trk_sube0_Mixing->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  // subleading jet -trk 
			  hsldGenJet_Trk_sube0_Mixing->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
			}
		      else
			{
			  // leading jet -trk
			  hldGenJet_Trk_sube1_Mixing->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  // subleading jet -trk
			  hsldGenJet_Trk_sube1_Mixing->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
			}
		    }
		  
		  if(ldJet_eta > sldJet_eta) // crucial for rapidity asymmetry
		    {
		      gn_ldjtTRkCorr_RapAsym++;
		      if(do_sube) // for sube
			{
			  if(trk_sube == 0) {gn_sube0_ldjtTRkCorr_RapAsym++;}
			  else {gn_sube1_ldjtTRkCorr_RapAsym++;}
			}

		      // leading jet -trk
		      hldGenJet_Trk_Mixing_RapAsym->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
		      
		      // subleading jet -trk
		      hsldGenJet_Trk_Mixing_RapAsym->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));

		      if(do_sube) // for sube
			{
			  if(trk_sube == 0)
			    {
			      // subleading jet -trk
			      hldGenJet_Trk_sube0_Mixing_RapAsym->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			      // subleading jet -trk 
			      hsldGenJet_Trk_sube0_Mixing_RapAsym->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
			    }
			  else
			    {
			      // leading jet -trk
			      hldGenJet_Trk_sube1_Mixing_RapAsym->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			      // subleading jet -trk
			      hsldGenJet_Trk_sube1_Mixing_RapAsym->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
			    }
			}
		      
		      if(Deta_ldsldJet < 0.5) // mid rapidity
			{
			  gn_ldjtTRkCorr_RapAsym_1++;
			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0) {gn_sube0_ldjtTRkCorr_RapAsym_1++;}
			      else {gn_sube1_ldjtTRkCorr_RapAsym_1++;}
			    }
			  
			  // leading jet -trk
			  hldGenJet_Trk_Mixing_RapAsym_1->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldGenJet_Trk_Mixing_RapAsym_1->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));

			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0)
				{
				  // subleading jet -trk
				  hldGenJet_Trk_sube0_Mixing_RapAsym_1->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk 
				  hsldGenJet_Trk_sube0_Mixing_RapAsym_1->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			      else
				{
				  // leading jet -trk
				  hldGenJet_Trk_sube1_Mixing_RapAsym_1->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk
				  hsldGenJet_Trk_sube1_Mixing_RapAsym_1->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			    }
			}
		      else if(Deta_ldsldJet > 0.5 && Deta_ldsldJet < 1.0) // intermediate rapidity
			{
			  gn_ldjtTRkCorr_RapAsym_2++;
			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0) {gn_sube0_ldjtTRkCorr_RapAsym_2++;}
			      else {gn_sube1_ldjtTRkCorr_RapAsym_2++;}
			    }
  
			  // leading jet -trk
			  hldGenJet_Trk_Mixing_RapAsym_2->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldGenJet_Trk_Mixing_RapAsym_2->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));

			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0)
				{
				  // subleading jet -trk
				  hldGenJet_Trk_sube0_Mixing_RapAsym_2->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk 
				  hsldGenJet_Trk_sube0_Mixing_RapAsym_2->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			      else
				{
				  // leading jet -trk
				  hldGenJet_Trk_sube1_Mixing_RapAsym_2->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk
				  hsldGenJet_Trk_sube1_Mixing_RapAsym_2->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			    }
			}
		      else if(Deta_ldsldJet > 1.0) // higher rapidity
			{
			  gn_ldjtTRkCorr_RapAsym_3++;
			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0) {gn_sube0_ldjtTRkCorr_RapAsym_3++;}
			      else {gn_sube1_ldjtTRkCorr_RapAsym_3++;}
			    }

			  // leading jet -trk
			  hldGenJet_Trk_Mixing_RapAsym_3->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldGenJet_Trk_Mixing_RapAsym_3->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
			  
			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0)
				{
				  // subleading jet -trk
				  hldGenJet_Trk_sube0_Mixing_RapAsym_3->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk 
				  hsldGenJet_Trk_sube0_Mixing_RapAsym_3->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			      else
				{
				  // leading jet -trk
				  hldGenJet_Trk_sube1_Mixing_RapAsym_3->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk
				  hsldGenJet_Trk_sube1_Mixing_RapAsym_3->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			    }
			}
		      if((ldJet_eta)*(sldJet_eta) > 0) // same hemisphere
			{
			  gn_ldjtTRkCorr_RapAsym_sh++;
			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0) {gn_sube0_ldjtTRkCorr_RapAsym_sh++;}
			      else {gn_sube1_ldjtTRkCorr_RapAsym_sh++;}
			    }

			  // leading jet -trk
			  hldGenJet_Trk_Mixing_RapAsym_sh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldGenJet_Trk_Mixing_RapAsym_sh->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));

			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0)
				{
				  // subleading jet -trk
				  hldGenJet_Trk_sube0_Mixing_RapAsym_sh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk 
				  hsldGenJet_Trk_sube0_Mixing_RapAsym_sh->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			      else
				{
				  // leading jet -trk
				  hldGenJet_Trk_sube1_Mixing_RapAsym_sh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk
				  hsldGenJet_Trk_sube1_Mixing_RapAsym_sh->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			    }
			}
		      else if((ldJet_eta)*(sldJet_eta) < 0) // opposite hemisphere
			{
			  gn_ldjtTRkCorr_RapAsym_oh++;
			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0) {gn_sube0_ldjtTRkCorr_RapAsym_oh++;}
			      else {gn_sube1_ldjtTRkCorr_RapAsym_oh++;}
			    }
			  
			  // leading jet -trk
			  hldGenJet_Trk_Mixing_RapAsym_oh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
			  
			  // subleading jet -trk
			  hsldGenJet_Trk_Mixing_RapAsym_oh->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
			  
			  if(do_sube) // for sube
			    {
			      if(trk_sube == 0)
				{
				  // subleading jet -trk
				  hldGenJet_Trk_sube0_Mixing_RapAsym_oh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk 
				  hsldGenJet_Trk_sube0_Mixing_RapAsym_oh->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			      else
				{
				  // leading jet -trk
				  hldGenJet_Trk_sube1_Mixing_RapAsym_oh->Fill(ldjet_trk_mixing, (evtw*ldJetW*trk_w));
				  // subleading jet -trk
				  hsldGenJet_Trk_sube1_Mixing_RapAsym_oh->Fill(sldjet_trk_mixing, (evtw*sldJetW*trk_w));
				}
			    }
			}
		    } // if(ldjetEta > sldjetEta)
		} // for gen
	    }// Track loop end
	  /*
	  //ntrk sum
	  if(colliding_system == "PbPb" && !isrc)
	    {
	      if(ctbin == 0) hntrk_Mixing_0->Fill(TotalNtrk);
              else if(ctbin == 1) hntrk_Mixing_1->Fill(TotalNtrk);
              else if(ctbin == 2) hntrk_Mixing_2->Fill(TotalNtrk);
              else if(ctbin == 3) hntrk_Mixing_3->Fill(TotalNtrk);
              else if(ctbin == 4) hntrk_Mixing_4->Fill(TotalNtrk);
	    }
	  */
	  
	  // for jet pair
	  double jet_pair[3] = {0, (double)ldJetflvor, (double)ctbin};
	  if(isrc)
	    {
	      if(ldjtTRkCorr > 0){hldsld_Jet_pair_Mixing->Fill(jet_pair, (evtw*ldJetW));}
	      
	      if(ldJet_eta > sldJet_eta) // crucial for rapidity asymmetry
		{
		  if(ldjtTRkCorr_RapAsym > 0){hldsld_Jet_pair_Mixing_RapAsym->Fill(jet_pair, (evtw*ldJetW));}
		  
		  if(Deta_ldsldJet < 0.5) // mid rapidity
		    {
		      if(ldjtTRkCorr_RapAsym_1 > 0){hldsld_Jet_pair_Mixing_RapAsym_1->Fill(jet_pair, (evtw*ldJetW));}
		    }
		  else if(Deta_ldsldJet > 0.5 && Deta_ldsldJet < 1.0) // intermediate rapidity
		    {
		      if(ldjtTRkCorr_RapAsym_2 > 0){hldsld_Jet_pair_Mixing_RapAsym_2->Fill(jet_pair, (evtw*ldJetW));}
		    }
		  else if(Deta_ldsldJet > 1.0) // higher rapidity
		    {
		      if(ldjtTRkCorr_RapAsym_3 > 0){hldsld_Jet_pair_Mixing_RapAsym_3->Fill(jet_pair, (evtw*ldJetW));}
		    }
		  if((ldJet_eta)*(sldJet_eta) > 0) // same hemisphere
		    {
		      if(ldjtTRkCorr_RapAsym_sh > 0){hldsld_Jet_pair_Mixing_RapAsym_sh->Fill(jet_pair, (evtw*ldJetW));}
		    }
		  else if((ldJet_eta)*(sldJet_eta) < 0) // opposite hemisphere
		    {
		      if(ldjtTRkCorr_RapAsym_oh > 0){hldsld_Jet_pair_Mixing_RapAsym_oh->Fill(jet_pair, (evtw*ldJetW));}
		    }
		}
	    }
	  else
	    {
	      if(gn_ldjtTRkCorr > 0){hldsld_GenJet_pair_Mixing->Fill(jet_pair, (evtw*ldJetW));}
	      if(do_sube) // for sube
		{
		  if(gn_sube0_ldjtTRkCorr > 0){hldsld_GenJet_pair_sube0_Mixing->Fill(jet_pair, (evtw*ldJetW));}
		  if(gn_sube1_ldjtTRkCorr > 0){hldsld_GenJet_pair_sube1_Mixing->Fill(jet_pair, (evtw*ldJetW));} 
		} 
	      if(ldJet_eta > sldJet_eta) // crucial for rapidity asymmetry
		{
		  if(gn_ldjtTRkCorr_RapAsym > 0){hldsld_GenJet_pair_Mixing_RapAsym->Fill(jet_pair, (evtw*ldJetW));}
		  if(do_sube) // for sube
		    {
		      if(gn_sube0_ldjtTRkCorr_RapAsym > 0){hldsld_GenJet_pair_sube0_Mixing_RapAsym->Fill(jet_pair, (evtw*ldJetW));}
		      if(gn_sube1_ldjtTRkCorr_RapAsym > 0){hldsld_GenJet_pair_sube1_Mixing_RapAsym->Fill(jet_pair, (evtw*ldJetW));} 
		    }
		  
		  if(Deta_ldsldJet < 0.5) // mid rapidity
		    {
		      if(gn_ldjtTRkCorr_RapAsym_1 > 0){hldsld_GenJet_pair_Mixing_RapAsym_1->Fill(jet_pair, (evtw*ldJetW));}
		      if(do_sube) // for sube
			{
			  if(gn_sube0_ldjtTRkCorr_RapAsym_1 > 0){hldsld_GenJet_pair_sube0_Mixing_RapAsym_1->Fill(jet_pair, (evtw*ldJetW));}
			  if(gn_sube1_ldjtTRkCorr_RapAsym_1 > 0){hldsld_GenJet_pair_sube1_Mixing_RapAsym_1->Fill(jet_pair, (evtw*ldJetW));} 
			}
		    }
		  else if(Deta_ldsldJet > 0.5 && Deta_ldsldJet < 1.0) // intermediate rapidity
		    {
		      if(gn_ldjtTRkCorr_RapAsym_2 > 0){hldsld_GenJet_pair_Mixing_RapAsym_2->Fill(jet_pair, (evtw*ldJetW));}
		      if(do_sube) // for sube
			{
			  if(gn_sube0_ldjtTRkCorr_RapAsym_2 > 0){hldsld_GenJet_pair_sube0_Mixing_RapAsym_2->Fill(jet_pair, (evtw*ldJetW));}
			  if(gn_sube1_ldjtTRkCorr_RapAsym_2 > 0){hldsld_GenJet_pair_sube1_Mixing_RapAsym_2->Fill(jet_pair, (evtw*ldJetW));} 
			}
		    }
		  else if(Deta_ldsldJet > 1.0) // higher rapidity
		    {
		      if(gn_ldjtTRkCorr_RapAsym_3 > 0){hldsld_GenJet_pair_Mixing_RapAsym_3->Fill(jet_pair, (evtw*ldJetW));}
		      if(do_sube) // for sube
			{
			  if(gn_sube0_ldjtTRkCorr_RapAsym_3 > 0){hldsld_GenJet_pair_sube0_Mixing_RapAsym_3->Fill(jet_pair, (evtw*ldJetW));}
			  if(gn_sube1_ldjtTRkCorr_RapAsym_3 > 0){hldsld_GenJet_pair_sube1_Mixing_RapAsym_3->Fill(jet_pair, (evtw*ldJetW));} 
			}
		    }
		  if((ldJet_eta)*(sldJet_eta) > 0) // same hemisphere
		    {
		      if(gn_ldjtTRkCorr_RapAsym_sh > 0){hldsld_GenJet_pair_Mixing_RapAsym_sh->Fill(jet_pair, (evtw*ldJetW));}
		      if(do_sube) // for sube
			{
			  if(gn_sube0_ldjtTRkCorr_RapAsym_sh > 0){hldsld_GenJet_pair_sube0_Mixing_RapAsym_sh->Fill(jet_pair, (evtw*ldJetW));}
			  if(gn_sube1_ldjtTRkCorr_RapAsym_sh > 0){hldsld_GenJet_pair_sube1_Mixing_RapAsym_sh->Fill(jet_pair, (evtw*ldJetW));} 
			}
		    }
		  else if((ldJet_eta)*(sldJet_eta) < 0) // opposite hemisphere
		    {
		      if(gn_ldjtTRkCorr_RapAsym_oh > 0){hldsld_GenJet_pair_Mixing_RapAsym_oh->Fill(jet_pair, (evtw*ldJetW));}
		      if(do_sube) // for sube
			{
			  if(gn_sube0_ldjtTRkCorr_RapAsym_oh > 0){hldsld_GenJet_pair_sube0_Mixing_RapAsym_oh->Fill(jet_pair, (evtw*ldJetW));}
			  if(gn_sube1_ldjtTRkCorr_RapAsym_oh > 0){hldsld_GenJet_pair_sube1_Mixing_RapAsym_oh->Fill(jet_pair, (evtw*ldJetW));} 
			}
		    }
		} // ld eta > sub ld eta
	    } // for gen
	} // 2nd event loop end
    } // event/leading subleadin jet loop end
  std::cout<<std::endl;
  std::cout<<"TotaldijetEvent in Mxing loop is: "<<TotaldijetEvent<<std::endl;
  std::cout<<endl;
  if(isrc) {std::cout<<"~~~~~end leading subleading jet tracks mixing correlation~~~~~~~~"<<std::endl;}
  else {std::cout<<"~~~~~end gen leading subleading jet tracks mixing correlation~~~~~~~~"<<std::endl;}
  std::cout<<endl;
}// function loop

// For centrality determination from Hihf
const Int_t nBins = 200; // table of bin edges

const Double_t binTable_Nom[nBins+1] = {0, 10.5072, 11.2099, 11.8364, 12.478, 13.1194, 13.7623, 14.4081, 15.0709, 15.7532, 16.4673, 17.1881, 17.923, 18.673, 19.4865, 20.3033, 21.1536, 22.0086, 22.9046, 23.8196, 24.7924, 25.8082, 26.8714, 27.9481, 29.0828, 30.2757, 31.5043, 32.8044, 34.1572, 35.6142, 37.1211, 38.6798, 40.3116, 42.0398, 43.8572, 45.6977, 47.6312, 49.6899, 51.815, 54.028, 56.3037, 58.7091, 61.2024, 63.8353, 66.5926, 69.3617, 72.2068, 75.2459, 78.3873, 81.5916, 84.9419, 88.498, 92.1789, 95.9582, 99.8431, 103.739, 107.78, 111.97, 116.312, 120.806, 125.46, 130.269, 135.247, 140.389, 145.713, 151.212, 156.871, 162.729, 168.762, 174.998, 181.424, 188.063, 194.907, 201.942, 209.19, 216.683, 224.37, 232.291, 240.43, 248.807, 257.416, 266.256, 275.348, 284.668, 294.216, 304.053, 314.142, 324.488, 335.101, 345.974, 357.116, 368.547, 380.283, 392.29, 404.564, 417.122, 429.968, 443.116, 456.577, 470.357, 484.422, 498.78, 513.473, 528.479, 543.813, 559.445, 575.411, 591.724, 608.352, 625.344, 642.686, 660.361, 678.371, 696.749, 715.485, 734.608, 754.068, 773.846, 794.046, 814.649, 835.608, 856.972, 878.719, 900.887, 923.409, 946.374, 969.674, 993.435, 1017.62, 1042.21, 1067.28, 1092.72, 1118.64, 1144.96, 1171.71, 1198.98, 1226.67, 1254.82, 1283.46, 1312.65, 1342.21, 1372.27, 1402.85, 1433.93, 1465.49, 1497.62, 1530.29, 1563.49, 1597.22, 1631.49, 1666.37, 1701.8, 1737.75, 1774.35, 1811.51, 1849.29, 1887.75, 1926.79, 1966.6, 2006.97, 2047.99, 2089.71, 2132.1, 2175.23, 2219.17, 2263.72, 2309.2, 2355.43, 2402.47, 2450.33, 2499.05, 2548.66, 2599.16, 2650.59, 2703.03, 2756.32, 2810.75, 2866.27, 2922.91, 2980.54, 3039.47, 3099.53, 3160.98, 3223.66, 3287.71, 3353.18, 3420.34, 3489.13, 3559.72, 3632.06, 3706.18, 3782.42, 3860.78, 3941.42, 4024.52, 4110.27, 4199.4, 4292.8, 4394.49, 4519.52, 5199.95};

const Double_t binTable_Down[nBins+1] = {0, 10.5071, 11.2094, 11.8357, 12.4763, 13.117, 13.7597, 14.4049, 15.0671, 15.7491, 16.4622, 17.1812, 17.9144, 18.6674, 19.4797, 20.2963, 21.1435, 21.9974, 22.8928, 23.8068, 24.7805, 25.7931, 26.8556, 27.9308, 29.0638, 30.2582, 31.4795, 32.7816, 34.1349, 35.5834, 37.0941, 38.6474, 40.2782, 42.0035, 43.8112, 45.6576, 47.5758, 49.6381, 51.6667, 53.7353, 55.8903, 58.1259, 60.4528, 62.8712, 65.3859, 67.9968, 70.7065, 73.5231, 76.4519, 79.4922, 82.6461, 85.9264, 89.3269, 92.8562, 96.5212, 100.322, 104.262, 108.344, 112.585, 116.971, 121.521, 126.225, 131.09, 136.127, 141.328, 146.721, 152.284, 158.014, 163.935, 170.054, 176.372, 182.878, 189.602, 196.532, 203.653, 211.017, 218.599, 226.387, 234.418, 242.667, 251.16, 259.886, 268.852, 278.071, 287.498, 297.2, 307.184, 317.409, 327.894, 338.66, 349.686, 360.996, 372.607, 384.508, 396.669, 409.133, 421.86, 434.906, 448.258, 461.916, 475.906, 490.16, 504.74, 519.663, 534.911, 550.453, 566.322, 582.525, 599.08, 615.968, 633.211, 650.805, 668.76, 687.048, 705.707, 724.774, 744.163, 763.9, 783.999, 804.528, 825.432, 846.746, 868.429, 890.523, 913.007, 935.952, 959.211, 982.919, 1007.08, 1031.63, 1056.62, 1082.08, 1107.96, 1134.24, 1160.99, 1188.22, 1215.91, 1244.06, 1272.69, 1301.85, 1331.45, 1361.51, 1392.07, 1423.18, 1454.77, 1486.93, 1519.57, 1552.81, 1586.55, 1620.87, 1655.79, 1691.26, 1727.27, 1763.93, 1801.12, 1838.97, 1877.47, 1916.61, 1956.45, 1996.89, 2038.04, 2079.84, 2122.35, 2165.52, 2209.53, 2254.24, 2299.83, 2346.19, 2393.31, 2441.28, 2490.16, 2539.86, 2590.57, 2642.16, 2694.74, 2748.23, 2802.81, 2858.47, 2915.33, 2973.2, 3032.28, 3092.56, 3154.24, 3217.19, 3281.45, 3347.18, 3414.6, 3483.65, 3554.56, 3627.2, 3701.66, 3778.25, 3856.97, 3937.98, 4021.48, 4107.62, 4197.21, 4291.05, 4393.19, 4518.6, 5199.95};

const Double_t binTable_Up[nBins+1] = {0, 10.5075, 11.2107, 11.838, 12.4797, 13.1213, 13.7641, 14.4124, 15.0745, 15.7577, 16.473, 17.1939, 17.9297, 18.6812, 19.4958, 20.3143, 21.1648, 22.0218, 22.9159, 23.8328, 24.8059, 25.8204, 26.89, 27.9702, 29.1042, 30.3022, 31.528, 32.8347, 34.1896, 35.6439, 37.1542, 38.7172, 40.3518, 42.091, 43.9053, 45.7415, 47.6853, 49.7457, 51.8755, 54.0983, 56.3594, 58.7848, 61.2861, 63.9228, 66.6825, 69.4421, 72.297, 75.3547, 78.4967, 81.6977, 85.0755, 88.6211, 92.3058, 96.1071, 99.9975, 104.065, 108.272, 112.512, 116.906, 121.601, 126.465, 131.482, 136.866, 142.229, 147.786, 153.546, 159.571, 165.586, 171.902, 178.419, 185.063, 191.856, 199.055, 206.261, 213.999, 221.719, 229.671, 237.84, 246.088, 254.828, 263.883, 272.907, 282.236, 291.925, 301.519, 311.477, 321.691, 332.153, 342.892, 353.878, 365.161, 376.742, 388.577, 400.684, 413.075, 425.746, 438.711, 451.989, 465.556, 479.45, 493.608, 508.077, 522.891, 538.003, 553.415, 569.151, 585.216, 601.601, 618.354, 635.422, 652.84, 670.599, 688.699, 707.161, 726.014, 745.185, 764.687, 784.557, 804.838, 825.489, 846.537, 867.951, 889.752, 911.955, 934.588, 957.52, 980.912, 1004.73, 1028.94, 1053.57, 1078.67, 1104.17, 1130.07, 1156.39, 1183.2, 1210.47, 1238.17, 1266.38, 1295.02, 1324.16, 1353.71, 1383.77, 1414.35, 1445.41, 1477, 1509.09, 1541.74, 1574.88, 1608.59, 1642.83, 1677.66, 1713.07, 1748.98, 1785.47, 1822.63, 1860.33, 1898.72, 1937.73, 1977.42, 2017.71, 2058.62, 2100.25, 2142.57, 2185.56, 2229.38, 2273.91, 2319.2, 2365.33, 2412.22, 2459.94, 2508.52, 2557.98, 2608.35, 2659.61, 2711.86, 2765, 2819.23, 2874.58, 2930.97, 2988.46, 3047.12, 3106.95, 3168.15, 3230.6, 3294.37, 3359.58, 3426.47, 3494.95, 3565.21, 3637.21, 3711.03, 3786.91, 3864.85, 3945.11, 4027.8, 4113.06, 4201.73, 4294.72, 4395.9, 4520.5, 5199.95};

Int_t getHiBinFromhiHF(const Double_t& hiHF, const Int_t& cent_cond)
{
  Int_t binPos = -1;
  for(int i = 0; i < nBins; ++i)
    {
      if(cent_cond == 0)
	{
	  if(hiHF >= binTable_Nom[i] && hiHF < binTable_Nom[i+1])
	    {
	      binPos = i;
	      break;
	    }
	}
      else if(cent_cond == 1)
	{
	  if(hiHF >= binTable_Up[i] && hiHF < binTable_Up[i+1])
	    {
	      binPos = i;
	      break;
	    }
	}
      else if(cent_cond == 2)
	{
	  if(hiHF >= binTable_Down[i] && hiHF < binTable_Down[i+1])
	    {
	      binPos = i;
	      break;
	    }
	}
    }
  
  binPos = nBins - 1 - binPos;
  return (Int_t)(200*((Double_t)binPos)/((Double_t)nBins));
}

// For JER scale factor 
const Int_t nBins_Eta = 5; // eta bins
// for PbPb taken from Autumn18_RunD_V7b_MC_SF_AK4PF.txt
const Double_t SF_Nom_PbPb[nBins_Eta] = {1.1742, 1.1930, 1.1451, 1.1618, 1.1455};
const Double_t SF_Down_PbPb[nBins_Eta] = {1.1415, 1.1559, 1.0812, 1.1086, 1.0838};
const Double_t SF_Up_PbPb[nBins_Eta] = {1.2096, 1.2302, 1.2089, 1.2150, 1.2072};
// for pp taken from Fall17_V3b_MC_SF_AK4PF.txt
const Double_t SF_Nom_pp[nBins_Eta] = {1.1432, 1.1815, 1.0989, 1.1137, 1.1307};
const Double_t SF_Down_pp[nBins_Eta] = {1.1210, 1.1332, 1.0533, 0.9740, 0.9837};
const Double_t SF_Up_pp[nBins_Eta] = {1.1654, 1.2299, 1.1444, 1.2533, 1.2778};

const Double_t Eta_Value[nBins_Eta+1] = {0., 0.522, 0.783, 1.131, 1.305, 1.740};

Double_t getJERSFFromEta(const Double_t& jetEta, const Int_t& SF_cond, const TString& colliding_system)
{
  Int_t Etabin = -1;
  for(int i = 0; i < nBins_Eta; ++i)
    {
      if(fabs(jetEta) >= Eta_Value[i] && fabs(jetEta) < Eta_Value[i+1])
	{
	  Etabin = i;
	  break;
	}
    }
  
  double SF = 1.;
  
  if(colliding_system == "PbPb")
    {
      if(SF_cond == 0) SF = SF_Nom_PbPb[Etabin]; // for nominal
      else if(SF_cond == 1) SF = SF_Up_PbPb[Etabin]; // for up
      else if(SF_cond == 2) SF = SF_Down_PbPb[Etabin]; // for down
      else {std::cout<<"Jet eta out of range, please check. I am returning SF = 1"<< std::endl; SF = 1.;}
    }
  else if(colliding_system == "pp")
    {
      if(SF_cond == 0) SF = SF_Nom_pp[Etabin]; // for nominal
      else if(SF_cond == 1) SF = SF_Up_pp[Etabin]; // for up
      else if(SF_cond == 2) SF = SF_Down_pp[Etabin]; // for down
      else {std::cout<<"Jet eta out of range, please check. I am returning SF = 1"<< std::endl; SF = 1.;}
    }
  return SF;
}

void print_start()
{
  cout << endl;
  time_t init = time(0);
  char* init_time = ctime(&init); // convert now to string form                                        
  cout << "Starting at : " << init_time << endl;
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << endl;
  cout << endl;
}

void print_stop()
{
  time_t end = time(0);
  char* end_time = ctime(&end); // convert now to string form                                                   
  cout << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << endl;
  cout << "Stopping at : " << end_time << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << endl;
}
