#include "WWWTree.h"
WWWTree wwwbaby;

void WWWTree::Init(TTree *tree) {
  lep_p4_branch = tree->GetBranch("lep_p4");
  if (lep_p4_branch) lep_p4_branch->SetAddress(&lep_p4_);
  lep_closest_jet_p4_branch = tree->GetBranch("lep_closest_jet_p4");
  if (lep_closest_jet_p4_branch) lep_closest_jet_p4_branch->SetAddress(&lep_closest_jet_p4_);
  gamma_p4_branch = tree->GetBranch("gamma_p4");
  if (gamma_p4_branch) gamma_p4_branch->SetAddress(&gamma_p4_);
  genPart_p4_branch = tree->GetBranch("genPart_p4");
  if (genPart_p4_branch) genPart_p4_branch->SetAddress(&genPart_p4_);
  jets_p4_branch = tree->GetBranch("jets_p4");
  if (jets_p4_branch) jets_p4_branch->SetAddress(&jets_p4_);
  removed_jets_p4_branch = tree->GetBranch("removed_jets_p4");
  if (removed_jets_p4_branch) removed_jets_p4_branch->SetAddress(&removed_jets_p4_);
  jets_medb_p4_branch = tree->GetBranch("jets_medb_p4");
  if (jets_medb_p4_branch) jets_medb_p4_branch->SetAddress(&jets_medb_p4_);
  jets_up_p4_branch = tree->GetBranch("jets_up_p4");
  if (jets_up_p4_branch) jets_up_p4_branch->SetAddress(&jets_up_p4_);
  jets_medb_up_p4_branch = tree->GetBranch("jets_medb_up_p4");
  if (jets_medb_up_p4_branch) jets_medb_up_p4_branch->SetAddress(&jets_medb_up_p4_);
  jets_dn_p4_branch = tree->GetBranch("jets_dn_p4");
  if (jets_dn_p4_branch) jets_dn_p4_branch->SetAddress(&jets_dn_p4_);
  jets_medb_dn_p4_branch = tree->GetBranch("jets_medb_dn_p4");
  if (jets_medb_dn_p4_branch) jets_medb_dn_p4_branch->SetAddress(&jets_medb_dn_p4_);
  highPtPFcands_p4_branch = tree->GetBranch("highPtPFcands_p4");
  if (highPtPFcands_p4_branch) highPtPFcands_p4_branch->SetAddress(&highPtPFcands_p4_);
  decayedphoton_lep1_p4_branch = tree->GetBranch("decayedphoton_lep1_p4");
  if (decayedphoton_lep1_p4_branch) decayedphoton_lep1_p4_branch->SetAddress(&decayedphoton_lep1_p4_);
  decayedphoton_lep2_p4_branch = tree->GetBranch("decayedphoton_lep2_p4");
  if (decayedphoton_lep2_p4_branch) decayedphoton_lep2_p4_branch->SetAddress(&decayedphoton_lep2_p4_);
  decayedphoton_bosn_p4_branch = tree->GetBranch("decayedphoton_bosn_p4");
  if (decayedphoton_bosn_p4_branch) decayedphoton_bosn_p4_branch->SetAddress(&decayedphoton_bosn_p4_);

  tree->SetMakeClass(1);

  run_branch = tree->GetBranch("run");
  if (run_branch) run_branch->SetAddress(&run_);
  lumi_branch = tree->GetBranch("lumi");
  if (lumi_branch) lumi_branch->SetAddress(&lumi_);
  evt_branch = tree->GetBranch("evt");
  if (evt_branch) evt_branch->SetAddress(&evt_);
  isData_branch = tree->GetBranch("isData");
  if (isData_branch) isData_branch->SetAddress(&isData_);
  evt_passgoodrunlist_branch = tree->GetBranch("evt_passgoodrunlist");
  if (evt_passgoodrunlist_branch) evt_passgoodrunlist_branch->SetAddress(&evt_passgoodrunlist_);
  evt_scale1fb_branch = tree->GetBranch("evt_scale1fb");
  if (evt_scale1fb_branch) evt_scale1fb_branch->SetAddress(&evt_scale1fb_);
  evt_xsec_branch = tree->GetBranch("evt_xsec");
  if (evt_xsec_branch) evt_xsec_branch->SetAddress(&evt_xsec_);
  evt_kfactor_branch = tree->GetBranch("evt_kfactor");
  if (evt_kfactor_branch) evt_kfactor_branch->SetAddress(&evt_kfactor_);
  evt_filter_branch = tree->GetBranch("evt_filter");
  if (evt_filter_branch) evt_filter_branch->SetAddress(&evt_filter_);
  evt_nEvts_branch = tree->GetBranch("evt_nEvts");
  if (evt_nEvts_branch) evt_nEvts_branch->SetAddress(&evt_nEvts_);
  evt_dataset_branch = tree->GetBranch("evt_dataset");
  if (evt_dataset_branch) evt_dataset_branch->SetAddress(&evt_dataset_);
  puWeight_branch = tree->GetBranch("puWeight");
  if (puWeight_branch) puWeight_branch->SetAddress(&puWeight_);
  nVert_branch = tree->GetBranch("nVert");
  if (nVert_branch) nVert_branch->SetAddress(&nVert_);
  firstgoodvertex_branch = tree->GetBranch("firstgoodvertex");
  if (firstgoodvertex_branch) firstgoodvertex_branch->SetAddress(&firstgoodvertex_);
  nTrueInt_branch = tree->GetBranch("nTrueInt");
  if (nTrueInt_branch) nTrueInt_branch->SetAddress(&nTrueInt_);
  rho_branch = tree->GetBranch("rho");
  if (rho_branch) rho_branch->SetAddress(&rho_);
  rho25_branch = tree->GetBranch("rho25");
  if (rho25_branch) rho25_branch->SetAddress(&rho25_);
  gen_ht_branch = tree->GetBranch("gen_ht");
  if (gen_ht_branch) gen_ht_branch->SetAddress(&gen_ht_);
  nBJetTight_branch = tree->GetBranch("nBJetTight");
  if (nBJetTight_branch) nBJetTight_branch->SetAddress(&nBJetTight_);
  nBJetMedium_branch = tree->GetBranch("nBJetMedium");
  if (nBJetMedium_branch) nBJetMedium_branch->SetAddress(&nBJetMedium_);
  nBJetLoose_branch = tree->GetBranch("nBJetLoose");
  if (nBJetLoose_branch) nBJetLoose_branch->SetAddress(&nBJetLoose_);
  nBJetTight_up_branch = tree->GetBranch("nBJetTight_up");
  if (nBJetTight_up_branch) nBJetTight_up_branch->SetAddress(&nBJetTight_up_);
  nBJetMedium_up_branch = tree->GetBranch("nBJetMedium_up");
  if (nBJetMedium_up_branch) nBJetMedium_up_branch->SetAddress(&nBJetMedium_up_);
  nBJetLoose_up_branch = tree->GetBranch("nBJetLoose_up");
  if (nBJetLoose_up_branch) nBJetLoose_up_branch->SetAddress(&nBJetLoose_up_);
  nBJetTight_dn_branch = tree->GetBranch("nBJetTight_dn");
  if (nBJetTight_dn_branch) nBJetTight_dn_branch->SetAddress(&nBJetTight_dn_);
  nBJetMedium_dn_branch = tree->GetBranch("nBJetMedium_dn");
  if (nBJetMedium_dn_branch) nBJetMedium_dn_branch->SetAddress(&nBJetMedium_dn_);
  nBJetLoose_dn_branch = tree->GetBranch("nBJetLoose_dn");
  if (nBJetLoose_dn_branch) nBJetLoose_dn_branch->SetAddress(&nBJetLoose_dn_);
  nJet200MuFrac50DphiMet_branch = tree->GetBranch("nJet200MuFrac50DphiMet");
  if (nJet200MuFrac50DphiMet_branch) nJet200MuFrac50DphiMet_branch->SetAddress(&nJet200MuFrac50DphiMet_);
  nMuons10_branch = tree->GetBranch("nMuons10");
  if (nMuons10_branch) nMuons10_branch->SetAddress(&nMuons10_);
  nBadMuons20_branch = tree->GetBranch("nBadMuons20");
  if (nBadMuons20_branch) nBadMuons20_branch->SetAddress(&nBadMuons20_);
  nElectrons10_branch = tree->GetBranch("nElectrons10");
  if (nElectrons10_branch) nElectrons10_branch->SetAddress(&nElectrons10_);
  nGammas20_branch = tree->GetBranch("nGammas20");
  if (nGammas20_branch) nGammas20_branch->SetAddress(&nGammas20_);
  nTaus20_branch = tree->GetBranch("nTaus20");
  if (nTaus20_branch) nTaus20_branch->SetAddress(&nTaus20_);
  met_pt_branch = tree->GetBranch("met_pt");
  if (met_pt_branch) met_pt_branch->SetAddress(&met_pt_);
  met_phi_branch = tree->GetBranch("met_phi");
  if (met_phi_branch) met_phi_branch->SetAddress(&met_phi_);
  met_calo_pt_branch = tree->GetBranch("met_calo_pt");
  if (met_calo_pt_branch) met_calo_pt_branch->SetAddress(&met_calo_pt_);
  met_calo_phi_branch = tree->GetBranch("met_calo_phi");
  if (met_calo_phi_branch) met_calo_phi_branch->SetAddress(&met_calo_phi_);
  met_miniaod_pt_branch = tree->GetBranch("met_miniaod_pt");
  if (met_miniaod_pt_branch) met_miniaod_pt_branch->SetAddress(&met_miniaod_pt_);
  met_miniaod_phi_branch = tree->GetBranch("met_miniaod_phi");
  if (met_miniaod_phi_branch) met_miniaod_phi_branch->SetAddress(&met_miniaod_phi_);
  met_muegclean_pt_branch = tree->GetBranch("met_muegclean_pt");
  if (met_muegclean_pt_branch) met_muegclean_pt_branch->SetAddress(&met_muegclean_pt_);
  met_muegclean_phi_branch = tree->GetBranch("met_muegclean_phi");
  if (met_muegclean_phi_branch) met_muegclean_phi_branch->SetAddress(&met_muegclean_phi_);
  met_rawPt_branch = tree->GetBranch("met_rawPt");
  if (met_rawPt_branch) met_rawPt_branch->SetAddress(&met_rawPt_);
  met_rawPhi_branch = tree->GetBranch("met_rawPhi");
  if (met_rawPhi_branch) met_rawPhi_branch->SetAddress(&met_rawPhi_);
  met_genPt_branch = tree->GetBranch("met_genPt");
  if (met_genPt_branch) met_genPt_branch->SetAddress(&met_genPt_);
  met_genPhi_branch = tree->GetBranch("met_genPhi");
  if (met_genPhi_branch) met_genPhi_branch->SetAddress(&met_genPhi_);
  sumet_raw_branch = tree->GetBranch("sumet_raw");
  if (sumet_raw_branch) sumet_raw_branch->SetAddress(&sumet_raw_);
  Flag_ecalLaserCorrFilter_branch = tree->GetBranch("Flag_ecalLaserCorrFilter");
  if (Flag_ecalLaserCorrFilter_branch) Flag_ecalLaserCorrFilter_branch->SetAddress(&Flag_ecalLaserCorrFilter_);
  Flag_hcalLaserEventFilter_branch = tree->GetBranch("Flag_hcalLaserEventFilter");
  if (Flag_hcalLaserEventFilter_branch) Flag_hcalLaserEventFilter_branch->SetAddress(&Flag_hcalLaserEventFilter_);
  Flag_trackingFailureFilter_branch = tree->GetBranch("Flag_trackingFailureFilter");
  if (Flag_trackingFailureFilter_branch) Flag_trackingFailureFilter_branch->SetAddress(&Flag_trackingFailureFilter_);
  Flag_CSCTightHaloFilter_branch = tree->GetBranch("Flag_CSCTightHaloFilter");
  if (Flag_CSCTightHaloFilter_branch) Flag_CSCTightHaloFilter_branch->SetAddress(&Flag_CSCTightHaloFilter_);
  Flag_HBHENoiseFilter_branch = tree->GetBranch("Flag_HBHENoiseFilter");
  if (Flag_HBHENoiseFilter_branch) Flag_HBHENoiseFilter_branch->SetAddress(&Flag_HBHENoiseFilter_);
  Flag_HBHEIsoNoiseFilter_branch = tree->GetBranch("Flag_HBHEIsoNoiseFilter");
  if (Flag_HBHEIsoNoiseFilter_branch) Flag_HBHEIsoNoiseFilter_branch->SetAddress(&Flag_HBHEIsoNoiseFilter_);
  Flag_CSCTightHalo2015Filter_branch = tree->GetBranch("Flag_CSCTightHalo2015Filter");
  if (Flag_CSCTightHalo2015Filter_branch) Flag_CSCTightHalo2015Filter_branch->SetAddress(&Flag_CSCTightHalo2015Filter_);
  Flag_EcalDeadCellTriggerPrimitiveFilter_branch = tree->GetBranch("Flag_EcalDeadCellTriggerPrimitiveFilter");
  if (Flag_EcalDeadCellTriggerPrimitiveFilter_branch) Flag_EcalDeadCellTriggerPrimitiveFilter_branch->SetAddress(&Flag_EcalDeadCellTriggerPrimitiveFilter_);
  Flag_goodVertices_branch = tree->GetBranch("Flag_goodVertices");
  if (Flag_goodVertices_branch) Flag_goodVertices_branch->SetAddress(&Flag_goodVertices_);
  Flag_eeBadScFilter_branch = tree->GetBranch("Flag_eeBadScFilter");
  if (Flag_eeBadScFilter_branch) Flag_eeBadScFilter_branch->SetAddress(&Flag_eeBadScFilter_);
  Flag_globalTightHalo2016_branch = tree->GetBranch("Flag_globalTightHalo2016");
  if (Flag_globalTightHalo2016_branch) Flag_globalTightHalo2016_branch->SetAddress(&Flag_globalTightHalo2016_);
  Flag_badMuonFilter_branch = tree->GetBranch("Flag_badMuonFilter");
  if (Flag_badMuonFilter_branch) Flag_badMuonFilter_branch->SetAddress(&Flag_badMuonFilter_);
  Flag_badChargedCandidateFilter_branch = tree->GetBranch("Flag_badChargedCandidateFilter");
  if (Flag_badChargedCandidateFilter_branch) Flag_badChargedCandidateFilter_branch->SetAddress(&Flag_badChargedCandidateFilter_);
  Flag_badMuonFilterv2_branch = tree->GetBranch("Flag_badMuonFilterv2");
  if (Flag_badMuonFilterv2_branch) Flag_badMuonFilterv2_branch->SetAddress(&Flag_badMuonFilterv2_);
  Flag_badChargedCandidateFilterv2_branch = tree->GetBranch("Flag_badChargedCandidateFilterv2");
  if (Flag_badChargedCandidateFilterv2_branch) Flag_badChargedCandidateFilterv2_branch->SetAddress(&Flag_badChargedCandidateFilterv2_);
  Flag_badMuons_branch = tree->GetBranch("Flag_badMuons");
  if (Flag_badMuons_branch) Flag_badMuons_branch->SetAddress(&Flag_badMuons_);
  Flag_duplicateMuons_branch = tree->GetBranch("Flag_duplicateMuons");
  if (Flag_duplicateMuons_branch) Flag_duplicateMuons_branch->SetAddress(&Flag_duplicateMuons_);
  Flag_noBadMuons_branch = tree->GetBranch("Flag_noBadMuons");
  if (Flag_noBadMuons_branch) Flag_noBadMuons_branch->SetAddress(&Flag_noBadMuons_);
  HLT_singleEl_branch = tree->GetBranch("HLT_singleEl");
  if (HLT_singleEl_branch) HLT_singleEl_branch->SetAddress(&HLT_singleEl_);
  HLT_singleMu_branch = tree->GetBranch("HLT_singleMu");
  if (HLT_singleMu_branch) HLT_singleMu_branch->SetAddress(&HLT_singleMu_);
  HLT_singleMu_noiso_branch = tree->GetBranch("HLT_singleMu_noiso");
  if (HLT_singleMu_noiso_branch) HLT_singleMu_noiso_branch->SetAddress(&HLT_singleMu_noiso_);
  HLT_DoubleEl_noiso_branch = tree->GetBranch("HLT_DoubleEl_noiso");
  if (HLT_DoubleEl_noiso_branch) HLT_DoubleEl_noiso_branch->SetAddress(&HLT_DoubleEl_noiso_);
  HLT_DoubleEl_branch = tree->GetBranch("HLT_DoubleEl");
  if (HLT_DoubleEl_branch) HLT_DoubleEl_branch->SetAddress(&HLT_DoubleEl_);
  HLT_DoubleEl_DZ_branch = tree->GetBranch("HLT_DoubleEl_DZ");
  if (HLT_DoubleEl_DZ_branch) HLT_DoubleEl_DZ_branch->SetAddress(&HLT_DoubleEl_DZ_);
  HLT_DoubleEl_DZ_2_branch = tree->GetBranch("HLT_DoubleEl_DZ_2");
  if (HLT_DoubleEl_DZ_2_branch) HLT_DoubleEl_DZ_2_branch->SetAddress(&HLT_DoubleEl_DZ_2_);
  HLT_MuEG_branch = tree->GetBranch("HLT_MuEG");
  if (HLT_MuEG_branch) HLT_MuEG_branch->SetAddress(&HLT_MuEG_);
  HLT_MuEG_2_branch = tree->GetBranch("HLT_MuEG_2");
  if (HLT_MuEG_2_branch) HLT_MuEG_2_branch->SetAddress(&HLT_MuEG_2_);
  HLT_MuEG_noiso_branch = tree->GetBranch("HLT_MuEG_noiso");
  if (HLT_MuEG_noiso_branch) HLT_MuEG_noiso_branch->SetAddress(&HLT_MuEG_noiso_);
  HLT_MuEG_noiso_2_branch = tree->GetBranch("HLT_MuEG_noiso_2");
  if (HLT_MuEG_noiso_2_branch) HLT_MuEG_noiso_2_branch->SetAddress(&HLT_MuEG_noiso_2_);
  HLT_Mu8_EG17_branch = tree->GetBranch("HLT_Mu8_EG17");
  if (HLT_Mu8_EG17_branch) HLT_Mu8_EG17_branch->SetAddress(&HLT_Mu8_EG17_);
  HLT_Mu8_EG23_branch = tree->GetBranch("HLT_Mu8_EG23");
  if (HLT_Mu8_EG23_branch) HLT_Mu8_EG23_branch->SetAddress(&HLT_Mu8_EG23_);
  HLT_Mu8_EG23_DZ_branch = tree->GetBranch("HLT_Mu8_EG23_DZ");
  if (HLT_Mu8_EG23_DZ_branch) HLT_Mu8_EG23_DZ_branch->SetAddress(&HLT_Mu8_EG23_DZ_);
  HLT_Mu12_EG23_DZ_branch = tree->GetBranch("HLT_Mu12_EG23_DZ");
  if (HLT_Mu12_EG23_DZ_branch) HLT_Mu12_EG23_DZ_branch->SetAddress(&HLT_Mu12_EG23_DZ_);
  HLT_Mu17_EG12_branch = tree->GetBranch("HLT_Mu17_EG12");
  if (HLT_Mu17_EG12_branch) HLT_Mu17_EG12_branch->SetAddress(&HLT_Mu17_EG12_);
  HLT_Mu23_EG8_branch = tree->GetBranch("HLT_Mu23_EG8");
  if (HLT_Mu23_EG8_branch) HLT_Mu23_EG8_branch->SetAddress(&HLT_Mu23_EG8_);
  HLT_Mu23_EG8_DZ_branch = tree->GetBranch("HLT_Mu23_EG8_DZ");
  if (HLT_Mu23_EG8_DZ_branch) HLT_Mu23_EG8_DZ_branch->SetAddress(&HLT_Mu23_EG8_DZ_);
  HLT_Mu23_EG12_branch = tree->GetBranch("HLT_Mu23_EG12");
  if (HLT_Mu23_EG12_branch) HLT_Mu23_EG12_branch->SetAddress(&HLT_Mu23_EG12_);
  HLT_Mu23_EG12_DZ_branch = tree->GetBranch("HLT_Mu23_EG12_DZ");
  if (HLT_Mu23_EG12_DZ_branch) HLT_Mu23_EG12_DZ_branch->SetAddress(&HLT_Mu23_EG12_DZ_);
  HLT_DoubleMu_noiso_branch = tree->GetBranch("HLT_DoubleMu_noiso");
  if (HLT_DoubleMu_noiso_branch) HLT_DoubleMu_noiso_branch->SetAddress(&HLT_DoubleMu_noiso_);
  HLT_DoubleMu_noiso_27_8_branch = tree->GetBranch("HLT_DoubleMu_noiso_27_8");
  if (HLT_DoubleMu_noiso_27_8_branch) HLT_DoubleMu_noiso_27_8_branch->SetAddress(&HLT_DoubleMu_noiso_27_8_);
  HLT_DoubleMu_noiso_30_11_branch = tree->GetBranch("HLT_DoubleMu_noiso_30_11");
  if (HLT_DoubleMu_noiso_30_11_branch) HLT_DoubleMu_noiso_30_11_branch->SetAddress(&HLT_DoubleMu_noiso_30_11_);
  HLT_DoubleMu_noiso_40_11_branch = tree->GetBranch("HLT_DoubleMu_noiso_40_11");
  if (HLT_DoubleMu_noiso_40_11_branch) HLT_DoubleMu_noiso_40_11_branch->SetAddress(&HLT_DoubleMu_noiso_40_11_);
  HLT_DoubleMu_branch = tree->GetBranch("HLT_DoubleMu");
  if (HLT_DoubleMu_branch) HLT_DoubleMu_branch->SetAddress(&HLT_DoubleMu_);
  HLT_DoubleMu_tk_branch = tree->GetBranch("HLT_DoubleMu_tk");
  if (HLT_DoubleMu_tk_branch) HLT_DoubleMu_tk_branch->SetAddress(&HLT_DoubleMu_tk_);
  HLT_DoubleMu_dbltk_branch = tree->GetBranch("HLT_DoubleMu_dbltk");
  if (HLT_DoubleMu_dbltk_branch) HLT_DoubleMu_dbltk_branch->SetAddress(&HLT_DoubleMu_dbltk_);
  HLT_DoubleMu_nonDZ_branch = tree->GetBranch("HLT_DoubleMu_nonDZ");
  if (HLT_DoubleMu_nonDZ_branch) HLT_DoubleMu_nonDZ_branch->SetAddress(&HLT_DoubleMu_nonDZ_);
  HLT_DoubleMu_tk_nonDZ_branch = tree->GetBranch("HLT_DoubleMu_tk_nonDZ");
  if (HLT_DoubleMu_tk_nonDZ_branch) HLT_DoubleMu_tk_nonDZ_branch->SetAddress(&HLT_DoubleMu_tk_nonDZ_);
  HLT_Photon22_R9Id90_HE10_IsoM_branch = tree->GetBranch("HLT_Photon22_R9Id90_HE10_IsoM");
  if (HLT_Photon22_R9Id90_HE10_IsoM_branch) HLT_Photon22_R9Id90_HE10_IsoM_branch->SetAddress(&HLT_Photon22_R9Id90_HE10_IsoM_);
  HLT_Photon30_R9Id90_HE10_IsoM_branch = tree->GetBranch("HLT_Photon30_R9Id90_HE10_IsoM");
  if (HLT_Photon30_R9Id90_HE10_IsoM_branch) HLT_Photon30_R9Id90_HE10_IsoM_branch->SetAddress(&HLT_Photon30_R9Id90_HE10_IsoM_);
  HLT_Photon36_R9Id90_HE10_IsoM_branch = tree->GetBranch("HLT_Photon36_R9Id90_HE10_IsoM");
  if (HLT_Photon36_R9Id90_HE10_IsoM_branch) HLT_Photon36_R9Id90_HE10_IsoM_branch->SetAddress(&HLT_Photon36_R9Id90_HE10_IsoM_);
  HLT_Photon50_R9Id90_HE10_IsoM_branch = tree->GetBranch("HLT_Photon50_R9Id90_HE10_IsoM");
  if (HLT_Photon50_R9Id90_HE10_IsoM_branch) HLT_Photon50_R9Id90_HE10_IsoM_branch->SetAddress(&HLT_Photon50_R9Id90_HE10_IsoM_);
  HLT_Photon75_R9Id90_HE10_IsoM_branch = tree->GetBranch("HLT_Photon75_R9Id90_HE10_IsoM");
  if (HLT_Photon75_R9Id90_HE10_IsoM_branch) HLT_Photon75_R9Id90_HE10_IsoM_branch->SetAddress(&HLT_Photon75_R9Id90_HE10_IsoM_);
  HLT_Photon90_R9Id90_HE10_IsoM_branch = tree->GetBranch("HLT_Photon90_R9Id90_HE10_IsoM");
  if (HLT_Photon90_R9Id90_HE10_IsoM_branch) HLT_Photon90_R9Id90_HE10_IsoM_branch->SetAddress(&HLT_Photon90_R9Id90_HE10_IsoM_);
  HLT_Photon120_R9Id90_HE10_IsoM_branch = tree->GetBranch("HLT_Photon120_R9Id90_HE10_IsoM");
  if (HLT_Photon120_R9Id90_HE10_IsoM_branch) HLT_Photon120_R9Id90_HE10_IsoM_branch->SetAddress(&HLT_Photon120_R9Id90_HE10_IsoM_);
  HLT_Photon165_R9Id90_HE10_IsoM_branch = tree->GetBranch("HLT_Photon165_R9Id90_HE10_IsoM");
  if (HLT_Photon165_R9Id90_HE10_IsoM_branch) HLT_Photon165_R9Id90_HE10_IsoM_branch->SetAddress(&HLT_Photon165_R9Id90_HE10_IsoM_);
  HLT_Photon165_HE10_branch = tree->GetBranch("HLT_Photon165_HE10");
  if (HLT_Photon165_HE10_branch) HLT_Photon165_HE10_branch->SetAddress(&HLT_Photon165_HE10_);
  HLT_CaloJet500_NoJetID_branch = tree->GetBranch("HLT_CaloJet500_NoJetID");
  if (HLT_CaloJet500_NoJetID_branch) HLT_CaloJet500_NoJetID_branch->SetAddress(&HLT_CaloJet500_NoJetID_);
  HLT_ECALHT800_NoJetID_branch = tree->GetBranch("HLT_ECALHT800_NoJetID");
  if (HLT_ECALHT800_NoJetID_branch) HLT_ECALHT800_NoJetID_branch->SetAddress(&HLT_ECALHT800_NoJetID_);
  HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton_branch = tree->GetBranch("HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton");
  if (HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton_branch) HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton_branch->SetAddress(&HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton_);
  HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton_branch = tree->GetBranch("HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton");
  if (HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton_branch) HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton_branch->SetAddress(&HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton_);
  HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton_branch = tree->GetBranch("HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton");
  if (HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton_branch) HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton_branch->SetAddress(&HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton_);
  HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton_branch = tree->GetBranch("HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton");
  if (HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton_branch) HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton_branch->SetAddress(&HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton_);
  HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton_branch = tree->GetBranch("HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton");
  if (HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton_branch) HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton_branch->SetAddress(&HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton_);
  HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton_branch = tree->GetBranch("HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton");
  if (HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton_branch) HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton_branch->SetAddress(&HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton_);
  HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton_branch = tree->GetBranch("HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton");
  if (HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton_branch) HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton_branch->SetAddress(&HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton_);
  HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton_branch = tree->GetBranch("HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton");
  if (HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton_branch) HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton_branch->SetAddress(&HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton_);
  HLT_Photon165_HE10_matchedtophoton_branch = tree->GetBranch("HLT_Photon165_HE10_matchedtophoton");
  if (HLT_Photon165_HE10_matchedtophoton_branch) HLT_Photon165_HE10_matchedtophoton_branch->SetAddress(&HLT_Photon165_HE10_matchedtophoton_);
  dilmass_branch = tree->GetBranch("dilmass");
  if (dilmass_branch) dilmass_branch->SetAddress(&dilmass_);
  dilpt_branch = tree->GetBranch("dilpt");
  if (dilpt_branch) dilpt_branch->SetAddress(&dilpt_);
  dRll_branch = tree->GetBranch("dRll");
  if (dRll_branch) dRll_branch->SetAddress(&dRll_);
  matched_neutralemf_branch = tree->GetBranch("matched_neutralemf");
  if (matched_neutralemf_branch) matched_neutralemf_branch->SetAddress(&matched_neutralemf_);
  matched_emf_branch = tree->GetBranch("matched_emf");
  if (matched_emf_branch) matched_emf_branch->SetAddress(&matched_emf_);
  elveto_branch = tree->GetBranch("elveto");
  if (elveto_branch) elveto_branch->SetAddress(&elveto_);
  nlep_branch = tree->GetBranch("nlep");
  if (nlep_branch) nlep_branch->SetAddress(&nlep_);
  nveto_leptons_branch = tree->GetBranch("nveto_leptons");
  if (nveto_leptons_branch) nveto_leptons_branch->SetAddress(&nveto_leptons_);
  nVetoEl_relIso03EAless01_branch = tree->GetBranch("nVetoEl_relIso03EAless01");
  if (nVetoEl_relIso03EAless01_branch) nVetoEl_relIso03EAless01_branch->SetAddress(&nVetoEl_relIso03EAless01_);
  nVetoEl_relIso03EAless02_branch = tree->GetBranch("nVetoEl_relIso03EAless02");
  if (nVetoEl_relIso03EAless02_branch) nVetoEl_relIso03EAless02_branch->SetAddress(&nVetoEl_relIso03EAless02_);
  nVetoEl_relIso03EAless03_branch = tree->GetBranch("nVetoEl_relIso03EAless03");
  if (nVetoEl_relIso03EAless03_branch) nVetoEl_relIso03EAless03_branch->SetAddress(&nVetoEl_relIso03EAless03_);
  nVetoEl_relIso03EAless04_branch = tree->GetBranch("nVetoEl_relIso03EAless04");
  if (nVetoEl_relIso03EAless04_branch) nVetoEl_relIso03EAless04_branch->SetAddress(&nVetoEl_relIso03EAless04_);
  nVetoMu_relIso03EAless01_branch = tree->GetBranch("nVetoMu_relIso03EAless01");
  if (nVetoMu_relIso03EAless01_branch) nVetoMu_relIso03EAless01_branch->SetAddress(&nVetoMu_relIso03EAless01_);
  nVetoMu_relIso03EAless02_branch = tree->GetBranch("nVetoMu_relIso03EAless02");
  if (nVetoMu_relIso03EAless02_branch) nVetoMu_relIso03EAless02_branch->SetAddress(&nVetoMu_relIso03EAless02_);
  nVetoMu_relIso03EAless03_branch = tree->GetBranch("nVetoMu_relIso03EAless03");
  if (nVetoMu_relIso03EAless03_branch) nVetoMu_relIso03EAless03_branch->SetAddress(&nVetoMu_relIso03EAless03_);
  nVetoMu_relIso03EAless04_branch = tree->GetBranch("nVetoMu_relIso03EAless04");
  if (nVetoMu_relIso03EAless04_branch) nVetoMu_relIso03EAless04_branch->SetAddress(&nVetoMu_relIso03EAless04_);
  lep_pt_branch = tree->GetBranch("lep_pt");
  if (lep_pt_branch) lep_pt_branch->SetAddress(&lep_pt_);
  lep_eta_branch = tree->GetBranch("lep_eta");
  if (lep_eta_branch) lep_eta_branch->SetAddress(&lep_eta_);
  lep_phi_branch = tree->GetBranch("lep_phi");
  if (lep_phi_branch) lep_phi_branch->SetAddress(&lep_phi_);
  lep_mass_branch = tree->GetBranch("lep_mass");
  if (lep_mass_branch) lep_mass_branch->SetAddress(&lep_mass_);
  lep_charge_branch = tree->GetBranch("lep_charge");
  if (lep_charge_branch) lep_charge_branch->SetAddress(&lep_charge_);
  lep_3ch_agree_branch = tree->GetBranch("lep_3ch_agree");
  if (lep_3ch_agree_branch) lep_3ch_agree_branch->SetAddress(&lep_3ch_agree_);
  lep_isFromW_branch = tree->GetBranch("lep_isFromW");
  if (lep_isFromW_branch) lep_isFromW_branch->SetAddress(&lep_isFromW_);
  lep_isFromZ_branch = tree->GetBranch("lep_isFromZ");
  if (lep_isFromZ_branch) lep_isFromZ_branch->SetAddress(&lep_isFromZ_);
  lep_isFromB_branch = tree->GetBranch("lep_isFromB");
  if (lep_isFromB_branch) lep_isFromB_branch->SetAddress(&lep_isFromB_);
  lep_isFromC_branch = tree->GetBranch("lep_isFromC");
  if (lep_isFromC_branch) lep_isFromC_branch->SetAddress(&lep_isFromC_);
  lep_isFromL_branch = tree->GetBranch("lep_isFromL");
  if (lep_isFromL_branch) lep_isFromL_branch->SetAddress(&lep_isFromL_);
  lep_isFromLF_branch = tree->GetBranch("lep_isFromLF");
  if (lep_isFromLF_branch) lep_isFromLF_branch->SetAddress(&lep_isFromLF_);
  lep_ptRatio_branch = tree->GetBranch("lep_ptRatio");
  if (lep_ptRatio_branch) lep_ptRatio_branch->SetAddress(&lep_ptRatio_);
  lep_motherIdSS_branch = tree->GetBranch("lep_motherIdSS");
  if (lep_motherIdSS_branch) lep_motherIdSS_branch->SetAddress(&lep_motherIdSS_);
  lep_genPart_index_branch = tree->GetBranch("lep_genPart_index");
  if (lep_genPart_index_branch) lep_genPart_index_branch->SetAddress(&lep_genPart_index_);
  lep_coneCorrPt_branch = tree->GetBranch("lep_coneCorrPt");
  if (lep_coneCorrPt_branch) lep_coneCorrPt_branch->SetAddress(&lep_coneCorrPt_);
  lep_ptRel_branch = tree->GetBranch("lep_ptRel");
  if (lep_ptRel_branch) lep_ptRel_branch->SetAddress(&lep_ptRel_);
  lep_relIso03_branch = tree->GetBranch("lep_relIso03");
  if (lep_relIso03_branch) lep_relIso03_branch->SetAddress(&lep_relIso03_);
  lep_relIso03DB_branch = tree->GetBranch("lep_relIso03DB");
  if (lep_relIso03DB_branch) lep_relIso03DB_branch->SetAddress(&lep_relIso03DB_);
  lep_relIso03EA_branch = tree->GetBranch("lep_relIso03EA");
  if (lep_relIso03EA_branch) lep_relIso03EA_branch->SetAddress(&lep_relIso03EA_);
  lep_relIso03EAv2_branch = tree->GetBranch("lep_relIso03EAv2");
  if (lep_relIso03EAv2_branch) lep_relIso03EAv2_branch->SetAddress(&lep_relIso03EAv2_);
  lep_relIso04DB_branch = tree->GetBranch("lep_relIso04DB");
  if (lep_relIso04DB_branch) lep_relIso04DB_branch->SetAddress(&lep_relIso04DB_);
  lep_relIso04EA_branch = tree->GetBranch("lep_relIso04EA");
  if (lep_relIso04EA_branch) lep_relIso04EA_branch->SetAddress(&lep_relIso04EA_);
  lep_relIso04EAv2_branch = tree->GetBranch("lep_relIso04EAv2");
  if (lep_relIso04EAv2_branch) lep_relIso04EAv2_branch->SetAddress(&lep_relIso04EAv2_);
  lep_miniRelIsoCMS3_EA_branch = tree->GetBranch("lep_miniRelIsoCMS3_EA");
  if (lep_miniRelIsoCMS3_EA_branch) lep_miniRelIsoCMS3_EA_branch->SetAddress(&lep_miniRelIsoCMS3_EA_);
  lep_miniRelIsoCMS3_EAv2_branch = tree->GetBranch("lep_miniRelIsoCMS3_EAv2");
  if (lep_miniRelIsoCMS3_EAv2_branch) lep_miniRelIsoCMS3_EAv2_branch->SetAddress(&lep_miniRelIsoCMS3_EAv2_);
  lep_miniRelIsoCMS3_DB_branch = tree->GetBranch("lep_miniRelIsoCMS3_DB");
  if (lep_miniRelIsoCMS3_DB_branch) lep_miniRelIsoCMS3_DB_branch->SetAddress(&lep_miniRelIsoCMS3_DB_);
  lep_isTriggerSafenoIso_v1_branch = tree->GetBranch("lep_isTriggerSafenoIso_v1");
  if (lep_isTriggerSafenoIso_v1_branch) lep_isTriggerSafenoIso_v1_branch->SetAddress(&lep_isTriggerSafenoIso_v1_);
  lep_isTriggerSafenoIso_v2_branch = tree->GetBranch("lep_isTriggerSafenoIso_v2");
  if (lep_isTriggerSafenoIso_v2_branch) lep_isTriggerSafenoIso_v2_branch->SetAddress(&lep_isTriggerSafenoIso_v2_);
  lep_isTriggerSafe_v1_branch = tree->GetBranch("lep_isTriggerSafe_v1");
  if (lep_isTriggerSafe_v1_branch) lep_isTriggerSafe_v1_branch->SetAddress(&lep_isTriggerSafe_v1_);
  lep_isTriggerSafe_v2_branch = tree->GetBranch("lep_isTriggerSafe_v2");
  if (lep_isTriggerSafe_v2_branch) lep_isTriggerSafe_v2_branch->SetAddress(&lep_isTriggerSafe_v2_);
  lep_pass_VVV_cutbased_veto_branch = tree->GetBranch("lep_pass_VVV_cutbased_veto");
  if (lep_pass_VVV_cutbased_veto_branch) lep_pass_VVV_cutbased_veto_branch->SetAddress(&lep_pass_VVV_cutbased_veto_);
  lep_pass_VVV_cutbased_veto_noiso_branch = tree->GetBranch("lep_pass_VVV_cutbased_veto_noiso");
  if (lep_pass_VVV_cutbased_veto_noiso_branch) lep_pass_VVV_cutbased_veto_noiso_branch->SetAddress(&lep_pass_VVV_cutbased_veto_noiso_);
  lep_pass_VVV_cutbased_veto_noiso_noip_branch = tree->GetBranch("lep_pass_VVV_cutbased_veto_noiso_noip");
  if (lep_pass_VVV_cutbased_veto_noiso_noip_branch) lep_pass_VVV_cutbased_veto_noiso_noip_branch->SetAddress(&lep_pass_VVV_cutbased_veto_noiso_noip_);
  lep_pass_VVV_cutbased_veto_bak1_branch = tree->GetBranch("lep_pass_VVV_cutbased_veto_bak1");
  if (lep_pass_VVV_cutbased_veto_bak1_branch) lep_pass_VVV_cutbased_veto_bak1_branch->SetAddress(&lep_pass_VVV_cutbased_veto_bak1_);
  lep_pass_VVV_cutbased_veto_noiso_bak1_branch = tree->GetBranch("lep_pass_VVV_cutbased_veto_noiso_bak1");
  if (lep_pass_VVV_cutbased_veto_noiso_bak1_branch) lep_pass_VVV_cutbased_veto_noiso_bak1_branch->SetAddress(&lep_pass_VVV_cutbased_veto_noiso_bak1_);
  lep_pass_VVV_cutbased_veto_noiso_noip_bak1_branch = tree->GetBranch("lep_pass_VVV_cutbased_veto_noiso_noip_bak1");
  if (lep_pass_VVV_cutbased_veto_noiso_noip_bak1_branch) lep_pass_VVV_cutbased_veto_noiso_noip_bak1_branch->SetAddress(&lep_pass_VVV_cutbased_veto_noiso_noip_bak1_);
  lep_pass_VVV_cutbased_fo_branch = tree->GetBranch("lep_pass_VVV_cutbased_fo");
  if (lep_pass_VVV_cutbased_fo_branch) lep_pass_VVV_cutbased_fo_branch->SetAddress(&lep_pass_VVV_cutbased_fo_);
  lep_pass_VVV_cutbased_fo_noiso_branch = tree->GetBranch("lep_pass_VVV_cutbased_fo_noiso");
  if (lep_pass_VVV_cutbased_fo_noiso_branch) lep_pass_VVV_cutbased_fo_noiso_branch->SetAddress(&lep_pass_VVV_cutbased_fo_noiso_);
  lep_pass_VVV_cutbased_tight_noiso_branch = tree->GetBranch("lep_pass_VVV_cutbased_tight_noiso");
  if (lep_pass_VVV_cutbased_tight_noiso_branch) lep_pass_VVV_cutbased_tight_noiso_branch->SetAddress(&lep_pass_VVV_cutbased_tight_noiso_);
  lep_pass_VVV_cutbased_tight_branch = tree->GetBranch("lep_pass_VVV_cutbased_tight");
  if (lep_pass_VVV_cutbased_tight_branch) lep_pass_VVV_cutbased_tight_branch->SetAddress(&lep_pass_VVV_cutbased_tight_);
  lep_pass_VVV_MVAbased_tight_noiso_branch = tree->GetBranch("lep_pass_VVV_MVAbased_tight_noiso");
  if (lep_pass_VVV_MVAbased_tight_noiso_branch) lep_pass_VVV_MVAbased_tight_noiso_branch->SetAddress(&lep_pass_VVV_MVAbased_tight_noiso_);
  lep_pass_VVV_MVAbased_tight_branch = tree->GetBranch("lep_pass_VVV_MVAbased_tight");
  if (lep_pass_VVV_MVAbased_tight_branch) lep_pass_VVV_MVAbased_tight_branch->SetAddress(&lep_pass_VVV_MVAbased_tight_);
  lep_pass_VVV_baseline_branch = tree->GetBranch("lep_pass_VVV_baseline");
  if (lep_pass_VVV_baseline_branch) lep_pass_VVV_baseline_branch->SetAddress(&lep_pass_VVV_baseline_);
  nlep_VVV_cutbased_veto_branch = tree->GetBranch("nlep_VVV_cutbased_veto");
  if (nlep_VVV_cutbased_veto_branch) nlep_VVV_cutbased_veto_branch->SetAddress(&nlep_VVV_cutbased_veto_);
  nlep_VVV_cutbased_veto_noiso_branch = tree->GetBranch("nlep_VVV_cutbased_veto_noiso");
  if (nlep_VVV_cutbased_veto_noiso_branch) nlep_VVV_cutbased_veto_noiso_branch->SetAddress(&nlep_VVV_cutbased_veto_noiso_);
  nlep_VVV_cutbased_veto_noiso_noip_branch = tree->GetBranch("nlep_VVV_cutbased_veto_noiso_noip");
  if (nlep_VVV_cutbased_veto_noiso_noip_branch) nlep_VVV_cutbased_veto_noiso_noip_branch->SetAddress(&nlep_VVV_cutbased_veto_noiso_noip_);
  nlep_VVV_cutbased_veto_bak1_branch = tree->GetBranch("nlep_VVV_cutbased_veto_bak1");
  if (nlep_VVV_cutbased_veto_bak1_branch) nlep_VVV_cutbased_veto_bak1_branch->SetAddress(&nlep_VVV_cutbased_veto_bak1_);
  nlep_VVV_cutbased_veto_noiso_bak1_branch = tree->GetBranch("nlep_VVV_cutbased_veto_noiso_bak1");
  if (nlep_VVV_cutbased_veto_noiso_bak1_branch) nlep_VVV_cutbased_veto_noiso_bak1_branch->SetAddress(&nlep_VVV_cutbased_veto_noiso_bak1_);
  nlep_VVV_cutbased_veto_noiso_noip_bak1_branch = tree->GetBranch("nlep_VVV_cutbased_veto_noiso_noip_bak1");
  if (nlep_VVV_cutbased_veto_noiso_noip_bak1_branch) nlep_VVV_cutbased_veto_noiso_noip_bak1_branch->SetAddress(&nlep_VVV_cutbased_veto_noiso_noip_bak1_);
  nlep_VVV_cutbased_fo_branch = tree->GetBranch("nlep_VVV_cutbased_fo");
  if (nlep_VVV_cutbased_fo_branch) nlep_VVV_cutbased_fo_branch->SetAddress(&nlep_VVV_cutbased_fo_);
  nlep_VVV_cutbased_fo_noiso_branch = tree->GetBranch("nlep_VVV_cutbased_fo_noiso");
  if (nlep_VVV_cutbased_fo_noiso_branch) nlep_VVV_cutbased_fo_noiso_branch->SetAddress(&nlep_VVV_cutbased_fo_noiso_);
  nlep_VVV_cutbased_tight_noiso_branch = tree->GetBranch("nlep_VVV_cutbased_tight_noiso");
  if (nlep_VVV_cutbased_tight_noiso_branch) nlep_VVV_cutbased_tight_noiso_branch->SetAddress(&nlep_VVV_cutbased_tight_noiso_);
  nlep_VVV_cutbased_tight_branch = tree->GetBranch("nlep_VVV_cutbased_tight");
  if (nlep_VVV_cutbased_tight_branch) nlep_VVV_cutbased_tight_branch->SetAddress(&nlep_VVV_cutbased_tight_);
  nlep_VVV_MVAbased_tight_noiso_branch = tree->GetBranch("nlep_VVV_MVAbased_tight_noiso");
  if (nlep_VVV_MVAbased_tight_noiso_branch) nlep_VVV_MVAbased_tight_noiso_branch->SetAddress(&nlep_VVV_MVAbased_tight_noiso_);
  nlep_VVV_MVAbased_tight_branch = tree->GetBranch("nlep_VVV_MVAbased_tight");
  if (nlep_VVV_MVAbased_tight_branch) nlep_VVV_MVAbased_tight_branch->SetAddress(&nlep_VVV_MVAbased_tight_);
  nlep_VVV_baseline_branch = tree->GetBranch("nlep_VVV_baseline");
  if (nlep_VVV_baseline_branch) nlep_VVV_baseline_branch->SetAddress(&nlep_VVV_baseline_);
  lep_pdgId_branch = tree->GetBranch("lep_pdgId");
  if (lep_pdgId_branch) lep_pdgId_branch->SetAddress(&lep_pdgId_);
  lep_mc_Id_branch = tree->GetBranch("lep_mc_Id");
  if (lep_mc_Id_branch) lep_mc_Id_branch->SetAddress(&lep_mc_Id_);
  lep_dxy_branch = tree->GetBranch("lep_dxy");
  if (lep_dxy_branch) lep_dxy_branch->SetAddress(&lep_dxy_);
  lep_ip3d_branch = tree->GetBranch("lep_ip3d");
  if (lep_ip3d_branch) lep_ip3d_branch->SetAddress(&lep_ip3d_);
  lep_ip3derr_branch = tree->GetBranch("lep_ip3derr");
  if (lep_ip3derr_branch) lep_ip3derr_branch->SetAddress(&lep_ip3derr_);
  lep_etaSC_branch = tree->GetBranch("lep_etaSC");
  if (lep_etaSC_branch) lep_etaSC_branch->SetAddress(&lep_etaSC_);
  lep_r9_branch = tree->GetBranch("lep_r9");
  if (lep_r9_branch) lep_r9_branch->SetAddress(&lep_r9_);
  lep_dz_branch = tree->GetBranch("lep_dz");
  if (lep_dz_branch) lep_dz_branch->SetAddress(&lep_dz_);
  lep_tightId_branch = tree->GetBranch("lep_tightId");
  if (lep_tightId_branch) lep_tightId_branch->SetAddress(&lep_tightId_);
  lep_mcMatchId_branch = tree->GetBranch("lep_mcMatchId");
  if (lep_mcMatchId_branch) lep_mcMatchId_branch->SetAddress(&lep_mcMatchId_);
  lep_lostHits_branch = tree->GetBranch("lep_lostHits");
  if (lep_lostHits_branch) lep_lostHits_branch->SetAddress(&lep_lostHits_);
  lep_convVeto_branch = tree->GetBranch("lep_convVeto");
  if (lep_convVeto_branch) lep_convVeto_branch->SetAddress(&lep_convVeto_);
  lep_tightCharge_branch = tree->GetBranch("lep_tightCharge");
  if (lep_tightCharge_branch) lep_tightCharge_branch->SetAddress(&lep_tightCharge_);
  lep_MVA_branch = tree->GetBranch("lep_MVA");
  if (lep_MVA_branch) lep_MVA_branch->SetAddress(&lep_MVA_);
  lep_validfraction_branch = tree->GetBranch("lep_validfraction");
  if (lep_validfraction_branch) lep_validfraction_branch->SetAddress(&lep_validfraction_);
  lep_pterr_branch = tree->GetBranch("lep_pterr");
  if (lep_pterr_branch) lep_pterr_branch->SetAddress(&lep_pterr_);
  lep_sta_pterrOpt_branch = tree->GetBranch("lep_sta_pterrOpt");
  if (lep_sta_pterrOpt_branch) lep_sta_pterrOpt_branch->SetAddress(&lep_sta_pterrOpt_);
  lep_glb_pterrOpt_branch = tree->GetBranch("lep_glb_pterrOpt");
  if (lep_glb_pterrOpt_branch) lep_glb_pterrOpt_branch->SetAddress(&lep_glb_pterrOpt_);
  lep_x2ondof_branch = tree->GetBranch("lep_x2ondof");
  if (lep_x2ondof_branch) lep_x2ondof_branch->SetAddress(&lep_x2ondof_);
  lep_sta_x2ondof_branch = tree->GetBranch("lep_sta_x2ondof");
  if (lep_sta_x2ondof_branch) lep_sta_x2ondof_branch->SetAddress(&lep_sta_x2ondof_);
  lep_glb_x2ondof_branch = tree->GetBranch("lep_glb_x2ondof");
  if (lep_glb_x2ondof_branch) lep_glb_x2ondof_branch->SetAddress(&lep_glb_x2ondof_);
  nisoTrack_5gev_branch = tree->GetBranch("nisoTrack_5gev");
  if (nisoTrack_5gev_branch) nisoTrack_5gev_branch->SetAddress(&nisoTrack_5gev_);
  nisoTrack_mt2_branch = tree->GetBranch("nisoTrack_mt2");
  if (nisoTrack_mt2_branch) nisoTrack_mt2_branch->SetAddress(&nisoTrack_mt2_);
  nisoTrack_mt2_cleaned_VVV_cutbased_veto_branch = tree->GetBranch("nisoTrack_mt2_cleaned_VVV_cutbased_veto");
  if (nisoTrack_mt2_cleaned_VVV_cutbased_veto_branch) nisoTrack_mt2_cleaned_VVV_cutbased_veto_branch->SetAddress(&nisoTrack_mt2_cleaned_VVV_cutbased_veto_);
  nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_branch = tree->GetBranch("nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso");
  if (nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_branch) nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_branch->SetAddress(&nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_);
  nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_branch = tree->GetBranch("nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip");
  if (nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_branch) nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_branch->SetAddress(&nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_);
  nisoTrack_mt2_cleaned_VVV_cutbased_veto_bak1_branch = tree->GetBranch("nisoTrack_mt2_cleaned_VVV_cutbased_veto_bak1");
  if (nisoTrack_mt2_cleaned_VVV_cutbased_veto_bak1_branch) nisoTrack_mt2_cleaned_VVV_cutbased_veto_bak1_branch->SetAddress(&nisoTrack_mt2_cleaned_VVV_cutbased_veto_bak1_);
  nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_bak1_branch = tree->GetBranch("nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_bak1");
  if (nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_bak1_branch) nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_bak1_branch->SetAddress(&nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_bak1_);
  nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_bak1_branch = tree->GetBranch("nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_bak1");
  if (nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_bak1_branch) nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_bak1_branch->SetAddress(&nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_bak1_);
  nisoTrack_mt2_cleaned_VVV_cutbased_fo_branch = tree->GetBranch("nisoTrack_mt2_cleaned_VVV_cutbased_fo");
  if (nisoTrack_mt2_cleaned_VVV_cutbased_fo_branch) nisoTrack_mt2_cleaned_VVV_cutbased_fo_branch->SetAddress(&nisoTrack_mt2_cleaned_VVV_cutbased_fo_);
  nisoTrack_mt2_cleaned_VVV_cutbased_fo_noiso_branch = tree->GetBranch("nisoTrack_mt2_cleaned_VVV_cutbased_fo_noiso");
  if (nisoTrack_mt2_cleaned_VVV_cutbased_fo_noiso_branch) nisoTrack_mt2_cleaned_VVV_cutbased_fo_noiso_branch->SetAddress(&nisoTrack_mt2_cleaned_VVV_cutbased_fo_noiso_);
  nisoTrack_mt2_cleaned_VVV_cutbased_tight_noiso_branch = tree->GetBranch("nisoTrack_mt2_cleaned_VVV_cutbased_tight_noiso");
  if (nisoTrack_mt2_cleaned_VVV_cutbased_tight_noiso_branch) nisoTrack_mt2_cleaned_VVV_cutbased_tight_noiso_branch->SetAddress(&nisoTrack_mt2_cleaned_VVV_cutbased_tight_noiso_);
  nisoTrack_mt2_cleaned_VVV_cutbased_tight_branch = tree->GetBranch("nisoTrack_mt2_cleaned_VVV_cutbased_tight");
  if (nisoTrack_mt2_cleaned_VVV_cutbased_tight_branch) nisoTrack_mt2_cleaned_VVV_cutbased_tight_branch->SetAddress(&nisoTrack_mt2_cleaned_VVV_cutbased_tight_);
  nisoTrack_mt2_cleaned_VVV_MVAbased_tight_noiso_branch = tree->GetBranch("nisoTrack_mt2_cleaned_VVV_MVAbased_tight_noiso");
  if (nisoTrack_mt2_cleaned_VVV_MVAbased_tight_noiso_branch) nisoTrack_mt2_cleaned_VVV_MVAbased_tight_noiso_branch->SetAddress(&nisoTrack_mt2_cleaned_VVV_MVAbased_tight_noiso_);
  nisoTrack_mt2_cleaned_VVV_MVAbased_tight_branch = tree->GetBranch("nisoTrack_mt2_cleaned_VVV_MVAbased_tight");
  if (nisoTrack_mt2_cleaned_VVV_MVAbased_tight_branch) nisoTrack_mt2_cleaned_VVV_MVAbased_tight_branch->SetAddress(&nisoTrack_mt2_cleaned_VVV_MVAbased_tight_);
  nisoTrack_mt2_cleaned_VVV_baseline_branch = tree->GetBranch("nisoTrack_mt2_cleaned_VVV_baseline");
  if (nisoTrack_mt2_cleaned_VVV_baseline_branch) nisoTrack_mt2_cleaned_VVV_baseline_branch->SetAddress(&nisoTrack_mt2_cleaned_VVV_baseline_);
  nisoTrack_stop_branch = tree->GetBranch("nisoTrack_stop");
  if (nisoTrack_stop_branch) nisoTrack_stop_branch->SetAddress(&nisoTrack_stop_);
  nisoTrack_PFLep5_woverlaps_branch = tree->GetBranch("nisoTrack_PFLep5_woverlaps");
  if (nisoTrack_PFLep5_woverlaps_branch) nisoTrack_PFLep5_woverlaps_branch->SetAddress(&nisoTrack_PFLep5_woverlaps_);
  nisoTrack_PFHad10_woverlaps_branch = tree->GetBranch("nisoTrack_PFHad10_woverlaps");
  if (nisoTrack_PFHad10_woverlaps_branch) nisoTrack_PFHad10_woverlaps_branch->SetAddress(&nisoTrack_PFHad10_woverlaps_);
  ngamma_branch = tree->GetBranch("ngamma");
  if (ngamma_branch) ngamma_branch->SetAddress(&ngamma_);
  gamma_pt_branch = tree->GetBranch("gamma_pt");
  if (gamma_pt_branch) gamma_pt_branch->SetAddress(&gamma_pt_);
  gamma_eta_branch = tree->GetBranch("gamma_eta");
  if (gamma_eta_branch) gamma_eta_branch->SetAddress(&gamma_eta_);
  gamma_phi_branch = tree->GetBranch("gamma_phi");
  if (gamma_phi_branch) gamma_phi_branch->SetAddress(&gamma_phi_);
  gamma_mass_branch = tree->GetBranch("gamma_mass");
  if (gamma_mass_branch) gamma_mass_branch->SetAddress(&gamma_mass_);
  gamma_mcMatchId_branch = tree->GetBranch("gamma_mcMatchId");
  if (gamma_mcMatchId_branch) gamma_mcMatchId_branch->SetAddress(&gamma_mcMatchId_);
  gamma_genPt_branch = tree->GetBranch("gamma_genPt");
  if (gamma_genPt_branch) gamma_genPt_branch->SetAddress(&gamma_genPt_);
  gamma_genIso_branch = tree->GetBranch("gamma_genIso");
  if (gamma_genIso_branch) gamma_genIso_branch->SetAddress(&gamma_genIso_);
  gamma_chHadIso_branch = tree->GetBranch("gamma_chHadIso");
  if (gamma_chHadIso_branch) gamma_chHadIso_branch->SetAddress(&gamma_chHadIso_);
  gamma_neuHadIso_branch = tree->GetBranch("gamma_neuHadIso");
  if (gamma_neuHadIso_branch) gamma_neuHadIso_branch->SetAddress(&gamma_neuHadIso_);
  gamma_phIso_branch = tree->GetBranch("gamma_phIso");
  if (gamma_phIso_branch) gamma_phIso_branch->SetAddress(&gamma_phIso_);
  gamma_sigmaIetaIeta_branch = tree->GetBranch("gamma_sigmaIetaIeta");
  if (gamma_sigmaIetaIeta_branch) gamma_sigmaIetaIeta_branch->SetAddress(&gamma_sigmaIetaIeta_);
  gamma_r9_branch = tree->GetBranch("gamma_r9");
  if (gamma_r9_branch) gamma_r9_branch->SetAddress(&gamma_r9_);
  gamma_hOverE_branch = tree->GetBranch("gamma_hOverE");
  if (gamma_hOverE_branch) gamma_hOverE_branch->SetAddress(&gamma_hOverE_);
  gamma_hOverE_online_branch = tree->GetBranch("gamma_hOverE_online");
  if (gamma_hOverE_online_branch) gamma_hOverE_online_branch->SetAddress(&gamma_hOverE_online_);
  gamma_idCutBased_branch = tree->GetBranch("gamma_idCutBased");
  if (gamma_idCutBased_branch) gamma_idCutBased_branch->SetAddress(&gamma_idCutBased_);
  gamma_ecpfclusiso_branch = tree->GetBranch("gamma_ecpfclusiso");
  if (gamma_ecpfclusiso_branch) gamma_ecpfclusiso_branch->SetAddress(&gamma_ecpfclusiso_);
  gamma_hcpfclusiso_branch = tree->GetBranch("gamma_hcpfclusiso");
  if (gamma_hcpfclusiso_branch) gamma_hcpfclusiso_branch->SetAddress(&gamma_hcpfclusiso_);
  gamma_hollowtkiso03_branch = tree->GetBranch("gamma_hollowtkiso03");
  if (gamma_hollowtkiso03_branch) gamma_hollowtkiso03_branch->SetAddress(&gamma_hollowtkiso03_);
  gamma_genIsPromptFinalState_branch = tree->GetBranch("gamma_genIsPromptFinalState");
  if (gamma_genIsPromptFinalState_branch) gamma_genIsPromptFinalState_branch->SetAddress(&gamma_genIsPromptFinalState_);
  gamma_drMinParton_branch = tree->GetBranch("gamma_drMinParton");
  if (gamma_drMinParton_branch) gamma_drMinParton_branch->SetAddress(&gamma_drMinParton_);
  ngenPart_branch = tree->GetBranch("ngenPart");
  if (ngenPart_branch) ngenPart_branch->SetAddress(&ngenPart_);
  genPart_pt_branch = tree->GetBranch("genPart_pt");
  if (genPart_pt_branch) genPart_pt_branch->SetAddress(&genPart_pt_);
  genPart_eta_branch = tree->GetBranch("genPart_eta");
  if (genPart_eta_branch) genPart_eta_branch->SetAddress(&genPart_eta_);
  genPart_phi_branch = tree->GetBranch("genPart_phi");
  if (genPart_phi_branch) genPart_phi_branch->SetAddress(&genPart_phi_);
  genPart_mass_branch = tree->GetBranch("genPart_mass");
  if (genPart_mass_branch) genPart_mass_branch->SetAddress(&genPart_mass_);
  genPart_pdgId_branch = tree->GetBranch("genPart_pdgId");
  if (genPart_pdgId_branch) genPart_pdgId_branch->SetAddress(&genPart_pdgId_);
  genPart_status_branch = tree->GetBranch("genPart_status");
  if (genPart_status_branch) genPart_status_branch->SetAddress(&genPart_status_);
  genPart_charge_branch = tree->GetBranch("genPart_charge");
  if (genPart_charge_branch) genPart_charge_branch->SetAddress(&genPart_charge_);
  genPart_motherId_branch = tree->GetBranch("genPart_motherId");
  if (genPart_motherId_branch) genPart_motherId_branch->SetAddress(&genPart_motherId_);
  genPart_grandmaId_branch = tree->GetBranch("genPart_grandmaId");
  if (genPart_grandmaId_branch) genPart_grandmaId_branch->SetAddress(&genPart_grandmaId_);
  genPart_isp6status3_branch = tree->GetBranch("genPart_isp6status3");
  if (genPart_isp6status3_branch) genPart_isp6status3_branch->SetAddress(&genPart_isp6status3_);
  ngenLep_branch = tree->GetBranch("ngenLep");
  if (ngenLep_branch) ngenLep_branch->SetAddress(&ngenLep_);
  genLep_pt_branch = tree->GetBranch("genLep_pt");
  if (genLep_pt_branch) genLep_pt_branch->SetAddress(&genLep_pt_);
  genLep_eta_branch = tree->GetBranch("genLep_eta");
  if (genLep_eta_branch) genLep_eta_branch->SetAddress(&genLep_eta_);
  genLep_phi_branch = tree->GetBranch("genLep_phi");
  if (genLep_phi_branch) genLep_phi_branch->SetAddress(&genLep_phi_);
  genLep_mass_branch = tree->GetBranch("genLep_mass");
  if (genLep_mass_branch) genLep_mass_branch->SetAddress(&genLep_mass_);
  genLep_pdgId_branch = tree->GetBranch("genLep_pdgId");
  if (genLep_pdgId_branch) genLep_pdgId_branch->SetAddress(&genLep_pdgId_);
  genLep_status_branch = tree->GetBranch("genLep_status");
  if (genLep_status_branch) genLep_status_branch->SetAddress(&genLep_status_);
  genLep_charge_branch = tree->GetBranch("genLep_charge");
  if (genLep_charge_branch) genLep_charge_branch->SetAddress(&genLep_charge_);
  genLep_sourceId_branch = tree->GetBranch("genLep_sourceId");
  if (genLep_sourceId_branch) genLep_sourceId_branch->SetAddress(&genLep_sourceId_);
  genLep_isp6status3_branch = tree->GetBranch("genLep_isp6status3");
  if (genLep_isp6status3_branch) genLep_isp6status3_branch->SetAddress(&genLep_isp6status3_);
  ngenTau_branch = tree->GetBranch("ngenTau");
  if (ngenTau_branch) ngenTau_branch->SetAddress(&ngenTau_);
  genTau_pt_branch = tree->GetBranch("genTau_pt");
  if (genTau_pt_branch) genTau_pt_branch->SetAddress(&genTau_pt_);
  genTau_eta_branch = tree->GetBranch("genTau_eta");
  if (genTau_eta_branch) genTau_eta_branch->SetAddress(&genTau_eta_);
  genTau_phi_branch = tree->GetBranch("genTau_phi");
  if (genTau_phi_branch) genTau_phi_branch->SetAddress(&genTau_phi_);
  genTau_mass_branch = tree->GetBranch("genTau_mass");
  if (genTau_mass_branch) genTau_mass_branch->SetAddress(&genTau_mass_);
  genTau_pdgId_branch = tree->GetBranch("genTau_pdgId");
  if (genTau_pdgId_branch) genTau_pdgId_branch->SetAddress(&genTau_pdgId_);
  genTau_status_branch = tree->GetBranch("genTau_status");
  if (genTau_status_branch) genTau_status_branch->SetAddress(&genTau_status_);
  genTau_charge_branch = tree->GetBranch("genTau_charge");
  if (genTau_charge_branch) genTau_charge_branch->SetAddress(&genTau_charge_);
  genTau_sourceId_branch = tree->GetBranch("genTau_sourceId");
  if (genTau_sourceId_branch) genTau_sourceId_branch->SetAddress(&genTau_sourceId_);
  ngenLepFromTau_branch = tree->GetBranch("ngenLepFromTau");
  if (ngenLepFromTau_branch) ngenLepFromTau_branch->SetAddress(&ngenLepFromTau_);
  genLepFromTau_pt_branch = tree->GetBranch("genLepFromTau_pt");
  if (genLepFromTau_pt_branch) genLepFromTau_pt_branch->SetAddress(&genLepFromTau_pt_);
  genLepFromTau_eta_branch = tree->GetBranch("genLepFromTau_eta");
  if (genLepFromTau_eta_branch) genLepFromTau_eta_branch->SetAddress(&genLepFromTau_eta_);
  genLepFromTau_phi_branch = tree->GetBranch("genLepFromTau_phi");
  if (genLepFromTau_phi_branch) genLepFromTau_phi_branch->SetAddress(&genLepFromTau_phi_);
  genLepFromTau_mass_branch = tree->GetBranch("genLepFromTau_mass");
  if (genLepFromTau_mass_branch) genLepFromTau_mass_branch->SetAddress(&genLepFromTau_mass_);
  genLepFromTau_pdgId_branch = tree->GetBranch("genLepFromTau_pdgId");
  if (genLepFromTau_pdgId_branch) genLepFromTau_pdgId_branch->SetAddress(&genLepFromTau_pdgId_);
  genLepFromTau_status_branch = tree->GetBranch("genLepFromTau_status");
  if (genLepFromTau_status_branch) genLepFromTau_status_branch->SetAddress(&genLepFromTau_status_);
  genLepFromTau_charge_branch = tree->GetBranch("genLepFromTau_charge");
  if (genLepFromTau_charge_branch) genLepFromTau_charge_branch->SetAddress(&genLepFromTau_charge_);
  genLepFromTau_sourceId_branch = tree->GetBranch("genLepFromTau_sourceId");
  if (genLepFromTau_sourceId_branch) genLepFromTau_sourceId_branch->SetAddress(&genLepFromTau_sourceId_);
  njets_branch = tree->GetBranch("njets");
  if (njets_branch) njets_branch->SetAddress(&njets_);
  removed_jets_csv_branch = tree->GetBranch("removed_jets_csv");
  if (removed_jets_csv_branch) removed_jets_csv_branch->SetAddress(&removed_jets_csv_);
  njets_up_branch = tree->GetBranch("njets_up");
  if (njets_up_branch) njets_up_branch->SetAddress(&njets_up_);
  jets_csv_branch = tree->GetBranch("jets_csv");
  if (jets_csv_branch) jets_csv_branch->SetAddress(&jets_csv_);
  jets_up_csv_branch = tree->GetBranch("jets_up_csv");
  if (jets_up_csv_branch) jets_up_csv_branch->SetAddress(&jets_up_csv_);
  njets_dn_branch = tree->GetBranch("njets_dn");
  if (njets_dn_branch) njets_dn_branch->SetAddress(&njets_dn_);
  jets_dn_csv_branch = tree->GetBranch("jets_dn_csv");
  if (jets_dn_csv_branch) jets_dn_csv_branch->SetAddress(&jets_dn_csv_);
  jets_muf_branch = tree->GetBranch("jets_muf");
  if (jets_muf_branch) jets_muf_branch->SetAddress(&jets_muf_);
  jets_mcFlavour_branch = tree->GetBranch("jets_mcFlavour");
  if (jets_mcFlavour_branch) jets_mcFlavour_branch->SetAddress(&jets_mcFlavour_);
  jets_mcHadronFlav_branch = tree->GetBranch("jets_mcHadronFlav");
  if (jets_mcHadronFlav_branch) jets_mcHadronFlav_branch->SetAddress(&jets_mcHadronFlav_);
  nhighPtPFcands_branch = tree->GetBranch("nhighPtPFcands");
  if (nhighPtPFcands_branch) nhighPtPFcands_branch->SetAddress(&nhighPtPFcands_);
  highPtPFcands_dz_branch = tree->GetBranch("highPtPFcands_dz");
  if (highPtPFcands_dz_branch) highPtPFcands_dz_branch->SetAddress(&highPtPFcands_dz_);
  highPtPFcands_pdgId_branch = tree->GetBranch("highPtPFcands_pdgId");
  if (highPtPFcands_pdgId_branch) highPtPFcands_pdgId_branch->SetAddress(&highPtPFcands_pdgId_);
  ht_branch = tree->GetBranch("ht");
  if (ht_branch) ht_branch->SetAddress(&ht_);
  ht_up_branch = tree->GetBranch("ht_up");
  if (ht_up_branch) ht_up_branch->SetAddress(&ht_up_);
  ht_dn_branch = tree->GetBranch("ht_dn");
  if (ht_dn_branch) ht_dn_branch->SetAddress(&ht_dn_);
  metsig_unofficial_branch = tree->GetBranch("metsig_unofficial");
  if (metsig_unofficial_branch) metsig_unofficial_branch->SetAddress(&metsig_unofficial_);
  metsig_unofficial_up_branch = tree->GetBranch("metsig_unofficial_up");
  if (metsig_unofficial_up_branch) metsig_unofficial_up_branch->SetAddress(&metsig_unofficial_up_);
  metsig_unofficial_dn_branch = tree->GetBranch("metsig_unofficial_dn");
  if (metsig_unofficial_dn_branch) metsig_unofficial_dn_branch->SetAddress(&metsig_unofficial_dn_);
  mt_lep1_branch = tree->GetBranch("mt_lep1");
  if (mt_lep1_branch) mt_lep1_branch->SetAddress(&mt_lep1_);
  mt2_branch = tree->GetBranch("mt2");
  if (mt2_branch) mt2_branch->SetAddress(&mt2_);
  mt2_up_branch = tree->GetBranch("mt2_up");
  if (mt2_up_branch) mt2_up_branch->SetAddress(&mt2_up_);
  mt2_dn_branch = tree->GetBranch("mt2_dn");
  if (mt2_dn_branch) mt2_dn_branch->SetAddress(&mt2_dn_);
  mt2j_branch = tree->GetBranch("mt2j");
  if (mt2j_branch) mt2j_branch->SetAddress(&mt2j_);
  mt2b_branch = tree->GetBranch("mt2b");
  if (mt2b_branch) mt2b_branch->SetAddress(&mt2b_);
  mt2b_up_branch = tree->GetBranch("mt2b_up");
  if (mt2b_up_branch) mt2b_up_branch->SetAddress(&mt2b_up_);
  mt2b_dn_branch = tree->GetBranch("mt2b_dn");
  if (mt2b_dn_branch) mt2b_dn_branch->SetAddress(&mt2b_dn_);
  mt2_genmet_branch = tree->GetBranch("mt2_genmet");
  if (mt2_genmet_branch) mt2_genmet_branch->SetAddress(&mt2_genmet_);
  mt2b_genmet_branch = tree->GetBranch("mt2b_genmet");
  if (mt2b_genmet_branch) mt2b_genmet_branch->SetAddress(&mt2b_genmet_);
  mjj_mindphi_branch = tree->GetBranch("mjj_mindphi");
  if (mjj_mindphi_branch) mjj_mindphi_branch->SetAddress(&mjj_mindphi_);
  mjj_branch = tree->GetBranch("mjj");
  if (mjj_branch) mjj_branch->SetAddress(&mjj_);
  mbb_csv_branch = tree->GetBranch("mbb_csv");
  if (mbb_csv_branch) mbb_csv_branch->SetAddress(&mbb_csv_);
  mbb_bpt_branch = tree->GetBranch("mbb_bpt");
  if (mbb_bpt_branch) mbb_bpt_branch->SetAddress(&mbb_bpt_);
  dphi_jj_branch = tree->GetBranch("dphi_jj");
  if (dphi_jj_branch) dphi_jj_branch->SetAddress(&dphi_jj_);
  dphi_ll_branch = tree->GetBranch("dphi_ll");
  if (dphi_ll_branch) dphi_ll_branch->SetAddress(&dphi_ll_);
  sum_mlb_branch = tree->GetBranch("sum_mlb");
  if (sum_mlb_branch) sum_mlb_branch->SetAddress(&sum_mlb_);
  deta_jj_branch = tree->GetBranch("deta_jj");
  if (deta_jj_branch) deta_jj_branch->SetAddress(&deta_jj_);
  dR_jj_branch = tree->GetBranch("dR_jj");
  if (dR_jj_branch) dR_jj_branch->SetAddress(&dR_jj_);
  dphi_metj1_branch = tree->GetBranch("dphi_metj1");
  if (dphi_metj1_branch) dphi_metj1_branch->SetAddress(&dphi_metj1_);
  dphi_metj2_branch = tree->GetBranch("dphi_metj2");
  if (dphi_metj2_branch) dphi_metj2_branch->SetAddress(&dphi_metj2_);
  dphi_genmetj1_branch = tree->GetBranch("dphi_genmetj1");
  if (dphi_genmetj1_branch) dphi_genmetj1_branch->SetAddress(&dphi_genmetj1_);
  dphi_genmetj2_branch = tree->GetBranch("dphi_genmetj2");
  if (dphi_genmetj2_branch) dphi_genmetj2_branch->SetAddress(&dphi_genmetj2_);
  mjj_mindphi_up_branch = tree->GetBranch("mjj_mindphi_up");
  if (mjj_mindphi_up_branch) mjj_mindphi_up_branch->SetAddress(&mjj_mindphi_up_);
  mjj_up_branch = tree->GetBranch("mjj_up");
  if (mjj_up_branch) mjj_up_branch->SetAddress(&mjj_up_);
  mbb_csv_up_branch = tree->GetBranch("mbb_csv_up");
  if (mbb_csv_up_branch) mbb_csv_up_branch->SetAddress(&mbb_csv_up_);
  mbb_bpt_up_branch = tree->GetBranch("mbb_bpt_up");
  if (mbb_bpt_up_branch) mbb_bpt_up_branch->SetAddress(&mbb_bpt_up_);
  dphi_jj_up_branch = tree->GetBranch("dphi_jj_up");
  if (dphi_jj_up_branch) dphi_jj_up_branch->SetAddress(&dphi_jj_up_);
  dphi_ll_up_branch = tree->GetBranch("dphi_ll_up");
  if (dphi_ll_up_branch) dphi_ll_up_branch->SetAddress(&dphi_ll_up_);
  sum_mlb_up_branch = tree->GetBranch("sum_mlb_up");
  if (sum_mlb_up_branch) sum_mlb_up_branch->SetAddress(&sum_mlb_up_);
  deta_jj_up_branch = tree->GetBranch("deta_jj_up");
  if (deta_jj_up_branch) deta_jj_up_branch->SetAddress(&deta_jj_up_);
  dR_jj_up_branch = tree->GetBranch("dR_jj_up");
  if (dR_jj_up_branch) dR_jj_up_branch->SetAddress(&dR_jj_up_);
  dphi_metj1_up_branch = tree->GetBranch("dphi_metj1_up");
  if (dphi_metj1_up_branch) dphi_metj1_up_branch->SetAddress(&dphi_metj1_up_);
  dphi_metj2_up_branch = tree->GetBranch("dphi_metj2_up");
  if (dphi_metj2_up_branch) dphi_metj2_up_branch->SetAddress(&dphi_metj2_up_);
  mjj_mindphi_dn_branch = tree->GetBranch("mjj_mindphi_dn");
  if (mjj_mindphi_dn_branch) mjj_mindphi_dn_branch->SetAddress(&mjj_mindphi_dn_);
  mjj_dn_branch = tree->GetBranch("mjj_dn");
  if (mjj_dn_branch) mjj_dn_branch->SetAddress(&mjj_dn_);
  mbb_csv_dn_branch = tree->GetBranch("mbb_csv_dn");
  if (mbb_csv_dn_branch) mbb_csv_dn_branch->SetAddress(&mbb_csv_dn_);
  mbb_bpt_dn_branch = tree->GetBranch("mbb_bpt_dn");
  if (mbb_bpt_dn_branch) mbb_bpt_dn_branch->SetAddress(&mbb_bpt_dn_);
  dphi_jj_dn_branch = tree->GetBranch("dphi_jj_dn");
  if (dphi_jj_dn_branch) dphi_jj_dn_branch->SetAddress(&dphi_jj_dn_);
  dphi_ll_dn_branch = tree->GetBranch("dphi_ll_dn");
  if (dphi_ll_dn_branch) dphi_ll_dn_branch->SetAddress(&dphi_ll_dn_);
  sum_mlb_dn_branch = tree->GetBranch("sum_mlb_dn");
  if (sum_mlb_dn_branch) sum_mlb_dn_branch->SetAddress(&sum_mlb_dn_);
  deta_jj_dn_branch = tree->GetBranch("deta_jj_dn");
  if (deta_jj_dn_branch) deta_jj_dn_branch->SetAddress(&deta_jj_dn_);
  dR_jj_dn_branch = tree->GetBranch("dR_jj_dn");
  if (dR_jj_dn_branch) dR_jj_dn_branch->SetAddress(&dR_jj_dn_);
  dphi_metj1_dn_branch = tree->GetBranch("dphi_metj1_dn");
  if (dphi_metj1_dn_branch) dphi_metj1_dn_branch->SetAddress(&dphi_metj1_dn_);
  dphi_metj2_dn_branch = tree->GetBranch("dphi_metj2_dn");
  if (dphi_metj2_dn_branch) dphi_metj2_dn_branch->SetAddress(&dphi_metj2_dn_);
  weight_btagsf_branch = tree->GetBranch("weight_btagsf");
  if (weight_btagsf_branch) weight_btagsf_branch->SetAddress(&weight_btagsf_);
  weight_btagsf_heavy_UP_branch = tree->GetBranch("weight_btagsf_heavy_UP");
  if (weight_btagsf_heavy_UP_branch) weight_btagsf_heavy_UP_branch->SetAddress(&weight_btagsf_heavy_UP_);
  weight_btagsf_light_UP_branch = tree->GetBranch("weight_btagsf_light_UP");
  if (weight_btagsf_light_UP_branch) weight_btagsf_light_UP_branch->SetAddress(&weight_btagsf_light_UP_);
  weight_btagsf_heavy_DN_branch = tree->GetBranch("weight_btagsf_heavy_DN");
  if (weight_btagsf_heavy_DN_branch) weight_btagsf_heavy_DN_branch->SetAddress(&weight_btagsf_heavy_DN_);
  weight_btagsf_light_DN_branch = tree->GetBranch("weight_btagsf_light_DN");
  if (weight_btagsf_light_DN_branch) weight_btagsf_light_DN_branch->SetAddress(&weight_btagsf_light_DN_);
  chpfcands_0013_pt_branch = tree->GetBranch("chpfcands_0013_pt");
  if (chpfcands_0013_pt_branch) chpfcands_0013_pt_branch->SetAddress(&chpfcands_0013_pt_);
  chpfcands_1316_pt_branch = tree->GetBranch("chpfcands_1316_pt");
  if (chpfcands_1316_pt_branch) chpfcands_1316_pt_branch->SetAddress(&chpfcands_1316_pt_);
  chpfcands_1624_pt_branch = tree->GetBranch("chpfcands_1624_pt");
  if (chpfcands_1624_pt_branch) chpfcands_1624_pt_branch->SetAddress(&chpfcands_1624_pt_);
  chpfcands_2430_pt_branch = tree->GetBranch("chpfcands_2430_pt");
  if (chpfcands_2430_pt_branch) chpfcands_2430_pt_branch->SetAddress(&chpfcands_2430_pt_);
  chpfcands_30in_pt_branch = tree->GetBranch("chpfcands_30in_pt");
  if (chpfcands_30in_pt_branch) chpfcands_30in_pt_branch->SetAddress(&chpfcands_30in_pt_);
  phpfcands_0013_pt_branch = tree->GetBranch("phpfcands_0013_pt");
  if (phpfcands_0013_pt_branch) phpfcands_0013_pt_branch->SetAddress(&phpfcands_0013_pt_);
  phpfcands_1316_pt_branch = tree->GetBranch("phpfcands_1316_pt");
  if (phpfcands_1316_pt_branch) phpfcands_1316_pt_branch->SetAddress(&phpfcands_1316_pt_);
  phpfcands_1624_pt_branch = tree->GetBranch("phpfcands_1624_pt");
  if (phpfcands_1624_pt_branch) phpfcands_1624_pt_branch->SetAddress(&phpfcands_1624_pt_);
  phpfcands_2430_pt_branch = tree->GetBranch("phpfcands_2430_pt");
  if (phpfcands_2430_pt_branch) phpfcands_2430_pt_branch->SetAddress(&phpfcands_2430_pt_);
  phpfcands_30in_pt_branch = tree->GetBranch("phpfcands_30in_pt");
  if (phpfcands_30in_pt_branch) phpfcands_30in_pt_branch->SetAddress(&phpfcands_30in_pt_);
  nupfcands_0013_pt_branch = tree->GetBranch("nupfcands_0013_pt");
  if (nupfcands_0013_pt_branch) nupfcands_0013_pt_branch->SetAddress(&nupfcands_0013_pt_);
  nupfcands_1316_pt_branch = tree->GetBranch("nupfcands_1316_pt");
  if (nupfcands_1316_pt_branch) nupfcands_1316_pt_branch->SetAddress(&nupfcands_1316_pt_);
  nupfcands_1624_pt_branch = tree->GetBranch("nupfcands_1624_pt");
  if (nupfcands_1624_pt_branch) nupfcands_1624_pt_branch->SetAddress(&nupfcands_1624_pt_);
  nupfcands_2430_pt_branch = tree->GetBranch("nupfcands_2430_pt");
  if (nupfcands_2430_pt_branch) nupfcands_2430_pt_branch->SetAddress(&nupfcands_2430_pt_);
  nupfcands_30in_pt_branch = tree->GetBranch("nupfcands_30in_pt");
  if (nupfcands_30in_pt_branch) nupfcands_30in_pt_branch->SetAddress(&nupfcands_30in_pt_);
  chpfcands_0013_sumet_branch = tree->GetBranch("chpfcands_0013_sumet");
  if (chpfcands_0013_sumet_branch) chpfcands_0013_sumet_branch->SetAddress(&chpfcands_0013_sumet_);
  chpfcands_1316_sumet_branch = tree->GetBranch("chpfcands_1316_sumet");
  if (chpfcands_1316_sumet_branch) chpfcands_1316_sumet_branch->SetAddress(&chpfcands_1316_sumet_);
  chpfcands_1624_sumet_branch = tree->GetBranch("chpfcands_1624_sumet");
  if (chpfcands_1624_sumet_branch) chpfcands_1624_sumet_branch->SetAddress(&chpfcands_1624_sumet_);
  chpfcands_2430_sumet_branch = tree->GetBranch("chpfcands_2430_sumet");
  if (chpfcands_2430_sumet_branch) chpfcands_2430_sumet_branch->SetAddress(&chpfcands_2430_sumet_);
  chpfcands_30in_sumet_branch = tree->GetBranch("chpfcands_30in_sumet");
  if (chpfcands_30in_sumet_branch) chpfcands_30in_sumet_branch->SetAddress(&chpfcands_30in_sumet_);
  phpfcands_0013_sumet_branch = tree->GetBranch("phpfcands_0013_sumet");
  if (phpfcands_0013_sumet_branch) phpfcands_0013_sumet_branch->SetAddress(&phpfcands_0013_sumet_);
  phpfcands_1316_sumet_branch = tree->GetBranch("phpfcands_1316_sumet");
  if (phpfcands_1316_sumet_branch) phpfcands_1316_sumet_branch->SetAddress(&phpfcands_1316_sumet_);
  phpfcands_1624_sumet_branch = tree->GetBranch("phpfcands_1624_sumet");
  if (phpfcands_1624_sumet_branch) phpfcands_1624_sumet_branch->SetAddress(&phpfcands_1624_sumet_);
  phpfcands_2430_sumet_branch = tree->GetBranch("phpfcands_2430_sumet");
  if (phpfcands_2430_sumet_branch) phpfcands_2430_sumet_branch->SetAddress(&phpfcands_2430_sumet_);
  phpfcands_30in_sumet_branch = tree->GetBranch("phpfcands_30in_sumet");
  if (phpfcands_30in_sumet_branch) phpfcands_30in_sumet_branch->SetAddress(&phpfcands_30in_sumet_);
  nupfcands_0013_sumet_branch = tree->GetBranch("nupfcands_0013_sumet");
  if (nupfcands_0013_sumet_branch) nupfcands_0013_sumet_branch->SetAddress(&nupfcands_0013_sumet_);
  nupfcands_1316_sumet_branch = tree->GetBranch("nupfcands_1316_sumet");
  if (nupfcands_1316_sumet_branch) nupfcands_1316_sumet_branch->SetAddress(&nupfcands_1316_sumet_);
  nupfcands_1624_sumet_branch = tree->GetBranch("nupfcands_1624_sumet");
  if (nupfcands_1624_sumet_branch) nupfcands_1624_sumet_branch->SetAddress(&nupfcands_1624_sumet_);
  nupfcands_2430_sumet_branch = tree->GetBranch("nupfcands_2430_sumet");
  if (nupfcands_2430_sumet_branch) nupfcands_2430_sumet_branch->SetAddress(&nupfcands_2430_sumet_);
  nupfcands_30in_sumet_branch = tree->GetBranch("nupfcands_30in_sumet");
  if (nupfcands_30in_sumet_branch) nupfcands_30in_sumet_branch->SetAddress(&nupfcands_30in_sumet_);
  chpfcands_0013_phi_branch = tree->GetBranch("chpfcands_0013_phi");
  if (chpfcands_0013_phi_branch) chpfcands_0013_phi_branch->SetAddress(&chpfcands_0013_phi_);
  chpfcands_1316_phi_branch = tree->GetBranch("chpfcands_1316_phi");
  if (chpfcands_1316_phi_branch) chpfcands_1316_phi_branch->SetAddress(&chpfcands_1316_phi_);
  chpfcands_1624_phi_branch = tree->GetBranch("chpfcands_1624_phi");
  if (chpfcands_1624_phi_branch) chpfcands_1624_phi_branch->SetAddress(&chpfcands_1624_phi_);
  chpfcands_2430_phi_branch = tree->GetBranch("chpfcands_2430_phi");
  if (chpfcands_2430_phi_branch) chpfcands_2430_phi_branch->SetAddress(&chpfcands_2430_phi_);
  chpfcands_30in_phi_branch = tree->GetBranch("chpfcands_30in_phi");
  if (chpfcands_30in_phi_branch) chpfcands_30in_phi_branch->SetAddress(&chpfcands_30in_phi_);
  phpfcands_0013_phi_branch = tree->GetBranch("phpfcands_0013_phi");
  if (phpfcands_0013_phi_branch) phpfcands_0013_phi_branch->SetAddress(&phpfcands_0013_phi_);
  phpfcands_1316_phi_branch = tree->GetBranch("phpfcands_1316_phi");
  if (phpfcands_1316_phi_branch) phpfcands_1316_phi_branch->SetAddress(&phpfcands_1316_phi_);
  phpfcands_1624_phi_branch = tree->GetBranch("phpfcands_1624_phi");
  if (phpfcands_1624_phi_branch) phpfcands_1624_phi_branch->SetAddress(&phpfcands_1624_phi_);
  phpfcands_2430_phi_branch = tree->GetBranch("phpfcands_2430_phi");
  if (phpfcands_2430_phi_branch) phpfcands_2430_phi_branch->SetAddress(&phpfcands_2430_phi_);
  phpfcands_30in_phi_branch = tree->GetBranch("phpfcands_30in_phi");
  if (phpfcands_30in_phi_branch) phpfcands_30in_phi_branch->SetAddress(&phpfcands_30in_phi_);
  nupfcands_0013_phi_branch = tree->GetBranch("nupfcands_0013_phi");
  if (nupfcands_0013_phi_branch) nupfcands_0013_phi_branch->SetAddress(&nupfcands_0013_phi_);
  nupfcands_1316_phi_branch = tree->GetBranch("nupfcands_1316_phi");
  if (nupfcands_1316_phi_branch) nupfcands_1316_phi_branch->SetAddress(&nupfcands_1316_phi_);
  nupfcands_1624_phi_branch = tree->GetBranch("nupfcands_1624_phi");
  if (nupfcands_1624_phi_branch) nupfcands_1624_phi_branch->SetAddress(&nupfcands_1624_phi_);
  nupfcands_2430_phi_branch = tree->GetBranch("nupfcands_2430_phi");
  if (nupfcands_2430_phi_branch) nupfcands_2430_phi_branch->SetAddress(&nupfcands_2430_phi_);
  nupfcands_30in_phi_branch = tree->GetBranch("nupfcands_30in_phi");
  if (nupfcands_30in_phi_branch) nupfcands_30in_phi_branch->SetAddress(&nupfcands_30in_phi_);
  met_T1CHS_pt_branch = tree->GetBranch("met_T1CHS_pt");
  if (met_T1CHS_pt_branch) met_T1CHS_pt_branch->SetAddress(&met_T1CHS_pt_);
  met_T1CHS_phi_branch = tree->GetBranch("met_T1CHS_phi");
  if (met_T1CHS_phi_branch) met_T1CHS_phi_branch->SetAddress(&met_T1CHS_phi_);
  met_T1CHS_fromCORE_pt_branch = tree->GetBranch("met_T1CHS_fromCORE_pt");
  if (met_T1CHS_fromCORE_pt_branch) met_T1CHS_fromCORE_pt_branch->SetAddress(&met_T1CHS_fromCORE_pt_);
  met_T1CHS_fromCORE_phi_branch = tree->GetBranch("met_T1CHS_fromCORE_phi");
  if (met_T1CHS_fromCORE_phi_branch) met_T1CHS_fromCORE_phi_branch->SetAddress(&met_T1CHS_fromCORE_phi_);
  met_T1CHS_miniAOD_CORE_pt_branch = tree->GetBranch("met_T1CHS_miniAOD_CORE_pt");
  if (met_T1CHS_miniAOD_CORE_pt_branch) met_T1CHS_miniAOD_CORE_pt_branch->SetAddress(&met_T1CHS_miniAOD_CORE_pt_);
  met_T1CHS_miniAOD_CORE_phi_branch = tree->GetBranch("met_T1CHS_miniAOD_CORE_phi");
  if (met_T1CHS_miniAOD_CORE_phi_branch) met_T1CHS_miniAOD_CORE_phi_branch->SetAddress(&met_T1CHS_miniAOD_CORE_phi_);
  met_T1CHS_miniAOD_CORE_up_pt_branch = tree->GetBranch("met_T1CHS_miniAOD_CORE_up_pt");
  if (met_T1CHS_miniAOD_CORE_up_pt_branch) met_T1CHS_miniAOD_CORE_up_pt_branch->SetAddress(&met_T1CHS_miniAOD_CORE_up_pt_);
  met_T1CHS_miniAOD_CORE_up_phi_branch = tree->GetBranch("met_T1CHS_miniAOD_CORE_up_phi");
  if (met_T1CHS_miniAOD_CORE_up_phi_branch) met_T1CHS_miniAOD_CORE_up_phi_branch->SetAddress(&met_T1CHS_miniAOD_CORE_up_phi_);
  met_T1CHS_miniAOD_CORE_dn_pt_branch = tree->GetBranch("met_T1CHS_miniAOD_CORE_dn_pt");
  if (met_T1CHS_miniAOD_CORE_dn_pt_branch) met_T1CHS_miniAOD_CORE_dn_pt_branch->SetAddress(&met_T1CHS_miniAOD_CORE_dn_pt_);
  met_T1CHS_miniAOD_CORE_dn_phi_branch = tree->GetBranch("met_T1CHS_miniAOD_CORE_dn_phi");
  if (met_T1CHS_miniAOD_CORE_dn_phi_branch) met_T1CHS_miniAOD_CORE_dn_phi_branch->SetAddress(&met_T1CHS_miniAOD_CORE_dn_phi_);
  decayedphoton_mt2_branch = tree->GetBranch("decayedphoton_mt2");
  if (decayedphoton_mt2_branch) decayedphoton_mt2_branch->SetAddress(&decayedphoton_mt2_);
  hyp_type_branch = tree->GetBranch("hyp_type");
  if (hyp_type_branch) hyp_type_branch->SetAddress(&hyp_type_);
  evt_type_branch = tree->GetBranch("evt_type");
  if (evt_type_branch) evt_type_branch->SetAddress(&evt_type_);
  mass_gluino_branch = tree->GetBranch("mass_gluino");
  if (mass_gluino_branch) mass_gluino_branch->SetAddress(&mass_gluino_);
  mass_LSP_branch = tree->GetBranch("mass_LSP");
  if (mass_LSP_branch) mass_LSP_branch->SetAddress(&mass_LSP_);
  mass_chi_branch = tree->GetBranch("mass_chi");
  if (mass_chi_branch) mass_chi_branch->SetAddress(&mass_chi_);
  isrboost_branch = tree->GetBranch("isrboost");
  if (isrboost_branch) isrboost_branch->SetAddress(&isrboost_);
  isr_njets_branch = tree->GetBranch("isr_njets");
  if (isr_njets_branch) isr_njets_branch->SetAddress(&isr_njets_);
  isr_weight_branch = tree->GetBranch("isr_weight");
  if (isr_weight_branch) isr_weight_branch->SetAddress(&isr_weight_);
  isr_unc_branch = tree->GetBranch("isr_unc");
  if (isr_unc_branch) isr_unc_branch->SetAddress(&isr_unc_);
  weightsf_lepid_branch = tree->GetBranch("weightsf_lepid");
  if (weightsf_lepid_branch) weightsf_lepid_branch->SetAddress(&weightsf_lepid_);
  weightsf_lepiso_branch = tree->GetBranch("weightsf_lepiso");
  if (weightsf_lepiso_branch) weightsf_lepiso_branch->SetAddress(&weightsf_lepiso_);
  weightsf_lepip_branch = tree->GetBranch("weightsf_lepip");
  if (weightsf_lepip_branch) weightsf_lepip_branch->SetAddress(&weightsf_lepip_);
  weightsf_lepreco_branch = tree->GetBranch("weightsf_lepreco");
  if (weightsf_lepreco_branch) weightsf_lepreco_branch->SetAddress(&weightsf_lepreco_);
  weightsf_lepconv_branch = tree->GetBranch("weightsf_lepconv");
  if (weightsf_lepconv_branch) weightsf_lepconv_branch->SetAddress(&weightsf_lepconv_);
  weightsf_lepid_FS_branch = tree->GetBranch("weightsf_lepid_FS");
  if (weightsf_lepid_FS_branch) weightsf_lepid_FS_branch->SetAddress(&weightsf_lepid_FS_);
  weightsf_lepiso_FS_branch = tree->GetBranch("weightsf_lepiso_FS");
  if (weightsf_lepiso_FS_branch) weightsf_lepiso_FS_branch->SetAddress(&weightsf_lepiso_FS_);
  weightsf_lepip_FS_branch = tree->GetBranch("weightsf_lepip_FS");
  if (weightsf_lepip_FS_branch) weightsf_lepip_FS_branch->SetAddress(&weightsf_lepip_FS_);
  weight_rn_r1_n1_branch = tree->GetBranch("weight_rn_r1_n1");
  if (weight_rn_r1_n1_branch) weight_rn_r1_n1_branch->SetAddress(&weight_rn_r1_n1_);
  weight_rn_r1_n2_branch = tree->GetBranch("weight_rn_r1_n2");
  if (weight_rn_r1_n2_branch) weight_rn_r1_n2_branch->SetAddress(&weight_rn_r1_n2_);
  weight_rn_r1_np05_branch = tree->GetBranch("weight_rn_r1_np05");
  if (weight_rn_r1_np05_branch) weight_rn_r1_np05_branch->SetAddress(&weight_rn_r1_np05_);
  weight_rn_r2_n1_branch = tree->GetBranch("weight_rn_r2_n1");
  if (weight_rn_r2_n1_branch) weight_rn_r2_n1_branch->SetAddress(&weight_rn_r2_n1_);
  weight_rn_r2_n2_branch = tree->GetBranch("weight_rn_r2_n2");
  if (weight_rn_r2_n2_branch) weight_rn_r2_n2_branch->SetAddress(&weight_rn_r2_n2_);
  weight_rn_r2_n0p5_branch = tree->GetBranch("weight_rn_r2_n0p5");
  if (weight_rn_r2_n0p5_branch) weight_rn_r2_n0p5_branch->SetAddress(&weight_rn_r2_n0p5_);
  weight_rn_r0p5_n1_branch = tree->GetBranch("weight_rn_r0p5_n1");
  if (weight_rn_r0p5_n1_branch) weight_rn_r0p5_n1_branch->SetAddress(&weight_rn_r0p5_n1_);
  weight_rn_r0p5_n2_branch = tree->GetBranch("weight_rn_r0p5_n2");
  if (weight_rn_r0p5_n2_branch) weight_rn_r0p5_n2_branch->SetAddress(&weight_rn_r0p5_n2_);
  weight_rn_r0p5_n0p5_branch = tree->GetBranch("weight_rn_r0p5_n0p5");
  if (weight_rn_r0p5_n0p5_branch) weight_rn_r0p5_n0p5_branch->SetAddress(&weight_rn_r0p5_n0p5_);
  weight_pdf_up_branch = tree->GetBranch("weight_pdf_up");
  if (weight_pdf_up_branch) weight_pdf_up_branch->SetAddress(&weight_pdf_up_);
  weight_pdf_down_branch = tree->GetBranch("weight_pdf_down");
  if (weight_pdf_down_branch) weight_pdf_down_branch->SetAddress(&weight_pdf_down_);
  weight_alphas_down_branch = tree->GetBranch("weight_alphas_down");
  if (weight_alphas_down_branch) weight_alphas_down_branch->SetAddress(&weight_alphas_down_);
  weight_alphas_up_branch = tree->GetBranch("weight_alphas_up");
  if (weight_alphas_up_branch) weight_alphas_up_branch->SetAddress(&weight_alphas_up_);
  LooseElec_branch = tree->GetBranch("LooseElec");
  if (LooseElec_branch) LooseElec_branch->SetAddress(&LooseElec_);
  LooseLepton_branch = tree->GetBranch("LooseLepton");
  if (LooseLepton_branch) LooseLepton_branch->SetAddress(&LooseLepton_);
  LooseMuon_branch = tree->GetBranch("LooseMuon");
  if (LooseMuon_branch) LooseMuon_branch->SetAddress(&LooseMuon_);
  SignalElec_branch = tree->GetBranch("SignalElec");
  if (SignalElec_branch) SignalElec_branch->SetAddress(&SignalElec_);
  SignalLepton_branch = tree->GetBranch("SignalLepton");
  if (SignalLepton_branch) SignalLepton_branch->SetAddress(&SignalLepton_);
  SignalMuon_branch = tree->GetBranch("SignalMuon");
  if (SignalMuon_branch) SignalMuon_branch->SetAddress(&SignalMuon_);
  Tight3lElec_branch = tree->GetBranch("Tight3lElec");
  if (Tight3lElec_branch) Tight3lElec_branch->SetAddress(&Tight3lElec_);
  Tight3lLepton_branch = tree->GetBranch("Tight3lLepton");
  if (Tight3lLepton_branch) Tight3lLepton_branch->SetAddress(&Tight3lLepton_);
  TightElec_branch = tree->GetBranch("TightElec");
  if (TightElec_branch) TightElec_branch->SetAddress(&TightElec_);
  TightLepton_branch = tree->GetBranch("TightLepton");
  if (TightLepton_branch) TightLepton_branch->SetAddress(&TightLepton_);
  TightMuon_branch = tree->GetBranch("TightMuon");
  if (TightMuon_branch) TightMuon_branch->SetAddress(&TightMuon_);
  VetoElec_branch = tree->GetBranch("VetoElec");
  if (VetoElec_branch) VetoElec_branch->SetAddress(&VetoElec_);
  VetoLepton_branch = tree->GetBranch("VetoLepton");
  if (VetoLepton_branch) VetoLepton_branch->SetAddress(&VetoLepton_);
  VetoMuon_branch = tree->GetBranch("VetoMuon");
  if (VetoMuon_branch) VetoMuon_branch->SetAddress(&VetoMuon_);
  Good3LJet_branch = tree->GetBranch("Good3LJet");
  if (Good3LJet_branch) Good3LJet_branch->SetAddress(&Good3LJet_);
  GoodSSJet_branch = tree->GetBranch("GoodSSJet");
  if (GoodSSJet_branch) GoodSSJet_branch->SetAddress(&GoodSSJet_);
  GoodSSWJet_branch = tree->GetBranch("GoodSSWJet");
  if (GoodSSWJet_branch) GoodSSWJet_branch->SetAddress(&GoodSSWJet_);
  LooseBJet_branch = tree->GetBranch("LooseBJet");
  if (LooseBJet_branch) LooseBJet_branch->SetAddress(&LooseBJet_);

  tree->SetMakeClass(0);
}

void WWWTree::GetEntry(unsigned int idx) {
  // this only marks branches as not loaded, saving a lot of time
  index = idx;
  run_isLoaded = false;
  lumi_isLoaded = false;
  evt_isLoaded = false;
  isData_isLoaded = false;
  evt_passgoodrunlist_isLoaded = false;
  evt_scale1fb_isLoaded = false;
  evt_xsec_isLoaded = false;
  evt_kfactor_isLoaded = false;
  evt_filter_isLoaded = false;
  evt_nEvts_isLoaded = false;
  evt_dataset_isLoaded = false;
  puWeight_isLoaded = false;
  nVert_isLoaded = false;
  firstgoodvertex_isLoaded = false;
  nTrueInt_isLoaded = false;
  rho_isLoaded = false;
  rho25_isLoaded = false;
  gen_ht_isLoaded = false;
  nBJetTight_isLoaded = false;
  nBJetMedium_isLoaded = false;
  nBJetLoose_isLoaded = false;
  nBJetTight_up_isLoaded = false;
  nBJetMedium_up_isLoaded = false;
  nBJetLoose_up_isLoaded = false;
  nBJetTight_dn_isLoaded = false;
  nBJetMedium_dn_isLoaded = false;
  nBJetLoose_dn_isLoaded = false;
  nJet200MuFrac50DphiMet_isLoaded = false;
  nMuons10_isLoaded = false;
  nBadMuons20_isLoaded = false;
  nElectrons10_isLoaded = false;
  nGammas20_isLoaded = false;
  nTaus20_isLoaded = false;
  met_pt_isLoaded = false;
  met_phi_isLoaded = false;
  met_calo_pt_isLoaded = false;
  met_calo_phi_isLoaded = false;
  met_miniaod_pt_isLoaded = false;
  met_miniaod_phi_isLoaded = false;
  met_muegclean_pt_isLoaded = false;
  met_muegclean_phi_isLoaded = false;
  met_rawPt_isLoaded = false;
  met_rawPhi_isLoaded = false;
  met_genPt_isLoaded = false;
  met_genPhi_isLoaded = false;
  sumet_raw_isLoaded = false;
  Flag_ecalLaserCorrFilter_isLoaded = false;
  Flag_hcalLaserEventFilter_isLoaded = false;
  Flag_trackingFailureFilter_isLoaded = false;
  Flag_CSCTightHaloFilter_isLoaded = false;
  Flag_HBHENoiseFilter_isLoaded = false;
  Flag_HBHEIsoNoiseFilter_isLoaded = false;
  Flag_CSCTightHalo2015Filter_isLoaded = false;
  Flag_EcalDeadCellTriggerPrimitiveFilter_isLoaded = false;
  Flag_goodVertices_isLoaded = false;
  Flag_eeBadScFilter_isLoaded = false;
  Flag_globalTightHalo2016_isLoaded = false;
  Flag_badMuonFilter_isLoaded = false;
  Flag_badChargedCandidateFilter_isLoaded = false;
  Flag_badMuonFilterv2_isLoaded = false;
  Flag_badChargedCandidateFilterv2_isLoaded = false;
  Flag_badMuons_isLoaded = false;
  Flag_duplicateMuons_isLoaded = false;
  Flag_noBadMuons_isLoaded = false;
  HLT_singleEl_isLoaded = false;
  HLT_singleMu_isLoaded = false;
  HLT_singleMu_noiso_isLoaded = false;
  HLT_DoubleEl_noiso_isLoaded = false;
  HLT_DoubleEl_isLoaded = false;
  HLT_DoubleEl_DZ_isLoaded = false;
  HLT_DoubleEl_DZ_2_isLoaded = false;
  HLT_MuEG_isLoaded = false;
  HLT_MuEG_2_isLoaded = false;
  HLT_MuEG_noiso_isLoaded = false;
  HLT_MuEG_noiso_2_isLoaded = false;
  HLT_Mu8_EG17_isLoaded = false;
  HLT_Mu8_EG23_isLoaded = false;
  HLT_Mu8_EG23_DZ_isLoaded = false;
  HLT_Mu12_EG23_DZ_isLoaded = false;
  HLT_Mu17_EG12_isLoaded = false;
  HLT_Mu23_EG8_isLoaded = false;
  HLT_Mu23_EG8_DZ_isLoaded = false;
  HLT_Mu23_EG12_isLoaded = false;
  HLT_Mu23_EG12_DZ_isLoaded = false;
  HLT_DoubleMu_noiso_isLoaded = false;
  HLT_DoubleMu_noiso_27_8_isLoaded = false;
  HLT_DoubleMu_noiso_30_11_isLoaded = false;
  HLT_DoubleMu_noiso_40_11_isLoaded = false;
  HLT_DoubleMu_isLoaded = false;
  HLT_DoubleMu_tk_isLoaded = false;
  HLT_DoubleMu_dbltk_isLoaded = false;
  HLT_DoubleMu_nonDZ_isLoaded = false;
  HLT_DoubleMu_tk_nonDZ_isLoaded = false;
  HLT_Photon22_R9Id90_HE10_IsoM_isLoaded = false;
  HLT_Photon30_R9Id90_HE10_IsoM_isLoaded = false;
  HLT_Photon36_R9Id90_HE10_IsoM_isLoaded = false;
  HLT_Photon50_R9Id90_HE10_IsoM_isLoaded = false;
  HLT_Photon75_R9Id90_HE10_IsoM_isLoaded = false;
  HLT_Photon90_R9Id90_HE10_IsoM_isLoaded = false;
  HLT_Photon120_R9Id90_HE10_IsoM_isLoaded = false;
  HLT_Photon165_R9Id90_HE10_IsoM_isLoaded = false;
  HLT_Photon165_HE10_isLoaded = false;
  HLT_CaloJet500_NoJetID_isLoaded = false;
  HLT_ECALHT800_NoJetID_isLoaded = false;
  HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton_isLoaded = false;
  HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton_isLoaded = false;
  HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton_isLoaded = false;
  HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton_isLoaded = false;
  HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton_isLoaded = false;
  HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton_isLoaded = false;
  HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton_isLoaded = false;
  HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton_isLoaded = false;
  HLT_Photon165_HE10_matchedtophoton_isLoaded = false;
  dilmass_isLoaded = false;
  dilpt_isLoaded = false;
  dRll_isLoaded = false;
  matched_neutralemf_isLoaded = false;
  matched_emf_isLoaded = false;
  elveto_isLoaded = false;
  nlep_isLoaded = false;
  nveto_leptons_isLoaded = false;
  nVetoEl_relIso03EAless01_isLoaded = false;
  nVetoEl_relIso03EAless02_isLoaded = false;
  nVetoEl_relIso03EAless03_isLoaded = false;
  nVetoEl_relIso03EAless04_isLoaded = false;
  nVetoMu_relIso03EAless01_isLoaded = false;
  nVetoMu_relIso03EAless02_isLoaded = false;
  nVetoMu_relIso03EAless03_isLoaded = false;
  nVetoMu_relIso03EAless04_isLoaded = false;
  lep_p4_isLoaded = false;
  lep_pt_isLoaded = false;
  lep_eta_isLoaded = false;
  lep_phi_isLoaded = false;
  lep_mass_isLoaded = false;
  lep_charge_isLoaded = false;
  lep_3ch_agree_isLoaded = false;
  lep_isFromW_isLoaded = false;
  lep_isFromZ_isLoaded = false;
  lep_isFromB_isLoaded = false;
  lep_isFromC_isLoaded = false;
  lep_isFromL_isLoaded = false;
  lep_isFromLF_isLoaded = false;
  lep_closest_jet_p4_isLoaded = false;
  lep_ptRatio_isLoaded = false;
  lep_motherIdSS_isLoaded = false;
  lep_genPart_index_isLoaded = false;
  lep_coneCorrPt_isLoaded = false;
  lep_ptRel_isLoaded = false;
  lep_relIso03_isLoaded = false;
  lep_relIso03DB_isLoaded = false;
  lep_relIso03EA_isLoaded = false;
  lep_relIso03EAv2_isLoaded = false;
  lep_relIso04DB_isLoaded = false;
  lep_relIso04EA_isLoaded = false;
  lep_relIso04EAv2_isLoaded = false;
  lep_miniRelIsoCMS3_EA_isLoaded = false;
  lep_miniRelIsoCMS3_EAv2_isLoaded = false;
  lep_miniRelIsoCMS3_DB_isLoaded = false;
  lep_isTriggerSafenoIso_v1_isLoaded = false;
  lep_isTriggerSafenoIso_v2_isLoaded = false;
  lep_isTriggerSafe_v1_isLoaded = false;
  lep_isTriggerSafe_v2_isLoaded = false;
  lep_pass_VVV_cutbased_veto_isLoaded = false;
  lep_pass_VVV_cutbased_veto_noiso_isLoaded = false;
  lep_pass_VVV_cutbased_veto_noiso_noip_isLoaded = false;
  lep_pass_VVV_cutbased_veto_bak1_isLoaded = false;
  lep_pass_VVV_cutbased_veto_noiso_bak1_isLoaded = false;
  lep_pass_VVV_cutbased_veto_noiso_noip_bak1_isLoaded = false;
  lep_pass_VVV_cutbased_fo_isLoaded = false;
  lep_pass_VVV_cutbased_fo_noiso_isLoaded = false;
  lep_pass_VVV_cutbased_tight_noiso_isLoaded = false;
  lep_pass_VVV_cutbased_tight_isLoaded = false;
  lep_pass_VVV_MVAbased_tight_noiso_isLoaded = false;
  lep_pass_VVV_MVAbased_tight_isLoaded = false;
  lep_pass_VVV_baseline_isLoaded = false;
  nlep_VVV_cutbased_veto_isLoaded = false;
  nlep_VVV_cutbased_veto_noiso_isLoaded = false;
  nlep_VVV_cutbased_veto_noiso_noip_isLoaded = false;
  nlep_VVV_cutbased_veto_bak1_isLoaded = false;
  nlep_VVV_cutbased_veto_noiso_bak1_isLoaded = false;
  nlep_VVV_cutbased_veto_noiso_noip_bak1_isLoaded = false;
  nlep_VVV_cutbased_fo_isLoaded = false;
  nlep_VVV_cutbased_fo_noiso_isLoaded = false;
  nlep_VVV_cutbased_tight_noiso_isLoaded = false;
  nlep_VVV_cutbased_tight_isLoaded = false;
  nlep_VVV_MVAbased_tight_noiso_isLoaded = false;
  nlep_VVV_MVAbased_tight_isLoaded = false;
  nlep_VVV_baseline_isLoaded = false;
  lep_pdgId_isLoaded = false;
  lep_mc_Id_isLoaded = false;
  lep_dxy_isLoaded = false;
  lep_ip3d_isLoaded = false;
  lep_ip3derr_isLoaded = false;
  lep_etaSC_isLoaded = false;
  lep_r9_isLoaded = false;
  lep_dz_isLoaded = false;
  lep_tightId_isLoaded = false;
  lep_mcMatchId_isLoaded = false;
  lep_lostHits_isLoaded = false;
  lep_convVeto_isLoaded = false;
  lep_tightCharge_isLoaded = false;
  lep_MVA_isLoaded = false;
  lep_validfraction_isLoaded = false;
  lep_pterr_isLoaded = false;
  lep_sta_pterrOpt_isLoaded = false;
  lep_glb_pterrOpt_isLoaded = false;
  lep_x2ondof_isLoaded = false;
  lep_sta_x2ondof_isLoaded = false;
  lep_glb_x2ondof_isLoaded = false;
  nisoTrack_5gev_isLoaded = false;
  nisoTrack_mt2_isLoaded = false;
  nisoTrack_mt2_cleaned_VVV_cutbased_veto_isLoaded = false;
  nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_isLoaded = false;
  nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_isLoaded = false;
  nisoTrack_mt2_cleaned_VVV_cutbased_veto_bak1_isLoaded = false;
  nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_bak1_isLoaded = false;
  nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_bak1_isLoaded = false;
  nisoTrack_mt2_cleaned_VVV_cutbased_fo_isLoaded = false;
  nisoTrack_mt2_cleaned_VVV_cutbased_fo_noiso_isLoaded = false;
  nisoTrack_mt2_cleaned_VVV_cutbased_tight_noiso_isLoaded = false;
  nisoTrack_mt2_cleaned_VVV_cutbased_tight_isLoaded = false;
  nisoTrack_mt2_cleaned_VVV_MVAbased_tight_noiso_isLoaded = false;
  nisoTrack_mt2_cleaned_VVV_MVAbased_tight_isLoaded = false;
  nisoTrack_mt2_cleaned_VVV_baseline_isLoaded = false;
  nisoTrack_stop_isLoaded = false;
  nisoTrack_PFLep5_woverlaps_isLoaded = false;
  nisoTrack_PFHad10_woverlaps_isLoaded = false;
  ngamma_isLoaded = false;
  gamma_p4_isLoaded = false;
  gamma_pt_isLoaded = false;
  gamma_eta_isLoaded = false;
  gamma_phi_isLoaded = false;
  gamma_mass_isLoaded = false;
  gamma_mcMatchId_isLoaded = false;
  gamma_genPt_isLoaded = false;
  gamma_genIso_isLoaded = false;
  gamma_chHadIso_isLoaded = false;
  gamma_neuHadIso_isLoaded = false;
  gamma_phIso_isLoaded = false;
  gamma_sigmaIetaIeta_isLoaded = false;
  gamma_r9_isLoaded = false;
  gamma_hOverE_isLoaded = false;
  gamma_hOverE_online_isLoaded = false;
  gamma_idCutBased_isLoaded = false;
  gamma_ecpfclusiso_isLoaded = false;
  gamma_hcpfclusiso_isLoaded = false;
  gamma_hollowtkiso03_isLoaded = false;
  gamma_genIsPromptFinalState_isLoaded = false;
  gamma_drMinParton_isLoaded = false;
  ngenPart_isLoaded = false;
  genPart_p4_isLoaded = false;
  genPart_pt_isLoaded = false;
  genPart_eta_isLoaded = false;
  genPart_phi_isLoaded = false;
  genPart_mass_isLoaded = false;
  genPart_pdgId_isLoaded = false;
  genPart_status_isLoaded = false;
  genPart_charge_isLoaded = false;
  genPart_motherId_isLoaded = false;
  genPart_grandmaId_isLoaded = false;
  genPart_isp6status3_isLoaded = false;
  ngenLep_isLoaded = false;
  genLep_pt_isLoaded = false;
  genLep_eta_isLoaded = false;
  genLep_phi_isLoaded = false;
  genLep_mass_isLoaded = false;
  genLep_pdgId_isLoaded = false;
  genLep_status_isLoaded = false;
  genLep_charge_isLoaded = false;
  genLep_sourceId_isLoaded = false;
  genLep_isp6status3_isLoaded = false;
  ngenTau_isLoaded = false;
  genTau_pt_isLoaded = false;
  genTau_eta_isLoaded = false;
  genTau_phi_isLoaded = false;
  genTau_mass_isLoaded = false;
  genTau_pdgId_isLoaded = false;
  genTau_status_isLoaded = false;
  genTau_charge_isLoaded = false;
  genTau_sourceId_isLoaded = false;
  ngenLepFromTau_isLoaded = false;
  genLepFromTau_pt_isLoaded = false;
  genLepFromTau_eta_isLoaded = false;
  genLepFromTau_phi_isLoaded = false;
  genLepFromTau_mass_isLoaded = false;
  genLepFromTau_pdgId_isLoaded = false;
  genLepFromTau_status_isLoaded = false;
  genLepFromTau_charge_isLoaded = false;
  genLepFromTau_sourceId_isLoaded = false;
  njets_isLoaded = false;
  jets_p4_isLoaded = false;
  removed_jets_p4_isLoaded = false;
  removed_jets_csv_isLoaded = false;
  jets_medb_p4_isLoaded = false;
  njets_up_isLoaded = false;
  jets_up_p4_isLoaded = false;
  jets_medb_up_p4_isLoaded = false;
  jets_csv_isLoaded = false;
  jets_up_csv_isLoaded = false;
  njets_dn_isLoaded = false;
  jets_dn_p4_isLoaded = false;
  jets_medb_dn_p4_isLoaded = false;
  jets_dn_csv_isLoaded = false;
  jets_muf_isLoaded = false;
  jets_mcFlavour_isLoaded = false;
  jets_mcHadronFlav_isLoaded = false;
  nhighPtPFcands_isLoaded = false;
  highPtPFcands_p4_isLoaded = false;
  highPtPFcands_dz_isLoaded = false;
  highPtPFcands_pdgId_isLoaded = false;
  ht_isLoaded = false;
  ht_up_isLoaded = false;
  ht_dn_isLoaded = false;
  metsig_unofficial_isLoaded = false;
  metsig_unofficial_up_isLoaded = false;
  metsig_unofficial_dn_isLoaded = false;
  mt_lep1_isLoaded = false;
  mt2_isLoaded = false;
  mt2_up_isLoaded = false;
  mt2_dn_isLoaded = false;
  mt2j_isLoaded = false;
  mt2b_isLoaded = false;
  mt2b_up_isLoaded = false;
  mt2b_dn_isLoaded = false;
  mt2_genmet_isLoaded = false;
  mt2b_genmet_isLoaded = false;
  mjj_mindphi_isLoaded = false;
  mjj_isLoaded = false;
  mbb_csv_isLoaded = false;
  mbb_bpt_isLoaded = false;
  dphi_jj_isLoaded = false;
  dphi_ll_isLoaded = false;
  sum_mlb_isLoaded = false;
  deta_jj_isLoaded = false;
  dR_jj_isLoaded = false;
  dphi_metj1_isLoaded = false;
  dphi_metj2_isLoaded = false;
  dphi_genmetj1_isLoaded = false;
  dphi_genmetj2_isLoaded = false;
  mjj_mindphi_up_isLoaded = false;
  mjj_up_isLoaded = false;
  mbb_csv_up_isLoaded = false;
  mbb_bpt_up_isLoaded = false;
  dphi_jj_up_isLoaded = false;
  dphi_ll_up_isLoaded = false;
  sum_mlb_up_isLoaded = false;
  deta_jj_up_isLoaded = false;
  dR_jj_up_isLoaded = false;
  dphi_metj1_up_isLoaded = false;
  dphi_metj2_up_isLoaded = false;
  mjj_mindphi_dn_isLoaded = false;
  mjj_dn_isLoaded = false;
  mbb_csv_dn_isLoaded = false;
  mbb_bpt_dn_isLoaded = false;
  dphi_jj_dn_isLoaded = false;
  dphi_ll_dn_isLoaded = false;
  sum_mlb_dn_isLoaded = false;
  deta_jj_dn_isLoaded = false;
  dR_jj_dn_isLoaded = false;
  dphi_metj1_dn_isLoaded = false;
  dphi_metj2_dn_isLoaded = false;
  weight_btagsf_isLoaded = false;
  weight_btagsf_heavy_UP_isLoaded = false;
  weight_btagsf_light_UP_isLoaded = false;
  weight_btagsf_heavy_DN_isLoaded = false;
  weight_btagsf_light_DN_isLoaded = false;
  chpfcands_0013_pt_isLoaded = false;
  chpfcands_1316_pt_isLoaded = false;
  chpfcands_1624_pt_isLoaded = false;
  chpfcands_2430_pt_isLoaded = false;
  chpfcands_30in_pt_isLoaded = false;
  phpfcands_0013_pt_isLoaded = false;
  phpfcands_1316_pt_isLoaded = false;
  phpfcands_1624_pt_isLoaded = false;
  phpfcands_2430_pt_isLoaded = false;
  phpfcands_30in_pt_isLoaded = false;
  nupfcands_0013_pt_isLoaded = false;
  nupfcands_1316_pt_isLoaded = false;
  nupfcands_1624_pt_isLoaded = false;
  nupfcands_2430_pt_isLoaded = false;
  nupfcands_30in_pt_isLoaded = false;
  chpfcands_0013_sumet_isLoaded = false;
  chpfcands_1316_sumet_isLoaded = false;
  chpfcands_1624_sumet_isLoaded = false;
  chpfcands_2430_sumet_isLoaded = false;
  chpfcands_30in_sumet_isLoaded = false;
  phpfcands_0013_sumet_isLoaded = false;
  phpfcands_1316_sumet_isLoaded = false;
  phpfcands_1624_sumet_isLoaded = false;
  phpfcands_2430_sumet_isLoaded = false;
  phpfcands_30in_sumet_isLoaded = false;
  nupfcands_0013_sumet_isLoaded = false;
  nupfcands_1316_sumet_isLoaded = false;
  nupfcands_1624_sumet_isLoaded = false;
  nupfcands_2430_sumet_isLoaded = false;
  nupfcands_30in_sumet_isLoaded = false;
  chpfcands_0013_phi_isLoaded = false;
  chpfcands_1316_phi_isLoaded = false;
  chpfcands_1624_phi_isLoaded = false;
  chpfcands_2430_phi_isLoaded = false;
  chpfcands_30in_phi_isLoaded = false;
  phpfcands_0013_phi_isLoaded = false;
  phpfcands_1316_phi_isLoaded = false;
  phpfcands_1624_phi_isLoaded = false;
  phpfcands_2430_phi_isLoaded = false;
  phpfcands_30in_phi_isLoaded = false;
  nupfcands_0013_phi_isLoaded = false;
  nupfcands_1316_phi_isLoaded = false;
  nupfcands_1624_phi_isLoaded = false;
  nupfcands_2430_phi_isLoaded = false;
  nupfcands_30in_phi_isLoaded = false;
  met_T1CHS_pt_isLoaded = false;
  met_T1CHS_phi_isLoaded = false;
  met_T1CHS_fromCORE_pt_isLoaded = false;
  met_T1CHS_fromCORE_phi_isLoaded = false;
  met_T1CHS_miniAOD_CORE_pt_isLoaded = false;
  met_T1CHS_miniAOD_CORE_phi_isLoaded = false;
  met_T1CHS_miniAOD_CORE_up_pt_isLoaded = false;
  met_T1CHS_miniAOD_CORE_up_phi_isLoaded = false;
  met_T1CHS_miniAOD_CORE_dn_pt_isLoaded = false;
  met_T1CHS_miniAOD_CORE_dn_phi_isLoaded = false;
  decayedphoton_lep1_p4_isLoaded = false;
  decayedphoton_lep2_p4_isLoaded = false;
  decayedphoton_bosn_p4_isLoaded = false;
  decayedphoton_mt2_isLoaded = false;
  hyp_type_isLoaded = false;
  evt_type_isLoaded = false;
  mass_gluino_isLoaded = false;
  mass_LSP_isLoaded = false;
  mass_chi_isLoaded = false;
  isrboost_isLoaded = false;
  isr_njets_isLoaded = false;
  isr_weight_isLoaded = false;
  isr_unc_isLoaded = false;
  weightsf_lepid_isLoaded = false;
  weightsf_lepiso_isLoaded = false;
  weightsf_lepip_isLoaded = false;
  weightsf_lepreco_isLoaded = false;
  weightsf_lepconv_isLoaded = false;
  weightsf_lepid_FS_isLoaded = false;
  weightsf_lepiso_FS_isLoaded = false;
  weightsf_lepip_FS_isLoaded = false;
  weight_rn_r1_n1_isLoaded = false;
  weight_rn_r1_n2_isLoaded = false;
  weight_rn_r1_np05_isLoaded = false;
  weight_rn_r2_n1_isLoaded = false;
  weight_rn_r2_n2_isLoaded = false;
  weight_rn_r2_n0p5_isLoaded = false;
  weight_rn_r0p5_n1_isLoaded = false;
  weight_rn_r0p5_n2_isLoaded = false;
  weight_rn_r0p5_n0p5_isLoaded = false;
  weight_pdf_up_isLoaded = false;
  weight_pdf_down_isLoaded = false;
  weight_alphas_down_isLoaded = false;
  weight_alphas_up_isLoaded = false;
  LooseElec_isLoaded = false;
  LooseLepton_isLoaded = false;
  LooseMuon_isLoaded = false;
  SignalElec_isLoaded = false;
  SignalLepton_isLoaded = false;
  SignalMuon_isLoaded = false;
  Tight3lElec_isLoaded = false;
  Tight3lLepton_isLoaded = false;
  TightElec_isLoaded = false;
  TightLepton_isLoaded = false;
  TightMuon_isLoaded = false;
  VetoElec_isLoaded = false;
  VetoLepton_isLoaded = false;
  VetoMuon_isLoaded = false;
  Good3LJet_isLoaded = false;
  GoodSSJet_isLoaded = false;
  GoodSSWJet_isLoaded = false;
  LooseBJet_isLoaded = false;
}

void WWWTree::LoadAllBranches() {
  // load all branches
  if (run_branch != 0) run();
  if (lumi_branch != 0) lumi();
  if (evt_branch != 0) evt();
  if (isData_branch != 0) isData();
  if (evt_passgoodrunlist_branch != 0) evt_passgoodrunlist();
  if (evt_scale1fb_branch != 0) evt_scale1fb();
  if (evt_xsec_branch != 0) evt_xsec();
  if (evt_kfactor_branch != 0) evt_kfactor();
  if (evt_filter_branch != 0) evt_filter();
  if (evt_nEvts_branch != 0) evt_nEvts();
  if (evt_dataset_branch != 0) evt_dataset();
  if (puWeight_branch != 0) puWeight();
  if (nVert_branch != 0) nVert();
  if (firstgoodvertex_branch != 0) firstgoodvertex();
  if (nTrueInt_branch != 0) nTrueInt();
  if (rho_branch != 0) rho();
  if (rho25_branch != 0) rho25();
  if (gen_ht_branch != 0) gen_ht();
  if (nBJetTight_branch != 0) nBJetTight();
  if (nBJetMedium_branch != 0) nBJetMedium();
  if (nBJetLoose_branch != 0) nBJetLoose();
  if (nBJetTight_up_branch != 0) nBJetTight_up();
  if (nBJetMedium_up_branch != 0) nBJetMedium_up();
  if (nBJetLoose_up_branch != 0) nBJetLoose_up();
  if (nBJetTight_dn_branch != 0) nBJetTight_dn();
  if (nBJetMedium_dn_branch != 0) nBJetMedium_dn();
  if (nBJetLoose_dn_branch != 0) nBJetLoose_dn();
  if (nJet200MuFrac50DphiMet_branch != 0) nJet200MuFrac50DphiMet();
  if (nMuons10_branch != 0) nMuons10();
  if (nBadMuons20_branch != 0) nBadMuons20();
  if (nElectrons10_branch != 0) nElectrons10();
  if (nGammas20_branch != 0) nGammas20();
  if (nTaus20_branch != 0) nTaus20();
  if (met_pt_branch != 0) met_pt();
  if (met_phi_branch != 0) met_phi();
  if (met_calo_pt_branch != 0) met_calo_pt();
  if (met_calo_phi_branch != 0) met_calo_phi();
  if (met_miniaod_pt_branch != 0) met_miniaod_pt();
  if (met_miniaod_phi_branch != 0) met_miniaod_phi();
  if (met_muegclean_pt_branch != 0) met_muegclean_pt();
  if (met_muegclean_phi_branch != 0) met_muegclean_phi();
  if (met_rawPt_branch != 0) met_rawPt();
  if (met_rawPhi_branch != 0) met_rawPhi();
  if (met_genPt_branch != 0) met_genPt();
  if (met_genPhi_branch != 0) met_genPhi();
  if (sumet_raw_branch != 0) sumet_raw();
  if (Flag_ecalLaserCorrFilter_branch != 0) Flag_ecalLaserCorrFilter();
  if (Flag_hcalLaserEventFilter_branch != 0) Flag_hcalLaserEventFilter();
  if (Flag_trackingFailureFilter_branch != 0) Flag_trackingFailureFilter();
  if (Flag_CSCTightHaloFilter_branch != 0) Flag_CSCTightHaloFilter();
  if (Flag_HBHENoiseFilter_branch != 0) Flag_HBHENoiseFilter();
  if (Flag_HBHEIsoNoiseFilter_branch != 0) Flag_HBHEIsoNoiseFilter();
  if (Flag_CSCTightHalo2015Filter_branch != 0) Flag_CSCTightHalo2015Filter();
  if (Flag_EcalDeadCellTriggerPrimitiveFilter_branch != 0) Flag_EcalDeadCellTriggerPrimitiveFilter();
  if (Flag_goodVertices_branch != 0) Flag_goodVertices();
  if (Flag_eeBadScFilter_branch != 0) Flag_eeBadScFilter();
  if (Flag_globalTightHalo2016_branch != 0) Flag_globalTightHalo2016();
  if (Flag_badMuonFilter_branch != 0) Flag_badMuonFilter();
  if (Flag_badChargedCandidateFilter_branch != 0) Flag_badChargedCandidateFilter();
  if (Flag_badMuonFilterv2_branch != 0) Flag_badMuonFilterv2();
  if (Flag_badChargedCandidateFilterv2_branch != 0) Flag_badChargedCandidateFilterv2();
  if (Flag_badMuons_branch != 0) Flag_badMuons();
  if (Flag_duplicateMuons_branch != 0) Flag_duplicateMuons();
  if (Flag_noBadMuons_branch != 0) Flag_noBadMuons();
  if (HLT_singleEl_branch != 0) HLT_singleEl();
  if (HLT_singleMu_branch != 0) HLT_singleMu();
  if (HLT_singleMu_noiso_branch != 0) HLT_singleMu_noiso();
  if (HLT_DoubleEl_noiso_branch != 0) HLT_DoubleEl_noiso();
  if (HLT_DoubleEl_branch != 0) HLT_DoubleEl();
  if (HLT_DoubleEl_DZ_branch != 0) HLT_DoubleEl_DZ();
  if (HLT_DoubleEl_DZ_2_branch != 0) HLT_DoubleEl_DZ_2();
  if (HLT_MuEG_branch != 0) HLT_MuEG();
  if (HLT_MuEG_2_branch != 0) HLT_MuEG_2();
  if (HLT_MuEG_noiso_branch != 0) HLT_MuEG_noiso();
  if (HLT_MuEG_noiso_2_branch != 0) HLT_MuEG_noiso_2();
  if (HLT_Mu8_EG17_branch != 0) HLT_Mu8_EG17();
  if (HLT_Mu8_EG23_branch != 0) HLT_Mu8_EG23();
  if (HLT_Mu8_EG23_DZ_branch != 0) HLT_Mu8_EG23_DZ();
  if (HLT_Mu12_EG23_DZ_branch != 0) HLT_Mu12_EG23_DZ();
  if (HLT_Mu17_EG12_branch != 0) HLT_Mu17_EG12();
  if (HLT_Mu23_EG8_branch != 0) HLT_Mu23_EG8();
  if (HLT_Mu23_EG8_DZ_branch != 0) HLT_Mu23_EG8_DZ();
  if (HLT_Mu23_EG12_branch != 0) HLT_Mu23_EG12();
  if (HLT_Mu23_EG12_DZ_branch != 0) HLT_Mu23_EG12_DZ();
  if (HLT_DoubleMu_noiso_branch != 0) HLT_DoubleMu_noiso();
  if (HLT_DoubleMu_noiso_27_8_branch != 0) HLT_DoubleMu_noiso_27_8();
  if (HLT_DoubleMu_noiso_30_11_branch != 0) HLT_DoubleMu_noiso_30_11();
  if (HLT_DoubleMu_noiso_40_11_branch != 0) HLT_DoubleMu_noiso_40_11();
  if (HLT_DoubleMu_branch != 0) HLT_DoubleMu();
  if (HLT_DoubleMu_tk_branch != 0) HLT_DoubleMu_tk();
  if (HLT_DoubleMu_dbltk_branch != 0) HLT_DoubleMu_dbltk();
  if (HLT_DoubleMu_nonDZ_branch != 0) HLT_DoubleMu_nonDZ();
  if (HLT_DoubleMu_tk_nonDZ_branch != 0) HLT_DoubleMu_tk_nonDZ();
  if (HLT_Photon22_R9Id90_HE10_IsoM_branch != 0) HLT_Photon22_R9Id90_HE10_IsoM();
  if (HLT_Photon30_R9Id90_HE10_IsoM_branch != 0) HLT_Photon30_R9Id90_HE10_IsoM();
  if (HLT_Photon36_R9Id90_HE10_IsoM_branch != 0) HLT_Photon36_R9Id90_HE10_IsoM();
  if (HLT_Photon50_R9Id90_HE10_IsoM_branch != 0) HLT_Photon50_R9Id90_HE10_IsoM();
  if (HLT_Photon75_R9Id90_HE10_IsoM_branch != 0) HLT_Photon75_R9Id90_HE10_IsoM();
  if (HLT_Photon90_R9Id90_HE10_IsoM_branch != 0) HLT_Photon90_R9Id90_HE10_IsoM();
  if (HLT_Photon120_R9Id90_HE10_IsoM_branch != 0) HLT_Photon120_R9Id90_HE10_IsoM();
  if (HLT_Photon165_R9Id90_HE10_IsoM_branch != 0) HLT_Photon165_R9Id90_HE10_IsoM();
  if (HLT_Photon165_HE10_branch != 0) HLT_Photon165_HE10();
  if (HLT_CaloJet500_NoJetID_branch != 0) HLT_CaloJet500_NoJetID();
  if (HLT_ECALHT800_NoJetID_branch != 0) HLT_ECALHT800_NoJetID();
  if (HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton_branch != 0) HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton();
  if (HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton_branch != 0) HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton();
  if (HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton_branch != 0) HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton();
  if (HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton_branch != 0) HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton();
  if (HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton_branch != 0) HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton();
  if (HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton_branch != 0) HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton();
  if (HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton_branch != 0) HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton();
  if (HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton_branch != 0) HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton();
  if (HLT_Photon165_HE10_matchedtophoton_branch != 0) HLT_Photon165_HE10_matchedtophoton();
  if (dilmass_branch != 0) dilmass();
  if (dilpt_branch != 0) dilpt();
  if (dRll_branch != 0) dRll();
  if (matched_neutralemf_branch != 0) matched_neutralemf();
  if (matched_emf_branch != 0) matched_emf();
  if (elveto_branch != 0) elveto();
  if (nlep_branch != 0) nlep();
  if (nveto_leptons_branch != 0) nveto_leptons();
  if (nVetoEl_relIso03EAless01_branch != 0) nVetoEl_relIso03EAless01();
  if (nVetoEl_relIso03EAless02_branch != 0) nVetoEl_relIso03EAless02();
  if (nVetoEl_relIso03EAless03_branch != 0) nVetoEl_relIso03EAless03();
  if (nVetoEl_relIso03EAless04_branch != 0) nVetoEl_relIso03EAless04();
  if (nVetoMu_relIso03EAless01_branch != 0) nVetoMu_relIso03EAless01();
  if (nVetoMu_relIso03EAless02_branch != 0) nVetoMu_relIso03EAless02();
  if (nVetoMu_relIso03EAless03_branch != 0) nVetoMu_relIso03EAless03();
  if (nVetoMu_relIso03EAless04_branch != 0) nVetoMu_relIso03EAless04();
  if (lep_p4_branch != 0) lep_p4();
  if (lep_pt_branch != 0) lep_pt();
  if (lep_eta_branch != 0) lep_eta();
  if (lep_phi_branch != 0) lep_phi();
  if (lep_mass_branch != 0) lep_mass();
  if (lep_charge_branch != 0) lep_charge();
  if (lep_3ch_agree_branch != 0) lep_3ch_agree();
  if (lep_isFromW_branch != 0) lep_isFromW();
  if (lep_isFromZ_branch != 0) lep_isFromZ();
  if (lep_isFromB_branch != 0) lep_isFromB();
  if (lep_isFromC_branch != 0) lep_isFromC();
  if (lep_isFromL_branch != 0) lep_isFromL();
  if (lep_isFromLF_branch != 0) lep_isFromLF();
  if (lep_closest_jet_p4_branch != 0) lep_closest_jet_p4();
  if (lep_ptRatio_branch != 0) lep_ptRatio();
  if (lep_motherIdSS_branch != 0) lep_motherIdSS();
  if (lep_genPart_index_branch != 0) lep_genPart_index();
  if (lep_coneCorrPt_branch != 0) lep_coneCorrPt();
  if (lep_ptRel_branch != 0) lep_ptRel();
  if (lep_relIso03_branch != 0) lep_relIso03();
  if (lep_relIso03DB_branch != 0) lep_relIso03DB();
  if (lep_relIso03EA_branch != 0) lep_relIso03EA();
  if (lep_relIso03EAv2_branch != 0) lep_relIso03EAv2();
  if (lep_relIso04DB_branch != 0) lep_relIso04DB();
  if (lep_relIso04EA_branch != 0) lep_relIso04EA();
  if (lep_relIso04EAv2_branch != 0) lep_relIso04EAv2();
  if (lep_miniRelIsoCMS3_EA_branch != 0) lep_miniRelIsoCMS3_EA();
  if (lep_miniRelIsoCMS3_EAv2_branch != 0) lep_miniRelIsoCMS3_EAv2();
  if (lep_miniRelIsoCMS3_DB_branch != 0) lep_miniRelIsoCMS3_DB();
  if (lep_isTriggerSafenoIso_v1_branch != 0) lep_isTriggerSafenoIso_v1();
  if (lep_isTriggerSafenoIso_v2_branch != 0) lep_isTriggerSafenoIso_v2();
  if (lep_isTriggerSafe_v1_branch != 0) lep_isTriggerSafe_v1();
  if (lep_isTriggerSafe_v2_branch != 0) lep_isTriggerSafe_v2();
  if (lep_pass_VVV_cutbased_veto_branch != 0) lep_pass_VVV_cutbased_veto();
  if (lep_pass_VVV_cutbased_veto_noiso_branch != 0) lep_pass_VVV_cutbased_veto_noiso();
  if (lep_pass_VVV_cutbased_veto_noiso_noip_branch != 0) lep_pass_VVV_cutbased_veto_noiso_noip();
  if (lep_pass_VVV_cutbased_veto_bak1_branch != 0) lep_pass_VVV_cutbased_veto_bak1();
  if (lep_pass_VVV_cutbased_veto_noiso_bak1_branch != 0) lep_pass_VVV_cutbased_veto_noiso_bak1();
  if (lep_pass_VVV_cutbased_veto_noiso_noip_bak1_branch != 0) lep_pass_VVV_cutbased_veto_noiso_noip_bak1();
  if (lep_pass_VVV_cutbased_fo_branch != 0) lep_pass_VVV_cutbased_fo();
  if (lep_pass_VVV_cutbased_fo_noiso_branch != 0) lep_pass_VVV_cutbased_fo_noiso();
  if (lep_pass_VVV_cutbased_tight_noiso_branch != 0) lep_pass_VVV_cutbased_tight_noiso();
  if (lep_pass_VVV_cutbased_tight_branch != 0) lep_pass_VVV_cutbased_tight();
  if (lep_pass_VVV_MVAbased_tight_noiso_branch != 0) lep_pass_VVV_MVAbased_tight_noiso();
  if (lep_pass_VVV_MVAbased_tight_branch != 0) lep_pass_VVV_MVAbased_tight();
  if (lep_pass_VVV_baseline_branch != 0) lep_pass_VVV_baseline();
  if (nlep_VVV_cutbased_veto_branch != 0) nlep_VVV_cutbased_veto();
  if (nlep_VVV_cutbased_veto_noiso_branch != 0) nlep_VVV_cutbased_veto_noiso();
  if (nlep_VVV_cutbased_veto_noiso_noip_branch != 0) nlep_VVV_cutbased_veto_noiso_noip();
  if (nlep_VVV_cutbased_veto_bak1_branch != 0) nlep_VVV_cutbased_veto_bak1();
  if (nlep_VVV_cutbased_veto_noiso_bak1_branch != 0) nlep_VVV_cutbased_veto_noiso_bak1();
  if (nlep_VVV_cutbased_veto_noiso_noip_bak1_branch != 0) nlep_VVV_cutbased_veto_noiso_noip_bak1();
  if (nlep_VVV_cutbased_fo_branch != 0) nlep_VVV_cutbased_fo();
  if (nlep_VVV_cutbased_fo_noiso_branch != 0) nlep_VVV_cutbased_fo_noiso();
  if (nlep_VVV_cutbased_tight_noiso_branch != 0) nlep_VVV_cutbased_tight_noiso();
  if (nlep_VVV_cutbased_tight_branch != 0) nlep_VVV_cutbased_tight();
  if (nlep_VVV_MVAbased_tight_noiso_branch != 0) nlep_VVV_MVAbased_tight_noiso();
  if (nlep_VVV_MVAbased_tight_branch != 0) nlep_VVV_MVAbased_tight();
  if (nlep_VVV_baseline_branch != 0) nlep_VVV_baseline();
  if (lep_pdgId_branch != 0) lep_pdgId();
  if (lep_mc_Id_branch != 0) lep_mc_Id();
  if (lep_dxy_branch != 0) lep_dxy();
  if (lep_ip3d_branch != 0) lep_ip3d();
  if (lep_ip3derr_branch != 0) lep_ip3derr();
  if (lep_etaSC_branch != 0) lep_etaSC();
  if (lep_r9_branch != 0) lep_r9();
  if (lep_dz_branch != 0) lep_dz();
  if (lep_tightId_branch != 0) lep_tightId();
  if (lep_mcMatchId_branch != 0) lep_mcMatchId();
  if (lep_lostHits_branch != 0) lep_lostHits();
  if (lep_convVeto_branch != 0) lep_convVeto();
  if (lep_tightCharge_branch != 0) lep_tightCharge();
  if (lep_MVA_branch != 0) lep_MVA();
  if (lep_validfraction_branch != 0) lep_validfraction();
  if (lep_pterr_branch != 0) lep_pterr();
  if (lep_sta_pterrOpt_branch != 0) lep_sta_pterrOpt();
  if (lep_glb_pterrOpt_branch != 0) lep_glb_pterrOpt();
  if (lep_x2ondof_branch != 0) lep_x2ondof();
  if (lep_sta_x2ondof_branch != 0) lep_sta_x2ondof();
  if (lep_glb_x2ondof_branch != 0) lep_glb_x2ondof();
  if (nisoTrack_5gev_branch != 0) nisoTrack_5gev();
  if (nisoTrack_mt2_branch != 0) nisoTrack_mt2();
  if (nisoTrack_mt2_cleaned_VVV_cutbased_veto_branch != 0) nisoTrack_mt2_cleaned_VVV_cutbased_veto();
  if (nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_branch != 0) nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso();
  if (nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_branch != 0) nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip();
  if (nisoTrack_mt2_cleaned_VVV_cutbased_veto_bak1_branch != 0) nisoTrack_mt2_cleaned_VVV_cutbased_veto_bak1();
  if (nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_bak1_branch != 0) nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_bak1();
  if (nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_bak1_branch != 0) nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_bak1();
  if (nisoTrack_mt2_cleaned_VVV_cutbased_fo_branch != 0) nisoTrack_mt2_cleaned_VVV_cutbased_fo();
  if (nisoTrack_mt2_cleaned_VVV_cutbased_fo_noiso_branch != 0) nisoTrack_mt2_cleaned_VVV_cutbased_fo_noiso();
  if (nisoTrack_mt2_cleaned_VVV_cutbased_tight_noiso_branch != 0) nisoTrack_mt2_cleaned_VVV_cutbased_tight_noiso();
  if (nisoTrack_mt2_cleaned_VVV_cutbased_tight_branch != 0) nisoTrack_mt2_cleaned_VVV_cutbased_tight();
  if (nisoTrack_mt2_cleaned_VVV_MVAbased_tight_noiso_branch != 0) nisoTrack_mt2_cleaned_VVV_MVAbased_tight_noiso();
  if (nisoTrack_mt2_cleaned_VVV_MVAbased_tight_branch != 0) nisoTrack_mt2_cleaned_VVV_MVAbased_tight();
  if (nisoTrack_mt2_cleaned_VVV_baseline_branch != 0) nisoTrack_mt2_cleaned_VVV_baseline();
  if (nisoTrack_stop_branch != 0) nisoTrack_stop();
  if (nisoTrack_PFLep5_woverlaps_branch != 0) nisoTrack_PFLep5_woverlaps();
  if (nisoTrack_PFHad10_woverlaps_branch != 0) nisoTrack_PFHad10_woverlaps();
  if (ngamma_branch != 0) ngamma();
  if (gamma_p4_branch != 0) gamma_p4();
  if (gamma_pt_branch != 0) gamma_pt();
  if (gamma_eta_branch != 0) gamma_eta();
  if (gamma_phi_branch != 0) gamma_phi();
  if (gamma_mass_branch != 0) gamma_mass();
  if (gamma_mcMatchId_branch != 0) gamma_mcMatchId();
  if (gamma_genPt_branch != 0) gamma_genPt();
  if (gamma_genIso_branch != 0) gamma_genIso();
  if (gamma_chHadIso_branch != 0) gamma_chHadIso();
  if (gamma_neuHadIso_branch != 0) gamma_neuHadIso();
  if (gamma_phIso_branch != 0) gamma_phIso();
  if (gamma_sigmaIetaIeta_branch != 0) gamma_sigmaIetaIeta();
  if (gamma_r9_branch != 0) gamma_r9();
  if (gamma_hOverE_branch != 0) gamma_hOverE();
  if (gamma_hOverE_online_branch != 0) gamma_hOverE_online();
  if (gamma_idCutBased_branch != 0) gamma_idCutBased();
  if (gamma_ecpfclusiso_branch != 0) gamma_ecpfclusiso();
  if (gamma_hcpfclusiso_branch != 0) gamma_hcpfclusiso();
  if (gamma_hollowtkiso03_branch != 0) gamma_hollowtkiso03();
  if (gamma_genIsPromptFinalState_branch != 0) gamma_genIsPromptFinalState();
  if (gamma_drMinParton_branch != 0) gamma_drMinParton();
  if (ngenPart_branch != 0) ngenPart();
  if (genPart_p4_branch != 0) genPart_p4();
  if (genPart_pt_branch != 0) genPart_pt();
  if (genPart_eta_branch != 0) genPart_eta();
  if (genPart_phi_branch != 0) genPart_phi();
  if (genPart_mass_branch != 0) genPart_mass();
  if (genPart_pdgId_branch != 0) genPart_pdgId();
  if (genPart_status_branch != 0) genPart_status();
  if (genPart_charge_branch != 0) genPart_charge();
  if (genPart_motherId_branch != 0) genPart_motherId();
  if (genPart_grandmaId_branch != 0) genPart_grandmaId();
  if (genPart_isp6status3_branch != 0) genPart_isp6status3();
  if (ngenLep_branch != 0) ngenLep();
  if (genLep_pt_branch != 0) genLep_pt();
  if (genLep_eta_branch != 0) genLep_eta();
  if (genLep_phi_branch != 0) genLep_phi();
  if (genLep_mass_branch != 0) genLep_mass();
  if (genLep_pdgId_branch != 0) genLep_pdgId();
  if (genLep_status_branch != 0) genLep_status();
  if (genLep_charge_branch != 0) genLep_charge();
  if (genLep_sourceId_branch != 0) genLep_sourceId();
  if (genLep_isp6status3_branch != 0) genLep_isp6status3();
  if (ngenTau_branch != 0) ngenTau();
  if (genTau_pt_branch != 0) genTau_pt();
  if (genTau_eta_branch != 0) genTau_eta();
  if (genTau_phi_branch != 0) genTau_phi();
  if (genTau_mass_branch != 0) genTau_mass();
  if (genTau_pdgId_branch != 0) genTau_pdgId();
  if (genTau_status_branch != 0) genTau_status();
  if (genTau_charge_branch != 0) genTau_charge();
  if (genTau_sourceId_branch != 0) genTau_sourceId();
  if (ngenLepFromTau_branch != 0) ngenLepFromTau();
  if (genLepFromTau_pt_branch != 0) genLepFromTau_pt();
  if (genLepFromTau_eta_branch != 0) genLepFromTau_eta();
  if (genLepFromTau_phi_branch != 0) genLepFromTau_phi();
  if (genLepFromTau_mass_branch != 0) genLepFromTau_mass();
  if (genLepFromTau_pdgId_branch != 0) genLepFromTau_pdgId();
  if (genLepFromTau_status_branch != 0) genLepFromTau_status();
  if (genLepFromTau_charge_branch != 0) genLepFromTau_charge();
  if (genLepFromTau_sourceId_branch != 0) genLepFromTau_sourceId();
  if (njets_branch != 0) njets();
  if (jets_p4_branch != 0) jets_p4();
  if (removed_jets_p4_branch != 0) removed_jets_p4();
  if (removed_jets_csv_branch != 0) removed_jets_csv();
  if (jets_medb_p4_branch != 0) jets_medb_p4();
  if (njets_up_branch != 0) njets_up();
  if (jets_up_p4_branch != 0) jets_up_p4();
  if (jets_medb_up_p4_branch != 0) jets_medb_up_p4();
  if (jets_csv_branch != 0) jets_csv();
  if (jets_up_csv_branch != 0) jets_up_csv();
  if (njets_dn_branch != 0) njets_dn();
  if (jets_dn_p4_branch != 0) jets_dn_p4();
  if (jets_medb_dn_p4_branch != 0) jets_medb_dn_p4();
  if (jets_dn_csv_branch != 0) jets_dn_csv();
  if (jets_muf_branch != 0) jets_muf();
  if (jets_mcFlavour_branch != 0) jets_mcFlavour();
  if (jets_mcHadronFlav_branch != 0) jets_mcHadronFlav();
  if (nhighPtPFcands_branch != 0) nhighPtPFcands();
  if (highPtPFcands_p4_branch != 0) highPtPFcands_p4();
  if (highPtPFcands_dz_branch != 0) highPtPFcands_dz();
  if (highPtPFcands_pdgId_branch != 0) highPtPFcands_pdgId();
  if (ht_branch != 0) ht();
  if (ht_up_branch != 0) ht_up();
  if (ht_dn_branch != 0) ht_dn();
  if (metsig_unofficial_branch != 0) metsig_unofficial();
  if (metsig_unofficial_up_branch != 0) metsig_unofficial_up();
  if (metsig_unofficial_dn_branch != 0) metsig_unofficial_dn();
  if (mt_lep1_branch != 0) mt_lep1();
  if (mt2_branch != 0) mt2();
  if (mt2_up_branch != 0) mt2_up();
  if (mt2_dn_branch != 0) mt2_dn();
  if (mt2j_branch != 0) mt2j();
  if (mt2b_branch != 0) mt2b();
  if (mt2b_up_branch != 0) mt2b_up();
  if (mt2b_dn_branch != 0) mt2b_dn();
  if (mt2_genmet_branch != 0) mt2_genmet();
  if (mt2b_genmet_branch != 0) mt2b_genmet();
  if (mjj_mindphi_branch != 0) mjj_mindphi();
  if (mjj_branch != 0) mjj();
  if (mbb_csv_branch != 0) mbb_csv();
  if (mbb_bpt_branch != 0) mbb_bpt();
  if (dphi_jj_branch != 0) dphi_jj();
  if (dphi_ll_branch != 0) dphi_ll();
  if (sum_mlb_branch != 0) sum_mlb();
  if (deta_jj_branch != 0) deta_jj();
  if (dR_jj_branch != 0) dR_jj();
  if (dphi_metj1_branch != 0) dphi_metj1();
  if (dphi_metj2_branch != 0) dphi_metj2();
  if (dphi_genmetj1_branch != 0) dphi_genmetj1();
  if (dphi_genmetj2_branch != 0) dphi_genmetj2();
  if (mjj_mindphi_up_branch != 0) mjj_mindphi_up();
  if (mjj_up_branch != 0) mjj_up();
  if (mbb_csv_up_branch != 0) mbb_csv_up();
  if (mbb_bpt_up_branch != 0) mbb_bpt_up();
  if (dphi_jj_up_branch != 0) dphi_jj_up();
  if (dphi_ll_up_branch != 0) dphi_ll_up();
  if (sum_mlb_up_branch != 0) sum_mlb_up();
  if (deta_jj_up_branch != 0) deta_jj_up();
  if (dR_jj_up_branch != 0) dR_jj_up();
  if (dphi_metj1_up_branch != 0) dphi_metj1_up();
  if (dphi_metj2_up_branch != 0) dphi_metj2_up();
  if (mjj_mindphi_dn_branch != 0) mjj_mindphi_dn();
  if (mjj_dn_branch != 0) mjj_dn();
  if (mbb_csv_dn_branch != 0) mbb_csv_dn();
  if (mbb_bpt_dn_branch != 0) mbb_bpt_dn();
  if (dphi_jj_dn_branch != 0) dphi_jj_dn();
  if (dphi_ll_dn_branch != 0) dphi_ll_dn();
  if (sum_mlb_dn_branch != 0) sum_mlb_dn();
  if (deta_jj_dn_branch != 0) deta_jj_dn();
  if (dR_jj_dn_branch != 0) dR_jj_dn();
  if (dphi_metj1_dn_branch != 0) dphi_metj1_dn();
  if (dphi_metj2_dn_branch != 0) dphi_metj2_dn();
  if (weight_btagsf_branch != 0) weight_btagsf();
  if (weight_btagsf_heavy_UP_branch != 0) weight_btagsf_heavy_UP();
  if (weight_btagsf_light_UP_branch != 0) weight_btagsf_light_UP();
  if (weight_btagsf_heavy_DN_branch != 0) weight_btagsf_heavy_DN();
  if (weight_btagsf_light_DN_branch != 0) weight_btagsf_light_DN();
  if (chpfcands_0013_pt_branch != 0) chpfcands_0013_pt();
  if (chpfcands_1316_pt_branch != 0) chpfcands_1316_pt();
  if (chpfcands_1624_pt_branch != 0) chpfcands_1624_pt();
  if (chpfcands_2430_pt_branch != 0) chpfcands_2430_pt();
  if (chpfcands_30in_pt_branch != 0) chpfcands_30in_pt();
  if (phpfcands_0013_pt_branch != 0) phpfcands_0013_pt();
  if (phpfcands_1316_pt_branch != 0) phpfcands_1316_pt();
  if (phpfcands_1624_pt_branch != 0) phpfcands_1624_pt();
  if (phpfcands_2430_pt_branch != 0) phpfcands_2430_pt();
  if (phpfcands_30in_pt_branch != 0) phpfcands_30in_pt();
  if (nupfcands_0013_pt_branch != 0) nupfcands_0013_pt();
  if (nupfcands_1316_pt_branch != 0) nupfcands_1316_pt();
  if (nupfcands_1624_pt_branch != 0) nupfcands_1624_pt();
  if (nupfcands_2430_pt_branch != 0) nupfcands_2430_pt();
  if (nupfcands_30in_pt_branch != 0) nupfcands_30in_pt();
  if (chpfcands_0013_sumet_branch != 0) chpfcands_0013_sumet();
  if (chpfcands_1316_sumet_branch != 0) chpfcands_1316_sumet();
  if (chpfcands_1624_sumet_branch != 0) chpfcands_1624_sumet();
  if (chpfcands_2430_sumet_branch != 0) chpfcands_2430_sumet();
  if (chpfcands_30in_sumet_branch != 0) chpfcands_30in_sumet();
  if (phpfcands_0013_sumet_branch != 0) phpfcands_0013_sumet();
  if (phpfcands_1316_sumet_branch != 0) phpfcands_1316_sumet();
  if (phpfcands_1624_sumet_branch != 0) phpfcands_1624_sumet();
  if (phpfcands_2430_sumet_branch != 0) phpfcands_2430_sumet();
  if (phpfcands_30in_sumet_branch != 0) phpfcands_30in_sumet();
  if (nupfcands_0013_sumet_branch != 0) nupfcands_0013_sumet();
  if (nupfcands_1316_sumet_branch != 0) nupfcands_1316_sumet();
  if (nupfcands_1624_sumet_branch != 0) nupfcands_1624_sumet();
  if (nupfcands_2430_sumet_branch != 0) nupfcands_2430_sumet();
  if (nupfcands_30in_sumet_branch != 0) nupfcands_30in_sumet();
  if (chpfcands_0013_phi_branch != 0) chpfcands_0013_phi();
  if (chpfcands_1316_phi_branch != 0) chpfcands_1316_phi();
  if (chpfcands_1624_phi_branch != 0) chpfcands_1624_phi();
  if (chpfcands_2430_phi_branch != 0) chpfcands_2430_phi();
  if (chpfcands_30in_phi_branch != 0) chpfcands_30in_phi();
  if (phpfcands_0013_phi_branch != 0) phpfcands_0013_phi();
  if (phpfcands_1316_phi_branch != 0) phpfcands_1316_phi();
  if (phpfcands_1624_phi_branch != 0) phpfcands_1624_phi();
  if (phpfcands_2430_phi_branch != 0) phpfcands_2430_phi();
  if (phpfcands_30in_phi_branch != 0) phpfcands_30in_phi();
  if (nupfcands_0013_phi_branch != 0) nupfcands_0013_phi();
  if (nupfcands_1316_phi_branch != 0) nupfcands_1316_phi();
  if (nupfcands_1624_phi_branch != 0) nupfcands_1624_phi();
  if (nupfcands_2430_phi_branch != 0) nupfcands_2430_phi();
  if (nupfcands_30in_phi_branch != 0) nupfcands_30in_phi();
  if (met_T1CHS_pt_branch != 0) met_T1CHS_pt();
  if (met_T1CHS_phi_branch != 0) met_T1CHS_phi();
  if (met_T1CHS_fromCORE_pt_branch != 0) met_T1CHS_fromCORE_pt();
  if (met_T1CHS_fromCORE_phi_branch != 0) met_T1CHS_fromCORE_phi();
  if (met_T1CHS_miniAOD_CORE_pt_branch != 0) met_T1CHS_miniAOD_CORE_pt();
  if (met_T1CHS_miniAOD_CORE_phi_branch != 0) met_T1CHS_miniAOD_CORE_phi();
  if (met_T1CHS_miniAOD_CORE_up_pt_branch != 0) met_T1CHS_miniAOD_CORE_up_pt();
  if (met_T1CHS_miniAOD_CORE_up_phi_branch != 0) met_T1CHS_miniAOD_CORE_up_phi();
  if (met_T1CHS_miniAOD_CORE_dn_pt_branch != 0) met_T1CHS_miniAOD_CORE_dn_pt();
  if (met_T1CHS_miniAOD_CORE_dn_phi_branch != 0) met_T1CHS_miniAOD_CORE_dn_phi();
  if (decayedphoton_lep1_p4_branch != 0) decayedphoton_lep1_p4();
  if (decayedphoton_lep2_p4_branch != 0) decayedphoton_lep2_p4();
  if (decayedphoton_bosn_p4_branch != 0) decayedphoton_bosn_p4();
  if (decayedphoton_mt2_branch != 0) decayedphoton_mt2();
  if (hyp_type_branch != 0) hyp_type();
  if (evt_type_branch != 0) evt_type();
  if (mass_gluino_branch != 0) mass_gluino();
  if (mass_LSP_branch != 0) mass_LSP();
  if (mass_chi_branch != 0) mass_chi();
  if (isrboost_branch != 0) isrboost();
  if (isr_njets_branch != 0) isr_njets();
  if (isr_weight_branch != 0) isr_weight();
  if (isr_unc_branch != 0) isr_unc();
  if (weightsf_lepid_branch != 0) weightsf_lepid();
  if (weightsf_lepiso_branch != 0) weightsf_lepiso();
  if (weightsf_lepip_branch != 0) weightsf_lepip();
  if (weightsf_lepreco_branch != 0) weightsf_lepreco();
  if (weightsf_lepconv_branch != 0) weightsf_lepconv();
  if (weightsf_lepid_FS_branch != 0) weightsf_lepid_FS();
  if (weightsf_lepiso_FS_branch != 0) weightsf_lepiso_FS();
  if (weightsf_lepip_FS_branch != 0) weightsf_lepip_FS();
  if (weight_rn_r1_n1_branch != 0) weight_rn_r1_n1();
  if (weight_rn_r1_n2_branch != 0) weight_rn_r1_n2();
  if (weight_rn_r1_np05_branch != 0) weight_rn_r1_np05();
  if (weight_rn_r2_n1_branch != 0) weight_rn_r2_n1();
  if (weight_rn_r2_n2_branch != 0) weight_rn_r2_n2();
  if (weight_rn_r2_n0p5_branch != 0) weight_rn_r2_n0p5();
  if (weight_rn_r0p5_n1_branch != 0) weight_rn_r0p5_n1();
  if (weight_rn_r0p5_n2_branch != 0) weight_rn_r0p5_n2();
  if (weight_rn_r0p5_n0p5_branch != 0) weight_rn_r0p5_n0p5();
  if (weight_pdf_up_branch != 0) weight_pdf_up();
  if (weight_pdf_down_branch != 0) weight_pdf_down();
  if (weight_alphas_down_branch != 0) weight_alphas_down();
  if (weight_alphas_up_branch != 0) weight_alphas_up();
  if (LooseElec_branch != 0) LooseElec();
  if (LooseLepton_branch != 0) LooseLepton();
  if (LooseMuon_branch != 0) LooseMuon();
  if (SignalElec_branch != 0) SignalElec();
  if (SignalLepton_branch != 0) SignalLepton();
  if (SignalMuon_branch != 0) SignalMuon();
  if (Tight3lElec_branch != 0) Tight3lElec();
  if (Tight3lLepton_branch != 0) Tight3lLepton();
  if (TightElec_branch != 0) TightElec();
  if (TightLepton_branch != 0) TightLepton();
  if (TightMuon_branch != 0) TightMuon();
  if (VetoElec_branch != 0) VetoElec();
  if (VetoLepton_branch != 0) VetoLepton();
  if (VetoMuon_branch != 0) VetoMuon();
  if (Good3LJet_branch != 0) Good3LJet();
  if (GoodSSJet_branch != 0) GoodSSJet();
  if (GoodSSWJet_branch != 0) GoodSSWJet();
  if (LooseBJet_branch != 0) LooseBJet();
}

const int &WWWTree::run() {
  if (not run_isLoaded) {
    if (run_branch != 0) {
      run_branch->GetEntry(index);
    } else {
      printf("branch run_branch does not exist!\n");
      exit(1);
    }
    run_isLoaded = true;
  }
  return run_;
}

const int &WWWTree::lumi() {
  if (not lumi_isLoaded) {
    if (lumi_branch != 0) {
      lumi_branch->GetEntry(index);
    } else {
      printf("branch lumi_branch does not exist!\n");
      exit(1);
    }
    lumi_isLoaded = true;
  }
  return lumi_;
}

const unsigned long long &WWWTree::evt() {
  if (not evt_isLoaded) {
    if (evt_branch != 0) {
      evt_branch->GetEntry(index);
    } else {
      printf("branch evt_branch does not exist!\n");
      exit(1);
    }
    evt_isLoaded = true;
  }
  return evt_;
}

const int &WWWTree::isData() {
  if (not isData_isLoaded) {
    if (isData_branch != 0) {
      isData_branch->GetEntry(index);
    } else {
      printf("branch isData_branch does not exist!\n");
      exit(1);
    }
    isData_isLoaded = true;
  }
  return isData_;
}

const bool &WWWTree::evt_passgoodrunlist() {
  if (not evt_passgoodrunlist_isLoaded) {
    if (evt_passgoodrunlist_branch != 0) {
      evt_passgoodrunlist_branch->GetEntry(index);
    } else {
      printf("branch evt_passgoodrunlist_branch does not exist!\n");
      exit(1);
    }
    evt_passgoodrunlist_isLoaded = true;
  }
  return evt_passgoodrunlist_;
}

const float &WWWTree::evt_scale1fb() {
  if (not evt_scale1fb_isLoaded) {
    if (evt_scale1fb_branch != 0) {
      evt_scale1fb_branch->GetEntry(index);
    } else {
      printf("branch evt_scale1fb_branch does not exist!\n");
      exit(1);
    }
    evt_scale1fb_isLoaded = true;
  }
  return evt_scale1fb_;
}

const float &WWWTree::evt_xsec() {
  if (not evt_xsec_isLoaded) {
    if (evt_xsec_branch != 0) {
      evt_xsec_branch->GetEntry(index);
    } else {
      printf("branch evt_xsec_branch does not exist!\n");
      exit(1);
    }
    evt_xsec_isLoaded = true;
  }
  return evt_xsec_;
}

const float &WWWTree::evt_kfactor() {
  if (not evt_kfactor_isLoaded) {
    if (evt_kfactor_branch != 0) {
      evt_kfactor_branch->GetEntry(index);
    } else {
      printf("branch evt_kfactor_branch does not exist!\n");
      exit(1);
    }
    evt_kfactor_isLoaded = true;
  }
  return evt_kfactor_;
}

const float &WWWTree::evt_filter() {
  if (not evt_filter_isLoaded) {
    if (evt_filter_branch != 0) {
      evt_filter_branch->GetEntry(index);
    } else {
      printf("branch evt_filter_branch does not exist!\n");
      exit(1);
    }
    evt_filter_isLoaded = true;
  }
  return evt_filter_;
}

const int &WWWTree::evt_nEvts() {
  if (not evt_nEvts_isLoaded) {
    if (evt_nEvts_branch != 0) {
      evt_nEvts_branch->GetEntry(index);
    } else {
      printf("branch evt_nEvts_branch does not exist!\n");
      exit(1);
    }
    evt_nEvts_isLoaded = true;
  }
  return evt_nEvts_;
}

const vector<TString> &WWWTree::evt_dataset() {
  if (not evt_dataset_isLoaded) {
    if (evt_dataset_branch != 0) {
      evt_dataset_branch->GetEntry(index);
    } else {
      printf("branch evt_dataset_branch does not exist!\n");
      exit(1);
    }
    evt_dataset_isLoaded = true;
  }
  return *evt_dataset_;
}

const float &WWWTree::puWeight() {
  if (not puWeight_isLoaded) {
    if (puWeight_branch != 0) {
      puWeight_branch->GetEntry(index);
    } else {
      printf("branch puWeight_branch does not exist!\n");
      exit(1);
    }
    puWeight_isLoaded = true;
  }
  return puWeight_;
}

const int &WWWTree::nVert() {
  if (not nVert_isLoaded) {
    if (nVert_branch != 0) {
      nVert_branch->GetEntry(index);
    } else {
      printf("branch nVert_branch does not exist!\n");
      exit(1);
    }
    nVert_isLoaded = true;
  }
  return nVert_;
}

const int &WWWTree::firstgoodvertex() {
  if (not firstgoodvertex_isLoaded) {
    if (firstgoodvertex_branch != 0) {
      firstgoodvertex_branch->GetEntry(index);
    } else {
      printf("branch firstgoodvertex_branch does not exist!\n");
      exit(1);
    }
    firstgoodvertex_isLoaded = true;
  }
  return firstgoodvertex_;
}

const int &WWWTree::nTrueInt() {
  if (not nTrueInt_isLoaded) {
    if (nTrueInt_branch != 0) {
      nTrueInt_branch->GetEntry(index);
    } else {
      printf("branch nTrueInt_branch does not exist!\n");
      exit(1);
    }
    nTrueInt_isLoaded = true;
  }
  return nTrueInt_;
}

const float &WWWTree::rho() {
  if (not rho_isLoaded) {
    if (rho_branch != 0) {
      rho_branch->GetEntry(index);
    } else {
      printf("branch rho_branch does not exist!\n");
      exit(1);
    }
    rho_isLoaded = true;
  }
  return rho_;
}

const float &WWWTree::rho25() {
  if (not rho25_isLoaded) {
    if (rho25_branch != 0) {
      rho25_branch->GetEntry(index);
    } else {
      printf("branch rho25_branch does not exist!\n");
      exit(1);
    }
    rho25_isLoaded = true;
  }
  return rho25_;
}

const float &WWWTree::gen_ht() {
  if (not gen_ht_isLoaded) {
    if (gen_ht_branch != 0) {
      gen_ht_branch->GetEntry(index);
    } else {
      printf("branch gen_ht_branch does not exist!\n");
      exit(1);
    }
    gen_ht_isLoaded = true;
  }
  return gen_ht_;
}

const int &WWWTree::nBJetTight() {
  if (not nBJetTight_isLoaded) {
    if (nBJetTight_branch != 0) {
      nBJetTight_branch->GetEntry(index);
    } else {
      printf("branch nBJetTight_branch does not exist!\n");
      exit(1);
    }
    nBJetTight_isLoaded = true;
  }
  return nBJetTight_;
}

const int &WWWTree::nBJetMedium() {
  if (not nBJetMedium_isLoaded) {
    if (nBJetMedium_branch != 0) {
      nBJetMedium_branch->GetEntry(index);
    } else {
      printf("branch nBJetMedium_branch does not exist!\n");
      exit(1);
    }
    nBJetMedium_isLoaded = true;
  }
  return nBJetMedium_;
}

const int &WWWTree::nBJetLoose() {
  if (not nBJetLoose_isLoaded) {
    if (nBJetLoose_branch != 0) {
      nBJetLoose_branch->GetEntry(index);
    } else {
      printf("branch nBJetLoose_branch does not exist!\n");
      exit(1);
    }
    nBJetLoose_isLoaded = true;
  }
  return nBJetLoose_;
}

const int &WWWTree::nBJetTight_up() {
  if (not nBJetTight_up_isLoaded) {
    if (nBJetTight_up_branch != 0) {
      nBJetTight_up_branch->GetEntry(index);
    } else {
      printf("branch nBJetTight_up_branch does not exist!\n");
      exit(1);
    }
    nBJetTight_up_isLoaded = true;
  }
  return nBJetTight_up_;
}

const int &WWWTree::nBJetMedium_up() {
  if (not nBJetMedium_up_isLoaded) {
    if (nBJetMedium_up_branch != 0) {
      nBJetMedium_up_branch->GetEntry(index);
    } else {
      printf("branch nBJetMedium_up_branch does not exist!\n");
      exit(1);
    }
    nBJetMedium_up_isLoaded = true;
  }
  return nBJetMedium_up_;
}

const int &WWWTree::nBJetLoose_up() {
  if (not nBJetLoose_up_isLoaded) {
    if (nBJetLoose_up_branch != 0) {
      nBJetLoose_up_branch->GetEntry(index);
    } else {
      printf("branch nBJetLoose_up_branch does not exist!\n");
      exit(1);
    }
    nBJetLoose_up_isLoaded = true;
  }
  return nBJetLoose_up_;
}

const int &WWWTree::nBJetTight_dn() {
  if (not nBJetTight_dn_isLoaded) {
    if (nBJetTight_dn_branch != 0) {
      nBJetTight_dn_branch->GetEntry(index);
    } else {
      printf("branch nBJetTight_dn_branch does not exist!\n");
      exit(1);
    }
    nBJetTight_dn_isLoaded = true;
  }
  return nBJetTight_dn_;
}

const int &WWWTree::nBJetMedium_dn() {
  if (not nBJetMedium_dn_isLoaded) {
    if (nBJetMedium_dn_branch != 0) {
      nBJetMedium_dn_branch->GetEntry(index);
    } else {
      printf("branch nBJetMedium_dn_branch does not exist!\n");
      exit(1);
    }
    nBJetMedium_dn_isLoaded = true;
  }
  return nBJetMedium_dn_;
}

const int &WWWTree::nBJetLoose_dn() {
  if (not nBJetLoose_dn_isLoaded) {
    if (nBJetLoose_dn_branch != 0) {
      nBJetLoose_dn_branch->GetEntry(index);
    } else {
      printf("branch nBJetLoose_dn_branch does not exist!\n");
      exit(1);
    }
    nBJetLoose_dn_isLoaded = true;
  }
  return nBJetLoose_dn_;
}

const int &WWWTree::nJet200MuFrac50DphiMet() {
  if (not nJet200MuFrac50DphiMet_isLoaded) {
    if (nJet200MuFrac50DphiMet_branch != 0) {
      nJet200MuFrac50DphiMet_branch->GetEntry(index);
    } else {
      printf("branch nJet200MuFrac50DphiMet_branch does not exist!\n");
      exit(1);
    }
    nJet200MuFrac50DphiMet_isLoaded = true;
  }
  return nJet200MuFrac50DphiMet_;
}

const int &WWWTree::nMuons10() {
  if (not nMuons10_isLoaded) {
    if (nMuons10_branch != 0) {
      nMuons10_branch->GetEntry(index);
    } else {
      printf("branch nMuons10_branch does not exist!\n");
      exit(1);
    }
    nMuons10_isLoaded = true;
  }
  return nMuons10_;
}

const int &WWWTree::nBadMuons20() {
  if (not nBadMuons20_isLoaded) {
    if (nBadMuons20_branch != 0) {
      nBadMuons20_branch->GetEntry(index);
    } else {
      printf("branch nBadMuons20_branch does not exist!\n");
      exit(1);
    }
    nBadMuons20_isLoaded = true;
  }
  return nBadMuons20_;
}

const int &WWWTree::nElectrons10() {
  if (not nElectrons10_isLoaded) {
    if (nElectrons10_branch != 0) {
      nElectrons10_branch->GetEntry(index);
    } else {
      printf("branch nElectrons10_branch does not exist!\n");
      exit(1);
    }
    nElectrons10_isLoaded = true;
  }
  return nElectrons10_;
}

const int &WWWTree::nGammas20() {
  if (not nGammas20_isLoaded) {
    if (nGammas20_branch != 0) {
      nGammas20_branch->GetEntry(index);
    } else {
      printf("branch nGammas20_branch does not exist!\n");
      exit(1);
    }
    nGammas20_isLoaded = true;
  }
  return nGammas20_;
}

const int &WWWTree::nTaus20() {
  if (not nTaus20_isLoaded) {
    if (nTaus20_branch != 0) {
      nTaus20_branch->GetEntry(index);
    } else {
      printf("branch nTaus20_branch does not exist!\n");
      exit(1);
    }
    nTaus20_isLoaded = true;
  }
  return nTaus20_;
}

const float &WWWTree::met_pt() {
  if (not met_pt_isLoaded) {
    if (met_pt_branch != 0) {
      met_pt_branch->GetEntry(index);
    } else {
      printf("branch met_pt_branch does not exist!\n");
      exit(1);
    }
    met_pt_isLoaded = true;
  }
  return met_pt_;
}

const float &WWWTree::met_phi() {
  if (not met_phi_isLoaded) {
    if (met_phi_branch != 0) {
      met_phi_branch->GetEntry(index);
    } else {
      printf("branch met_phi_branch does not exist!\n");
      exit(1);
    }
    met_phi_isLoaded = true;
  }
  return met_phi_;
}

const float &WWWTree::met_calo_pt() {
  if (not met_calo_pt_isLoaded) {
    if (met_calo_pt_branch != 0) {
      met_calo_pt_branch->GetEntry(index);
    } else {
      printf("branch met_calo_pt_branch does not exist!\n");
      exit(1);
    }
    met_calo_pt_isLoaded = true;
  }
  return met_calo_pt_;
}

const float &WWWTree::met_calo_phi() {
  if (not met_calo_phi_isLoaded) {
    if (met_calo_phi_branch != 0) {
      met_calo_phi_branch->GetEntry(index);
    } else {
      printf("branch met_calo_phi_branch does not exist!\n");
      exit(1);
    }
    met_calo_phi_isLoaded = true;
  }
  return met_calo_phi_;
}

const float &WWWTree::met_miniaod_pt() {
  if (not met_miniaod_pt_isLoaded) {
    if (met_miniaod_pt_branch != 0) {
      met_miniaod_pt_branch->GetEntry(index);
    } else {
      printf("branch met_miniaod_pt_branch does not exist!\n");
      exit(1);
    }
    met_miniaod_pt_isLoaded = true;
  }
  return met_miniaod_pt_;
}

const float &WWWTree::met_miniaod_phi() {
  if (not met_miniaod_phi_isLoaded) {
    if (met_miniaod_phi_branch != 0) {
      met_miniaod_phi_branch->GetEntry(index);
    } else {
      printf("branch met_miniaod_phi_branch does not exist!\n");
      exit(1);
    }
    met_miniaod_phi_isLoaded = true;
  }
  return met_miniaod_phi_;
}

const float &WWWTree::met_muegclean_pt() {
  if (not met_muegclean_pt_isLoaded) {
    if (met_muegclean_pt_branch != 0) {
      met_muegclean_pt_branch->GetEntry(index);
    } else {
      printf("branch met_muegclean_pt_branch does not exist!\n");
      exit(1);
    }
    met_muegclean_pt_isLoaded = true;
  }
  return met_muegclean_pt_;
}

const float &WWWTree::met_muegclean_phi() {
  if (not met_muegclean_phi_isLoaded) {
    if (met_muegclean_phi_branch != 0) {
      met_muegclean_phi_branch->GetEntry(index);
    } else {
      printf("branch met_muegclean_phi_branch does not exist!\n");
      exit(1);
    }
    met_muegclean_phi_isLoaded = true;
  }
  return met_muegclean_phi_;
}

const float &WWWTree::met_rawPt() {
  if (not met_rawPt_isLoaded) {
    if (met_rawPt_branch != 0) {
      met_rawPt_branch->GetEntry(index);
    } else {
      printf("branch met_rawPt_branch does not exist!\n");
      exit(1);
    }
    met_rawPt_isLoaded = true;
  }
  return met_rawPt_;
}

const float &WWWTree::met_rawPhi() {
  if (not met_rawPhi_isLoaded) {
    if (met_rawPhi_branch != 0) {
      met_rawPhi_branch->GetEntry(index);
    } else {
      printf("branch met_rawPhi_branch does not exist!\n");
      exit(1);
    }
    met_rawPhi_isLoaded = true;
  }
  return met_rawPhi_;
}

const float &WWWTree::met_genPt() {
  if (not met_genPt_isLoaded) {
    if (met_genPt_branch != 0) {
      met_genPt_branch->GetEntry(index);
    } else {
      printf("branch met_genPt_branch does not exist!\n");
      exit(1);
    }
    met_genPt_isLoaded = true;
  }
  return met_genPt_;
}

const float &WWWTree::met_genPhi() {
  if (not met_genPhi_isLoaded) {
    if (met_genPhi_branch != 0) {
      met_genPhi_branch->GetEntry(index);
    } else {
      printf("branch met_genPhi_branch does not exist!\n");
      exit(1);
    }
    met_genPhi_isLoaded = true;
  }
  return met_genPhi_;
}

const float &WWWTree::sumet_raw() {
  if (not sumet_raw_isLoaded) {
    if (sumet_raw_branch != 0) {
      sumet_raw_branch->GetEntry(index);
    } else {
      printf("branch sumet_raw_branch does not exist!\n");
      exit(1);
    }
    sumet_raw_isLoaded = true;
  }
  return sumet_raw_;
}

const int &WWWTree::Flag_ecalLaserCorrFilter() {
  if (not Flag_ecalLaserCorrFilter_isLoaded) {
    if (Flag_ecalLaserCorrFilter_branch != 0) {
      Flag_ecalLaserCorrFilter_branch->GetEntry(index);
    } else {
      printf("branch Flag_ecalLaserCorrFilter_branch does not exist!\n");
      exit(1);
    }
    Flag_ecalLaserCorrFilter_isLoaded = true;
  }
  return Flag_ecalLaserCorrFilter_;
}

const int &WWWTree::Flag_hcalLaserEventFilter() {
  if (not Flag_hcalLaserEventFilter_isLoaded) {
    if (Flag_hcalLaserEventFilter_branch != 0) {
      Flag_hcalLaserEventFilter_branch->GetEntry(index);
    } else {
      printf("branch Flag_hcalLaserEventFilter_branch does not exist!\n");
      exit(1);
    }
    Flag_hcalLaserEventFilter_isLoaded = true;
  }
  return Flag_hcalLaserEventFilter_;
}

const int &WWWTree::Flag_trackingFailureFilter() {
  if (not Flag_trackingFailureFilter_isLoaded) {
    if (Flag_trackingFailureFilter_branch != 0) {
      Flag_trackingFailureFilter_branch->GetEntry(index);
    } else {
      printf("branch Flag_trackingFailureFilter_branch does not exist!\n");
      exit(1);
    }
    Flag_trackingFailureFilter_isLoaded = true;
  }
  return Flag_trackingFailureFilter_;
}

const int &WWWTree::Flag_CSCTightHaloFilter() {
  if (not Flag_CSCTightHaloFilter_isLoaded) {
    if (Flag_CSCTightHaloFilter_branch != 0) {
      Flag_CSCTightHaloFilter_branch->GetEntry(index);
    } else {
      printf("branch Flag_CSCTightHaloFilter_branch does not exist!\n");
      exit(1);
    }
    Flag_CSCTightHaloFilter_isLoaded = true;
  }
  return Flag_CSCTightHaloFilter_;
}

const int &WWWTree::Flag_HBHENoiseFilter() {
  if (not Flag_HBHENoiseFilter_isLoaded) {
    if (Flag_HBHENoiseFilter_branch != 0) {
      Flag_HBHENoiseFilter_branch->GetEntry(index);
    } else {
      printf("branch Flag_HBHENoiseFilter_branch does not exist!\n");
      exit(1);
    }
    Flag_HBHENoiseFilter_isLoaded = true;
  }
  return Flag_HBHENoiseFilter_;
}

const int &WWWTree::Flag_HBHEIsoNoiseFilter() {
  if (not Flag_HBHEIsoNoiseFilter_isLoaded) {
    if (Flag_HBHEIsoNoiseFilter_branch != 0) {
      Flag_HBHEIsoNoiseFilter_branch->GetEntry(index);
    } else {
      printf("branch Flag_HBHEIsoNoiseFilter_branch does not exist!\n");
      exit(1);
    }
    Flag_HBHEIsoNoiseFilter_isLoaded = true;
  }
  return Flag_HBHEIsoNoiseFilter_;
}

const int &WWWTree::Flag_CSCTightHalo2015Filter() {
  if (not Flag_CSCTightHalo2015Filter_isLoaded) {
    if (Flag_CSCTightHalo2015Filter_branch != 0) {
      Flag_CSCTightHalo2015Filter_branch->GetEntry(index);
    } else {
      printf("branch Flag_CSCTightHalo2015Filter_branch does not exist!\n");
      exit(1);
    }
    Flag_CSCTightHalo2015Filter_isLoaded = true;
  }
  return Flag_CSCTightHalo2015Filter_;
}

const int &WWWTree::Flag_EcalDeadCellTriggerPrimitiveFilter() {
  if (not Flag_EcalDeadCellTriggerPrimitiveFilter_isLoaded) {
    if (Flag_EcalDeadCellTriggerPrimitiveFilter_branch != 0) {
      Flag_EcalDeadCellTriggerPrimitiveFilter_branch->GetEntry(index);
    } else {
      printf("branch Flag_EcalDeadCellTriggerPrimitiveFilter_branch does not exist!\n");
      exit(1);
    }
    Flag_EcalDeadCellTriggerPrimitiveFilter_isLoaded = true;
  }
  return Flag_EcalDeadCellTriggerPrimitiveFilter_;
}

const int &WWWTree::Flag_goodVertices() {
  if (not Flag_goodVertices_isLoaded) {
    if (Flag_goodVertices_branch != 0) {
      Flag_goodVertices_branch->GetEntry(index);
    } else {
      printf("branch Flag_goodVertices_branch does not exist!\n");
      exit(1);
    }
    Flag_goodVertices_isLoaded = true;
  }
  return Flag_goodVertices_;
}

const int &WWWTree::Flag_eeBadScFilter() {
  if (not Flag_eeBadScFilter_isLoaded) {
    if (Flag_eeBadScFilter_branch != 0) {
      Flag_eeBadScFilter_branch->GetEntry(index);
    } else {
      printf("branch Flag_eeBadScFilter_branch does not exist!\n");
      exit(1);
    }
    Flag_eeBadScFilter_isLoaded = true;
  }
  return Flag_eeBadScFilter_;
}

const int &WWWTree::Flag_globalTightHalo2016() {
  if (not Flag_globalTightHalo2016_isLoaded) {
    if (Flag_globalTightHalo2016_branch != 0) {
      Flag_globalTightHalo2016_branch->GetEntry(index);
    } else {
      printf("branch Flag_globalTightHalo2016_branch does not exist!\n");
      exit(1);
    }
    Flag_globalTightHalo2016_isLoaded = true;
  }
  return Flag_globalTightHalo2016_;
}

const int &WWWTree::Flag_badMuonFilter() {
  if (not Flag_badMuonFilter_isLoaded) {
    if (Flag_badMuonFilter_branch != 0) {
      Flag_badMuonFilter_branch->GetEntry(index);
    } else {
      printf("branch Flag_badMuonFilter_branch does not exist!\n");
      exit(1);
    }
    Flag_badMuonFilter_isLoaded = true;
  }
  return Flag_badMuonFilter_;
}

const int &WWWTree::Flag_badChargedCandidateFilter() {
  if (not Flag_badChargedCandidateFilter_isLoaded) {
    if (Flag_badChargedCandidateFilter_branch != 0) {
      Flag_badChargedCandidateFilter_branch->GetEntry(index);
    } else {
      printf("branch Flag_badChargedCandidateFilter_branch does not exist!\n");
      exit(1);
    }
    Flag_badChargedCandidateFilter_isLoaded = true;
  }
  return Flag_badChargedCandidateFilter_;
}

const int &WWWTree::Flag_badMuonFilterv2() {
  if (not Flag_badMuonFilterv2_isLoaded) {
    if (Flag_badMuonFilterv2_branch != 0) {
      Flag_badMuonFilterv2_branch->GetEntry(index);
    } else {
      printf("branch Flag_badMuonFilterv2_branch does not exist!\n");
      exit(1);
    }
    Flag_badMuonFilterv2_isLoaded = true;
  }
  return Flag_badMuonFilterv2_;
}

const int &WWWTree::Flag_badChargedCandidateFilterv2() {
  if (not Flag_badChargedCandidateFilterv2_isLoaded) {
    if (Flag_badChargedCandidateFilterv2_branch != 0) {
      Flag_badChargedCandidateFilterv2_branch->GetEntry(index);
    } else {
      printf("branch Flag_badChargedCandidateFilterv2_branch does not exist!\n");
      exit(1);
    }
    Flag_badChargedCandidateFilterv2_isLoaded = true;
  }
  return Flag_badChargedCandidateFilterv2_;
}

const int &WWWTree::Flag_badMuons() {
  if (not Flag_badMuons_isLoaded) {
    if (Flag_badMuons_branch != 0) {
      Flag_badMuons_branch->GetEntry(index);
    } else {
      printf("branch Flag_badMuons_branch does not exist!\n");
      exit(1);
    }
    Flag_badMuons_isLoaded = true;
  }
  return Flag_badMuons_;
}

const int &WWWTree::Flag_duplicateMuons() {
  if (not Flag_duplicateMuons_isLoaded) {
    if (Flag_duplicateMuons_branch != 0) {
      Flag_duplicateMuons_branch->GetEntry(index);
    } else {
      printf("branch Flag_duplicateMuons_branch does not exist!\n");
      exit(1);
    }
    Flag_duplicateMuons_isLoaded = true;
  }
  return Flag_duplicateMuons_;
}

const int &WWWTree::Flag_noBadMuons() {
  if (not Flag_noBadMuons_isLoaded) {
    if (Flag_noBadMuons_branch != 0) {
      Flag_noBadMuons_branch->GetEntry(index);
    } else {
      printf("branch Flag_noBadMuons_branch does not exist!\n");
      exit(1);
    }
    Flag_noBadMuons_isLoaded = true;
  }
  return Flag_noBadMuons_;
}

const int &WWWTree::HLT_singleEl() {
  if (not HLT_singleEl_isLoaded) {
    if (HLT_singleEl_branch != 0) {
      HLT_singleEl_branch->GetEntry(index);
    } else {
      printf("branch HLT_singleEl_branch does not exist!\n");
      exit(1);
    }
    HLT_singleEl_isLoaded = true;
  }
  return HLT_singleEl_;
}

const int &WWWTree::HLT_singleMu() {
  if (not HLT_singleMu_isLoaded) {
    if (HLT_singleMu_branch != 0) {
      HLT_singleMu_branch->GetEntry(index);
    } else {
      printf("branch HLT_singleMu_branch does not exist!\n");
      exit(1);
    }
    HLT_singleMu_isLoaded = true;
  }
  return HLT_singleMu_;
}

const int &WWWTree::HLT_singleMu_noiso() {
  if (not HLT_singleMu_noiso_isLoaded) {
    if (HLT_singleMu_noiso_branch != 0) {
      HLT_singleMu_noiso_branch->GetEntry(index);
    } else {
      printf("branch HLT_singleMu_noiso_branch does not exist!\n");
      exit(1);
    }
    HLT_singleMu_noiso_isLoaded = true;
  }
  return HLT_singleMu_noiso_;
}

const int &WWWTree::HLT_DoubleEl_noiso() {
  if (not HLT_DoubleEl_noiso_isLoaded) {
    if (HLT_DoubleEl_noiso_branch != 0) {
      HLT_DoubleEl_noiso_branch->GetEntry(index);
    } else {
      printf("branch HLT_DoubleEl_noiso_branch does not exist!\n");
      exit(1);
    }
    HLT_DoubleEl_noiso_isLoaded = true;
  }
  return HLT_DoubleEl_noiso_;
}

const int &WWWTree::HLT_DoubleEl() {
  if (not HLT_DoubleEl_isLoaded) {
    if (HLT_DoubleEl_branch != 0) {
      HLT_DoubleEl_branch->GetEntry(index);
    } else {
      printf("branch HLT_DoubleEl_branch does not exist!\n");
      exit(1);
    }
    HLT_DoubleEl_isLoaded = true;
  }
  return HLT_DoubleEl_;
}

const int &WWWTree::HLT_DoubleEl_DZ() {
  if (not HLT_DoubleEl_DZ_isLoaded) {
    if (HLT_DoubleEl_DZ_branch != 0) {
      HLT_DoubleEl_DZ_branch->GetEntry(index);
    } else {
      printf("branch HLT_DoubleEl_DZ_branch does not exist!\n");
      exit(1);
    }
    HLT_DoubleEl_DZ_isLoaded = true;
  }
  return HLT_DoubleEl_DZ_;
}

const int &WWWTree::HLT_DoubleEl_DZ_2() {
  if (not HLT_DoubleEl_DZ_2_isLoaded) {
    if (HLT_DoubleEl_DZ_2_branch != 0) {
      HLT_DoubleEl_DZ_2_branch->GetEntry(index);
    } else {
      printf("branch HLT_DoubleEl_DZ_2_branch does not exist!\n");
      exit(1);
    }
    HLT_DoubleEl_DZ_2_isLoaded = true;
  }
  return HLT_DoubleEl_DZ_2_;
}

const int &WWWTree::HLT_MuEG() {
  if (not HLT_MuEG_isLoaded) {
    if (HLT_MuEG_branch != 0) {
      HLT_MuEG_branch->GetEntry(index);
    } else {
      printf("branch HLT_MuEG_branch does not exist!\n");
      exit(1);
    }
    HLT_MuEG_isLoaded = true;
  }
  return HLT_MuEG_;
}

const int &WWWTree::HLT_MuEG_2() {
  if (not HLT_MuEG_2_isLoaded) {
    if (HLT_MuEG_2_branch != 0) {
      HLT_MuEG_2_branch->GetEntry(index);
    } else {
      printf("branch HLT_MuEG_2_branch does not exist!\n");
      exit(1);
    }
    HLT_MuEG_2_isLoaded = true;
  }
  return HLT_MuEG_2_;
}

const int &WWWTree::HLT_MuEG_noiso() {
  if (not HLT_MuEG_noiso_isLoaded) {
    if (HLT_MuEG_noiso_branch != 0) {
      HLT_MuEG_noiso_branch->GetEntry(index);
    } else {
      printf("branch HLT_MuEG_noiso_branch does not exist!\n");
      exit(1);
    }
    HLT_MuEG_noiso_isLoaded = true;
  }
  return HLT_MuEG_noiso_;
}

const int &WWWTree::HLT_MuEG_noiso_2() {
  if (not HLT_MuEG_noiso_2_isLoaded) {
    if (HLT_MuEG_noiso_2_branch != 0) {
      HLT_MuEG_noiso_2_branch->GetEntry(index);
    } else {
      printf("branch HLT_MuEG_noiso_2_branch does not exist!\n");
      exit(1);
    }
    HLT_MuEG_noiso_2_isLoaded = true;
  }
  return HLT_MuEG_noiso_2_;
}

const int &WWWTree::HLT_Mu8_EG17() {
  if (not HLT_Mu8_EG17_isLoaded) {
    if (HLT_Mu8_EG17_branch != 0) {
      HLT_Mu8_EG17_branch->GetEntry(index);
    } else {
      printf("branch HLT_Mu8_EG17_branch does not exist!\n");
      exit(1);
    }
    HLT_Mu8_EG17_isLoaded = true;
  }
  return HLT_Mu8_EG17_;
}

const int &WWWTree::HLT_Mu8_EG23() {
  if (not HLT_Mu8_EG23_isLoaded) {
    if (HLT_Mu8_EG23_branch != 0) {
      HLT_Mu8_EG23_branch->GetEntry(index);
    } else {
      printf("branch HLT_Mu8_EG23_branch does not exist!\n");
      exit(1);
    }
    HLT_Mu8_EG23_isLoaded = true;
  }
  return HLT_Mu8_EG23_;
}

const int &WWWTree::HLT_Mu8_EG23_DZ() {
  if (not HLT_Mu8_EG23_DZ_isLoaded) {
    if (HLT_Mu8_EG23_DZ_branch != 0) {
      HLT_Mu8_EG23_DZ_branch->GetEntry(index);
    } else {
      printf("branch HLT_Mu8_EG23_DZ_branch does not exist!\n");
      exit(1);
    }
    HLT_Mu8_EG23_DZ_isLoaded = true;
  }
  return HLT_Mu8_EG23_DZ_;
}

const int &WWWTree::HLT_Mu12_EG23_DZ() {
  if (not HLT_Mu12_EG23_DZ_isLoaded) {
    if (HLT_Mu12_EG23_DZ_branch != 0) {
      HLT_Mu12_EG23_DZ_branch->GetEntry(index);
    } else {
      printf("branch HLT_Mu12_EG23_DZ_branch does not exist!\n");
      exit(1);
    }
    HLT_Mu12_EG23_DZ_isLoaded = true;
  }
  return HLT_Mu12_EG23_DZ_;
}

const int &WWWTree::HLT_Mu17_EG12() {
  if (not HLT_Mu17_EG12_isLoaded) {
    if (HLT_Mu17_EG12_branch != 0) {
      HLT_Mu17_EG12_branch->GetEntry(index);
    } else {
      printf("branch HLT_Mu17_EG12_branch does not exist!\n");
      exit(1);
    }
    HLT_Mu17_EG12_isLoaded = true;
  }
  return HLT_Mu17_EG12_;
}

const int &WWWTree::HLT_Mu23_EG8() {
  if (not HLT_Mu23_EG8_isLoaded) {
    if (HLT_Mu23_EG8_branch != 0) {
      HLT_Mu23_EG8_branch->GetEntry(index);
    } else {
      printf("branch HLT_Mu23_EG8_branch does not exist!\n");
      exit(1);
    }
    HLT_Mu23_EG8_isLoaded = true;
  }
  return HLT_Mu23_EG8_;
}

const int &WWWTree::HLT_Mu23_EG8_DZ() {
  if (not HLT_Mu23_EG8_DZ_isLoaded) {
    if (HLT_Mu23_EG8_DZ_branch != 0) {
      HLT_Mu23_EG8_DZ_branch->GetEntry(index);
    } else {
      printf("branch HLT_Mu23_EG8_DZ_branch does not exist!\n");
      exit(1);
    }
    HLT_Mu23_EG8_DZ_isLoaded = true;
  }
  return HLT_Mu23_EG8_DZ_;
}

const int &WWWTree::HLT_Mu23_EG12() {
  if (not HLT_Mu23_EG12_isLoaded) {
    if (HLT_Mu23_EG12_branch != 0) {
      HLT_Mu23_EG12_branch->GetEntry(index);
    } else {
      printf("branch HLT_Mu23_EG12_branch does not exist!\n");
      exit(1);
    }
    HLT_Mu23_EG12_isLoaded = true;
  }
  return HLT_Mu23_EG12_;
}

const int &WWWTree::HLT_Mu23_EG12_DZ() {
  if (not HLT_Mu23_EG12_DZ_isLoaded) {
    if (HLT_Mu23_EG12_DZ_branch != 0) {
      HLT_Mu23_EG12_DZ_branch->GetEntry(index);
    } else {
      printf("branch HLT_Mu23_EG12_DZ_branch does not exist!\n");
      exit(1);
    }
    HLT_Mu23_EG12_DZ_isLoaded = true;
  }
  return HLT_Mu23_EG12_DZ_;
}

const int &WWWTree::HLT_DoubleMu_noiso() {
  if (not HLT_DoubleMu_noiso_isLoaded) {
    if (HLT_DoubleMu_noiso_branch != 0) {
      HLT_DoubleMu_noiso_branch->GetEntry(index);
    } else {
      printf("branch HLT_DoubleMu_noiso_branch does not exist!\n");
      exit(1);
    }
    HLT_DoubleMu_noiso_isLoaded = true;
  }
  return HLT_DoubleMu_noiso_;
}

const int &WWWTree::HLT_DoubleMu_noiso_27_8() {
  if (not HLT_DoubleMu_noiso_27_8_isLoaded) {
    if (HLT_DoubleMu_noiso_27_8_branch != 0) {
      HLT_DoubleMu_noiso_27_8_branch->GetEntry(index);
    } else {
      printf("branch HLT_DoubleMu_noiso_27_8_branch does not exist!\n");
      exit(1);
    }
    HLT_DoubleMu_noiso_27_8_isLoaded = true;
  }
  return HLT_DoubleMu_noiso_27_8_;
}

const int &WWWTree::HLT_DoubleMu_noiso_30_11() {
  if (not HLT_DoubleMu_noiso_30_11_isLoaded) {
    if (HLT_DoubleMu_noiso_30_11_branch != 0) {
      HLT_DoubleMu_noiso_30_11_branch->GetEntry(index);
    } else {
      printf("branch HLT_DoubleMu_noiso_30_11_branch does not exist!\n");
      exit(1);
    }
    HLT_DoubleMu_noiso_30_11_isLoaded = true;
  }
  return HLT_DoubleMu_noiso_30_11_;
}

const int &WWWTree::HLT_DoubleMu_noiso_40_11() {
  if (not HLT_DoubleMu_noiso_40_11_isLoaded) {
    if (HLT_DoubleMu_noiso_40_11_branch != 0) {
      HLT_DoubleMu_noiso_40_11_branch->GetEntry(index);
    } else {
      printf("branch HLT_DoubleMu_noiso_40_11_branch does not exist!\n");
      exit(1);
    }
    HLT_DoubleMu_noiso_40_11_isLoaded = true;
  }
  return HLT_DoubleMu_noiso_40_11_;
}

const int &WWWTree::HLT_DoubleMu() {
  if (not HLT_DoubleMu_isLoaded) {
    if (HLT_DoubleMu_branch != 0) {
      HLT_DoubleMu_branch->GetEntry(index);
    } else {
      printf("branch HLT_DoubleMu_branch does not exist!\n");
      exit(1);
    }
    HLT_DoubleMu_isLoaded = true;
  }
  return HLT_DoubleMu_;
}

const int &WWWTree::HLT_DoubleMu_tk() {
  if (not HLT_DoubleMu_tk_isLoaded) {
    if (HLT_DoubleMu_tk_branch != 0) {
      HLT_DoubleMu_tk_branch->GetEntry(index);
    } else {
      printf("branch HLT_DoubleMu_tk_branch does not exist!\n");
      exit(1);
    }
    HLT_DoubleMu_tk_isLoaded = true;
  }
  return HLT_DoubleMu_tk_;
}

const int &WWWTree::HLT_DoubleMu_dbltk() {
  if (not HLT_DoubleMu_dbltk_isLoaded) {
    if (HLT_DoubleMu_dbltk_branch != 0) {
      HLT_DoubleMu_dbltk_branch->GetEntry(index);
    } else {
      printf("branch HLT_DoubleMu_dbltk_branch does not exist!\n");
      exit(1);
    }
    HLT_DoubleMu_dbltk_isLoaded = true;
  }
  return HLT_DoubleMu_dbltk_;
}

const int &WWWTree::HLT_DoubleMu_nonDZ() {
  if (not HLT_DoubleMu_nonDZ_isLoaded) {
    if (HLT_DoubleMu_nonDZ_branch != 0) {
      HLT_DoubleMu_nonDZ_branch->GetEntry(index);
    } else {
      printf("branch HLT_DoubleMu_nonDZ_branch does not exist!\n");
      exit(1);
    }
    HLT_DoubleMu_nonDZ_isLoaded = true;
  }
  return HLT_DoubleMu_nonDZ_;
}

const int &WWWTree::HLT_DoubleMu_tk_nonDZ() {
  if (not HLT_DoubleMu_tk_nonDZ_isLoaded) {
    if (HLT_DoubleMu_tk_nonDZ_branch != 0) {
      HLT_DoubleMu_tk_nonDZ_branch->GetEntry(index);
    } else {
      printf("branch HLT_DoubleMu_tk_nonDZ_branch does not exist!\n");
      exit(1);
    }
    HLT_DoubleMu_tk_nonDZ_isLoaded = true;
  }
  return HLT_DoubleMu_tk_nonDZ_;
}

const int &WWWTree::HLT_Photon22_R9Id90_HE10_IsoM() {
  if (not HLT_Photon22_R9Id90_HE10_IsoM_isLoaded) {
    if (HLT_Photon22_R9Id90_HE10_IsoM_branch != 0) {
      HLT_Photon22_R9Id90_HE10_IsoM_branch->GetEntry(index);
    } else {
      printf("branch HLT_Photon22_R9Id90_HE10_IsoM_branch does not exist!\n");
      exit(1);
    }
    HLT_Photon22_R9Id90_HE10_IsoM_isLoaded = true;
  }
  return HLT_Photon22_R9Id90_HE10_IsoM_;
}

const int &WWWTree::HLT_Photon30_R9Id90_HE10_IsoM() {
  if (not HLT_Photon30_R9Id90_HE10_IsoM_isLoaded) {
    if (HLT_Photon30_R9Id90_HE10_IsoM_branch != 0) {
      HLT_Photon30_R9Id90_HE10_IsoM_branch->GetEntry(index);
    } else {
      printf("branch HLT_Photon30_R9Id90_HE10_IsoM_branch does not exist!\n");
      exit(1);
    }
    HLT_Photon30_R9Id90_HE10_IsoM_isLoaded = true;
  }
  return HLT_Photon30_R9Id90_HE10_IsoM_;
}

const int &WWWTree::HLT_Photon36_R9Id90_HE10_IsoM() {
  if (not HLT_Photon36_R9Id90_HE10_IsoM_isLoaded) {
    if (HLT_Photon36_R9Id90_HE10_IsoM_branch != 0) {
      HLT_Photon36_R9Id90_HE10_IsoM_branch->GetEntry(index);
    } else {
      printf("branch HLT_Photon36_R9Id90_HE10_IsoM_branch does not exist!\n");
      exit(1);
    }
    HLT_Photon36_R9Id90_HE10_IsoM_isLoaded = true;
  }
  return HLT_Photon36_R9Id90_HE10_IsoM_;
}

const int &WWWTree::HLT_Photon50_R9Id90_HE10_IsoM() {
  if (not HLT_Photon50_R9Id90_HE10_IsoM_isLoaded) {
    if (HLT_Photon50_R9Id90_HE10_IsoM_branch != 0) {
      HLT_Photon50_R9Id90_HE10_IsoM_branch->GetEntry(index);
    } else {
      printf("branch HLT_Photon50_R9Id90_HE10_IsoM_branch does not exist!\n");
      exit(1);
    }
    HLT_Photon50_R9Id90_HE10_IsoM_isLoaded = true;
  }
  return HLT_Photon50_R9Id90_HE10_IsoM_;
}

const int &WWWTree::HLT_Photon75_R9Id90_HE10_IsoM() {
  if (not HLT_Photon75_R9Id90_HE10_IsoM_isLoaded) {
    if (HLT_Photon75_R9Id90_HE10_IsoM_branch != 0) {
      HLT_Photon75_R9Id90_HE10_IsoM_branch->GetEntry(index);
    } else {
      printf("branch HLT_Photon75_R9Id90_HE10_IsoM_branch does not exist!\n");
      exit(1);
    }
    HLT_Photon75_R9Id90_HE10_IsoM_isLoaded = true;
  }
  return HLT_Photon75_R9Id90_HE10_IsoM_;
}

const int &WWWTree::HLT_Photon90_R9Id90_HE10_IsoM() {
  if (not HLT_Photon90_R9Id90_HE10_IsoM_isLoaded) {
    if (HLT_Photon90_R9Id90_HE10_IsoM_branch != 0) {
      HLT_Photon90_R9Id90_HE10_IsoM_branch->GetEntry(index);
    } else {
      printf("branch HLT_Photon90_R9Id90_HE10_IsoM_branch does not exist!\n");
      exit(1);
    }
    HLT_Photon90_R9Id90_HE10_IsoM_isLoaded = true;
  }
  return HLT_Photon90_R9Id90_HE10_IsoM_;
}

const int &WWWTree::HLT_Photon120_R9Id90_HE10_IsoM() {
  if (not HLT_Photon120_R9Id90_HE10_IsoM_isLoaded) {
    if (HLT_Photon120_R9Id90_HE10_IsoM_branch != 0) {
      HLT_Photon120_R9Id90_HE10_IsoM_branch->GetEntry(index);
    } else {
      printf("branch HLT_Photon120_R9Id90_HE10_IsoM_branch does not exist!\n");
      exit(1);
    }
    HLT_Photon120_R9Id90_HE10_IsoM_isLoaded = true;
  }
  return HLT_Photon120_R9Id90_HE10_IsoM_;
}

const int &WWWTree::HLT_Photon165_R9Id90_HE10_IsoM() {
  if (not HLT_Photon165_R9Id90_HE10_IsoM_isLoaded) {
    if (HLT_Photon165_R9Id90_HE10_IsoM_branch != 0) {
      HLT_Photon165_R9Id90_HE10_IsoM_branch->GetEntry(index);
    } else {
      printf("branch HLT_Photon165_R9Id90_HE10_IsoM_branch does not exist!\n");
      exit(1);
    }
    HLT_Photon165_R9Id90_HE10_IsoM_isLoaded = true;
  }
  return HLT_Photon165_R9Id90_HE10_IsoM_;
}

const int &WWWTree::HLT_Photon165_HE10() {
  if (not HLT_Photon165_HE10_isLoaded) {
    if (HLT_Photon165_HE10_branch != 0) {
      HLT_Photon165_HE10_branch->GetEntry(index);
    } else {
      printf("branch HLT_Photon165_HE10_branch does not exist!\n");
      exit(1);
    }
    HLT_Photon165_HE10_isLoaded = true;
  }
  return HLT_Photon165_HE10_;
}

const int &WWWTree::HLT_CaloJet500_NoJetID() {
  if (not HLT_CaloJet500_NoJetID_isLoaded) {
    if (HLT_CaloJet500_NoJetID_branch != 0) {
      HLT_CaloJet500_NoJetID_branch->GetEntry(index);
    } else {
      printf("branch HLT_CaloJet500_NoJetID_branch does not exist!\n");
      exit(1);
    }
    HLT_CaloJet500_NoJetID_isLoaded = true;
  }
  return HLT_CaloJet500_NoJetID_;
}

const int &WWWTree::HLT_ECALHT800_NoJetID() {
  if (not HLT_ECALHT800_NoJetID_isLoaded) {
    if (HLT_ECALHT800_NoJetID_branch != 0) {
      HLT_ECALHT800_NoJetID_branch->GetEntry(index);
    } else {
      printf("branch HLT_ECALHT800_NoJetID_branch does not exist!\n");
      exit(1);
    }
    HLT_ECALHT800_NoJetID_isLoaded = true;
  }
  return HLT_ECALHT800_NoJetID_;
}

const bool &WWWTree::HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton() {
  if (not HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton_isLoaded) {
    if (HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton_branch != 0) {
      HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton_branch->GetEntry(index);
    } else {
      printf("branch HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton_branch does not exist!\n");
      exit(1);
    }
    HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton_isLoaded = true;
  }
  return HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton_;
}

const bool &WWWTree::HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton() {
  if (not HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton_isLoaded) {
    if (HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton_branch != 0) {
      HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton_branch->GetEntry(index);
    } else {
      printf("branch HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton_branch does not exist!\n");
      exit(1);
    }
    HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton_isLoaded = true;
  }
  return HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton_;
}

const bool &WWWTree::HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton() {
  if (not HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton_isLoaded) {
    if (HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton_branch != 0) {
      HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton_branch->GetEntry(index);
    } else {
      printf("branch HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton_branch does not exist!\n");
      exit(1);
    }
    HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton_isLoaded = true;
  }
  return HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton_;
}

const bool &WWWTree::HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton() {
  if (not HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton_isLoaded) {
    if (HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton_branch != 0) {
      HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton_branch->GetEntry(index);
    } else {
      printf("branch HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton_branch does not exist!\n");
      exit(1);
    }
    HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton_isLoaded = true;
  }
  return HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton_;
}

const bool &WWWTree::HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton() {
  if (not HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton_isLoaded) {
    if (HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton_branch != 0) {
      HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton_branch->GetEntry(index);
    } else {
      printf("branch HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton_branch does not exist!\n");
      exit(1);
    }
    HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton_isLoaded = true;
  }
  return HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton_;
}

const bool &WWWTree::HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton() {
  if (not HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton_isLoaded) {
    if (HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton_branch != 0) {
      HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton_branch->GetEntry(index);
    } else {
      printf("branch HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton_branch does not exist!\n");
      exit(1);
    }
    HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton_isLoaded = true;
  }
  return HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton_;
}

const bool &WWWTree::HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton() {
  if (not HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton_isLoaded) {
    if (HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton_branch != 0) {
      HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton_branch->GetEntry(index);
    } else {
      printf("branch HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton_branch does not exist!\n");
      exit(1);
    }
    HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton_isLoaded = true;
  }
  return HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton_;
}

const bool &WWWTree::HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton() {
  if (not HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton_isLoaded) {
    if (HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton_branch != 0) {
      HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton_branch->GetEntry(index);
    } else {
      printf("branch HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton_branch does not exist!\n");
      exit(1);
    }
    HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton_isLoaded = true;
  }
  return HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton_;
}

const bool &WWWTree::HLT_Photon165_HE10_matchedtophoton() {
  if (not HLT_Photon165_HE10_matchedtophoton_isLoaded) {
    if (HLT_Photon165_HE10_matchedtophoton_branch != 0) {
      HLT_Photon165_HE10_matchedtophoton_branch->GetEntry(index);
    } else {
      printf("branch HLT_Photon165_HE10_matchedtophoton_branch does not exist!\n");
      exit(1);
    }
    HLT_Photon165_HE10_matchedtophoton_isLoaded = true;
  }
  return HLT_Photon165_HE10_matchedtophoton_;
}

const float &WWWTree::dilmass() {
  if (not dilmass_isLoaded) {
    if (dilmass_branch != 0) {
      dilmass_branch->GetEntry(index);
    } else {
      printf("branch dilmass_branch does not exist!\n");
      exit(1);
    }
    dilmass_isLoaded = true;
  }
  return dilmass_;
}

const float &WWWTree::dilpt() {
  if (not dilpt_isLoaded) {
    if (dilpt_branch != 0) {
      dilpt_branch->GetEntry(index);
    } else {
      printf("branch dilpt_branch does not exist!\n");
      exit(1);
    }
    dilpt_isLoaded = true;
  }
  return dilpt_;
}

const float &WWWTree::dRll() {
  if (not dRll_isLoaded) {
    if (dRll_branch != 0) {
      dRll_branch->GetEntry(index);
    } else {
      printf("branch dRll_branch does not exist!\n");
      exit(1);
    }
    dRll_isLoaded = true;
  }
  return dRll_;
}

const float &WWWTree::matched_neutralemf() {
  if (not matched_neutralemf_isLoaded) {
    if (matched_neutralemf_branch != 0) {
      matched_neutralemf_branch->GetEntry(index);
    } else {
      printf("branch matched_neutralemf_branch does not exist!\n");
      exit(1);
    }
    matched_neutralemf_isLoaded = true;
  }
  return matched_neutralemf_;
}

const float &WWWTree::matched_emf() {
  if (not matched_emf_isLoaded) {
    if (matched_emf_branch != 0) {
      matched_emf_branch->GetEntry(index);
    } else {
      printf("branch matched_emf_branch does not exist!\n");
      exit(1);
    }
    matched_emf_isLoaded = true;
  }
  return matched_emf_;
}

const bool &WWWTree::elveto() {
  if (not elveto_isLoaded) {
    if (elveto_branch != 0) {
      elveto_branch->GetEntry(index);
    } else {
      printf("branch elveto_branch does not exist!\n");
      exit(1);
    }
    elveto_isLoaded = true;
  }
  return elveto_;
}

const int &WWWTree::nlep() {
  if (not nlep_isLoaded) {
    if (nlep_branch != 0) {
      nlep_branch->GetEntry(index);
    } else {
      printf("branch nlep_branch does not exist!\n");
      exit(1);
    }
    nlep_isLoaded = true;
  }
  return nlep_;
}

const int &WWWTree::nveto_leptons() {
  if (not nveto_leptons_isLoaded) {
    if (nveto_leptons_branch != 0) {
      nveto_leptons_branch->GetEntry(index);
    } else {
      printf("branch nveto_leptons_branch does not exist!\n");
      exit(1);
    }
    nveto_leptons_isLoaded = true;
  }
  return nveto_leptons_;
}

const int &WWWTree::nVetoEl_relIso03EAless01() {
  if (not nVetoEl_relIso03EAless01_isLoaded) {
    if (nVetoEl_relIso03EAless01_branch != 0) {
      nVetoEl_relIso03EAless01_branch->GetEntry(index);
    } else {
      printf("branch nVetoEl_relIso03EAless01_branch does not exist!\n");
      exit(1);
    }
    nVetoEl_relIso03EAless01_isLoaded = true;
  }
  return nVetoEl_relIso03EAless01_;
}

const int &WWWTree::nVetoEl_relIso03EAless02() {
  if (not nVetoEl_relIso03EAless02_isLoaded) {
    if (nVetoEl_relIso03EAless02_branch != 0) {
      nVetoEl_relIso03EAless02_branch->GetEntry(index);
    } else {
      printf("branch nVetoEl_relIso03EAless02_branch does not exist!\n");
      exit(1);
    }
    nVetoEl_relIso03EAless02_isLoaded = true;
  }
  return nVetoEl_relIso03EAless02_;
}

const int &WWWTree::nVetoEl_relIso03EAless03() {
  if (not nVetoEl_relIso03EAless03_isLoaded) {
    if (nVetoEl_relIso03EAless03_branch != 0) {
      nVetoEl_relIso03EAless03_branch->GetEntry(index);
    } else {
      printf("branch nVetoEl_relIso03EAless03_branch does not exist!\n");
      exit(1);
    }
    nVetoEl_relIso03EAless03_isLoaded = true;
  }
  return nVetoEl_relIso03EAless03_;
}

const int &WWWTree::nVetoEl_relIso03EAless04() {
  if (not nVetoEl_relIso03EAless04_isLoaded) {
    if (nVetoEl_relIso03EAless04_branch != 0) {
      nVetoEl_relIso03EAless04_branch->GetEntry(index);
    } else {
      printf("branch nVetoEl_relIso03EAless04_branch does not exist!\n");
      exit(1);
    }
    nVetoEl_relIso03EAless04_isLoaded = true;
  }
  return nVetoEl_relIso03EAless04_;
}

const int &WWWTree::nVetoMu_relIso03EAless01() {
  if (not nVetoMu_relIso03EAless01_isLoaded) {
    if (nVetoMu_relIso03EAless01_branch != 0) {
      nVetoMu_relIso03EAless01_branch->GetEntry(index);
    } else {
      printf("branch nVetoMu_relIso03EAless01_branch does not exist!\n");
      exit(1);
    }
    nVetoMu_relIso03EAless01_isLoaded = true;
  }
  return nVetoMu_relIso03EAless01_;
}

const int &WWWTree::nVetoMu_relIso03EAless02() {
  if (not nVetoMu_relIso03EAless02_isLoaded) {
    if (nVetoMu_relIso03EAless02_branch != 0) {
      nVetoMu_relIso03EAless02_branch->GetEntry(index);
    } else {
      printf("branch nVetoMu_relIso03EAless02_branch does not exist!\n");
      exit(1);
    }
    nVetoMu_relIso03EAless02_isLoaded = true;
  }
  return nVetoMu_relIso03EAless02_;
}

const int &WWWTree::nVetoMu_relIso03EAless03() {
  if (not nVetoMu_relIso03EAless03_isLoaded) {
    if (nVetoMu_relIso03EAless03_branch != 0) {
      nVetoMu_relIso03EAless03_branch->GetEntry(index);
    } else {
      printf("branch nVetoMu_relIso03EAless03_branch does not exist!\n");
      exit(1);
    }
    nVetoMu_relIso03EAless03_isLoaded = true;
  }
  return nVetoMu_relIso03EAless03_;
}

const int &WWWTree::nVetoMu_relIso03EAless04() {
  if (not nVetoMu_relIso03EAless04_isLoaded) {
    if (nVetoMu_relIso03EAless04_branch != 0) {
      nVetoMu_relIso03EAless04_branch->GetEntry(index);
    } else {
      printf("branch nVetoMu_relIso03EAless04_branch does not exist!\n");
      exit(1);
    }
    nVetoMu_relIso03EAless04_isLoaded = true;
  }
  return nVetoMu_relIso03EAless04_;
}

const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &WWWTree::lep_p4() {
  if (not lep_p4_isLoaded) {
    if (lep_p4_branch != 0) {
      lep_p4_branch->GetEntry(index);
    } else {
      printf("branch lep_p4_branch does not exist!\n");
      exit(1);
    }
    lep_p4_isLoaded = true;
  }
  return *lep_p4_;
}

const vector<float> &WWWTree::lep_pt() {
  if (not lep_pt_isLoaded) {
    if (lep_pt_branch != 0) {
      lep_pt_branch->GetEntry(index);
    } else {
      printf("branch lep_pt_branch does not exist!\n");
      exit(1);
    }
    lep_pt_isLoaded = true;
  }
  return *lep_pt_;
}

const vector<float> &WWWTree::lep_eta() {
  if (not lep_eta_isLoaded) {
    if (lep_eta_branch != 0) {
      lep_eta_branch->GetEntry(index);
    } else {
      printf("branch lep_eta_branch does not exist!\n");
      exit(1);
    }
    lep_eta_isLoaded = true;
  }
  return *lep_eta_;
}

const vector<float> &WWWTree::lep_phi() {
  if (not lep_phi_isLoaded) {
    if (lep_phi_branch != 0) {
      lep_phi_branch->GetEntry(index);
    } else {
      printf("branch lep_phi_branch does not exist!\n");
      exit(1);
    }
    lep_phi_isLoaded = true;
  }
  return *lep_phi_;
}

const vector<float> &WWWTree::lep_mass() {
  if (not lep_mass_isLoaded) {
    if (lep_mass_branch != 0) {
      lep_mass_branch->GetEntry(index);
    } else {
      printf("branch lep_mass_branch does not exist!\n");
      exit(1);
    }
    lep_mass_isLoaded = true;
  }
  return *lep_mass_;
}

const vector<int> &WWWTree::lep_charge() {
  if (not lep_charge_isLoaded) {
    if (lep_charge_branch != 0) {
      lep_charge_branch->GetEntry(index);
    } else {
      printf("branch lep_charge_branch does not exist!\n");
      exit(1);
    }
    lep_charge_isLoaded = true;
  }
  return *lep_charge_;
}

const vector<bool> &WWWTree::lep_3ch_agree() {
  if (not lep_3ch_agree_isLoaded) {
    if (lep_3ch_agree_branch != 0) {
      lep_3ch_agree_branch->GetEntry(index);
    } else {
      printf("branch lep_3ch_agree_branch does not exist!\n");
      exit(1);
    }
    lep_3ch_agree_isLoaded = true;
  }
  return *lep_3ch_agree_;
}

const vector<bool> &WWWTree::lep_isFromW() {
  if (not lep_isFromW_isLoaded) {
    if (lep_isFromW_branch != 0) {
      lep_isFromW_branch->GetEntry(index);
    } else {
      printf("branch lep_isFromW_branch does not exist!\n");
      exit(1);
    }
    lep_isFromW_isLoaded = true;
  }
  return *lep_isFromW_;
}

const vector<bool> &WWWTree::lep_isFromZ() {
  if (not lep_isFromZ_isLoaded) {
    if (lep_isFromZ_branch != 0) {
      lep_isFromZ_branch->GetEntry(index);
    } else {
      printf("branch lep_isFromZ_branch does not exist!\n");
      exit(1);
    }
    lep_isFromZ_isLoaded = true;
  }
  return *lep_isFromZ_;
}

const vector<bool> &WWWTree::lep_isFromB() {
  if (not lep_isFromB_isLoaded) {
    if (lep_isFromB_branch != 0) {
      lep_isFromB_branch->GetEntry(index);
    } else {
      printf("branch lep_isFromB_branch does not exist!\n");
      exit(1);
    }
    lep_isFromB_isLoaded = true;
  }
  return *lep_isFromB_;
}

const vector<bool> &WWWTree::lep_isFromC() {
  if (not lep_isFromC_isLoaded) {
    if (lep_isFromC_branch != 0) {
      lep_isFromC_branch->GetEntry(index);
    } else {
      printf("branch lep_isFromC_branch does not exist!\n");
      exit(1);
    }
    lep_isFromC_isLoaded = true;
  }
  return *lep_isFromC_;
}

const vector<bool> &WWWTree::lep_isFromL() {
  if (not lep_isFromL_isLoaded) {
    if (lep_isFromL_branch != 0) {
      lep_isFromL_branch->GetEntry(index);
    } else {
      printf("branch lep_isFromL_branch does not exist!\n");
      exit(1);
    }
    lep_isFromL_isLoaded = true;
  }
  return *lep_isFromL_;
}

const vector<bool> &WWWTree::lep_isFromLF() {
  if (not lep_isFromLF_isLoaded) {
    if (lep_isFromLF_branch != 0) {
      lep_isFromLF_branch->GetEntry(index);
    } else {
      printf("branch lep_isFromLF_branch does not exist!\n");
      exit(1);
    }
    lep_isFromLF_isLoaded = true;
  }
  return *lep_isFromLF_;
}

const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &WWWTree::lep_closest_jet_p4() {
  if (not lep_closest_jet_p4_isLoaded) {
    if (lep_closest_jet_p4_branch != 0) {
      lep_closest_jet_p4_branch->GetEntry(index);
    } else {
      printf("branch lep_closest_jet_p4_branch does not exist!\n");
      exit(1);
    }
    lep_closest_jet_p4_isLoaded = true;
  }
  return *lep_closest_jet_p4_;
}

const vector<double> &WWWTree::lep_ptRatio() {
  if (not lep_ptRatio_isLoaded) {
    if (lep_ptRatio_branch != 0) {
      lep_ptRatio_branch->GetEntry(index);
    } else {
      printf("branch lep_ptRatio_branch does not exist!\n");
      exit(1);
    }
    lep_ptRatio_isLoaded = true;
  }
  return *lep_ptRatio_;
}

const vector<int> &WWWTree::lep_motherIdSS() {
  if (not lep_motherIdSS_isLoaded) {
    if (lep_motherIdSS_branch != 0) {
      lep_motherIdSS_branch->GetEntry(index);
    } else {
      printf("branch lep_motherIdSS_branch does not exist!\n");
      exit(1);
    }
    lep_motherIdSS_isLoaded = true;
  }
  return *lep_motherIdSS_;
}

const vector<int> &WWWTree::lep_genPart_index() {
  if (not lep_genPart_index_isLoaded) {
    if (lep_genPart_index_branch != 0) {
      lep_genPart_index_branch->GetEntry(index);
    } else {
      printf("branch lep_genPart_index_branch does not exist!\n");
      exit(1);
    }
    lep_genPart_index_isLoaded = true;
  }
  return *lep_genPart_index_;
}

const vector<double> &WWWTree::lep_coneCorrPt() {
  if (not lep_coneCorrPt_isLoaded) {
    if (lep_coneCorrPt_branch != 0) {
      lep_coneCorrPt_branch->GetEntry(index);
    } else {
      printf("branch lep_coneCorrPt_branch does not exist!\n");
      exit(1);
    }
    lep_coneCorrPt_isLoaded = true;
  }
  return *lep_coneCorrPt_;
}

const vector<double> &WWWTree::lep_ptRel() {
  if (not lep_ptRel_isLoaded) {
    if (lep_ptRel_branch != 0) {
      lep_ptRel_branch->GetEntry(index);
    } else {
      printf("branch lep_ptRel_branch does not exist!\n");
      exit(1);
    }
    lep_ptRel_isLoaded = true;
  }
  return *lep_ptRel_;
}

const vector<double> &WWWTree::lep_relIso03() {
  if (not lep_relIso03_isLoaded) {
    if (lep_relIso03_branch != 0) {
      lep_relIso03_branch->GetEntry(index);
    } else {
      printf("branch lep_relIso03_branch does not exist!\n");
      exit(1);
    }
    lep_relIso03_isLoaded = true;
  }
  return *lep_relIso03_;
}

const vector<double> &WWWTree::lep_relIso03DB() {
  if (not lep_relIso03DB_isLoaded) {
    if (lep_relIso03DB_branch != 0) {
      lep_relIso03DB_branch->GetEntry(index);
    } else {
      printf("branch lep_relIso03DB_branch does not exist!\n");
      exit(1);
    }
    lep_relIso03DB_isLoaded = true;
  }
  return *lep_relIso03DB_;
}

const vector<double> &WWWTree::lep_relIso03EA() {
  if (not lep_relIso03EA_isLoaded) {
    if (lep_relIso03EA_branch != 0) {
      lep_relIso03EA_branch->GetEntry(index);
    } else {
      printf("branch lep_relIso03EA_branch does not exist!\n");
      exit(1);
    }
    lep_relIso03EA_isLoaded = true;
  }
  return *lep_relIso03EA_;
}

const vector<double> &WWWTree::lep_relIso03EAv2() {
  if (not lep_relIso03EAv2_isLoaded) {
    if (lep_relIso03EAv2_branch != 0) {
      lep_relIso03EAv2_branch->GetEntry(index);
    } else {
      printf("branch lep_relIso03EAv2_branch does not exist!\n");
      exit(1);
    }
    lep_relIso03EAv2_isLoaded = true;
  }
  return *lep_relIso03EAv2_;
}

const vector<double> &WWWTree::lep_relIso04DB() {
  if (not lep_relIso04DB_isLoaded) {
    if (lep_relIso04DB_branch != 0) {
      lep_relIso04DB_branch->GetEntry(index);
    } else {
      printf("branch lep_relIso04DB_branch does not exist!\n");
      exit(1);
    }
    lep_relIso04DB_isLoaded = true;
  }
  return *lep_relIso04DB_;
}

const vector<double> &WWWTree::lep_relIso04EA() {
  if (not lep_relIso04EA_isLoaded) {
    if (lep_relIso04EA_branch != 0) {
      lep_relIso04EA_branch->GetEntry(index);
    } else {
      printf("branch lep_relIso04EA_branch does not exist!\n");
      exit(1);
    }
    lep_relIso04EA_isLoaded = true;
  }
  return *lep_relIso04EA_;
}

const vector<double> &WWWTree::lep_relIso04EAv2() {
  if (not lep_relIso04EAv2_isLoaded) {
    if (lep_relIso04EAv2_branch != 0) {
      lep_relIso04EAv2_branch->GetEntry(index);
    } else {
      printf("branch lep_relIso04EAv2_branch does not exist!\n");
      exit(1);
    }
    lep_relIso04EAv2_isLoaded = true;
  }
  return *lep_relIso04EAv2_;
}

const vector<double> &WWWTree::lep_miniRelIsoCMS3_EA() {
  if (not lep_miniRelIsoCMS3_EA_isLoaded) {
    if (lep_miniRelIsoCMS3_EA_branch != 0) {
      lep_miniRelIsoCMS3_EA_branch->GetEntry(index);
    } else {
      printf("branch lep_miniRelIsoCMS3_EA_branch does not exist!\n");
      exit(1);
    }
    lep_miniRelIsoCMS3_EA_isLoaded = true;
  }
  return *lep_miniRelIsoCMS3_EA_;
}

const vector<double> &WWWTree::lep_miniRelIsoCMS3_EAv2() {
  if (not lep_miniRelIsoCMS3_EAv2_isLoaded) {
    if (lep_miniRelIsoCMS3_EAv2_branch != 0) {
      lep_miniRelIsoCMS3_EAv2_branch->GetEntry(index);
    } else {
      printf("branch lep_miniRelIsoCMS3_EAv2_branch does not exist!\n");
      exit(1);
    }
    lep_miniRelIsoCMS3_EAv2_isLoaded = true;
  }
  return *lep_miniRelIsoCMS3_EAv2_;
}

const vector<double> &WWWTree::lep_miniRelIsoCMS3_DB() {
  if (not lep_miniRelIsoCMS3_DB_isLoaded) {
    if (lep_miniRelIsoCMS3_DB_branch != 0) {
      lep_miniRelIsoCMS3_DB_branch->GetEntry(index);
    } else {
      printf("branch lep_miniRelIsoCMS3_DB_branch does not exist!\n");
      exit(1);
    }
    lep_miniRelIsoCMS3_DB_isLoaded = true;
  }
  return *lep_miniRelIsoCMS3_DB_;
}

const vector<bool> &WWWTree::lep_isTriggerSafenoIso_v1() {
  if (not lep_isTriggerSafenoIso_v1_isLoaded) {
    if (lep_isTriggerSafenoIso_v1_branch != 0) {
      lep_isTriggerSafenoIso_v1_branch->GetEntry(index);
    } else {
      printf("branch lep_isTriggerSafenoIso_v1_branch does not exist!\n");
      exit(1);
    }
    lep_isTriggerSafenoIso_v1_isLoaded = true;
  }
  return *lep_isTriggerSafenoIso_v1_;
}

const vector<bool> &WWWTree::lep_isTriggerSafenoIso_v2() {
  if (not lep_isTriggerSafenoIso_v2_isLoaded) {
    if (lep_isTriggerSafenoIso_v2_branch != 0) {
      lep_isTriggerSafenoIso_v2_branch->GetEntry(index);
    } else {
      printf("branch lep_isTriggerSafenoIso_v2_branch does not exist!\n");
      exit(1);
    }
    lep_isTriggerSafenoIso_v2_isLoaded = true;
  }
  return *lep_isTriggerSafenoIso_v2_;
}

const vector<bool> &WWWTree::lep_isTriggerSafe_v1() {
  if (not lep_isTriggerSafe_v1_isLoaded) {
    if (lep_isTriggerSafe_v1_branch != 0) {
      lep_isTriggerSafe_v1_branch->GetEntry(index);
    } else {
      printf("branch lep_isTriggerSafe_v1_branch does not exist!\n");
      exit(1);
    }
    lep_isTriggerSafe_v1_isLoaded = true;
  }
  return *lep_isTriggerSafe_v1_;
}

const vector<bool> &WWWTree::lep_isTriggerSafe_v2() {
  if (not lep_isTriggerSafe_v2_isLoaded) {
    if (lep_isTriggerSafe_v2_branch != 0) {
      lep_isTriggerSafe_v2_branch->GetEntry(index);
    } else {
      printf("branch lep_isTriggerSafe_v2_branch does not exist!\n");
      exit(1);
    }
    lep_isTriggerSafe_v2_isLoaded = true;
  }
  return *lep_isTriggerSafe_v2_;
}

const vector<bool> &WWWTree::lep_pass_VVV_cutbased_veto() {
  if (not lep_pass_VVV_cutbased_veto_isLoaded) {
    if (lep_pass_VVV_cutbased_veto_branch != 0) {
      lep_pass_VVV_cutbased_veto_branch->GetEntry(index);
    } else {
      printf("branch lep_pass_VVV_cutbased_veto_branch does not exist!\n");
      exit(1);
    }
    lep_pass_VVV_cutbased_veto_isLoaded = true;
  }
  return *lep_pass_VVV_cutbased_veto_;
}

const vector<bool> &WWWTree::lep_pass_VVV_cutbased_veto_noiso() {
  if (not lep_pass_VVV_cutbased_veto_noiso_isLoaded) {
    if (lep_pass_VVV_cutbased_veto_noiso_branch != 0) {
      lep_pass_VVV_cutbased_veto_noiso_branch->GetEntry(index);
    } else {
      printf("branch lep_pass_VVV_cutbased_veto_noiso_branch does not exist!\n");
      exit(1);
    }
    lep_pass_VVV_cutbased_veto_noiso_isLoaded = true;
  }
  return *lep_pass_VVV_cutbased_veto_noiso_;
}

const vector<bool> &WWWTree::lep_pass_VVV_cutbased_veto_noiso_noip() {
  if (not lep_pass_VVV_cutbased_veto_noiso_noip_isLoaded) {
    if (lep_pass_VVV_cutbased_veto_noiso_noip_branch != 0) {
      lep_pass_VVV_cutbased_veto_noiso_noip_branch->GetEntry(index);
    } else {
      printf("branch lep_pass_VVV_cutbased_veto_noiso_noip_branch does not exist!\n");
      exit(1);
    }
    lep_pass_VVV_cutbased_veto_noiso_noip_isLoaded = true;
  }
  return *lep_pass_VVV_cutbased_veto_noiso_noip_;
}

const vector<bool> &WWWTree::lep_pass_VVV_cutbased_veto_bak1() {
  if (not lep_pass_VVV_cutbased_veto_bak1_isLoaded) {
    if (lep_pass_VVV_cutbased_veto_bak1_branch != 0) {
      lep_pass_VVV_cutbased_veto_bak1_branch->GetEntry(index);
    } else {
      printf("branch lep_pass_VVV_cutbased_veto_bak1_branch does not exist!\n");
      exit(1);
    }
    lep_pass_VVV_cutbased_veto_bak1_isLoaded = true;
  }
  return *lep_pass_VVV_cutbased_veto_bak1_;
}

const vector<bool> &WWWTree::lep_pass_VVV_cutbased_veto_noiso_bak1() {
  if (not lep_pass_VVV_cutbased_veto_noiso_bak1_isLoaded) {
    if (lep_pass_VVV_cutbased_veto_noiso_bak1_branch != 0) {
      lep_pass_VVV_cutbased_veto_noiso_bak1_branch->GetEntry(index);
    } else {
      printf("branch lep_pass_VVV_cutbased_veto_noiso_bak1_branch does not exist!\n");
      exit(1);
    }
    lep_pass_VVV_cutbased_veto_noiso_bak1_isLoaded = true;
  }
  return *lep_pass_VVV_cutbased_veto_noiso_bak1_;
}

const vector<bool> &WWWTree::lep_pass_VVV_cutbased_veto_noiso_noip_bak1() {
  if (not lep_pass_VVV_cutbased_veto_noiso_noip_bak1_isLoaded) {
    if (lep_pass_VVV_cutbased_veto_noiso_noip_bak1_branch != 0) {
      lep_pass_VVV_cutbased_veto_noiso_noip_bak1_branch->GetEntry(index);
    } else {
      printf("branch lep_pass_VVV_cutbased_veto_noiso_noip_bak1_branch does not exist!\n");
      exit(1);
    }
    lep_pass_VVV_cutbased_veto_noiso_noip_bak1_isLoaded = true;
  }
  return *lep_pass_VVV_cutbased_veto_noiso_noip_bak1_;
}

const vector<bool> &WWWTree::lep_pass_VVV_cutbased_fo() {
  if (not lep_pass_VVV_cutbased_fo_isLoaded) {
    if (lep_pass_VVV_cutbased_fo_branch != 0) {
      lep_pass_VVV_cutbased_fo_branch->GetEntry(index);
    } else {
      printf("branch lep_pass_VVV_cutbased_fo_branch does not exist!\n");
      exit(1);
    }
    lep_pass_VVV_cutbased_fo_isLoaded = true;
  }
  return *lep_pass_VVV_cutbased_fo_;
}

const vector<bool> &WWWTree::lep_pass_VVV_cutbased_fo_noiso() {
  if (not lep_pass_VVV_cutbased_fo_noiso_isLoaded) {
    if (lep_pass_VVV_cutbased_fo_noiso_branch != 0) {
      lep_pass_VVV_cutbased_fo_noiso_branch->GetEntry(index);
    } else {
      printf("branch lep_pass_VVV_cutbased_fo_noiso_branch does not exist!\n");
      exit(1);
    }
    lep_pass_VVV_cutbased_fo_noiso_isLoaded = true;
  }
  return *lep_pass_VVV_cutbased_fo_noiso_;
}

const vector<bool> &WWWTree::lep_pass_VVV_cutbased_tight_noiso() {
  if (not lep_pass_VVV_cutbased_tight_noiso_isLoaded) {
    if (lep_pass_VVV_cutbased_tight_noiso_branch != 0) {
      lep_pass_VVV_cutbased_tight_noiso_branch->GetEntry(index);
    } else {
      printf("branch lep_pass_VVV_cutbased_tight_noiso_branch does not exist!\n");
      exit(1);
    }
    lep_pass_VVV_cutbased_tight_noiso_isLoaded = true;
  }
  return *lep_pass_VVV_cutbased_tight_noiso_;
}

const vector<bool> &WWWTree::lep_pass_VVV_cutbased_tight() {
  if (not lep_pass_VVV_cutbased_tight_isLoaded) {
    if (lep_pass_VVV_cutbased_tight_branch != 0) {
      lep_pass_VVV_cutbased_tight_branch->GetEntry(index);
    } else {
      printf("branch lep_pass_VVV_cutbased_tight_branch does not exist!\n");
      exit(1);
    }
    lep_pass_VVV_cutbased_tight_isLoaded = true;
  }
  return *lep_pass_VVV_cutbased_tight_;
}

const vector<bool> &WWWTree::lep_pass_VVV_MVAbased_tight_noiso() {
  if (not lep_pass_VVV_MVAbased_tight_noiso_isLoaded) {
    if (lep_pass_VVV_MVAbased_tight_noiso_branch != 0) {
      lep_pass_VVV_MVAbased_tight_noiso_branch->GetEntry(index);
    } else {
      printf("branch lep_pass_VVV_MVAbased_tight_noiso_branch does not exist!\n");
      exit(1);
    }
    lep_pass_VVV_MVAbased_tight_noiso_isLoaded = true;
  }
  return *lep_pass_VVV_MVAbased_tight_noiso_;
}

const vector<bool> &WWWTree::lep_pass_VVV_MVAbased_tight() {
  if (not lep_pass_VVV_MVAbased_tight_isLoaded) {
    if (lep_pass_VVV_MVAbased_tight_branch != 0) {
      lep_pass_VVV_MVAbased_tight_branch->GetEntry(index);
    } else {
      printf("branch lep_pass_VVV_MVAbased_tight_branch does not exist!\n");
      exit(1);
    }
    lep_pass_VVV_MVAbased_tight_isLoaded = true;
  }
  return *lep_pass_VVV_MVAbased_tight_;
}

const vector<bool> &WWWTree::lep_pass_VVV_baseline() {
  if (not lep_pass_VVV_baseline_isLoaded) {
    if (lep_pass_VVV_baseline_branch != 0) {
      lep_pass_VVV_baseline_branch->GetEntry(index);
    } else {
      printf("branch lep_pass_VVV_baseline_branch does not exist!\n");
      exit(1);
    }
    lep_pass_VVV_baseline_isLoaded = true;
  }
  return *lep_pass_VVV_baseline_;
}

const int &WWWTree::nlep_VVV_cutbased_veto() {
  if (not nlep_VVV_cutbased_veto_isLoaded) {
    if (nlep_VVV_cutbased_veto_branch != 0) {
      nlep_VVV_cutbased_veto_branch->GetEntry(index);
    } else {
      printf("branch nlep_VVV_cutbased_veto_branch does not exist!\n");
      exit(1);
    }
    nlep_VVV_cutbased_veto_isLoaded = true;
  }
  return nlep_VVV_cutbased_veto_;
}

const int &WWWTree::nlep_VVV_cutbased_veto_noiso() {
  if (not nlep_VVV_cutbased_veto_noiso_isLoaded) {
    if (nlep_VVV_cutbased_veto_noiso_branch != 0) {
      nlep_VVV_cutbased_veto_noiso_branch->GetEntry(index);
    } else {
      printf("branch nlep_VVV_cutbased_veto_noiso_branch does not exist!\n");
      exit(1);
    }
    nlep_VVV_cutbased_veto_noiso_isLoaded = true;
  }
  return nlep_VVV_cutbased_veto_noiso_;
}

const int &WWWTree::nlep_VVV_cutbased_veto_noiso_noip() {
  if (not nlep_VVV_cutbased_veto_noiso_noip_isLoaded) {
    if (nlep_VVV_cutbased_veto_noiso_noip_branch != 0) {
      nlep_VVV_cutbased_veto_noiso_noip_branch->GetEntry(index);
    } else {
      printf("branch nlep_VVV_cutbased_veto_noiso_noip_branch does not exist!\n");
      exit(1);
    }
    nlep_VVV_cutbased_veto_noiso_noip_isLoaded = true;
  }
  return nlep_VVV_cutbased_veto_noiso_noip_;
}

const int &WWWTree::nlep_VVV_cutbased_veto_bak1() {
  if (not nlep_VVV_cutbased_veto_bak1_isLoaded) {
    if (nlep_VVV_cutbased_veto_bak1_branch != 0) {
      nlep_VVV_cutbased_veto_bak1_branch->GetEntry(index);
    } else {
      printf("branch nlep_VVV_cutbased_veto_bak1_branch does not exist!\n");
      exit(1);
    }
    nlep_VVV_cutbased_veto_bak1_isLoaded = true;
  }
  return nlep_VVV_cutbased_veto_bak1_;
}

const int &WWWTree::nlep_VVV_cutbased_veto_noiso_bak1() {
  if (not nlep_VVV_cutbased_veto_noiso_bak1_isLoaded) {
    if (nlep_VVV_cutbased_veto_noiso_bak1_branch != 0) {
      nlep_VVV_cutbased_veto_noiso_bak1_branch->GetEntry(index);
    } else {
      printf("branch nlep_VVV_cutbased_veto_noiso_bak1_branch does not exist!\n");
      exit(1);
    }
    nlep_VVV_cutbased_veto_noiso_bak1_isLoaded = true;
  }
  return nlep_VVV_cutbased_veto_noiso_bak1_;
}

const int &WWWTree::nlep_VVV_cutbased_veto_noiso_noip_bak1() {
  if (not nlep_VVV_cutbased_veto_noiso_noip_bak1_isLoaded) {
    if (nlep_VVV_cutbased_veto_noiso_noip_bak1_branch != 0) {
      nlep_VVV_cutbased_veto_noiso_noip_bak1_branch->GetEntry(index);
    } else {
      printf("branch nlep_VVV_cutbased_veto_noiso_noip_bak1_branch does not exist!\n");
      exit(1);
    }
    nlep_VVV_cutbased_veto_noiso_noip_bak1_isLoaded = true;
  }
  return nlep_VVV_cutbased_veto_noiso_noip_bak1_;
}

const int &WWWTree::nlep_VVV_cutbased_fo() {
  if (not nlep_VVV_cutbased_fo_isLoaded) {
    if (nlep_VVV_cutbased_fo_branch != 0) {
      nlep_VVV_cutbased_fo_branch->GetEntry(index);
    } else {
      printf("branch nlep_VVV_cutbased_fo_branch does not exist!\n");
      exit(1);
    }
    nlep_VVV_cutbased_fo_isLoaded = true;
  }
  return nlep_VVV_cutbased_fo_;
}

const int &WWWTree::nlep_VVV_cutbased_fo_noiso() {
  if (not nlep_VVV_cutbased_fo_noiso_isLoaded) {
    if (nlep_VVV_cutbased_fo_noiso_branch != 0) {
      nlep_VVV_cutbased_fo_noiso_branch->GetEntry(index);
    } else {
      printf("branch nlep_VVV_cutbased_fo_noiso_branch does not exist!\n");
      exit(1);
    }
    nlep_VVV_cutbased_fo_noiso_isLoaded = true;
  }
  return nlep_VVV_cutbased_fo_noiso_;
}

const int &WWWTree::nlep_VVV_cutbased_tight_noiso() {
  if (not nlep_VVV_cutbased_tight_noiso_isLoaded) {
    if (nlep_VVV_cutbased_tight_noiso_branch != 0) {
      nlep_VVV_cutbased_tight_noiso_branch->GetEntry(index);
    } else {
      printf("branch nlep_VVV_cutbased_tight_noiso_branch does not exist!\n");
      exit(1);
    }
    nlep_VVV_cutbased_tight_noiso_isLoaded = true;
  }
  return nlep_VVV_cutbased_tight_noiso_;
}

const int &WWWTree::nlep_VVV_cutbased_tight() {
  if (not nlep_VVV_cutbased_tight_isLoaded) {
    if (nlep_VVV_cutbased_tight_branch != 0) {
      nlep_VVV_cutbased_tight_branch->GetEntry(index);
    } else {
      printf("branch nlep_VVV_cutbased_tight_branch does not exist!\n");
      exit(1);
    }
    nlep_VVV_cutbased_tight_isLoaded = true;
  }
  return nlep_VVV_cutbased_tight_;
}

const int &WWWTree::nlep_VVV_MVAbased_tight_noiso() {
  if (not nlep_VVV_MVAbased_tight_noiso_isLoaded) {
    if (nlep_VVV_MVAbased_tight_noiso_branch != 0) {
      nlep_VVV_MVAbased_tight_noiso_branch->GetEntry(index);
    } else {
      printf("branch nlep_VVV_MVAbased_tight_noiso_branch does not exist!\n");
      exit(1);
    }
    nlep_VVV_MVAbased_tight_noiso_isLoaded = true;
  }
  return nlep_VVV_MVAbased_tight_noiso_;
}

const int &WWWTree::nlep_VVV_MVAbased_tight() {
  if (not nlep_VVV_MVAbased_tight_isLoaded) {
    if (nlep_VVV_MVAbased_tight_branch != 0) {
      nlep_VVV_MVAbased_tight_branch->GetEntry(index);
    } else {
      printf("branch nlep_VVV_MVAbased_tight_branch does not exist!\n");
      exit(1);
    }
    nlep_VVV_MVAbased_tight_isLoaded = true;
  }
  return nlep_VVV_MVAbased_tight_;
}

const int &WWWTree::nlep_VVV_baseline() {
  if (not nlep_VVV_baseline_isLoaded) {
    if (nlep_VVV_baseline_branch != 0) {
      nlep_VVV_baseline_branch->GetEntry(index);
    } else {
      printf("branch nlep_VVV_baseline_branch does not exist!\n");
      exit(1);
    }
    nlep_VVV_baseline_isLoaded = true;
  }
  return nlep_VVV_baseline_;
}

const vector<int> &WWWTree::lep_pdgId() {
  if (not lep_pdgId_isLoaded) {
    if (lep_pdgId_branch != 0) {
      lep_pdgId_branch->GetEntry(index);
    } else {
      printf("branch lep_pdgId_branch does not exist!\n");
      exit(1);
    }
    lep_pdgId_isLoaded = true;
  }
  return *lep_pdgId_;
}

const vector<int> &WWWTree::lep_mc_Id() {
  if (not lep_mc_Id_isLoaded) {
    if (lep_mc_Id_branch != 0) {
      lep_mc_Id_branch->GetEntry(index);
    } else {
      printf("branch lep_mc_Id_branch does not exist!\n");
      exit(1);
    }
    lep_mc_Id_isLoaded = true;
  }
  return *lep_mc_Id_;
}

const vector<float> &WWWTree::lep_dxy() {
  if (not lep_dxy_isLoaded) {
    if (lep_dxy_branch != 0) {
      lep_dxy_branch->GetEntry(index);
    } else {
      printf("branch lep_dxy_branch does not exist!\n");
      exit(1);
    }
    lep_dxy_isLoaded = true;
  }
  return *lep_dxy_;
}

const vector<float> &WWWTree::lep_ip3d() {
  if (not lep_ip3d_isLoaded) {
    if (lep_ip3d_branch != 0) {
      lep_ip3d_branch->GetEntry(index);
    } else {
      printf("branch lep_ip3d_branch does not exist!\n");
      exit(1);
    }
    lep_ip3d_isLoaded = true;
  }
  return *lep_ip3d_;
}

const vector<float> &WWWTree::lep_ip3derr() {
  if (not lep_ip3derr_isLoaded) {
    if (lep_ip3derr_branch != 0) {
      lep_ip3derr_branch->GetEntry(index);
    } else {
      printf("branch lep_ip3derr_branch does not exist!\n");
      exit(1);
    }
    lep_ip3derr_isLoaded = true;
  }
  return *lep_ip3derr_;
}

const vector<float> &WWWTree::lep_etaSC() {
  if (not lep_etaSC_isLoaded) {
    if (lep_etaSC_branch != 0) {
      lep_etaSC_branch->GetEntry(index);
    } else {
      printf("branch lep_etaSC_branch does not exist!\n");
      exit(1);
    }
    lep_etaSC_isLoaded = true;
  }
  return *lep_etaSC_;
}

const vector<float> &WWWTree::lep_r9() {
  if (not lep_r9_isLoaded) {
    if (lep_r9_branch != 0) {
      lep_r9_branch->GetEntry(index);
    } else {
      printf("branch lep_r9_branch does not exist!\n");
      exit(1);
    }
    lep_r9_isLoaded = true;
  }
  return *lep_r9_;
}

const vector<float> &WWWTree::lep_dz() {
  if (not lep_dz_isLoaded) {
    if (lep_dz_branch != 0) {
      lep_dz_branch->GetEntry(index);
    } else {
      printf("branch lep_dz_branch does not exist!\n");
      exit(1);
    }
    lep_dz_isLoaded = true;
  }
  return *lep_dz_;
}

const vector<int> &WWWTree::lep_tightId() {
  if (not lep_tightId_isLoaded) {
    if (lep_tightId_branch != 0) {
      lep_tightId_branch->GetEntry(index);
    } else {
      printf("branch lep_tightId_branch does not exist!\n");
      exit(1);
    }
    lep_tightId_isLoaded = true;
  }
  return *lep_tightId_;
}

const vector<int> &WWWTree::lep_mcMatchId() {
  if (not lep_mcMatchId_isLoaded) {
    if (lep_mcMatchId_branch != 0) {
      lep_mcMatchId_branch->GetEntry(index);
    } else {
      printf("branch lep_mcMatchId_branch does not exist!\n");
      exit(1);
    }
    lep_mcMatchId_isLoaded = true;
  }
  return *lep_mcMatchId_;
}

const vector<int> &WWWTree::lep_lostHits() {
  if (not lep_lostHits_isLoaded) {
    if (lep_lostHits_branch != 0) {
      lep_lostHits_branch->GetEntry(index);
    } else {
      printf("branch lep_lostHits_branch does not exist!\n");
      exit(1);
    }
    lep_lostHits_isLoaded = true;
  }
  return *lep_lostHits_;
}

const vector<int> &WWWTree::lep_convVeto() {
  if (not lep_convVeto_isLoaded) {
    if (lep_convVeto_branch != 0) {
      lep_convVeto_branch->GetEntry(index);
    } else {
      printf("branch lep_convVeto_branch does not exist!\n");
      exit(1);
    }
    lep_convVeto_isLoaded = true;
  }
  return *lep_convVeto_;
}

const vector<int> &WWWTree::lep_tightCharge() {
  if (not lep_tightCharge_isLoaded) {
    if (lep_tightCharge_branch != 0) {
      lep_tightCharge_branch->GetEntry(index);
    } else {
      printf("branch lep_tightCharge_branch does not exist!\n");
      exit(1);
    }
    lep_tightCharge_isLoaded = true;
  }
  return *lep_tightCharge_;
}

const vector<float> &WWWTree::lep_MVA() {
  if (not lep_MVA_isLoaded) {
    if (lep_MVA_branch != 0) {
      lep_MVA_branch->GetEntry(index);
    } else {
      printf("branch lep_MVA_branch does not exist!\n");
      exit(1);
    }
    lep_MVA_isLoaded = true;
  }
  return *lep_MVA_;
}

const vector<float> &WWWTree::lep_validfraction() {
  if (not lep_validfraction_isLoaded) {
    if (lep_validfraction_branch != 0) {
      lep_validfraction_branch->GetEntry(index);
    } else {
      printf("branch lep_validfraction_branch does not exist!\n");
      exit(1);
    }
    lep_validfraction_isLoaded = true;
  }
  return *lep_validfraction_;
}

const vector<float> &WWWTree::lep_pterr() {
  if (not lep_pterr_isLoaded) {
    if (lep_pterr_branch != 0) {
      lep_pterr_branch->GetEntry(index);
    } else {
      printf("branch lep_pterr_branch does not exist!\n");
      exit(1);
    }
    lep_pterr_isLoaded = true;
  }
  return *lep_pterr_;
}

const vector<float> &WWWTree::lep_sta_pterrOpt() {
  if (not lep_sta_pterrOpt_isLoaded) {
    if (lep_sta_pterrOpt_branch != 0) {
      lep_sta_pterrOpt_branch->GetEntry(index);
    } else {
      printf("branch lep_sta_pterrOpt_branch does not exist!\n");
      exit(1);
    }
    lep_sta_pterrOpt_isLoaded = true;
  }
  return *lep_sta_pterrOpt_;
}

const vector<float> &WWWTree::lep_glb_pterrOpt() {
  if (not lep_glb_pterrOpt_isLoaded) {
    if (lep_glb_pterrOpt_branch != 0) {
      lep_glb_pterrOpt_branch->GetEntry(index);
    } else {
      printf("branch lep_glb_pterrOpt_branch does not exist!\n");
      exit(1);
    }
    lep_glb_pterrOpt_isLoaded = true;
  }
  return *lep_glb_pterrOpt_;
}

const vector<float> &WWWTree::lep_x2ondof() {
  if (not lep_x2ondof_isLoaded) {
    if (lep_x2ondof_branch != 0) {
      lep_x2ondof_branch->GetEntry(index);
    } else {
      printf("branch lep_x2ondof_branch does not exist!\n");
      exit(1);
    }
    lep_x2ondof_isLoaded = true;
  }
  return *lep_x2ondof_;
}

const vector<float> &WWWTree::lep_sta_x2ondof() {
  if (not lep_sta_x2ondof_isLoaded) {
    if (lep_sta_x2ondof_branch != 0) {
      lep_sta_x2ondof_branch->GetEntry(index);
    } else {
      printf("branch lep_sta_x2ondof_branch does not exist!\n");
      exit(1);
    }
    lep_sta_x2ondof_isLoaded = true;
  }
  return *lep_sta_x2ondof_;
}

const vector<float> &WWWTree::lep_glb_x2ondof() {
  if (not lep_glb_x2ondof_isLoaded) {
    if (lep_glb_x2ondof_branch != 0) {
      lep_glb_x2ondof_branch->GetEntry(index);
    } else {
      printf("branch lep_glb_x2ondof_branch does not exist!\n");
      exit(1);
    }
    lep_glb_x2ondof_isLoaded = true;
  }
  return *lep_glb_x2ondof_;
}

const int &WWWTree::nisoTrack_5gev() {
  if (not nisoTrack_5gev_isLoaded) {
    if (nisoTrack_5gev_branch != 0) {
      nisoTrack_5gev_branch->GetEntry(index);
    } else {
      printf("branch nisoTrack_5gev_branch does not exist!\n");
      exit(1);
    }
    nisoTrack_5gev_isLoaded = true;
  }
  return nisoTrack_5gev_;
}

const int &WWWTree::nisoTrack_mt2() {
  if (not nisoTrack_mt2_isLoaded) {
    if (nisoTrack_mt2_branch != 0) {
      nisoTrack_mt2_branch->GetEntry(index);
    } else {
      printf("branch nisoTrack_mt2_branch does not exist!\n");
      exit(1);
    }
    nisoTrack_mt2_isLoaded = true;
  }
  return nisoTrack_mt2_;
}

const int &WWWTree::nisoTrack_mt2_cleaned_VVV_cutbased_veto() {
  if (not nisoTrack_mt2_cleaned_VVV_cutbased_veto_isLoaded) {
    if (nisoTrack_mt2_cleaned_VVV_cutbased_veto_branch != 0) {
      nisoTrack_mt2_cleaned_VVV_cutbased_veto_branch->GetEntry(index);
    } else {
      printf("branch nisoTrack_mt2_cleaned_VVV_cutbased_veto_branch does not exist!\n");
      exit(1);
    }
    nisoTrack_mt2_cleaned_VVV_cutbased_veto_isLoaded = true;
  }
  return nisoTrack_mt2_cleaned_VVV_cutbased_veto_;
}

const int &WWWTree::nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso() {
  if (not nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_isLoaded) {
    if (nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_branch != 0) {
      nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_branch->GetEntry(index);
    } else {
      printf("branch nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_branch does not exist!\n");
      exit(1);
    }
    nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_isLoaded = true;
  }
  return nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_;
}

const int &WWWTree::nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip() {
  if (not nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_isLoaded) {
    if (nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_branch != 0) {
      nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_branch->GetEntry(index);
    } else {
      printf("branch nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_branch does not exist!\n");
      exit(1);
    }
    nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_isLoaded = true;
  }
  return nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_;
}

const int &WWWTree::nisoTrack_mt2_cleaned_VVV_cutbased_veto_bak1() {
  if (not nisoTrack_mt2_cleaned_VVV_cutbased_veto_bak1_isLoaded) {
    if (nisoTrack_mt2_cleaned_VVV_cutbased_veto_bak1_branch != 0) {
      nisoTrack_mt2_cleaned_VVV_cutbased_veto_bak1_branch->GetEntry(index);
    } else {
      printf("branch nisoTrack_mt2_cleaned_VVV_cutbased_veto_bak1_branch does not exist!\n");
      exit(1);
    }
    nisoTrack_mt2_cleaned_VVV_cutbased_veto_bak1_isLoaded = true;
  }
  return nisoTrack_mt2_cleaned_VVV_cutbased_veto_bak1_;
}

const int &WWWTree::nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_bak1() {
  if (not nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_bak1_isLoaded) {
    if (nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_bak1_branch != 0) {
      nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_bak1_branch->GetEntry(index);
    } else {
      printf("branch nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_bak1_branch does not exist!\n");
      exit(1);
    }
    nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_bak1_isLoaded = true;
  }
  return nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_bak1_;
}

const int &WWWTree::nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_bak1() {
  if (not nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_bak1_isLoaded) {
    if (nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_bak1_branch != 0) {
      nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_bak1_branch->GetEntry(index);
    } else {
      printf("branch nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_bak1_branch does not exist!\n");
      exit(1);
    }
    nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_bak1_isLoaded = true;
  }
  return nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_bak1_;
}

const int &WWWTree::nisoTrack_mt2_cleaned_VVV_cutbased_fo() {
  if (not nisoTrack_mt2_cleaned_VVV_cutbased_fo_isLoaded) {
    if (nisoTrack_mt2_cleaned_VVV_cutbased_fo_branch != 0) {
      nisoTrack_mt2_cleaned_VVV_cutbased_fo_branch->GetEntry(index);
    } else {
      printf("branch nisoTrack_mt2_cleaned_VVV_cutbased_fo_branch does not exist!\n");
      exit(1);
    }
    nisoTrack_mt2_cleaned_VVV_cutbased_fo_isLoaded = true;
  }
  return nisoTrack_mt2_cleaned_VVV_cutbased_fo_;
}

const int &WWWTree::nisoTrack_mt2_cleaned_VVV_cutbased_fo_noiso() {
  if (not nisoTrack_mt2_cleaned_VVV_cutbased_fo_noiso_isLoaded) {
    if (nisoTrack_mt2_cleaned_VVV_cutbased_fo_noiso_branch != 0) {
      nisoTrack_mt2_cleaned_VVV_cutbased_fo_noiso_branch->GetEntry(index);
    } else {
      printf("branch nisoTrack_mt2_cleaned_VVV_cutbased_fo_noiso_branch does not exist!\n");
      exit(1);
    }
    nisoTrack_mt2_cleaned_VVV_cutbased_fo_noiso_isLoaded = true;
  }
  return nisoTrack_mt2_cleaned_VVV_cutbased_fo_noiso_;
}

const int &WWWTree::nisoTrack_mt2_cleaned_VVV_cutbased_tight_noiso() {
  if (not nisoTrack_mt2_cleaned_VVV_cutbased_tight_noiso_isLoaded) {
    if (nisoTrack_mt2_cleaned_VVV_cutbased_tight_noiso_branch != 0) {
      nisoTrack_mt2_cleaned_VVV_cutbased_tight_noiso_branch->GetEntry(index);
    } else {
      printf("branch nisoTrack_mt2_cleaned_VVV_cutbased_tight_noiso_branch does not exist!\n");
      exit(1);
    }
    nisoTrack_mt2_cleaned_VVV_cutbased_tight_noiso_isLoaded = true;
  }
  return nisoTrack_mt2_cleaned_VVV_cutbased_tight_noiso_;
}

const int &WWWTree::nisoTrack_mt2_cleaned_VVV_cutbased_tight() {
  if (not nisoTrack_mt2_cleaned_VVV_cutbased_tight_isLoaded) {
    if (nisoTrack_mt2_cleaned_VVV_cutbased_tight_branch != 0) {
      nisoTrack_mt2_cleaned_VVV_cutbased_tight_branch->GetEntry(index);
    } else {
      printf("branch nisoTrack_mt2_cleaned_VVV_cutbased_tight_branch does not exist!\n");
      exit(1);
    }
    nisoTrack_mt2_cleaned_VVV_cutbased_tight_isLoaded = true;
  }
  return nisoTrack_mt2_cleaned_VVV_cutbased_tight_;
}

const int &WWWTree::nisoTrack_mt2_cleaned_VVV_MVAbased_tight_noiso() {
  if (not nisoTrack_mt2_cleaned_VVV_MVAbased_tight_noiso_isLoaded) {
    if (nisoTrack_mt2_cleaned_VVV_MVAbased_tight_noiso_branch != 0) {
      nisoTrack_mt2_cleaned_VVV_MVAbased_tight_noiso_branch->GetEntry(index);
    } else {
      printf("branch nisoTrack_mt2_cleaned_VVV_MVAbased_tight_noiso_branch does not exist!\n");
      exit(1);
    }
    nisoTrack_mt2_cleaned_VVV_MVAbased_tight_noiso_isLoaded = true;
  }
  return nisoTrack_mt2_cleaned_VVV_MVAbased_tight_noiso_;
}

const int &WWWTree::nisoTrack_mt2_cleaned_VVV_MVAbased_tight() {
  if (not nisoTrack_mt2_cleaned_VVV_MVAbased_tight_isLoaded) {
    if (nisoTrack_mt2_cleaned_VVV_MVAbased_tight_branch != 0) {
      nisoTrack_mt2_cleaned_VVV_MVAbased_tight_branch->GetEntry(index);
    } else {
      printf("branch nisoTrack_mt2_cleaned_VVV_MVAbased_tight_branch does not exist!\n");
      exit(1);
    }
    nisoTrack_mt2_cleaned_VVV_MVAbased_tight_isLoaded = true;
  }
  return nisoTrack_mt2_cleaned_VVV_MVAbased_tight_;
}

const int &WWWTree::nisoTrack_mt2_cleaned_VVV_baseline() {
  if (not nisoTrack_mt2_cleaned_VVV_baseline_isLoaded) {
    if (nisoTrack_mt2_cleaned_VVV_baseline_branch != 0) {
      nisoTrack_mt2_cleaned_VVV_baseline_branch->GetEntry(index);
    } else {
      printf("branch nisoTrack_mt2_cleaned_VVV_baseline_branch does not exist!\n");
      exit(1);
    }
    nisoTrack_mt2_cleaned_VVV_baseline_isLoaded = true;
  }
  return nisoTrack_mt2_cleaned_VVV_baseline_;
}

const int &WWWTree::nisoTrack_stop() {
  if (not nisoTrack_stop_isLoaded) {
    if (nisoTrack_stop_branch != 0) {
      nisoTrack_stop_branch->GetEntry(index);
    } else {
      printf("branch nisoTrack_stop_branch does not exist!\n");
      exit(1);
    }
    nisoTrack_stop_isLoaded = true;
  }
  return nisoTrack_stop_;
}

const int &WWWTree::nisoTrack_PFLep5_woverlaps() {
  if (not nisoTrack_PFLep5_woverlaps_isLoaded) {
    if (nisoTrack_PFLep5_woverlaps_branch != 0) {
      nisoTrack_PFLep5_woverlaps_branch->GetEntry(index);
    } else {
      printf("branch nisoTrack_PFLep5_woverlaps_branch does not exist!\n");
      exit(1);
    }
    nisoTrack_PFLep5_woverlaps_isLoaded = true;
  }
  return nisoTrack_PFLep5_woverlaps_;
}

const int &WWWTree::nisoTrack_PFHad10_woverlaps() {
  if (not nisoTrack_PFHad10_woverlaps_isLoaded) {
    if (nisoTrack_PFHad10_woverlaps_branch != 0) {
      nisoTrack_PFHad10_woverlaps_branch->GetEntry(index);
    } else {
      printf("branch nisoTrack_PFHad10_woverlaps_branch does not exist!\n");
      exit(1);
    }
    nisoTrack_PFHad10_woverlaps_isLoaded = true;
  }
  return nisoTrack_PFHad10_woverlaps_;
}

const int &WWWTree::ngamma() {
  if (not ngamma_isLoaded) {
    if (ngamma_branch != 0) {
      ngamma_branch->GetEntry(index);
    } else {
      printf("branch ngamma_branch does not exist!\n");
      exit(1);
    }
    ngamma_isLoaded = true;
  }
  return ngamma_;
}

const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &WWWTree::gamma_p4() {
  if (not gamma_p4_isLoaded) {
    if (gamma_p4_branch != 0) {
      gamma_p4_branch->GetEntry(index);
    } else {
      printf("branch gamma_p4_branch does not exist!\n");
      exit(1);
    }
    gamma_p4_isLoaded = true;
  }
  return *gamma_p4_;
}

const vector<float> &WWWTree::gamma_pt() {
  if (not gamma_pt_isLoaded) {
    if (gamma_pt_branch != 0) {
      gamma_pt_branch->GetEntry(index);
    } else {
      printf("branch gamma_pt_branch does not exist!\n");
      exit(1);
    }
    gamma_pt_isLoaded = true;
  }
  return *gamma_pt_;
}

const vector<float> &WWWTree::gamma_eta() {
  if (not gamma_eta_isLoaded) {
    if (gamma_eta_branch != 0) {
      gamma_eta_branch->GetEntry(index);
    } else {
      printf("branch gamma_eta_branch does not exist!\n");
      exit(1);
    }
    gamma_eta_isLoaded = true;
  }
  return *gamma_eta_;
}

const vector<float> &WWWTree::gamma_phi() {
  if (not gamma_phi_isLoaded) {
    if (gamma_phi_branch != 0) {
      gamma_phi_branch->GetEntry(index);
    } else {
      printf("branch gamma_phi_branch does not exist!\n");
      exit(1);
    }
    gamma_phi_isLoaded = true;
  }
  return *gamma_phi_;
}

const vector<float> &WWWTree::gamma_mass() {
  if (not gamma_mass_isLoaded) {
    if (gamma_mass_branch != 0) {
      gamma_mass_branch->GetEntry(index);
    } else {
      printf("branch gamma_mass_branch does not exist!\n");
      exit(1);
    }
    gamma_mass_isLoaded = true;
  }
  return *gamma_mass_;
}

const vector<int> &WWWTree::gamma_mcMatchId() {
  if (not gamma_mcMatchId_isLoaded) {
    if (gamma_mcMatchId_branch != 0) {
      gamma_mcMatchId_branch->GetEntry(index);
    } else {
      printf("branch gamma_mcMatchId_branch does not exist!\n");
      exit(1);
    }
    gamma_mcMatchId_isLoaded = true;
  }
  return *gamma_mcMatchId_;
}

const vector<float> &WWWTree::gamma_genPt() {
  if (not gamma_genPt_isLoaded) {
    if (gamma_genPt_branch != 0) {
      gamma_genPt_branch->GetEntry(index);
    } else {
      printf("branch gamma_genPt_branch does not exist!\n");
      exit(1);
    }
    gamma_genPt_isLoaded = true;
  }
  return *gamma_genPt_;
}

const vector<float> &WWWTree::gamma_genIso() {
  if (not gamma_genIso_isLoaded) {
    if (gamma_genIso_branch != 0) {
      gamma_genIso_branch->GetEntry(index);
    } else {
      printf("branch gamma_genIso_branch does not exist!\n");
      exit(1);
    }
    gamma_genIso_isLoaded = true;
  }
  return *gamma_genIso_;
}

const vector<float> &WWWTree::gamma_chHadIso() {
  if (not gamma_chHadIso_isLoaded) {
    if (gamma_chHadIso_branch != 0) {
      gamma_chHadIso_branch->GetEntry(index);
    } else {
      printf("branch gamma_chHadIso_branch does not exist!\n");
      exit(1);
    }
    gamma_chHadIso_isLoaded = true;
  }
  return *gamma_chHadIso_;
}

const vector<float> &WWWTree::gamma_neuHadIso() {
  if (not gamma_neuHadIso_isLoaded) {
    if (gamma_neuHadIso_branch != 0) {
      gamma_neuHadIso_branch->GetEntry(index);
    } else {
      printf("branch gamma_neuHadIso_branch does not exist!\n");
      exit(1);
    }
    gamma_neuHadIso_isLoaded = true;
  }
  return *gamma_neuHadIso_;
}

const vector<float> &WWWTree::gamma_phIso() {
  if (not gamma_phIso_isLoaded) {
    if (gamma_phIso_branch != 0) {
      gamma_phIso_branch->GetEntry(index);
    } else {
      printf("branch gamma_phIso_branch does not exist!\n");
      exit(1);
    }
    gamma_phIso_isLoaded = true;
  }
  return *gamma_phIso_;
}

const vector<float> &WWWTree::gamma_sigmaIetaIeta() {
  if (not gamma_sigmaIetaIeta_isLoaded) {
    if (gamma_sigmaIetaIeta_branch != 0) {
      gamma_sigmaIetaIeta_branch->GetEntry(index);
    } else {
      printf("branch gamma_sigmaIetaIeta_branch does not exist!\n");
      exit(1);
    }
    gamma_sigmaIetaIeta_isLoaded = true;
  }
  return *gamma_sigmaIetaIeta_;
}

const vector<float> &WWWTree::gamma_r9() {
  if (not gamma_r9_isLoaded) {
    if (gamma_r9_branch != 0) {
      gamma_r9_branch->GetEntry(index);
    } else {
      printf("branch gamma_r9_branch does not exist!\n");
      exit(1);
    }
    gamma_r9_isLoaded = true;
  }
  return *gamma_r9_;
}

const vector<float> &WWWTree::gamma_hOverE() {
  if (not gamma_hOverE_isLoaded) {
    if (gamma_hOverE_branch != 0) {
      gamma_hOverE_branch->GetEntry(index);
    } else {
      printf("branch gamma_hOverE_branch does not exist!\n");
      exit(1);
    }
    gamma_hOverE_isLoaded = true;
  }
  return *gamma_hOverE_;
}

const vector<float> &WWWTree::gamma_hOverE_online() {
  if (not gamma_hOverE_online_isLoaded) {
    if (gamma_hOverE_online_branch != 0) {
      gamma_hOverE_online_branch->GetEntry(index);
    } else {
      printf("branch gamma_hOverE_online_branch does not exist!\n");
      exit(1);
    }
    gamma_hOverE_online_isLoaded = true;
  }
  return *gamma_hOverE_online_;
}

const vector<int> &WWWTree::gamma_idCutBased() {
  if (not gamma_idCutBased_isLoaded) {
    if (gamma_idCutBased_branch != 0) {
      gamma_idCutBased_branch->GetEntry(index);
    } else {
      printf("branch gamma_idCutBased_branch does not exist!\n");
      exit(1);
    }
    gamma_idCutBased_isLoaded = true;
  }
  return *gamma_idCutBased_;
}

const vector<float> &WWWTree::gamma_ecpfclusiso() {
  if (not gamma_ecpfclusiso_isLoaded) {
    if (gamma_ecpfclusiso_branch != 0) {
      gamma_ecpfclusiso_branch->GetEntry(index);
    } else {
      printf("branch gamma_ecpfclusiso_branch does not exist!\n");
      exit(1);
    }
    gamma_ecpfclusiso_isLoaded = true;
  }
  return *gamma_ecpfclusiso_;
}

const vector<float> &WWWTree::gamma_hcpfclusiso() {
  if (not gamma_hcpfclusiso_isLoaded) {
    if (gamma_hcpfclusiso_branch != 0) {
      gamma_hcpfclusiso_branch->GetEntry(index);
    } else {
      printf("branch gamma_hcpfclusiso_branch does not exist!\n");
      exit(1);
    }
    gamma_hcpfclusiso_isLoaded = true;
  }
  return *gamma_hcpfclusiso_;
}

const vector<float> &WWWTree::gamma_hollowtkiso03() {
  if (not gamma_hollowtkiso03_isLoaded) {
    if (gamma_hollowtkiso03_branch != 0) {
      gamma_hollowtkiso03_branch->GetEntry(index);
    } else {
      printf("branch gamma_hollowtkiso03_branch does not exist!\n");
      exit(1);
    }
    gamma_hollowtkiso03_isLoaded = true;
  }
  return *gamma_hollowtkiso03_;
}

const vector<int> &WWWTree::gamma_genIsPromptFinalState() {
  if (not gamma_genIsPromptFinalState_isLoaded) {
    if (gamma_genIsPromptFinalState_branch != 0) {
      gamma_genIsPromptFinalState_branch->GetEntry(index);
    } else {
      printf("branch gamma_genIsPromptFinalState_branch does not exist!\n");
      exit(1);
    }
    gamma_genIsPromptFinalState_isLoaded = true;
  }
  return *gamma_genIsPromptFinalState_;
}

const vector<float> &WWWTree::gamma_drMinParton() {
  if (not gamma_drMinParton_isLoaded) {
    if (gamma_drMinParton_branch != 0) {
      gamma_drMinParton_branch->GetEntry(index);
    } else {
      printf("branch gamma_drMinParton_branch does not exist!\n");
      exit(1);
    }
    gamma_drMinParton_isLoaded = true;
  }
  return *gamma_drMinParton_;
}

const int &WWWTree::ngenPart() {
  if (not ngenPart_isLoaded) {
    if (ngenPart_branch != 0) {
      ngenPart_branch->GetEntry(index);
    } else {
      printf("branch ngenPart_branch does not exist!\n");
      exit(1);
    }
    ngenPart_isLoaded = true;
  }
  return ngenPart_;
}

const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &WWWTree::genPart_p4() {
  if (not genPart_p4_isLoaded) {
    if (genPart_p4_branch != 0) {
      genPart_p4_branch->GetEntry(index);
    } else {
      printf("branch genPart_p4_branch does not exist!\n");
      exit(1);
    }
    genPart_p4_isLoaded = true;
  }
  return *genPart_p4_;
}

const vector<float> &WWWTree::genPart_pt() {
  if (not genPart_pt_isLoaded) {
    if (genPart_pt_branch != 0) {
      genPart_pt_branch->GetEntry(index);
    } else {
      printf("branch genPart_pt_branch does not exist!\n");
      exit(1);
    }
    genPart_pt_isLoaded = true;
  }
  return *genPart_pt_;
}

const vector<float> &WWWTree::genPart_eta() {
  if (not genPart_eta_isLoaded) {
    if (genPart_eta_branch != 0) {
      genPart_eta_branch->GetEntry(index);
    } else {
      printf("branch genPart_eta_branch does not exist!\n");
      exit(1);
    }
    genPart_eta_isLoaded = true;
  }
  return *genPart_eta_;
}

const vector<float> &WWWTree::genPart_phi() {
  if (not genPart_phi_isLoaded) {
    if (genPart_phi_branch != 0) {
      genPart_phi_branch->GetEntry(index);
    } else {
      printf("branch genPart_phi_branch does not exist!\n");
      exit(1);
    }
    genPart_phi_isLoaded = true;
  }
  return *genPart_phi_;
}

const vector<float> &WWWTree::genPart_mass() {
  if (not genPart_mass_isLoaded) {
    if (genPart_mass_branch != 0) {
      genPart_mass_branch->GetEntry(index);
    } else {
      printf("branch genPart_mass_branch does not exist!\n");
      exit(1);
    }
    genPart_mass_isLoaded = true;
  }
  return *genPart_mass_;
}

const vector<int> &WWWTree::genPart_pdgId() {
  if (not genPart_pdgId_isLoaded) {
    if (genPart_pdgId_branch != 0) {
      genPart_pdgId_branch->GetEntry(index);
    } else {
      printf("branch genPart_pdgId_branch does not exist!\n");
      exit(1);
    }
    genPart_pdgId_isLoaded = true;
  }
  return *genPart_pdgId_;
}

const vector<int> &WWWTree::genPart_status() {
  if (not genPart_status_isLoaded) {
    if (genPart_status_branch != 0) {
      genPart_status_branch->GetEntry(index);
    } else {
      printf("branch genPart_status_branch does not exist!\n");
      exit(1);
    }
    genPart_status_isLoaded = true;
  }
  return *genPart_status_;
}

const vector<float> &WWWTree::genPart_charge() {
  if (not genPart_charge_isLoaded) {
    if (genPart_charge_branch != 0) {
      genPart_charge_branch->GetEntry(index);
    } else {
      printf("branch genPart_charge_branch does not exist!\n");
      exit(1);
    }
    genPart_charge_isLoaded = true;
  }
  return *genPart_charge_;
}

const vector<int> &WWWTree::genPart_motherId() {
  if (not genPart_motherId_isLoaded) {
    if (genPart_motherId_branch != 0) {
      genPart_motherId_branch->GetEntry(index);
    } else {
      printf("branch genPart_motherId_branch does not exist!\n");
      exit(1);
    }
    genPart_motherId_isLoaded = true;
  }
  return *genPart_motherId_;
}

const vector<int> &WWWTree::genPart_grandmaId() {
  if (not genPart_grandmaId_isLoaded) {
    if (genPart_grandmaId_branch != 0) {
      genPart_grandmaId_branch->GetEntry(index);
    } else {
      printf("branch genPart_grandmaId_branch does not exist!\n");
      exit(1);
    }
    genPart_grandmaId_isLoaded = true;
  }
  return *genPart_grandmaId_;
}

const vector<bool> &WWWTree::genPart_isp6status3() {
  if (not genPart_isp6status3_isLoaded) {
    if (genPart_isp6status3_branch != 0) {
      genPart_isp6status3_branch->GetEntry(index);
    } else {
      printf("branch genPart_isp6status3_branch does not exist!\n");
      exit(1);
    }
    genPart_isp6status3_isLoaded = true;
  }
  return *genPart_isp6status3_;
}

const int &WWWTree::ngenLep() {
  if (not ngenLep_isLoaded) {
    if (ngenLep_branch != 0) {
      ngenLep_branch->GetEntry(index);
    } else {
      printf("branch ngenLep_branch does not exist!\n");
      exit(1);
    }
    ngenLep_isLoaded = true;
  }
  return ngenLep_;
}

const vector<float> &WWWTree::genLep_pt() {
  if (not genLep_pt_isLoaded) {
    if (genLep_pt_branch != 0) {
      genLep_pt_branch->GetEntry(index);
    } else {
      printf("branch genLep_pt_branch does not exist!\n");
      exit(1);
    }
    genLep_pt_isLoaded = true;
  }
  return *genLep_pt_;
}

const vector<float> &WWWTree::genLep_eta() {
  if (not genLep_eta_isLoaded) {
    if (genLep_eta_branch != 0) {
      genLep_eta_branch->GetEntry(index);
    } else {
      printf("branch genLep_eta_branch does not exist!\n");
      exit(1);
    }
    genLep_eta_isLoaded = true;
  }
  return *genLep_eta_;
}

const vector<float> &WWWTree::genLep_phi() {
  if (not genLep_phi_isLoaded) {
    if (genLep_phi_branch != 0) {
      genLep_phi_branch->GetEntry(index);
    } else {
      printf("branch genLep_phi_branch does not exist!\n");
      exit(1);
    }
    genLep_phi_isLoaded = true;
  }
  return *genLep_phi_;
}

const vector<float> &WWWTree::genLep_mass() {
  if (not genLep_mass_isLoaded) {
    if (genLep_mass_branch != 0) {
      genLep_mass_branch->GetEntry(index);
    } else {
      printf("branch genLep_mass_branch does not exist!\n");
      exit(1);
    }
    genLep_mass_isLoaded = true;
  }
  return *genLep_mass_;
}

const vector<int> &WWWTree::genLep_pdgId() {
  if (not genLep_pdgId_isLoaded) {
    if (genLep_pdgId_branch != 0) {
      genLep_pdgId_branch->GetEntry(index);
    } else {
      printf("branch genLep_pdgId_branch does not exist!\n");
      exit(1);
    }
    genLep_pdgId_isLoaded = true;
  }
  return *genLep_pdgId_;
}

const vector<int> &WWWTree::genLep_status() {
  if (not genLep_status_isLoaded) {
    if (genLep_status_branch != 0) {
      genLep_status_branch->GetEntry(index);
    } else {
      printf("branch genLep_status_branch does not exist!\n");
      exit(1);
    }
    genLep_status_isLoaded = true;
  }
  return *genLep_status_;
}

const vector<float> &WWWTree::genLep_charge() {
  if (not genLep_charge_isLoaded) {
    if (genLep_charge_branch != 0) {
      genLep_charge_branch->GetEntry(index);
    } else {
      printf("branch genLep_charge_branch does not exist!\n");
      exit(1);
    }
    genLep_charge_isLoaded = true;
  }
  return *genLep_charge_;
}

const vector<int> &WWWTree::genLep_sourceId() {
  if (not genLep_sourceId_isLoaded) {
    if (genLep_sourceId_branch != 0) {
      genLep_sourceId_branch->GetEntry(index);
    } else {
      printf("branch genLep_sourceId_branch does not exist!\n");
      exit(1);
    }
    genLep_sourceId_isLoaded = true;
  }
  return *genLep_sourceId_;
}

const vector<bool> &WWWTree::genLep_isp6status3() {
  if (not genLep_isp6status3_isLoaded) {
    if (genLep_isp6status3_branch != 0) {
      genLep_isp6status3_branch->GetEntry(index);
    } else {
      printf("branch genLep_isp6status3_branch does not exist!\n");
      exit(1);
    }
    genLep_isp6status3_isLoaded = true;
  }
  return *genLep_isp6status3_;
}

const int &WWWTree::ngenTau() {
  if (not ngenTau_isLoaded) {
    if (ngenTau_branch != 0) {
      ngenTau_branch->GetEntry(index);
    } else {
      printf("branch ngenTau_branch does not exist!\n");
      exit(1);
    }
    ngenTau_isLoaded = true;
  }
  return ngenTau_;
}

const vector<float> &WWWTree::genTau_pt() {
  if (not genTau_pt_isLoaded) {
    if (genTau_pt_branch != 0) {
      genTau_pt_branch->GetEntry(index);
    } else {
      printf("branch genTau_pt_branch does not exist!\n");
      exit(1);
    }
    genTau_pt_isLoaded = true;
  }
  return *genTau_pt_;
}

const vector<float> &WWWTree::genTau_eta() {
  if (not genTau_eta_isLoaded) {
    if (genTau_eta_branch != 0) {
      genTau_eta_branch->GetEntry(index);
    } else {
      printf("branch genTau_eta_branch does not exist!\n");
      exit(1);
    }
    genTau_eta_isLoaded = true;
  }
  return *genTau_eta_;
}

const vector<float> &WWWTree::genTau_phi() {
  if (not genTau_phi_isLoaded) {
    if (genTau_phi_branch != 0) {
      genTau_phi_branch->GetEntry(index);
    } else {
      printf("branch genTau_phi_branch does not exist!\n");
      exit(1);
    }
    genTau_phi_isLoaded = true;
  }
  return *genTau_phi_;
}

const vector<float> &WWWTree::genTau_mass() {
  if (not genTau_mass_isLoaded) {
    if (genTau_mass_branch != 0) {
      genTau_mass_branch->GetEntry(index);
    } else {
      printf("branch genTau_mass_branch does not exist!\n");
      exit(1);
    }
    genTau_mass_isLoaded = true;
  }
  return *genTau_mass_;
}

const vector<int> &WWWTree::genTau_pdgId() {
  if (not genTau_pdgId_isLoaded) {
    if (genTau_pdgId_branch != 0) {
      genTau_pdgId_branch->GetEntry(index);
    } else {
      printf("branch genTau_pdgId_branch does not exist!\n");
      exit(1);
    }
    genTau_pdgId_isLoaded = true;
  }
  return *genTau_pdgId_;
}

const vector<int> &WWWTree::genTau_status() {
  if (not genTau_status_isLoaded) {
    if (genTau_status_branch != 0) {
      genTau_status_branch->GetEntry(index);
    } else {
      printf("branch genTau_status_branch does not exist!\n");
      exit(1);
    }
    genTau_status_isLoaded = true;
  }
  return *genTau_status_;
}

const vector<float> &WWWTree::genTau_charge() {
  if (not genTau_charge_isLoaded) {
    if (genTau_charge_branch != 0) {
      genTau_charge_branch->GetEntry(index);
    } else {
      printf("branch genTau_charge_branch does not exist!\n");
      exit(1);
    }
    genTau_charge_isLoaded = true;
  }
  return *genTau_charge_;
}

const vector<int> &WWWTree::genTau_sourceId() {
  if (not genTau_sourceId_isLoaded) {
    if (genTau_sourceId_branch != 0) {
      genTau_sourceId_branch->GetEntry(index);
    } else {
      printf("branch genTau_sourceId_branch does not exist!\n");
      exit(1);
    }
    genTau_sourceId_isLoaded = true;
  }
  return *genTau_sourceId_;
}

const int &WWWTree::ngenLepFromTau() {
  if (not ngenLepFromTau_isLoaded) {
    if (ngenLepFromTau_branch != 0) {
      ngenLepFromTau_branch->GetEntry(index);
    } else {
      printf("branch ngenLepFromTau_branch does not exist!\n");
      exit(1);
    }
    ngenLepFromTau_isLoaded = true;
  }
  return ngenLepFromTau_;
}

const vector<float> &WWWTree::genLepFromTau_pt() {
  if (not genLepFromTau_pt_isLoaded) {
    if (genLepFromTau_pt_branch != 0) {
      genLepFromTau_pt_branch->GetEntry(index);
    } else {
      printf("branch genLepFromTau_pt_branch does not exist!\n");
      exit(1);
    }
    genLepFromTau_pt_isLoaded = true;
  }
  return *genLepFromTau_pt_;
}

const vector<float> &WWWTree::genLepFromTau_eta() {
  if (not genLepFromTau_eta_isLoaded) {
    if (genLepFromTau_eta_branch != 0) {
      genLepFromTau_eta_branch->GetEntry(index);
    } else {
      printf("branch genLepFromTau_eta_branch does not exist!\n");
      exit(1);
    }
    genLepFromTau_eta_isLoaded = true;
  }
  return *genLepFromTau_eta_;
}

const vector<float> &WWWTree::genLepFromTau_phi() {
  if (not genLepFromTau_phi_isLoaded) {
    if (genLepFromTau_phi_branch != 0) {
      genLepFromTau_phi_branch->GetEntry(index);
    } else {
      printf("branch genLepFromTau_phi_branch does not exist!\n");
      exit(1);
    }
    genLepFromTau_phi_isLoaded = true;
  }
  return *genLepFromTau_phi_;
}

const vector<float> &WWWTree::genLepFromTau_mass() {
  if (not genLepFromTau_mass_isLoaded) {
    if (genLepFromTau_mass_branch != 0) {
      genLepFromTau_mass_branch->GetEntry(index);
    } else {
      printf("branch genLepFromTau_mass_branch does not exist!\n");
      exit(1);
    }
    genLepFromTau_mass_isLoaded = true;
  }
  return *genLepFromTau_mass_;
}

const vector<int> &WWWTree::genLepFromTau_pdgId() {
  if (not genLepFromTau_pdgId_isLoaded) {
    if (genLepFromTau_pdgId_branch != 0) {
      genLepFromTau_pdgId_branch->GetEntry(index);
    } else {
      printf("branch genLepFromTau_pdgId_branch does not exist!\n");
      exit(1);
    }
    genLepFromTau_pdgId_isLoaded = true;
  }
  return *genLepFromTau_pdgId_;
}

const vector<int> &WWWTree::genLepFromTau_status() {
  if (not genLepFromTau_status_isLoaded) {
    if (genLepFromTau_status_branch != 0) {
      genLepFromTau_status_branch->GetEntry(index);
    } else {
      printf("branch genLepFromTau_status_branch does not exist!\n");
      exit(1);
    }
    genLepFromTau_status_isLoaded = true;
  }
  return *genLepFromTau_status_;
}

const vector<float> &WWWTree::genLepFromTau_charge() {
  if (not genLepFromTau_charge_isLoaded) {
    if (genLepFromTau_charge_branch != 0) {
      genLepFromTau_charge_branch->GetEntry(index);
    } else {
      printf("branch genLepFromTau_charge_branch does not exist!\n");
      exit(1);
    }
    genLepFromTau_charge_isLoaded = true;
  }
  return *genLepFromTau_charge_;
}

const vector<int> &WWWTree::genLepFromTau_sourceId() {
  if (not genLepFromTau_sourceId_isLoaded) {
    if (genLepFromTau_sourceId_branch != 0) {
      genLepFromTau_sourceId_branch->GetEntry(index);
    } else {
      printf("branch genLepFromTau_sourceId_branch does not exist!\n");
      exit(1);
    }
    genLepFromTau_sourceId_isLoaded = true;
  }
  return *genLepFromTau_sourceId_;
}

const int &WWWTree::njets() {
  if (not njets_isLoaded) {
    if (njets_branch != 0) {
      njets_branch->GetEntry(index);
    } else {
      printf("branch njets_branch does not exist!\n");
      exit(1);
    }
    njets_isLoaded = true;
  }
  return njets_;
}

const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &WWWTree::jets_p4() {
  if (not jets_p4_isLoaded) {
    if (jets_p4_branch != 0) {
      jets_p4_branch->GetEntry(index);
    } else {
      printf("branch jets_p4_branch does not exist!\n");
      exit(1);
    }
    jets_p4_isLoaded = true;
  }
  return *jets_p4_;
}

const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &WWWTree::removed_jets_p4() {
  if (not removed_jets_p4_isLoaded) {
    if (removed_jets_p4_branch != 0) {
      removed_jets_p4_branch->GetEntry(index);
    } else {
      printf("branch removed_jets_p4_branch does not exist!\n");
      exit(1);
    }
    removed_jets_p4_isLoaded = true;
  }
  return *removed_jets_p4_;
}

const vector<float> &WWWTree::removed_jets_csv() {
  if (not removed_jets_csv_isLoaded) {
    if (removed_jets_csv_branch != 0) {
      removed_jets_csv_branch->GetEntry(index);
    } else {
      printf("branch removed_jets_csv_branch does not exist!\n");
      exit(1);
    }
    removed_jets_csv_isLoaded = true;
  }
  return *removed_jets_csv_;
}

const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &WWWTree::jets_medb_p4() {
  if (not jets_medb_p4_isLoaded) {
    if (jets_medb_p4_branch != 0) {
      jets_medb_p4_branch->GetEntry(index);
    } else {
      printf("branch jets_medb_p4_branch does not exist!\n");
      exit(1);
    }
    jets_medb_p4_isLoaded = true;
  }
  return *jets_medb_p4_;
}

const int &WWWTree::njets_up() {
  if (not njets_up_isLoaded) {
    if (njets_up_branch != 0) {
      njets_up_branch->GetEntry(index);
    } else {
      printf("branch njets_up_branch does not exist!\n");
      exit(1);
    }
    njets_up_isLoaded = true;
  }
  return njets_up_;
}

const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &WWWTree::jets_up_p4() {
  if (not jets_up_p4_isLoaded) {
    if (jets_up_p4_branch != 0) {
      jets_up_p4_branch->GetEntry(index);
    } else {
      printf("branch jets_up_p4_branch does not exist!\n");
      exit(1);
    }
    jets_up_p4_isLoaded = true;
  }
  return *jets_up_p4_;
}

const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &WWWTree::jets_medb_up_p4() {
  if (not jets_medb_up_p4_isLoaded) {
    if (jets_medb_up_p4_branch != 0) {
      jets_medb_up_p4_branch->GetEntry(index);
    } else {
      printf("branch jets_medb_up_p4_branch does not exist!\n");
      exit(1);
    }
    jets_medb_up_p4_isLoaded = true;
  }
  return *jets_medb_up_p4_;
}

const vector<float> &WWWTree::jets_csv() {
  if (not jets_csv_isLoaded) {
    if (jets_csv_branch != 0) {
      jets_csv_branch->GetEntry(index);
    } else {
      printf("branch jets_csv_branch does not exist!\n");
      exit(1);
    }
    jets_csv_isLoaded = true;
  }
  return *jets_csv_;
}

const vector<float> &WWWTree::jets_up_csv() {
  if (not jets_up_csv_isLoaded) {
    if (jets_up_csv_branch != 0) {
      jets_up_csv_branch->GetEntry(index);
    } else {
      printf("branch jets_up_csv_branch does not exist!\n");
      exit(1);
    }
    jets_up_csv_isLoaded = true;
  }
  return *jets_up_csv_;
}

const int &WWWTree::njets_dn() {
  if (not njets_dn_isLoaded) {
    if (njets_dn_branch != 0) {
      njets_dn_branch->GetEntry(index);
    } else {
      printf("branch njets_dn_branch does not exist!\n");
      exit(1);
    }
    njets_dn_isLoaded = true;
  }
  return njets_dn_;
}

const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &WWWTree::jets_dn_p4() {
  if (not jets_dn_p4_isLoaded) {
    if (jets_dn_p4_branch != 0) {
      jets_dn_p4_branch->GetEntry(index);
    } else {
      printf("branch jets_dn_p4_branch does not exist!\n");
      exit(1);
    }
    jets_dn_p4_isLoaded = true;
  }
  return *jets_dn_p4_;
}

const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &WWWTree::jets_medb_dn_p4() {
  if (not jets_medb_dn_p4_isLoaded) {
    if (jets_medb_dn_p4_branch != 0) {
      jets_medb_dn_p4_branch->GetEntry(index);
    } else {
      printf("branch jets_medb_dn_p4_branch does not exist!\n");
      exit(1);
    }
    jets_medb_dn_p4_isLoaded = true;
  }
  return *jets_medb_dn_p4_;
}

const vector<float> &WWWTree::jets_dn_csv() {
  if (not jets_dn_csv_isLoaded) {
    if (jets_dn_csv_branch != 0) {
      jets_dn_csv_branch->GetEntry(index);
    } else {
      printf("branch jets_dn_csv_branch does not exist!\n");
      exit(1);
    }
    jets_dn_csv_isLoaded = true;
  }
  return *jets_dn_csv_;
}

const vector<float> &WWWTree::jets_muf() {
  if (not jets_muf_isLoaded) {
    if (jets_muf_branch != 0) {
      jets_muf_branch->GetEntry(index);
    } else {
      printf("branch jets_muf_branch does not exist!\n");
      exit(1);
    }
    jets_muf_isLoaded = true;
  }
  return *jets_muf_;
}

const vector<int> &WWWTree::jets_mcFlavour() {
  if (not jets_mcFlavour_isLoaded) {
    if (jets_mcFlavour_branch != 0) {
      jets_mcFlavour_branch->GetEntry(index);
    } else {
      printf("branch jets_mcFlavour_branch does not exist!\n");
      exit(1);
    }
    jets_mcFlavour_isLoaded = true;
  }
  return *jets_mcFlavour_;
}

const vector<int> &WWWTree::jets_mcHadronFlav() {
  if (not jets_mcHadronFlav_isLoaded) {
    if (jets_mcHadronFlav_branch != 0) {
      jets_mcHadronFlav_branch->GetEntry(index);
    } else {
      printf("branch jets_mcHadronFlav_branch does not exist!\n");
      exit(1);
    }
    jets_mcHadronFlav_isLoaded = true;
  }
  return *jets_mcHadronFlav_;
}

const int &WWWTree::nhighPtPFcands() {
  if (not nhighPtPFcands_isLoaded) {
    if (nhighPtPFcands_branch != 0) {
      nhighPtPFcands_branch->GetEntry(index);
    } else {
      printf("branch nhighPtPFcands_branch does not exist!\n");
      exit(1);
    }
    nhighPtPFcands_isLoaded = true;
  }
  return nhighPtPFcands_;
}

const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &WWWTree::highPtPFcands_p4() {
  if (not highPtPFcands_p4_isLoaded) {
    if (highPtPFcands_p4_branch != 0) {
      highPtPFcands_p4_branch->GetEntry(index);
    } else {
      printf("branch highPtPFcands_p4_branch does not exist!\n");
      exit(1);
    }
    highPtPFcands_p4_isLoaded = true;
  }
  return *highPtPFcands_p4_;
}

const vector<float> &WWWTree::highPtPFcands_dz() {
  if (not highPtPFcands_dz_isLoaded) {
    if (highPtPFcands_dz_branch != 0) {
      highPtPFcands_dz_branch->GetEntry(index);
    } else {
      printf("branch highPtPFcands_dz_branch does not exist!\n");
      exit(1);
    }
    highPtPFcands_dz_isLoaded = true;
  }
  return *highPtPFcands_dz_;
}

const vector<int> &WWWTree::highPtPFcands_pdgId() {
  if (not highPtPFcands_pdgId_isLoaded) {
    if (highPtPFcands_pdgId_branch != 0) {
      highPtPFcands_pdgId_branch->GetEntry(index);
    } else {
      printf("branch highPtPFcands_pdgId_branch does not exist!\n");
      exit(1);
    }
    highPtPFcands_pdgId_isLoaded = true;
  }
  return *highPtPFcands_pdgId_;
}

const float &WWWTree::ht() {
  if (not ht_isLoaded) {
    if (ht_branch != 0) {
      ht_branch->GetEntry(index);
    } else {
      printf("branch ht_branch does not exist!\n");
      exit(1);
    }
    ht_isLoaded = true;
  }
  return ht_;
}

const float &WWWTree::ht_up() {
  if (not ht_up_isLoaded) {
    if (ht_up_branch != 0) {
      ht_up_branch->GetEntry(index);
    } else {
      printf("branch ht_up_branch does not exist!\n");
      exit(1);
    }
    ht_up_isLoaded = true;
  }
  return ht_up_;
}

const float &WWWTree::ht_dn() {
  if (not ht_dn_isLoaded) {
    if (ht_dn_branch != 0) {
      ht_dn_branch->GetEntry(index);
    } else {
      printf("branch ht_dn_branch does not exist!\n");
      exit(1);
    }
    ht_dn_isLoaded = true;
  }
  return ht_dn_;
}

const float &WWWTree::metsig_unofficial() {
  if (not metsig_unofficial_isLoaded) {
    if (metsig_unofficial_branch != 0) {
      metsig_unofficial_branch->GetEntry(index);
    } else {
      printf("branch metsig_unofficial_branch does not exist!\n");
      exit(1);
    }
    metsig_unofficial_isLoaded = true;
  }
  return metsig_unofficial_;
}

const float &WWWTree::metsig_unofficial_up() {
  if (not metsig_unofficial_up_isLoaded) {
    if (metsig_unofficial_up_branch != 0) {
      metsig_unofficial_up_branch->GetEntry(index);
    } else {
      printf("branch metsig_unofficial_up_branch does not exist!\n");
      exit(1);
    }
    metsig_unofficial_up_isLoaded = true;
  }
  return metsig_unofficial_up_;
}

const float &WWWTree::metsig_unofficial_dn() {
  if (not metsig_unofficial_dn_isLoaded) {
    if (metsig_unofficial_dn_branch != 0) {
      metsig_unofficial_dn_branch->GetEntry(index);
    } else {
      printf("branch metsig_unofficial_dn_branch does not exist!\n");
      exit(1);
    }
    metsig_unofficial_dn_isLoaded = true;
  }
  return metsig_unofficial_dn_;
}

const float &WWWTree::mt_lep1() {
  if (not mt_lep1_isLoaded) {
    if (mt_lep1_branch != 0) {
      mt_lep1_branch->GetEntry(index);
    } else {
      printf("branch mt_lep1_branch does not exist!\n");
      exit(1);
    }
    mt_lep1_isLoaded = true;
  }
  return mt_lep1_;
}

const float &WWWTree::mt2() {
  if (not mt2_isLoaded) {
    if (mt2_branch != 0) {
      mt2_branch->GetEntry(index);
    } else {
      printf("branch mt2_branch does not exist!\n");
      exit(1);
    }
    mt2_isLoaded = true;
  }
  return mt2_;
}

const float &WWWTree::mt2_up() {
  if (not mt2_up_isLoaded) {
    if (mt2_up_branch != 0) {
      mt2_up_branch->GetEntry(index);
    } else {
      printf("branch mt2_up_branch does not exist!\n");
      exit(1);
    }
    mt2_up_isLoaded = true;
  }
  return mt2_up_;
}

const float &WWWTree::mt2_dn() {
  if (not mt2_dn_isLoaded) {
    if (mt2_dn_branch != 0) {
      mt2_dn_branch->GetEntry(index);
    } else {
      printf("branch mt2_dn_branch does not exist!\n");
      exit(1);
    }
    mt2_dn_isLoaded = true;
  }
  return mt2_dn_;
}

const float &WWWTree::mt2j() {
  if (not mt2j_isLoaded) {
    if (mt2j_branch != 0) {
      mt2j_branch->GetEntry(index);
    } else {
      printf("branch mt2j_branch does not exist!\n");
      exit(1);
    }
    mt2j_isLoaded = true;
  }
  return mt2j_;
}

const float &WWWTree::mt2b() {
  if (not mt2b_isLoaded) {
    if (mt2b_branch != 0) {
      mt2b_branch->GetEntry(index);
    } else {
      printf("branch mt2b_branch does not exist!\n");
      exit(1);
    }
    mt2b_isLoaded = true;
  }
  return mt2b_;
}

const float &WWWTree::mt2b_up() {
  if (not mt2b_up_isLoaded) {
    if (mt2b_up_branch != 0) {
      mt2b_up_branch->GetEntry(index);
    } else {
      printf("branch mt2b_up_branch does not exist!\n");
      exit(1);
    }
    mt2b_up_isLoaded = true;
  }
  return mt2b_up_;
}

const float &WWWTree::mt2b_dn() {
  if (not mt2b_dn_isLoaded) {
    if (mt2b_dn_branch != 0) {
      mt2b_dn_branch->GetEntry(index);
    } else {
      printf("branch mt2b_dn_branch does not exist!\n");
      exit(1);
    }
    mt2b_dn_isLoaded = true;
  }
  return mt2b_dn_;
}

const float &WWWTree::mt2_genmet() {
  if (not mt2_genmet_isLoaded) {
    if (mt2_genmet_branch != 0) {
      mt2_genmet_branch->GetEntry(index);
    } else {
      printf("branch mt2_genmet_branch does not exist!\n");
      exit(1);
    }
    mt2_genmet_isLoaded = true;
  }
  return mt2_genmet_;
}

const float &WWWTree::mt2b_genmet() {
  if (not mt2b_genmet_isLoaded) {
    if (mt2b_genmet_branch != 0) {
      mt2b_genmet_branch->GetEntry(index);
    } else {
      printf("branch mt2b_genmet_branch does not exist!\n");
      exit(1);
    }
    mt2b_genmet_isLoaded = true;
  }
  return mt2b_genmet_;
}

const float &WWWTree::mjj_mindphi() {
  if (not mjj_mindphi_isLoaded) {
    if (mjj_mindphi_branch != 0) {
      mjj_mindphi_branch->GetEntry(index);
    } else {
      printf("branch mjj_mindphi_branch does not exist!\n");
      exit(1);
    }
    mjj_mindphi_isLoaded = true;
  }
  return mjj_mindphi_;
}

const float &WWWTree::mjj() {
  if (not mjj_isLoaded) {
    if (mjj_branch != 0) {
      mjj_branch->GetEntry(index);
    } else {
      printf("branch mjj_branch does not exist!\n");
      exit(1);
    }
    mjj_isLoaded = true;
  }
  return mjj_;
}

const float &WWWTree::mbb_csv() {
  if (not mbb_csv_isLoaded) {
    if (mbb_csv_branch != 0) {
      mbb_csv_branch->GetEntry(index);
    } else {
      printf("branch mbb_csv_branch does not exist!\n");
      exit(1);
    }
    mbb_csv_isLoaded = true;
  }
  return mbb_csv_;
}

const float &WWWTree::mbb_bpt() {
  if (not mbb_bpt_isLoaded) {
    if (mbb_bpt_branch != 0) {
      mbb_bpt_branch->GetEntry(index);
    } else {
      printf("branch mbb_bpt_branch does not exist!\n");
      exit(1);
    }
    mbb_bpt_isLoaded = true;
  }
  return mbb_bpt_;
}

const float &WWWTree::dphi_jj() {
  if (not dphi_jj_isLoaded) {
    if (dphi_jj_branch != 0) {
      dphi_jj_branch->GetEntry(index);
    } else {
      printf("branch dphi_jj_branch does not exist!\n");
      exit(1);
    }
    dphi_jj_isLoaded = true;
  }
  return dphi_jj_;
}

const float &WWWTree::dphi_ll() {
  if (not dphi_ll_isLoaded) {
    if (dphi_ll_branch != 0) {
      dphi_ll_branch->GetEntry(index);
    } else {
      printf("branch dphi_ll_branch does not exist!\n");
      exit(1);
    }
    dphi_ll_isLoaded = true;
  }
  return dphi_ll_;
}

const float &WWWTree::sum_mlb() {
  if (not sum_mlb_isLoaded) {
    if (sum_mlb_branch != 0) {
      sum_mlb_branch->GetEntry(index);
    } else {
      printf("branch sum_mlb_branch does not exist!\n");
      exit(1);
    }
    sum_mlb_isLoaded = true;
  }
  return sum_mlb_;
}

const float &WWWTree::deta_jj() {
  if (not deta_jj_isLoaded) {
    if (deta_jj_branch != 0) {
      deta_jj_branch->GetEntry(index);
    } else {
      printf("branch deta_jj_branch does not exist!\n");
      exit(1);
    }
    deta_jj_isLoaded = true;
  }
  return deta_jj_;
}

const float &WWWTree::dR_jj() {
  if (not dR_jj_isLoaded) {
    if (dR_jj_branch != 0) {
      dR_jj_branch->GetEntry(index);
    } else {
      printf("branch dR_jj_branch does not exist!\n");
      exit(1);
    }
    dR_jj_isLoaded = true;
  }
  return dR_jj_;
}

const float &WWWTree::dphi_metj1() {
  if (not dphi_metj1_isLoaded) {
    if (dphi_metj1_branch != 0) {
      dphi_metj1_branch->GetEntry(index);
    } else {
      printf("branch dphi_metj1_branch does not exist!\n");
      exit(1);
    }
    dphi_metj1_isLoaded = true;
  }
  return dphi_metj1_;
}

const float &WWWTree::dphi_metj2() {
  if (not dphi_metj2_isLoaded) {
    if (dphi_metj2_branch != 0) {
      dphi_metj2_branch->GetEntry(index);
    } else {
      printf("branch dphi_metj2_branch does not exist!\n");
      exit(1);
    }
    dphi_metj2_isLoaded = true;
  }
  return dphi_metj2_;
}

const float &WWWTree::dphi_genmetj1() {
  if (not dphi_genmetj1_isLoaded) {
    if (dphi_genmetj1_branch != 0) {
      dphi_genmetj1_branch->GetEntry(index);
    } else {
      printf("branch dphi_genmetj1_branch does not exist!\n");
      exit(1);
    }
    dphi_genmetj1_isLoaded = true;
  }
  return dphi_genmetj1_;
}

const float &WWWTree::dphi_genmetj2() {
  if (not dphi_genmetj2_isLoaded) {
    if (dphi_genmetj2_branch != 0) {
      dphi_genmetj2_branch->GetEntry(index);
    } else {
      printf("branch dphi_genmetj2_branch does not exist!\n");
      exit(1);
    }
    dphi_genmetj2_isLoaded = true;
  }
  return dphi_genmetj2_;
}

const float &WWWTree::mjj_mindphi_up() {
  if (not mjj_mindphi_up_isLoaded) {
    if (mjj_mindphi_up_branch != 0) {
      mjj_mindphi_up_branch->GetEntry(index);
    } else {
      printf("branch mjj_mindphi_up_branch does not exist!\n");
      exit(1);
    }
    mjj_mindphi_up_isLoaded = true;
  }
  return mjj_mindphi_up_;
}

const float &WWWTree::mjj_up() {
  if (not mjj_up_isLoaded) {
    if (mjj_up_branch != 0) {
      mjj_up_branch->GetEntry(index);
    } else {
      printf("branch mjj_up_branch does not exist!\n");
      exit(1);
    }
    mjj_up_isLoaded = true;
  }
  return mjj_up_;
}

const float &WWWTree::mbb_csv_up() {
  if (not mbb_csv_up_isLoaded) {
    if (mbb_csv_up_branch != 0) {
      mbb_csv_up_branch->GetEntry(index);
    } else {
      printf("branch mbb_csv_up_branch does not exist!\n");
      exit(1);
    }
    mbb_csv_up_isLoaded = true;
  }
  return mbb_csv_up_;
}

const float &WWWTree::mbb_bpt_up() {
  if (not mbb_bpt_up_isLoaded) {
    if (mbb_bpt_up_branch != 0) {
      mbb_bpt_up_branch->GetEntry(index);
    } else {
      printf("branch mbb_bpt_up_branch does not exist!\n");
      exit(1);
    }
    mbb_bpt_up_isLoaded = true;
  }
  return mbb_bpt_up_;
}

const float &WWWTree::dphi_jj_up() {
  if (not dphi_jj_up_isLoaded) {
    if (dphi_jj_up_branch != 0) {
      dphi_jj_up_branch->GetEntry(index);
    } else {
      printf("branch dphi_jj_up_branch does not exist!\n");
      exit(1);
    }
    dphi_jj_up_isLoaded = true;
  }
  return dphi_jj_up_;
}

const float &WWWTree::dphi_ll_up() {
  if (not dphi_ll_up_isLoaded) {
    if (dphi_ll_up_branch != 0) {
      dphi_ll_up_branch->GetEntry(index);
    } else {
      printf("branch dphi_ll_up_branch does not exist!\n");
      exit(1);
    }
    dphi_ll_up_isLoaded = true;
  }
  return dphi_ll_up_;
}

const float &WWWTree::sum_mlb_up() {
  if (not sum_mlb_up_isLoaded) {
    if (sum_mlb_up_branch != 0) {
      sum_mlb_up_branch->GetEntry(index);
    } else {
      printf("branch sum_mlb_up_branch does not exist!\n");
      exit(1);
    }
    sum_mlb_up_isLoaded = true;
  }
  return sum_mlb_up_;
}

const float &WWWTree::deta_jj_up() {
  if (not deta_jj_up_isLoaded) {
    if (deta_jj_up_branch != 0) {
      deta_jj_up_branch->GetEntry(index);
    } else {
      printf("branch deta_jj_up_branch does not exist!\n");
      exit(1);
    }
    deta_jj_up_isLoaded = true;
  }
  return deta_jj_up_;
}

const float &WWWTree::dR_jj_up() {
  if (not dR_jj_up_isLoaded) {
    if (dR_jj_up_branch != 0) {
      dR_jj_up_branch->GetEntry(index);
    } else {
      printf("branch dR_jj_up_branch does not exist!\n");
      exit(1);
    }
    dR_jj_up_isLoaded = true;
  }
  return dR_jj_up_;
}

const float &WWWTree::dphi_metj1_up() {
  if (not dphi_metj1_up_isLoaded) {
    if (dphi_metj1_up_branch != 0) {
      dphi_metj1_up_branch->GetEntry(index);
    } else {
      printf("branch dphi_metj1_up_branch does not exist!\n");
      exit(1);
    }
    dphi_metj1_up_isLoaded = true;
  }
  return dphi_metj1_up_;
}

const float &WWWTree::dphi_metj2_up() {
  if (not dphi_metj2_up_isLoaded) {
    if (dphi_metj2_up_branch != 0) {
      dphi_metj2_up_branch->GetEntry(index);
    } else {
      printf("branch dphi_metj2_up_branch does not exist!\n");
      exit(1);
    }
    dphi_metj2_up_isLoaded = true;
  }
  return dphi_metj2_up_;
}

const float &WWWTree::mjj_mindphi_dn() {
  if (not mjj_mindphi_dn_isLoaded) {
    if (mjj_mindphi_dn_branch != 0) {
      mjj_mindphi_dn_branch->GetEntry(index);
    } else {
      printf("branch mjj_mindphi_dn_branch does not exist!\n");
      exit(1);
    }
    mjj_mindphi_dn_isLoaded = true;
  }
  return mjj_mindphi_dn_;
}

const float &WWWTree::mjj_dn() {
  if (not mjj_dn_isLoaded) {
    if (mjj_dn_branch != 0) {
      mjj_dn_branch->GetEntry(index);
    } else {
      printf("branch mjj_dn_branch does not exist!\n");
      exit(1);
    }
    mjj_dn_isLoaded = true;
  }
  return mjj_dn_;
}

const float &WWWTree::mbb_csv_dn() {
  if (not mbb_csv_dn_isLoaded) {
    if (mbb_csv_dn_branch != 0) {
      mbb_csv_dn_branch->GetEntry(index);
    } else {
      printf("branch mbb_csv_dn_branch does not exist!\n");
      exit(1);
    }
    mbb_csv_dn_isLoaded = true;
  }
  return mbb_csv_dn_;
}

const float &WWWTree::mbb_bpt_dn() {
  if (not mbb_bpt_dn_isLoaded) {
    if (mbb_bpt_dn_branch != 0) {
      mbb_bpt_dn_branch->GetEntry(index);
    } else {
      printf("branch mbb_bpt_dn_branch does not exist!\n");
      exit(1);
    }
    mbb_bpt_dn_isLoaded = true;
  }
  return mbb_bpt_dn_;
}

const float &WWWTree::dphi_jj_dn() {
  if (not dphi_jj_dn_isLoaded) {
    if (dphi_jj_dn_branch != 0) {
      dphi_jj_dn_branch->GetEntry(index);
    } else {
      printf("branch dphi_jj_dn_branch does not exist!\n");
      exit(1);
    }
    dphi_jj_dn_isLoaded = true;
  }
  return dphi_jj_dn_;
}

const float &WWWTree::dphi_ll_dn() {
  if (not dphi_ll_dn_isLoaded) {
    if (dphi_ll_dn_branch != 0) {
      dphi_ll_dn_branch->GetEntry(index);
    } else {
      printf("branch dphi_ll_dn_branch does not exist!\n");
      exit(1);
    }
    dphi_ll_dn_isLoaded = true;
  }
  return dphi_ll_dn_;
}

const float &WWWTree::sum_mlb_dn() {
  if (not sum_mlb_dn_isLoaded) {
    if (sum_mlb_dn_branch != 0) {
      sum_mlb_dn_branch->GetEntry(index);
    } else {
      printf("branch sum_mlb_dn_branch does not exist!\n");
      exit(1);
    }
    sum_mlb_dn_isLoaded = true;
  }
  return sum_mlb_dn_;
}

const float &WWWTree::deta_jj_dn() {
  if (not deta_jj_dn_isLoaded) {
    if (deta_jj_dn_branch != 0) {
      deta_jj_dn_branch->GetEntry(index);
    } else {
      printf("branch deta_jj_dn_branch does not exist!\n");
      exit(1);
    }
    deta_jj_dn_isLoaded = true;
  }
  return deta_jj_dn_;
}

const float &WWWTree::dR_jj_dn() {
  if (not dR_jj_dn_isLoaded) {
    if (dR_jj_dn_branch != 0) {
      dR_jj_dn_branch->GetEntry(index);
    } else {
      printf("branch dR_jj_dn_branch does not exist!\n");
      exit(1);
    }
    dR_jj_dn_isLoaded = true;
  }
  return dR_jj_dn_;
}

const float &WWWTree::dphi_metj1_dn() {
  if (not dphi_metj1_dn_isLoaded) {
    if (dphi_metj1_dn_branch != 0) {
      dphi_metj1_dn_branch->GetEntry(index);
    } else {
      printf("branch dphi_metj1_dn_branch does not exist!\n");
      exit(1);
    }
    dphi_metj1_dn_isLoaded = true;
  }
  return dphi_metj1_dn_;
}

const float &WWWTree::dphi_metj2_dn() {
  if (not dphi_metj2_dn_isLoaded) {
    if (dphi_metj2_dn_branch != 0) {
      dphi_metj2_dn_branch->GetEntry(index);
    } else {
      printf("branch dphi_metj2_dn_branch does not exist!\n");
      exit(1);
    }
    dphi_metj2_dn_isLoaded = true;
  }
  return dphi_metj2_dn_;
}

const float &WWWTree::weight_btagsf() {
  if (not weight_btagsf_isLoaded) {
    if (weight_btagsf_branch != 0) {
      weight_btagsf_branch->GetEntry(index);
    } else {
      printf("branch weight_btagsf_branch does not exist!\n");
      exit(1);
    }
    weight_btagsf_isLoaded = true;
  }
  return weight_btagsf_;
}

const float &WWWTree::weight_btagsf_heavy_UP() {
  if (not weight_btagsf_heavy_UP_isLoaded) {
    if (weight_btagsf_heavy_UP_branch != 0) {
      weight_btagsf_heavy_UP_branch->GetEntry(index);
    } else {
      printf("branch weight_btagsf_heavy_UP_branch does not exist!\n");
      exit(1);
    }
    weight_btagsf_heavy_UP_isLoaded = true;
  }
  return weight_btagsf_heavy_UP_;
}

const float &WWWTree::weight_btagsf_light_UP() {
  if (not weight_btagsf_light_UP_isLoaded) {
    if (weight_btagsf_light_UP_branch != 0) {
      weight_btagsf_light_UP_branch->GetEntry(index);
    } else {
      printf("branch weight_btagsf_light_UP_branch does not exist!\n");
      exit(1);
    }
    weight_btagsf_light_UP_isLoaded = true;
  }
  return weight_btagsf_light_UP_;
}

const float &WWWTree::weight_btagsf_heavy_DN() {
  if (not weight_btagsf_heavy_DN_isLoaded) {
    if (weight_btagsf_heavy_DN_branch != 0) {
      weight_btagsf_heavy_DN_branch->GetEntry(index);
    } else {
      printf("branch weight_btagsf_heavy_DN_branch does not exist!\n");
      exit(1);
    }
    weight_btagsf_heavy_DN_isLoaded = true;
  }
  return weight_btagsf_heavy_DN_;
}

const float &WWWTree::weight_btagsf_light_DN() {
  if (not weight_btagsf_light_DN_isLoaded) {
    if (weight_btagsf_light_DN_branch != 0) {
      weight_btagsf_light_DN_branch->GetEntry(index);
    } else {
      printf("branch weight_btagsf_light_DN_branch does not exist!\n");
      exit(1);
    }
    weight_btagsf_light_DN_isLoaded = true;
  }
  return weight_btagsf_light_DN_;
}

const float &WWWTree::chpfcands_0013_pt() {
  if (not chpfcands_0013_pt_isLoaded) {
    if (chpfcands_0013_pt_branch != 0) {
      chpfcands_0013_pt_branch->GetEntry(index);
    } else {
      printf("branch chpfcands_0013_pt_branch does not exist!\n");
      exit(1);
    }
    chpfcands_0013_pt_isLoaded = true;
  }
  return chpfcands_0013_pt_;
}

const float &WWWTree::chpfcands_1316_pt() {
  if (not chpfcands_1316_pt_isLoaded) {
    if (chpfcands_1316_pt_branch != 0) {
      chpfcands_1316_pt_branch->GetEntry(index);
    } else {
      printf("branch chpfcands_1316_pt_branch does not exist!\n");
      exit(1);
    }
    chpfcands_1316_pt_isLoaded = true;
  }
  return chpfcands_1316_pt_;
}

const float &WWWTree::chpfcands_1624_pt() {
  if (not chpfcands_1624_pt_isLoaded) {
    if (chpfcands_1624_pt_branch != 0) {
      chpfcands_1624_pt_branch->GetEntry(index);
    } else {
      printf("branch chpfcands_1624_pt_branch does not exist!\n");
      exit(1);
    }
    chpfcands_1624_pt_isLoaded = true;
  }
  return chpfcands_1624_pt_;
}

const float &WWWTree::chpfcands_2430_pt() {
  if (not chpfcands_2430_pt_isLoaded) {
    if (chpfcands_2430_pt_branch != 0) {
      chpfcands_2430_pt_branch->GetEntry(index);
    } else {
      printf("branch chpfcands_2430_pt_branch does not exist!\n");
      exit(1);
    }
    chpfcands_2430_pt_isLoaded = true;
  }
  return chpfcands_2430_pt_;
}

const float &WWWTree::chpfcands_30in_pt() {
  if (not chpfcands_30in_pt_isLoaded) {
    if (chpfcands_30in_pt_branch != 0) {
      chpfcands_30in_pt_branch->GetEntry(index);
    } else {
      printf("branch chpfcands_30in_pt_branch does not exist!\n");
      exit(1);
    }
    chpfcands_30in_pt_isLoaded = true;
  }
  return chpfcands_30in_pt_;
}

const float &WWWTree::phpfcands_0013_pt() {
  if (not phpfcands_0013_pt_isLoaded) {
    if (phpfcands_0013_pt_branch != 0) {
      phpfcands_0013_pt_branch->GetEntry(index);
    } else {
      printf("branch phpfcands_0013_pt_branch does not exist!\n");
      exit(1);
    }
    phpfcands_0013_pt_isLoaded = true;
  }
  return phpfcands_0013_pt_;
}

const float &WWWTree::phpfcands_1316_pt() {
  if (not phpfcands_1316_pt_isLoaded) {
    if (phpfcands_1316_pt_branch != 0) {
      phpfcands_1316_pt_branch->GetEntry(index);
    } else {
      printf("branch phpfcands_1316_pt_branch does not exist!\n");
      exit(1);
    }
    phpfcands_1316_pt_isLoaded = true;
  }
  return phpfcands_1316_pt_;
}

const float &WWWTree::phpfcands_1624_pt() {
  if (not phpfcands_1624_pt_isLoaded) {
    if (phpfcands_1624_pt_branch != 0) {
      phpfcands_1624_pt_branch->GetEntry(index);
    } else {
      printf("branch phpfcands_1624_pt_branch does not exist!\n");
      exit(1);
    }
    phpfcands_1624_pt_isLoaded = true;
  }
  return phpfcands_1624_pt_;
}

const float &WWWTree::phpfcands_2430_pt() {
  if (not phpfcands_2430_pt_isLoaded) {
    if (phpfcands_2430_pt_branch != 0) {
      phpfcands_2430_pt_branch->GetEntry(index);
    } else {
      printf("branch phpfcands_2430_pt_branch does not exist!\n");
      exit(1);
    }
    phpfcands_2430_pt_isLoaded = true;
  }
  return phpfcands_2430_pt_;
}

const float &WWWTree::phpfcands_30in_pt() {
  if (not phpfcands_30in_pt_isLoaded) {
    if (phpfcands_30in_pt_branch != 0) {
      phpfcands_30in_pt_branch->GetEntry(index);
    } else {
      printf("branch phpfcands_30in_pt_branch does not exist!\n");
      exit(1);
    }
    phpfcands_30in_pt_isLoaded = true;
  }
  return phpfcands_30in_pt_;
}

const float &WWWTree::nupfcands_0013_pt() {
  if (not nupfcands_0013_pt_isLoaded) {
    if (nupfcands_0013_pt_branch != 0) {
      nupfcands_0013_pt_branch->GetEntry(index);
    } else {
      printf("branch nupfcands_0013_pt_branch does not exist!\n");
      exit(1);
    }
    nupfcands_0013_pt_isLoaded = true;
  }
  return nupfcands_0013_pt_;
}

const float &WWWTree::nupfcands_1316_pt() {
  if (not nupfcands_1316_pt_isLoaded) {
    if (nupfcands_1316_pt_branch != 0) {
      nupfcands_1316_pt_branch->GetEntry(index);
    } else {
      printf("branch nupfcands_1316_pt_branch does not exist!\n");
      exit(1);
    }
    nupfcands_1316_pt_isLoaded = true;
  }
  return nupfcands_1316_pt_;
}

const float &WWWTree::nupfcands_1624_pt() {
  if (not nupfcands_1624_pt_isLoaded) {
    if (nupfcands_1624_pt_branch != 0) {
      nupfcands_1624_pt_branch->GetEntry(index);
    } else {
      printf("branch nupfcands_1624_pt_branch does not exist!\n");
      exit(1);
    }
    nupfcands_1624_pt_isLoaded = true;
  }
  return nupfcands_1624_pt_;
}

const float &WWWTree::nupfcands_2430_pt() {
  if (not nupfcands_2430_pt_isLoaded) {
    if (nupfcands_2430_pt_branch != 0) {
      nupfcands_2430_pt_branch->GetEntry(index);
    } else {
      printf("branch nupfcands_2430_pt_branch does not exist!\n");
      exit(1);
    }
    nupfcands_2430_pt_isLoaded = true;
  }
  return nupfcands_2430_pt_;
}

const float &WWWTree::nupfcands_30in_pt() {
  if (not nupfcands_30in_pt_isLoaded) {
    if (nupfcands_30in_pt_branch != 0) {
      nupfcands_30in_pt_branch->GetEntry(index);
    } else {
      printf("branch nupfcands_30in_pt_branch does not exist!\n");
      exit(1);
    }
    nupfcands_30in_pt_isLoaded = true;
  }
  return nupfcands_30in_pt_;
}

const float &WWWTree::chpfcands_0013_sumet() {
  if (not chpfcands_0013_sumet_isLoaded) {
    if (chpfcands_0013_sumet_branch != 0) {
      chpfcands_0013_sumet_branch->GetEntry(index);
    } else {
      printf("branch chpfcands_0013_sumet_branch does not exist!\n");
      exit(1);
    }
    chpfcands_0013_sumet_isLoaded = true;
  }
  return chpfcands_0013_sumet_;
}

const float &WWWTree::chpfcands_1316_sumet() {
  if (not chpfcands_1316_sumet_isLoaded) {
    if (chpfcands_1316_sumet_branch != 0) {
      chpfcands_1316_sumet_branch->GetEntry(index);
    } else {
      printf("branch chpfcands_1316_sumet_branch does not exist!\n");
      exit(1);
    }
    chpfcands_1316_sumet_isLoaded = true;
  }
  return chpfcands_1316_sumet_;
}

const float &WWWTree::chpfcands_1624_sumet() {
  if (not chpfcands_1624_sumet_isLoaded) {
    if (chpfcands_1624_sumet_branch != 0) {
      chpfcands_1624_sumet_branch->GetEntry(index);
    } else {
      printf("branch chpfcands_1624_sumet_branch does not exist!\n");
      exit(1);
    }
    chpfcands_1624_sumet_isLoaded = true;
  }
  return chpfcands_1624_sumet_;
}

const float &WWWTree::chpfcands_2430_sumet() {
  if (not chpfcands_2430_sumet_isLoaded) {
    if (chpfcands_2430_sumet_branch != 0) {
      chpfcands_2430_sumet_branch->GetEntry(index);
    } else {
      printf("branch chpfcands_2430_sumet_branch does not exist!\n");
      exit(1);
    }
    chpfcands_2430_sumet_isLoaded = true;
  }
  return chpfcands_2430_sumet_;
}

const float &WWWTree::chpfcands_30in_sumet() {
  if (not chpfcands_30in_sumet_isLoaded) {
    if (chpfcands_30in_sumet_branch != 0) {
      chpfcands_30in_sumet_branch->GetEntry(index);
    } else {
      printf("branch chpfcands_30in_sumet_branch does not exist!\n");
      exit(1);
    }
    chpfcands_30in_sumet_isLoaded = true;
  }
  return chpfcands_30in_sumet_;
}

const float &WWWTree::phpfcands_0013_sumet() {
  if (not phpfcands_0013_sumet_isLoaded) {
    if (phpfcands_0013_sumet_branch != 0) {
      phpfcands_0013_sumet_branch->GetEntry(index);
    } else {
      printf("branch phpfcands_0013_sumet_branch does not exist!\n");
      exit(1);
    }
    phpfcands_0013_sumet_isLoaded = true;
  }
  return phpfcands_0013_sumet_;
}

const float &WWWTree::phpfcands_1316_sumet() {
  if (not phpfcands_1316_sumet_isLoaded) {
    if (phpfcands_1316_sumet_branch != 0) {
      phpfcands_1316_sumet_branch->GetEntry(index);
    } else {
      printf("branch phpfcands_1316_sumet_branch does not exist!\n");
      exit(1);
    }
    phpfcands_1316_sumet_isLoaded = true;
  }
  return phpfcands_1316_sumet_;
}

const float &WWWTree::phpfcands_1624_sumet() {
  if (not phpfcands_1624_sumet_isLoaded) {
    if (phpfcands_1624_sumet_branch != 0) {
      phpfcands_1624_sumet_branch->GetEntry(index);
    } else {
      printf("branch phpfcands_1624_sumet_branch does not exist!\n");
      exit(1);
    }
    phpfcands_1624_sumet_isLoaded = true;
  }
  return phpfcands_1624_sumet_;
}

const float &WWWTree::phpfcands_2430_sumet() {
  if (not phpfcands_2430_sumet_isLoaded) {
    if (phpfcands_2430_sumet_branch != 0) {
      phpfcands_2430_sumet_branch->GetEntry(index);
    } else {
      printf("branch phpfcands_2430_sumet_branch does not exist!\n");
      exit(1);
    }
    phpfcands_2430_sumet_isLoaded = true;
  }
  return phpfcands_2430_sumet_;
}

const float &WWWTree::phpfcands_30in_sumet() {
  if (not phpfcands_30in_sumet_isLoaded) {
    if (phpfcands_30in_sumet_branch != 0) {
      phpfcands_30in_sumet_branch->GetEntry(index);
    } else {
      printf("branch phpfcands_30in_sumet_branch does not exist!\n");
      exit(1);
    }
    phpfcands_30in_sumet_isLoaded = true;
  }
  return phpfcands_30in_sumet_;
}

const float &WWWTree::nupfcands_0013_sumet() {
  if (not nupfcands_0013_sumet_isLoaded) {
    if (nupfcands_0013_sumet_branch != 0) {
      nupfcands_0013_sumet_branch->GetEntry(index);
    } else {
      printf("branch nupfcands_0013_sumet_branch does not exist!\n");
      exit(1);
    }
    nupfcands_0013_sumet_isLoaded = true;
  }
  return nupfcands_0013_sumet_;
}

const float &WWWTree::nupfcands_1316_sumet() {
  if (not nupfcands_1316_sumet_isLoaded) {
    if (nupfcands_1316_sumet_branch != 0) {
      nupfcands_1316_sumet_branch->GetEntry(index);
    } else {
      printf("branch nupfcands_1316_sumet_branch does not exist!\n");
      exit(1);
    }
    nupfcands_1316_sumet_isLoaded = true;
  }
  return nupfcands_1316_sumet_;
}

const float &WWWTree::nupfcands_1624_sumet() {
  if (not nupfcands_1624_sumet_isLoaded) {
    if (nupfcands_1624_sumet_branch != 0) {
      nupfcands_1624_sumet_branch->GetEntry(index);
    } else {
      printf("branch nupfcands_1624_sumet_branch does not exist!\n");
      exit(1);
    }
    nupfcands_1624_sumet_isLoaded = true;
  }
  return nupfcands_1624_sumet_;
}

const float &WWWTree::nupfcands_2430_sumet() {
  if (not nupfcands_2430_sumet_isLoaded) {
    if (nupfcands_2430_sumet_branch != 0) {
      nupfcands_2430_sumet_branch->GetEntry(index);
    } else {
      printf("branch nupfcands_2430_sumet_branch does not exist!\n");
      exit(1);
    }
    nupfcands_2430_sumet_isLoaded = true;
  }
  return nupfcands_2430_sumet_;
}

const float &WWWTree::nupfcands_30in_sumet() {
  if (not nupfcands_30in_sumet_isLoaded) {
    if (nupfcands_30in_sumet_branch != 0) {
      nupfcands_30in_sumet_branch->GetEntry(index);
    } else {
      printf("branch nupfcands_30in_sumet_branch does not exist!\n");
      exit(1);
    }
    nupfcands_30in_sumet_isLoaded = true;
  }
  return nupfcands_30in_sumet_;
}

const float &WWWTree::chpfcands_0013_phi() {
  if (not chpfcands_0013_phi_isLoaded) {
    if (chpfcands_0013_phi_branch != 0) {
      chpfcands_0013_phi_branch->GetEntry(index);
    } else {
      printf("branch chpfcands_0013_phi_branch does not exist!\n");
      exit(1);
    }
    chpfcands_0013_phi_isLoaded = true;
  }
  return chpfcands_0013_phi_;
}

const float &WWWTree::chpfcands_1316_phi() {
  if (not chpfcands_1316_phi_isLoaded) {
    if (chpfcands_1316_phi_branch != 0) {
      chpfcands_1316_phi_branch->GetEntry(index);
    } else {
      printf("branch chpfcands_1316_phi_branch does not exist!\n");
      exit(1);
    }
    chpfcands_1316_phi_isLoaded = true;
  }
  return chpfcands_1316_phi_;
}

const float &WWWTree::chpfcands_1624_phi() {
  if (not chpfcands_1624_phi_isLoaded) {
    if (chpfcands_1624_phi_branch != 0) {
      chpfcands_1624_phi_branch->GetEntry(index);
    } else {
      printf("branch chpfcands_1624_phi_branch does not exist!\n");
      exit(1);
    }
    chpfcands_1624_phi_isLoaded = true;
  }
  return chpfcands_1624_phi_;
}

const float &WWWTree::chpfcands_2430_phi() {
  if (not chpfcands_2430_phi_isLoaded) {
    if (chpfcands_2430_phi_branch != 0) {
      chpfcands_2430_phi_branch->GetEntry(index);
    } else {
      printf("branch chpfcands_2430_phi_branch does not exist!\n");
      exit(1);
    }
    chpfcands_2430_phi_isLoaded = true;
  }
  return chpfcands_2430_phi_;
}

const float &WWWTree::chpfcands_30in_phi() {
  if (not chpfcands_30in_phi_isLoaded) {
    if (chpfcands_30in_phi_branch != 0) {
      chpfcands_30in_phi_branch->GetEntry(index);
    } else {
      printf("branch chpfcands_30in_phi_branch does not exist!\n");
      exit(1);
    }
    chpfcands_30in_phi_isLoaded = true;
  }
  return chpfcands_30in_phi_;
}

const float &WWWTree::phpfcands_0013_phi() {
  if (not phpfcands_0013_phi_isLoaded) {
    if (phpfcands_0013_phi_branch != 0) {
      phpfcands_0013_phi_branch->GetEntry(index);
    } else {
      printf("branch phpfcands_0013_phi_branch does not exist!\n");
      exit(1);
    }
    phpfcands_0013_phi_isLoaded = true;
  }
  return phpfcands_0013_phi_;
}

const float &WWWTree::phpfcands_1316_phi() {
  if (not phpfcands_1316_phi_isLoaded) {
    if (phpfcands_1316_phi_branch != 0) {
      phpfcands_1316_phi_branch->GetEntry(index);
    } else {
      printf("branch phpfcands_1316_phi_branch does not exist!\n");
      exit(1);
    }
    phpfcands_1316_phi_isLoaded = true;
  }
  return phpfcands_1316_phi_;
}

const float &WWWTree::phpfcands_1624_phi() {
  if (not phpfcands_1624_phi_isLoaded) {
    if (phpfcands_1624_phi_branch != 0) {
      phpfcands_1624_phi_branch->GetEntry(index);
    } else {
      printf("branch phpfcands_1624_phi_branch does not exist!\n");
      exit(1);
    }
    phpfcands_1624_phi_isLoaded = true;
  }
  return phpfcands_1624_phi_;
}

const float &WWWTree::phpfcands_2430_phi() {
  if (not phpfcands_2430_phi_isLoaded) {
    if (phpfcands_2430_phi_branch != 0) {
      phpfcands_2430_phi_branch->GetEntry(index);
    } else {
      printf("branch phpfcands_2430_phi_branch does not exist!\n");
      exit(1);
    }
    phpfcands_2430_phi_isLoaded = true;
  }
  return phpfcands_2430_phi_;
}

const float &WWWTree::phpfcands_30in_phi() {
  if (not phpfcands_30in_phi_isLoaded) {
    if (phpfcands_30in_phi_branch != 0) {
      phpfcands_30in_phi_branch->GetEntry(index);
    } else {
      printf("branch phpfcands_30in_phi_branch does not exist!\n");
      exit(1);
    }
    phpfcands_30in_phi_isLoaded = true;
  }
  return phpfcands_30in_phi_;
}

const float &WWWTree::nupfcands_0013_phi() {
  if (not nupfcands_0013_phi_isLoaded) {
    if (nupfcands_0013_phi_branch != 0) {
      nupfcands_0013_phi_branch->GetEntry(index);
    } else {
      printf("branch nupfcands_0013_phi_branch does not exist!\n");
      exit(1);
    }
    nupfcands_0013_phi_isLoaded = true;
  }
  return nupfcands_0013_phi_;
}

const float &WWWTree::nupfcands_1316_phi() {
  if (not nupfcands_1316_phi_isLoaded) {
    if (nupfcands_1316_phi_branch != 0) {
      nupfcands_1316_phi_branch->GetEntry(index);
    } else {
      printf("branch nupfcands_1316_phi_branch does not exist!\n");
      exit(1);
    }
    nupfcands_1316_phi_isLoaded = true;
  }
  return nupfcands_1316_phi_;
}

const float &WWWTree::nupfcands_1624_phi() {
  if (not nupfcands_1624_phi_isLoaded) {
    if (nupfcands_1624_phi_branch != 0) {
      nupfcands_1624_phi_branch->GetEntry(index);
    } else {
      printf("branch nupfcands_1624_phi_branch does not exist!\n");
      exit(1);
    }
    nupfcands_1624_phi_isLoaded = true;
  }
  return nupfcands_1624_phi_;
}

const float &WWWTree::nupfcands_2430_phi() {
  if (not nupfcands_2430_phi_isLoaded) {
    if (nupfcands_2430_phi_branch != 0) {
      nupfcands_2430_phi_branch->GetEntry(index);
    } else {
      printf("branch nupfcands_2430_phi_branch does not exist!\n");
      exit(1);
    }
    nupfcands_2430_phi_isLoaded = true;
  }
  return nupfcands_2430_phi_;
}

const float &WWWTree::nupfcands_30in_phi() {
  if (not nupfcands_30in_phi_isLoaded) {
    if (nupfcands_30in_phi_branch != 0) {
      nupfcands_30in_phi_branch->GetEntry(index);
    } else {
      printf("branch nupfcands_30in_phi_branch does not exist!\n");
      exit(1);
    }
    nupfcands_30in_phi_isLoaded = true;
  }
  return nupfcands_30in_phi_;
}

const float &WWWTree::met_T1CHS_pt() {
  if (not met_T1CHS_pt_isLoaded) {
    if (met_T1CHS_pt_branch != 0) {
      met_T1CHS_pt_branch->GetEntry(index);
    } else {
      printf("branch met_T1CHS_pt_branch does not exist!\n");
      exit(1);
    }
    met_T1CHS_pt_isLoaded = true;
  }
  return met_T1CHS_pt_;
}

const float &WWWTree::met_T1CHS_phi() {
  if (not met_T1CHS_phi_isLoaded) {
    if (met_T1CHS_phi_branch != 0) {
      met_T1CHS_phi_branch->GetEntry(index);
    } else {
      printf("branch met_T1CHS_phi_branch does not exist!\n");
      exit(1);
    }
    met_T1CHS_phi_isLoaded = true;
  }
  return met_T1CHS_phi_;
}

const float &WWWTree::met_T1CHS_fromCORE_pt() {
  if (not met_T1CHS_fromCORE_pt_isLoaded) {
    if (met_T1CHS_fromCORE_pt_branch != 0) {
      met_T1CHS_fromCORE_pt_branch->GetEntry(index);
    } else {
      printf("branch met_T1CHS_fromCORE_pt_branch does not exist!\n");
      exit(1);
    }
    met_T1CHS_fromCORE_pt_isLoaded = true;
  }
  return met_T1CHS_fromCORE_pt_;
}

const float &WWWTree::met_T1CHS_fromCORE_phi() {
  if (not met_T1CHS_fromCORE_phi_isLoaded) {
    if (met_T1CHS_fromCORE_phi_branch != 0) {
      met_T1CHS_fromCORE_phi_branch->GetEntry(index);
    } else {
      printf("branch met_T1CHS_fromCORE_phi_branch does not exist!\n");
      exit(1);
    }
    met_T1CHS_fromCORE_phi_isLoaded = true;
  }
  return met_T1CHS_fromCORE_phi_;
}

const float &WWWTree::met_T1CHS_miniAOD_CORE_pt() {
  if (not met_T1CHS_miniAOD_CORE_pt_isLoaded) {
    if (met_T1CHS_miniAOD_CORE_pt_branch != 0) {
      met_T1CHS_miniAOD_CORE_pt_branch->GetEntry(index);
    } else {
      printf("branch met_T1CHS_miniAOD_CORE_pt_branch does not exist!\n");
      exit(1);
    }
    met_T1CHS_miniAOD_CORE_pt_isLoaded = true;
  }
  return met_T1CHS_miniAOD_CORE_pt_;
}

const float &WWWTree::met_T1CHS_miniAOD_CORE_phi() {
  if (not met_T1CHS_miniAOD_CORE_phi_isLoaded) {
    if (met_T1CHS_miniAOD_CORE_phi_branch != 0) {
      met_T1CHS_miniAOD_CORE_phi_branch->GetEntry(index);
    } else {
      printf("branch met_T1CHS_miniAOD_CORE_phi_branch does not exist!\n");
      exit(1);
    }
    met_T1CHS_miniAOD_CORE_phi_isLoaded = true;
  }
  return met_T1CHS_miniAOD_CORE_phi_;
}

const float &WWWTree::met_T1CHS_miniAOD_CORE_up_pt() {
  if (not met_T1CHS_miniAOD_CORE_up_pt_isLoaded) {
    if (met_T1CHS_miniAOD_CORE_up_pt_branch != 0) {
      met_T1CHS_miniAOD_CORE_up_pt_branch->GetEntry(index);
    } else {
      printf("branch met_T1CHS_miniAOD_CORE_up_pt_branch does not exist!\n");
      exit(1);
    }
    met_T1CHS_miniAOD_CORE_up_pt_isLoaded = true;
  }
  return met_T1CHS_miniAOD_CORE_up_pt_;
}

const float &WWWTree::met_T1CHS_miniAOD_CORE_up_phi() {
  if (not met_T1CHS_miniAOD_CORE_up_phi_isLoaded) {
    if (met_T1CHS_miniAOD_CORE_up_phi_branch != 0) {
      met_T1CHS_miniAOD_CORE_up_phi_branch->GetEntry(index);
    } else {
      printf("branch met_T1CHS_miniAOD_CORE_up_phi_branch does not exist!\n");
      exit(1);
    }
    met_T1CHS_miniAOD_CORE_up_phi_isLoaded = true;
  }
  return met_T1CHS_miniAOD_CORE_up_phi_;
}

const float &WWWTree::met_T1CHS_miniAOD_CORE_dn_pt() {
  if (not met_T1CHS_miniAOD_CORE_dn_pt_isLoaded) {
    if (met_T1CHS_miniAOD_CORE_dn_pt_branch != 0) {
      met_T1CHS_miniAOD_CORE_dn_pt_branch->GetEntry(index);
    } else {
      printf("branch met_T1CHS_miniAOD_CORE_dn_pt_branch does not exist!\n");
      exit(1);
    }
    met_T1CHS_miniAOD_CORE_dn_pt_isLoaded = true;
  }
  return met_T1CHS_miniAOD_CORE_dn_pt_;
}

const float &WWWTree::met_T1CHS_miniAOD_CORE_dn_phi() {
  if (not met_T1CHS_miniAOD_CORE_dn_phi_isLoaded) {
    if (met_T1CHS_miniAOD_CORE_dn_phi_branch != 0) {
      met_T1CHS_miniAOD_CORE_dn_phi_branch->GetEntry(index);
    } else {
      printf("branch met_T1CHS_miniAOD_CORE_dn_phi_branch does not exist!\n");
      exit(1);
    }
    met_T1CHS_miniAOD_CORE_dn_phi_isLoaded = true;
  }
  return met_T1CHS_miniAOD_CORE_dn_phi_;
}

const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &WWWTree::decayedphoton_lep1_p4() {
  if (not decayedphoton_lep1_p4_isLoaded) {
    if (decayedphoton_lep1_p4_branch != 0) {
      decayedphoton_lep1_p4_branch->GetEntry(index);
    } else {
      printf("branch decayedphoton_lep1_p4_branch does not exist!\n");
      exit(1);
    }
    decayedphoton_lep1_p4_isLoaded = true;
  }
  return *decayedphoton_lep1_p4_;
}

const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &WWWTree::decayedphoton_lep2_p4() {
  if (not decayedphoton_lep2_p4_isLoaded) {
    if (decayedphoton_lep2_p4_branch != 0) {
      decayedphoton_lep2_p4_branch->GetEntry(index);
    } else {
      printf("branch decayedphoton_lep2_p4_branch does not exist!\n");
      exit(1);
    }
    decayedphoton_lep2_p4_isLoaded = true;
  }
  return *decayedphoton_lep2_p4_;
}

const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &WWWTree::decayedphoton_bosn_p4() {
  if (not decayedphoton_bosn_p4_isLoaded) {
    if (decayedphoton_bosn_p4_branch != 0) {
      decayedphoton_bosn_p4_branch->GetEntry(index);
    } else {
      printf("branch decayedphoton_bosn_p4_branch does not exist!\n");
      exit(1);
    }
    decayedphoton_bosn_p4_isLoaded = true;
  }
  return *decayedphoton_bosn_p4_;
}

const float &WWWTree::decayedphoton_mt2() {
  if (not decayedphoton_mt2_isLoaded) {
    if (decayedphoton_mt2_branch != 0) {
      decayedphoton_mt2_branch->GetEntry(index);
    } else {
      printf("branch decayedphoton_mt2_branch does not exist!\n");
      exit(1);
    }
    decayedphoton_mt2_isLoaded = true;
  }
  return decayedphoton_mt2_;
}

const int &WWWTree::hyp_type() {
  if (not hyp_type_isLoaded) {
    if (hyp_type_branch != 0) {
      hyp_type_branch->GetEntry(index);
    } else {
      printf("branch hyp_type_branch does not exist!\n");
      exit(1);
    }
    hyp_type_isLoaded = true;
  }
  return hyp_type_;
}

const int &WWWTree::evt_type() {
  if (not evt_type_isLoaded) {
    if (evt_type_branch != 0) {
      evt_type_branch->GetEntry(index);
    } else {
      printf("branch evt_type_branch does not exist!\n");
      exit(1);
    }
    evt_type_isLoaded = true;
  }
  return evt_type_;
}

const int &WWWTree::mass_gluino() {
  if (not mass_gluino_isLoaded) {
    if (mass_gluino_branch != 0) {
      mass_gluino_branch->GetEntry(index);
    } else {
      printf("branch mass_gluino_branch does not exist!\n");
      exit(1);
    }
    mass_gluino_isLoaded = true;
  }
  return mass_gluino_;
}

const int &WWWTree::mass_LSP() {
  if (not mass_LSP_isLoaded) {
    if (mass_LSP_branch != 0) {
      mass_LSP_branch->GetEntry(index);
    } else {
      printf("branch mass_LSP_branch does not exist!\n");
      exit(1);
    }
    mass_LSP_isLoaded = true;
  }
  return mass_LSP_;
}

const int &WWWTree::mass_chi() {
  if (not mass_chi_isLoaded) {
    if (mass_chi_branch != 0) {
      mass_chi_branch->GetEntry(index);
    } else {
      printf("branch mass_chi_branch does not exist!\n");
      exit(1);
    }
    mass_chi_isLoaded = true;
  }
  return mass_chi_;
}

const float &WWWTree::isrboost() {
  if (not isrboost_isLoaded) {
    if (isrboost_branch != 0) {
      isrboost_branch->GetEntry(index);
    } else {
      printf("branch isrboost_branch does not exist!\n");
      exit(1);
    }
    isrboost_isLoaded = true;
  }
  return isrboost_;
}

const int &WWWTree::isr_njets() {
  if (not isr_njets_isLoaded) {
    if (isr_njets_branch != 0) {
      isr_njets_branch->GetEntry(index);
    } else {
      printf("branch isr_njets_branch does not exist!\n");
      exit(1);
    }
    isr_njets_isLoaded = true;
  }
  return isr_njets_;
}

const float &WWWTree::isr_weight() {
  if (not isr_weight_isLoaded) {
    if (isr_weight_branch != 0) {
      isr_weight_branch->GetEntry(index);
    } else {
      printf("branch isr_weight_branch does not exist!\n");
      exit(1);
    }
    isr_weight_isLoaded = true;
  }
  return isr_weight_;
}

const float &WWWTree::isr_unc() {
  if (not isr_unc_isLoaded) {
    if (isr_unc_branch != 0) {
      isr_unc_branch->GetEntry(index);
    } else {
      printf("branch isr_unc_branch does not exist!\n");
      exit(1);
    }
    isr_unc_isLoaded = true;
  }
  return isr_unc_;
}

const vector<float> &WWWTree::weightsf_lepid() {
  if (not weightsf_lepid_isLoaded) {
    if (weightsf_lepid_branch != 0) {
      weightsf_lepid_branch->GetEntry(index);
    } else {
      printf("branch weightsf_lepid_branch does not exist!\n");
      exit(1);
    }
    weightsf_lepid_isLoaded = true;
  }
  return *weightsf_lepid_;
}

const vector<float> &WWWTree::weightsf_lepiso() {
  if (not weightsf_lepiso_isLoaded) {
    if (weightsf_lepiso_branch != 0) {
      weightsf_lepiso_branch->GetEntry(index);
    } else {
      printf("branch weightsf_lepiso_branch does not exist!\n");
      exit(1);
    }
    weightsf_lepiso_isLoaded = true;
  }
  return *weightsf_lepiso_;
}

const vector<float> &WWWTree::weightsf_lepip() {
  if (not weightsf_lepip_isLoaded) {
    if (weightsf_lepip_branch != 0) {
      weightsf_lepip_branch->GetEntry(index);
    } else {
      printf("branch weightsf_lepip_branch does not exist!\n");
      exit(1);
    }
    weightsf_lepip_isLoaded = true;
  }
  return *weightsf_lepip_;
}

const vector<float> &WWWTree::weightsf_lepreco() {
  if (not weightsf_lepreco_isLoaded) {
    if (weightsf_lepreco_branch != 0) {
      weightsf_lepreco_branch->GetEntry(index);
    } else {
      printf("branch weightsf_lepreco_branch does not exist!\n");
      exit(1);
    }
    weightsf_lepreco_isLoaded = true;
  }
  return *weightsf_lepreco_;
}

const vector<float> &WWWTree::weightsf_lepconv() {
  if (not weightsf_lepconv_isLoaded) {
    if (weightsf_lepconv_branch != 0) {
      weightsf_lepconv_branch->GetEntry(index);
    } else {
      printf("branch weightsf_lepconv_branch does not exist!\n");
      exit(1);
    }
    weightsf_lepconv_isLoaded = true;
  }
  return *weightsf_lepconv_;
}

const vector<float> &WWWTree::weightsf_lepid_FS() {
  if (not weightsf_lepid_FS_isLoaded) {
    if (weightsf_lepid_FS_branch != 0) {
      weightsf_lepid_FS_branch->GetEntry(index);
    } else {
      printf("branch weightsf_lepid_FS_branch does not exist!\n");
      exit(1);
    }
    weightsf_lepid_FS_isLoaded = true;
  }
  return *weightsf_lepid_FS_;
}

const vector<float> &WWWTree::weightsf_lepiso_FS() {
  if (not weightsf_lepiso_FS_isLoaded) {
    if (weightsf_lepiso_FS_branch != 0) {
      weightsf_lepiso_FS_branch->GetEntry(index);
    } else {
      printf("branch weightsf_lepiso_FS_branch does not exist!\n");
      exit(1);
    }
    weightsf_lepiso_FS_isLoaded = true;
  }
  return *weightsf_lepiso_FS_;
}

const vector<float> &WWWTree::weightsf_lepip_FS() {
  if (not weightsf_lepip_FS_isLoaded) {
    if (weightsf_lepip_FS_branch != 0) {
      weightsf_lepip_FS_branch->GetEntry(index);
    } else {
      printf("branch weightsf_lepip_FS_branch does not exist!\n");
      exit(1);
    }
    weightsf_lepip_FS_isLoaded = true;
  }
  return *weightsf_lepip_FS_;
}

const double &WWWTree::weight_rn_r1_n1() {
  if (not weight_rn_r1_n1_isLoaded) {
    if (weight_rn_r1_n1_branch != 0) {
      weight_rn_r1_n1_branch->GetEntry(index);
    } else {
      printf("branch weight_rn_r1_n1_branch does not exist!\n");
      exit(1);
    }
    weight_rn_r1_n1_isLoaded = true;
  }
  return weight_rn_r1_n1_;
}

const double &WWWTree::weight_rn_r1_n2() {
  if (not weight_rn_r1_n2_isLoaded) {
    if (weight_rn_r1_n2_branch != 0) {
      weight_rn_r1_n2_branch->GetEntry(index);
    } else {
      printf("branch weight_rn_r1_n2_branch does not exist!\n");
      exit(1);
    }
    weight_rn_r1_n2_isLoaded = true;
  }
  return weight_rn_r1_n2_;
}

const double &WWWTree::weight_rn_r1_np05() {
  if (not weight_rn_r1_np05_isLoaded) {
    if (weight_rn_r1_np05_branch != 0) {
      weight_rn_r1_np05_branch->GetEntry(index);
    } else {
      printf("branch weight_rn_r1_np05_branch does not exist!\n");
      exit(1);
    }
    weight_rn_r1_np05_isLoaded = true;
  }
  return weight_rn_r1_np05_;
}

const double &WWWTree::weight_rn_r2_n1() {
  if (not weight_rn_r2_n1_isLoaded) {
    if (weight_rn_r2_n1_branch != 0) {
      weight_rn_r2_n1_branch->GetEntry(index);
    } else {
      printf("branch weight_rn_r2_n1_branch does not exist!\n");
      exit(1);
    }
    weight_rn_r2_n1_isLoaded = true;
  }
  return weight_rn_r2_n1_;
}

const double &WWWTree::weight_rn_r2_n2() {
  if (not weight_rn_r2_n2_isLoaded) {
    if (weight_rn_r2_n2_branch != 0) {
      weight_rn_r2_n2_branch->GetEntry(index);
    } else {
      printf("branch weight_rn_r2_n2_branch does not exist!\n");
      exit(1);
    }
    weight_rn_r2_n2_isLoaded = true;
  }
  return weight_rn_r2_n2_;
}

const double &WWWTree::weight_rn_r2_n0p5() {
  if (not weight_rn_r2_n0p5_isLoaded) {
    if (weight_rn_r2_n0p5_branch != 0) {
      weight_rn_r2_n0p5_branch->GetEntry(index);
    } else {
      printf("branch weight_rn_r2_n0p5_branch does not exist!\n");
      exit(1);
    }
    weight_rn_r2_n0p5_isLoaded = true;
  }
  return weight_rn_r2_n0p5_;
}

const double &WWWTree::weight_rn_r0p5_n1() {
  if (not weight_rn_r0p5_n1_isLoaded) {
    if (weight_rn_r0p5_n1_branch != 0) {
      weight_rn_r0p5_n1_branch->GetEntry(index);
    } else {
      printf("branch weight_rn_r0p5_n1_branch does not exist!\n");
      exit(1);
    }
    weight_rn_r0p5_n1_isLoaded = true;
  }
  return weight_rn_r0p5_n1_;
}

const double &WWWTree::weight_rn_r0p5_n2() {
  if (not weight_rn_r0p5_n2_isLoaded) {
    if (weight_rn_r0p5_n2_branch != 0) {
      weight_rn_r0p5_n2_branch->GetEntry(index);
    } else {
      printf("branch weight_rn_r0p5_n2_branch does not exist!\n");
      exit(1);
    }
    weight_rn_r0p5_n2_isLoaded = true;
  }
  return weight_rn_r0p5_n2_;
}

const double &WWWTree::weight_rn_r0p5_n0p5() {
  if (not weight_rn_r0p5_n0p5_isLoaded) {
    if (weight_rn_r0p5_n0p5_branch != 0) {
      weight_rn_r0p5_n0p5_branch->GetEntry(index);
    } else {
      printf("branch weight_rn_r0p5_n0p5_branch does not exist!\n");
      exit(1);
    }
    weight_rn_r0p5_n0p5_isLoaded = true;
  }
  return weight_rn_r0p5_n0p5_;
}

const double &WWWTree::weight_pdf_up() {
  if (not weight_pdf_up_isLoaded) {
    if (weight_pdf_up_branch != 0) {
      weight_pdf_up_branch->GetEntry(index);
    } else {
      printf("branch weight_pdf_up_branch does not exist!\n");
      exit(1);
    }
    weight_pdf_up_isLoaded = true;
  }
  return weight_pdf_up_;
}

const double &WWWTree::weight_pdf_down() {
  if (not weight_pdf_down_isLoaded) {
    if (weight_pdf_down_branch != 0) {
      weight_pdf_down_branch->GetEntry(index);
    } else {
      printf("branch weight_pdf_down_branch does not exist!\n");
      exit(1);
    }
    weight_pdf_down_isLoaded = true;
  }
  return weight_pdf_down_;
}

const double &WWWTree::weight_alphas_down() {
  if (not weight_alphas_down_isLoaded) {
    if (weight_alphas_down_branch != 0) {
      weight_alphas_down_branch->GetEntry(index);
    } else {
      printf("branch weight_alphas_down_branch does not exist!\n");
      exit(1);
    }
    weight_alphas_down_isLoaded = true;
  }
  return weight_alphas_down_;
}

const double &WWWTree::weight_alphas_up() {
  if (not weight_alphas_up_isLoaded) {
    if (weight_alphas_up_branch != 0) {
      weight_alphas_up_branch->GetEntry(index);
    } else {
      printf("branch weight_alphas_up_branch does not exist!\n");
      exit(1);
    }
    weight_alphas_up_isLoaded = true;
  }
  return weight_alphas_up_;
}

const vector<int> &WWWTree::LooseElec() {
  if (not LooseElec_isLoaded) {
    if (LooseElec_branch != 0) {
      LooseElec_branch->GetEntry(index);
    } else {
      printf("branch LooseElec_branch does not exist!\n");
      exit(1);
    }
    LooseElec_isLoaded = true;
  }
  return *LooseElec_;
}

const vector<int> &WWWTree::LooseLepton() {
  if (not LooseLepton_isLoaded) {
    if (LooseLepton_branch != 0) {
      LooseLepton_branch->GetEntry(index);
    } else {
      printf("branch LooseLepton_branch does not exist!\n");
      exit(1);
    }
    LooseLepton_isLoaded = true;
  }
  return *LooseLepton_;
}

const vector<int> &WWWTree::LooseMuon() {
  if (not LooseMuon_isLoaded) {
    if (LooseMuon_branch != 0) {
      LooseMuon_branch->GetEntry(index);
    } else {
      printf("branch LooseMuon_branch does not exist!\n");
      exit(1);
    }
    LooseMuon_isLoaded = true;
  }
  return *LooseMuon_;
}

const vector<int> &WWWTree::SignalElec() {
  if (not SignalElec_isLoaded) {
    if (SignalElec_branch != 0) {
      SignalElec_branch->GetEntry(index);
    } else {
      printf("branch SignalElec_branch does not exist!\n");
      exit(1);
    }
    SignalElec_isLoaded = true;
  }
  return *SignalElec_;
}

const vector<int> &WWWTree::SignalLepton() {
  if (not SignalLepton_isLoaded) {
    if (SignalLepton_branch != 0) {
      SignalLepton_branch->GetEntry(index);
    } else {
      printf("branch SignalLepton_branch does not exist!\n");
      exit(1);
    }
    SignalLepton_isLoaded = true;
  }
  return *SignalLepton_;
}

const vector<int> &WWWTree::SignalMuon() {
  if (not SignalMuon_isLoaded) {
    if (SignalMuon_branch != 0) {
      SignalMuon_branch->GetEntry(index);
    } else {
      printf("branch SignalMuon_branch does not exist!\n");
      exit(1);
    }
    SignalMuon_isLoaded = true;
  }
  return *SignalMuon_;
}

const vector<int> &WWWTree::Tight3lElec() {
  if (not Tight3lElec_isLoaded) {
    if (Tight3lElec_branch != 0) {
      Tight3lElec_branch->GetEntry(index);
    } else {
      printf("branch Tight3lElec_branch does not exist!\n");
      exit(1);
    }
    Tight3lElec_isLoaded = true;
  }
  return *Tight3lElec_;
}

const vector<int> &WWWTree::Tight3lLepton() {
  if (not Tight3lLepton_isLoaded) {
    if (Tight3lLepton_branch != 0) {
      Tight3lLepton_branch->GetEntry(index);
    } else {
      printf("branch Tight3lLepton_branch does not exist!\n");
      exit(1);
    }
    Tight3lLepton_isLoaded = true;
  }
  return *Tight3lLepton_;
}

const vector<int> &WWWTree::TightElec() {
  if (not TightElec_isLoaded) {
    if (TightElec_branch != 0) {
      TightElec_branch->GetEntry(index);
    } else {
      printf("branch TightElec_branch does not exist!\n");
      exit(1);
    }
    TightElec_isLoaded = true;
  }
  return *TightElec_;
}

const vector<int> &WWWTree::TightLepton() {
  if (not TightLepton_isLoaded) {
    if (TightLepton_branch != 0) {
      TightLepton_branch->GetEntry(index);
    } else {
      printf("branch TightLepton_branch does not exist!\n");
      exit(1);
    }
    TightLepton_isLoaded = true;
  }
  return *TightLepton_;
}

const vector<int> &WWWTree::TightMuon() {
  if (not TightMuon_isLoaded) {
    if (TightMuon_branch != 0) {
      TightMuon_branch->GetEntry(index);
    } else {
      printf("branch TightMuon_branch does not exist!\n");
      exit(1);
    }
    TightMuon_isLoaded = true;
  }
  return *TightMuon_;
}

const vector<int> &WWWTree::VetoElec() {
  if (not VetoElec_isLoaded) {
    if (VetoElec_branch != 0) {
      VetoElec_branch->GetEntry(index);
    } else {
      printf("branch VetoElec_branch does not exist!\n");
      exit(1);
    }
    VetoElec_isLoaded = true;
  }
  return *VetoElec_;
}

const vector<int> &WWWTree::VetoLepton() {
  if (not VetoLepton_isLoaded) {
    if (VetoLepton_branch != 0) {
      VetoLepton_branch->GetEntry(index);
    } else {
      printf("branch VetoLepton_branch does not exist!\n");
      exit(1);
    }
    VetoLepton_isLoaded = true;
  }
  return *VetoLepton_;
}

const vector<int> &WWWTree::VetoMuon() {
  if (not VetoMuon_isLoaded) {
    if (VetoMuon_branch != 0) {
      VetoMuon_branch->GetEntry(index);
    } else {
      printf("branch VetoMuon_branch does not exist!\n");
      exit(1);
    }
    VetoMuon_isLoaded = true;
  }
  return *VetoMuon_;
}

const vector<int> &WWWTree::Good3LJet() {
  if (not Good3LJet_isLoaded) {
    if (Good3LJet_branch != 0) {
      Good3LJet_branch->GetEntry(index);
    } else {
      printf("branch Good3LJet_branch does not exist!\n");
      exit(1);
    }
    Good3LJet_isLoaded = true;
  }
  return *Good3LJet_;
}

const vector<int> &WWWTree::GoodSSJet() {
  if (not GoodSSJet_isLoaded) {
    if (GoodSSJet_branch != 0) {
      GoodSSJet_branch->GetEntry(index);
    } else {
      printf("branch GoodSSJet_branch does not exist!\n");
      exit(1);
    }
    GoodSSJet_isLoaded = true;
  }
  return *GoodSSJet_;
}

const vector<int> &WWWTree::GoodSSWJet() {
  if (not GoodSSWJet_isLoaded) {
    if (GoodSSWJet_branch != 0) {
      GoodSSWJet_branch->GetEntry(index);
    } else {
      printf("branch GoodSSWJet_branch does not exist!\n");
      exit(1);
    }
    GoodSSWJet_isLoaded = true;
  }
  return *GoodSSWJet_;
}

const vector<int> &WWWTree::LooseBJet() {
  if (not LooseBJet_isLoaded) {
    if (LooseBJet_branch != 0) {
      LooseBJet_branch->GetEntry(index);
    } else {
      printf("branch LooseBJet_branch does not exist!\n");
      exit(1);
    }
    LooseBJet_isLoaded = true;
  }
  return *LooseBJet_;
}


void WWWTree::progress( int nEventsTotal, int nEventsChain ){
  int period = 1000;
  if (nEventsTotal%1000 == 0) {
    // xterm magic from L. Vacavant and A. Cerri
    if (isatty(1)) {
      if ((nEventsChain - nEventsTotal) > period) {
        float frac = (float)nEventsTotal/(nEventsChain*0.01);
        printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
             "\033[0m\033[32m <---\033[0m\015", frac);
        fflush(stdout);
      }
      else {
        printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
               "\033[0m\033[32m <---\033[0m\015", 100.);
        cout << endl;
      }
    }
  }
}

namespace www {

const int &run() { return wwwbaby.run(); }
const int &lumi() { return wwwbaby.lumi(); }
const unsigned long long &evt() { return wwwbaby.evt(); }
const int &isData() { return wwwbaby.isData(); }
const bool &evt_passgoodrunlist() { return wwwbaby.evt_passgoodrunlist(); }
const float &evt_scale1fb() { return wwwbaby.evt_scale1fb(); }
const float &evt_xsec() { return wwwbaby.evt_xsec(); }
const float &evt_kfactor() { return wwwbaby.evt_kfactor(); }
const float &evt_filter() { return wwwbaby.evt_filter(); }
const int &evt_nEvts() { return wwwbaby.evt_nEvts(); }
const vector<TString> &evt_dataset() { return wwwbaby.evt_dataset(); }
const float &puWeight() { return wwwbaby.puWeight(); }
const int &nVert() { return wwwbaby.nVert(); }
const int &firstgoodvertex() { return wwwbaby.firstgoodvertex(); }
const int &nTrueInt() { return wwwbaby.nTrueInt(); }
const float &rho() { return wwwbaby.rho(); }
const float &rho25() { return wwwbaby.rho25(); }
const float &gen_ht() { return wwwbaby.gen_ht(); }
const int &nBJetTight() { return wwwbaby.nBJetTight(); }
const int &nBJetMedium() { return wwwbaby.nBJetMedium(); }
const int &nBJetLoose() { return wwwbaby.nBJetLoose(); }
const int &nBJetTight_up() { return wwwbaby.nBJetTight_up(); }
const int &nBJetMedium_up() { return wwwbaby.nBJetMedium_up(); }
const int &nBJetLoose_up() { return wwwbaby.nBJetLoose_up(); }
const int &nBJetTight_dn() { return wwwbaby.nBJetTight_dn(); }
const int &nBJetMedium_dn() { return wwwbaby.nBJetMedium_dn(); }
const int &nBJetLoose_dn() { return wwwbaby.nBJetLoose_dn(); }
const int &nJet200MuFrac50DphiMet() { return wwwbaby.nJet200MuFrac50DphiMet(); }
const int &nMuons10() { return wwwbaby.nMuons10(); }
const int &nBadMuons20() { return wwwbaby.nBadMuons20(); }
const int &nElectrons10() { return wwwbaby.nElectrons10(); }
const int &nGammas20() { return wwwbaby.nGammas20(); }
const int &nTaus20() { return wwwbaby.nTaus20(); }
const float &met_pt() { return wwwbaby.met_pt(); }
const float &met_phi() { return wwwbaby.met_phi(); }
const float &met_calo_pt() { return wwwbaby.met_calo_pt(); }
const float &met_calo_phi() { return wwwbaby.met_calo_phi(); }
const float &met_miniaod_pt() { return wwwbaby.met_miniaod_pt(); }
const float &met_miniaod_phi() { return wwwbaby.met_miniaod_phi(); }
const float &met_muegclean_pt() { return wwwbaby.met_muegclean_pt(); }
const float &met_muegclean_phi() { return wwwbaby.met_muegclean_phi(); }
const float &met_rawPt() { return wwwbaby.met_rawPt(); }
const float &met_rawPhi() { return wwwbaby.met_rawPhi(); }
const float &met_genPt() { return wwwbaby.met_genPt(); }
const float &met_genPhi() { return wwwbaby.met_genPhi(); }
const float &sumet_raw() { return wwwbaby.sumet_raw(); }
const int &Flag_ecalLaserCorrFilter() { return wwwbaby.Flag_ecalLaserCorrFilter(); }
const int &Flag_hcalLaserEventFilter() { return wwwbaby.Flag_hcalLaserEventFilter(); }
const int &Flag_trackingFailureFilter() { return wwwbaby.Flag_trackingFailureFilter(); }
const int &Flag_CSCTightHaloFilter() { return wwwbaby.Flag_CSCTightHaloFilter(); }
const int &Flag_HBHENoiseFilter() { return wwwbaby.Flag_HBHENoiseFilter(); }
const int &Flag_HBHEIsoNoiseFilter() { return wwwbaby.Flag_HBHEIsoNoiseFilter(); }
const int &Flag_CSCTightHalo2015Filter() { return wwwbaby.Flag_CSCTightHalo2015Filter(); }
const int &Flag_EcalDeadCellTriggerPrimitiveFilter() { return wwwbaby.Flag_EcalDeadCellTriggerPrimitiveFilter(); }
const int &Flag_goodVertices() { return wwwbaby.Flag_goodVertices(); }
const int &Flag_eeBadScFilter() { return wwwbaby.Flag_eeBadScFilter(); }
const int &Flag_globalTightHalo2016() { return wwwbaby.Flag_globalTightHalo2016(); }
const int &Flag_badMuonFilter() { return wwwbaby.Flag_badMuonFilter(); }
const int &Flag_badChargedCandidateFilter() { return wwwbaby.Flag_badChargedCandidateFilter(); }
const int &Flag_badMuonFilterv2() { return wwwbaby.Flag_badMuonFilterv2(); }
const int &Flag_badChargedCandidateFilterv2() { return wwwbaby.Flag_badChargedCandidateFilterv2(); }
const int &Flag_badMuons() { return wwwbaby.Flag_badMuons(); }
const int &Flag_duplicateMuons() { return wwwbaby.Flag_duplicateMuons(); }
const int &Flag_noBadMuons() { return wwwbaby.Flag_noBadMuons(); }
const int &HLT_singleEl() { return wwwbaby.HLT_singleEl(); }
const int &HLT_singleMu() { return wwwbaby.HLT_singleMu(); }
const int &HLT_singleMu_noiso() { return wwwbaby.HLT_singleMu_noiso(); }
const int &HLT_DoubleEl_noiso() { return wwwbaby.HLT_DoubleEl_noiso(); }
const int &HLT_DoubleEl() { return wwwbaby.HLT_DoubleEl(); }
const int &HLT_DoubleEl_DZ() { return wwwbaby.HLT_DoubleEl_DZ(); }
const int &HLT_DoubleEl_DZ_2() { return wwwbaby.HLT_DoubleEl_DZ_2(); }
const int &HLT_MuEG() { return wwwbaby.HLT_MuEG(); }
const int &HLT_MuEG_2() { return wwwbaby.HLT_MuEG_2(); }
const int &HLT_MuEG_noiso() { return wwwbaby.HLT_MuEG_noiso(); }
const int &HLT_MuEG_noiso_2() { return wwwbaby.HLT_MuEG_noiso_2(); }
const int &HLT_Mu8_EG17() { return wwwbaby.HLT_Mu8_EG17(); }
const int &HLT_Mu8_EG23() { return wwwbaby.HLT_Mu8_EG23(); }
const int &HLT_Mu8_EG23_DZ() { return wwwbaby.HLT_Mu8_EG23_DZ(); }
const int &HLT_Mu12_EG23_DZ() { return wwwbaby.HLT_Mu12_EG23_DZ(); }
const int &HLT_Mu17_EG12() { return wwwbaby.HLT_Mu17_EG12(); }
const int &HLT_Mu23_EG8() { return wwwbaby.HLT_Mu23_EG8(); }
const int &HLT_Mu23_EG8_DZ() { return wwwbaby.HLT_Mu23_EG8_DZ(); }
const int &HLT_Mu23_EG12() { return wwwbaby.HLT_Mu23_EG12(); }
const int &HLT_Mu23_EG12_DZ() { return wwwbaby.HLT_Mu23_EG12_DZ(); }
const int &HLT_DoubleMu_noiso() { return wwwbaby.HLT_DoubleMu_noiso(); }
const int &HLT_DoubleMu_noiso_27_8() { return wwwbaby.HLT_DoubleMu_noiso_27_8(); }
const int &HLT_DoubleMu_noiso_30_11() { return wwwbaby.HLT_DoubleMu_noiso_30_11(); }
const int &HLT_DoubleMu_noiso_40_11() { return wwwbaby.HLT_DoubleMu_noiso_40_11(); }
const int &HLT_DoubleMu() { return wwwbaby.HLT_DoubleMu(); }
const int &HLT_DoubleMu_tk() { return wwwbaby.HLT_DoubleMu_tk(); }
const int &HLT_DoubleMu_dbltk() { return wwwbaby.HLT_DoubleMu_dbltk(); }
const int &HLT_DoubleMu_nonDZ() { return wwwbaby.HLT_DoubleMu_nonDZ(); }
const int &HLT_DoubleMu_tk_nonDZ() { return wwwbaby.HLT_DoubleMu_tk_nonDZ(); }
const int &HLT_Photon22_R9Id90_HE10_IsoM() { return wwwbaby.HLT_Photon22_R9Id90_HE10_IsoM(); }
const int &HLT_Photon30_R9Id90_HE10_IsoM() { return wwwbaby.HLT_Photon30_R9Id90_HE10_IsoM(); }
const int &HLT_Photon36_R9Id90_HE10_IsoM() { return wwwbaby.HLT_Photon36_R9Id90_HE10_IsoM(); }
const int &HLT_Photon50_R9Id90_HE10_IsoM() { return wwwbaby.HLT_Photon50_R9Id90_HE10_IsoM(); }
const int &HLT_Photon75_R9Id90_HE10_IsoM() { return wwwbaby.HLT_Photon75_R9Id90_HE10_IsoM(); }
const int &HLT_Photon90_R9Id90_HE10_IsoM() { return wwwbaby.HLT_Photon90_R9Id90_HE10_IsoM(); }
const int &HLT_Photon120_R9Id90_HE10_IsoM() { return wwwbaby.HLT_Photon120_R9Id90_HE10_IsoM(); }
const int &HLT_Photon165_R9Id90_HE10_IsoM() { return wwwbaby.HLT_Photon165_R9Id90_HE10_IsoM(); }
const int &HLT_Photon165_HE10() { return wwwbaby.HLT_Photon165_HE10(); }
const int &HLT_CaloJet500_NoJetID() { return wwwbaby.HLT_CaloJet500_NoJetID(); }
const int &HLT_ECALHT800_NoJetID() { return wwwbaby.HLT_ECALHT800_NoJetID(); }
const bool &HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton() { return wwwbaby.HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton(); }
const bool &HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton() { return wwwbaby.HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton(); }
const bool &HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton() { return wwwbaby.HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton(); }
const bool &HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton() { return wwwbaby.HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton(); }
const bool &HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton() { return wwwbaby.HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton(); }
const bool &HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton() { return wwwbaby.HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton(); }
const bool &HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton() { return wwwbaby.HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton(); }
const bool &HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton() { return wwwbaby.HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton(); }
const bool &HLT_Photon165_HE10_matchedtophoton() { return wwwbaby.HLT_Photon165_HE10_matchedtophoton(); }
const float &dilmass() { return wwwbaby.dilmass(); }
const float &dilpt() { return wwwbaby.dilpt(); }
const float &dRll() { return wwwbaby.dRll(); }
const float &matched_neutralemf() { return wwwbaby.matched_neutralemf(); }
const float &matched_emf() { return wwwbaby.matched_emf(); }
const bool &elveto() { return wwwbaby.elveto(); }
const int &nlep() { return wwwbaby.nlep(); }
const int &nveto_leptons() { return wwwbaby.nveto_leptons(); }
const int &nVetoEl_relIso03EAless01() { return wwwbaby.nVetoEl_relIso03EAless01(); }
const int &nVetoEl_relIso03EAless02() { return wwwbaby.nVetoEl_relIso03EAless02(); }
const int &nVetoEl_relIso03EAless03() { return wwwbaby.nVetoEl_relIso03EAless03(); }
const int &nVetoEl_relIso03EAless04() { return wwwbaby.nVetoEl_relIso03EAless04(); }
const int &nVetoMu_relIso03EAless01() { return wwwbaby.nVetoMu_relIso03EAless01(); }
const int &nVetoMu_relIso03EAless02() { return wwwbaby.nVetoMu_relIso03EAless02(); }
const int &nVetoMu_relIso03EAless03() { return wwwbaby.nVetoMu_relIso03EAless03(); }
const int &nVetoMu_relIso03EAless04() { return wwwbaby.nVetoMu_relIso03EAless04(); }
const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &lep_p4() { return wwwbaby.lep_p4(); }
const vector<float> &lep_pt() { return wwwbaby.lep_pt(); }
const vector<float> &lep_eta() { return wwwbaby.lep_eta(); }
const vector<float> &lep_phi() { return wwwbaby.lep_phi(); }
const vector<float> &lep_mass() { return wwwbaby.lep_mass(); }
const vector<int> &lep_charge() { return wwwbaby.lep_charge(); }
const vector<bool> &lep_3ch_agree() { return wwwbaby.lep_3ch_agree(); }
const vector<bool> &lep_isFromW() { return wwwbaby.lep_isFromW(); }
const vector<bool> &lep_isFromZ() { return wwwbaby.lep_isFromZ(); }
const vector<bool> &lep_isFromB() { return wwwbaby.lep_isFromB(); }
const vector<bool> &lep_isFromC() { return wwwbaby.lep_isFromC(); }
const vector<bool> &lep_isFromL() { return wwwbaby.lep_isFromL(); }
const vector<bool> &lep_isFromLF() { return wwwbaby.lep_isFromLF(); }
const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &lep_closest_jet_p4() { return wwwbaby.lep_closest_jet_p4(); }
const vector<double> &lep_ptRatio() { return wwwbaby.lep_ptRatio(); }
const vector<int> &lep_motherIdSS() { return wwwbaby.lep_motherIdSS(); }
const vector<int> &lep_genPart_index() { return wwwbaby.lep_genPart_index(); }
const vector<double> &lep_coneCorrPt() { return wwwbaby.lep_coneCorrPt(); }
const vector<double> &lep_ptRel() { return wwwbaby.lep_ptRel(); }
const vector<double> &lep_relIso03() { return wwwbaby.lep_relIso03(); }
const vector<double> &lep_relIso03DB() { return wwwbaby.lep_relIso03DB(); }
const vector<double> &lep_relIso03EA() { return wwwbaby.lep_relIso03EA(); }
const vector<double> &lep_relIso03EAv2() { return wwwbaby.lep_relIso03EAv2(); }
const vector<double> &lep_relIso04DB() { return wwwbaby.lep_relIso04DB(); }
const vector<double> &lep_relIso04EA() { return wwwbaby.lep_relIso04EA(); }
const vector<double> &lep_relIso04EAv2() { return wwwbaby.lep_relIso04EAv2(); }
const vector<double> &lep_miniRelIsoCMS3_EA() { return wwwbaby.lep_miniRelIsoCMS3_EA(); }
const vector<double> &lep_miniRelIsoCMS3_EAv2() { return wwwbaby.lep_miniRelIsoCMS3_EAv2(); }
const vector<double> &lep_miniRelIsoCMS3_DB() { return wwwbaby.lep_miniRelIsoCMS3_DB(); }
const vector<bool> &lep_isTriggerSafenoIso_v1() { return wwwbaby.lep_isTriggerSafenoIso_v1(); }
const vector<bool> &lep_isTriggerSafenoIso_v2() { return wwwbaby.lep_isTriggerSafenoIso_v2(); }
const vector<bool> &lep_isTriggerSafe_v1() { return wwwbaby.lep_isTriggerSafe_v1(); }
const vector<bool> &lep_isTriggerSafe_v2() { return wwwbaby.lep_isTriggerSafe_v2(); }
const vector<bool> &lep_pass_VVV_cutbased_veto() { return wwwbaby.lep_pass_VVV_cutbased_veto(); }
const vector<bool> &lep_pass_VVV_cutbased_veto_noiso() { return wwwbaby.lep_pass_VVV_cutbased_veto_noiso(); }
const vector<bool> &lep_pass_VVV_cutbased_veto_noiso_noip() { return wwwbaby.lep_pass_VVV_cutbased_veto_noiso_noip(); }
const vector<bool> &lep_pass_VVV_cutbased_veto_bak1() { return wwwbaby.lep_pass_VVV_cutbased_veto_bak1(); }
const vector<bool> &lep_pass_VVV_cutbased_veto_noiso_bak1() { return wwwbaby.lep_pass_VVV_cutbased_veto_noiso_bak1(); }
const vector<bool> &lep_pass_VVV_cutbased_veto_noiso_noip_bak1() { return wwwbaby.lep_pass_VVV_cutbased_veto_noiso_noip_bak1(); }
const vector<bool> &lep_pass_VVV_cutbased_fo() { return wwwbaby.lep_pass_VVV_cutbased_fo(); }
const vector<bool> &lep_pass_VVV_cutbased_fo_noiso() { return wwwbaby.lep_pass_VVV_cutbased_fo_noiso(); }
const vector<bool> &lep_pass_VVV_cutbased_tight_noiso() { return wwwbaby.lep_pass_VVV_cutbased_tight_noiso(); }
const vector<bool> &lep_pass_VVV_cutbased_tight() { return wwwbaby.lep_pass_VVV_cutbased_tight(); }
const vector<bool> &lep_pass_VVV_MVAbased_tight_noiso() { return wwwbaby.lep_pass_VVV_MVAbased_tight_noiso(); }
const vector<bool> &lep_pass_VVV_MVAbased_tight() { return wwwbaby.lep_pass_VVV_MVAbased_tight(); }
const vector<bool> &lep_pass_VVV_baseline() { return wwwbaby.lep_pass_VVV_baseline(); }
const int &nlep_VVV_cutbased_veto() { return wwwbaby.nlep_VVV_cutbased_veto(); }
const int &nlep_VVV_cutbased_veto_noiso() { return wwwbaby.nlep_VVV_cutbased_veto_noiso(); }
const int &nlep_VVV_cutbased_veto_noiso_noip() { return wwwbaby.nlep_VVV_cutbased_veto_noiso_noip(); }
const int &nlep_VVV_cutbased_veto_bak1() { return wwwbaby.nlep_VVV_cutbased_veto_bak1(); }
const int &nlep_VVV_cutbased_veto_noiso_bak1() { return wwwbaby.nlep_VVV_cutbased_veto_noiso_bak1(); }
const int &nlep_VVV_cutbased_veto_noiso_noip_bak1() { return wwwbaby.nlep_VVV_cutbased_veto_noiso_noip_bak1(); }
const int &nlep_VVV_cutbased_fo() { return wwwbaby.nlep_VVV_cutbased_fo(); }
const int &nlep_VVV_cutbased_fo_noiso() { return wwwbaby.nlep_VVV_cutbased_fo_noiso(); }
const int &nlep_VVV_cutbased_tight_noiso() { return wwwbaby.nlep_VVV_cutbased_tight_noiso(); }
const int &nlep_VVV_cutbased_tight() { return wwwbaby.nlep_VVV_cutbased_tight(); }
const int &nlep_VVV_MVAbased_tight_noiso() { return wwwbaby.nlep_VVV_MVAbased_tight_noiso(); }
const int &nlep_VVV_MVAbased_tight() { return wwwbaby.nlep_VVV_MVAbased_tight(); }
const int &nlep_VVV_baseline() { return wwwbaby.nlep_VVV_baseline(); }
const vector<int> &lep_pdgId() { return wwwbaby.lep_pdgId(); }
const vector<int> &lep_mc_Id() { return wwwbaby.lep_mc_Id(); }
const vector<float> &lep_dxy() { return wwwbaby.lep_dxy(); }
const vector<float> &lep_ip3d() { return wwwbaby.lep_ip3d(); }
const vector<float> &lep_ip3derr() { return wwwbaby.lep_ip3derr(); }
const vector<float> &lep_etaSC() { return wwwbaby.lep_etaSC(); }
const vector<float> &lep_r9() { return wwwbaby.lep_r9(); }
const vector<float> &lep_dz() { return wwwbaby.lep_dz(); }
const vector<int> &lep_tightId() { return wwwbaby.lep_tightId(); }
const vector<int> &lep_mcMatchId() { return wwwbaby.lep_mcMatchId(); }
const vector<int> &lep_lostHits() { return wwwbaby.lep_lostHits(); }
const vector<int> &lep_convVeto() { return wwwbaby.lep_convVeto(); }
const vector<int> &lep_tightCharge() { return wwwbaby.lep_tightCharge(); }
const vector<float> &lep_MVA() { return wwwbaby.lep_MVA(); }
const vector<float> &lep_validfraction() { return wwwbaby.lep_validfraction(); }
const vector<float> &lep_pterr() { return wwwbaby.lep_pterr(); }
const vector<float> &lep_sta_pterrOpt() { return wwwbaby.lep_sta_pterrOpt(); }
const vector<float> &lep_glb_pterrOpt() { return wwwbaby.lep_glb_pterrOpt(); }
const vector<float> &lep_x2ondof() { return wwwbaby.lep_x2ondof(); }
const vector<float> &lep_sta_x2ondof() { return wwwbaby.lep_sta_x2ondof(); }
const vector<float> &lep_glb_x2ondof() { return wwwbaby.lep_glb_x2ondof(); }
const int &nisoTrack_5gev() { return wwwbaby.nisoTrack_5gev(); }
const int &nisoTrack_mt2() { return wwwbaby.nisoTrack_mt2(); }
const int &nisoTrack_mt2_cleaned_VVV_cutbased_veto() { return wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto(); }
const int &nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso() { return wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso(); }
const int &nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip() { return wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip(); }
const int &nisoTrack_mt2_cleaned_VVV_cutbased_veto_bak1() { return wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto_bak1(); }
const int &nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_bak1() { return wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_bak1(); }
const int &nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_bak1() { return wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_veto_noiso_noip_bak1(); }
const int &nisoTrack_mt2_cleaned_VVV_cutbased_fo() { return wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_fo(); }
const int &nisoTrack_mt2_cleaned_VVV_cutbased_fo_noiso() { return wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_fo_noiso(); }
const int &nisoTrack_mt2_cleaned_VVV_cutbased_tight_noiso() { return wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_tight_noiso(); }
const int &nisoTrack_mt2_cleaned_VVV_cutbased_tight() { return wwwbaby.nisoTrack_mt2_cleaned_VVV_cutbased_tight(); }
const int &nisoTrack_mt2_cleaned_VVV_MVAbased_tight_noiso() { return wwwbaby.nisoTrack_mt2_cleaned_VVV_MVAbased_tight_noiso(); }
const int &nisoTrack_mt2_cleaned_VVV_MVAbased_tight() { return wwwbaby.nisoTrack_mt2_cleaned_VVV_MVAbased_tight(); }
const int &nisoTrack_mt2_cleaned_VVV_baseline() { return wwwbaby.nisoTrack_mt2_cleaned_VVV_baseline(); }
const int &nisoTrack_stop() { return wwwbaby.nisoTrack_stop(); }
const int &nisoTrack_PFLep5_woverlaps() { return wwwbaby.nisoTrack_PFLep5_woverlaps(); }
const int &nisoTrack_PFHad10_woverlaps() { return wwwbaby.nisoTrack_PFHad10_woverlaps(); }
const int &ngamma() { return wwwbaby.ngamma(); }
const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &gamma_p4() { return wwwbaby.gamma_p4(); }
const vector<float> &gamma_pt() { return wwwbaby.gamma_pt(); }
const vector<float> &gamma_eta() { return wwwbaby.gamma_eta(); }
const vector<float> &gamma_phi() { return wwwbaby.gamma_phi(); }
const vector<float> &gamma_mass() { return wwwbaby.gamma_mass(); }
const vector<int> &gamma_mcMatchId() { return wwwbaby.gamma_mcMatchId(); }
const vector<float> &gamma_genPt() { return wwwbaby.gamma_genPt(); }
const vector<float> &gamma_genIso() { return wwwbaby.gamma_genIso(); }
const vector<float> &gamma_chHadIso() { return wwwbaby.gamma_chHadIso(); }
const vector<float> &gamma_neuHadIso() { return wwwbaby.gamma_neuHadIso(); }
const vector<float> &gamma_phIso() { return wwwbaby.gamma_phIso(); }
const vector<float> &gamma_sigmaIetaIeta() { return wwwbaby.gamma_sigmaIetaIeta(); }
const vector<float> &gamma_r9() { return wwwbaby.gamma_r9(); }
const vector<float> &gamma_hOverE() { return wwwbaby.gamma_hOverE(); }
const vector<float> &gamma_hOverE_online() { return wwwbaby.gamma_hOverE_online(); }
const vector<int> &gamma_idCutBased() { return wwwbaby.gamma_idCutBased(); }
const vector<float> &gamma_ecpfclusiso() { return wwwbaby.gamma_ecpfclusiso(); }
const vector<float> &gamma_hcpfclusiso() { return wwwbaby.gamma_hcpfclusiso(); }
const vector<float> &gamma_hollowtkiso03() { return wwwbaby.gamma_hollowtkiso03(); }
const vector<int> &gamma_genIsPromptFinalState() { return wwwbaby.gamma_genIsPromptFinalState(); }
const vector<float> &gamma_drMinParton() { return wwwbaby.gamma_drMinParton(); }
const int &ngenPart() { return wwwbaby.ngenPart(); }
const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genPart_p4() { return wwwbaby.genPart_p4(); }
const vector<float> &genPart_pt() { return wwwbaby.genPart_pt(); }
const vector<float> &genPart_eta() { return wwwbaby.genPart_eta(); }
const vector<float> &genPart_phi() { return wwwbaby.genPart_phi(); }
const vector<float> &genPart_mass() { return wwwbaby.genPart_mass(); }
const vector<int> &genPart_pdgId() { return wwwbaby.genPart_pdgId(); }
const vector<int> &genPart_status() { return wwwbaby.genPart_status(); }
const vector<float> &genPart_charge() { return wwwbaby.genPart_charge(); }
const vector<int> &genPart_motherId() { return wwwbaby.genPart_motherId(); }
const vector<int> &genPart_grandmaId() { return wwwbaby.genPart_grandmaId(); }
const vector<bool> &genPart_isp6status3() { return wwwbaby.genPart_isp6status3(); }
const int &ngenLep() { return wwwbaby.ngenLep(); }
const vector<float> &genLep_pt() { return wwwbaby.genLep_pt(); }
const vector<float> &genLep_eta() { return wwwbaby.genLep_eta(); }
const vector<float> &genLep_phi() { return wwwbaby.genLep_phi(); }
const vector<float> &genLep_mass() { return wwwbaby.genLep_mass(); }
const vector<int> &genLep_pdgId() { return wwwbaby.genLep_pdgId(); }
const vector<int> &genLep_status() { return wwwbaby.genLep_status(); }
const vector<float> &genLep_charge() { return wwwbaby.genLep_charge(); }
const vector<int> &genLep_sourceId() { return wwwbaby.genLep_sourceId(); }
const vector<bool> &genLep_isp6status3() { return wwwbaby.genLep_isp6status3(); }
const int &ngenTau() { return wwwbaby.ngenTau(); }
const vector<float> &genTau_pt() { return wwwbaby.genTau_pt(); }
const vector<float> &genTau_eta() { return wwwbaby.genTau_eta(); }
const vector<float> &genTau_phi() { return wwwbaby.genTau_phi(); }
const vector<float> &genTau_mass() { return wwwbaby.genTau_mass(); }
const vector<int> &genTau_pdgId() { return wwwbaby.genTau_pdgId(); }
const vector<int> &genTau_status() { return wwwbaby.genTau_status(); }
const vector<float> &genTau_charge() { return wwwbaby.genTau_charge(); }
const vector<int> &genTau_sourceId() { return wwwbaby.genTau_sourceId(); }
const int &ngenLepFromTau() { return wwwbaby.ngenLepFromTau(); }
const vector<float> &genLepFromTau_pt() { return wwwbaby.genLepFromTau_pt(); }
const vector<float> &genLepFromTau_eta() { return wwwbaby.genLepFromTau_eta(); }
const vector<float> &genLepFromTau_phi() { return wwwbaby.genLepFromTau_phi(); }
const vector<float> &genLepFromTau_mass() { return wwwbaby.genLepFromTau_mass(); }
const vector<int> &genLepFromTau_pdgId() { return wwwbaby.genLepFromTau_pdgId(); }
const vector<int> &genLepFromTau_status() { return wwwbaby.genLepFromTau_status(); }
const vector<float> &genLepFromTau_charge() { return wwwbaby.genLepFromTau_charge(); }
const vector<int> &genLepFromTau_sourceId() { return wwwbaby.genLepFromTau_sourceId(); }
const int &njets() { return wwwbaby.njets(); }
const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_p4() { return wwwbaby.jets_p4(); }
const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &removed_jets_p4() { return wwwbaby.removed_jets_p4(); }
const vector<float> &removed_jets_csv() { return wwwbaby.removed_jets_csv(); }
const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_medb_p4() { return wwwbaby.jets_medb_p4(); }
const int &njets_up() { return wwwbaby.njets_up(); }
const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_up_p4() { return wwwbaby.jets_up_p4(); }
const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_medb_up_p4() { return wwwbaby.jets_medb_up_p4(); }
const vector<float> &jets_csv() { return wwwbaby.jets_csv(); }
const vector<float> &jets_up_csv() { return wwwbaby.jets_up_csv(); }
const int &njets_dn() { return wwwbaby.njets_dn(); }
const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_dn_p4() { return wwwbaby.jets_dn_p4(); }
const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_medb_dn_p4() { return wwwbaby.jets_medb_dn_p4(); }
const vector<float> &jets_dn_csv() { return wwwbaby.jets_dn_csv(); }
const vector<float> &jets_muf() { return wwwbaby.jets_muf(); }
const vector<int> &jets_mcFlavour() { return wwwbaby.jets_mcFlavour(); }
const vector<int> &jets_mcHadronFlav() { return wwwbaby.jets_mcHadronFlav(); }
const int &nhighPtPFcands() { return wwwbaby.nhighPtPFcands(); }
const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &highPtPFcands_p4() { return wwwbaby.highPtPFcands_p4(); }
const vector<float> &highPtPFcands_dz() { return wwwbaby.highPtPFcands_dz(); }
const vector<int> &highPtPFcands_pdgId() { return wwwbaby.highPtPFcands_pdgId(); }
const float &ht() { return wwwbaby.ht(); }
const float &ht_up() { return wwwbaby.ht_up(); }
const float &ht_dn() { return wwwbaby.ht_dn(); }
const float &metsig_unofficial() { return wwwbaby.metsig_unofficial(); }
const float &metsig_unofficial_up() { return wwwbaby.metsig_unofficial_up(); }
const float &metsig_unofficial_dn() { return wwwbaby.metsig_unofficial_dn(); }
const float &mt_lep1() { return wwwbaby.mt_lep1(); }
const float &mt2() { return wwwbaby.mt2(); }
const float &mt2_up() { return wwwbaby.mt2_up(); }
const float &mt2_dn() { return wwwbaby.mt2_dn(); }
const float &mt2j() { return wwwbaby.mt2j(); }
const float &mt2b() { return wwwbaby.mt2b(); }
const float &mt2b_up() { return wwwbaby.mt2b_up(); }
const float &mt2b_dn() { return wwwbaby.mt2b_dn(); }
const float &mt2_genmet() { return wwwbaby.mt2_genmet(); }
const float &mt2b_genmet() { return wwwbaby.mt2b_genmet(); }
const float &mjj_mindphi() { return wwwbaby.mjj_mindphi(); }
const float &mjj() { return wwwbaby.mjj(); }
const float &mbb_csv() { return wwwbaby.mbb_csv(); }
const float &mbb_bpt() { return wwwbaby.mbb_bpt(); }
const float &dphi_jj() { return wwwbaby.dphi_jj(); }
const float &dphi_ll() { return wwwbaby.dphi_ll(); }
const float &sum_mlb() { return wwwbaby.sum_mlb(); }
const float &deta_jj() { return wwwbaby.deta_jj(); }
const float &dR_jj() { return wwwbaby.dR_jj(); }
const float &dphi_metj1() { return wwwbaby.dphi_metj1(); }
const float &dphi_metj2() { return wwwbaby.dphi_metj2(); }
const float &dphi_genmetj1() { return wwwbaby.dphi_genmetj1(); }
const float &dphi_genmetj2() { return wwwbaby.dphi_genmetj2(); }
const float &mjj_mindphi_up() { return wwwbaby.mjj_mindphi_up(); }
const float &mjj_up() { return wwwbaby.mjj_up(); }
const float &mbb_csv_up() { return wwwbaby.mbb_csv_up(); }
const float &mbb_bpt_up() { return wwwbaby.mbb_bpt_up(); }
const float &dphi_jj_up() { return wwwbaby.dphi_jj_up(); }
const float &dphi_ll_up() { return wwwbaby.dphi_ll_up(); }
const float &sum_mlb_up() { return wwwbaby.sum_mlb_up(); }
const float &deta_jj_up() { return wwwbaby.deta_jj_up(); }
const float &dR_jj_up() { return wwwbaby.dR_jj_up(); }
const float &dphi_metj1_up() { return wwwbaby.dphi_metj1_up(); }
const float &dphi_metj2_up() { return wwwbaby.dphi_metj2_up(); }
const float &mjj_mindphi_dn() { return wwwbaby.mjj_mindphi_dn(); }
const float &mjj_dn() { return wwwbaby.mjj_dn(); }
const float &mbb_csv_dn() { return wwwbaby.mbb_csv_dn(); }
const float &mbb_bpt_dn() { return wwwbaby.mbb_bpt_dn(); }
const float &dphi_jj_dn() { return wwwbaby.dphi_jj_dn(); }
const float &dphi_ll_dn() { return wwwbaby.dphi_ll_dn(); }
const float &sum_mlb_dn() { return wwwbaby.sum_mlb_dn(); }
const float &deta_jj_dn() { return wwwbaby.deta_jj_dn(); }
const float &dR_jj_dn() { return wwwbaby.dR_jj_dn(); }
const float &dphi_metj1_dn() { return wwwbaby.dphi_metj1_dn(); }
const float &dphi_metj2_dn() { return wwwbaby.dphi_metj2_dn(); }
const float &weight_btagsf() { return wwwbaby.weight_btagsf(); }
const float &weight_btagsf_heavy_UP() { return wwwbaby.weight_btagsf_heavy_UP(); }
const float &weight_btagsf_light_UP() { return wwwbaby.weight_btagsf_light_UP(); }
const float &weight_btagsf_heavy_DN() { return wwwbaby.weight_btagsf_heavy_DN(); }
const float &weight_btagsf_light_DN() { return wwwbaby.weight_btagsf_light_DN(); }
const float &chpfcands_0013_pt() { return wwwbaby.chpfcands_0013_pt(); }
const float &chpfcands_1316_pt() { return wwwbaby.chpfcands_1316_pt(); }
const float &chpfcands_1624_pt() { return wwwbaby.chpfcands_1624_pt(); }
const float &chpfcands_2430_pt() { return wwwbaby.chpfcands_2430_pt(); }
const float &chpfcands_30in_pt() { return wwwbaby.chpfcands_30in_pt(); }
const float &phpfcands_0013_pt() { return wwwbaby.phpfcands_0013_pt(); }
const float &phpfcands_1316_pt() { return wwwbaby.phpfcands_1316_pt(); }
const float &phpfcands_1624_pt() { return wwwbaby.phpfcands_1624_pt(); }
const float &phpfcands_2430_pt() { return wwwbaby.phpfcands_2430_pt(); }
const float &phpfcands_30in_pt() { return wwwbaby.phpfcands_30in_pt(); }
const float &nupfcands_0013_pt() { return wwwbaby.nupfcands_0013_pt(); }
const float &nupfcands_1316_pt() { return wwwbaby.nupfcands_1316_pt(); }
const float &nupfcands_1624_pt() { return wwwbaby.nupfcands_1624_pt(); }
const float &nupfcands_2430_pt() { return wwwbaby.nupfcands_2430_pt(); }
const float &nupfcands_30in_pt() { return wwwbaby.nupfcands_30in_pt(); }
const float &chpfcands_0013_sumet() { return wwwbaby.chpfcands_0013_sumet(); }
const float &chpfcands_1316_sumet() { return wwwbaby.chpfcands_1316_sumet(); }
const float &chpfcands_1624_sumet() { return wwwbaby.chpfcands_1624_sumet(); }
const float &chpfcands_2430_sumet() { return wwwbaby.chpfcands_2430_sumet(); }
const float &chpfcands_30in_sumet() { return wwwbaby.chpfcands_30in_sumet(); }
const float &phpfcands_0013_sumet() { return wwwbaby.phpfcands_0013_sumet(); }
const float &phpfcands_1316_sumet() { return wwwbaby.phpfcands_1316_sumet(); }
const float &phpfcands_1624_sumet() { return wwwbaby.phpfcands_1624_sumet(); }
const float &phpfcands_2430_sumet() { return wwwbaby.phpfcands_2430_sumet(); }
const float &phpfcands_30in_sumet() { return wwwbaby.phpfcands_30in_sumet(); }
const float &nupfcands_0013_sumet() { return wwwbaby.nupfcands_0013_sumet(); }
const float &nupfcands_1316_sumet() { return wwwbaby.nupfcands_1316_sumet(); }
const float &nupfcands_1624_sumet() { return wwwbaby.nupfcands_1624_sumet(); }
const float &nupfcands_2430_sumet() { return wwwbaby.nupfcands_2430_sumet(); }
const float &nupfcands_30in_sumet() { return wwwbaby.nupfcands_30in_sumet(); }
const float &chpfcands_0013_phi() { return wwwbaby.chpfcands_0013_phi(); }
const float &chpfcands_1316_phi() { return wwwbaby.chpfcands_1316_phi(); }
const float &chpfcands_1624_phi() { return wwwbaby.chpfcands_1624_phi(); }
const float &chpfcands_2430_phi() { return wwwbaby.chpfcands_2430_phi(); }
const float &chpfcands_30in_phi() { return wwwbaby.chpfcands_30in_phi(); }
const float &phpfcands_0013_phi() { return wwwbaby.phpfcands_0013_phi(); }
const float &phpfcands_1316_phi() { return wwwbaby.phpfcands_1316_phi(); }
const float &phpfcands_1624_phi() { return wwwbaby.phpfcands_1624_phi(); }
const float &phpfcands_2430_phi() { return wwwbaby.phpfcands_2430_phi(); }
const float &phpfcands_30in_phi() { return wwwbaby.phpfcands_30in_phi(); }
const float &nupfcands_0013_phi() { return wwwbaby.nupfcands_0013_phi(); }
const float &nupfcands_1316_phi() { return wwwbaby.nupfcands_1316_phi(); }
const float &nupfcands_1624_phi() { return wwwbaby.nupfcands_1624_phi(); }
const float &nupfcands_2430_phi() { return wwwbaby.nupfcands_2430_phi(); }
const float &nupfcands_30in_phi() { return wwwbaby.nupfcands_30in_phi(); }
const float &met_T1CHS_pt() { return wwwbaby.met_T1CHS_pt(); }
const float &met_T1CHS_phi() { return wwwbaby.met_T1CHS_phi(); }
const float &met_T1CHS_fromCORE_pt() { return wwwbaby.met_T1CHS_fromCORE_pt(); }
const float &met_T1CHS_fromCORE_phi() { return wwwbaby.met_T1CHS_fromCORE_phi(); }
const float &met_T1CHS_miniAOD_CORE_pt() { return wwwbaby.met_T1CHS_miniAOD_CORE_pt(); }
const float &met_T1CHS_miniAOD_CORE_phi() { return wwwbaby.met_T1CHS_miniAOD_CORE_phi(); }
const float &met_T1CHS_miniAOD_CORE_up_pt() { return wwwbaby.met_T1CHS_miniAOD_CORE_up_pt(); }
const float &met_T1CHS_miniAOD_CORE_up_phi() { return wwwbaby.met_T1CHS_miniAOD_CORE_up_phi(); }
const float &met_T1CHS_miniAOD_CORE_dn_pt() { return wwwbaby.met_T1CHS_miniAOD_CORE_dn_pt(); }
const float &met_T1CHS_miniAOD_CORE_dn_phi() { return wwwbaby.met_T1CHS_miniAOD_CORE_dn_phi(); }
const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &decayedphoton_lep1_p4() { return wwwbaby.decayedphoton_lep1_p4(); }
const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &decayedphoton_lep2_p4() { return wwwbaby.decayedphoton_lep2_p4(); }
const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &decayedphoton_bosn_p4() { return wwwbaby.decayedphoton_bosn_p4(); }
const float &decayedphoton_mt2() { return wwwbaby.decayedphoton_mt2(); }
const int &hyp_type() { return wwwbaby.hyp_type(); }
const int &evt_type() { return wwwbaby.evt_type(); }
const int &mass_gluino() { return wwwbaby.mass_gluino(); }
const int &mass_LSP() { return wwwbaby.mass_LSP(); }
const int &mass_chi() { return wwwbaby.mass_chi(); }
const float &isrboost() { return wwwbaby.isrboost(); }
const int &isr_njets() { return wwwbaby.isr_njets(); }
const float &isr_weight() { return wwwbaby.isr_weight(); }
const float &isr_unc() { return wwwbaby.isr_unc(); }
const vector<float> &weightsf_lepid() { return wwwbaby.weightsf_lepid(); }
const vector<float> &weightsf_lepiso() { return wwwbaby.weightsf_lepiso(); }
const vector<float> &weightsf_lepip() { return wwwbaby.weightsf_lepip(); }
const vector<float> &weightsf_lepreco() { return wwwbaby.weightsf_lepreco(); }
const vector<float> &weightsf_lepconv() { return wwwbaby.weightsf_lepconv(); }
const vector<float> &weightsf_lepid_FS() { return wwwbaby.weightsf_lepid_FS(); }
const vector<float> &weightsf_lepiso_FS() { return wwwbaby.weightsf_lepiso_FS(); }
const vector<float> &weightsf_lepip_FS() { return wwwbaby.weightsf_lepip_FS(); }
const double &weight_rn_r1_n1() { return wwwbaby.weight_rn_r1_n1(); }
const double &weight_rn_r1_n2() { return wwwbaby.weight_rn_r1_n2(); }
const double &weight_rn_r1_np05() { return wwwbaby.weight_rn_r1_np05(); }
const double &weight_rn_r2_n1() { return wwwbaby.weight_rn_r2_n1(); }
const double &weight_rn_r2_n2() { return wwwbaby.weight_rn_r2_n2(); }
const double &weight_rn_r2_n0p5() { return wwwbaby.weight_rn_r2_n0p5(); }
const double &weight_rn_r0p5_n1() { return wwwbaby.weight_rn_r0p5_n1(); }
const double &weight_rn_r0p5_n2() { return wwwbaby.weight_rn_r0p5_n2(); }
const double &weight_rn_r0p5_n0p5() { return wwwbaby.weight_rn_r0p5_n0p5(); }
const double &weight_pdf_up() { return wwwbaby.weight_pdf_up(); }
const double &weight_pdf_down() { return wwwbaby.weight_pdf_down(); }
const double &weight_alphas_down() { return wwwbaby.weight_alphas_down(); }
const double &weight_alphas_up() { return wwwbaby.weight_alphas_up(); }
const vector<int> &LooseElec() { return wwwbaby.LooseElec(); }
const vector<int> &LooseLepton() { return wwwbaby.LooseLepton(); }
const vector<int> &LooseMuon() { return wwwbaby.LooseMuon(); }
const vector<int> &SignalElec() { return wwwbaby.SignalElec(); }
const vector<int> &SignalLepton() { return wwwbaby.SignalLepton(); }
const vector<int> &SignalMuon() { return wwwbaby.SignalMuon(); }
const vector<int> &Tight3lElec() { return wwwbaby.Tight3lElec(); }
const vector<int> &Tight3lLepton() { return wwwbaby.Tight3lLepton(); }
const vector<int> &TightElec() { return wwwbaby.TightElec(); }
const vector<int> &TightLepton() { return wwwbaby.TightLepton(); }
const vector<int> &TightMuon() { return wwwbaby.TightMuon(); }
const vector<int> &VetoElec() { return wwwbaby.VetoElec(); }
const vector<int> &VetoLepton() { return wwwbaby.VetoLepton(); }
const vector<int> &VetoMuon() { return wwwbaby.VetoMuon(); }
const vector<int> &Good3LJet() { return wwwbaby.Good3LJet(); }
const vector<int> &GoodSSJet() { return wwwbaby.GoodSSJet(); }
const vector<int> &GoodSSWJet() { return wwwbaby.GoodSSWJet(); }
const vector<int> &LooseBJet() { return wwwbaby.LooseBJet(); }

}
