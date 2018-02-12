// -*- C++ -*-
//
// Package:    EleNtupler/EleNtupler
// Class:      EleNtupler
// 
/**\class EleNtupler EleNtupler.cc EleNtupler/EleNtupler/plugins/EleNtupler.cc

 Description: A class to make simple ntuples with enough data to do trigger efficiency for electron triggers

 Implementation:
     []
*/
//
// Original Author:  Colin James Jacob
//         Created:  Tue, 06 Feb 2018 11:29:41 GMT
//
//

#include "EleNtupler/EleNtupler/interface/EleNtupler.h"

using edm::InputTag;
using edm::View;

EleNtupler::EleNtupler(const edm::ParameterSet& iConfig)
{
  //usesResource("TFileService");

  isMC_ = iConfig.getParameter<bool>("isMC");

  trigFilterDeltaRCut_ = iConfig.getParameter<double>("trigFilterDeltaRCut");

  vtxLabel_ = consumes<reco::VertexCollection>(iConfig.getParameter<InputTag>("VtxLabel"));
  vtxBSLabel_ = consumes<reco::VertexCollection>(iConfig.getParameter<InputTag>("VtxBSLabel"));
  trgEventLabel_            = consumes<trigger::TriggerEvent>         (iConfig.getParameter<InputTag>("triggerEvent"));
  triggerObjectsLabel_      = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerEvent"));
  trgResultsLabel_          = consumes<edm::TriggerResults>           (iConfig.getParameter<InputTag>("triggerResults"));
  trgResultsProcess_        =                                          iConfig.getParameter<InputTag>("triggerResults").process();
  puCollection_             = consumes<vector<PileupSummaryInfo> >    (iConfig.getParameter<InputTag>("pileupCollection"));
  electronCollection_       = consumes<View<pat::Electron> >          (iConfig.getParameter<InputTag>("electronSrc"));
  pfAllParticles_           = consumes<reco::PFCandidateCollection>   (iConfig.getParameter<InputTag>("PFAllCandidates"));
  pckPFCandidateCollection_ = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<InputTag>("packedPFCands"));
  pckPFCdsLabel_            = consumes<vector<pat::PackedCandidate>>  (iConfig.getParameter<InputTag>("packedPFCands"));

  generatorLabel_           = consumes<GenEventInfoProduct>           (iConfig.getParameter<InputTag>("generatorLabel"));
  lheEventLabel_            = consumes<LHEEventProduct>               (iConfig.getParameter<InputTag>("LHEEventLabel"));
  genParticlesCollection_   = consumes<vector<reco::GenParticle>>     (iConfig.getParameter<InputTag>("genParticleSrc"));

  // electron ID 
  eleVetoIdMapToken_       = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"));
  eleLooseIdMapToken_      = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"));
  eleMediumIdMapToken_     = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"));
  eleTightIdMapToken_      = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"));
  eleHLTIdMapToken_        = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHLTIdMap"));
  eleHEEPIdMapToken_       = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap"));
  eleMVAValuesMapToken_    = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("eleMVAValuesMap"));
  elePFClusEcalIsoToken_   = mayConsume<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("elePFClusEcalIsoProducer"));
  elePFClusHcalIsoToken_   = mayConsume<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("elePFClusHcalIsoProducer"));

  // global event data
  rhoLabel_        = consumes<double>(iConfig.getParameter<InputTag>("rhoLabel"));
  rhoCentralLabel_ = consumes<double>(iConfig.getParameter<InputTag>("rhoCentralLabel"));

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("EventTree", "Event data");
  hEvents_ = fs->make<TH1F>("hEvents", "total processed and skimmed events", 2, 0, 2);

  // global data
  tree_->Branch("run",             &run_);
  tree_->Branch("event",           &event_);
  tree_->Branch("lumis",           &lumis_);
  tree_->Branch("nVtx",            &nVtx_);
  tree_->Branch("nGoodVtx",        &nGoodVtx_);
  tree_->Branch("nTracksPV",       &nTracksPV_);
  tree_->Branch("vtx",             &vtx_);
  tree_->Branch("vty",             &vty_);
  tree_->Branch("vtz",             &vtz_);
  tree_->Branch("rho",             &rho_);
  tree_->Branch("rhoCentral",      &rhoCentral_);
  tree_->Branch("eleHLTs",         &HLTEle_);
  tree_->Branch("eleHLTprescales", &HLTElePrescaled_);

  // electron variables
  tree_->Branch("nEle",                    &nEle_);
  tree_->Branch("eleCharge",               &eleCharge_);
  tree_->Branch("eleEnergy",               &eleEnergy_);
  tree_->Branch("eleD0",                   &eleD0_);
  tree_->Branch("eleDz",                   &eleDz_);
  tree_->Branch("elePt",                   &elePt_);
  tree_->Branch("eleEta",                  &eleEta_);
  tree_->Branch("elePhi",                  &elePhi_);
  tree_->Branch("eleZ",                    &eleZ_);
  tree_->Branch("eleMatchedObjPt",         &eleMatchedObjPt_);
  tree_->Branch("eleMatchedObjEta",        &eleMatchedObjEta_);
  tree_->Branch("eleMatchedObjPhi",        &eleMatchedObjPhi_);
  tree_->Branch("eleMatchedObjZ",          &eleMatchedObjZ_);
  tree_->Branch("eleMatchedObjDz",         &eleMatchedObjDz_);
  tree_->Branch("eleMatchedObjDR",         &eleMatchedObjDR_);
  tree_->Branch("eleR9",                   &eleR9_);
  //tree_->Branch("eleCalibPt",              &eleCalibPt_);
  //tree_->Branch("eleCalibEn",              &eleCalibEn_);
  tree_->Branch("eleSCEta",                &eleSCEta_);
  tree_->Branch("eleSCPhi",                &eleSCPhi_);
  tree_->Branch("eleHoverE",               &eleHoverE_);
  tree_->Branch("eleSigmaIEtaIEta",        &eleSigmaIEtaIEta_);
  tree_->Branch("eleSigmaIEtaIPhi",        &eleSigmaIEtaIPhi_);
  tree_->Branch("eleSigmaIPhiIPhi",        &eleSigmaIPhiIPhi_);
  tree_->Branch("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5_);
  tree_->Branch("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5_);
  tree_->Branch("eleConvVeto",             &eleConvVeto_);
  tree_->Branch("eleMissHits",             &eleMissHits_);
  tree_->Branch("elePFChIso",              &elePFChIso_);
  tree_->Branch("elePFPhoIso",             &elePFPhoIso_);
  tree_->Branch("elePFNeuIso",             &elePFNeuIso_);
  tree_->Branch("elePFPUIso",              &elePFPUIso_);
  tree_->Branch("elePFClusEcalIso",        &elePFClusEcalIso_);
  tree_->Branch("elePFClusHcalIso",        &elePFClusHcalIso_);
  tree_->Branch("elePFMiniIso",            &elePFMiniIso_);
  tree_->Branch("eleIDMVA",                &eleIDMVA_);
  tree_->Branch("eleFiredSingleTrgs",      &eleFiredSingleTrgs_);
  tree_->Branch("eleSingleTrigNames",      &eleSingleTrigNames_);
  tree_->Branch("eleFiredDoubleTrgs",      &eleFiredDoubleTrgs_);
  tree_->Branch("eleDoubleTrigNames",      &eleDoubleTrigNames_);
  tree_->Branch("eleFiredHLTFilters",      &eleFiredHLTFilters_);
  tree_->Branch("eleFilterNames",          &eleFilterNames_);
  tree_->Branch("eleIDbit",                &eleIDbit_);

  eleSingleTrigNames_ = {""};
  eleDoubleTrigNames_ = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"};
  eleFilterNames_ = {"hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter","hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter","hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter","hltEle27WPTightGsfTrackIsoFilter"};

}


EleNtupler::~EleNtupler()
{

}


//
// member functions
//

double EleNtupler::dDeltaPhi(const double& phi1, const double& phi2)
{
  double dPhi = phi1 - phi2;
  if ( dPhi > TMath::Pi() ) { dPhi -= 2.*TMath::Pi(); }
  if ( dPhi < -TMath::Pi() ) { dPhi += 2.*TMath::Pi(); }
  return dPhi;
}

double EleNtupler::dDeltaR(const double& eta1, const double& phi1, const double& eta2, const double& phi2)
{
  double dPhi = dDeltaPhi(phi1, phi2);
  double dEta = eta1 - eta2;
  return sqrt(dEta*dEta + dPhi*dPhi);
}

// ------------ method called for each event  ------------
void
EleNtupler::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  hEvents_->Fill(0.5); // opened event

  nEle_ = 0;
  eleCharge_.clear();
  eleEnergy_.clear();
  eleD0_.clear();
  eleDz_.clear();
  elePt_.clear();
  eleEta_.clear();
  elePhi_.clear();
  eleZ_.clear();
  eleMatchedObjPt_.clear();
  eleMatchedObjEta_.clear();
  eleMatchedObjPhi_.clear();
  eleMatchedObjZ_.clear();
  eleMatchedObjDz_.clear();
  eleMatchedObjDR_.clear();
  eleR9_.clear();
  //eleCalibPt_.clear();
  //eleCalibEn_.clear();
  eleSCEta_.clear();
  eleSCPhi_.clear();
  eleHoverE_.clear();
  eleSigmaIEtaIEta_.clear();
  eleSigmaIEtaIPhi_.clear();
  eleSigmaIPhiIPhi_.clear();
  eleSigmaIEtaIEtaFull5x5_.clear();
  eleSigmaIPhiIPhiFull5x5_.clear();
  eleConvVeto_.clear();
  eleMissHits_.clear();
  elePFChIso_.clear();
  elePFPhoIso_.clear();
  elePFNeuIso_.clear();
  elePFPUIso_.clear();
  elePFClusEcalIso_.clear();
  elePFClusHcalIso_.clear();
  elePFMiniIso_.clear();
  eleIDMVA_.clear();
  eleFiredSingleTrgs_.clear();
  eleFiredDoubleTrgs_.clear();
  eleFiredHLTFilters_.clear();
  eleIDbit_.clear();

  edm::Handle<edm::View<pat::Electron> > electronHandle;
  e.getByToken(electronCollection_, electronHandle);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjHandle;
  e.getByToken(triggerObjectsLabel_, triggerObjHandle);

  edm::Handle<edm::TriggerResults> triggerResultsHandle;
  e.getByToken(trgResultsLabel_, triggerResultsHandle);

  bool cfg_changed = true;
  HLTConfigProvider hltCfg;
  hltCfg.init(e.getRun(), es, trgResultsProcess_, cfg_changed);

  const edm::TriggerNames& names = e.triggerNames(*triggerResultsHandle);

  HLTEle_          = 0U;
  HLTElePrescaled_ = 0U;

  for ( size_t iName = 0; iName < names.size(); ++iName ) {
    const string& name = names.triggerName(iName);
    
    int bitEle = -1;
    if      (name.find("HLT_Ele25_eta2p1_WPTight_Gsf_v")                      != string::npos) bitEle =  0;
    else if (name.find("HLT_Ele27_eta2p1_WPTight_Gsf_v")                      != string::npos) bitEle =  1; 
    else if (name.find("HLT_Ele27_eta2p1_WPLoose_Gsf_v")                      != string::npos) bitEle =  2;
    else if (name.find("HLT_Ele32_eta2p1_WPTight_Gsf_v")                      != string::npos) bitEle =  3; 
    else if (name.find("HLT_Ele27_WPTight_Gsf_v")                             != string::npos) bitEle =  4; 
    else if (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")         != string::npos) bitEle =  5; 
    else if (name.find("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v")             != string::npos) bitEle =  6; 
    else if (name.find("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v")           != string::npos) bitEle =  7; 
    else if (name.find("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v")             != string::npos) bitEle =  8;
    else if (name.find("HLT_DoubleEle33_CaloIdL_MW_v")                        != string::npos) bitEle =  9;
    else if (name.find("HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL_v")          != string::npos) bitEle = 10;
    else if (name.find("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v")             != string::npos) bitEle = 11;
    else if (name.find("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v")      != string::npos) bitEle = 12;
    else if (name.find("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30_v")      != string::npos) bitEle = 13;
    else if (name.find("HLT_Ele17_Ele12_CaloId_TrackId_Iso_DZ_v")             != string::npos) bitEle = 14;
    else if (name.find("HLT_DoubleEle33_CaloId_GsfTrackIdVL_v")               != string::npos) bitEle = 15;
    else if (name.find("HLT_Ele27_WPTight_Gsf_L1JetTauSeeded_v")              != string::npos) bitEle = 16; 
    else if (name.find("HLT_Ele30_WPTight_Gsf_L1JetTauSeeded_v")              != string::npos) bitEle = 17; 
    else if (name.find("HLT_Ele32_WPTight_Gsf_L1JetTauSeeded_v")              != string::npos) bitEle = 18; 
    else if (name.find("HLT_Ele115_CaloIdVT_GsfTrkIdT_v")                     != string::npos) bitEle = 19;
    else if (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded_v") != string::npos) bitEle = 20;
    else if (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v")            != string::npos) bitEle = 21;

    UInt_t isFired     = (triggerResultsHandle->accept(iName)) ? 1 : 0;
    UInt_t isPrescaled = (hltCfg.prescaleValue(0, name) != 1)  ? 1 : 0;

    if ( bitEle >= 0 ) {
      HLTEle_          |= ( isFired << bitEle );
      HLTElePrescaled_ |= ( isPrescaled << bitEle );
    }
  }

  //edm::Handle<edm::View<pat::Electron> > calibelectronHandle;
  //e.getByToken(calibelectronCollection_, calibelectronHandle);

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  e.getByToken(pckPFCandidateCollection_, pfcands);

  if (!electronHandle.isValid()) {
    edm::LogWarning("EleNtupler") << "no electrons in event";
    return;
  }

  edm::Handle<edm::ValueMap<bool> >  veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  tight_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  hlt_id_decisions; 
  edm::Handle<edm::ValueMap<bool> >  heep_id_decisions;
  edm::Handle<edm::ValueMap<float> > eleMVAValues;
  edm::Handle<edm::ValueMap<float> > elePFClusEcalIsoValues;
  edm::Handle<edm::ValueMap<float> > elePFClusHcalIsoValues;

  e.getByToken(eleVetoIdMapToken_ ,       veto_id_decisions);
  e.getByToken(eleLooseIdMapToken_ ,      loose_id_decisions);
  e.getByToken(eleMediumIdMapToken_,      medium_id_decisions);
  e.getByToken(eleTightIdMapToken_,       tight_id_decisions);
  e.getByToken(eleHLTIdMapToken_,         hlt_id_decisions);
  e.getByToken(eleHEEPIdMapToken_ ,       heep_id_decisions);
  e.getByToken(eleMVAValuesMapToken_,     eleMVAValues);
  e.getByToken(elePFClusEcalIsoToken_,    elePFClusEcalIsoValues);
  e.getByToken(elePFClusHcalIsoToken_,    elePFClusHcalIsoValues);

  edm::Handle<reco::VertexCollection> recVtxs;
  e.getByToken(vtxLabel_, recVtxs);

  reco::Vertex vtx;
  math::XYZPoint pv(0, 0, 0);

  nVtx_     = -1;
  nGoodVtx_ = -1;
  if ( recVtxs.isValid() ) {
    nVtx_     = 0;
    nGoodVtx_ = 0;

    // best-known primary vertex coordinates 
    for (vector<reco::Vertex>::const_iterator v = recVtxs->begin(); v != recVtxs->end(); ++v) {

      bool isFake = (v->chi2() == 0 && v->ndof() == 0);

      if ( nVtx_ == 0 ) {
	nTracksPV_ = v->nTracks();
	vtx_ = v->x();
	vty_ = v->y();
	vtz_ = v->z();
      }

      if (!isFake) {
	pv.SetXYZ(v->x(), v->y(), v->z());
	vtx = *v;
	break;
      }

      if ( !v->isFake() && v->ndof() > 4. && fabs(v->z()) <= 24. && fabs(v->position().rho()) <= 2. ) {
	++nGoodVtx_;
      }
      ++nVtx_;

    } // vtx loop

  } // valid vtx handle
  else {
    edm::LogWarning("EleNtupler") << "Primary vertices info not available";
  }

  edm::Handle<double> rhoHandle;
  e.getByToken(rhoLabel_, rhoHandle);

  edm::Handle<double> rhoCentralHandle;
  e.getByToken(rhoCentralLabel_, rhoCentralHandle);

  run_    = e.id().run();
  event_  = e.id().event();
  lumis_  = e.luminosityBlock();
  rho_    = *(rhoHandle.product());
  if ( rhoCentralHandle.isValid() ) { rhoCentral_ = *(rhoCentralHandle.product()); }
  else { rhoCentral_ = -99.; }

  for ( edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle ) {

    eleCharge_       .push_back(iEle->charge());
    eleEnergy_       .push_back(iEle->energy());
    eleD0_           .push_back(iEle->gsfTrack()->dxy(pv));
    eleDz_           .push_back(iEle->gsfTrack()->dz(pv));
    elePt_           .push_back(iEle->pt());
    eleEta_          .push_back(iEle->eta());
    elePhi_          .push_back(iEle->phi());
    eleZ_            .push_back(iEle->vz());
    eleR9_           .push_back(iEle->r9());
    eleSCEta_        .push_back(iEle->superCluster()->eta());
    eleSCPhi_        .push_back(iEle->superCluster()->phi());
    eleHoverE_       .push_back(iEle->hcalOverEcal());
    eleSigmaIEtaIEta_.push_back(iEle->sigmaIetaIeta());
    eleSigmaIEtaIPhi_.push_back(iEle->sigmaIetaIphi());
    eleSigmaIPhiIPhi_.push_back(iEle->sigmaIphiIphi());
    eleConvVeto_     .push_back((Int_t)iEle->passConversionVeto());
    eleMissHits_     .push_back(iEle->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS));

    reco::GsfElectron::PflowIsolationVariables pfIso = iEle->pfIsolationVariables();
    elePFChIso_ .push_back(pfIso.sumChargedHadronPt);
    elePFPhoIso_.push_back(pfIso.sumPhotonEt);
    elePFNeuIso_.push_back(pfIso.sumNeutralHadronEt);
    elePFPUIso_ .push_back(pfIso.sumPUPt);
    //elePFMiniIso_.push_back(getMiniIsolation(pfcands, dynamic_cast<const reco::Candidate*>(&(*iEle)), 0.05, 0.2, 10., false));

    elePFClusEcalIso_.push_back(iEle->ecalPFClusterIso());
    elePFClusHcalIso_.push_back(iEle->hcalPFClusterIso());

    eleSigmaIEtaIEtaFull5x5_.push_back(iEle->full5x5_sigmaIetaIeta());
    eleSigmaIPhiIPhiFull5x5_.push_back(iEle->full5x5_sigmaIphiIphi());

    const auto el = electronHandle->ptrAt(nEle_);

    eleIDbit_.push_back(0U);

    bool passVeto = (*veto_id_decisions)[el];
    bool passLoose = (*loose_id_decisions)[el];
    bool passMedium = (*medium_id_decisions)[el];
    bool passTight = (*tight_id_decisions)[el];
    bool passHEEP = (*heep_id_decisions)[el];
    bool passHLT = (*hlt_id_decisions)[el];

    if ( passVeto ) { eleIDbit_.at(nEle_) |= ( 0b1 << 0); }
    if ( passLoose ) { eleIDbit_.at(nEle_) |= ( 0b1 << 1 ); }
    if ( passMedium ) { eleIDbit_.at(nEle_) |= ( 0b1 << 2 ); }
    if ( passTight ) { eleIDbit_.at(nEle_) |= ( 0b1 << 3 ); }
    if ( passHEEP ) { eleIDbit_.at(nEle_) |= ( 0b1 << 4 ); }
    if ( passHLT ) { eleIDbit_.at(nEle_) |= ( 0b1 << 5 ); }

    eleIDMVA_.push_back((*eleMVAValues)[el]);

    eleFiredHLTFilters_.push_back(0U);

    // in case there is no matched obj, we'll fill the data with zeros
    // then, if there is a matched obj, we'll overwrite the zeros at position nEle_
    eleMatchedObjPt_ .push_back(0.);
    eleMatchedObjEta_.push_back(0.);
    eleMatchedObjPhi_.push_back(0.);
    eleMatchedObjZ_  .push_back(0.);
    eleMatchedObjDz_ .push_back(0.);
    eleMatchedObjDR_ .push_back(0.);

    for ( pat::TriggerObjectStandAlone obj : *triggerObjHandle ) {
      obj.unpackPathNames(names);

      vector<bool> hasFilters(eleFilterNames_.size(),false);
      for ( string iFilter : obj.filterLabels() ) {
	auto it = std::find(eleFilterNames_.begin(), eleFilterNames_.end(), iFilter);
	if ( it != eleFilterNames_.end() ) {
	  hasFilters.at(it - eleFilterNames_.begin()) = true;
	}
      } // loop on filters

      if ( std::any_of(hasFilters.begin(), hasFilters.end(), [](bool b){return b;}) ) {
	double dR = dDeltaR(iEle->eta(), iEle->phi(), obj.eta(), obj.phi());
	if ( dR < trigFilterDeltaRCut_ ) {
	  // as described above, we want the matched obj to be aligned with the electron it matches
	  // we could be more space efficient with an additional variable, vector<size_t> eleMatchedObjWhichEle_
	  // then, however, we'd have to loop over the objects rather than the electrons
	  eleMatchedObjPt_ .at(nEle_) = obj.pt();
	  eleMatchedObjEta_.at(nEle_) = obj.eta();
	  eleMatchedObjPhi_.at(nEle_) = obj.phi();
	  eleMatchedObjZ_  .at(nEle_) = obj.vz();
	  //eleMatchedObjDz_ .at(nEle_) = obj.bestTrack()->dz(pv);
	  eleMatchedObjDR_ .at(nEle_) = dR;
	  auto track = obj.bestTrack();
	  if ( track ) {
	    eleMatchedObjDz_.at(nEle_) = track->dz(pv);
	  }
	  else {
	    eleMatchedObjDz_.at(nEle_) = pv.Z() - obj.vz();
	  }
	  for ( size_t i = 0; i < hasFilters.size(); ++i ) {
	    if ( hasFilters.at(i) ) {
	      eleFiredHLTFilters_.at(nEle_) |= ( 0b1<<i );
	    }
	  } // pushing filter info into tree variable
	} // trig obj matches electron
      } // is one of the filters of interest

      /*
      // probably faster than above method. N_objFilt log(N_objFilt) + min(N_objFilt,N_filt)
      // whereas the above is N_objFilt * N_filt, but requires more memory
      vector<string> filters = obj.filterLabels();
      vector<string> ourFilts = eleFilterNames_;
      vector<string> intersect;
      std::sort(filters.begin(), filters.end());
      std::sort(eleFilterNames_.begin(), eleFilterNames_.end());
      auto it = std::set_intersect(filters.begin(), filters.end(), eleFilterNames_.begin(), eleFilterNames_.end(), intersect.begin());
      if ( it != intersect.begin() ) {
        for ( string iF : intersect ) {
	  for ( size_t iFN = 0; iFN < eleFilterNames_.size(); ++iFN ) {
	    if ( iF == eleFilterNames_.at(iFN) ) {
  	      eleFiredHLTFilters_.at(nEle_) |= ( 0b1 << iFN );
	    }
	  }
	}
      }
      */

    } // loop on trigger objs

    ++nEle_;

  } // loop on electrons

  tree_->Fill();
  hEvents_->Fill(1.5); // processed event with electrons
}


// ------------ method called once each job just before starting event loop  ------------
//void 
//EleNtupler::beginJob()
//{
//}

// ------------ method called once each job just after ending the event loop  ------------
//void 
//EleNtupler::endJob() 
//{
//}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
/*
void
EleNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
*/

//define this as a plug-in
DEFINE_FWK_MODULE(EleNtupler);
