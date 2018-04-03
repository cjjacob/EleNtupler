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
  vtxBSLabel_ = consumes<reco::BeamSpot>(iConfig.getParameter<InputTag>("VtxBSLabel"));
  trgEventLabel_            = consumes<trigger::TriggerEvent>         (iConfig.getParameter<InputTag>("triggerEvent"));
  triggerObjectsLabel_      = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerEvent"));
  trgResultsLabel_          = consumes<edm::TriggerResults>           (iConfig.getParameter<InputTag>("triggerResults"));
  trgResultsProcess_        =                                          iConfig.getParameter<InputTag>("triggerResults").process();
  // puCollection_             = consumes<vector<PileupSummaryInfo> >    (iConfig.getParameter<InputTag>("pileupCollection"));
  electronCollection_       = consumes<View<pat::Electron>>           (iConfig.getParameter<InputTag>("electronSrc"));
  photonCollection_         = consumes<View<pat::Photon>>             (iConfig.getParameter<InputTag>("photonSrc"));
  superClusterCollection_   = consumes<View<reco::SuperCluster>>      (iConfig.getParameter<InputTag>("superClusterSrc"));
  // pfAllParticles_           = consumes<reco::PFCandidateCollection>   (iConfig.getParameter<InputTag>("PFAllCandidates"));
  // pckPFCandidateCollection_ = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<InputTag>("packedPFCands"));
  // pckPFCdsLabel_            = consumes<vector<pat::PackedCandidate>>  (iConfig.getParameter<InputTag>("packedPFCands"));

  // generatorLabel_           = consumes<GenEventInfoProduct>           (iConfig.getParameter<InputTag>("generatorLabel"));
  // lheEventLabel_            = consumes<LHEEventProduct>               (iConfig.getParameter<InputTag>("LHEEventLabel"));
  // genParticlesCollection_   = consumes<vector<reco::GenParticle>>     (iConfig.getParameter<InputTag>("genParticleSrc"));

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
  // rhoLabel_        = consumes<double>(iConfig.getParameter<InputTag>("rhoLabel"));
  // rhoCentralLabel_ = consumes<double>(iConfig.getParameter<InputTag>("rhoCentralLabel"));

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("EventTree", "Event data");
  hEvents_ = fs->make<TH1D>("hEvents", "total processed and skimmed events", 2, 0, 2);

  // Z data
  tree_->Branch("nZ",     &nZ_);
  tree_->Branch("ZPt",    &ZPt_);
  tree_->Branch("ZPz",    &ZPz_);
  tree_->Branch("ZEta",   &ZEta_);
  tree_->Branch("ZPhi",   &ZPhi_);
  tree_->Branch("ZM",     &ZM_);
  tree_->Branch("Ze0",    &Ze0_);
  tree_->Branch("Ze1",    &Ze1_);

  // global data
  tree_->Branch("run",             &run_);
  tree_->Branch("event",           &event_);
  tree_->Branch("lumis",           &lumis_);
  tree_->Branch("nVtx",            &nVtx_);
  // tree_->Branch("nGoodVtx",        &nGoodVtx_);
  // tree_->Branch("nTracksPV",       &nTracksPV_);
  tree_->Branch("vtx",             &vtx_);
  tree_->Branch("vty",             &vty_);
  tree_->Branch("vtz",             &vtz_);
  tree_->Branch("bsx",             &bsx_);
  tree_->Branch("bsy",             &bsy_);
  tree_->Branch("bsz",             &bsz_);
  // tree_->Branch("rho",             &rho_);
  // tree_->Branch("rhoCentral",      &rhoCentral_);
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
  tree_->Branch("eleX",                    &eleX_);
  tree_->Branch("eleY",                    &eleY_);
  tree_->Branch("eleZ",                    &eleZ_);
  tree_->Branch("gsfTrackX",               &gsfTrackX_);
  tree_->Branch("gsfTrackY",               &gsfTrackY_);
  tree_->Branch("gsfTrackZ",               &gsfTrackZ_);
  // tree_->Branch("eleR9",                   &eleR9_);
  // tree_->Branch("eleCalibPt",              &eleCalibPt_);
  // tree_->Branch("eleCalibEn",              &eleCalibEn_);
  tree_->Branch("eleSCEta",                &eleSCEta_);
  tree_->Branch("eleSCPhi",                &eleSCPhi_);
  // tree_->Branch("eleHoverE",               &eleHoverE_);
  // tree_->Branch("eleSigmaIEtaIEta",        &eleSigmaIEtaIEta_);
  // tree_->Branch("eleSigmaIEtaIPhi",        &eleSigmaIEtaIPhi_);
  // tree_->Branch("eleSigmaIPhiIPhi",        &eleSigmaIPhiIPhi_);
  // tree_->Branch("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5_);
  // tree_->Branch("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5_);
  // tree_->Branch("eleConvVeto",             &eleConvVeto_);
  tree_->Branch("eleHits",                 &eleHits_);
  tree_->Branch("eleValidPixHits",         &eleValidPixHits_);
  tree_->Branch("eleMissHits",             &eleMissHits_);
  // tree_->Branch("elePFChIso",              &elePFChIso_);
  // tree_->Branch("elePFPhoIso",             &elePFPhoIso_);
  // tree_->Branch("elePFNeuIso",             &elePFNeuIso_);
  // tree_->Branch("elePFPUIso",              &elePFPUIso_);
  // tree_->Branch("elePFClusEcalIso",        &elePFClusEcalIso_);
  // tree_->Branch("elePFClusHcalIso",        &elePFClusHcalIso_);
  // tree_->Branch("elePFMiniIso",            &elePFMiniIso_);
  tree_->Branch("eleIDMVA",                &eleIDMVA_);
  tree_->Branch("eleFiredHLTFilters",      &eleFiredHLTFilters_);
  // tree_->Branch("eleFilterNames",          &eleFilterNames_);
  tree_->Branch("eleIDbit",                &eleIDbit_);

  // TrigObj
  tree_->Branch("nTO",                 &nTO_);
  tree_->Branch("TrigObjPt",           &TrigObjPt_);
  tree_->Branch("TrigObjEta",          &TrigObjEta_);
  tree_->Branch("TrigObjPhi",          &TrigObjPhi_);
  tree_->Branch("TrigObjEnergy",       &TrigObjEnergy_);
  tree_->Branch("TrigObjDR",           &TrigObjDR_);
  tree_->Branch("TrigObjMatchedEle",   &TrigObjMatchedEle_);
  tree_->Branch("TrigObjType",         &TrigObjType_);
  tree_->Branch("TrigObjFiredFilters", &TrigObjFiredFilters_);

  // photon
  tree_->Branch("nPho",      &nPho_);
  tree_->Branch("phoPt",     &phoPt_);
  tree_->Branch("phoEta",    &phoEta_);
  tree_->Branch("phoPhi",    &phoPhi_);
  tree_->Branch("phoEnergy", &phoEnergy_);

  // supercluster
  tree_->Branch("nSC",      &nSC_);
  // tree_->Branch("scPt",     &scPt_);
  tree_->Branch("scEta",    &scEta_);
  tree_->Branch("scPhi",    &scPhi_);
  tree_->Branch("scEnergy", &scEnergy_);

  eleFilterNames_ = {"hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter","hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter","hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter","hltEle27WPTightGsfTrackIsoFilter","hltEgammaEleGsfTrackIso","hltEGL1SingleAndDoubleEGOrPairFilter","hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg1Filter","hltEle23Ele12CaloIdLTrackIdLIsoVLEtLeg2Filter"};

}


EleNtupler::~EleNtupler()
{

}


//
// member functions
//

EleNtupler::ZWithEles::ZWithEles(const pat::Electron& e0, const size_t& ei0, const pat::Electron& e1, const size_t& ei1, const vector<UShort_t>& eleIDs)
{
  // auto start = chrono::high_resolution_clock::now();

  if ( rand() % 2 ) { eles_.first = ei0; eles_.second = ei1; }
  else { eles_.first = ei1; eles_.second = ei0; }

  TLorentzVector ve0, ve1;
  ve0.SetPtEtaPhiM(e0.pt(), e0.eta(), e0.phi(), 0.000511);
  ve1.SetPtEtaPhiM(e1.pt(), e1.eta(), e1.phi(), 0.000511);
  TLorentzVector vz = ve0 + ve1;

  Z_ = vz;
  
  bool neutral = (e0.charge() * e1.charge()) < 0.0;
  bool goodMass = fabs( vz.M() - 91.19 ) <= 20.0;
  
  bool loose0, loose1;
  if ( (eleIDs.at(ei0) & 2) > 0 ) { loose0 = true; }
  if ( (eleIDs.at(ei1) & 2) > 0 ) { loose1 = true; }

  if ( neutral && goodMass && loose0 && loose1 ) {
    isGoodZ_ = true;
  } // opp charge and within 20 GeV of Z mass
  else {
    isGoodZ_ = false;
  }

  // auto end = chrono::high_resolution_clock::now();
  // cout << " |- Constructed ZWithEles in " << chrono::duration_cast<chrono::microseconds>(end-start).count() << " us\n";
}

EleNtupler::ZWithEles::~ZWithEles() {}

const std::pair<size_t,size_t>& EleNtupler::ZWithEles::Eles() const
{
  return this->eles_;
}

const double EleNtupler::ZWithEles::M() const
{
  return this->Z_.M();
}

const TLorentzVector EleNtupler::ZWithEles::Z() const
{
  return this->Z_;
}

const bool EleNtupler::ZWithEles::IsGoodZ() const
{
  return this->isGoodZ_;
}

void EleNtupler::MakeZList(const edm::View<pat::Electron>& eleView, const vector<UShort_t>& eleIDs)
{
  // auto start = chrono::high_resolution_clock::now();

  for ( size_t iEle = 0; iEle < eleView.size()-1; ++iEle ) {
    for ( size_t jEle = iEle+1; jEle < eleView.size(); ++jEle ) {
      ZWithEles z(eleView[iEle], iEle, eleView[jEle], jEle, eleIDs);
      if ( z.IsGoodZ() ) {
	// keep only good Z's
	ZList.push_back(z);
      }
    }
  }

  // auto listed = chrono::high_resolution_clock::now();
  // cout << "--- Made ZList in " << chrono::duration_cast<chrono::microseconds>(listed-start).count() << " us\n";

  // sort by closeness to Z mass
  std::sort( ZList.begin(), ZList.end(), 
	     [](const ZWithEles& a, const ZWithEles& b){ return fabs(a.M()-91.19) < fabs(b.M()-91.19); } 
	     );

  // auto sorted = chrono::high_resolution_clock::now();
  // cout << "--- Sorted ZList in " << chrono::duration_cast<chrono::microseconds>(sorted-listed).count() << " us\n";

  return;
}

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
  // cout << "Analyze\n";

  // auto start = chrono::high_resolution_clock::now();

  hEvents_->Fill(0.5); // opened event

  nEle_ = 0;
  eleCharge_.clear();
  eleEnergy_.clear();
  eleD0_.clear();
  eleDz_.clear();
  elePt_.clear();
  eleEta_.clear();
  elePhi_.clear();
  eleX_.clear();
  eleY_.clear();
  eleZ_.clear();
  gsfTrackX_.clear();
  gsfTrackY_.clear();
  gsfTrackZ_.clear();
  // eleR9_.clear();
  // eleCalibPt_.clear();
  // eleCalibEn_.clear();
  eleSCEta_.clear();
  eleSCPhi_.clear();
  // eleHoverE_.clear();
  // eleSigmaIEtaIEta_.clear();
  // eleSigmaIEtaIPhi_.clear();
  // eleSigmaIPhiIPhi_.clear();
  // eleSigmaIEtaIEtaFull5x5_.clear();
  // eleSigmaIPhiIPhiFull5x5_.clear();
  // eleConvVeto_.clear();
  eleHits_.clear();
  eleValidPixHits_.clear();
  eleMissHits_.clear();
  // elePFChIso_.clear();
  // elePFPhoIso_.clear();
  // elePFNeuIso_.clear();
  // elePFPUIso_.clear();
  // elePFClusEcalIso_.clear();
  // elePFClusHcalIso_.clear();
  // elePFMiniIso_.clear();
  eleIDMVA_.clear();
  eleFiredHLTFilters_.clear();
  eleIDbit_.clear();

  // TrigObj
  nTO_ = 0;
  TrigObjPt_.clear();
  TrigObjEta_.clear();
  TrigObjPhi_.clear();
  TrigObjEnergy_.clear();
  TrigObjDR_.clear();
  TrigObjMatchedEle_.clear();
  TrigObjType_.clear();
  TrigObjFiredFilters_.clear();

  // photon
  nPho_ = 0;
  phoPt_.clear();
  phoEta_.clear();
  phoPhi_.clear();
  phoEnergy_.clear();

  // auto cleared = chrono::high_resolution_clock::now();
  // cout << "Cleared vecs in " << chrono::duration_cast<chrono::microseconds>(cleared-start).count() << " us\n";

  edm::Handle<edm::View<pat::Electron>> electronHandle;
  e.getByToken(electronCollection_, electronHandle);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjHandle;
  e.getByToken(triggerObjectsLabel_, triggerObjHandle);

  edm::Handle<edm::TriggerResults> triggerResultsHandle;
  e.getByToken(trgResultsLabel_, triggerResultsHandle);

  edm::Handle<edm::View<pat::Photon>> photonHandle;
  e.getByToken(photonCollection_, photonHandle);

  edm::Handle<edm::View<reco::SuperCluster>> superClusterHandle;
  e.getByToken(superClusterCollection_, superClusterHandle);

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

  // auto hlts = chrono::high_resolution_clock::now();
  // cout << "Found HLTs in " << chrono::duration_cast<chrono::microseconds>(hlts-cleared).count() << " us\n";

  //edm::Handle<edm::View<pat::Electron> > calibelectronHandle;
  //e.getByToken(calibelectronCollection_, calibelectronHandle);

  // edm::Handle<pat::PackedCandidateCollection> pfcands;
  // e.getByToken(pckPFCandidateCollection_, pfcands);

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

  // auto ids = chrono::high_resolution_clock::now();
  // cout << "Found IDs in " << chrono::duration_cast<chrono::microseconds>(ids-hlts).count() << " us\n";

  edm::Handle<reco::VertexCollection> recVtxs;
  e.getByToken(vtxLabel_, recVtxs);

  edm::Handle<reco::BeamSpot> bsHandle;
  e.getByToken(vtxBSLabel_, bsHandle);
  const reco::BeamSpot& beamspot = *bsHandle.product();

  bsx_ = beamspot.position().x();
  bsy_ = beamspot.position().y();
  bsz_ = beamspot.position().z();

  reco::Vertex vtx;
  math::XYZPoint pv(0, 0, 0);

  nVtx_     = -1;
  // nGoodVtx_ = -1;
  if ( recVtxs.isValid() ) {
    nVtx_     = 0;
    // nGoodVtx_ = 0;

    // best-known primary vertex coordinates 
    for (vector<reco::Vertex>::const_iterator v = recVtxs->begin(); v != recVtxs->end(); ++v) {

      bool isFake = (v->chi2() == 0 && v->ndof() == 0);

      if ( v == recVtxs->begin() && !isFake ) {
	// nTracksPV_ = v->nTracks();
	vtx_ = v->x();
	vty_ = v->y();
	vtz_ = v->z();
	pv.SetXYZ(v->x(), v->y(), v->z());
	vtx = *v;
      }

      // if ( !v->isFake() && v->ndof() > 4. && fabs(v->z()) <= 24. && fabs(v->position().rho()) <= 2. ) {
      // 	++nGoodVtx_;
      // }
      ++nVtx_;

    } // vtx loop

  } // valid vtx handle
  else {
    edm::LogWarning("EleNtupler") << "Primary vertices info not available";
  }

  // auto vtxs = chrono::high_resolution_clock::now();
  // cout << "Found Vtxs in " << chrono::duration_cast<chrono::microseconds>(vtxs-ids).count() << " us\n";

  // edm::Handle<double> rhoHandle;
  // e.getByToken(rhoLabel_, rhoHandle);

  // edm::Handle<double> rhoCentralHandle;
  // e.getByToken(rhoCentralLabel_, rhoCentralHandle);

  run_    = e.id().run();
  event_  = e.id().event();
  lumis_  = e.luminosityBlock();
  // rho_    = *(rhoHandle.product());
  // if ( rhoCentralHandle.isValid() ) { rhoCentral_ = *(rhoCentralHandle.product()); }
  // else { rhoCentral_ = -99.; }

  for ( edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle ) {

    eleCharge_       .push_back(iEle->charge());
    eleEnergy_       .push_back(iEle->energy());
    eleD0_           .push_back(iEle->gsfTrack()->dxy(pv));
    eleDz_           .push_back(iEle->gsfTrack()->dz(pv));
    elePt_           .push_back(iEle->pt());
    eleEta_          .push_back(iEle->eta());
    elePhi_          .push_back(iEle->phi());
    eleX_            .push_back(iEle->vx());
    eleY_            .push_back(iEle->vy());
    eleZ_            .push_back(iEle->vz());
    gsfTrackX_       .push_back(iEle->gsfTrack()->vx());
    gsfTrackY_       .push_back(iEle->gsfTrack()->vy());
    gsfTrackZ_       .push_back(iEle->gsfTrack()->vz());
    // eleR9_           .push_back(iEle->r9());
    eleSCEta_        .push_back(iEle->superCluster()->eta());
    eleSCPhi_        .push_back(iEle->superCluster()->phi());
    // eleHoverE_       .push_back(iEle->hcalOverEcal());
    // eleSigmaIEtaIEta_.push_back(iEle->sigmaIetaIeta());
    // eleSigmaIEtaIPhi_.push_back(iEle->sigmaIetaIphi());
    // eleSigmaIPhiIPhi_.push_back(iEle->sigmaIphiIphi());
    // eleConvVeto_     .push_back((Int_t)iEle->passConversionVeto());
    eleHits_         .push_back(iEle->gsfTrack()->hitPattern().trackerLayersWithMeasurement());
    eleValidPixHits_ .push_back(iEle->gsfTrack()->hitPattern().numberOfValidPixelHits());
    eleMissHits_     .push_back(iEle->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS));

    // reco::GsfElectron::PflowIsolationVariables pfIso = iEle->pfIsolationVariables();
    // elePFChIso_ .push_back(pfIso.sumChargedHadronPt);
    // elePFPhoIso_.push_back(pfIso.sumPhotonEt);
    // elePFNeuIso_.push_back(pfIso.sumNeutralHadronEt);
    // elePFPUIso_ .push_back(pfIso.sumPUPt);
    // //elePFMiniIso_.push_back(getMiniIsolation(pfcands, dynamic_cast<const reco::Candidate*>(&(*iEle)), 0.05, 0.2, 10., false));

    // elePFClusEcalIso_.push_back(iEle->ecalPFClusterIso());
    // elePFClusHcalIso_.push_back(iEle->hcalPFClusterIso());

    // eleSigmaIEtaIEtaFull5x5_.push_back(iEle->full5x5_sigmaIetaIeta());
    // eleSigmaIPhiIPhiFull5x5_.push_back(iEle->full5x5_sigmaIphiIphi());

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

    //// moved trigger objects to their own loop

    // in case there is no matched obj, we'll fill the data with zeros
    // then, if there is a matched obj, we'll overwrite the zeros at position nEle_
    // eleMatchedObjPt_ .push_back(0.);
    // eleMatchedObjEta_.push_back(0.);
    // eleMatchedObjPhi_.push_back(0.);
    // eleMatchedObjDR_ .push_back(0.);

    // for ( pat::TriggerObjectStandAlone obj : *triggerObjHandle ) {
    //   obj.unpackPathNames(names);

    //   vector<bool> hasFilters(eleFilterNames_.size(),false);
    //   for ( string iFilter : obj.filterLabels() ) {
    // 	auto it = std::find(eleFilterNames_.begin(), eleFilterNames_.end(), iFilter);
    // 	if ( it != eleFilterNames_.end() ) {
    // 	  hasFilters.at(it - eleFilterNames_.begin()) = true;
    // 	}
    //   } // loop on filters

    //   if ( std::any_of(hasFilters.begin(), hasFilters.end(), [](bool b){return b;}) ) {
    // 	double dR = dDeltaR(iEle->eta(), iEle->phi(), obj.eta(), obj.phi());
    // 	if ( dR < trigFilterDeltaRCut_ ) {
    // 	  // as described above, we want the matched obj to be aligned with the electron it matches
    // 	  // we could be more space efficient with an additional variable, vector<size_t> eleMatchedObjWhichEle_
    // 	  // then, however, we'd have to loop over the objects rather than the electrons
    // 	  eleMatchedObjPt_ .at(nEle_) = obj.pt();
    // 	  eleMatchedObjEta_.at(nEle_) = obj.eta();
    // 	  eleMatchedObjPhi_.at(nEle_) = obj.phi();
    // 	  eleMatchedObjDR_ .at(nEle_) = dR;
    // 	  for ( size_t i = 0; i < hasFilters.size(); ++i ) {
    // 	    if ( hasFilters.at(i) ) {
    // 	      eleFiredHLTFilters_.at(nEle_) |= ( 0b1<<i );
    // 	    }
    // 	  } // pushing filter info into tree variable
    // 	} // trig obj matches electron
    //   } // is one of the filters of interest

    //   /*
    //   // probably faster than above method. N_objFilt log(N_objFilt) + min(N_objFilt,N_filt)
    //   // whereas the above is N_objFilt * N_filt, but requires more memory
    //   vector<string> filters = obj.filterLabels();
    //   vector<string> ourFilts = eleFilterNames_;
    //   vector<string> intersect;
    //   std::sort(filters.begin(), filters.end());
    //   std::sort(eleFilterNames_.begin(), eleFilterNames_.end());
    //   auto it = std::set_intersect(filters.begin(), filters.end(), eleFilterNames_.begin(), eleFilterNames_.end(), intersect.begin());
    //   if ( it != intersect.begin() ) {
    //     for ( string iF : intersect ) {
    // 	  for ( size_t iFN = 0; iFN < eleFilterNames_.size(); ++iFN ) {
    // 	    if ( iF == eleFilterNames_.at(iFN) ) {
    // 	      eleFiredHLTFilters_.at(nEle_) |= ( 0b1 << iFN );
    // 	    }
    // 	  }
    // 	}
    //   }
    //   */

    // } // loop on trigger objs

    ++nEle_;

  } // loop on electrons

  // auto elecs = chrono::high_resolution_clock::now();
  // cout << "Looped over eles in " << chrono::duration_cast<chrono::microseconds>(elecs-vtxs).count() << " us\n";

  ZList.clear();
  nZ_ = 0;
  ZPt_.clear();
  ZPz_.clear();
  ZEta_.clear();
  ZPhi_.clear();
  ZM_.clear();
  Ze0_.clear();
  Ze1_.clear();
  if ( nEle_ > 1 ) {
    MakeZList(*electronHandle, eleIDbit_);
    nZ_ = static_cast<int>(ZList.size());
    for ( const auto& z : ZList ) {
      ZPt_ .push_back(z.Z().Pt());
      ZPz_ .push_back(z.Z().Pz());
      ZEta_.push_back(z.Z().Eta());
      ZPhi_.push_back(z.Z().Phi());
      ZM_  .push_back(z.M());
      Ze0_ .push_back( static_cast<int>(z.Eles().first) );
      Ze1_ .push_back( static_cast<int>(z.Eles().second) );
    }
  }

  // auto zeds = chrono::high_resolution_clock::now();
  // cout << "Looped over Zs in " << chrono::duration_cast<chrono::microseconds>(zeds-elecs).count() << " us\n";

  for ( edm::View<pat::Photon>::const_iterator iPho = photonHandle->begin(); iPho != photonHandle->end(); ++iPho ) {
    phoPt_    .push_back(iPho->pt());
    phoEta_   .push_back(iPho->eta());
    phoPhi_   .push_back(iPho->phi());
    phoEnergy_.push_back(iPho->energy());
    ++nPho_;
  } // photon loop

  for ( edm::View<reco::SuperCluster>::const_iterator iSC = superClusterHandle->begin(); iSC != superClusterHandle->end(); ++iSC ) {
    // scPt_    .push_back(iSC->pt());
    scEta_   .push_back(iSC->eta());
    scPhi_   .push_back(iSC->phi());
    scEnergy_.push_back(iSC->energy());
    ++nSC_;
  } // super cluster loop



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

      TrigObjPt_    .push_back(obj.pt());
      TrigObjEta_   .push_back(obj.eta());
      TrigObjPhi_   .push_back(obj.phi());
      TrigObjEnergy_.push_back(obj.energy());

      // TriggerPhoton   = +81,
      // TriggerElectron = +82,
      // TriggerMuon     = +83,
      // TriggerTau      = +84,
      // TriggerTrack    = +91,
      // TriggerCluster  = +92,
      TrigObjType_.push_back( obj.triggerObjectTypes() );

      TrigObjFiredFilters_.push_back(0U);
      for ( size_t i = 0; i < hasFilters.size(); ++i ) {
	if ( hasFilters.at(i) ) {
	  TrigObjFiredFilters_.back() |= ( 0b1 << i );
	}
      }

      TrigObjMatchedEle_.push_back(-1);
      size_t ele_ind = 0U;
      for ( edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle ) {
	double dR = dDeltaR(iEle->eta(), iEle->phi(), obj.eta(), obj.phi());
	if ( dR < trigFilterDeltaRCut_ ) {
	  TrigObjMatchedEle_.back() = ele_ind;
	  eleFiredHLTFilters_.at(ele_ind) = TrigObjFiredFilters_.back();
	} // trig obj matches electron
	++ele_ind;
      } // loop on electrons

      ++nTO_;

    } // has one of the HLT_Ele23_Ele12_... filters
  } // trigger obj loop

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
