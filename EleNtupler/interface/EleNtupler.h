#ifndef EleNtupler_h
#define EleNtupler_h

#include "TTree.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include <memory>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;
using namespace reco;

class EleNtupler : public edm::EDAnalyzer {
 public:

  explicit EleNtupler(const edm::ParameterSet&);
  ~EleNtupler();

  //static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
  TTree   *tree_;
  TH1F    *hEvents_;
  TH1F    *hPU_;
  TH1F    *hPUTrue_;

  //virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  //virtual void endJob();

  double dDeltaPhi(const double& phi1, const double& phi2);
  double dDeltaR(const double& eta1, const double& phi1, const double& eta2, const double& phi2);

  edm::EDGetTokenT<reco::VertexCollection>         vtxLabel_;
  edm::EDGetTokenT<reco::VertexCollection>         vtxBSLabel_;
  edm::EDGetTokenT<trigger::TriggerEvent>          trgEventLabel_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsLabel_;
  edm::EDGetTokenT<edm::TriggerResults>            trgResultsLabel_;
  string                                           trgResultsProcess_;
  edm::EDGetTokenT<vector<PileupSummaryInfo> >     puCollection_;
  edm::EDGetTokenT<edm::View<pat::Electron> >      electronCollection_;
  edm::EDGetTokenT<reco::PFCandidateCollection>    pfAllParticles_;
  edm::EDGetTokenT<vector<pat::PackedCandidate> >  pckPFCdsLabel_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pckPFCandidateCollection_;
  edm::EDGetTokenT<GenEventInfoProduct>            generatorLabel_;
  edm::EDGetTokenT<LHEEventProduct>                lheEventLabel_;
  edm::EDGetTokenT<vector<reco::GenParticle> >     genParticlesCollection_;

  // electron ID decisions objects
  edm::EDGetTokenT<edm::ValueMap<bool> >  eleVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> >  eleLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> >  eleMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> >  eleTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> >  eleHLTIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> >  eleHEEPIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > eleMVAValuesMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > elePFClusEcalIsoToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > elePFClusHcalIsoToken_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidateCollection_;

  // global event data
  edm::EDGetTokenT<double> rhoLabel_;
  edm::EDGetTokenT<double> rhoCentralLabel_;

  bool isMC_;
  double trigFilterDeltaRCut_;

  // global event data
  Int_t    run_;
  Long64_t event_;
  Int_t    lumis_;
  Int_t    nVtx_;
  Int_t    nGoodVtx_;
  Int_t    nTracksPV_;
  float    vtx_;
  float    vty_;
  float    vtz_;
  float    rho_;
  float    rhoCentral_;
  UInt_t   HLTEle_;
  UInt_t   HLTElePrescaled_;

  // ntuple tree electron variables
  Int_t          nEle_;
  vector<int>    eleCharge_;
  vector<float>  eleEnergy_;
  vector<float>  eleD0_;
  vector<float>  eleDz_;
  vector<float>  elePt_;
  vector<float>  eleEta_;
  vector<float>  elePhi_;
  vector<float>  eleZ_;
  vector<float>  eleMatchedObjPt_;
  vector<float>  eleMatchedObjEta_;
  vector<float>  eleMatchedObjPhi_;
  vector<float>  eleMatchedObjZ_;
  vector<float>  eleMatchedObjDz_;
  vector<float>  eleMatchedObjDR_;
  vector<float>  eleR9_;
  vector<float>  eleCalibPt_;
  vector<float>  eleCalibEn_;
  vector<float>  eleSCEta_;
  vector<float>  eleSCPhi_;
  vector<float>  eleHoverE_;
  vector<float>  eleSigmaIEtaIEta_;
  vector<float>  eleSigmaIEtaIPhi_;
  vector<float>  eleSigmaIPhiIPhi_;
  vector<float>  eleSigmaIEtaIEtaFull5x5_;
  vector<float>  eleSigmaIPhiIPhiFull5x5_;
  vector<int>    eleConvVeto_;
  vector<int>    eleMissHits_;
  vector<float>  elePFChIso_;
  vector<float>  elePFPhoIso_;
  vector<float>  elePFNeuIso_;
  vector<float>  elePFPUIso_;
  vector<float>  elePFClusEcalIso_;
  vector<float>  elePFClusHcalIso_;
  vector<float>  elePFMiniIso_;
  vector<float>  eleIDMVA_;
  vector<UInt_t> eleFiredSingleTrgs_;
  vector<string> eleSingleTrigNames_;
  vector<UInt_t> eleFiredHLTFilters_;
  vector<string> eleFilterNames_;
  vector<UInt_t> eleFiredDoubleTrgs_;
  vector<string> eleDoubleTrigNames_;
  vector<UShort_t> eleIDbit_;

};

#endif // EleNtupler_h
