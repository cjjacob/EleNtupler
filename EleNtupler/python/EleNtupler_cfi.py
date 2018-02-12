import FWCore.ParameterSet.Config as cms

EleNtupler = cms.EDAnalyzer(
  "EleNtupler",

  isMC = cms.bool(False),

  GlobalTag = cms.string('80X_dataRun2_2016SeptRepro_v7'),

  trigFilterDeltaRCut = cms.double(0.2),

  triggerEvent      = cms.InputTag("selectedPatTrigger", "", ""),
  triggerResults    = cms.InputTag("TriggerResults", "", "HLT"),
  patTriggerResults = cms.InputTag("TriggerResults", "", "RECO"),

  rhoLabel          = cms.InputTag("fixedGridRhoFastjetAll"),
  rhoCentralLabel   = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
  pileupCollection  = cms.InputTag("slimmedAddPileupInfo"),
  VtxLabel          = cms.InputTag("offlineSlimmedPrimaryVertices"),
  VtxBSLabel        = cms.InputTag("offlinePrimaryVerticesWithBS"),
  electronSrc       = cms.InputTag("slimmedElectrons"),
  PFAllCandidates   = cms.InputTag("particleFlow"),
  packedPFCands     = cms.InputTag("packedPFCandidates"),
  genParticleSrc    = cms.InputTag("genParticles"),
  generatorLabel    = cms.InputTag("generator"),
  LHEEventLabel     = cms.InputTag("externalLHEProducer"),

  eleVetoIdMap    = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
  eleLooseIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
  eleMediumIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
  eleTightIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
  eleHLTIdMap     = cms.InputTag("egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1"),
  eleHEEPIdMap    = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),

  eleMVAValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),

  elePFClusEcalIsoProducer = cms.InputTag("electronEcalPFClusterIsolationProducer"),
  elePFClusHcalIsoProducer = cms.InputTag("electronHcalPFClusterIsolationProducer"),

)
