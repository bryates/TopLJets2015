import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer("MiniAnalyzer",
  triggerBits     = cms.InputTag("TriggerResults","","HLT"),
  prescales       = cms.InputTag("patTrigger"),
  rho             = cms.InputTag("fixedGridRhoFastjetAll"),
  vertices        = cms.InputTag("offlineSlimmedPrimaryVertices"),
  muons           = cms.InputTag("slimmedMuons"),
  electrons       = cms.InputTag("slimmedElectrons"),
  jets            = cms.InputTag("slimmedJets"),
  mets            = cms.InputTag("slimmedMETs"),
  pfCands         = cms.InputTag("packedPFCandidates"),
  eleVetoIdMap    = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
  eleLooseIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
  eleMediumIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
  eleTightIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
  
  eleMediumMVAIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp80"),
  eleTightMVAIdMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp90"),
  mvaValuesMap      = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigValues"),
  mvaCategoriesMap  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigCategories"),
  DEBUG             = cms.untracked.bool(False)
  )

