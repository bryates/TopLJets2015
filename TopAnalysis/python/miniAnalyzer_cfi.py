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
                      eleVetoIdMap    = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                      eleTightIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
                      eleMediumIdFullInfoMap = cms.InputTag(),
                      muTriggersToUse = cms.vstring(),
                      mlTriggersToUse = cms.vstring()
  )

