import FWCore.ParameterSet.Config as cms

analysis = cms.EDAnalyzer("MiniAnalyzer",
                          saveTree        = cms.bool(True),
                          triggerBits     = cms.InputTag("TriggerResults","","HLT"),
                          prescales       = cms.InputTag("patTrigger"),
                          rho             = cms.InputTag("fixedGridRhoFastjetAll"),
                          vertices        = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          muons           = cms.InputTag("slimmedMuons"),
                          electrons       = cms.InputTag("slimmedElectrons"),
                          jets            = cms.InputTag("slimmedJetsReapplyJEC"),
                          mets            = cms.InputTag("slimmedMETs"),                          
                          puppimets       = cms.InputTag('slimmedMETsPuppi'),
                          pfCands         = cms.InputTag('packedPFCandidates'),                          
                          eleVetoIdMap    = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
                          eleTightIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                          eleTightIdFullInfoMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                          muTriggersToUse = cms.vstring('HLT_IsoMu24_v',           'HLT_IsoMu22_eta2p1_v'),
                          elTriggersToUse = cms.vstring('HLT_Ele23_WPLoose_Gsf_v', 'HLT_Ele25_WPTight_Gsf_v')
                          )
