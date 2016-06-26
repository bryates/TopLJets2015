import FWCore.ParameterSet.Config as cms

analysis = cms.EDAnalyzer("MiniAnalyzer",
                          saveTree        = cms.bool(True),
                          savePF          = cms.bool(True),
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
                          muTriggersToUse = cms.vstring('HLT_IsoMu20_v', 'HLT_IsoTkMu20_v', 'HLT_IsoMu22_v', 'HLT_IsoTkMu22_v',
                                                        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v'),
                          elTriggersToUse = cms.vstring('HLT_Ele27_WPTight_Gsf_v', 'HLT_Ele25_WPTight_Gsf_v',
                                                        'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                                        'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v', 'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v')                          
                          )
