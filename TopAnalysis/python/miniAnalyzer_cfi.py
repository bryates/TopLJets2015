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
                          eleVetoIdMap    = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
                          eleTightIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
                          eleTightIdFullInfoMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
                          muTriggersToUse = cms.vstring('HLT_IsoMu24_v', 'HLT_IsoTkMu24_v','HLT_IsoMu22_v', 'HLT_IsoTkMu22_v',
                                                        'HLT_IsoMu22_eta2p1_v', 'HLT_IsoTkMu22_eta2p1_v',
                                                        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v'),
                          elTriggersToUse = cms.vstring('HLT_Ele27_WPTight_Gsf_v', 'HLT_Ele32_WPTight_Gsf_v',
                                                        'HLT_Ele32_eta2p1_WPTight_Gsf_v', 'HLT_Ele25_eta2p1_WPTight_Gsf_v'
                                                        'HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                                        'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
                                                        'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v', 
                                                        'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v', 'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v'
                                                        'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v', 'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v')
                          )
