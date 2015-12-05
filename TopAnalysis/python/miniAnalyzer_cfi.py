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
                          pfCands         = cms.InputTag("packedPFCandidates"),
                          genTtbarId      = cms.InputTag("categorizeGenTtbar", "genTtbarId"),
                          eleVetoIdMap    = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
                          eleTightIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                          eleTightIdFullInfoMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                          muTriggersToUse = cms.vstring('IsoMu20_eta2p1_v','IsoTkMu20_v', 'IsoMu20_eta2p1_TriCentralPFJet30_v','IsoMu20_eta2p1_TriCentralPFJet50_40_30_v''IsoMu24_eta2p1_v','IsoMu24_eta2p1_TriCentralPFJet30_v','IsoMu24_eta2p1_TriCentralPFJet50_40_30_v'),
                          elTriggersToUse = cms.vstring('Ele22_eta2p1_WP75_Gsf_v', 'Ele22_eta2p1_WP75_Gsf_TriCentralPFJet30_v','Ele27_eta2p1_WP75_Gsf_v', 'Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v','Ele23_CaloIdL_TrackIdL_IsoVL_v')
                          )

