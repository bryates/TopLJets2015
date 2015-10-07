import FWCore.ParameterSet.Config as cms

analysis = cms.EDAnalyzer("MiniAnalyzer",
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
                          eleMediumIdFullInfoMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
                          muTriggersToUse = cms.vstring('IsoMu20_eta2p1_v','IsoMu20_eta2p1_TriCentralPFJet30_v','IsoMu20_eta2p1_TriCentralPFJet50_40_30_v''IsoMu24_eta2p1_v','IsoMu24_eta2p1_TriCentralPFJet30_v','IsoMu24_eta2p1_TriCentralPFJet50_40_30_v'),
                          elTriggersToUse = cms.vstring('Ele22_eta2p1_WP75_Gsf_v', 'Ele22_eta2p1_WP75_Gsf_TriCentralPFJet30_v','Ele27_eta2p1_WP75_Gsf_v', 'Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v')
                          )

