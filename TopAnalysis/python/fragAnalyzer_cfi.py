import FWCore.ParameterSet.Config as cms

fragAnalyzer = cms.EDAnalyzer("FragmentationAnalyzer",
                              genJets = cms.InputTag("slimmedGenJets")
                             )
