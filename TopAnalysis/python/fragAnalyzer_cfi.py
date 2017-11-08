import FWCore.ParameterSet.Config as cms

fragAnalyzer = cms.EDAnalyzer("FragmentationAnalyzer",
                              genJets = cms.InputTag("slimmedGenJets"),
                              cfg = cms.FileInPath('TopLJets2015/TopAnalysis/data/era2016/bfragweights.root')
                             )
