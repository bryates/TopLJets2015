import FWCore.ParameterSet.Config as cms

bfragAnalysis = cms.EDAnalyzer("FragmentationAnalyzer",
                               hadronList = cms.vint32(511,521,531,5122),
                               numEntries = cms.vint32(1)
                               )
