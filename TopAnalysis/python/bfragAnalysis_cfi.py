import FWCore.ParameterSet.Config as cms

bfragAnalysis = cms.EDAnalyzer("FragmentationAnalyzer",
                               hadronList = cms.vint32(511,521,531,5122,413,423,431),
                               hadronUncDzb = cms.vdouble(2.8,4,0.13,1.1,0.5,0.9,0.4),
                               hadronUncDz = cms.vdouble(1.5,0.7,0.5,1.1,0.5,0.9,0.4),
                               numEntries = cms.vint32(500000)
                               )
