import FWCore.ParameterSet.Config as cms

fragWgtProducer = cms.EDAnalyzer("FragmentationWeightProducer",
                                 genJets = cms.InputTag("pseudoTop:jets"),
                                 cfg = cms.FileInPath('TopLJets2015/TopAnalysis/data/era2016/bfragweights.root'),
                                 fragModel = cms.string("PetersonModel")
                                )
