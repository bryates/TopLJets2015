import FWCore.ParameterSet.Config as cms

bfragWgtProducer = cms.EDProducer('BFragmentationProducer',
                                  cfg = cms.FileInPath('TopLJets2015/TopAnalysis/data/era2016/bfragweights.root')
                                  )
