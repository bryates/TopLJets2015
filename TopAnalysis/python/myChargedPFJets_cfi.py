import FWCore.ParameterSet.Config as cms
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets as ak4chargedPFJets
from RecoBTag.Configuration.RecoBTag_cff import *
from RecoJets.Configuration.RecoJetAssociations_cff import *
from PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi import *


#cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#Packed_ParticleFlow_Candidates
#selected charged candidates not entering the fit of other PVs
selectedChargedPFs = cms.EDFilter("CandPtrSelector",
                                  src = cms.InputTag("packedPFCandidates"),
                                  cut = cms.string("fromPV>0 && charge!=0")
                                  )

#jet algorithm
ak4chargedPFJets.src=cms.InputTag('selectedChargedPFs')
ak4chargedPFJets.jetPtMin=cms.double(1.0)

#b-tagging
#ak4chargedPFJetsTracksAssociatorAtVertexPF=ak5JetTracksAssociatorAtVertexPF.clone( jets = cms.InputTag("ak4chargedPFJets"),
#                                                                                   tracks = cms.InputTag("unpackedTracksAndVertices"),
#                                                                                   primaryVertex = cms.InputTag("unpackedTracksAndVertices") )
#inclusiveSecondaryVertexFinderTagInfos.extSVCollection = cms.InputTag("unpackedTracksAndVertices","secondary",""),
#combinedSecondaryVertex.trackMultiplicityMin = 1

#the sequence
myChargedPFJets = cms.Sequence(selectedChargedPFs*ak4chargedPFJets) #*ak4chargedPFJetsTracksAssociatorAtVertexPF)

