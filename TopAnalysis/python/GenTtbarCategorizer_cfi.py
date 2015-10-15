import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
from PhysicsTools.JetMCAlgos.sequences.GenHFHadronMatching_cff import *
from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons

def defineGenTtbarCategorizerSequence(process):
    
    process.ak4GenJetsCustom = ak4GenJets.clone(src = 'packedGenParticles',
                                                rParam = cms.double(0.4),
                                                jetAlgorithm = cms.string("AntiKt")
                                                )
    process.genJetFlavourPlusLeptonInfos = genJetFlavourPlusLeptonInfos.clone( jets = 'ak4GenJetsCustom',
                                                                               rParam = cms.double(0.4),
                                                                               jetAlgorithm = cms.string("AntiKt")
                                                                               )
    process.matchGenBHadron = matchGenBHadron.clone(genParticles = 'prunedGenParticles')
    process.matchGenCHadron = matchGenCHadron.clone(genParticles = 'prunedGenParticles')
    process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(particles = 'prunedGenParticles')
    process.categorizeGenTtbar = cms.EDProducer("GenTtbarCategorizer",
                                                # Phase space of additional jets
                                                genJetPtMin = cms.double(20.),
                                                genJetAbsEtaMax = cms.double(2.4),
                                                # Input tags holding information about b/c hadron matching
                                                genJets = cms.InputTag("ak4GenJetsCustom"),
                                                genBHadJetIndex = cms.InputTag("matchGenBHadron", "genBHadJetIndex"),
                                                genBHadFlavour = cms.InputTag("matchGenBHadron", "genBHadFlavour"),
                                                genBHadFromTopWeakDecay = cms.InputTag("matchGenBHadron", "genBHadFromTopWeakDecay"),
                                                genBHadPlusMothers = cms.InputTag("matchGenBHadron", "genBHadPlusMothers"),
                                                genBHadPlusMothersIndices = cms.InputTag("matchGenBHadron", "genBHadPlusMothersIndices"),
                                                genBHadIndex = cms.InputTag("matchGenBHadron", "genBHadIndex"),
                                                genBHadLeptonHadronIndex = cms.InputTag("matchGenBHadron", "genBHadLeptonHadronIndex"),
                                                genBHadLeptonViaTau = cms.InputTag("matchGenBHadron", "genBHadLeptonViaTau"),
                                                genCHadJetIndex = cms.InputTag("matchGenCHadron", "genCHadJetIndex"),
                                                genCHadFlavour = cms.InputTag("matchGenCHadron", "genCHadFlavour"),
                                                genCHadFromTopWeakDecay = cms.InputTag("matchGenCHadron", "genCHadFromTopWeakDecay"),
                                                genCHadBHadronId = cms.InputTag("matchGenCHadron", "genCHadBHadronId"),
                                                )
    process.genTtbarCategorizerSequence=cms.Sequence(process.ak4GenJetsCustom
                                                     *process.genJetFlavourPlusLeptonInfos
                                                     *process.matchGenBHadron
                                                     *process.matchGenCHadron
                                                     *process.categorizeGenTtbar)



