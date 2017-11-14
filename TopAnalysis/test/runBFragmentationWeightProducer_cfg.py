import FWCore.ParameterSet.Config as cms

#process = cms.Process('GEN')
process = cms.Process("Analysis")

#import sys
#model=sys.argv[2]
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('inputFile', 
                 '/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/36CDAE89-B3BE-E611-B022-0025905B8604.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "input file to process"
                 )
options.register('fragModel', 'PetersonFrag',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Run on this model"
                 )
options.parseArguments()

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 5000

process.maxEvents = cms.untracked.PSet (
    input = cms.untracked.int32(-1)
)

process.load('TopLJets2015.TopAnalysis.bfragWgtProducer_cfi')
print "Running on: %s" % options.fragModel
process.fragWgtProducer.fragModel=options.fragModel

process.TFileService = cms.Service("TFileService", fileName = cms.string("FragmentationDist_%s.root"%options.fragModel))

#from FWCore.ParameterSet.VarParsing import VarParsing
#options = VarParsing ('python')
#options.register('inputFile', 
                 #'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/36CDAE89-B3BE-E611-B022-0025905B8604.root',
                 #VarParsing.multiplicity.singleton,
                 #VarParsing.varType.string,
                 #"input file to process"
                 #)
#options.parseArguments()
process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/36CDAE89-B3BE-E611-B022-0025905B8604.root'),
                            fileNames = cms.untracked.vstring(options.inputFile),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )

# pseudo-top
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
    inputPruned = cms.InputTag("prunedGenParticles"),
    inputPacked = cms.InputTag("packedGenParticles"),
)
process.load('GeneratorInterface.RivetInterface.genParticles2HepMC_cfi')
process.genParticles2HepMC.genParticles = cms.InputTag("mergedGenParticles")
process.genParticles2HepMC.genEventInfo = cms.InputTag("generator")
process.load('TopQuarkAnalysis.TopEventProducers.producers.pseudoTop_cfi')
#process.pseudoTop.leptonMinPt=cms.double(20)
#process.pseudoTop.leptonMaxEta=cms.double(2.5)
#process.pseudoTop.jetMaxEta=cms.double(5.0)

process.analysis = cms.Path(process.mergedGenParticles*process.genParticles2HepMC*process.pseudoTop*process.fragWgtProducer)
