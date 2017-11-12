import FWCore.ParameterSet.Config as cms

process = cms.Process('GEN')

#import sys
#model=sys.argv[2]
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
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

# pseudo-top
process.load('TopQuarkAnalysis.TopEventProducers.producers.pseudoTop_cfi')
#process.pseudoTop.leptonMinPt=cms.double(20)
#process.pseudoTop.leptonMaxEta=cms.double(2.5)
#process.pseudoTop.jetMaxEta=cms.double(5.0)

process.load('TopLJets2015.TopAnalysis.fragAnalyzer_cfi')
print options.fragModel
process.fragAnalyzer.fragModel=options.fragModel

#process.TFileService = cms.Service("TFileService", fileName = cms.string("FragmentationDist.root"))
#process.TFileService = cms.Service("TFileService", filename = cms.string("FragmentationDist_%s.root" % options.fragModel))
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
                            fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/36CDAE89-B3BE-E611-B022-0025905B8604.root'),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )

process.analysis = cms.Path(process.pseudoTop*process.fragAnalyzer)

#process.option = cms.untracked.PSet()

#from UserCode.TopMassSecVtx.MarkusSherpaSamples_cfi import getMarkusSherpaSamplesFor
#process.source = cms.Source("PoolSource",
                            #filename = getMarkusSherpaSamplesFor(model),
                            #duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            #)
