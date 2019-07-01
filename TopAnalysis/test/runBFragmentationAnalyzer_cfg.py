import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('python')
options.register('inputFile', 
                 '/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/36CDAE89-B3BE-E611-B022-0025905B8604.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "input file to process"
                 )
options.register('frag', 
		 'BL', 
		 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
		 "Fragmentation formula to use: BL (Bowler-Lund) P (Peterson)")
options.register('param', 
		 0.855, 
		 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
		 "StringZ:rFactB (BL) or StringZ:epsilonB (P)")
options.parseArguments()

#process = cms.Process('GEN')
process = cms.Process("Analysis")

#message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(options.maxEvents)
#)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.options = cms.untracked.PSet(
 wantSummary = cms.untracked.bool(True)
)

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/36CDAE89-B3BE-E611-B022-0025905B8604.root'),
                            fileNames = cms.untracked.vstring(options.inputFile),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Configuration/GenProduction/python/TOP-RunIIWinter15GS-00001-fragment.py nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

#pythia 8 config
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
from Configuration.Generator.Pythia8PowhegEmissionVetoSettings_cfi import *
process.generator = cms.EDFilter("Pythia8HadronizerFilter",
                                 maxEventsToPrint = cms.untracked.int32(0),
                                 pythiaPylistVerbosity = cms.untracked.int32(0),
                                 filterEfficiency = cms.untracked.double(1.0),
                                 pythiaHepMCVerbosity = cms.untracked.bool(False),
                                 comEnergy = cms.double(13000.),
                                 PythiaParameters = cms.PSet( pythia8CommonSettingsBlock,
                                                              pythia8CUEP8M1SettingsBlock,
                                                              pythia8PowhegEmissionVetoSettingsBlock,
                                                              processParameters = cms.vstring( 'POWHEG:nFinal = 2',
                                                                                               'TimeShower:mMaxGamma = 1.0',
                                                                                               'TimeShower:renormMultFac   = 1',
                                                                                               'TimeShower:factorMultFac   = 1',
                                                                                               'TimeShower:MEcorrections   = on'											      
                                                                                               ),
                                                              parameterSets = cms.vstring('pythia8CommonSettings',
                                                                                          'pythia8CUEP8M1Settings',
                                                                                          'pythia8PowhegEmissionVetoSettings',
                                                                                          'processParameters'
                                                                                          )
                                                              )
                                 )

if options.frag=='P':
	process.generator.PythiaParameters.processParameters.append('StringZ:usePetersonB = on')     
        print "Running Peterson with epsilonB = %s" % options.param
	process.generator.PythiaParameters.processParameters.append('StringZ:epsilonB = %f'%options.param)
if options.frag=='BL':
	process.generator.PythiaParameters.processParameters.append('StringZ:rFactB = %f'%options.param)
        print "Running Bowler-Lund with rFactB = %s" % options.param
print "max events = %s" % options.maxEvents

#pseudo-top config
from GeneratorInterface.RivetInterface.genParticles2HepMC_cfi import genParticles2HepMC
process.genParticles2HepMC = genParticles2HepMC.clone( genParticles = cms.InputTag("genParticles") )
process.load("TopQuarkAnalysis.TopEventProducers.producers.pseudoTop_cfi")
process.pseudoTop.jetMaxEta=cms.double(5.0)

#analysis config
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outputFile)
                                   )
process.load('TopLJets2015.TopAnalysis.bfragAnalysis_cfi')

# Path and EndPath definitions
process.ProductionFilterSequence = cms.Sequence(process.generator)
process.generation_step = cms.Path(process.pgen)
process.AnalysisSequence = cms.Path(process.genParticles2HepMC*process.pseudoTop*process.bfragAnalysis)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.AnalysisSequence,process.genfiltersummary_step,process.endjob_step)

# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 


