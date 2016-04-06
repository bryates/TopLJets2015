import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('runOnData', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Run this on real data"
                 )
options.register('outFilename', 'MiniEvents.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Output file name"
                 )
options.register('inputDir', '',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "input directory with files to process"
                 )
options.register('saveTree', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "save summary tree"
                 )
options.parseArguments()

process = cms.Process("MiniAnalysis")

# Load the standard set of configuration modules
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#jet energy corrections
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
from TopLJets2015.TopAnalysis.customizeJetTools_cff import *
jecLevels=['L1FastJet','L2Relative','L3Absolute']
jecFile='Fall15_25nsV2_MC.db'
jecTag='Fall15_25nsV2_MC_AK4PFchs'
if options.runOnData : 
    jecLevels.append( 'L2L3Residual' )
    jecFile='Fall15_25nsV2_DATA.db'
    jecTag='Fall15_25nsV2_DATA_AK4PFchs'
customizeJetTools(process=process,jecLevels=jecLevels,jecFile=jecFile,jecTag=jecTag)


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = '76X_dataRun2_v15' if options.runOnData else '76X_mcRun2_asymptotic_v12'


# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
  allowUnscheduled = cms.untracked.bool(True),
  #  wantSummary = cms.untracked.bool(True)
)

#Number of events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

from TopLJets2015.TopAnalysis.storeTools import getEOSlslist
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        #getEOSlslist(directory='/store/mc/RunIIFall15DR76/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/25nsFlat10to25TSG_76X_mcRun2_asymptotic_v11_ext3-v1/30000/')
        getEOSlslist(directory='/store/mc/RunIIFall15MiniAODv2/TT_TuneEE5C_13TeV-powheg-herwigpp/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/10000')
        )
                            )
if options.inputDir!='':  
    print 'Will process files from',options.inputDir
    process.source.fileNames=cms.untracked.vstring(getEOSlslist(directory=options.inputDir))
                                
#reduce verbosity
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#Tfileservice
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outFilename))

# Set up electron ID (VID framework)
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff'] 
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


#PUPPI
process.load('CommonTools/PileupAlgos/Puppi_cff')
process.puppi.candName = cms.InputTag('packedPFCandidates')
process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')

#analyzer
process.load('TopLJets2015.TopAnalysis.miniAnalyzer_cfi')
if options.runOnData:
    print 'Adapting to run on data'
    process.source.fileNames = cms.untracked.vstring('/store/data/Run2015D/DoubleMuon/MINIAOD/16Dec2015-v1/10000/F8146472-37A8-E511-9C12-0CC47A4C8E70.root')

if not options.saveTree:
    print 'Summary tree won\'t be saved'
    process.analysis.saveTree=cms.bool(False)



if options.runOnData:
    process.p = cms.Path( process.puppi
                          *process.egmGsfElectronIDSequence
                          *process.customizeJetToolsSequence
                          *process.analysis
                          )
else:
    #pseudo-top
    process.load('TopQuarkAnalysis.TopEventProducers.producers.pseudoTop_cfi')
    process.pseudoTop.leptonMinPt=cms.double(20)
    process.pseudoTop.leptonMaxEta=cms.double(2.5)
    process.pseudoTop.jetMaxEta=cms.double(5.0)

    process.p = cms.Path( process.puppi
                          *process.egmGsfElectronIDSequence
                          *process.customizeJetToolsSequence
                          *process.pseudoTop
                          *process.analysis
                          )

#save in edm format
#process.Out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string("/tmp/psilva/test.root")
#                               )    

#process.end = cms.EndPath(process.Out)

print process.p
