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
options.register('savePF', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 'save PF candidates'
                 )
options.parseArguments()

process = cms.Process("MiniAnalysis")


# Load the standard set of configuration modules
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')


# global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v8' if options.runOnData else '80X_mcRun2_asymptotic_2016_v3')

#message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


# set input to process
from TopLJets2015.TopAnalysis.Compressed_T2tt_cfi import *
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root'),
                            #fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/F8B2494E-1ABF-E611-8FB5-70106F4A9248.root'),
                            #fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/50000/F254D922-4BE3-E611-AE37-0CC47A78A426.root'),
                            #fileNames = cms.untracked.vstring(Compressed_T2tt_2),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck') 
                            )
if options.runOnData:
    process.source.fileNames = cms.untracked.vstring('/store/data/Run2016G/MuonEG/MINIAOD/23Sep2016-v1/50000/A2AAED1F-F18F-E611-A51D-002590D0AFBA.root')
    #/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/06277EC1-181A-E611-870F-02163E0145E5.root')

#this make the process crash ?!
#if options.inputDir!='': 
#    from TopLJets2015.TopAnalysis.storeTools import getEOSlslist 
#    print 'Will process files from',options.inputDir
#    process.source.fileNames=cms.untracked.vstring(getEOSlslist(directory=options.inputDir))


######### Skim Filter
process.selectedMuons = cms.EDFilter("CandPtrSelector",
                                     src = cms.InputTag("slimmedMuons"),
                                     cut = cms.string("pt>9.8 && abs(eta)<2.4"))

process.selectedElectrons = cms.EDFilter("CandPtrSelector",
                                         src = cms.InputTag("slimmedElectrons"),
                                         cut = cms.string("pt>9.8 && abs(eta)<2.5"))

process.selectedJets = cms.EDFilter("CandPtrSelector",
                                         src = cms.InputTag("slimmedJets"),
                                         cut = cms.string("pt>20 && abs(eta)<2.5"))

process.allLeps = cms.EDProducer("CandViewMerger",
                                 src = cms.VInputTag(
                                                        cms.InputTag("selectedElectrons"),
                                                        cms.InputTag("selectedMuons")))

process.countLeps = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("allLeps"),
                                 minNumber = cms.uint32(1))

process.countJets = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("selectedJets"),
                                 minNumber = cms.uint32(2))

process.preYieldFilter = cms.Sequence(process.selectedMuons+process.selectedElectrons+process.allLeps+process.countLeps+process.selectedJets+process.countJets)

#analysis
process.load('TopLJets2015.TopAnalysis.miniAnalyzer_cfi')
if not options.saveTree:
    print 'Summary tree won\'t be saved'
    process.analysis.saveTree=cms.bool(False)
if not options.savePF:
    print 'Summary PF info won\'t be saved'
    process.analysis.savePF=cms.bool(False)

#pseudo-top
if not options.runOnData:
    process.load('TopQuarkAnalysis.TopEventProducers.producers.pseudoTop_cfi')
    process.pseudoTop.leptonMinPt=cms.double(20)
    process.pseudoTop.leptonMaxEta=cms.double(2.5)
    process.pseudoTop.jetMaxEta=cms.double(5.0)

# Set up electron ID (VID framework)
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff'] 
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#jet energy corrections
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
from TopLJets2015.TopAnalysis.customizeJetTools_cff import *
jecLevels=['L1FastJet','L2Relative','L3Absolute']
jecFile='Summer16_23Sep2016V4_MC.db'
jecTag='Summer16_23Sep2016V4_MC_AK4PFchs'
if options.runOnData : 
    print 'Warning we\'re still using Spring16 MC corrections for data - to be updated'
    jecLevels.append( 'L2L3Residual' )
    jecFile='Summer16_23Sep2016AllV4_DATA.db'
    jecTag='Summer16_23Sep2016AllV4_DATA_AK4PFchs'
customizeJetTools(process=process,jecLevels=jecLevels,jecFile=jecFile,jecTag=jecTag)

#tfile service
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('histo.root')
                                   )

if options.runOnData:
    process.p = cms.Path(process.preYieldFilter*process.egmGsfElectronIDSequence*process.customizeJetToolsSequence*process.analysis)
else:
    process.p = cms.Path(process.preYieldFilter*process.egmGsfElectronIDSequence*process.customizeJetToolsSequence*process.pseudoTop*process.analysis)
