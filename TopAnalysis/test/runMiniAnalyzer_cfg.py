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
jecFile='Summer15_25nsV6_MC.db'
jecTag='Summer15_25nsV6_MC_AK4PFchs'
if options.runOnData : 
    jecLevels.append( 'L2L3Residual' )
    jecFile='Summer15_25nsV6_DATA.db'
    jecTag='Summer15_25nsV6_DATA_AK4PFchs'
customizeJetTools(process=process,jecLevels=jecLevels,jecFile=jecFile,jecTag=jecTag)


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = '74X_dataRun2_v2' if options.runOnData else '74X_mcRun2_asymptotic_v4'


# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
  allowUnscheduled = cms.untracked.bool(True),
  #  wantSummary = cms.untracked.bool(True)
)

#Number of events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

from TopLJets2015.TopAnalysis.storeTools import getEOSlslist
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(getEOSlslist(directory='/store/mc/RunIISpring15MiniAODv2/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v3/60000/')
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
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


#analyzer
process.load('TopLJets2015.TopAnalysis.miniAnalyzer_cfi')
if options.runOnData:
    print 'Adapting to run on data'
    process.source.fileNames = cms.untracked.vstring('/store/data/Run2015D/DoubleMuon/MINIAOD/PromptReco-v4/000/258/177/00000/00C9FA0D-576D-E511-B810-02163E011D21.root')
    process.analysis.muTriggersToUse = cms.vstring('IsoMu22_v','IsoTkMu20_v','IsoMu18_v', 'IsoMu22_v', 'IsoMu18_TriCentralPFJet50_40_30_v', 'IsoMu22_TriCentralPFJet50_40_30_v', )
    process.analysis.elTriggersToUse = cms.vstring('Ele23_WPLoose_Gsf_v','Ele23_WPLoose_Gsf_TriCentralPFJet50_40_30','Ele27_WPLoose_Gsf_v','Ele27_WPLoose_Gsf_TriCentralPFJet50_40_30')

if not options.saveTree:
    print 'Summary tree won\'t be saved'
    process.analysis.saveTree=cms.bool(False)

if options.runOnData:
    process.p = cms.Path( process.egmGsfElectronIDSequence
                          *process.analysis)
else:
    from TopLJets2015.TopAnalysis.GenTtbarCategorizer_cfi import defineGenTtbarCategorizerSequence
    defineGenTtbarCategorizerSequence(process)
    process.p = cms.Path( process.genTtbarCategorizerSequence
                          *process.egmGsfElectronIDSequence
                          *process.customizeJetToolsSequence
                          *process.analysis)
    
print process.p
