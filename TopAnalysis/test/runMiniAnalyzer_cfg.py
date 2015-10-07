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

process = cms.Process("MiniAnalysis")

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = '74X_dataRun2_v2' if options.runOnData else '74X_mcRun2_asymptotic_v2'


# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
  allowUnscheduled = cms.untracked.bool(True),
  #  wantSummary = cms.untracked.bool(True)
)

#Number of events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/mc/RunIISpring15MiniAODv2/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v3/60000/FEEAA420-1E6A-E511-8E6D-00261894393D.root')
                            )


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
    process.source.fileNames = cms.untracked.vstring('/store/data/Run2015D/SingleMuon/MINIAOD/PromptReco-v3/000/258/155/00000/72A3374B-A76A-E511-AC17-02163E01388A.root')
    process.analysis.muTriggersToUse = cms.vstring('IsoMu18_v', 'IsoMu18_TriCentralPFJet50_40_30_v', 'IsoMu22_v', 'IsoMu22_TriCentralPFJet50_40_30_v', )
    process.analysis.elTriggersToUse = cms.vstring('Ele23_WPLoose_Gsf_v','Ele23_WPLoose_Gsf_TriCentralPFJet50_40_30','Ele27_WPLoose_Gsf_v','Ele27_WPLoose_Gsf_TriCentralPFJet50_40_30')

process.p = cms.Path(  process.egmGsfElectronIDSequence*process.analysis)

