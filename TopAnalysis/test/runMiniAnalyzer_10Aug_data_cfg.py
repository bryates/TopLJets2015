import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_P_V56::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))

#from UserCode.TopAnalysis.sp15.TT_TuneCUETP8M1_13TeV_Powheg_cfi import source as events_source
#from UserCode.TopAnalysis.sp15.TT_new_cfi import source as events_source
#process.source=events_source

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
# 'root://cms-xrd-global.cern.ch//store/data/Run2015B/JetHT/MINIAOD/PromptReco-v1/000/251/163/00000/F05CF208-A026-E511-85F4-02163E011B15.root'
 'root://cms-xrd-global.cern.ch//store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/252/116/00000/7E03BD9B-7730-E511-BA13-02163E011976.root'
 )
)

#reduce verbosity
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 10

#tfileservice
process.TFileService = cms.Service("TFileService",
  fileName = cms.string('Data_Run2015B.root')
)

# Set up electron ID (VID framework)
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

# turn on VID producer, indicate data format  to be  DataFormat.AOD or DataFormat.MiniAOD, as appropriate 

dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
# define which IDs we want to produce
my_id_modules_el = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff']

#add them to the VID producer
for idmod in my_id_modules_el:
  setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#running sequence
process.load('UserCode.TopAnalysis.miniAnalyzer_cfi')
process.p = cms.Path(process.egmGsfElectronIDSequence*process.demo)


