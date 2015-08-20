import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

runOnData=False #data/MC switch

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

if runOnData:
  #process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v1'
  process.GlobalTag.globaltag = 'GR_P_V56'
else:
  process.GlobalTag.globaltag = 'MCRUN2_74_V9'

#Number of events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))

#Input files
from UserCode.TopAnalysis.sp15.TT_TuneCUETP8M1_13TeV_powheg_pythia8_cfi import source as mc_events_source
from UserCode.TopAnalysis.Data_Mu_2015.data_mu_cfi import source as data_events_source

process.source=mc_events_source
process.source=data_events_source

if runOnData:
  process.source=data_events_source
else:
  process.source=mc_events_source

# Define the input source
#if runOnData:
#  fname = 'root://eoscms.cern.ch//store/data/Run2015B/JetHT/MINIAOD/PromptReco-v1/000/251/252/00000/263D331F-AF27-E511-969B-02163E012627.root'
#else:
#  fname = 'root://eoscms.cern.ch//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/60000/001C7571-0511-E511-9B8E-549F35AE4FAF.root'
#process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring([ fname ]))

#luminosity
import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList

if runOnData:
  process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt').getVLuminosityBlockRange()

#reduce verbosity
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 10

#tfileservice
process.TFileService = cms.Service("TFileService", fileName = cms.string('TT_TuneCUETP8M1_13TeV_Powheg.root'))

# Set up electron ID (VID framework)
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

# turn on VID producer, indicate data format  to be  DataFormat.AOD or DataFormat.MiniAOD, as appropriate 

dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
# define which IDs we want to produce
my_id_modules_el = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff','RecoEgamma.ElectronIdentification.Identification.mvaElectronID_PHYS14_PU20bx25_nonTrig_V1_cff']
#my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_PHYS14_PU20bx25_nonTrig_V1_cff']

#add them to the VID producer
for idmod in my_id_modules_el:
  setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#running sequence
process.load('UserCode.TopAnalysis.miniAnalyzer_cfi')
process.p = cms.Path(process.egmGsfElectronIDSequence*process.demo)

