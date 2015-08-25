import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')


from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

runOnData=True  #data/MC switch

use_nonTrigV1wp80=True #cut-based/mva-based elec_id switch
use_nonTrigV1wp90=False #cut-based/mva-based elec_id switch

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

if runOnData:
  process.GlobalTag.globaltag = 'GR_P_V56'
else:
  process.GlobalTag.globaltag = 'MCRUN2_74_V9'

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
  allowUnscheduled = cms.untracked.bool(True),
  wantSummary = cms.untracked.bool(True)
)

#Number of events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))

#Input files
from UserCode.TopAnalysis.sp15.TT_TuneCUETP8M1_13TeV_powheg_pythia8_cfi import source as mc_events_source
from UserCode.TopAnalysis.DataMu15.data_mu_cfi import source as data_events_source

process.source=mc_events_source
process.source=data_events_source

if runOnData:
  process.source=data_events_source
  outfilename='data_minitree.root'
else:
  process.source=mc_events_source
  outfilename='mc_minitree.root'

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

#Tfileservice
process.TFileService = cms.Service("TFileService",fileName = cms.string(outfilename))

#JEC: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JecGlobalTag
process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
process.patJetCorrFactorsReapplyJEC = process.patJetCorrFactorsUpdated.clone(
src = cms.InputTag("slimmedJets"),
  levels = ['L1FastJet','L2Relative', 'L3Absolute'],
  payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
process.patJetsReapplyJEC = process.patJetsUpdated.clone(
  jetSource = cms.InputTag("slimmedJets"),
  jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
  )
process.reapplyJEC = cms.Sequence( process.patJetCorrFactorsReapplyJEC + process.patJetsReapplyJEC)

# Set up electron ID (VID framework)
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

# turn on VID producer, indicate data format  to be  DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
# define which IDs we want to produce
my_id_modules_el = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff','RecoEgamma.ElectronIdentification.Identification.mvaElectronID_PHYS14_PU20bx25_nonTrig_V1_cff']

#add them to the VID producer
for idmod in my_id_modules_el:
  setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#running sequence
process.load('UserCode.TopAnalysis.miniAnalyzer_cfi')
process.p = cms.Path(
  process.reapplyJEC*
  process.egmGsfElectronIDSequence*
  process.demo
  )

