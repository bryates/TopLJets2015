import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#from Configuration.AlCa.GlobalTag import GlobalTag
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'MCRUN2_74_V9A::All' #for Simulation

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000))

#from UserCode.TopAnalysis.sp15.TT_TuneCUETP8M1_13TeV_Powheg_cfi import source as events_source
from UserCode.TopAnalysis.sp15.TT_new_cfi import source as events_source

process.source=events_source

#reduce verbosity
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 10

#tfileservice
process.TFileService = cms.Service("TFileService",
  fileName = cms.string('TT_TuneCUETP8M1_13TeV_Powheg.root')
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

#============JEC
process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",
  DBParameters = cms.PSet(
  messageLevel = cms.untracked.int32(1)
  ),
  timetype = cms.string('runnumber'),
  toGet = cms.VPSet(
  cms.PSet(
    record = cms.string('JetCorrectionsRecord'),
    tag    = cms.string('JetCorrectorParametersCollection_Summer15_50nsV2_MC_AK4PFchs'),
    # tag    = cms.string('JetCorrectorParametersCollection_Summer12_V3_MC_AK5PF'),
    label  = cms.untracked.string('AK4PFchs')
    ),
   ## here you add as many jet types as you need
   ## note that the tag name is specific for the particular sqlite file 
   ), 
   connect = cms.string('sqlite:Summer15_50nsV2_MC.db')
   # uncomment above tag lines and this comment to use MC JEC
   # connect = cms.string('sqlite:Summer12_V7_MC.db')
)
## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
process.patJetCorrFactorsReapplyJEC = process.patJetCorrFactorsUpdated.clone(
  src = cms.InputTag("slimmedJets"),
  levels = ['L1FastJet', 'L2Relative', 'L3Absolute'],
  payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
process.patJetsReapplyJEC = process.patJetsUpdated.clone(
  jetSource = cms.InputTag("slimmedJets"),
  jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
  )
process.JEC = cms.Sequence( process.patJetCorrFactorsReapplyJEC + process. patJetsReapplyJEC )
#============JEC

#running sequence
process.load('UserCode.TopAnalysis.miniAnalyzer_cfi')
process.p = cms.Path(process.egmGsfElectronIDSequence*process.demo)
process.p += cms.Sequence( process.patJetCorrFactorsReapplyJEC + process. patJetsReapplyJEC )

