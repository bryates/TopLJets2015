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

runOnData= True #data/MC switch

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

if runOnData:
  process.GlobalTag.globaltag = 'GR_P_V56'
else:
  process.GlobalTag.globaltag = 'MCRUN2_74_V9'

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
  allowUnscheduled = cms.untracked.bool(True),
#  wantSummary = cms.untracked.bool(True)
)

#Number of events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

#Input files

# Signal & BKG
#from UserCode.TopAnalysis.sp15.TT_TuneCUETP8M1_13TeV_poWheg_pythia8_cfi import source as mc_events_source
#from UserCode.TopAnalysis.sp15.DYJetsToLL_M50_TuneCUETP8M1_13TeV_amcatnloFXFX_pythia8_cfi import source as mc_events_source
#from UserCode.TopAnalysis.sp15.WJetsToLNu_TuneCUETP8M1_13TeV_amcatnloFXFX_pythia8_cfi import source as mc_events_source
#from UserCode.TopAnalysis.sp15.ST_s_channel_4f_leptonDecays_13TeV_amcatnlo_pythia8_TuneCUETP8M1_cfi import source as mc_events_source
#from UserCode.TopAnalysis.sp15.ST_t_channel_4f_leptonDecays_13TeV_amcatnlo_pythia8_TuneCUETP8M1_cfi import source as mc_events_source
#from UserCode.TopAnalysis.sp15.ST_t_channel_antitop_4f_leptonDecays_13TeV_amcatnlo_pythia8_TuneCUETP8M1_cfi import source as mc_events_source
#from UserCode.TopAnalysis.sp15.ST_t_channel_top_4f_leptonDecays_13TeV_amcatnlo_pythia8_TuneCUETP8M1_cfi import source as mc_events_source
#from UserCode.TopAnalysis.sp15.ST_tW_antitop_5f_DS_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M1_cfi import source as mc_events_source
#from UserCode.TopAnalysis.sp15.ST_tW_top_5f_DS_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M1_cfi import source as mc_events_source
#from UserCode.TopAnalysis.sp15.QCD_Pt_80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8_cfi import source as mc_events_source
#from UserCode.TopAnalysis.sp15.QCD_Pt_120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8_cfi import source as mc_events_source
#from UserCode.TopAnalysis.sp15.QCD_Pt_170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8_cfi import source as mc_events_source
#from UserCode.TopAnalysis.sp15.QCD_Pt_300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8_cfi import source as mc_events_source

# Systematics
#from UserCode.TopAnalysis.sp15.TTJets_TuneCUETP8M1_13TeV_amcatnloFXFX_pythia8_cfi import source as mc_events_source
#from UserCode.TopAnalysis.sp15.TTJets_TuneCUETP8M1_13TeV_madgraphMLM_pythia8_cfi import source as mc_events_source
#from UserCode.TopAnalysis.sp15.TT_TuneCUETP8M1_13TeV_powheg_scaleup_pythia8_cfi import source as mc_events_source
#from UserCode.TopAnalysis.sp15.TT_TuneCUETP8M1_13TeV_powheg_scaledown_pythia8_cfi import source as mc_events_source
#from UserCode.TopAnalysis.sp15.TT_TuneCUETP8M1_mtop1695_13TeV_powheg_pythia8_cfi import source as mc_events_source
from UserCode.TopAnalysis.sp15.TT_TuneCUETP8M1_mtop1755_13TeV_powheg_pythia8_cfi import source as mc_events_source

# DATA
#from UserCode.TopAnalysis.DataMu15.singlemu_2015C_cfi import source as data_events_source
from UserCode.TopAnalysis.DataMu15.singlemu_2015D_cfi import source as data_events_source
#from UserCode.TopAnalysis.DataMu15.singleEle_28Aug2015_cfi import source as data_events_source
#from UserCode.TopAnalysis.DataMu15.singleEle_PromptReco2015C_cfi import source as data_events_source
#from UserCode.TopAnalysis.DataMu15.singleEle_2015D_cfi import source as data_events_source

process.source=mc_events_source
process.source=data_events_source

if runOnData:
  process.source=data_events_source
  outfilename='data_minitree.root'
else:
  process.source=mc_events_source
  outfilename='mc_minitree.root'

#luminosity
import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
if runOnData:
#  process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt').getVLuminosityBlockRange()
  process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-256869_13TeV_PromptReco_Collisions15_25ns_JSON.txt').getVLuminosityBlockRange()

#reduce verbosity
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#Tfileservice
process.TFileService = cms.Service("TFileService",

# Signal & BKG
#fileName = cms.string("TT_TuneCUETP8M1_SEP22.root")
#fileName = cms.string("DYJetsToLL_M50_SEP22.root")
#fileName = cms.string("WJetsToLNu_TuneCUETP8M1_SEP22.root")
#fileName = cms.string("ST_s_channel_SEP22.root")
#fileName = cms.string("ST_t_channel_SEP22.root")
#fileName = cms.string("ST_t_antitop_channel_SEP22.root")
#fileName = cms.string("ST_t_top_channel_SEP22.root")
#fileName = cms.string("ST_tW_antitop_SEP22.root")
#fileName = cms.string("ST_tW_top_SEP22.root")
#fileName = cms.string("QCD_Pt_80to120_EMEnriched_SEP22.root")
#fileName = cms.string("QCD_Pt_120to170_EMEnriched_SEP22.root")
#fileName = cms.string("QCD_Pt_170to300_EMEnriched_SEP22.root")
#fileName = cms.string("QCD_Pt_300toInf_EMEnriched_SEP22.root")

# Systematics
#fileName = cms.string("TTJets_TuneCUETP8M1_13TeV_amcatnloFXFX_SEP22.root")
#fileName = cms.string("TTJets_TuneCUETP8M1_13TeV_madgraphMLM_SEP22.root")
#fileName = cms.string("TT_TuneCUETP8M1_13TeV_powheg_scaleup_SEP22.root")
#fileName = cms.string("TT_TuneCUETP8M1_13TeV_powheg_scaledown_SEP22.root")
#fileName = cms.string("TT_TuneCUETP8M1_mtop1695_SEP22.root")
#fileName = cms.string("TT_TuneCUETP8M1_mtop1755_SEP17.root")

# DATA
#fileName = cms.string("singlemu_2015C_SEP22.root")
#fileName = cms.string("singlemu_2015C_SEP25.root")
fileName = cms.string("singlemu_2015D_SEP25.root")
#fileName = cms.string("singleEle_28Aug2015_SEP25.root")
#fileName = cms.string("singleEle_28Aug2015_SEP22.root")
#fileName = cms.string("singleEle_PromptReco2015C_SEP25.root")
#fileName = cms.string("singleEle_PromptReco2015C_SEP22.root")
#fileName = cms.string("singleEle_2015D_SEP25.root")
)

#JEC: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JecGlobalTag
from CondCore.DBCommon.CondDBSetup_cfi import *
import os
if runOnData:
  jecfile="Summer15_50nsV4_DATA"
else:
  jecfile="Summer15_50nsV4_MC"

#dBFile = os.path.expandvars("$CMSSW_BASE/src/UserCode/TopAnalysis/jec/"+jecfile+".db")
dBFile = os.path.expandvars(jecfile+".db")
process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
  #connect = cms.string( "sqlite_file://"+dBFile ),
  connect = cms.string( "sqlite:"+dBFile),
  toGet =  cms.VPSet(
  cms.PSet(
    record = cms.string("JetCorrectionsRecord"),
    tag = cms.string("JetCorrectorParametersCollection_"+jecfile+"_AK4PF"),
    label= cms.untracked.string("AK4PF")
    ),
  cms.PSet(
    record = cms.string("JetCorrectionsRecord"),
    tag = cms.string("JetCorrectorParametersCollection_"+jecfile+"_AK4PF"),
    label= cms.untracked.string("AK4PFchs")
    ),
  )
  )

process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')

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

