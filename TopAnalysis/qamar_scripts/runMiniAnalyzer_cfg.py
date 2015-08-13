import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

#from UserCode.TopAnalysis.csa14.TT_PU_v2_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.TT_PU_S14_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.TT_PU_S14_V5_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.TT_PU_S14_V7_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.TTJets_MG_PU20bx25_POSTLS170_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.TT_wjets_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.TT_w1234jets_V5_v1_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.TT_w1234jets_V5_v2_cfi import source as events_source 
#from UserCode.TopAnalysis.csa14.TT_DY_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_80_120_MuEnriched_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_120_170_MuEnriched_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_170_300_MuEnriched_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_300_470_MuEnriched_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_470_600_MuEnriched_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_600_800_MuEnriched_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_800_1000_MuEnriched_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_1000_MuEnriched_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_1000_1400_Tune4C_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_1400_1800_Tune4C_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_1800_Tune4C_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_1800_2400_Tune4C_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_2400_3200_Tune4C_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_3200_Tune4C_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_15_3000_Tune4C_flat_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_170_300_Tune4C_pythia8_cfi import source as events_source 
#from UserCode.TopAnalysis.csa14.QCD_300_470_Tune4C_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_470_600_Tune4C_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_600_800_Tune4C_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_800_1000_Tune4C_pythia8_cfi import source as events_source 
#from UserCode.TopAnalysis.csa14.QCD_80_120_Tune4C_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_120_170_Tune4C_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.TT_TBarToLeptons_s_channel_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.TT_TToLeptons_s_channel_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.TT_TBarToLeptons_t_channel_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.TT_TToLeptons_t_channel_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.T_tW_channel_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.Tbar_tW_channel_cfi import source as events_source

process.source=events_source

#reduce verbosity
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#tfileservice
process.TFileService = cms.Service("TFileService",
#        fileName = cms.string("miniAOD_MC_v2.root")
#        fileName = cms.string("miniAOD_norm_TT_PU_S14_V5.root")
#        fileName = cms.string("miniAOD_pre_TT_PU_S14_V7.root")
#        fileName = cms.string("pf_combinedleptons_TT_PU20bx25_v2.root")
#        fileName = cms.string("miniAOD_TT_PU20bx25_v2.root")
#        fileName = cms.string("miniAOD_TT_S14_PU40bx50.root")
#        fileName = cms.string('TTJets_MG_PU20bx25_POSTLS170.root')
#        fileName = cms.string("miniAOD_wjets_PU20bx25_new.root")
#        fileName = cms.string("w1234jets_v5_v1.root")
#        fileName = cms.string("w1234jets_v5_v2.root") 
#        fileName = cms.string("miniAOD_DY_PU20bx25_new.root")
#        fileName = cms.string("chargediso_QCD_80_120_MuEnriched_PU20bx25.root")
#        fileName = cms.string("chargediso_QCD_120_170_MuEnriched_PU20bx25.root")
#        fileName = cms.string("chargediso_QCD_170_300_MuEnriched_PU20bx25.root")
#        fileName = cms.string("chargediso_QCD_300_470_MuEnriched_PU20bx25.root")
#        fileName = cms.string("chargediso_QCD_470_600_MuEnriched_PU20bx25.root")
#        fileName = cms.string("chargediso_QCD_600_800_MuEnriched_PU20bx25.root")
#        fileName = cms.string("chargediso_QCD_800_1000_MuEnriched_PU20bx25.root")
#        fileName = cms.string("chargediso_QCD_1000_MuEnriched_PU20bx25.root")
#        fileName = cms.string("miniAOD_QCD_1000_1400_Tune4C.root")
#        fileName = cms.string("miniAOD_QCD_1400_1800_Tune4C.root")
#        fileName = cms.string("miniAOD_QCD_1800_Tune4C.root")
#        fileName = cms.string("miniAOD_QCD_1800_2400_Tune4C.root")
#        fileName = cms.string("miniAOD_QCD_2400_3200_Tune4C.root")
#        fileName = cms.string("miniAOD_QCD_3200_Tune4C.root")
#        fileName = cms.string("QCD_15_3000_Tune4C.root")
#        fileName = cms.string("miniAOD_QCD_170_300_Tune4C.root") 
#        fileName = cms.string("miniAOD_QCD_300_470_Tune4C.root")
#        fileName = cms.string("miniAOD_QCD_470_600_Tune4C.root")
#        fileName = cms.string("miniAOD_QCD_600_800_Tune4C.root")
#        fileName = cms.string("miniAOD_QCD_800_1000_Tune4C.root")
#        fileName = cms.string("miniAOD_QCD_80_120_Tune4C.root")
#        fileName = cms.string("miniAOD_QCD_120_170_Tune4C.root")
#        fileName = cms.string("miniAOD_TBarToLeptons_s_channel_new.root")
#        fileName = cms.string("miniAOD_TToLeptons_s_channel_new.root")
#        fileName = cms.string("miniAOD_TBarToLeptons_t_channel_new.root")
#        fileName = cms.string("miniAOD_TToLeptons_t_channel_new.root")
#        fileName = cms.string("miniAOD_T_tW_channel_new.root")
#        fileName = cms.string("miniAOD_Tbar_tW_channel_new.root")
)

#running sequence
#process.load('UserCode.TopAnalysis.myChargedPFJets_cfi')
#process.p = cms.Path(process.myChargedPFJets*process.demo)
process.load('UserCode.TopAnalysis.miniAnalyzer_cfi')
process.p = cms.Path(process.demo)


