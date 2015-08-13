from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'TTBar'
config.General.workArea = 'crab_projects'

config.General.transferOutputs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/w/wajid/qamar_ttbar/CMSSW_7_4_5/src/UserCode/TopAnalysis/test/runMiniAnalyzer_10Aug_mc_cfg.py'

config.section_("Data")
config.Data.inputDataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v4/MINIAODSIM'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = 'T2_IT_Legnaro'
