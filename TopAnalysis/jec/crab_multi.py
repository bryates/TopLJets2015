from WMCore.Configuration import Configuration
config = Configuration()

# Usage: python crabConfig.py (to create jobs)
# Usage: ./multicrab -c status -w crab_projects to (check jobs status)
# Usage: ./multicrab -c status -w crab_projects/ -o "--long --sort=site"
# Usage: ./multicrab -c report -w crab_projects/ -o "--dbs=yes"
# Usage: ./multicrab -c resubmit -w crab_projects

# Common configuration
config.section_("General")
config.General.workArea        = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs    = True

config.section_("JobType")
config.JobType.pluginName   = 'Analysis' # PrivateMC
config.JobType.psetName     = 'UserCode/TopAnalysis/test/runMiniAnalyzer_cfg.py'
config.JobType.inputFiles   = ['Summer15_50nsV4_MC.db']
config.JobType.outputFiles  = ['minitree.root']

config.section_("Data")
config.Data.inputDBS        = 'global'    
config.Data.splitting       = 'FileBased' # EventBased, FileBased, LumiBased (1 lumi ~= 300 events)
config.Data.publication     = False


config.section_("Site")
config.Site.storageSite     = 'T2_IT_Legnaro'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException
    from multiprocessing import Process

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print hte.headers

    # dataset dependent configuration
    config.General.requestName  = 's13'
    config.Data.totalUnits      = 1
    config.Data.unitsPerJob     = 1
    config.Data.inputDataset = '/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

#    config.General.requestName  = 's11'
#    config.Data.totalUnits      = 1
#    config.Data.unitsPerJob     = 1
#    config.Data.inputDataset = '/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/MINIAODSIM'
#    p = Process(target=submit, args=(config,))
#    p.start()
#    p.join()
#
#    config.General.requestName  = 's12'
#    config.Data.totalUnits      = 1
#    config.Data.unitsPerJob     = 1
#    config.Data.inputDataset = '/ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/MINIAODSIM'
#    p = Process(target=submit, args=(config,))
#    p.start()
#    p.join()



