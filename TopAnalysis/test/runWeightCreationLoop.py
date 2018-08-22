#!/usr/bin/env python

from WMCore.Configuration import Configuration
import os

#TUNES=[ ('up','BL',0.955), ('cuetp8m2t4','BL',0.855), ('down','BL',0.955) ]
#TUNES=[ ('sup','BL',0.865), ('scentral','BL',0.854), ('sdown','BL',0.845), ('cuetp8m2t4','BL',0.855) ]
#TUNES=[ ('scentral','BL',0.825), ('sdown','BL',0.655) ]
#TUNES=[ ('sdown','BL',0.655), ('up','BL',1.055),  ('uup','BL',1.000), ('uuup','BL',0.975), ('central','BL',0.955), ('ccentral','BL',0.900), ('cccentral','BL',0.875), ('cuetp8m2t4','BL',0.855), ('scentral','BL',0.825), ('ddown','BL',0.800), ('dddown','BL',0.775), ('down','BL',0.755)]
#TUNES=[ ('725','BL',0.725), ('700','BL',0.700) ]
#TUNES=[ ('1025','BL',1.025),  ('1075','BL',1.075), ('675','BL',0.675), ('625','BL',0.625), ('sdown','BL',0.655), ('up','BL',1.055),  ('uup','BL',1.000), ('uuup','BL',0.975), ('central','BL',0.955), ('ccentral','BL',0.900), ('cccentral','BL',0.875), ('scentral','BL',0.825), ('ddown','BL',0.800), ('dddown','BL',0.775), ('down','BL',0.755), ('725','BL',0.725), ('700','BL',0.700) ]
#TUNES=[ ('1025','BL',1.025), ('1075','BL',1.075), ('up','BL',1.055), ('uup','BL',1.000), ('uuup','BL',0.975), ('central','BL',0.955), ('925','BL',0.925), ('ccentral','BL',0.900), ('cccentral','BL',0.875), ('cuetp8m2t4','BL',0.855), ('scentral','BL',0.825), ('ddown','BL',0.800), ('dddown','BL',0.775), ('down','BL',0.755), ('725','BL',0.725), ('700','BL',0.700), ('675','BL',0.675), ('sdown','BL',0.655), ('625','BL',0.625) ]
#chi^2 fit values +/- 1 sigma
#TUNES=[ ('d0','BL',0.765), ('d0up','BL',0.788), ('d0down','BL',0.742), ('d0mu','BL',0.867), ('d0muup','BL',0.926), ('d0mudown','BL',0.808), ('jpsi','BL',0.798), ('jpsiup','BL',0.851), ('jpsidown','BL',0.745) ]
#chi^2 best fit +/- 1 sigma
#TUNES=[ ('fit','BL',0.80), ('fitup','BL',0.82), ('fitdown','BL',0.78) ]
TUNES=[ ('d0','BL',0.806), ('d0up','BL',0.781), ('d0down','BL',0.831), ('d0mu','BL',0.867), ('d0muup','BL',0.926), ('d0mudown','BL',0.808), ('jpsi','BL',0.798), ('jpsiup','BL',0.851), ('jpsidown','BL',0.745), ('fit','BL',0.80), ('fitup','BL',0.82), ('fitdown','BL',0.78), ('1025','BL',1.025), ('1075','BL',1.075), ('up','BL',1.055), ('uup','BL',1.000), ('uuup','BL',0.975), ('central','BL',0.955), ('925','BL',0.925), ('ccentral','BL',0.900), ('cccentral','BL',0.875), ('cuetp8m2t4','BL',0.855), ('scentral','BL',0.825), ('ddown','BL',0.800), ('dddown','BL',0.775), ('down','BL',0.755), ('725','BL',0.725), ('700','BL',0.700), ('675','BL',0.675), ('sdown','BL',0.655), ('625','BL',0.625), ('600','BL',0.600), ('555','BL',0.555) ]
#TUNES=[ ('sup','BL',0.900), ('scentral','BL',0.875), ('cuetp8m2t4','BL',0.855), ('sdown','BL',0.825)]
#TUNES=[ ('up','BL',1.055), ('central','BL',0.955), ('cuetp8m2t4','BL',0.855), ('down','BL',0.755)]
#TUNES=[ ('up','BL',1.079), ('central','BL',0.8949), ('cuetp8m2t4','BL',0.855), ('down','BL',0.6981)]
#TUNES=[ ('up','BL',1.079), ('central','BL',0.8949), ('cuetp8m2t4','BL',0.855), ('down','BL',0.6981),('Peterson','P',0.003271)]
#MAXEVENTS=500
#MAXEVENTS=500000
MAXEVENTS=-1
CFG='${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runBFragmentationAnalyzer_cfg.py'

def main():

    cmsswBase=os.environ['CMSSW_BASE']

    #run the analyzer
    for tag,frag,param in TUNES:
        #os.system('cmsRun %s maxEvents=%d frag=%s param=%f outputFile=xb_%s.root'%(CFG,MAXEVENTS,frag,param,tag))
        #target = '%s/src/TopLJets2015/TopAnalysis/%s' % (cmsswBase,tag)
        #condorFile = open(target,'w')
        #condorFile.write('universe              = vanilla\n')
        #condorFile.write('executable            = condor/cond_frag.sh\n')
        #condorFile.write('arguments             = $(ClusterID) $(ProcId) %s %s %s %s\n' % (tag, MAXEVENTS, frag, param))
        #condorFile.write('output                = condor/log/%s_$(ProcId).out\n' % tag)
        #condorFile.write('error                 = condor/log/%s_$(ProcId).err\n' % tag)
        #condorFile.write('log                   = condor/log/%s.log\n' % tag)
        #condorFile.write('Rank                  = Memory >= 64\n')
        #condorFile.write('Request_Memory        = 32 Mb\n')
        #condorFile.write('+JobFlavour           = "workday"\n')
        #condorFile.write('Should_Transfer_Files = NO\n')
        #condorFile.write('queue')
        #condorFile.close()
        #os.system('condor_submit %s -batch-name %s' % (target,tag))
        #os.system('rm %s' % (tag))

        workDir='grid'
        crabConfigFile=workDir+'/'+tag+'_cfg.py'
        config_file=open(crabConfigFile,'w')
        config_file.write('from WMCore.Configuration import Configuration\n')
        config_file.write('import os\n')
        config_file.write('config = Configuration()\n')
        config_file.write('\n')
        
        config_file.write('config.section_("General")\n')
        config_file.write('config.General.requestName = "%s"\n' % tag)
        config_file.write('config.General.workArea = "grid"\n')
        config_file.write('config.General.transferOutputs=True\n')
        config_file.write('\n')
        config_file.write('config.section_("JobType")\n')
        config_file.write('config.JobType.pluginName = "Analysis"\n')
        config_file.write('config.JobType.psetName = "/afs/cern.ch/user/b/byates/CMSSW_8_0_26/src/TopLJets2015/TopAnalysis/test/runBFragmentationAnalyzer_cfg.py"\n')
        config_file.write('config.JobType.disableAutomaticOutputCollection = False\n')
        config_file.write('config.JobType.pyCfgParams = [\'frag=%s\',\'param=%f\',\'maxEvents=%d\']' % (frag,param,MAXEVENTS))
        config_file.write('\n')
        config_file.write('\n')
        config_file.write('config.section_("Data")\n')
        config_file.write('config.Data.inputDataset = "/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM"\n')
        config_file.write('config.Data.inputDBS = "global"\n')
        config_file.write('config.Data.splitting = "Automatic"\n')
        #config_file.write('config.Data.splitting = "FileBased"\n')
        #config_file.write('config.Data.unitsPerJob = 4\n')
        config_file.write('config.Data.publication = True\n')
        config_file.write('config.Data.ignoreLocality = False\n')
        config_file.write('config.Data.outLFNDirBase = \'/store/group/phys_top/byates/bfrag/%s\'\n' % tag)
        config_file.write('\n')
        config_file.write('config.section_("Site")\n')
        config_file.write('config.Site.storageSite = "T2_CH_CERN"')
        config_file.close()

        os.system('alias crab=\'/cvmfs/cms.cern.ch/crab3/crab-env-bootstrap.sh\' && crab submit -c %s' % crabConfigFile )

if __name__ == "__main__":
    main()




