#!/usr/bin/env python

import os

#TUNES=[ ('up','BL',0.955), ('cuetp8m2t4','BL',0.855), ('down','BL',0.955) ]
#TUNES=[ ('sup','BL',0.865), ('scentral','BL',0.854), ('sdown','BL',0.845), ('cuetp8m2t4','BL',0.855) ]
#TUNES=[ ('scentral','BL',0.825), ('sdown','BL',0.655) ]
#TUNES=[ ('sdown','BL',0.655), ('up','BL',1.055),  ('uup','BL',1.000), ('uuup','BL',0.975), ('central','BL',0.955), ('ccentral','BL',0.900), ('cccentral','BL',0.875), ('cuetp8m2t4','BL',0.855), ('scentral','BL',0.825), ('ddown','BL',0.800), ('dddown','BL',0.775), ('down','BL',0.755)]
TUNES=[ ('1025','BL',1.025),  ('1075','BL',1.075), ('675','BL',0.675), ('625','BL',0.625), ('sdown','BL',0.655), ('up','BL',1.055),  ('uup','BL',1.000), ('uuup','BL',0.975), ('central','BL',0.955), ('ccentral','BL',0.900), ('cccentral','BL',0.875), ('cuetp8m2t4','BL',0.855), ('scentral','BL',0.825), ('ddown','BL',0.800), ('dddown','BL',0.775), ('down','BL',0.755) ]
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
        target = '%s/src/TopLJets2015/TopAnalysis/%s' % (cmsswBase,tag)
        condorFile = open(target,'w')
        condorFile.write('universe              = vanilla\n')
        condorFile.write('executable            = condor/cond_frag.sh\n')
        condorFile.write('arguments             = $(ClusterID) $(ProcId) %s %s %s %s\n' % (tag, MAXEVENTS, frag, param))
        condorFile.write('output                = condor/log/%s_$(ProcId).out\n' % tag)
        condorFile.write('error                 = condor/log/%s_$(ProcId).err\n' % tag)
        condorFile.write('log                   = condor/log/%s.log\n' % tag)
        condorFile.write('Rank                  = Memory >= 64\n')
        condorFile.write('Request_Memory        = 32 Mb\n')
        condorFile.write('+JobFlavour           = "workday"\n')
        condorFile.write('Should_Transfer_Files = NO\n')
        condorFile.write('queue')
        condorFile.close()
        os.system('condor_submit %s -batch-name %s' % (target,tag))
        os.system('rm %s' % (tag))

if __name__ == "__main__":
    main()




