import optparse
import os,sys
import json
import commands

"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--json',      dest='inJson'  ,      help='json file with processed runs',      default=None,    type='string')
    parser.add_option('--puJson',    dest='puJson'  ,      help='pileup json file',      default='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_JSON_1.2_PRE_10-09-2015.txt',    type='string')
    (opt, args) = parser.parse_args()


    MINBIASXSEC={'nom':88000,'up':97680,'down':79200}
    PUBINS=50    
    for scenario in MINBIASXSEC:
        print scenario, 'xsec=',MINBIASXSEC[scenario]
        cmd='pileupCalc.py -i %s --inputLumiJSON %s --calcMode true --minBiasXsec %f --maxPileupBin %d --numPileupBins %s ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/PileupData_%s.root'%(opt.inJson,opt.puJson,MINBIASXSEC[scenario],PUBINS,PUBINS,scenario)
        commands.getstatusoutput(cmd)

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
