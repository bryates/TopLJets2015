import os
import sys
import optparse
from TopLJets2015.TopAnalysis.storeTools import getEOSlslist
from subprocess import Popen, PIPE

"""
steer the script
"""
def main():

    cmsswBase=os.environ['CMSSW_BASE']
    eos_cmd = '/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select'

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inDir',      dest='inDir',       help='input directory with files',               default=None,   type='string')
    parser.add_option('-o', '--outDir',     dest='outDir',      help='output directory with files',              default=None,   type='string')
    parser.add_option('-c', '--cleanup',    dest='cleanup',     help='removes original crab directory',          default =False, action='store_true')
    parser.add_option('-q', '--queue',      dest='queue',       help='batch queue',                              default='2nd',  type='string')
    parser.add_option(      '--only',       dest='only',        help='only this tag',                            default=None,   type='string')
    (opt, args) = parser.parse_args()

    Popen([eos_cmd, ' -b fuse mount', 'eos'],stdout=PIPE).communicate()

    #prepare output directory
    if opt.outDir is None: opt.outDir=opt.inDir
    Popen([eos_cmd, 'mkdir', '/eos/cms/'+opt.outDir],stdout=PIPE).communicate()

    dset_list=getEOSlslist(directory=opt.inDir,prepend='')
    for dset in dset_list:
        dsetname=dset.split('/')[-1]

        pub_list=getEOSlslist(directory=dset,prepend='')
        for pubDir in pub_list:

            if not 'crab' in pubDir:
                print 'Ambiguity found @ <publication-name> for <primary-dataset>=%s , bailing out'%dsetname
                continue
            pub=pubDir.split('/crab_')[-1]
            onlyList=[]

            if opt.only:
                if pub!=opt.only: #or pub not in opt.only:
                    continue


            #if 'Data13TeV' in pub : continue

            localMerge='python scripts/checkProductionIntegrity.py -i %s -o %s --nocheck --only %s'%(opt.inDir,opt.outDir,pub)
            if opt.cleanup :
                localMerge += ' --cleanup'
            cmd='bsub -q %s %s/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh \"%s\"' % (opt.queue,cmsswBase,localMerge)
            os.system(cmd)

    #Popen([eos_cmd, ' -b fuse umount', 'eos'],stdout=PIPE).communicate()

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

